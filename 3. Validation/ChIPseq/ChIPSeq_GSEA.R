# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(GEOquery)
library(rtracklayer)
library(tidyverse)
library(biomaRt)
library(Gviz)
library(readxl)
library(chipseq)
library(BSgenome.Mmusculus.UCSC.mm9)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/ChIPseq"))

# Binning function
BinChIPseq = function(reads, bins){
  mcols(bins)$score = countOverlaps(bins, reads) 
  return(bins) 
}

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Get sequence info
genome = BSgenome.Mmusculus.UCSC.mm9
si = seqinfo(genome)
si = si[ paste0('chr', c(1:19))]

# Import BW files
MeCP2_rep1 <- import.bw("Data/GSE71126_RAW/GSM1827604_MeCP2_ChIP_WT_rep1.extend200.coverage.mm9.bw")
MeCP2_rep2 <- import.bw("Data/GSE71126_RAW/GSM1827605_MeCP2_ChIP_WT_rep2.extend200.coverage.mm9.bw")
Input_rep1 <- import.bw("Data/GSE71126_RAW/GSM1827606_Input_WT_rep1.extend200.coverage.mm9.bw")
Input_rep2 <- import.bw("Data/GSE71126_RAW/GSM1827607_Input_WT_rep2.extend200.coverage.mm9.bw")

# Select chr1 - chr19
MeCP2_rep1_fil <- MeCP2_rep1[!(seqnames(MeCP2_rep1) %in% c("chrM", "chrY", "chrX"))]
rm(MeCP2_rep1)
MeCP2_rep2_fil <- MeCP2_rep2[!(seqnames(MeCP2_rep2) %in% c("chrM", "chrY", "chrX"))]
rm(MeCP2_rep2)
Input_rep1_fil <- Input_rep1[!(seqnames(Input_rep1) %in% c("chrM", "chrY", "chrX"))]
rm(Input_rep1)
Input_rep2_fil <- Input_rep2[!(seqnames(Input_rep2) %in% c("chrM", "chrY", "chrX"))]
rm(Input_rep2)


#******************************************************************************#
# MeCP2 enrichment for the GO terms
#******************************************************************************#

# Load data
load("Data/egs.RData")
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_ENSEMBL_Mm.RData"))
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

pvalues <- rep(NA, length(terms_ordered))
meanDiff <- rep(NA, length(terms_ordered))
Nsel <- rep(NA, length(terms_ordered))
Nnonsel <- rep(NA, length(terms_ordered))

for (terms in 1:length(terms_ordered)){
  
  # Get selected and non-selected genes
  ID <- GOannotation$ID[GOannotation$Name %in% terms_ordered[terms]]
  sel_genes <- GOgenes[[ID]]
  egs_sel <- egs[egs$ensembl_gene_id %in% sel_genes,]
  egs_nonsel <- egs[!(egs$ensembl_gene_id %in% sel_genes),]
  
  # Settings
  range <- 3000 # Range around TSS
  lo <- 60      # length out (bin size)
  
  # Get counts around TSS of selected genes
  tiles_sel = sapply(1:nrow(egs_sel), function(i)
    if(egs_sel$strand[i] == "1" )
      egs_sel$TSS[i] + seq(-1*range, range-100, length.out=lo)
    else
      egs_sel$TSS[i] + seq(range-100, -1*range, length.out=lo))
  
  tiles_sel = GRanges(tilename = paste(rep(egs_sel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                      seqnames = Rle(rep(paste0('chr', egs_sel$chromosome_name), each=lo)), 
                      ranges = IRanges(start = as.vector(tiles_sel),
                                       width = 100),
                      strand = Rle(rep("*", length(as.vector(tiles_sel)))),
                      seqinfo=si)
  
  selLoci = countOverlaps(tiles_sel, MeCP2_rep1_fil) + 
    countOverlaps(tiles_sel, MeCP2_rep2_fil)
  
  selLoci_matrix = matrix(selLoci, nrow=nrow(egs_sel), 
                          ncol=lo, byrow=TRUE)
  
  
  
  # Get counts around TSS of non-selected genes
  tiles_nonsel = sapply(1:nrow(egs_nonsel), function(i)
    if(egs_nonsel$strand[i] == "1" )
      egs_nonsel$TSS[i] + seq(-1*range, range-100, length.out=lo)
    else
      egs_nonsel$TSS[i] + seq(range-100, -1*range, length.out=lo))
  
  tiles_nonsel = GRanges(tilename = paste(rep(egs_nonsel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                         seqnames = Rle(rep(paste0('chr', egs_nonsel$chromosome_name), each=lo)), 
                         ranges = IRanges(start = as.vector(tiles_nonsel),
                                          width = 100),
                         strand = Rle(rep("*", length(as.vector(tiles_nonsel)))),
                         seqinfo=si)
  
  nonselLoci = countOverlaps(tiles_nonsel, MeCP2_rep1_fil) + 
    countOverlaps(tiles_nonsel, MeCP2_rep2_fil)
  
  nonselLoci_matrix = matrix(nonselLoci, nrow=nrow(egs_nonsel), 
                             ncol=lo, byrow=TRUE)
  
  
  # Test whether the mean enrichment is different for imprinted vs non-imprinted genes
  statistics <- t.test(rowMeans(nonselLoci_matrix), rowMeans(selLoci_matrix))
  pvalues[terms] <- statistics$p.value
  meanDiff[terms] <- statistics$estimate[2] - statistics$estimate[1]
  Nsel[terms] <- nrow(selLoci_matrix)
  Nnonsel[terms] <- nrow(nonselLoci_matrix)
}

# Combine results into data frame for plotting
plotDF <- data.frame(TermName = factor(terms_ordered, levels = terms_ordered),
                     Pvalue = pvalues,
                     meanDiff = meanDiff,
                     Nsel = Nsel,
                     Nnonsel = Nnonsel)

# Prepare data for plotting
plotDF <- inner_join(plotDF,GOannotation, by = c("TermName" = "Name"))
plotDF$TermName <- factor(plotDF$TermName, levels = terms_ordered)
plotDF <- arrange(plotDF, by = TermName)
plotDF$Description <- firstup(plotDF$Description)
plotDF$Description <- factor(plotDF$Description, levels = unique(plotDF$Description))

plotDF$Group <- NA
plotDF$Group[1:6] <- "  "
plotDF$Group[7:18] <- " "
plotDF$Group[19:26] <- ""

# Make plot
p <- ggplot(plotDF) +
  geom_bar(aes(x = -log10(Pvalue) * sign(meanDiff), y = Description),
           stat  ="identity", position = position_dodge(), 
           fill = "grey", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = -log10(0.05), linewidth = 1, linetype = "dashed") +
  ylab(NULL) +
  xlab(expression("signed"~-log[10]~"P value")) +
  facet_grid(rows = vars(Group), scale = "free", space = "free") +
  theme_bw()

panel_colors <- c("#FEC44F","#74C476","#4292C6")

# convert to grob
gp <- ggplotGrob(p) # where p is the original ggplot object

# assign the colors to the panels
for(i in 1:3){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}
grid::grid.draw(gp)

# Save plot
ggsave(gp, file = "ChIPSeq_GO.png",
       width = 6, height = 6)



#******************************************************************************#
# Ameboidal-type cell migration
#******************************************************************************#

#==============================================================================#
# Heatmap plot
#==============================================================================#

# Get selected and non-selected genes
ID <- GOannotation$ID[GOannotation$Name %in% "ameboidal-type cell migration (GO:0001667)"]
sel_genes <- GOgenes[[ID]]
egs_sel <- egs[egs$ensembl_gene_id %in% sel_genes,]
egs_nonsel <- egs[!(egs$ensembl_gene_id %in% sel_genes),]

# Settings
range <- 3000
lo <- 60

# Get counts around TSS of selected genes
tiles_sel = sapply(1:nrow(egs_sel), function(i)
  if(egs_sel$strand[i] == "1" )
    egs_sel$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_sel$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_sel = GRanges(tilename = paste(rep(egs_sel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                    seqnames = Rle(rep(paste0('chr', egs_sel$chromosome_name), each=lo)), 
                    ranges = IRanges(start = as.vector(tiles_sel),
                                     width = 100),
                    strand = Rle(rep("*", length(as.vector(tiles_sel)))),
                    seqinfo=si)

# Get bin count
selLoci = countOverlaps(tiles_sel, MeCP2_rep1_fil) + 
  countOverlaps(tiles_sel, MeCP2_rep2_fil)

selLoci_matrix = matrix(selLoci, nrow=nrow(egs_sel), 
                        ncol=lo, byrow=TRUE)

rownames(selLoci_matrix) <- egs_sel$ensembl_gene_id
colnames(selLoci_matrix) <- seq(-1*range, range, length.out=lo)

# Make data frame for plotting
plotMatrix <- gather(as.data.frame(selLoci_matrix))
plotMatrix$Gene <- rep(rownames(selLoci_matrix), ncol(selLoci_matrix))
plotMatrix$key <- as.numeric(plotMatrix$key)
plotMatrix$Gene <- factor(plotMatrix$Gene, levels = names(sort(rowSums(selLoci_matrix))))

# Make heatmap plot
p <- ggplot(plotMatrix) +
  geom_tile(aes(x = key, y = Gene, fill = value)) +
  scale_fill_viridis_c(option = "plasma", limits = c(0,70),
                       na.value = "yellow") +
  xlab("Distance from TSS") +
  ggtitle("Ameboidal-type cell migration") +
  labs(fill = "Enrichment") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# save plot
ggsave(p, file = "Heatmap_Ameboidal-type cell migration.png",
       width = 6, height = 6)


#==============================================================================#
# Density plot
#==============================================================================#

# Get counts around TSS of non-selected genes
tiles_nonsel = sapply(1:nrow(egs_nonsel), function(i)
  if(egs_nonsel$strand[i] == "1" )
    egs_nonsel$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_nonsel$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_nonsel = GRanges(tilename = paste(rep(egs_nonsel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                       seqnames = Rle(rep(paste0('chr', egs_nonsel$chromosome_name), each=lo)), 
                       ranges = IRanges(start = as.vector(tiles_nonsel),
                                        width = 100),
                       strand = Rle(rep("*", length(as.vector(tiles_nonsel)))),
                       seqinfo=si)

# Get bin counts
nonselLoci = countOverlaps(tiles_nonsel, MeCP2_rep1_fil) + 
  countOverlaps(tiles_nonsel, MeCP2_rep2_fil)
nonselLoci_matrix = matrix(nonselLoci, nrow=nrow(egs_nonsel), 
                           ncol=lo, byrow=TRUE)

# prepare data for plot
plotDF <- data.frame(value = c(colMeans(selLoci_matrix),
                               colMeans(nonselLoci_matrix)),
                     location = rep(seq(-1*range, range, length.out=lo),2),
                     group = c(rep("Selected",lo),
                               rep("Non-selected",lo))
)


# Get mean counts around TSS for different permutations
nPerm <- 100
plotDF_permutation  <- NULL
set.seed(123)
for (p in 1:nPerm){
  
  # Select random genes
  egs_permutation <- egs[sample(1:nrow(egs),nrow(selLoci_matrix)),]
  
  # Get tiles around TSS
  range <- 3000
  lo <- 60
  tiles_permutation = sapply(1:nrow(egs_permutation), function(i)
    if(egs_permutation$strand[i] == "1" )
      egs_permutation$TSS[i] + seq(-1*range, range-100, length.out=lo)
    else
      egs_permutation$TSS[i] + seq(range-100, -1*range, length.out=lo))
  
  tiles_permutation= GRanges(tilename = paste(rep(egs_permutation$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                             seqnames = Rle(rep(paste0('chr', egs_permutation$chromosome_name), each=lo)), 
                             ranges = IRanges(start = as.vector(tiles_permutation),
                                              width = 100),
                             strand = Rle(rep("*", length(as.vector(tiles_permutation)))),
                             seqinfo=si)
  
  # Get bin counts
  permutationLoci = countOverlaps(tiles_permutation, MeCP2_rep1_fil) + 
    countOverlaps(tiles_permutation, MeCP2_rep2_fil)
  permutationLoci_matrix = matrix(permutationLoci, nrow=nrow(egs_permutation), 
                                  ncol=lo, byrow=TRUE)
  
  # Put results into dataframe
  temp <- data.frame(value = colMeans(permutationLoci_matrix),
                     location = rep(seq(-1*range, range, length.out=lo),2),
                     group = p)
  
  # Combine with results from previous permutation
  plotDF_permutation <- rbind.data.frame(plotDF_permutation,temp)
}

# Save plotting data frames
plotDF$group <- ifelse(plotDF$group == "Selected", "Ameboidal-type cell migration", "Other")
save(plotDF, plotDF_permutation, 
     file = "Data/ChIPSeqData_Ameboidal-type cell migration.RData")

# Load plotting data frames (if needed)
load("Data/ChIPSeqData_Ameboidal-type cell migration.RData")

# Make density plot
colors <- c("#4292C6", "#737373")
p <- ggplot() +
  geom_line(data = plotDF_permutation, 
            aes(x = location, y = value, group = group),
            linewidth = 0.5, color = "#D9D9D9") +
  geom_line(data = plotDF,
            aes(x = location, y = value, group = group, color = group),
            linewidth = 1.5) +
  ylab("Mean enrichment") +
  xlab("Distance from the TSS") +
  ylim((c(10,40)))+
  labs(color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_color_manual(values = colors)

# Save plot
ggsave(p, file = "Density_Ameboidal-type cell migration.png",
       width = 6, height = 4)



#******************************************************************************#
# RNA splicing, via transesterification reactions (GO:0000375)
#******************************************************************************#

#==============================================================================#
# Heatmap plot
#==============================================================================#

# Get selected and non-selected genes
ID <- GOannotation$ID[GOannotation$Name %in% "RNA splicing, via transesterification reactions (GO:0000375)"]
sel_genes <- GOgenes[[ID]]
egs_sel <- egs[egs$ensembl_gene_id %in% sel_genes,]
egs_nonsel <- egs[!(egs$ensembl_gene_id %in% sel_genes),]

# Settings
range <- 3000
lo <- 60

# Get counts around TSS of selected genes
tiles_sel = sapply(1:nrow(egs_sel), function(i)
  if(egs_sel$strand[i] == "1" )
    egs_sel$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_sel$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_sel = GRanges(tilename = paste(rep(egs_sel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                    seqnames = Rle(rep(paste0('chr', egs_sel$chromosome_name), each=lo)), 
                    ranges = IRanges(start = as.vector(tiles_sel),
                                     width = 100),
                    strand = Rle(rep("*", length(as.vector(tiles_sel)))),
                    seqinfo=si)

# Get bin counts
selLoci = countOverlaps(tiles_sel, MeCP2_rep1_fil) + 
  countOverlaps(tiles_sel, MeCP2_rep2_fil)

selLoci_matrix = matrix(selLoci, nrow=nrow(egs_sel), 
                        ncol=lo, byrow=TRUE)

rownames(selLoci_matrix) <- egs_sel$ensembl_gene_id
colnames(selLoci_matrix) <- seq(-1*range, range, length.out=lo)

# Make data frame for plotting
plotMatrix <- gather(as.data.frame(selLoci_matrix))
plotMatrix$Gene <- rep(rownames(selLoci_matrix), ncol(selLoci_matrix))
plotMatrix$key <- as.numeric(plotMatrix$key)
plotMatrix$Gene <- factor(plotMatrix$Gene, levels = names(sort(rowSums(selLoci_matrix))))

# Make heatmap plot
p <- ggplot(plotMatrix) +
  geom_tile(aes(x = key, y = Gene, fill = value)) +
  scale_fill_viridis_c(option = "plasma", limits = c(0,70),
                       na.value = "yellow") +
  xlab("Distance from TSS") +
  ggtitle("RNA splicing, via transesterification reactions") +
  labs(fill = "Enrichment") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# Save plot
ggsave(p, file = "Heatmap_RNA splicing.png",
       width = 6, height = 6)


#==============================================================================#
# Density plot
#==============================================================================#

# Get counts around TSS of non-selected genes
tiles_nonsel = sapply(1:nrow(egs_nonsel), function(i)
  if(egs_nonsel$strand[i] == "1" )
    egs_nonsel$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_nonsel$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_nonsel = GRanges(tilename = paste(rep(egs_nonsel$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                       seqnames = Rle(rep(paste0('chr', egs_nonsel$chromosome_name), each=lo)), 
                       ranges = IRanges(start = as.vector(tiles_nonsel),
                                        width = 100),
                       strand = Rle(rep("*", length(as.vector(tiles_nonsel)))),
                       seqinfo=si)

# Get bin counts
nonselLoci = countOverlaps(tiles_nonsel, MeCP2_rep1_fil) + 
  countOverlaps(tiles_nonsel, MeCP2_rep2_fil)

nonselLoci_matrix = matrix(nonselLoci, nrow=nrow(egs_nonsel), 
                           ncol=lo, byrow=TRUE)

# prepare data for plotting
plotDF <- data.frame(value = c(colMeans(selLoci_matrix),
                               colMeans(nonselLoci_matrix)),
                     location = rep(seq(-1*range, range, length.out=lo),2),
                     group = c(rep("Selected",lo),
                               rep("Non-selected",lo))
)

# Get mean counts around TSS for different permutations
nPerm <- 100
plotDF_permutation  <- NULL
set.seed(123)
for (p in 1:nPerm){
  egs_permutation <- egs[sample(1:nrow(egs),nrow(selLoci_matrix)),]
  
  # Get tiles around TSS
  range <- 3000
  lo <- 60
  tiles_permutation = sapply(1:nrow(egs_permutation), function(i)
    if(egs_permutation$strand[i] == "1" )
      egs_permutation$TSS[i] + seq(-1*range, range-100, length.out=lo)
    else
      egs_permutation$TSS[i] + seq(range-100, -1*range, length.out=lo))
  
  tiles_permutation= GRanges(tilename = paste(rep(egs_permutation$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                             seqnames = Rle(rep(paste0('chr', egs_permutation$chromosome_name), each=lo)), 
                             ranges = IRanges(start = as.vector(tiles_permutation),
                                              width = 100),
                             strand = Rle(rep("*", length(as.vector(tiles_permutation)))),
                             seqinfo=si)
  
  # Get bin counts
  permutationLoci = countOverlaps(tiles_permutation, MeCP2_rep1_fil) + 
    countOverlaps(tiles_permutation, MeCP2_rep2_fil)
  
  permutationLoci_matrix = matrix(permutationLoci, nrow=nrow(egs_permutation), 
                                  ncol=lo, byrow=TRUE)
  
  # Put results into data frame
  temp <- data.frame(value = colMeans(permutationLoci_matrix),
                     location = rep(seq(-1*range, range, length.out=lo),2),
                     group = p)
  
  # Combine with the results from previous permutation
  plotDF_permutation <- rbind.data.frame(plotDF_permutation,temp)
}

# Save plotting data frames
plotDF$group <- ifelse(plotDF$group == "Selected", "RNA splicing, via transesterification reactions", "Other")
save(plotDF, plotDF_permutation, 
     file = "Data/ChIPSeqData_RNA splicing.RData")

# Load plotting data frames (if needed)
load("Data/ChIPSeqData_RNA splicing.RData")

# Make plot
colors <- rev(c("#EF3B2C", "#737373"))
p <- ggplot() +
  geom_line(data = plotDF_permutation, 
            aes(x = location, y = value, group = group),
            linewidth = 0.5, color = "#D9D9D9") +
  geom_line(data = plotDF,
            aes(x = location, y = value, group = group, color = group),
            linewidth = 1.5) +
  ylab("Mean enrichment") +
  xlab("Distance from the TSS") +
  ylim((c(10,40)))+
  labs(color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_color_manual(values = colors)

# Save plot
ggsave(p, file = "Density_RNA splicing.png",
       width = 6, height = 4)
