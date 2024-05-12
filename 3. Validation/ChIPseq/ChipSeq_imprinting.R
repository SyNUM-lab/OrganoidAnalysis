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

# Get sequence info
genome = BSgenome.Mmusculus.UCSC.mm9
si = seqinfo(genome)
si = si[ paste0('chr', c(1:19))]

# Import BW files (GSE71126, GEO)
MeCP2_rep1 <- import.bw("Data/GSE71126_RAW/GSM1827604_MeCP2_ChIP_WT_rep1.extend200.coverage.mm9.bw")
MeCP2_rep2 <- import.bw("Data/GSE71126_RAW/GSM1827605_MeCP2_ChIP_WT_rep2.extend200.coverage.mm9.bw")
Input_rep1 <- import.bw("Data/GSE71126_RAW/GSM1827606_Input_WT_rep1.extend200.coverage.mm9.bw")
Input_rep2 <- import.bw("Data/GSE71126_RAW/GSM1827607_Input_WT_rep2.extend200.coverage.mm9.bw")

# Select chr1 - chr19 (tremove sex and mitochrondrial genes)
MeCP2_rep1_fil <- MeCP2_rep1[!(seqnames(MeCP2_rep1) %in% c("chrM", "chrY", "chrX"))]
rm(MeCP2_rep1)
MeCP2_rep2_fil <- MeCP2_rep2[!(seqnames(MeCP2_rep2) %in% c("chrM", "chrY", "chrX"))]
rm(MeCP2_rep2)
Input_rep1_fil <- Input_rep1[!(seqnames(Input_rep1) %in% c("chrM", "chrY", "chrX"))]
rm(Input_rep1)
Input_rep2_fil <- Input_rep2[!(seqnames(Input_rep2) %in% c("chrM", "chrY", "chrX"))]
rm(Input_rep2)


#******************************************************************************#
# Get imprinted genes
#******************************************************************************#

mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl',
               host="may2009.archive.ensembl.org")

# Load genome annotations for each chromosome
egs = getBM(attributes = c('ensembl_gene_id','external_gene_id',
                           'chromosome_name','start_position',
                           'end_position','strand'), 
            filters='chromosome_name',
            values=1:19,
            mart=mart)

# Get transcription start site
egs$TSS = ifelse(egs$strand == "1", egs$start_position, egs$end_position)

# Save genome annotations
save(egs, file = "Data/egs.RData")

# load data (if needed)
load("Data/egs.RData")

# All genes in data
all_genes <- egs$external_gene_id

# Imprinted genes in data (https://www.geneimprint.com/site/genes-by-species.Mus+musculus)
imprintDF <- read_xlsx("Data/imprinted_genes_mouse.xlsx")
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,4)],]
imprinted_genes <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[1])]
imprinted_genes <- imprinted_genes[!is.na(imprinted_genes)]
imprinted_genes <- intersect(imprinted_genes, all_genes)

# Imprinted vs non-imprinted genes
egs_imprinted <- egs[egs$external_gene_id %in% imprinted_genes,]
egs_nonimprinted <- egs[!(egs$external_gene_id %in% imprinted_genes),]

# Settings:
range <- 3000 # Range around TSS
lo <- 60      # length out (bin size)

#==============================================================================#
# 1) Get bin counts around TSS of imprinted genes
#==============================================================================#

# Set tiles around TSS of imprinted genes
tiles_imprinted = sapply(1:nrow(egs_imprinted), function(i)
  if(egs_imprinted$strand[i] == "1" )
    egs_imprinted$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_imprinted$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_imprinted = GRanges(tilename = paste(rep(egs_imprinted$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                seqnames = Rle(rep(paste0('chr', egs_imprinted$chromosome_name), each=lo)), 
                ranges = IRanges(start = as.vector(tiles_imprinted),
                                 width = 100),
                strand = Rle(rep("*", length(as.vector(tiles_imprinted)))),
                seqinfo=si)

# Get bin counts
imprintedLoci = countOverlaps(tiles_imprinted, MeCP2_rep1_fil) + 
  countOverlaps(tiles_imprinted, MeCP2_rep2_fil)

imprintedLoci_matrix = matrix(imprintedLoci, nrow=nrow(egs_imprinted), 
                           ncol=lo, byrow=TRUE)


#==============================================================================#
# 2) Get bin counts around TSS of non-imprinted genes
#==============================================================================#

# Set tiles around the TSS of non-imprinted genes
tiles_nonimprinted = sapply(1:nrow(egs_nonimprinted), function(i)
  if(egs_nonimprinted$strand[i] == "1" )
    egs_nonimprinted$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_nonimprinted$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_nonimprinted = GRanges(tilename = paste(rep(egs_nonimprinted$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                          seqnames = Rle(rep(paste0('chr', egs_nonimprinted$chromosome_name), each=lo)), 
                          ranges = IRanges(start = as.vector(tiles_nonimprinted),
                                           width = 100),
                          strand = Rle(rep("*", length(as.vector(tiles_nonimprinted)))),
                          seqinfo=si)

# Get bin counts
nonimprintedLoci = countOverlaps(tiles_nonimprinted, MeCP2_rep1_fil) + 
  countOverlaps(tiles_nonimprinted, MeCP2_rep2_fil)

nonimprintedLoci_matrix = matrix(nonimprintedLoci, nrow=nrow(egs_nonimprinted), 
                              ncol=lo, byrow=TRUE)

# prepare data for plotting
plotDF <- data.frame(value = c(colMeans(imprintedLoci_matrix),
                               colMeans(nonimprintedLoci_matrix)),
                     location = rep(seq(-1*range, range, length.out=lo),2),
                     group = c(rep("Imprinted (n = 79)",lo),
                               rep("Non-imprinted (n = 25,900)",lo))
)


#==============================================================================#
# 3) Get bin counts around TSS for random permutations of genes
#==============================================================================#

nPerm <- 100 # Number of permutations
plotDF_permutation  <- NULL
set.seed(123)
for (p in 1:nPerm){
  
  # Select 79 random genes (same number of imprinted genes)
  egs_permutation <- egs[sample(1:nrow(egs),79),]
  
  # Set tiles around the TSS of these 79 genes
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
  
  
  # Put result into data frame
  temp <- data.frame(value = colMeans(permutationLoci_matrix),
                     location = rep(seq(-1*range, range, length.out=lo),2),
                     group = p)
  
  # Combine with results from previous permutation
  plotDF_permutation <- rbind.data.frame(plotDF_permutation,temp)
}


# Save plotting data
save(plotDF, plotDF_permutation, 
     file = "Data/plotDF_ChIPSeq.RData")

#==============================================================================#
# 4) Make plot
#==============================================================================#

# Load data for plotting
load("Data/plotDF_ChIPSeq.RData")

# Make plot
colors <- c("#EF4040", "#FFA732")
p <- ggplot() +
  geom_line(data = plotDF_permutation, 
            aes(x = location, y = value, group = group),
            linewidth = 0.5, color = "#FFE4B5") +
  geom_line(data = plotDF,
            aes(x = location, y = value, group = group, color = group),
            linewidth = 1.5) +
  ylab("Mean tag count") +
  xlab("Distance from the TSS") +
  labs(color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_color_manual(values = colors)

# Save plot
ggsave(p, 
       file = "ChipSeqPlot.png",
       width = 6, height = 4.5)
