
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readxl)
library(biomaRt)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/8. Imprinting/ASE"))

# Select SNP IDs
pos <- c("chr1:228321726", "chr14:100734523", "chr19:31276811", "chr19:58571238",  
         "chr5:95755497", "chr7:78018005")

# Get annotation data for the selected SNPs
granges <- TxDb.Hsapiens.UCSC.hg38.knownGene
granges <- transcriptsBy(granges, by = "gene")

gr_pos <- makeGRangesFromDataFrame(
  data.frame(
    chr=c("chr1", "chr14", "chr19", "chr19","chr5", "chr7"),
    start=c(228321726, 100734523, 31276811, 58571238,  
            95755497, 78018005)-500000,
    end=c(228321726, 100734523, 31276811, 58571238,  
          95755497, 78018005)+500000),
  keep.extra.columns=TRUE)

names(gr_pos) <- pos
overlaps <- findOverlaps(granges, gr_pos, type = 'any')

idx_query <- queryHits(overlaps)
genes1 <- names(granges)[idx_query] 

idx_subject <- subjectHits(overlaps)
pos_info1 <- names(gr_pos)[idx_subject]

genomeRegion <- data.frame(Gene = genes1,
                           Pos = pos_info1)


#==============================================================================#
# Plot signed -log P value of genes in the region
#==============================================================================#

# Load expression data
load(paste0(homeDir,"/sampleInfo.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_norm.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/DEresults_RTTvsIC_gx.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/geneAnnotation.RData"))
geneAnnotation$entrezgene_id <- as.character(geneAnnotation$entrezgene_id)
geneAnnotation <- geneAnnotation[geneAnnotation$chromosome_name != "X",]
geneAnnotation <- unique(geneAnnotation[,c(1,2,3,7)])

# MEG8 did not have an associated ENTREZ gene id, so I add it manually
geneAnnotation$entrezgene_id[geneAnnotation$GeneName == "MEG8"] <- "79104"

# Get logFCs and P values
genes_all <- rownames(topList[[1]])
logFCs <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))

spvalue <- -log10(pvalues)*sign(logFCs)

rownames(spvalue) <- genes_all
colnames(spvalue) <- names(topList)
spvalue_df <- data.frame(spvalue,
                       ID = rownames(spvalue))

genomeRegion <- inner_join(genomeRegion, geneAnnotation, by = c("Gene" = "entrezgene_id"))
genomeRegion <- inner_join(genomeRegion, spvalue_df, by = c("gene_id" = "ID"))

plotDF <- NULL
for (i in 6:12){
  temp <- genomeRegion[,c(1:5,i)]
  temp$TimeRegion <- colnames(temp)[6]
  colnames(temp) <- c(colnames(temp)[1:5],"spvalue", "TimeRegion")
  plotDF <- rbind.data.frame(plotDF,temp)
}


plotDF$Time <- unlist(lapply(str_split(plotDF$TimeRegion, "_"), function(x){x[2]}))
plotDF$Region <- unlist(lapply(str_split(plotDF$TimeRegion, "_"), function(x){x[1]}))
plotDF$Region[plotDF$Region == "Cell"] <- "iPSC"


plotDF$Region <- factor(plotDF$Region, levels = c("iPSC", "Dorsal", "Ventral"))

colors <- RColorBrewer::brewer.pal(n = 5, name = "Reds")[2:5]

plotDF$Pos <- factor(plotDF$Pos, levels = c("chr1:228321726",
                                            "chr19:31276811",
                                            "chr19:58571238",
                                            "chr5:95755497",   
                                            "chr7:78018005",
                                            "chr14:100734523"))

# Get start location of genes: this will be used to sort the genes in the plot
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
annotations <- getBM(attributes=c("start_position", "entrezgene_id"), 
                     filters = "entrezgene_id",
                     values = unique(plotDF$Gene),
                     mart = ensembl)
annotations$entrezgene_id <- as.character(annotations$entrezgene_id)

# Add annotations to the data frame
plotDF <- left_join(plotDF,annotations, by = c("Gene" = "entrezgene_id"))
plotDF$start_position[plotDF$GeneName == "MEG8"] <- 100889649

# Sort based on the start location
plotDF <- arrange(plotDF, by = start_position)
plotDF$GeneName <- factor(plotDF$GeneName, levels = unique(plotDF$GeneName))
plotDF <- unique(plotDF[,1:10])

# Panel names should be empty
empty_names <- setNames(rep("",6), levels(plotDF$Pos))

# Make plot
p <- ggplot(plotDF) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x = GeneName,y = spvalue, color = Time, shape = Region)) +
  facet_grid(cols = vars(Pos), scale = "free", space = "free", labeller = labeller(Pos = empty_names)) +
  xlab(NULL) +
  ylab(expression("Signed"~-log[10]~"P value")) +
  guides(shape = guide_legend(override.aes = list(size = 4)),
         color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(angle = 90, size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")) +
  scale_color_manual(values = colors)

panel_colors <- c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC")

# convert to grob
gp <- ggplotGrob(p)
for(i in 1:6){
  grob.i <- grep("strip-t", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}
#gp <- grid::grid.draw(gp)

# Save plot
ggsave(gp, file = paste0("RegionASE.png"), width = 9, height = 4)
