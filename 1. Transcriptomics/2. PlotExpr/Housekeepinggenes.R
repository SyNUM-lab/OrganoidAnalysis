# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/2. PlotExpr/"))

# Load data
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(preprocessing_dir,"FPKM.RData"))
load(paste0(homeDir,"/SampleInfo.RData"))

geneAnnotation1 <- unique(geneAnnotation[,c(1,6)])
rownames(geneAnnotation1) <- geneAnnotation1$gene_id

#==============================================================================#
# select genes
test <- geneAnnotation[which(str_detect(geneAnnotation$GeneName, "GAPDH")),]
selectedGene <- c("ENSG00000142330_CAPN10", 
                  "ENSG00000176396_EID2",
                  "ENSG00000150991_UBC",
                  "ENSG00000112592_TBP")
#==============================================================================#

# Get HGNC symbol of selected gene
#nameGene <- geneAnnotation$GeneName[geneAnnotation$gene_id == selectedGene][1]

# Select expression data of selected protein
dataMatrix_sel <-  gather(data.frame(FPKM[rownames(FPKM) %in% selectedGene,]))
dataMatrix_sel$gene_id <- rep(rownames(FPKM[rownames(FPKM) %in% selectedGene,]), ncol(FPKM))

# Prepare data for plotting
plotDF <- inner_join(sampleInfo, dataMatrix_sel, by = c("SampleID" = "key"))
plotDF$Tissue[plotDF$Tissue == "Cell"] <- "iPSC"
plotDF$TimeRegion <- paste0(plotDF$Tissue, "_",plotDF$Time)
plotDF$TimeRegion <- factor(plotDF$TimeRegion,levels = c("Dorsal_D75",
                                                         "Dorsal_D40",
                                                         "Dorsal_D13",
                                                         "iPSC_D0",
                                                         "Ventral_D13",
                                                         "Ventral_D40",
                                                         "Ventral_D75"))
plotDF$ID <-  paste0(plotDF$TimeRegion, "_", plotDF$Gene)

plotDF <- inner_join(plotDF, geneAnnotation, by = c("gene_id" = "gene_id"))

# Main plot: expression over time
mainPlot <- ggplot(plotDF) +
  geom_point(aes(x = TimeRegion, y = value, color = GeneName, shape = Group), size = 2) +
  geom_line(aes(x = TimeRegion, y = value, color = GeneName, group = paste0(Replicate, Group, gene_id))) +
  ylab(expression(log[2]~"FPKM")) +
  xlab(NULL) +
  #ggtitle(paste0(nameGene, " (", selectedGene, ")")) +
  #ggtitle(nameGene) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) +
  scale_color_brewer(palette = "Dark2")

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotDF$TimeRegion,
  time = plotDF$Time,
  region = plotDF$Tissue))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# Tissue/region
colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = region)) +
  geom_text(data = colSideColor[(str_detect(colSideColor$sample, "D40") |
                                   str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = region)) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots into a single figure
finalPlot <- ggarrange(mainPlot,
                       colSideColorPlot_time,
                       colSideColorPlot_tissue,
                       heights = c(6,0.5,0.5),nrow = 3,ncol = 1,
                       common.legend = TRUE,
                       legend = "right",
                       align = "v")
# Print plot
finalPlot

# Save plot
ggsave(finalPlot, file = "Housekeepinggenes/Housekeepinggenes.png", width = 7, height = 4)

