# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readxl)
library(ggpubr)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/8. Imprinting/ASE"))

# Load data
load("outputASE.RData")
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/geneAnnotation.RData"))

# Filter for imprinted genes
ref_imp <- ref[rownames(ref) %in% pos_info_imprinted$id,]
alt_imp <- alt[rownames(alt) %in% pos_info_imprinted$id,]
total_imp <- total[rownames(total) %in% pos_info_imprinted$id,]
score_imp <- score[rownames(score) %in% pos_info_imprinted$id,]

# Make sample information data frame
sampleInfo <- data.frame(
  SampleID = colnames(alt),
  Region = "iPSC",
  Time = "D0",
  Group = "IC",
  Replicate = "1"
)

sampleInfo$Region[str_detect(sampleInfo$SampleID, "Ventral")] <- "Ventral"
sampleInfo$Region[str_detect(sampleInfo$SampleID, "Dorsal")] <- "Dorsal"

sampleInfo$Time[str_detect(sampleInfo$SampleID, "D13")] <- "D13"
sampleInfo$Time[str_detect(sampleInfo$SampleID, "D40")] <- "D40"
sampleInfo$Time[str_detect(sampleInfo$SampleID, "D75")] <- "D75"

sampleInfo$Group[str_detect(sampleInfo$SampleID, "RTT")] <- "RTT"

sampleInfo$TimeRegion <- paste0(sampleInfo$Time, "_", sampleInfo$Region)
sampleInfo$TimeRegionGroup <- paste0(sampleInfo$Time, "_", 
                                     sampleInfo$Region, "_",
                                     sampleInfo$Group)

sampleInfo$Replicate[str_detect(sampleInfo$SampleID, "_2")] <- "2"
sampleInfo$Replicate[str_detect(sampleInfo$SampleID, "_3")] <- "3"

# Select SNP ID for plotting
selectedID <- "chr1:228321726"
geneName <- pos_info_imprinted$hgnc_symbol[pos_info_imprinted$id == selectedID]

# Select expression data of selected protein
dataMatrix_sel <- data.frame(Expr = score_imp[selectedID,],
                             SampleID = colnames(score_imp))

# Prepare data for plotting
plotDF <- inner_join(sampleInfo, dataMatrix_sel, by = c("SampleID" = "SampleID"))
plotDF$TimeRegion <- factor(plotDF$TimeRegion,levels = c("D75_Dorsal",
                                                         "D40_Dorsal",
                                                         "D13_Dorsal",
                                                         "D0_iPSC",
                                                         "D13_Ventral",
                                                         "D40_Ventral",
                                                         "D75_Ventral"))



# Main plot: expression over time
mainPlot <- ggplot(plotDF) +
  geom_point(aes(x = TimeRegion, y = Expr, color = Group), size = 2) +
  geom_line(aes(x = TimeRegion, y = Expr, color = Group, group = paste0(Replicate, Group))) +
  ylab("ASE") +
  xlab(NULL) +
  ggtitle(geneName) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) +
  scale_color_manual(values = c("#377EB8","#E41A1C"))

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotDF$TimeRegion,
  time = plotDF$Time,
  region = plotDF$Region))

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
