# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(readxl)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/8. Imprinting"))

# Load data
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_norm.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/DEresults_RTTvsIC_gx.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/geneAnnotation.RData"))
load(paste0(homeDir,"/sampleInfo.RData"))
genes_all <- rownames(topList[[1]])

# Get p-values
pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues) <- genes_all
colnames(pvalues) <- names(topList)

# get adjusted p-values
adj_pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"FDR"]})), ncol = length(topList))
rownames(adj_pvalues) <- genes_all
colnames(adj_pvalues) <- names(topList)

# Get logFCs
logFCs <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
rownames(logFCs) <- genes_all
colnames(logFCs) <- names(topList)

# Collect logFCs of the MEG members
logFCs_MEGs <- logFCs[c("ENSG00000214548_MEG3",
                        "ENSG00000225746_MEG8",
                        "ENSG00000223403_MEG9"),]

# Put data in data.frame
plotDF <- gather(as.data.frame(logFCs_MEGs))
plotDF$GeneID <- rep(rownames(logFCs_MEGs), ncol(logFCs_MEGs))

# Add Gene name
plotDF$GeneName <- unlist(lapply(str_split(plotDF$GeneID, "_"), function(x){x[2]}))

# Add region
plotDF$Region <- "iPSC"
plotDF$Region[str_detect(plotDF$key, "Ventral")] <- "Ventral"
plotDF$Region[str_detect(plotDF$key, "Dorsal")] <- "Dorsal"

# Add time
plotDF$Time <- "D0"
plotDF$Time[str_detect(plotDF$key, "D13")] <- "D13"
plotDF$Time[str_detect(plotDF$key, "D40")] <- "D40"
plotDF$Time[str_detect(plotDF$key, "D75")] <- "D75"

# Order the time points
plotDF$key <- factor(plotDF$key,
                     levels = c("Dorsal_D75_RTTvsIC",  
                                "Dorsal_D40_RTTvsIC",
                                "Dorsal_D13_RTTvsIC",
                                "Cell_D0_RTTvsIC",  
                                "Ventral_D13_RTTvsIC",
                                "Ventral_D40_RTTvsIC",
                                "Ventral_D75_RTTvsIC"))

# Make main plot
mainPlot <- ggplot(plotDF) +
  geom_line(aes(x = key, y = value, color = GeneName, group = GeneName),
            linewidth = 1.5) +
  #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5), fill= "grey", alpha = 0.05) +
  geom_hline(yintercept = 1, color = "black", linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "black", linewidth = 1, linetype = "solid") +
  geom_text(x = "Cell_D0_RTTvsIC", y = 1.4, label = expression(log[2]~"FC = 1")) +
  ylab(expression(log[2]~"FC")) +
  xlab(NULL) +
  ylim(c(0,10)) +
  theme_bw() +
  scale_color_brewer(palette = "Accent") +
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotDF$key,
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

ggsave(finalPlot, file = paste0("MEG_logFCs.png"), width = 6, height = 4)
