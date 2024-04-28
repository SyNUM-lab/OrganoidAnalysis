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

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/4. DEAnalysis")

# Load data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/SampleInfo.RData")

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

# Set colors
color_time <- c("#FEE5D9","#FCAE91","#FB6A4A","#CB181D","#FCAE91","#FB6A4A","#CB181D")

# Make volcano plot for each of the seven comparisons
for (i in 1:ncol(pvalues)){
  
  # Prepare data for plotting
  plotDF <- data.frame(Pvalue = pvalues[,i],
                       adjPvalue = adj_pvalues[,i],
                       logFC = logFCs[,i],
                       GeneID = rownames(pvalues))
  
  # make plot
  p <- ggplot() +
    geom_point(data = plotDF, 
               aes(x = logFC, y = -log10(Pvalue)),
               color = "#F0F0F0") +
    geom_point(data = plotDF[(plotDF$adjPvalue < 0.05) & (abs(plotDF$logFC) > 1),], 
               aes(x = logFC, y = -log10(Pvalue)),
               color = color_time[i]) +
    xlab(expression(log[2]~"FC")) +
    ylab(expression(-log[10]~"P value")) +
    labs(shape = NULL) +
    xlim(c(-17,17)) +
    ylim(c(0,50)) +
    scale_shape_manual(values = c(21,22,24)) +
    theme_void() +
    theme(axis.line.x = element_line(),
          axis.ticks.x = element_line(linewidth = 0.5, lineend = "butt"),
          axis.text.x = element_text(),
          axis.title.x = element_text(size = 12),
          axis.line.y = element_line(),
          axis.ticks.y = element_line(linewidth = 0.5, lineend = "butt"),
          axis.text.y = element_text(hjust = 1, vjust = 0.5,
                                     margin = margin(t = 0, r = 2, b = 0, l = 0)),
          axis.title.y = element_text(angle = 90, vjust = 1,
                                      margin = margin(t = 0, r = 5, b = 0, l = 0),
                                      size = 12),
          axis.ticks.length  = unit(0.1, "cm"),
          plot.margin = margin(1,1,1,1, "cm"),
          legend.position = "bottom")
  
  # Save plot
  ggsave(p, 
         file = paste0("Plots/",
                       colnames(pvalues)[i],
                       "_volcano.png"),
         width = 6, height = 4.5)
}
