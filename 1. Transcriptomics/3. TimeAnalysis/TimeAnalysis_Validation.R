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

# Set working directory
homeDir <- "E:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/3. TimeAnalysis"))

# Data from https://stemcell.libd.org/scb/
load("libd_stemcell_timecourse_rseGene_n157.rda")
metaData <- as.data.frame(rse_gene@colData)
exprData <- getRPKM(rse_gene)
rownames(exprData) <- str_remove(rownames(exprData), "\\.[^.]*")

# Own data data:
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(homeDir,"/SampleInfo.RData"))
all(sampleInfo$SampleID == colnames(gxMatrix_norm))

dorsal <- sampleInfo[sampleInfo$Tissue == "Cell",]
dorsal$Tissue <- "Dorsal"
ventral <- sampleInfo[sampleInfo$Tissue == "Cell",]
ventral$Tissue <- "Ventral"
sampleInfo <- rbind.data.frame(sampleInfo[sampleInfo$Tissue != "Cell",],
                               rbind.data.frame(dorsal,ventral))


# GO terms 
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData"))
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
GOterms <- c("GO:0140053",
             "GO:0009451",
             "GO:0007059",
             "GO:0006360",
             "GO:0000959",
             "GO:0050890",
             "GO:0050803",
             "GO:0008038",
             "GO:0099601",
             "GO:0035249")

GOnames <- c("Mitochondrial gene expression",
             "RNA modification",
             "Chromosome segregation",
             "Transcription by RNA polymerase I",
             "Mitochondrial RNA metabolic process",
             "Cognition",
             "Regulation of synapse structure or activity",
             "Neuron recognition",
             "Regulation of neurotransmitter receptor activity",
             "Synaptic transmission, glutamatergic"
             )
colors <- c("#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#4A1486",
            "#FDAE6B","#FD8D3C","#F16913","#D94801","#8C2D04")
names(colors) <- c("Mitochondrial gene expression",
                   "Transcription by RNA polymerase I",
                   "Mitochondrial RNA metabolic process",
                   "RNA modification",
                   "Chromosome segregation",
                   "Neuron recognition",
                   "Regulation of synapse structure or activity",
                   "Synaptic transmission, glutamatergic",
                   "Cognition",
                   "Regulation of neurotransmitter receptor activity"
)

for (g in 1:length(GOterms)){
  
  # Select genes
  genes <- GOgenes[[GOterms[g]]]
  
  
  #============================================================================#
  # LIBD data 
  #============================================================================#
  
  # Filter expression data
  exprData_fil <- exprData[rownames(exprData) %in% genes,]
  exprData_fil <- exprData_fil[!(apply(exprData_fil,1,var)==0),]
  
  # Scale the data
  exprData_scaled <- (exprData_fil - rowMeans(exprData_fil))/apply(exprData_fil,1,sd)
  
  # Get the median expression per sample
  plotDF <- data.frame(value = colMedians(exprData_scaled),
                       sampleID = colnames(exprData_scaled))
  
  # Add the meta data
  plotDF <- inner_join(plotDF, metaData, by = c("sampleID" = "SampleID"))
  plotDF <- plotDF[plotDF$DX == "CNT",]
  plotDF$Group <- paste0(plotDF$LINE, "_", plotDF$BIO_REP)
  plotDF$GroupDay <- paste0(plotDF$Group, "_", plotDF$DAY)
  plotDF <- plotDF %>%
    group_by(GroupDay) %>%
    reframe(meanValue = mean(value),
              DAY = DAY,
              Group = Group)
  plotDF <- plotDF[!duplicated(plotDF$GroupDay),]
  
  # Make plot
  p1 <- ggplot(plotDF) +
    geom_line(aes(x = as.numeric(DAY), y = meanValue, group = Group, color = GOnames[g])) +
    geom_point(aes(x = as.numeric(DAY), y = meanValue, color = GOnames[g])) +
    ggtitle("LIBD") +
    ylab("Median scaled expression") +
    xlab("Days in vitro") +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      legend.position = "none"
    )
  
  #============================================================================#
  # Our data 
  #============================================================================#
  
  # Filter expression data
  exprData_fil <- gxMatrix_norm[rownames(gxMatrix_norm) %in% geneAnnotation$gene_id[geneAnnotation$EnsemblID %in% genes],]
  exprData_fil <- exprData_fil[!(apply(exprData_fil,1,var)==0),]
  
  # Scale the data
  exprData_scaled <- (exprData_fil - rowMeans(exprData_fil))/apply(exprData_fil,1,sd)
  
  # Get the median expression per sample
  plotDF <- data.frame(value = colMedians(exprData_scaled),
                       sampleID = colnames(exprData_scaled))
  
  plotDF <- inner_join(plotDF, sampleInfo, by = c("sampleID" = "SampleID"))
  plotDF$Linegroup <- paste0(plotDF$Group, "_", plotDF$Tissue, "_", plotDF$Replicate)
  plotDF$Time1 <- str_remove(plotDF$Time, "D")
  
  # Make plot
  p2 <- ggplot(plotDF) +
    geom_line(aes(x = as.numeric(Time1), y = value, group = Linegroup, color = GOnames[g])) +
    geom_point(aes(x = as.numeric(Time1), y = value, color = GOnames[g])) +
    ggtitle("Our study") +
    ylab(NULL) +
    xlab("Days in vitro") +
    #scale_color_manual(values = c("#1B9E77", "#7570B3")) +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      legend.position = "none"
    )
  
  p_all <- p1 + p2 + patchwork::plot_layout(nrow = 1, ncol = 2)
  
  ggsave(p_all, file = paste0("Validation/", GOnames[g], ".png"),
         width = 7, height = 4)
}









