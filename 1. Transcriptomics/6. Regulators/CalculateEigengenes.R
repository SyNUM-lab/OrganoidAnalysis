# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(igraph)
library(tidyverse)
library(ggpubr)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/6. Regulators"))

# Load data:

# GO terms
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

# GO annotation
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData"))
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))

# Gene annotation
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"geneAnnotation.RData"))

# Expression matrix
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))

# Sample info
load(paste0(homeDir,"/SampleInfo.RData"))
sampleInfo$Tissue[sampleInfo$Tissue == "Cell"] <- "iPSC"


# Prepare empty matrix for eigengenes
eigengenes <- matrix(NA, nrow = length(terms_ordered), ncol = ncol(gxMatrix_norm))
rownames(eigengenes) <- terms_ordered
colnames(eigengenes) <- colnames(gxMatrix_norm)

# Calculate eigengenes for each process
for (i in 1:length(terms_ordered)){
  # Get GO id
  GOid <- GOannotation$ID[GOannotation$Name == terms_ordered[i]]
  
  # Retrieve genes
  geneList <- GOgenes[[GOid]]
  
  # filter gx matrix
  gxMatrix_fil <- gxMatrix_norm[rownames(gxMatrix_norm) %in% geneAnnotation$gene_id[geneAnnotation$EnsemblID %in% geneList],]
  
  # Run PCA
  pca <- prcomp(t(gxMatrix_fil), 
                retx = TRUE, # Give rotated variables (PCA loadings)
                center = TRUE,
                scale = TRUE) # Variables are scaled to unit variance
  
  # Get eigengene
  eigengenes[i,] <- pca$x[,1]
  
}

# Save data
save(eigengenes, file = "Data/eigengenes_all.RData")

