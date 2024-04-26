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

#==============================================================================#
# Gene expression
#==============================================================================#

# Load data:

# GO terms
load("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# GO annotation
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOannotation.RData"))
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOgenes_ENSEMBL.RData"))

# Gene annotation
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/geneAnnotation.RData")

# Expression matrix
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_norm.RData")

# Sample info
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")
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

save(eigengenes, file = "D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Eigengene/eigengenes_all.RData")

