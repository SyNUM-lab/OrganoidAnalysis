
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(igraph)
library(RCy3)

# FUNCTION: Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# FUNCTION: calculate Jaccard Index
JI <- function(x,y){length(intersect(x,y))/length(union(x,y))}

# Set working directory
setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/5. GSEA")

# Load necessary data
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData")
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOannotation.RData")
load("Data/terms_ordered1.RData")
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/SampleInfo.RData")

# Calculate pairwise Jaccard Index of selected GO terms
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
graph_matrix <- matrix(NA, nrow = length(GOgenes_fil), ncol = length(GOgenes_fil))
for (i in 1:length(GOgenes_fil)){
  graph_matrix[i,] <- unlist(lapply(GOgenes_fil,function(x){JI(x,GOgenes_fil[[i]])}))
}
rownames(graph_matrix) <- names(GOgenes_fil)
colnames(graph_matrix) <- names(GOgenes_fil)

# Set names of matrix
rownames(GOannotation) <- GOannotation$ID
rownames(graph_matrix) <- firstup(GOannotation[rownames(graph_matrix),"Description"])
colnames(graph_matrix) <- firstup(GOannotation[colnames(graph_matrix),"Description"])

# Make graph
g <- graph_from_adjacency_matrix(graph_matrix, 
                                 mode = "lower", 
                                 weighted = "weight", 
                                 diag = FALSE)
# Export graph into Cytoscape
createNetworkFromIgraph(
  g,
  title = "Jaccard Index",
  collection = "GO similarity"
)




