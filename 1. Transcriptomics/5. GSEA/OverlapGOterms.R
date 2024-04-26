# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(igraph)
library(RCy3)

# Load data
load("Data/reducedTerms_RTTvsIC_BP1.RData")
load("Data/simMatrix_RTTvsIC_BP1.RData")
load("Data/terms_ordered1.RData")
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOannotation.RData"))

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Get GO similarity of terms
terms <- GOannotation$ID[GOannotation$Name %in% terms_ordered]
graph_matrix <- simMatrix[terms,terms]

# Set names of matrix
rownames(GOannotation) <- GOannotation$ID
rownames(graph_matrix) <- firstup(GOannotation[rownames(graph_matrix),"Description"])
colnames(graph_matrix) <- firstup(GOannotation[colnames(graph_matrix),"Description"])




g <- graph_from_adjacency_matrix(graph_matrix, 
                                 mode = "lower", 
                                 weighted = "weight", 
                                 diag = FALSE)

createNetworkFromIgraph(
  g,
  title = "GO relations1",
  collection = "test"
)





# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(igraph)
library(RCy3)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

JI <- function(x,y){length(intersect(x,y))/length(union(x,y))}

# Load GO annotation data
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOgenes_ENSEMBL.RData"))
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOannotation.RData"))
load("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_norm.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/DEresults_RTTvsIC_gx.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/geneAnnotation.RData")
load("D:/RTTproject/CellAnalysis/Data/SampleInfo.RData")


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




g <- graph_from_adjacency_matrix(graph_matrix, 
                                 mode = "lower", 
                                 weighted = "weight", 
                                 diag = FALSE)

createNetworkFromIgraph(
  g,
  title = "GO relations1",
  collection = "test"
)




