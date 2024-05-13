
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/6. Regulators"))

# load packages
library(data.table)
library(tidyverse)
library(RCy3)
library(igraph)

# Read GWAS files
GWASdata <- NULL
for (file in list.files("GWAS", full.names = TRUE)){
  temp <- fread(file)
  gene <- str_remove(str_remove(file, ".*/"),"_.*")
  temp$Gene <- gene
  
  GWASdata <- rbind.data.frame(GWASdata, temp)
}

# Write output
write.csv(GWASdata, file = "Data/GWASdata.csv")

# Add categories manual
plotDF <- fread("GWASdata.csv")
plotDF <- plotDF[plotDF$Category != "Other",]
plotDF <- plotDF %>%
  group_by(Gene,efoTraits_abbreviated) %>%
  summarise(efoTraits_abbreviated = efoTraits_abbreviated,
            Gene = Gene,
            Category = Category,
            Count = sum(associationCount))

plotDF <- unique(plotDF)
plotDF <- plotDF[plotDF$Category != "Anthropometry",]
plotDF <- plotDF[plotDF$Gene != "BASP1-AS1",]


# Make graph:

# 1) set edges
edges <- plotDF[,c("Gene", "efoTraits_abbreviated", "Count")]
colnames(edges) <- c("from", "to", "count")

# 2) set nodes
nodes <- unique(data.frame(name = c(plotDF$efoTraits_abbreviated, plotDF$Gene),
                           group = c(plotDF$Category, rep("Gene",nrow(plotDF)))))

# 3) Make igraph object
g <- graph_from_data_frame(edges, directed=TRUE, vertices=nodes)

# 4) export graph to cytoscape
createNetworkFromIgraph(
  g,
  title = "GWASnetwork2",
  collection = "GWAS"
)


