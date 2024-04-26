
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Eigengene")

# load packages
library(data.table)
library(tidyverse)

# Read GWAS files
GWASdata <- NULL
for (file in list.files("GWAS", full.names = TRUE)){
  temp <- fread(file)
  gene <- str_remove(str_remove(file, ".*/"),"_.*")
  temp$Gene <- gene
  
  GWASdata <- rbind.data.frame(GWASdata, temp)
}

write.csv(GWASdata, file = "GWASdata.csv")

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

library(RCy3)
library(igraph)

edges <- plotDF[,c("Gene", "efoTraits_abbreviated", "Count")]
colnames(edges) <- c("from", "to", "count")

nodes <- unique(data.frame(name = c(plotDF$efoTraits_abbreviated, plotDF$Gene),
                           group = c(plotDF$Category, rep("Gene",nrow(plotDF)))))


g <- graph_from_data_frame(edges, directed=TRUE, vertices=nodes)

createNetworkFromIgraph(
  g,
  title = "GWASnetwork2",
  collection = "GWAS"
)


