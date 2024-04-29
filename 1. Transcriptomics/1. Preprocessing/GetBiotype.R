# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(biomaRt)

# Set working directory
setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing")

# Load data
load("geneAnnotation.RData")

# Set biomaRt datatset
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# Get annotations
annotations <- getBM(attributes=c('ensembl_gene_id',
                                  "gene_biotype"),
                     filters = 'ensembl_gene_id',
                     values = geneAnnotation$EnsemblID,
                     mart = ensembl)

annotations <- annotations[!duplicated(annotations),]
table(annotations$gene_biotype)

# save annotations
save(annotations, file = "Gene_biotype.RData")
