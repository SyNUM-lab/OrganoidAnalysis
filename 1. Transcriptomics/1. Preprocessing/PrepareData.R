# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(data.table)
library(biomaRt)

# Directory with RSEM output
dataDir <- "D:/RTTproject/GeneticFiles/RSEM/a_Part"

# Working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing"))

# Get raw gene expression values:
for (i in 1:7){
  
  # Get gene expression files
  files <- list.files(path = paste0(dataDir,i), pattern = ".genes.results")
  
  # Get sample names
  sampleNames <- str_remove(files, "_Clean_Data_unaligned.genes.results")
  
  for (j in 1:length(files)){
    
    # Read gene expression file
    GeneExpr <- fread(paste0(dataDir,i,"/",files[j]))
    print(paste0(i,",",j,": ", nrow(GeneExpr)))
    
    # Retrieve expression value from file
    GeneExpr <- GeneExpr[,c("gene_id", "expected_count")]
    
    # Set column name to sample name
    colnames(GeneExpr) <- c("gene_id",sampleNames[j])
    
    # Combine with the already-collected gene expression values
    if ((i == 1)&(j==1)){
      GeneExpr_all <- GeneExpr
    } else{
      GeneExpr_all <- inner_join(GeneExpr_all,GeneExpr, by = c("gene_id" = "gene_id"))
    }
  }
}

# Format data
gxMatrix_raw <- as.matrix(GeneExpr_all[,-1])
rownames(gxMatrix_raw) <- GeneExpr_all$gene_id

# Save expression matrix
save(gxMatrix_raw, file = "gxMatrix_raw1.RData")

# Make gene annotation file:

# Read gene expression file
GeneExpr <- fread(paste0(dataDir,i,"/",files[j]))

# Retrieve relevant columns
geneAnnotation <- GeneExpr[,c("gene_id", "transcript_id(s)", "length", "effective_length")]

# Add gene name and Ensembl ID to the file
geneAnnotation$EnsemblID <- unlist(lapply(str_split(geneAnnotation$gene_id, "_"), function(x){x[1]}))
geneAnnotation$GeneName <- unlist(lapply(str_split(geneAnnotation$gene_id, "_"), function(x){x[2]}))

# Change order of columns
geneAnnotation <- geneAnnotation[,c("gene_id", "EnsemblID", "GeneName", "transcript_id(s)", "length", "effective_length")]

# Add entrez gene id and chromosome name to the gene annotation
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "entrezgene_id",
                                  "chromosome_name"), 
                     filters = 'ensembl_gene_id',
                     values = geneAnnotation$EnsemblID,
                     mart = ensembl)

geneAnnotation <- inner_join(geneAnnotation, annotations, by = c("EnsemblID" = "ensembl_gene_id"))

# Save file
save(geneAnnotation, file = "geneAnnotation.RData")


