# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(data.table)
library(biomaRt)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/7. Splicing"))

# Data directory (contains RSEM output)
data_dir <- paste0("D:/RTTproject/GeneticFiles/RSEM/a_Part", 1:7)

# Get raw transcript expression values:
for (i in 1:length(data_dir)){
  
  # Get gene expression files
  files <- list.files(path = data_dir[i], pattern = ".isoforms.results")
  
  # Get sample names
  sampleNames <- str_remove(files, "_Clean_Data_unaligned.isoforms.results")
  
  for (j in 1:length(files)){
    
    # Read gene expression file
    GeneExpr <- fread(paste0(data_dir[i],"/",files[j]))
    print(paste0(i,",",j,": ", nrow(GeneExpr)))
    # Retrieve expression value from file
    GeneExpr <- GeneExpr[,c("transcript_id", "expected_count")]
    
    # Set column name to sample name
    colnames(GeneExpr) <- c("transcript_id",sampleNames[j])
    
    # Combine with the already-collected gene expression values
    if ((i == 1)&(j==1)){
      GeneExpr_all <- GeneExpr
    } else{
      GeneExpr_all <- inner_join(GeneExpr_all,GeneExpr, by = c("transcript_id" = "transcript_id"))
    }
  }
}

# Format data
txMatrix_raw <- as.matrix(GeneExpr_all[,-1])
rownames(txMatrix_raw) <- GeneExpr_all$transcript_id
save(txMatrix_raw, file = "Data/txMatrix_raw.RData")


# Get raw transcript percentage values:
for (i in 1:length(data_dir)){
  
  # Get gene expression files
  files <- list.files(path = data_dir[i], pattern = ".isoforms.results")
  
  # Get sample names
  sampleNames <- str_remove(files, "_Clean_Data_unaligned.isoforms.results")
  
  for (j in 1:length(files)){
    
    # Read gene expression file
    GeneExpr <- fread(paste0(data_dir[i],"/",files[j]))
    print(paste0(i,",",j,": ", nrow(GeneExpr)))
    # Retrieve expression value from file
    GeneExpr <- GeneExpr[,c("transcript_id", "IsoPct")]
    
    # Set column name to sample name
    colnames(GeneExpr) <- c("transcript_id",sampleNames[j])
    
    # Combine with the already-collected gene expression values
    if ((i == 1)&(j==1)){
      GeneExpr_all <- GeneExpr
    } else{
      GeneExpr_all <- inner_join(GeneExpr_all,GeneExpr, by = c("transcript_id" = "transcript_id"))
    }
  }
}

# Format data
IsoPctMatrix <- as.matrix(GeneExpr_all[,-1])
rownames(IsoPctMatrix) <- GeneExpr_all$transcript_id
save(IsoPctMatrix, file = "Data/IsoPctMatrix.RData")

# Transcript annotation
i = 1
j = 1
files <- list.files(path = data_dir[i], pattern = ".isoforms.results")
txAnnotation <- fread(paste0(data_dir[i],"/",files[j]))
txAnnotation <- txAnnotation[,c("transcript_id","gene_id","length","effective_length")]
save(txAnnotation, file = "Data/txAnnotation.RData")

