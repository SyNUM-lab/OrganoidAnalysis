# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(data.table)
library(biomaRt)
library(edgeR)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/7. Splicing"))

# Load data
load("Data/IsoPctMatrix.RData")
load("Data/txMatrix_raw.RData")
load("Data/txAnnotation.RData")
load(paste0(homeDir,"/sampleInfo.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_raw1.RData"))

all_transcripts <- unlist(lapply(str_split(rownames(txMatrix_raw), "_"), function(x)x[[1]]))
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# Load genome annotations for each chromosome
egs = getBM(attributes = c('ensembl_transcript_id',
                           'ensembl_gene_id',
                           "hgnc_symbol",
                           'transcription_start_site'), 
            filters='ensembl_transcript_id',
            values= all_transcripts,
            mart=ensembl)

# Get transcription start site
egs$TSS = ifelse(egs$strand == "1", egs$start_position, egs$end_position)