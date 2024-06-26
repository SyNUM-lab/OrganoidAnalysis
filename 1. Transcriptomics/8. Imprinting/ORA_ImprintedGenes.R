# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(readxl)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/8. Imprinting"))

load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_norm.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/geneAnnotation.RData"))
all_genes <- unique(geneAnnotation$GeneName[geneAnnotation$gene_id %in% rownames(gxMatrix_norm)])
all_genes <- all_genes[!is.na(all_genes)]

# Get imprinted genes: https://www.geneimprint.com/site/genes-by-species
imprintDF <- read_xlsx("ImprintedGenes.xlsx")
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,3)],]
imprinted_genes <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[3])]
imprinted_genes <- imprinted_genes[!is.na(imprinted_genes)]
imprinted_genes <- intersect(imprinted_genes, all_genes)

# get non-imprinted genes
nimprinted_genes <- setdiff(all_genes, imprinted_genes)
nimprinted_genes <- nimprinted_genes[!is.na(nimprinted_genes)]


load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_SYMBOL_Hs.RData"))
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")


test <- enricher(
  gene = imprinted_genes,
  pvalueCutoff = Inf,
  pAdjustMethod = "fdr",
  universe = all_genes,
  minGSSize = -Inf,
  maxGSSize = Inf,
  qvalueCutoff = Inf,
  gson = NULL,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)

results <- test@result
