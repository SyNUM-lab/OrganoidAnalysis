
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# FUNCTION: capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/2. Proteomics/4. GSEA"))

# Load data
preprocessing_dir <- paste0(homeDir, "/2. Proteomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"pxData_imp.RData"))
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load(paste0(preprocessing_dir,"DEresults_px.RData"))
load(paste0(preprocessing_dir,"proteinAnnotation.RData"))
all(annotations$uniprot_gn_id %in% str_remove(rownames(pxMatrix_imp), "-.*"))
load(paste0(homeDir, "/sampleInfo.RData"))

# Get p-values
pvalues <- DEresults_px[,3:9]
rownames(pvalues) <- DEresults_px$name

# get adjusted p-values
adjpvalues <- pvalues
for (i in 1:ncol(pvalues)){
  adjpvalues[,i] <- p.adjust(pvalues[,i], "fdr")
}

# Get logFCs
logFCs <- DEresults_px[,str_detect(colnames(DEresults_px), "ratio")]
rownames(logFCs) <- DEresults_px$name

###############################################################################

# GSEA

###############################################################################

# Perform GSEA for each time point/region
for (i in 1:ncol(adjpvalues)){
  gsea_input <- -log(pvalues[,i])*sign(logFCs[,i])
  names(gsea_input) <- str_remove(rownames(pvalues), "-.*")
  keep <- sort(abs(gsea_input), decreasing = TRUE)
  keep <- names(keep[!duplicated(names(keep))])
  gsea_input <- gsea_input[keep]
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  
  set.seed(123)
  GOtest <- gseGO(
    geneList = gsea_input,
    keyType = "UNIPROT",
    OrgDb = org.Hs.eg.db,
    ont =  "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500)
  
  if (i == 1){
    GOresults <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                          GOtest@result$ID, ")"),
                            pvalue = GOtest@result$pvalue)
    GOresults_NES <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                          GOtest@result$ID, ")"),
                            NES = GOtest@result$NES)
  }
  if (i > 1){
    temp <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                     GOtest@result$ID, ")"),
                       pvalue = GOtest@result$pvalue)
    GOresults <- full_join(GOresults,temp, by = c("Name" = "Name"))
    
    temp_NES <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                     GOtest@result$ID, ")"),
                       NES = GOtest@result$NES)
    GOresults_NES <- full_join(GOresults_NES,temp_NES, by = c("Name" = "Name"))
  }
  
}

# Save GSEA results
save(GOresults, file = "GOresults_GSEA.RData")
save(GOresults_NES, file = "GOresults_NES_GSEA.RData")
