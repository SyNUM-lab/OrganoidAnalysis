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
library(ggrepel)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/2. Proteomics/1. Preprocessing"))

# Load sample information
load(paste0(homeDir,"/sampleInfo.RData"))

#******************************************************************************#
# Collect proteomics data
#******************************************************************************#

# Load data
load("pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load("DEresults_px.RData")
load("proteinAnnotation.RData")
all(annotations$uniprot_gn_id %in% str_remove(rownames(pxMatrix_imp), "-.*"))

# Get p-values
pvalues_px <- DEresults_px[,3:9]
rownames(pvalues_px) <- DEresults_px$name
colnames(pvalues_px) <- c("Cell_D0",
                          "Dorsal_D13",
                          "Dorsal_D40",
                          "Dorsal_D75",
                          "Ventral_D13",
                          "Ventral_D40",
                          "Ventral_D75")

# get adjusted p-values
adjpvalues_px <- pvalues_px
for (i in 1:ncol(pvalues_px)){
  adjpvalues_px[,i] <- p.adjust(pvalues_px[,i], "fdr")
}

# Get logFCs
logFCs_px <- DEresults_px[,str_detect(colnames(DEresults_px), "ratio")]
rownames(logFCs_px) <- DEresults_px$name
colnames(logFCs_px) <- c("Cell_D0",
                          "Dorsal_D13",
                          "Dorsal_D40",
                          "Dorsal_D75",
                          "Ventral_D13",
                          "Ventral_D40",
                          "Ventral_D75")

#******************************************************************************#
# Collect transcriptomics data
#******************************************************************************#

# Load data
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
genes_all <- rownames(topList[[1]])

# Get p-values
pvalues_gx <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues_gx) <- genes_all
colnames(pvalues_gx) <- names(topList)

# get adjusted p-values
adj_pvalues_gx <- matrix(unlist(lapply(topList, function(x){x[genes_all,"FDR"]})), ncol = length(topList))
rownames(adj_pvalues_gx) <- genes_all
colnames(adj_pvalues_gx) <- names(topList)

# Get logFCs
logFCs_gx <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
rownames(logFCs_gx) <- genes_all
colnames(logFCs_gx) <- names(topList)


#******************************************************************************#
# Combine proteomics and transcriptomics
#******************************************************************************#

# uniprot to ensembl IDs
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
px_gx_annotations <- getBM(attributes=c("uniprot_gn_id",
                                   "ensembl_gene_id"), 
                      filters = 'uniprot_gn_id',
                      values = str_remove(rownames(pxMatrix_imp), "-.*"),
                      mart = ensembl)

# Add gene annotation
geneAnnotation_fil <- unique(geneAnnotation[geneAnnotation$gene_id %in% rownames(gxMatrix_norm),1:3])
px_gx_annotations <- full_join(px_gx_annotations, geneAnnotation_fil,
                               by = c("ensembl_gene_id" = "EnsemblID"))

# Add protein annotation
proteinAnnotation_fil <- unique(annotations[annotations$ID %in% rownames(pxMatrix_imp),c(1,2,5)])
px_gx_annotations <- full_join(px_gx_annotations, proteinAnnotation_fil,
                                by = c("uniprot_gn_id" = "uniprot_gn_id"))

# Remove empty rows
px_gx_annotations <- px_gx_annotations[(px_gx_annotations$gene_id %in% rownames(gxMatrix_norm)) |
                                           is.na(px_gx_annotations$gene_id),]

px_gx_annotations <- px_gx_annotations[(px_gx_annotations$ID %in% rownames(pxMatrix_imp)) |
                                           is.na(px_gx_annotations$ID),]

# Some proteins are linked to NA genes

# duplicated proteins with at least one gene id available
duplicatedProteinIDs <- unique(intersect(px_gx_annotations$ID[(duplicated(px_gx_annotations$ID))], 
                                         px_gx_annotations$ID[(!is.na(px_gx_annotations$gene_id))]))
duplicatedProteinIDs <- duplicatedProteinIDs[!is.na(duplicatedProteinIDs)]

# Remove these duplicated proteins without gene ID
px_gx_annotations <- px_gx_annotations[!((is.na(px_gx_annotations$gene_id)) &
                                         (px_gx_annotations$ID %in% duplicatedProteinIDs)),]

# Save annotation file
save(px_gx_annotations, file = "px_gx_annotations1.RData")


#******************************************************************************#
# Get statistics
#******************************************************************************#

# Load data
load("px_gx_annotations1.RData")

# Total number of genes that pass QC
length(unique(px_gx_annotations$gene_id[!is.na(px_gx_annotations$gene_id)]))

# Total number of proteins that pass QC
length(unique(px_gx_annotations$ID[!is.na(px_gx_annotations$ID)]))

# Total number of genes with protein ID
length(unique(px_gx_annotations$gene_id[(!is.na(px_gx_annotations$ID)) & (!is.na(px_gx_annotations$gene_id))]))

# Total number of genes w/o protein ID
length(unique(px_gx_annotations$gene_id[is.na(px_gx_annotations$ID) & (!is.na(px_gx_annotations$gene_id))]))

# Total number of proteins with gene ID
length(unique(px_gx_annotations$ID[(!is.na(px_gx_annotations$gene_id)) & (!is.na(px_gx_annotations$ID))]))

# Total number of proteins w/o gene ID
length(unique(px_gx_annotations$ID[(is.na(px_gx_annotations$gene_id)) & (!is.na(px_gx_annotations$ID))]))
