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
setwd("D:/RTTproject/CellAnalysis")
load("D:/RTTproject/CellAnalysis/Data/SampleInfo.RData")

#******************************************************************************#
# Collect proteomics data
#******************************************************************************#

# Load data
load("Proteins/Preprocessing/pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load("Proteins/Preprocessing/DEresults_px.RData")
load("Proteins/Preprocessing/proteinAnnotation.RData")
all(annotations$uniprot_gn_id %in% str_remove(rownames(pxMatrix_imp), "-.*"))
load("Data/sampleInfo.RData")

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
load("Genes/Preprocessing/gxMatrix_norm.RData")
load("Genes/Preprocessing/DEresults_RTTvsIC_gx.RData")
load("Genes/Preprocessing/geneAnnotation.RData")

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

#biomaRt::listAttributes(ensembl)
px_gx_annotations <- getBM(attributes=c("uniprot_gn_id",
                                   "ensembl_gene_id"), 
                      filters = 'uniprot_gn_id',
                      values = str_remove(rownames(pxMatrix_imp), "-.*"),
                      mart = ensembl)

geneAnnotation_fil <- unique(geneAnnotation[geneAnnotation$gene_id %in% rownames(gxMatrix_norm),1:3])
proteinAnnotation_fil <- unique(annotations[annotations$ID %in% rownames(pxMatrix_imp),c(1,2,5)])

px_gx_annotations <- full_join(px_gx_annotations, geneAnnotation_fil,
                                by = c("ensembl_gene_id" = "EnsemblID"))

px_gx_annotations <- full_join(px_gx_annotations, proteinAnnotation_fil,
                                by = c("uniprot_gn_id" = "uniprot_gn_id"))


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

save(px_gx_annotations, file = "Proteins/Preprocessing/px_gx_annotations1.RData")

load("Proteins/Preprocessing/px_gx_annotations.RData")
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



#******************************************************************************#
# Make venn diagram
#******************************************************************************#
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")
load("Proteins/Preprocessing/px_gx_annotations.RData")

px_gx_annotations <- unique(px_gx_annotations[,c("gene_id", "ID")])

sum(is.na(px_gx_annotations$gene_id))
length(unique(px_gx_annotations$ID[is.na(px_gx_annotations$gene_id)]))

sum(is.na(px_gx_annotations$ID))
length(unique(px_gx_annotations$gene_id[is.na(px_gx_annotations$ID)]))

sum(!is.na(px_gx_annotations$gene_id) & !is.na(px_gx_annotations$ID))



overlap <- px_gx_annotations[!is.na(px_gx_annotations$gene_id) & !is.na(px_gx_annotations$ID),]

overlap_px_expr <- pxMatrix_imp[overlap$ID,]
overlap_gx_expr <- gxMatrix_norm[overlap$gene_id,]

corMatrix <- matrix(NA, nrow = 42, ncol = 42)
for (i in 1:ncol(overlap_gx_expr)){
  corMatrix[,i] <- apply(overlap_px_expr,2,function(x){cor(x,overlap_gx_expr[,i])})
}

colnames(corMatrix) <- colnames(overlap_gx_expr)
rownames(corMatrix) <- colnames(overlap_gx_expr)

plotDF <- gather(data.frame(corMatrix))
plotDF$sample <- rep(rownames(corMatrix), ncol(corMatrix))
colnames(plotDF) <- c("transcriptomics", "value", "proteomics")


plotDF <- arrange(plotDF, by = transcriptomics)

plotDF <- inner_join(plotDF, sampleInfo, by = c("transcriptomics" = "SampleID"))

# Set colors
group_colors <- setNames(c("#377EB8","#E41A1C"),
                         c("IC","RTT"))

time_colors <- setNames(c("#FEE5D9","#FCAE91","#FB6A4A","#CB181D"),
                         c("D0", "D13", "D40", "D75"))

tissue_colors <- setNames(c("#1B9E77","#D95F02","#7570B3"),
                          c("Dorsal", "Cell", "Ventral"))


main <- ggplot(plotDF) +
  geom_tile(aes(x = transcriptomics, y = proteomics, fill = value), width = 0.95, height = 0.95) +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "none")

tissue_hor <- ggplot(plotDF) +
  geom_tile(aes(x = transcriptomics, y = 1, fill = Tissue)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = tissue_colors)

tissue_ver <- ggplot(plotDF) +
  geom_tile(aes(y = transcriptomics, x = 1, fill = Tissue)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = tissue_colors)

time_hor <- ggplot(plotDF) +
  geom_tile(aes(x = transcriptomics, y = 1, fill = Time)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = time_colors)

time_ver <- ggplot(plotDF) +
  geom_tile(aes(y = transcriptomics, x = 1, fill = Time)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = time_colors)

group_hor <- ggplot(plotDF) +
  geom_tile(aes(x = transcriptomics, y = 1, fill = Group)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = group_colors)

group_ver <- ggplot(plotDF) +
  geom_tile(aes(y = transcriptomics, x = 1, fill = Group)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = group_colors)


library(ggpubr)

p <- ggarrange(NULL, NULL, NULL, group_hor,
          NULL, NULL, NULL, tissue_hor,
          NULL, NULL, NULL, time_hor,
          group_ver, tissue_ver, time_ver, main,
          ncol = 4, nrow = 4, widths = c(0.5,0.5,0.5,9),
          heights = c(0.5,0.5,0.5,9)
          )

ggsave(p, file = "Proteins/Preprocessing/px_gx_correlations.png",
       width = 8, height =  8)


# Get legends:
legendPlot <- ggplot(plotDF) +
  geom_tile(aes(x = transcriptomics, y = proteomics, fill = value), width = 0.95, height = 0.95) +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "right",
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))


legend <- cowplot::get_legend(legendPlot)

ggsave(legend, file = "Proteins/Preprocessing/legend_px_gx_cor.png", width = 8, height = 8)

#plotDF <- plotDF[!duplicated(plotDF$corID),]
plotDF$key <- factor(plotDF$transcriptomics, levels = sort(unique(plotDF$transcriptomics)))
plotDF$sample <- factor(plotDF$sample, levels = sort(unique(plotDF$sample)))



plotDF$corID <- NA
for (j in 1:nrow(plotDF)){
  plotDF[j, "corID"] <-  paste(sort(c(plotDF[j,"key"],plotDF[j,"sample"]))[1],
                               sort(c(plotDF[j,"key"],plotDF[j,"sample"]))[2], sep = "_")
}
