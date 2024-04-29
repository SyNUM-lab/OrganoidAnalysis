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

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/5. GSEA")

# Load necessary data
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData")
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOannotation.RData")

load("Data/terms_ordered1.RData")
load("Data/GOresults_GSEA_gx1.RData")
load("Data/GOresults_NES_GSEA_gx1.RData")
load("Data/reducedTerms_RTTvsIC_BP1.RData")

preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/SampleInfo.RData")
genes_all <- rownames(topList[[1]])

# Get p-values
pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues) <- unlist(lapply(str_split(genes_all, "_"), function(x){x[1]}))
colnames(pvalues) <- names(topList)

# get adjusted p-values
adj_pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"FDR"]})), ncol = length(topList))
rownames(adj_pvalues) <- unlist(lapply(str_split(genes_all, "_"), function(x){x[1]}))
colnames(adj_pvalues) <- names(topList)

# Get logFCs
logFCs <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
rownames(logFCs) <- unlist(lapply(str_split(genes_all, "_"), function(x){x[1]}))
colnames(logFCs) <- names(topList)

# Get adjusted p-values
GOresults_adj <- apply(GOresults[,2:8],2,function(x){p.adjust(x, method = "fdr")})
rownames(GOresults_adj) <- GOresults$Name

# Get GO terms
rownames(GOannotation) <- GOannotation$Name
GOterms_all <- GOannotation[GOresults$Name, "ID"]

# Reduce number of GO terms based on previous selection
GOannotation_fil <- GOannotation[GOannotation$ID %in% unique(reducedTerms$parent),]
GOresults_fil <- GOresults[GOresults$Name %in% GOannotation_fil$Name,]
GOresults_NES_fil <- GOresults_NES[GOresults_NES$Name %in% GOannotation_fil$Name,]
GOresults_adj_fil <- GOresults_adj[rownames(GOresults_adj) %in% GOannotation_fil$Name,]

# Get signed -log10 pvalue
values <- -log10(GOresults[,2:8]) * sign(GOresults_NES[,2:8])
rownames(values) <- GOresults$Name


# Get time- and region-specific markers
GOresults_NES_fil <- GOresults_NES_fil[,-1]
GOresults_fil <- GOresults_fil[,-1]
sigMatrix_pos <- matrix((GOresults_NES_fil > 0) & (GOresults_adj_fil < 0.05), ncol = ncol(GOresults_NES_fil))
sigMatrix_neg <- matrix((GOresults_NES_fil < 0) & (GOresults_adj_fil < 0.05), ncol = ncol(GOresults_NES_fil))
sigMatrix_all <- matrix((GOresults_adj_fil < 0.05), ncol = ncol(GOresults_adj_fil))

# Early-time markers
early <- c(rownames(GOresults_adj_fil)[(rowSums(sigMatrix_pos[,c(1,2,5)]) == 3) &
                                  (rowSums(sigMatrix_all[,c(3,4,6,7)]) ==0)],
           rownames(GOresults_adj_fil)[(rowSums(sigMatrix_neg[,c(1,2,5)]) == 3) &
                                  (rowSums(sigMatrix_all[,c(3,4,6,7)]) ==0)]
)

# Mid-time markers
mid <- c(rownames(GOresults_adj_fil)[(rowSums(sigMatrix_pos[,c(2,3,5,6)]) >= 3) &
                                (rowSums(sigMatrix_all[,c(1,4,7)]) ==0)],
         rownames(GOresults_adj_fil)[(rowSums(sigMatrix_neg[,c(2,3,5,6)]) >= 3) &
                                (rowSums(sigMatrix_all[,c(1,4,7)]) ==0)]
)

# Late-time markers
late <- c(rownames(GOresults_adj_fil)[(rowSums(sigMatrix_pos[,c(3,4,6,7)]) >= 3) &
                                 (rowSums(sigMatrix_all[,c(1,2,5)]) ==0)],
          rownames(GOresults_adj_fil)[(rowSums(sigMatrix_neg[,c(3,4,6,7)]) >= 3) &
                                 (rowSums(sigMatrix_all[,c(1,2,5)]) ==0)]
)

# Vental-region markers
ventral <- c(rownames(GOresults_adj_fil)[(rowSums(sigMatrix_pos[,c(1,5,6,7)]) >= 3) &
                                    (rowSums(sigMatrix_all[,c(2,3,4)]) == 0)],
             rownames(GOresults_adj_fil)[(rowSums(sigMatrix_neg[,c(1,5,6,7)]) >= 3) &
                                    (rowSums(sigMatrix_all[,c(2,3,4)]) == 0)]
)

# Dorsal-region markers
dorsal <- c(rownames(GOresults_adj_fil)[(rowSums(sigMatrix_pos[,c(1,2,3,4)]) >= 3) &
                                   (rowSums(sigMatrix_all[,c(5,6,7)]) ==0)],
            rownames(GOresults_adj_fil)[(rowSums(sigMatrix_neg[,c(1,2,3,4)]) >= 3) &
                                   (rowSums(sigMatrix_all[,c(5,6,7)]) ==0)]
)


# Collect all selected GO terms
selTerms <- rev(c(dorsal, ventral, early, mid, late))
values_fil <- values[selTerms,]
GOresults_adj_fil <- GOresults_adj[selTerms,]

# Prepare data for plotting
plotData <- gather(values_fil)
plotData$Name <- rep(rownames(values_fil),7)

# combine with GO annotaiton
plotData <- inner_join(plotData, GOannotation, by = c('Name' = 'Name'))

# Combine with adjusted pvalues
GOresults_adj_fil <- as.data.frame(GOresults_adj_fil)
GOresults_adj_fil <- gather(GOresults_adj_fil[,1:7])
plotData <- cbind.data.frame(plotData, GOresults_adj_fil$value)
plotData$Sig <- ifelse(plotData$`GOresults_adj_fil$value` < 0.05, "Yes", "No")

# Add sample info
plotData$Region <- NA
plotData$Region[str_detect(plotData$key, "Cell")] <- "iPSC"
plotData$Region[str_detect(plotData$key, "Dorsal")] <- "Dorsal"
plotData$Region[str_detect(plotData$key, "Ventral")] <- "Ventral"

plotData$Time <- "D0"
plotData$Time[str_detect(plotData$key, "D13")] <- "D13"
plotData$Time[str_detect(plotData$key, "D40")] <- "D40"
plotData$Time[str_detect(plotData$key, "D75")] <- "D75"

plotData$key <- factor(plotData$key,
                       levels = c("Dorsal_D75_RTTvsIC",
                                  "Dorsal_D40_RTTvsIC",
                                  "Dorsal_D13_RTTvsIC",
                                  "Cell_D0_RTTvsIC",
                                  "Ventral_D13_RTTvsIC",
                                  "Ventral_D40_RTTvsIC",
                                  "Ventral_D75_RTTvsIC"))

plotData$Description <- factor(firstup(plotData$Description), 
                               levels = firstup(GOannotation[selTerms, "Description"]))
plotData$Name <- factor(plotData$Name, levels = selTerms)


# Make main plot: Heatmap
mainPlot <- ggplot(data = plotData, aes(x = key, y = Description, fill = value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(.~Region, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "Signed\n-log10 p-value") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotData$key,
  time = plotData$Time,
  tissue = plotData$Region))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# Tissue/region
colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[(str_detect(colSideColor$sample, "D40") |
                                   str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots into a single figure
finalPlot <- ggarrange(colSideColorPlot_tissue,
                       mainPlot,
                       colSideColorPlot_time,
                       heights = c(0.5,9,0.5), widths = c(2,8),nrow = 3,ncol = 1,
                       align = "hv",
                       common.legend = FALSE)

# Save plot
ggsave(finalPlot, file = "Plots/GOEnrichment_TimePointSpecific.png", width = 8, height = 6)


# Get legends:
legendPlot <- ggplot(data = plotData, aes(x = key, y = Description, fill = value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(.~Region, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "") +
  guides(color = "none") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))


legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "Plots/legend_GOEnrichment_TimePointSpecific.png", width = 8, height = 8)

