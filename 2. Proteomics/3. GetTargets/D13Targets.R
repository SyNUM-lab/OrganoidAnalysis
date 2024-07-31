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
library(grid)
library(grid)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/2. Proteomics/3. GetTargets"))

# Load data
preprocessing_dir <- paste0(homeDir,"/2. Proteomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"pxData_imp.RData"))
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load(paste0(preprocessing_dir,"DEresults_px.RData"))
load(paste0(preprocessing_dir,"proteinAnnotation.RData"))
load(paste0(homeDir,"/sampleInfo.RData"))

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


sigMatrix_logFC_p <- matrix((abs(logFCs) > 1) & (adjpvalues < 0.01), ncol = ncol(logFCs))
sigMatrix_p <- matrix(pvalues < 0.1, ncol = ncol(logFCs))

proteinTargets <- rownames(adjpvalues)[(rowSums(sigMatrix_p[,c(3,4,6,7)]) == 0) & (rowSums(sigMatrix_logFC_p[,c(2,5)]) > 0)]
proteinTargetNames <- annotations$hgnc_symbol[annotations$ID %in% proteinTargets]

# Make plot
pvalues_sel <- adjpvalues[proteinTargets,]
plotP <- gather(pvalues_sel, value = "FDR")
plotP$ProteinID <- rep(rownames(pvalues_sel), ncol(pvalues_sel))
plotP$key <- str_remove(plotP$key, "_p.val")
plotP$ID <- paste0(plotP$key, "_",plotP$ProteinID)

logfc_sel <- logFCs[proteinTargets,]
plotFC <- gather(logfc_sel, value = "log2FC")
plotFC$ProteinID <- rep(rownames(logfc_sel), ncol(logfc_sel))
plotFC$key <- str_remove(plotFC$key, "_ratio")
plotFC$ID <- paste0(plotFC$key, "_",plotFC$ProteinID)

# combine data
plotData <- inner_join(plotFC, plotP[,c(2,4)], by = c("ID" = "ID"))

# Add a column that indicates statistical signficance (FDR < 0.05)
plotData$Sig <- ifelse((plotData$FDR < 0.05) & (abs(plotData$log2FC) > 0.58), "Yes", "No")

# Combine with feature info
#plotData <- inner_join(plotData, selection, by = c("ProteinID" = "ID"))
plotData$Time <- NA
plotData$Time[str_detect(plotData$key, "D0")] <- "D0"
plotData$Time[str_detect(plotData$key, "D13")] <- "D13"
plotData$Time[str_detect(plotData$key, "D40")] <- "D40"
plotData$Time[str_detect(plotData$key, "D75")] <- "D75"

plotData$Region <- NA
plotData$Region[str_detect(plotData$key, "Dorsal")] <- "Dorsal"
plotData$Region[str_detect(plotData$key, "Ventral")] <- "Ventral"
plotData$Region[str_detect(plotData$key, "Cell")] <- "iPSC"
plotData$Region <- factor(plotData$Region, levels = c("Dorsal", "iPSC", "Ventral"))

plotData$Group <- paste0(plotData$Region, "_",plotData$Time)
plotData$Group <- factor(plotData$Group,levels = c("Dorsal_D75",
                                                   "Dorsal_D40",
                                                   "Dorsal_D13",
                                                   "iPSC_D0",
                                                   "Ventral_D13",
                                                   "Ventral_D40",
                                                   "Ventral_D75"))

plotData <- left_join(plotData, annotations, by = c("ProteinID" = "ID"))

# Perform clustering
model <- hclust(dist(logfc_sel), "ward.D2")
proteins_ordered <- model$labels[model$order]
plotData$ProteinID <- factor(plotData$ProteinID, 
                             levels = proteins_ordered)
annotations_fil <- annotations[annotations$ID %in% unique(plotData$ProteinID),]
rownames(annotations_fil) <- annotations_fil$ID
names_ordered <- annotations_fil[proteins_ordered, "hgnc_symbol"]
plotData$hgnc_symbol <- factor(plotData$hgnc_symbol, levels = names_ordered)

plotData$CompleteName <- factor(paste0(plotData$hgnc_symbol, " (", plotData$ProteinID, ")"),
                                levels = paste0(names_ordered, " (", proteins_ordered, ")"))


# Make main plot: Heatmap
mainPlot <- ggplot(data = plotData, aes(x = Group, y = CompleteName, fill = log2FC,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(cols = vars(Region), scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-6,6),trans = "pseudo_log", oob = scales::squish) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "logFC") +
  ylab("Protein Name\n") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, hjust = 1),
        #axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())



# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotData$Group,
  time = plotData$Time,
  region = plotData$Region))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
  facet_grid(.~region, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# Tissue/region
colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = region)) +
  geom_text(data = colSideColor[(str_detect(colSideColor$sample, "D40") |
                                   str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = region)) +
  facet_grid(.~region, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")



# Combine plots
finalPlot <-  
  colSideColorPlot_tissue +
  mainPlot +
  colSideColorPlot_time +
  plot_layout(ncol = 1, nrow = 3, 
              heights = c(0.5,9,0.5))


# Save plot
ggsave(finalPlot, file = "D13Targets.png", width = 7.5, height = 10)


# Make legend
legendPlot <- ggplot(data = plotData, aes(x = Group, y = CompleteName, fill = log2FC,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(cols = vars(Region), scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-6,6),trans = "pseudo_log", oob = scales::squish) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "logFC") +
  ylab("Protein Name\n") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, hjust = 1),
        #axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))



legend <- cowplot::get_legend(legendPlot)
ggsave(legend, file = "D13Targets_legend.png", width = 8, height = 8)





# Working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics//4. DEAnalysis"))

# Load data
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(homeDir,"/SampleInfo.RData"))

genes_all <- rownames(topList[[1]])

# Get p-values
pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues) <- genes_all
colnames(pvalues) <- names(topList)


# get adjusted p-values
adjpvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"FDR"]})), ncol = length(topList))
rownames(adjpvalues) <- genes_all
colnames(adjpvalues) <- names(topList)

# Get logFCs
logFCs <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
rownames(logFCs) <- genes_all
colnames(logFCs) <- names(topList)



sigMatrix_logFC_p <- matrix((abs(logFCs) > 0) & (pvalues < 0.05), ncol = ncol(logFCs))
sigMatrix_p <- matrix(pvalues < 0.1, ncol = ncol(logFCs))

geneTargets <- rownames(adjpvalues)[(rowSums(sigMatrix_p[,c(3,4,6,7)]) == 0) & (rowSums(sigMatrix_logFC_p[,c(2,5)]) > 0)]
geneTargetNames <- geneAnnotation$GeneName[geneAnnotation$gene_id %in% geneTargets]


targets_all <- intersect(proteinTargetNames, geneTargetNames)


#EIF2B1 
