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

# Upregulated proteins
sigMatrix_pos <- matrix((logFCs > 1) & (adjpvalues < 0.05), ncol = ncol(logFCs))

# Downregulated proteins
sigMatrix_neg <- matrix((logFCs < -1) & (adjpvalues < 0.05), ncol = ncol(logFCs))

# Early-time markers
early <- c(rownames(adjpvalues)[(rowSums(sigMatrix_pos[,c(1,2,5)]) == 3) &
                                  (rowSums(sigMatrix_pos[,c(3,4,6,7)]) ==0)],
           rownames(adjpvalues)[(rowSums(sigMatrix_neg[,c(1,2,5)]) == 3) &
                                  (rowSums(sigMatrix_neg[,c(3,4,6,7)]) ==0)]
)

# Mid-time markers
mid <- c(rownames(adjpvalues)[(rowSums(sigMatrix_pos[,c(2,3,5,6)]) >= 3) &
                                  (rowSums(sigMatrix_pos[,c(1,4,7)]) ==0)],
           rownames(adjpvalues)[(rowSums(sigMatrix_neg[,c(2,3,5,6)]) >= 3) &
                                  (rowSums(sigMatrix_neg[,c(1,4,7)]) ==0)]
)

# Late-time markers
late <- c(rownames(adjpvalues)[(rowSums(sigMatrix_pos[,c(3,4,6,7)]) >= 3) &
                                (rowSums(sigMatrix_pos[,c(1,2,5)]) ==0)],
         rownames(adjpvalues)[(rowSums(sigMatrix_neg[,c(3,4,6,7)]) >= 3) &
                                (rowSums(sigMatrix_neg[,c(1,2,5)]) ==0)]
)

# Ventral-region markers
ventral <- c(rownames(adjpvalues)[(rowSums(sigMatrix_pos[,c(1,5,6,7)]) >= 3) &
                                 (rowSums(sigMatrix_pos[,c(2,3,4)]) ==0)],
          rownames(adjpvalues)[(rowSums(sigMatrix_neg[,c(1,5,6,7)]) >= 3) &
                                 (rowSums(sigMatrix_neg[,c(2,3,4)]) ==0)]
)

# Dorsal-region markers
dorsal <- c(rownames(adjpvalues)[(rowSums(sigMatrix_pos[,c(1,2,3,4)]) >= 3) &
                                    (rowSums(sigMatrix_pos[,c(5,6,7)]) ==0)],
             rownames(adjpvalues)[(rowSums(sigMatrix_neg[,c(1,2,3,4)]) >= 3) &
                                    (rowSums(sigMatrix_neg[,c(5,6,7)]) ==0)]
)

# Overall markers
overall <- c(rownames(adjpvalues)[rowSums((logFCs > 1) & (adjpvalues < 0.05)) >= 6],
             rownames(adjpvalues)[rowSums((logFCs < -1) & (adjpvalues < 0.05)) >= 6])



#*****************************************************************************#
#   Make plot
#*****************************************************************************#

# Selection of markers in plot
selection <- data.frame(ID = c(mid, late, dorsal, ventral, overall),
                        ProteinGroup = c(rep("Mid", length(mid)),
                                         rep("Late", length(late)),
                                         rep("Early", length(early)),
                                         rep("Dorsal", length(dorsal)),
                                         rep("Ventral", length(ventral)),
                                         rep("Overall", length(overall)))
                        )
selection$ProteinGroup <- factor(selection$ProteinGroup,
                                 levels = c("Dorsal","Ventral","Early","Mid","Late", "Overall"))

# Combine with protein annotation
selection <- left_join(selection, unique(annotations[,c(1,2,5)]), by = c("ID" = "ID"))

# Save the selection of protein markers
save(selection, file = "proteinSelection.RData")

# Get P-values
pvalues_sel <- adjpvalues[selection$ID,]
plotP <- gather(pvalues_sel, value = "FDR")
plotP$ProteinID <- rep(rownames(pvalues_sel), ncol(pvalues_sel))
plotP$key <- str_remove(plotP$key, "_p.val")
plotP$ID <- paste0(plotP$key, "_",plotP$ProteinID)

# Get logFCs
logfc_sel <- logFCs[selection$ID,]
plotFC <- gather(logfc_sel, value = "log2FC")
plotFC$ProteinID <- rep(rownames(logfc_sel), ncol(logfc_sel))
plotFC$key <- str_remove(plotFC$key, "_ratio")
plotFC$ID <- paste0(plotFC$key, "_",plotFC$ProteinID)

# combine data
plotData <- inner_join(plotFC, plotP[,c(2,4)], by = c("ID" = "ID"))

# Add a column that indicates statistical signficance (FDR < 0.05)
plotData$Sig <- ifelse((plotData$FDR < 0.05) & (abs(plotData$log2FC) > 1), "Yes", "No")

# Combine with feature info
plotData <- inner_join(plotData, selection, by = c("ProteinID" = "ID"))
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

# The order of proteins is determined by clustering
model <- hclust(dist(logfc_sel), "ward.D2")
proteins_ordered <- model$labels[model$order]
plotData$ProteinID <- factor(plotData$ProteinID, 
                                levels = proteins_ordered)
rownames(selection) <- selection$ID
names_ordered <- selection[proteins_ordered, "hgnc_symbol"]
plotData$hgnc_symbol <- factor(plotData$hgnc_symbol, levels = names_ordered)

plotData$CompleteName <- factor(paste0(plotData$hgnc_symbol, " (", plotData$ProteinID, ")"),
                                levels = paste0(names_ordered, " (", proteins_ordered, ")"))
# Main plot
mainPlot <- ggplot() +
  geom_bar(data = plotData, aes(x = abs(log2FC), y = CompleteName, fill = log2FC,  color = Sig),
           linewidth = 0.5, stat = "identity") +
  facet_grid(cols = vars(Group), rows = vars(ProteinGroup),scales = "free", space = "free_y") +
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
  facet_grid(.~sample, scales = "free", space = "free") +
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
  facet_grid(.~sample, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")

rowSideColorPlot <- ggplot(plotData[plotData$Time == "D0",]) +
  geom_tile(aes(x = 1, y = ProteinID), fill = "lightgrey") +
  facet_grid(rows = vars(ProteinGroup), scales = "free", space = "free") +
  theme_void() +
  theme(strip.text = element_blank())


# Combine plots
finalPlot <-  
  colSideColorPlot_tissue + plot_spacer() +
  mainPlot + rowSideColorPlot +
  colSideColorPlot_time + plot_spacer() +
  plot_layout(ncol = 2, nrow = 3, 
              heights = c(0.5,9,0.5),
              widths = c(9,0.1))

# Save plot
ggsave(finalPlot, file = "RTTvsIC_ProteinTargets.png", width = 7.5, height = 10)

# Make legend
legendPlot <- ggplot(data = plotData, aes(x = Group, y = CompleteName, fill = log2FC,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(cols = vars(Region), rows = vars(ProteinGroup), scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-6,6),trans = "pseudo_log", oob = scales::squish) + 
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "") +
  guides(color  ="none") +
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

# Save legend
ggsave(legend, file = "legend_RTTvsIC_ProteinTargets.png", width = 8, height = 8)


#*****************************************************************************#
#   Check targets in transcriptomics data
#*****************************************************************************#

# Load data
load("proteinSelection.RData")
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
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

# Get genes that corresponding to selected proteins
load(paste0(homeDir,"/2. Proteomics/1. Preprocessing/px_gx_annotations1.RData"))
px_gx_annotations <- unique(px_gx_annotations[,c("gene_id", "ID")])
selection_gx <- inner_join(selection, px_gx_annotations, by = c("ID" = "ID"))

# Format P values
pvalues_sel <- data.frame(adj_pvalues_gx[selection_gx$gene_id[!is.na(selection_gx$gene_id)],])
plotP <- gather(pvalues_sel, value = "pvalue")
plotP$GeneID <- rep(rownames(pvalues_sel), ncol(pvalues_sel))
plotP$key <- str_remove(plotP$key, "_p.val")
plotP$ID <- paste0(plotP$key, "_",plotP$GeneID)

# Format logFCs
logfc_sel <- data.frame(logFCs_gx[selection_gx$gene_id[!is.na(selection_gx$gene_id)],])
plotFC <- gather(logfc_sel, value = "log2FC")
plotFC$GeneID <- rep(rownames(logfc_sel), ncol(logfc_sel))
plotFC$key <- str_remove(plotFC$key, "_ratio")
plotFC$ID <- paste0(plotFC$key, "_",plotFC$GeneID)

# combine data
plotData <- inner_join(plotFC, plotP[,c(2,4)], by = c("ID" = "ID"))

# Add a column that indicates statistical signficance (FDR < 0.05)
plotData$Sig <- ifelse((plotData$pvalue < 0.05) & (abs(plotData$log2FC) >1), "Yes", "No")

# Combine with feature info
plotData <- inner_join(plotData, selection_gx, by = c("GeneID" = "gene_id"))

# Format data for plotting
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

# Set column order
plotData$ProteinID <- factor(plotData$ID.y, 
                             levels = proteins_ordered)
rownames(selection_gx) <- selection_gx$ID
names_ordered <- selection_gx[proteins_ordered, "hgnc_symbol"]

# Set row order
plotData$hgnc_symbol <- factor(plotData$hgnc_symbol, levels = names_ordered)
plotData$CompleteName <- factor(paste0(plotData$hgnc_symbol, " (", plotData$ProteinID, ")"),
                                levels = paste0(names_ordered, " (", proteins_ordered, ")"))

# Main plot
mainPlot <- ggplot() +
  geom_bar(data = plotData, aes(x = abs(log2FC), y = CompleteName, fill = log2FC,  color = Sig),
           linewidth = 0.5, stat = "identity") +
  facet_grid(cols = vars(Group), rows = vars(ProteinGroup),scales = "free", space = "free_y") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-6,6),trans = "pseudo_log", oob = scales::squish) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "logFC") +
  ylab("Protein Name\n") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, hjust = 1),
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
  facet_grid(.~sample, scales = "free", space = "free") +
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
  facet_grid(.~sample, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")

rowSideColorPlot <- ggplot(plotData[plotData$Time == "D0",]) +
  geom_tile(aes(x = 1, y = ProteinID), fill = "lightgrey") +
  facet_grid(rows = vars(ProteinGroup), scales = "free", space = "free") +
  theme_void() +
  theme(strip.text = element_blank())


# Combine plots
finalPlot <-  
  colSideColorPlot_tissue + plot_spacer() +
  mainPlot + rowSideColorPlot +
  colSideColorPlot_time + plot_spacer() +
  plot_layout(ncol = 2, nrow = 3, 
              heights = c(0.5,9,0.5),
              widths = c(9,0.1))

# Save plot
ggsave(finalPlot, file = "RTTvsIC_GeneTargets.png", width = 7.5, height = 10)

# Make legend
legendPlot <- 
  ggplot(data = plotData, aes(x = Group, y = CompleteName, fill = log2FC,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(cols = vars(Region), rows = vars(ProteinGroup), scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-6,6), trans = "pseudo_log", oob = scales::squish) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "") +
  guides(color = "none")+
  ylab("Protein Name\n") +
  scale_y_discrete(position = "right") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y.right = element_text(size = 10, hjust = 0),
        #axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))

legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "legend_RTTvsIC_GeneTargets.png", width = 8, height = 8)

