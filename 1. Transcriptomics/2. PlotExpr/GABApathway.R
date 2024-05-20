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
library(ggnewscale)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/2. PlotExpr/"))

# Load data
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(homeDir,"/SampleInfo.RData"))
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


# Make data frame for plotting
GABAsig <- data.frame(
  GeneName = c("GAD1", "GAD2", 
               "SLC32A1", "SLC6A1", "SLC6A11", 
               "GABRA1", "GABRA2", "GABRA4", "GABRA5", "GABRD", "GABRP",
               "GABRB1", "GABRB2", "GABRB3", "GABRG1", "GABRG2", "GABRG3",
               "GABBR1", "GABBR2"),
  ProcessName = c(rep("Decarboxylase",2), 
                  rep("Transporter", 3),
                  rep("GABA-A receptor", 12),
                  rep("GABA-B receptor", 2)
                  
  )
  
)

# Combine with gene annotation
GABAsig <- left_join(GABAsig, unique(geneAnnotation[,1:3]),
                     by = c("GeneName" = "GeneName"))

# Re-format data for plotting
GABAsig_D0 <- GABAsig
GABAsig_D0$Value <- -log10(pvalues[GABAsig$EnsemblID,1]) * sign(logFCs[GABAsig$EnsemblID,1])
GABAsig_D0$Sig <- ifelse(adj_pvalues[GABAsig$EnsemblID,1]<0.05,"Yes","No")
GABAsig_D0$Time <- 1

GABAsig_D13 <- GABAsig
GABAsig_D13$Value <- -log10(pvalues[GABAsig$EnsemblID,5]) * sign(logFCs[GABAsig$EnsemblID,5])
GABAsig_D13$Sig <- ifelse(adj_pvalues[GABAsig$EnsemblID,5]<0.05,"Yes","No")
GABAsig_D13$Time <- 2

GABAsig_D40 <- GABAsig
GABAsig_D40$Value <- -log10(pvalues[GABAsig$EnsemblID,6]) * sign(logFCs[GABAsig$EnsemblID,6])
GABAsig_D40$Sig <- ifelse(adj_pvalues[GABAsig$EnsemblID,6]<0.05,"Yes","No")
GABAsig_D40$Time <- 3

GABAsig_D75 <- GABAsig
GABAsig_D75$Value <- -log10(pvalues[GABAsig$EnsemblID,7]) * sign(logFCs[GABAsig$EnsemblID,7])
GABAsig_D75$Sig <- ifelse(adj_pvalues[GABAsig$EnsemblID,7]<0.05,"Yes","No")
GABAsig_D75$Time <- 4

GABAsig_all <- rbind.data.frame(
  GABAsig_D0,
  GABAsig_D13,
  GABAsig_D40,
  GABAsig_D75
)

GABAsig_all$GeneName <- factor(GABAsig_all$GeneName,
                               levels = c("GAD1", "GAD2", 
                                          "SLC32A1", "SLC6A1", "SLC6A11", 
                                          "GABRA1", "GABRA2", "GABRA4", "GABRA5", "GABRD", "GABRP",
                                          "GABRB1", "GABRB2", "GABRB3", "GABRG1", "GABRG2", "GABRG3",
                                          "GABBR1", "GABBR2"))


# Make plot
p <- ggplot(GABAsig_all) +
  geom_tile(aes(x = GeneName, y = Time, fill = Value, color = Sig), 
            width = 0.9, height = 0.9, size = 1) +
  #facet_grid(cols = vars(ProcessName),  scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red",
                       midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  
  new_scale_fill() +
  geom_tile(aes(x = GeneName, y = -0.5, fill = ProcessName), alpha = 0.2) +
  scale_fill_manual(values = c("#A6761D", "#E7298A", "#66A61E", "#E6AB02")) +
  
  ggtitle("GABA receptor signaling (ventral region)") +
  ylim(c(-3,4.5))+
  annotate(x="",y=1:4,label=c("D0", "D13", "D40", "D75"),size=4,geom="text") +
  coord_polar(start=0.15, clip = "off") + 
  theme_void() + 
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# Save plot
ggsave(p, file = "GABAmarkers/GABAsignaling.png", width = 10, height = 8)


# Get legends:
legendPlot <- ggplot(GABAsig_all) +
  geom_tile(aes(x = GeneName, y = Time, fill = Value), 
            width = 0.9, height = 0.9) +
  #facet_grid(cols = vars(ProcessName),  scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red",
                       midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  
  ggtitle("GABA receptor signaling (ventral region)") +
  ylim(c(-3,4.5))+
  annotate(x="",y=1:4,label=c("D0", "D13", "D40", "D75"),size=4,geom="text") +
  coord_polar(start=0.15, clip = "off") + 
  theme_void() + 
  labs(fill = "") +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))


legend <- cowplot::get_legend(legendPlot)

# Save plot
ggsave(legend, file = "GABAmarkers/GABAsignaling_legend.png", width = 8, height = 8)

