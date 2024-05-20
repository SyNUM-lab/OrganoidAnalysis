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
Glutamatesig <- data.frame(
  GeneName = c("SLC38A1", "SLC38A2","SLC38A3",
               "GLS", "GLS2",
               "SLC17A7", "SLC17A6", "SLC1A6", "SLC1A1",
               "GRIN1", "GRIN2A", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B",
               "GRIA1", "GRIA2", "GRIA4",
               "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
               "GRM1", "GRM2", "GRM4", "GRM5", "GRM7", "GRM8"),
  ProcessName = c(rep("Glutamine transport",3),
                  rep("Glutamate generation",2), 
                  rep("Glutamate transport", 4), 
                  rep("NMDA receptor", 6),
                  rep("AMPA receptor", 3),
                  rep("Kainate receptor",5),
                  rep("Metabotropic receptor",6)
  )
  
)

# Combine with gene annotation
Glutamatesig <- left_join(Glutamatesig, unique(geneAnnotation[,1:3]),
                     by = c("GeneName" = "GeneName"))

# Re-format data for plotting
GABAsig_D0 <- Glutamatesig
GABAsig_D0$Value <- -log10(pvalues[Glutamatesig$EnsemblID,1]) * sign(logFCs[Glutamatesig$EnsemblID,1])
GABAsig_D0$Sig <- ifelse(adj_pvalues[Glutamatesig$EnsemblID,1]<0.05,"Yes","No")
GABAsig_D0$Time <- 1

GABAsig_D13 <- Glutamatesig
GABAsig_D13$Value <- -log10(pvalues[Glutamatesig$EnsemblID,2]) * sign(logFCs[Glutamatesig$EnsemblID,2])
GABAsig_D13$Sig <- ifelse(adj_pvalues[Glutamatesig$EnsemblID,2]<0.05,"Yes","No")
GABAsig_D13$Time <- 2

GABAsig_D40 <- Glutamatesig
GABAsig_D40$Value <- -log10(pvalues[Glutamatesig$EnsemblID,3]) * sign(logFCs[Glutamatesig$EnsemblID,3])
GABAsig_D40$Sig <- ifelse(adj_pvalues[Glutamatesig$EnsemblID,3]<0.05,"Yes","No")
GABAsig_D40$Time <- 3

GABAsig_D75 <- Glutamatesig
GABAsig_D75$Value <- -log10(pvalues[Glutamatesig$EnsemblID,4]) * sign(logFCs[Glutamatesig$EnsemblID,4])
GABAsig_D75$Sig <- ifelse(adj_pvalues[Glutamatesig$EnsemblID,4]<0.05,"Yes","No")
GABAsig_D75$Time <- 4

GABAsig_all <- rbind.data.frame(
  GABAsig_D0,
  GABAsig_D13,
  GABAsig_D40,
  GABAsig_D75
)

GABAsig_all$GeneName <- factor(GABAsig_all$GeneName,
                               levels = c("SLC38A1", "SLC38A2","SLC38A3",
                                          "GLS", "GLS2",
                                          "SLC17A7", "SLC17A6", "SLC1A6", "SLC1A1",
                                          "GRIN1", "GRIN2A", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B",
                                          "GRIA1", "GRIA2", "GRIA4",
                                          "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
                                          "GRM1", "GRM2", "GRM4", "GRM5", "GRM7", "GRM8"))

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
  #scale_fill_manual(values = c("#A6761D", "#E7298A", "#66A61E", "#E6AB02")) +
  
  ggtitle("Glutamate receptor signaling (dorsal region)") +
  ylim(c(-3,4.5))+
  annotate(x="",y=1:4,label=c("D0", "D13", "D40", "D75"),size=4,geom="text") +
  coord_polar(start=0.1, clip = "off") + 
  theme_void() + 
  theme(axis.text.x = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# Save plot
ggsave(p, file = "GABAmarkers/Glutamatesignaling.png", width = 10, height = 8)


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

ggsave(legend, file = "GABAmarkers/Glutamatesignaling_legend.png", width = 8, height = 8)

