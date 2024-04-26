#*****************************************************************************#
#   Preparation
#*****************************************************************************#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
setwd("D:/RTTproject/CellAnalysis/Genes/TimeAnalysis")

# Load packages
library(tidyverse)
library(ggrepel)
library(patchwork)

# Load data:
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_norm.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/geneAnnotation.RData")
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")
all(sampleInfo$SampleID == colnames(gxMatrix_norm))

#*****************************************************************************#
#   PCA
#*****************************************************************************#

# Run PCA
pca <- prcomp(t(gxMatrix_norm), 
              retx = TRUE, # Give rotated variables (PCA loadings)
              center = TRUE,
              scale = TRUE) # Variables are scaled to unit variance

# Combine sample info in PCA scores
plotPCA <- as.data.frame(pca$x)

# Add sample information to dataframe
if (all(rownames(plotPCA) == sampleInfo$SampleID)){
  plotPCA <- cbind(plotPCA, sampleInfo)
}

# Make plot
plotPCA$Colour <- paste0(plotPCA$Group, ": ", plotPCA$Time)
plotPCA$Tissue[plotPCA$Tissue == "Cell"] <- "iPSC"
plotPCA$Tissue <- factor(plotPCA$Tissue, levels = c("iPSC", "Dorsal", "Ventral"))


# Get expained variances
expl_var <- round((pca$sdev^2)/sum(pca$sdev^2),3)
loadings <- as.data.frame(pca$rotation)
loadings$GeneID <- rownames(loadings)
loadings <- inner_join(loadings, geneAnnotation, by = c("GeneID" = "gene_id"))

markers <- c("ENSG00000111704_NANOG", 
             "ENSG00000176165_FOXG1",
             "ENSG00000106852_LHX6",
             "ENSG00000136352_NKX2-1",
             "ENSG00000115844_DLX2",
             "ENSG00000007372_PAX6",
             "ENSG00000136535_TBR1",
             "ENSG00000132688_NES")

loadings_fil <- loadings[loadings$GeneID %in% markers,]


fills <- c("#F0F0F0","#D9D9D9","#BDBDBD","#969696")

p <- ggplot() +
  geom_point(data = plotPCA, 
             aes(x = PC1, y = PC2, fill = Time), 
             color = "grey", alpha = 1, shape = 21, size = 3) +
  geom_segment(data = loadings_fil,
               aes(x = 0, y = 0, xend = PC1*8000, yend = PC2*8000),
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.5) +
  ggrepel::geom_text_repel(data = loadings_fil,
               aes(x = PC1*8000, y = PC2*8000, label = GeneName),
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.5) +
  #scale_color_manual(values = colors) +
  scale_fill_manual(values = fills) +
  xlab(paste0("PC1 (", expl_var[1]*100,"%)")) +
  ylab(paste0("PC2 (", expl_var[2]*100,"%)")) +
  labs(fill = "Time", color = "GO term") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 18,
                                  face = "bold"),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.1,
                                        linetype = 1),
        panel.grid.minor= element_line(color = "lightgrey",
                                       size = 0.05,
                                       linetype = 1))
#theme(legend.title = element_text(face = "bold"))

# Make plot
plotPCA$Colour <- paste0(plotPCA$Group, ": ", plotPCA$Time)
plotPCA$Tissue[plotPCA$Tissue == "Cell"] <- "iPSC"
plotPCA$Tissue <- factor(plotPCA$Tissue, levels = c("iPSC", "Dorsal", "Ventral"))

colors <- c("#6BAED6","#4292C6","#2171B5","#084594",
            "#FB6A4A","#EF3B2C","#CB181D","#99000D")
p <- ggplot() +
  geom_segment(data = loadings_fil,
               aes(x = 0, y = 0, xend = PC1*7000, yend = PC2*7000),
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.5, color = "grey") +
  geom_point(data = plotPCA, 
             aes(x = PC1, y = PC2, color = Colour, shape = Tissue),size = 3, alpha = 1) +
  geom_text(data = loadings_fil,
                           aes(x = PC1*7000, y = PC2*7000, label = GeneName),
                           arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.5) +
  xlab(paste0("PC1 (", expl_var[1]*100,"%)")) +
  ylab(paste0("PC2 (", expl_var[2]*100,"%)")) +
  #ggtitle("Gene Expression Data") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 18,
                                  face = "bold"),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.1,
                                        linetype = 1),
        panel.grid.minor= element_line(color = "lightgrey",
                                       size = 0.05,
                                       linetype = 1)) +
  guides(color = guide_legend("Time"), shape = guide_legend("Region"))

ggsave(filename = "Plots/PCA_Scores_QCmarkers.png", p, width = 7, height = 5)



#*****************************************************************************#
#   Heatmap
#*****************************************************************************#

# Select markers
markers <- c("ENSG00000111704_NANOG", 
             "ENSG00000176165_FOXG1",
             "ENSG00000106852_LHX6",
             "ENSG00000136352_NKX2-1",
             "ENSG00000115844_DLX2",
             "ENSG00000007372_PAX6",
             "ENSG00000136535_TBR1",
             "ENSG00000132688_NES")

gxMatrix_fil <- gxMatrix_norm[markers,]

# Standardize expression
gxMatrix_fil <- (gxMatrix_fil - rowMeans(gxMatrix_fil))/apply(gxMatrix_fil,1,sd)


plotData <- gather(data.frame(gxMatrix_fil))
plotData$GeneID <- rep(rownames(gxMatrix_fil), ncol(gxMatrix_fil))

plotData <- inner_join(plotData, geneAnnotation, by = c("GeneID" = "gene_id"))

plotData <- inner_join(plotData, sampleInfo, by = c("key" = "SampleID"))

plotData$Tissue[plotData$Tissue == "Cell"] <- "iPSC"
plotData$TimeRegion <- paste0(plotData$Time, "_", plotData$Tissue)
plotData$TimeRegion <- factor(plotData$TimeRegion,
                              levels = c("D75_Dorsal", "D40_Dorsal", "D13_Dorsal",
                                         "D0_iPSC",
                                         "D13_Ventral", "D40_Ventral", "D75_Ventral"))
plotData$Sample <- paste0(plotData$Group, "_", plotData$Replicate)


plotData$key <- factor(plotData$key,
                       levels = c(rev(sampleInfo$SampleID[sampleInfo$Tissue == "Dorsal"]),
                                  sampleInfo$SampleID[sampleInfo$Tissue == "Cell"],
                                  sampleInfo$SampleID[sampleInfo$Tissue == "Ventral"]))

plotData$GeneName <- factor(plotData$GeneName, 
                            levels = c("NANOG",
                                       "NES", 
                                       "FOXG1",
                                       "LHX6",
                                       "NKX2-1",
                                       "DLX2",
                                       "PAX6",
                                       "TBR1"))


main <- ggplot(plotData) +
  #geom_rect(xmax = "D75_Dorsal", xmin = 0, ymin = Inf, ymax = 3.5, fill = "#377EB8", color = "white") +
  #geom_rect(xmax = "D75_Dorsal", xmin = 0, ymax = -Inf, ymin = 3.5, fill = "#E41A1C", color = "white") +
  geom_tile(aes(x = TimeRegion, y = Sample, fill = value), color = "white", size = 0.1) +
  facet_grid(cols = vars(Tissue), rows = vars(GeneName), scale = "free", space = "free") +
  scale_fill_viridis_c() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        #strip.background.y = element_rect(color="black", fill="grey", size=0.5, linetype="solid"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 9, color = "black", angle = 270))



#Column side color (Time and tissue/region)
colSideColor<- unique(data.frame(
  sample = factor(plotData$TimeRegion,
                  levels = c("D75_Dorsal", "D40_Dorsal", "D13_Dorsal",
                             "D0_iPSC",
                             "D13_Ventral", "D40_Ventral", "D75_Ventral")),
  time = plotData$Time,
  tissue = plotData$Tissue,
  Test = "1"))


# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
  facet_grid(cols = vars(tissue), rows = vars(Test), scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# Tissue/region
colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[(colSideColor$time == "D40") | (colSideColor$time == "D0"),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(cols = vars(tissue), rows = vars(Test), scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15)) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots
p <- ggpubr::ggarrange(main,
                       colSideColorPlot_time,
                       colSideColorPlot_tissue,
                       nrow = 3,
                       ncol = 1,
                       heights = c(9,0.5,0.5),
                       align = "v",
                       common.legend = FALSE)

ggsave(p, file = "Plots/Heatmap_QCmarkers1.png", width = 5, height = 7)


# Get legends:
legendPlot <- ggplot(plotData) +
  geom_tile(aes(x = TimeRegion, y = Sample, fill = value), color = "white", size = 0.1) +
  facet_grid(cols = vars(Tissue), rows = vars(GeneName), scale = "free", space = "free") +
  scale_fill_viridis_c() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        #strip.background.y = element_rect(color="black", fill="grey", size=0.5, linetype="solid"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 9, color = "black", angle = 270),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))


legend <- cowplot::get_legend(legendPlot)

ggsave(legend, file = "Plots/legend_QCmarkers.png", width = 8, height = 8)




