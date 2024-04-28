#*****************************************************************************#
#   Preparation
#*****************************************************************************#

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/3. TimeAnalysis")

# Load packages
library(tidyverse)
library(ggrepel)
library(patchwork)

# Load data:
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/SampleInfo.RData")
all(sampleInfo$SampleID == colnames(gxMatrix_norm))


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

# Prepare data for plotting
plotData <- gather(data.frame(gxMatrix_fil))
plotData$GeneID <- rep(rownames(gxMatrix_fil), ncol(gxMatrix_fil))

# Add gene annotation to data frame
plotData <- inner_join(plotData, geneAnnotation, by = c("GeneID" = "gene_id"))

# Add sample information to data frame
plotData <- inner_join(plotData, sampleInfo, by = c("key" = "SampleID"))

# Set order of columns
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

# Set order of rows
plotData$GeneName <- factor(plotData$GeneName, 
                            levels = c("NANOG",
                                       "NES", 
                                       "FOXG1",
                                       "LHX6",
                                       "NKX2-1",
                                       "DLX2",
                                       "PAX6",
                                       "TBR1"))

# Make plot
main <- ggplot(plotData) +
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

# Save plot
ggsave(p, file = "Heatmap_QCmarkers.png", width = 5, height = 7)


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

# Save plot
ggsave(legend, file = "legend_QCmarkers.png", width = 8, height = 8)




