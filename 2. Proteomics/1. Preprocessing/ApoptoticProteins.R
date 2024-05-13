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

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/2. Proteomics/1. Preprocessing"))

# Load sample information
load(paste0(homeDir,"/sampleInfo.RData"))

# Load data
load("pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load("proteinAnnotation.RData")
all(annotations$uniprot_gn_id %in% str_remove(rownames(pxMatrix_imp), "-.*"))

#******************************************************************************#
# Make plot
#******************************************************************************#

# Exclude RTT samples
pxMatrix_IC <- pxMatrix_imp[,!str_detect(colnames(pxMatrix_imp), "RTT")]

# Get apoptosis genes
apoptosis_genes <- annotations[annotations$hgnc_symbol %in% 
                                    c("BAX", "CASP3", "CASP6",
                                      "BAG1"),]

# Gather data
plotDF <- gather(data.frame(pxMatrix_IC))
plotDF$Protein <- rep(rownames(pxMatrix_IC, ncol(pxMatrix_IC)))

# Prepare data for plotting
sampleInfo$SampleID <- paste0(sampleInfo$Group,"_",
                              sampleInfo$Tissue, "_",
                              sampleInfo$Time, "_",
                              sampleInfo$Replicate
                              )
plotDF <- inner_join(sampleInfo, plotDF, by = c("SampleID" = "key"))
plotDF$Tissue[plotDF$Tissue == "Cell"] <- "iPSC"
plotDF$TimeRegion <- paste0(plotDF$Tissue, "_",plotDF$Time)
plotDF$TimeRegion <- factor(plotDF$TimeRegion,levels = c("Dorsal_D75",
                                                         "Dorsal_D40",
                                                         "Dorsal_D13",
                                                         "iPSC_D0",
                                                         "Ventral_D13",
                                                         "Ventral_D40",
                                                         "Ventral_D75"))
plotDF$ID <-  paste0(plotDF$TimeRegion, "_", plotDF$Gene)

# Get mean expression per gene and time point/region
plotDF$Sample_Protein <- paste0(plotDF$ReplicateGroup, "_", plotDF$Protein)
plotDF <- plotDF %>%
  group_by(Sample_Protein) %>%
  mutate(MeanExpr = mean(value))

# Add gene annotation
plotDF <- inner_join(plotDF, unique(annotations[,c(1,2,5)]), by = c("Protein" = "ID"))

# Set colors
colors <- c("#1B9E77",
            "#D95F02",
            "#7570B3",
            "#E7298A",
            "#66A61E",
            "#E6AB02",
            "#A6761D",
            "#666666",
            "#CB181D")


# Make main plot
mainPlot <- ggplot() +
  geom_boxplot(data = plotDF,
               aes(x = TimeRegion, y = MeanExpr), color = "#525252", fill = "#F0F0F0") +
  geom_point(data = plotDF[plotDF$Protein %in% apoptosis_genes$ID,],
             aes(x = TimeRegion, y = MeanExpr, color = hgnc_symbol), size = 3) +
  geom_line(data = plotDF[plotDF$Protein %in% apoptosis_genes$ID,],
            aes(x = TimeRegion, y = MeanExpr, group = Protein, color = hgnc_symbol), linewidth = 1) +
  ylab(expression(log[2]~"intensity")) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) +
  scale_color_manual(values = colors)



# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotDF$TimeRegion,
  time = plotDF$Time,
  region = plotDF$Tissue))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
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
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots into a single figure
finalPlot <- ggarrange(mainPlot,
                       colSideColorPlot_time,
                       colSideColorPlot_tissue,
                       heights = c(6,0.5,0.5),nrow = 3,ncol = 1,
                       common.legend = TRUE,
                       legend = "right",
                       align = "v")
# Print plot
finalPlot

# Save plot
ggsave(finalPlot, file = paste0("QCplots/ApoptoticProteins.png"), width = 7, height = 5)
