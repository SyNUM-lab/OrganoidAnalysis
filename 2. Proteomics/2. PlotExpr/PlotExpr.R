
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)

# Set working directory
setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/2. Proteomics/2. PlotExpr")

# Load data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/2. Proteomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"pxData_imp.RData"))
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load(paste0(preprocessing_dir,"DEresults_px.RData"))
load(paste0(preprocessing_dir,"proteinAnnotation.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/sampleInfo.RData")

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

#******************************************************************************
#******************************************************************************

annotations[which(annotations$hgnc_symbol == "LY6H"),]
selectedProtein <- "O94772-2"

#******************************************************************************
#******************************************************************************

# Chromosome
annotations$chromosome_name[annotations$uniprot_gn_id %in% str_remove(selectedProtein, "-.*")][1]

# Get HGNC symbol of selected protein
nameProtein <- annotations$hgnc_symbol[annotations$uniprot_gn_id %in% str_remove(selectedProtein, "-.*")][1]

# Select expression data of selected protein
dataMatrix_sel <- data.frame(Value = pxMatrix_imp[str_detect(rownames(pxMatrix_imp),selectedProtein),],
                             Protein = selectedProtein,
                             SampleID = pxData_imp@colData@listData$label)

# Prepare data for plotting
plotDF <- inner_join(sampleInfo, dataMatrix_sel, by = c("SampleID" = "SampleID"))
plotDF$Tissue[plotDF$Tissue == "Cell"] <- "iPSC"
plotDF$TimeRegion <- paste0(plotDF$Tissue, "_",plotDF$Time)
plotDF$TimeRegion <- factor(plotDF$TimeRegion,levels = c("Dorsal_D75",
                                                   "Dorsal_D40",
                                                   "Dorsal_D13",
                                                   "iPSC_D0",
                                                   "Ventral_D13",
                                                   "Ventral_D40",
                                                   "Ventral_D75"))

# Main plot: expression over time
mainPlot <- ggplot(plotDF) +
  geom_point(aes(x = TimeRegion, y = Value, color = Group), size = 2) +
  geom_line(aes(x = TimeRegion, y = Value, color = Group, group = paste0(Replicate, Group))) +
  ylab(expression(log[2]~"intensity")) +
  xlab(NULL) +
  ggtitle(paste0(nameProtein, " (", selectedProtein, ")")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) +
  scale_color_manual(values = c("#377EB8","#E41A1C"))

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
ggsave(finalPlot, file = "D:/RTTproject/CellAnalysis/Proteins/RTTvsIC/Plots/LY6H.png", width = 7, height = 4)



