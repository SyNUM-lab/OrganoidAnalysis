# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)

# Set working directory
setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/2. PlotExpr")

# Load data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/SampleInfo.RData")

# Get p-values
genes_all <- rownames(topList[[1]])
pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues) <- genes_all
colnames(pvalues) <- names(topList)

# Gather p-values
p_DF <- gather(data.frame(pvalues))
p_DF$Gene <- rep(rownames(pvalues), ncol(pvalues))
p_DF$Sig <- ifelse(p_DF$value < 0.05,"*","")
p_DF$Sig[p_DF$value < 0.01] <- "**"
p_DF$Sig[p_DF$value < 0.001] <- "***"

# Add region to data frame
p_DF$Region <- NA
p_DF$Region[str_detect(p_DF$key, "Cell")] <- "iPSC"
p_DF$Region[str_detect(p_DF$key, "Dorsal")] <- "Dorsal"
p_DF$Region[str_detect(p_DF$key, "Ventral")] <- "Ventral"

# Add time to data frame
p_DF$Time <- "D0"
p_DF$Time[str_detect(p_DF$key, "D13")] <- "D13"
p_DF$Time[str_detect(p_DF$key, "D40")] <- "D40"
p_DF$Time[str_detect(p_DF$key, "D75")] <- "D75"

# Make IDs
p_DF$TimeRegion <- paste0(p_DF$Region,"_", p_DF$Time)
p_DF$ID <- paste0(p_DF$TimeRegion, "_", p_DF$Gene)

#==============================================================================#
# select genes
test <- geneAnnotation[which(str_detect(geneAnnotation$GeneName, "PCSK6")),]
selectedGene <- "ENSG00000140479_PCSK6"
#==============================================================================#

# Get HGNC symbol of selected gene
nameGene <- geneAnnotation$GeneName[geneAnnotation$gene_id == selectedGene][1]

# Select expression data of selected protein
dataMatrix_sel <- data.frame(Expr = gxMatrix_norm[selectedGene,],
                             Gene = selectedGene,
                             SampleID = colnames(gxMatrix_norm))

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
plotDF$ID <-  paste0(plotDF$TimeRegion, "_", plotDF$Gene)
plotDF <- inner_join(plotDF, p_DF[,c("ID", "value", "Sig")], by = c("ID" = "ID"))

plotDF1 <- plotDF %>%
  group_by(ID) %>%
  mutate(testValue = max(Expr) + 0.3)


# Main plot: expression over time
mainPlot <- ggplot(plotDF1) +
  geom_point(aes(x = TimeRegion, y = Expr, color = Group), size = 2) +
  geom_line(aes(x = TimeRegion, y = Expr, color = Group, group = paste0(Replicate, Group))) +
  geom_text(aes(x = TimeRegion, y = testValue, label = Sig), color = "#595959", size = 6) +
  ylab(expression(log[2]~"CPM")) +
  xlab(NULL) +
  #ggtitle(paste0(nameGene, " (", selectedGene, ")")) +
  ggtitle(nameGene) +
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
ggsave(finalPlot, file = paste0(nameGene,".png"), width = 6, height = 4)
