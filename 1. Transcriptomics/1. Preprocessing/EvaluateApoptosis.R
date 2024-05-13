# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/"))

# Load data
load("gxMatrix_norm.RData")
load("geneAnnotation.RData")
load(paste0(homeDir,"/SampleInfo.RData"))
load("DEresults_Time_gx.RData")

# Load packages
library(ggpubr)
library(tidyverse)

#******************************************************************************#
# Get FPKM values
#******************************************************************************#

# Convert log CPM to CPM
gxMatrix_count <- 2^(gxMatrix_norm)-1
  
# Filter annotation file
geneAnnotation1 <- unique(geneAnnotation[,c(1,6)])
rownames(geneAnnotation1) <- geneAnnotation1$gene_id
  
# Calculate FPKM
RPM <- t(gxMatrix_count)
FPKM <- matrix(NA, nrow = nrow(gxMatrix_count), ncol = ncol(gxMatrix_count))
for (i in 1:nrow(gxMatrix_count)){
  FPKM[i,] <- t(RPM[,rownames(gxMatrix_count)[i]])/(as.numeric(geneAnnotation1[geneAnnotation1$gene_id %in% rownames(gxMatrix_count)[i],"effective_length"])/1000)
}
FPKM <- log2(FPKM+1)
colnames(FPKM) <- colnames(gxMatrix_count)
rownames(FPKM) <- rownames(gxMatrix_count)

# Save FPKM values
save(FPKM, file = "FPKM.RData")


#******************************************************************************#
# Make plot with expression values of apoptotic genes
#******************************************************************************#

# Load FPKM values
load("FPKM.RData")

# Exclude RTT samples
FPKM_IC <- FPKM[,sampleInfo$SampleID[sampleInfo$Group == "IC"]]

# Select expression of apoptotic genes:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7169302/#:~:text=Genetic%20regulators%20(mostly%20pro%E2%80%90apoptotic,DRs%20and%20the%20caspase%20family.

apoptosis_genes <- geneAnnotation[geneAnnotation$GeneName %in% 
                                    c("BCL2","TP53", "BAX", "CASP3",
                                      "FAS"),]

# Get statistics
topList$Dorsal_D75_vs_Dorsal_D13[rownames(topList$Dorsal_D75_vs_Dorsal_D13) %in%
                                   apoptosis_genes$gene_id,]

topList$Ventral_D75_vs_Ventral_D13[rownames(topList$Ventral_D75_vs_Ventral_D13) %in%
                                   apoptosis_genes$gene_id,]


# Gather data
plotDF <- gather(data.frame(FPKM_IC))
plotDF$Gene <- rep(rownames(FPKM_IC), ncol(FPKM_IC))

# Prepare data for plotting
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
plotDF$Sample_Gene <- paste0(plotDF$ReplicateGroup, "_", plotDF$Gene)
plotDF <- plotDF %>%
  group_by(Sample_Gene) %>%
  mutate(MeanExpr = mean(value))

# Add gene annotation
plotDF <- inner_join(plotDF, geneAnnotation, by = c("Gene" = "gene_id"))

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
  geom_point(data = plotDF[plotDF$Gene %in% apoptosis_genes$gene_id,],
                           aes(x = TimeRegion, y = MeanExpr, color = GeneName), size = 3) +
  geom_line(data = plotDF[plotDF$Gene %in% apoptosis_genes$gene_id,],
            aes(x = TimeRegion, y = MeanExpr, group = Gene, color = GeneName), linewidth = 1) +
  ylab(expression(log[2]~"FPKM")) +
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
ggsave(finalPlot, file = "QCplots/ApoptoticGenes.png", width = 7, height = 5)
