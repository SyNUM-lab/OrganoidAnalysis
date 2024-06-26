
#*****************************************************************************#
#   Preparation
#*****************************************************************************#

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/2. Proteomics/1. Preprocessing"))

# Load packages
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(rrvgo)
library(ggrepel)
library(rrvgo)
library(BaseSet)

# FUNCTION: Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Load data:
load("pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
colnames(pxMatrix_imp) <- pxData_imp@colData$label
load("proteinAnnotation.RData")
load(paste0(homeDir,"/SampleInfo.RData"))


#*****************************************************************************#
#   PCA
#*****************************************************************************#

# Run PCA
pca <- prcomp(t(pxMatrix_imp), 
              retx = TRUE, # Give rotated variables (PCA loadings)
              center = TRUE,
              scale = TRUE) # Variables are scaled to unit variance

# Get expained variances
expl_var <- round((pca$sdev^2)/sum(pca$sdev^2),3)


#*****************************************************************************#
#   Score plot
#*****************************************************************************#

# Combine sample info in PCA scores
plotPCA <- as.data.frame(pca$x)

# Add sample information to dataframe
if (all(rownames(plotPCA) == sampleInfo$SampleID)){
  plotPCA <- cbind(plotPCA, sampleInfo)
}

# Get loadings
loadingsPCA <- data.frame(pca$rotation)
loadingsPCA$ID <- rownames(loadingsPCA)
loadingsPCA <- inner_join(loadingsPCA, unique(annotations[,c(1,2,5)]),
                          by = c("ID" = "ID"))

# Select proteins with maximum/minimum loading
selectedProteins <- loadingsPCA[loadingsPCA$ID %in% c(names(head(sort(pca$rotation[,1]),5)),
                                                      names(tail(sort(pca$rotation[,1]),5))),]
selectedProteins$hgnc_symbol <- factor(selectedProteins$hgnc_symbol,
                                       levels = arrange(selectedProteins, PC1)[,"hgnc_symbol"])


# Make plot
plotPCA$Colour <- paste0(plotPCA$Group, ": ", plotPCA$Time)
plotPCA$Tissue[plotPCA$Tissue == "Cell"] <- "iPSC"
plotPCA$Tissue <- factor(plotPCA$Tissue, levels = c("iPSC", "Dorsal", "Ventral"))

colors <- c("#6BAED6","#4292C6","#2171B5","#084594",
            "#FB6A4A","#EF3B2C","#CB181D","#99000D")

colors <- c("#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#4A1486",
            "#FDAE6B","#FD8D3C","#F16913","#D94801","#8C2D04")
fills <- c("#F0F0F0","#D9D9D9","#BDBDBD","#969696")


p <- ggplot() +
  geom_point(data = plotPCA, 
             aes(x = PC1, y = PC2, fill = Time), 
             color = "grey", alpha = 1, shape = 21, size = 3) +
  geom_segment(data = selectedProteins,
               aes(x = 0, y = 0, xend = PC1*1500, yend = PC2*1500, color = hgnc_symbol),
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.5) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = fills) +
  xlab(paste0("PC1 (", expl_var[1]*100,"%)")) +
  ylab(paste0("PC2 (", expl_var[2]*100,"%)")) +
  ggtitle("Proteomics") +
  labs(fill = "Time", color = "Protein") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 18,
                                  face = "bold"),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.1,
                                        linetype = 1),
        panel.grid.minor= element_line(color = "lightgrey",
                                       size = 0.05,
                                       linetype = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#theme(legend.title = element_text(face = "bold"))

ggsave(p, file = "TimeAnalysis_loadings.png", width = 7, height = 5)





