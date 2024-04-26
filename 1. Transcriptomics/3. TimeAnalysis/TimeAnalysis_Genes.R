
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
library(org.Hs.eg.db)
library(clusterProfiler)
library(rrvgo)
library(ggrepel)
library(rrvgo)
library(BaseSet)

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

# Make plot
plotPCA$Colour <- paste0(plotPCA$Group, ": ", plotPCA$Time)
plotPCA$Tissue[plotPCA$Tissue == "Cell"] <- "iPSC"
plotPCA$Tissue <- factor(plotPCA$Tissue, levels = c("iPSC", "Dorsal", "Ventral"))

colors <- c("#6BAED6","#4292C6","#2171B5","#084594",
            "#FB6A4A","#EF3B2C","#CB181D","#99000D")
p <- ggplot(plotPCA, aes(x = PC1, y = PC2, color = Colour, shape = Tissue)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", expl_var[1]*100,"%)")) +
  ylab(paste0("PC2 (", expl_var[2]*100,"%)")) +
  ggtitle("Gene Expression Data") +
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

ggsave(filename = "Plots/PCA_Scores_gx.png", p, width = 7, height = 5)



#*****************************************************************************#
#   GSEA
#*****************************************************************************#

# Get loadings
loadings <- as.data.frame(pca$rotation)
loadings$GeneID <- rownames(loadings)

# Loadings of annotated genes
geneList <- inner_join(geneAnnotation, loadings, by = c("gene_id" = "GeneID"))

# GSEA input: ordered vector
gsea_input <- geneList$PC1
names(gsea_input) <- geneList$EnsemblID
gsea_input <- sort(gsea_input, decreasing = TRUE)

# Biological process
set.seed(123)
BPresults <- gseGO(
  geneList = gsea_input,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont =  "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500)

save(BPresults, file = "Data/BP_GSEA1.RData")


#*****************************************************************************#
#   Make plots
#*****************************************************************************#
# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Load data
load("Data/BP_GSEA1.RData")

# Get BP results
BP <- BPresults@result

# Make similarity matrix of GO terms
simMatrix <- calculateSimMatrix(
  BP$ID,
  orgdb="org.Hs.eg.db",
  ont = "BP",
  method = "Resnik"
)
scores <- setNames(-log10(BP$pvalue), BP$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.85,
                                orgdb="org.Hs.eg.db")

save(reducedTerms, file = "Data/reducedTerms_PC1_BP1.RData")

load("Data/reducedTerms_PC1_BP1.RData")

# Reduce number of GO terms
BP_fil <- BP[unique(reducedTerms$parent),]

# Get median PC1 and PC2 for each GO term
averages <- as.data.frame(rep(0, nrow(BP_fil)))
colnames(averages) <- "PC1"
averages$PC2 <- rep(0, nrow(BP_fil))
averages$Description <- BP_fil$Description
averages$ID <- BP_fil$ID
averages$NES <- BP_fil$NES
averages$pvalue <- BP_fil$pvalue
for (i in 1:nrow(BP_fil)){
  genes <- BPresults@geneSets[[BP_fil$ID[i]]]
  geneIDs <- geneList$gene_id[geneList$EnsemblID %in% genes]
  averages$PC1[i] <- median(loadings[geneIDs,"PC1"])
  averages$PC2[i] <- median(loadings[geneIDs, "PC2"])
}


# Select relevant terms
averages <- arrange(averages, by = NES)
selectedTerms <- rbind.data.frame(head(averages,5), tail(averages,5))
selectedTerms$Description <- firstup(selectedTerms$Description)
selectedTerms$Description <- factor(selectedTerms$Description, 
                                    levels = c("Mitochondrial gene expression",
                                               "Transcription by RNA polymerase I",
                                               "Mitochondrial RNA metabolic process",
                                               "RNA modification",
                                               "Chromosome segregation",
                                               "Neuron recognition",
                                               "Regulation of synapse structure or activity",
                                               "Synaptic transmission, glutamatergic",
                                               "Cognition",
                                               "Regulation of neurotransmitter receptor activity"
                                              ))

colors <- c("#9ECAE1","#6BAED6","#4292C6","#2171B5","#084594",
            "#FC9272","#FB6A4A","#EF3B2C","#CB181D","#99000D")
colors <- c("#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#4A1486",
            "#FDAE6B","#FD8D3C","#F16913","#D94801","#8C2D04")
names(colors) <- c("Mitochondrial gene expression",
                   "Transcription by RNA polymerase I",
                   "Mitochondrial RNA metabolic process",
                   "RNA modification",
                   "Chromosome segregation",
                   "Neuron recognition",
                   "Regulation of synapse structure or activity",
                   "Synaptic transmission, glutamatergic",
                   "Cognition",
                   "Regulation of neurotransmitter receptor activity"
)

fills <- c("#F0F0F0","#D9D9D9","#BDBDBD","#969696")

p <- ggplot() +
  geom_point(data = plotPCA, 
             aes(x = PC1, y = PC2, fill = Time), 
             color = "grey", alpha = 1, shape = 21, size = 3) +
  geom_segment(data = selectedTerms,
               aes(x = 0, y = 0, xend = PC1*9000, yend = PC2*9000, color = Description),
               arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.5) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = fills) +
  xlab(paste0("PC1 (", expl_var[1]*100,"%)")) +
  ylab(paste0("PC2 (", expl_var[2]*100,"%)")) +
  ggtitle("Transcriptomics") +
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
                                       linetype = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  #theme(legend.title = element_text(face = "bold"))

ggsave(p, file = "Plots/TimeAnalysis_loadings_alt1.png", width = 9, height = 5)


