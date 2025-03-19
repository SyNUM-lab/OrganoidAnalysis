
################################################################################

# sample-wise correlations

################################################################################


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "E:/RTTproject/CellAnalysis/OrganoidAnalysis"
load(paste0(homeDir,"/SampleInfo.RData"))

# Load packages
library(ggpubr)
library(tidyverse)

# Load gene data
setwd(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/"))
load("geneAnnotation.RData")
load("FPKM.RData")

# load protein data
setwd(paste0(homeDir,"/2. Proteomics/1. Preprocessing/"))
load("px_gx_annotations1.RData")
load("pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]

# prepare protein-gene annotations
px_gx_annotations <- unique(px_gx_annotations[,c("gene_id", "ID")])
px_gx_annotations_fil <- px_gx_annotations[(!is.na(px_gx_annotations$gene_id)) & 
                                             (!is.na(px_gx_annotations$ID)),]
colnames(px_gx_annotations_fil) <- c("Gene", "Protein")

# Check if all gene and proteins IDs have expression data available
all(px_gx_annotations_fil$Gene %in% rownames(FPKM))
all(px_gx_annotations_fil$Protein %in% rownames(pxMatrix_imp))


# Calculate gene-protein correlation for each sample
corValue <- rep(NA, ncol(FPKM))
for (i in 1:ncol(FPKM)){
  gene_temp <- data.frame(GeneExpr = FPKM[,i],
                          Gene = rownames(FPKM))
  
  protein_temp <- data.frame(ProteinExpr = pxMatrix_imp[,i],
                             Protein = rownames(pxMatrix_imp))
  
  all_temp <- inner_join(gene_temp, px_gx_annotations_fil, by = c("Gene" = "Gene"))
  all_temp <- inner_join(protein_temp, all_temp, by = c("Protein" = "Protein"))
  corValue[i] <- cor.test(all_temp$ProteinExpr, all_temp$GeneExpr, method = "spearman")$estimate
}

# Prepare data for plotting
corDF <- data.frame(corValue = corValue,
                    SampleID = colnames(FPKM))
corDF <- inner_join(corDF, sampleInfo, by = c("SampleID" = "SampleID"))
corDF$Tissue[corDF$Tissue == "Cell"] <- "iPSC"
corDF$TimeRegion <- paste0(corDF$Tissue, "_",corDF$Time)
corDF$TimeRegion <- factor(corDF$TimeRegion,levels = c("Dorsal_D75",
                                                         "Dorsal_D40",
                                                         "Dorsal_D13",
                                                         "iPSC_D0",
                                                         "Ventral_D13",
                                                         "Ventral_D40",
                                                         "Ventral_D75"))
corDF$ID <-  paste0(corDF$TimeRegion, "_", corDF$Gene)


# Main plot: expression over time
mainPlot <- ggplot(corDF) +
  geom_point(aes(x = TimeRegion, y = corValue, color = Group), size = 2) +
  geom_line(aes(x = TimeRegion, y = corValue, color = Group, group = paste0(Replicate, Group))) +
  ylab("Gene-protein correlation (\u03c1)") +
  xlab(NULL) +
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
  sample = corDF$TimeRegion,
  time = corDF$Time,
  region = corDF$Tissue))

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
ggsave(finalPlot, width = 6, height = 4,
       file = paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/GeneProteinCorrelations.png"))


################################################################################

# gene-wise correlations

################################################################################


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "E:/RTTproject/CellAnalysis/OrganoidAnalysis"
load(paste0(homeDir,"/SampleInfo.RData"))

# Load packages
library(ggpubr)
library(tidyverse)

# Load gene data
setwd(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/"))
load("geneAnnotation.RData")
load("FPKM.RData")

# load protein data
setwd(paste0(homeDir,"/2. Proteomics/1. Preprocessing/"))
load("px_gx_annotations1.RData")
load("pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]

# prepare protein-gene annotations
px_gx_annotations <- unique(px_gx_annotations[,c("gene_id", "ID")])
px_gx_annotations_fil <- px_gx_annotations[(!is.na(px_gx_annotations$gene_id)) & 
                                             (!is.na(px_gx_annotations$ID)),]
colnames(px_gx_annotations_fil) <- c("Gene", "Protein")

# Check if all gene and proteins IDs have expression data available
all(px_gx_annotations_fil$Gene %in% rownames(FPKM))
all(px_gx_annotations_fil$Protein %in% rownames(pxMatrix_imp))


# Calculate gene-protein correlation for each sample
corValue <- rep(NA, nrow(px_gx_annotations_fil))
for (i in 1:nrow(px_gx_annotations_fil)){
  all_temp <- data.frame(GeneExpr = FPKM[rownames(FPKM) == px_gx_annotations_fil$Gene[i],],
                          ProteinExpr = pxMatrix_imp[rownames(pxMatrix_imp) == px_gx_annotations_fil$Protein[i],],
                          Gene = colnames(FPKM))
  corValue[i] <- cor.test(all_temp$ProteinExpr, all_temp$GeneExpr, method = "spearman")$estimate
}

sum(corValue < 0)/length(corValue)
plotDF <- data.frame(value = corValue,
                     gene = px_gx_annotations_fil$Gene,
                     protein = px_gx_annotations_fil$Protein)

p <- ggplot(plotDF) +
  geom_histogram(aes(x = value), bins = 50, color = "white", fill = "#737373") +
  xlab("Gene-protein correlation (\u03c1)") +
  ylab("Count") +
  theme_bw()

ggsave(p, width = 6, height = 4,
       file = paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/GeneProteinCorrelations2.png"))

################################################################################

# logFC correlations

################################################################################


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "E:/RTTproject/CellAnalysis/OrganoidAnalysis"
load(paste0(homeDir,"/SampleInfo.RData"))

# Load packages
library(ggpubr)
library(tidyverse)

# Load gene data
setwd(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/"))
load("geneAnnotation.RData")
load("DEresults_RTTvsIC_gx.RData")

# load protein data
setwd(paste0(homeDir,"/2. Proteomics/1. Preprocessing/"))
load("px_gx_annotations1.RData")
load("DEresults_px.RData")
logFCs_px <- DEresults_px[,str_detect(colnames(DEresults_px), "ratio")]
pvalues_px <- DEresults_px[,str_detect(colnames(DEresults_px), "p.val")]


# prepare protein-gene annotations
px_gx_annotations <- unique(px_gx_annotations[,c("gene_id", "ID")])
px_gx_annotations_fil <- px_gx_annotations[(!is.na(px_gx_annotations$gene_id)) & 
                                             (!is.na(px_gx_annotations$ID)),]
colnames(px_gx_annotations_fil) <- c("Gene", "Protein")

# Check if all gene and proteins IDs have expression data available
all(px_gx_annotations_fil$Protein %in% DEresults_px$name)
all(px_gx_annotations_fil$Gene %in% rownames(topList[[1]]))



# Calculate gene-protein correlation for each sample
corValue <- rep(NA, length(topList))
for (i in 1:length(topList)){
  gene_temp <- data.frame(GeneExpr = topList[[i]]$logFC,
                          Gene = rownames(topList[[i]]))
  
  protein_temp <- data.frame(ProteinExpr = logFCs_px[,i],
                             Protein = DEresults_px$name)
  
  sig_genes <- rownames(topList[[i]])[p.adjust(topList[[i]]$PValue, method = "fdr") < 0.05]
  sig_proteins <- DEresults_px$name[p.adjust(pvalues_px[,i], method = "fdr")< 0.05]
  sig_annotations <- px_gx_annotations_fil[(px_gx_annotations_fil$Gene %in% sig_genes) | (px_gx_annotations_fil$Protein %in% sig_proteins),]
  all_temp <- inner_join(gene_temp, sig_annotations, by = c("Gene" = "Gene"))
  all_temp <- inner_join(protein_temp, all_temp, by = c("Protein" = "Protein"))
  corValue[i] <- cor.test(all_temp$ProteinExpr, all_temp$GeneExpr)$estimate
}



