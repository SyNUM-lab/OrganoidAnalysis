# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(stringr)
library(edgeR)
library(patchwork)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/"))

#*****************************************************************************#
# Normalization
#*****************************************************************************#

# Load data
load("gxMatrix_raw1.RData")                                           # Raw count data
load("geneAnnotation.RData")                                          # Gene annotation data
load(paste0(homeDir, "/sampleInfo.RData"))                            # Sample Information
all(sampleInfo$SampleID == colnames(gxMatrix_raw))

# Remove sex chromosomal genes
gxMatrix_fil <- gxMatrix_raw[rownames(gxMatrix_raw) %in% geneAnnotation$gene_id[(geneAnnotation$chromosome_name != "X") & (geneAnnotation$chromosome_name != "Y")],]

# Combine region/tissue (iPSC, Dorsal, Ventral), time (D0, D13, D40, and D75), and group (IC and RTT)
sampleInfo$Time_Group <- paste0(sampleInfo$Tissue, "_", sampleInfo$Time, "_", sampleInfo$Group)
time_group <- factor(sampleInfo$Time_Group)

# Make DGE list object
y <- DGEList(counts = gxMatrix_fil,
             group = time_group)

# remove lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)
gxMatrix_norm <- log2(cpm(y) + 1)
save(gxMatrix_norm, file = "gxMatrix_norm.RData")

#*****************************************************************************#
# DE analysis: RTT vs IC
#*****************************************************************************#

# Make design matrix
design <- model.matrix(~0 + time_group)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(y, design)

# Make Contrasts: RTT vs IC at every time point
all_contrasts <- makeContrasts(
  Cell_D0_RTTvsIC = time_groupCell_D0_RTT-time_groupCell_D0_IC,
  Dorsal_D13_RTTvsIC = time_groupDorsal_D13_RTT-time_groupDorsal_D13_IC,
  Dorsal_D40_RTTvsIC = time_groupDorsal_D40_RTT-time_groupDorsal_D40_IC,
  Dorsal_D75_RTTvsIC = time_groupDorsal_D75_RTT-time_groupDorsal_D75_IC,
  Ventral_D13_RTTvsIC = time_groupVentral_D13_RTT-time_groupVentral_D13_IC,
  Ventral_D40_RTTvsIC = time_groupVentral_D40_RTT-time_groupVentral_D40_IC,
  Ventral_D75_RTTvsIC = time_groupVentral_D75_RTT-time_groupVentral_D75_IC,
  levels = design
)

# Perform statistical analysis for each comparison and collect
# top table in a list object
topList <- list()
for (i in 1:7){
  # Test for significance
  test <- glmQLFTest(fit, contrast = all_contrasts[,i])
  
  #Get top table
  topTable <- topTags(test, n = Inf)$table
  topTable$ID <- rownames(topTable)
  
  # Collect top table into list
  topList[[i]] <- topTable
}

names(topList) <- colnames(all_contrasts)

# Save differential expression analysis results
save(topList, file = "DEresults_RTTvsIC_gx.RData")


#*****************************************************************************#
# DE analysis: Time points
#*****************************************************************************#

# Make design matrix
design <- model.matrix(~0 + time_group)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(y, design)

# Make Contrasts: RTT vs IC at every time point
all_contrasts <- makeContrasts(
  Dorsal_D13_vs_Cell_D0 = time_groupDorsal_D13_IC-time_groupCell_D0_IC,
  Dorsal_D40_vs_Dorsal_D13 = time_groupDorsal_D40_IC-time_groupDorsal_D13_IC,
  Dorsal_D40_vs_Cell_D0 = time_groupDorsal_D40_IC-time_groupCell_D0_IC,
  Dorsal_D75_vs_Dorsal_D40 = time_groupDorsal_D75_IC-time_groupDorsal_D40_IC,
  Dorsal_D75_vs_Dorsal_D13 = time_groupDorsal_D75_IC-time_groupDorsal_D13_IC,
  Dorsal_D75_vs_Cell_D0 = time_groupDorsal_D75_IC-time_groupCell_D0_IC,
  
  Ventral_D13_vs_Cell_D0 = time_groupVentral_D13_IC-time_groupCell_D0_IC,
  Ventral_D40_vs_Ventral_D13 = time_groupVentral_D40_IC-time_groupVentral_D13_IC,
  Ventral_D40_vs_Cell_D0 = time_groupVentral_D40_IC-time_groupCell_D0_IC,
  Ventral_D75_vs_Ventral_D40 = time_groupVentral_D75_IC-time_groupVentral_D40_IC,
  Ventral_D75_vs_Ventral_D13 = time_groupVentral_D75_IC-time_groupVentral_D13_IC,
  Ventral_D75_vs_Cell_D0 = time_groupVentral_D75_IC-time_groupCell_D0_IC,
  levels = design
)

# Perform statistical analysis for each comparison and collect
# top table in a list object
topList <- list()
for (i in 1:12){
  # Test for significance
  test <- glmQLFTest(fit, contrast = all_contrasts[,i])
  
  #Get top table
  topTable <- topTags(test, n = Inf)$table
  topTable$ID <- rownames(topTable)
  
  # Collect top table into list
  topList[[i]] <- topTable
}

names(topList) <- colnames(all_contrasts)

# Save differential expression analysis results
save(topList, file = "DEresults_Time_gx.RData")



#*****************************************************************************#
# PCA
#*****************************************************************************#

# Load data
load("gxMatrix_norm.RData")
load("geneAnnotation.RData")
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/sampleInfo.RData")

# Run PCA
pca <- prcomp(t(gxMatrix_norm), 
              retx = TRUE, # Give rotated variables (PCA loadings)
              center = TRUE,
              scale = TRUE) # Variables are scaled to unit variance

# Get expained variances
expl_var <- round((pca$sdev^2)/sum(pca$sdev^2),3)

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
  labs(color = "Group: Time", shape = "Region") +
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
                                       linetype = 1))

# Save PCA plot
ggsave(filename = "QCplots/PCA_Scores_gx.png", p, width = 7, height = 5)


#*****************************************************************************#
#   Boxplots
#*****************************************************************************#

# Format data for plotting
exprPlot <- gather(as.data.frame(gxMatrix_norm))
exprPlot <- inner_join(exprPlot, sampleInfo, by = c("key" = "SampleID"))

# Prepare data
exprPlot$Colour <- paste0(exprPlot$Group, ": ", exprPlot$Time)
exprPlot$Tissue[exprPlot$Tissue == "Cell"] <- "iPSC"
exprPlot$Tissue <- factor(exprPlot$Tissue, levels = c("Dorsal","iPSC",  "Ventral"))

# Order the samples
order <- arrange(exprPlot, by = Time)
orderSamples <- unique(order$key)
exprPlot$key <- factor(exprPlot$key,
                       levels = c(rev(orderSamples[str_detect(orderSamples, "Dorsal")]),
                                  orderSamples[str_detect(orderSamples, "D0")],
                                  orderSamples[str_detect(orderSamples, "Ventral")]
                                  ))


# Set colors
colors <- c("#6BAED6","#4292C6","#2171B5","#084594",
            "#FB6A4A","#EF3B2C","#CB181D","#99000D")

# Make boxplot
p <- ggplot(exprPlot, aes(x = key, y = value, fill = Colour)) +
  facet_grid(cols = vars(Tissue), scale = "free", space = "free") +
  geom_boxplot(alpha = 0.8) +
  ylab(expression(log[2]~"CPM")) +
  xlab("")+
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        legend.position = "none", 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"),
        strip.text = element_blank())

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = exprPlot$key[(exprPlot$Group == "IC") & (exprPlot$Replicate == 2)],
  time = exprPlot$Time[(exprPlot$Group == "IC") & (exprPlot$Replicate == 2)],
  tissue = exprPlot$Tissue[(exprPlot$Group == "IC") & (exprPlot$Replicate == 2)]))


# Tissue/region
colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[(str_detect(colSideColor$sample, "D40") |
                                   str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# Combine plots
finalPlot <- p + colSideColorPlot_time + colSideColorPlot_tissue +
  plot_layout(nrow = 3, ncol = 1, heights = c(8.6,0.7,0.7))

# Save plot
ggsave(finalPlot, file = "QCplots/Boxplots_gx.png", width = 8.5, height = 6)

# Get legends:
legendPlot <- ggplot(exprPlot, aes(x = key, y = value, fill = Colour)) +
  facet_grid(cols = vars(Tissue), scale = "free", space = "free") +
  geom_boxplot(alpha = 0.8) +
  ylab(expression(log[2]~"CPM")) +
  xlab("")+
  labs(fill = "Group: Time") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "right", 
        legend.title = element_text(size = 30),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"),
        strip.text = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))

legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "QCplots/Boxplots_gx_legend.png", width = 8, height = 8)

