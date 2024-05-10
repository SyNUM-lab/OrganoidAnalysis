# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(stringr)
library(edgeR)
library(limma)

# Set working directory
setwd("D:/RTTproject/PublicDatasets")

# Load data
load("GSE128380/gxMatrix_raw.RData")
load("GSE128380/sampleInfo.RData")
all(sampleInfo$SampleID == colnames(gxMatrix))

sampleInfo$Genotype <- ifelse(str_detect(sampleInfo$Title,"RTT"), "RTT", "WT")
sampleInfo$Tissue <- ifelse(sampleInfo$Tissue == "Temporal cortex", "TC", "CC")

#==============================================================================#
# Statistical analysis
#==============================================================================#

# Make DGE list object
d <- DGEList(counts = gxMatrix,
             group = paste0(sampleInfo$Tissue,"_", sampleInfo$Genotype))

# remove lowly expressed miRNAs
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]

# Normalize the data
d <- calcNormFactors(d)
gxMatrix_norm <- log2(cpm(d) + 1)
save(gxMatrix_norm, file = "GSE117511/gxMatrix_norm.RData")

# Make design matrix
group <- paste0(sampleInfo$Genotype, "_", sampleInfo$Tissue)
age <- sampleInfo$Age
design <- model.matrix(~ 0 + group + age)

# Estimate dispersion
d <- estimateDisp(d, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(d, design)


# Compare RTT vs WT
contr <- makeContrasts(groupRTT_CC- groupWT_CC, 
                       groupRTT_TC- groupWT_TC, 
                       levels = colnames(coef(fit)))

# Test for significance
test <- glmQLFTest(fit, contrast = contr[,1])
top.table_CC <- topTags(test, n = Inf)$table
save(top.table_CC, file = "GSE128380/top.table_CC.RData")

# Test for significance
test <- glmQLFTest(fit, contrast = contr[,2])
top.table_TC <- topTags(test, n = Inf)$table
save(top.table_TC, file = "GSE128380/top.table_TC.RData")



#==============================================================================#
# Quality control
#==============================================================================#

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
colors <- rev(c("#084594","#99000D"))
p <- ggplot(plotPCA, aes(x = PC1, y = PC2, color = Genotype, shape = Time)) +
  geom_point(size = 3) +
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
  guides(color = guide_legend("Genotype"), shape = guide_legend("Cell type"))

ggsave(filename = "GSE117511/PCA_Scores.png", p, width = 7, height = 5)


#==============================================================================#
# Plotting
#==============================================================================#


