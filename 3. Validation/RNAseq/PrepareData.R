
# Load packages
library(tidyverse)
library(GEOquery)
library(data.table)
library(tidyverse)
library(stringr)
library(edgeR)
library(limma)
library(statmod)

################################################################################
#
#   GSE117511
#
################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/RNAseq"))

#=============================================================================#
# Load raw data
#=============================================================================#

accessID <- "GSE117511"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Get gene annotation
annotation <- fread("Human.GRCh38.p13.annot.tsv")
rownames(annotation) <- as.character(annotation$GeneID)

gene_names <- rep(NA, nrow(tbl))
for (i in 1:nrow(tbl)){
  gene_names[i] <- annotation[which(rownames(annotation) %in% rownames(tbl)[i])[1], "Symbol"]
}
rownames(tbl) <- gene_names

# Get meta data
test <- getGEO(accessID)
sampleInfo <- test$GSE117511_series_matrix.txt.gz@phenoData@data[,c(1:2,11:14)]
colnames(sampleInfo) <- c("Title", "SampleID", "Genotype", "Stimulation", "Treatment", "Time")
sampleInfo$Genotype <- str_remove(sampleInfo$Genotype, "genotype/variation: ")
sampleInfo$Stimulation <- str_remove(sampleInfo$Stimulation, "stimulation: ")
sampleInfo$Treatment <- str_remove(sampleInfo$Treatment, "treatment: ")
sampleInfo$Time <- str_remove(sampleInfo$Time, "time: ")

# Filter for no treatment and stimulation
sampleInfo <- sampleInfo[(sampleInfo$Stimulation == "None") & (sampleInfo$Treatment == "None"),]
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = "GSE117511/gxMatrix_raw.RData")
save(sampleInfo, file = "GSE117511/sampleInfo.RData")

#=============================================================================#
# Pre-processing and statistical analysis
#=============================================================================#

# Load data
load("GSE117511/gxMatrix_raw.RData")
load("GSE117511/sampleInfo.RData")
all(sampleInfo$SampleID == colnames(gxMatrix))
sampleInfo$Genotype <- ifelse(sampleInfo$Genotype == "WT", "WT", "RTT")

# Make DGE list object
d <- DGEList(counts = gxMatrix,
             group = paste0(sampleInfo$Time,"_", sampleInfo$Genotype))

# remove lowly expressed miRNAs
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]

# Normalize the data
d <- calcNormFactors(d)
gxMatrix_norm <- log2(cpm(d) + 1)

# Save normalized data
save(gxMatrix_norm, file = "GSE117511/gxMatrix_norm.RData")

# Make design matrix
time <- sampleInfo$Time
genotype <- sampleInfo$Genotype
design <- model.matrix(~0+genotype+time)

# Estimate dispersion
d <- estimateDisp(d, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(d, design)

# Compare RTT vs WT
contr <- makeContrasts(genotypeRTT- genotypeWT,
                       levels = colnames(coef(fit)))

# Test for significance
test <- glmQLFTest(fit, contrast = contr[,1])
top.table <- topTags(test, n = Inf)$table

# Save statistics
save(top.table, file = "GSE117511/top.table.RData")

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

# Save plot
ggsave(filename = "GSE117511/PCA_Scores.png", p, width = 7, height = 5)



################################################################################
#
#   GSE107399
#
################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/RNAseq"))

#=============================================================================#
# Load raw data
#=============================================================================#

accessID <- "GSE107399"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Load gene annotation data
annotation <- fread("Human.GRCh38.p13.annot.tsv")
rownames(annotation) <- as.character(annotation$GeneID)
gene_names <- rep(NA, nrow(tbl))
for (i in 1:nrow(tbl)){
  gene_names[i] <- annotation[which(rownames(annotation) %in% rownames(tbl)[i])[1], "Symbol"]
}
rownames(tbl) <- gene_names

# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE107399_series_matrix.txt.gz@phenoData@data[5:30,c(1:2,8,12)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Genotype")
sampleInfo$Genotype <- str_remove(sampleInfo$Genotype, "phenotype: ")

samples <- intersect(colnames(tbl),sampleInfo$SampleID)
sampleInfo <- sampleInfo[samples,]

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))

#=============================================================================#
# Preprocessing and statistical analysis
#=============================================================================#

# Load data
load("GSE107399/gxMatrix_raw.RData")
load("GSE107399/sampleInfo.RData")
all(sampleInfo$SampleID == colnames(gxMatrix))

# Add replicate to sample information
sampleInfo$Subject <- NA
sampleInfo$Subject[str_detect(sampleInfo$Title, "Normal")] <- "Normal"
sampleInfo$Subject[str_detect(sampleInfo$Title, "Patient982")] <- "Patient982"
sampleInfo$Subject[str_detect(sampleInfo$Title, "Patient567")] <- "Patient567"

# Make DGE list object
d <- DGEList(counts = gxMatrix,
             group = paste0(sampleInfo$Tissue,"_", sampleInfo$Genotype))

# remove lowly expressed miRNAs
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]

# Normalize the data
d <- calcNormFactors(d)
gxMatrix_norm <- log2(cpm(d) + 1)
save(gxMatrix_norm, file = "GSE107399/gxMatrix_norm.RData")
save(sampleInfo, file = "GSE107399/sampleInfo.RData")

# Make design matrix
group <- factor(paste(sampleInfo$Genotype,sampleInfo$Tissue,sep="."))
design <- model.matrix(~0+group+sampleInfo$Subject)
colnames(design) <- c(levels(group), "Patient567", "Patient982")

# Estimate dispersion
d <- estimateDisp(d, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(d, design)

# Compare RTT vs WT
contr <- makeContrasts(RTT.iPSC-WT.iPSC, 
                       RTT.NPC - WT.NPC,
                       RTT.Neuron - WT.Neuron,
                       levels = colnames(coef(fit)))

# Test for significance
test <- glmQLFTest(fit, contrast = contr[,1])
top.table_iPSC <- topTags(test, n = Inf)$table
head(top.table_iPSC , 20)
save(top.table_iPSC, file = "GSE107399/top.table_iPSC.RData")

test <- glmQLFTest(fit, contrast = contr[,2])
top.table_NPC <- topTags(test, n = Inf)$table
head(top.table_NPC , 20)
save(top.table_NPC, file = "GSE107399/top.table_NPC.RData")

test <- glmQLFTest(fit, contrast = contr[,3])
top.table_Neuron <- topTags(test, n = Inf)$table
head(top.table_Neuron , 20)
save(top.table_Neuron, file = "GSE107399/top.table_Neuron.RData")

#==============================================================================#
# Quality control
#==============================================================================#

# Load data
load("GSE107399/gxMatrix_norm.RData")
load("GSE107399/sampleInfo.RData")

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
plotPCA$Colour <- paste0(plotPCA$Genotype, ": ", plotPCA$Tissue)

colors <- rev(c("#084594","#99000D"))
p <- ggplot(plotPCA, aes(x = PC1, y = PC2, color = Genotype, shape = Tissue)) +
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
  guides(color = guide_legend("Time"), shape = guide_legend("Tissue"))

# Save plot
ggsave(filename = "GSE107399/PCA_Scores.png", p, width = 7, height = 5)


################################################################################
#
#   GSE123753
#
################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/RNAseq"))

#=============================================================================#
# Load raw data
#=============================================================================#
accessID <- "GSE123753"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE123753_series_matrix.txt.gz@phenoData@data[,c(1:2,8,11)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Genotype")
sampleInfo$Genotype <- str_remove(sampleInfo$Genotype, "genotype: ")
sampleInfo <- sampleInfo[!str_detect(sampleInfo$Title, "TRAP"),]
sampleInfo <- sampleInfo[sampleInfo$Tissue != "Induced pluripotent stem cells",]

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))

#=============================================================================#
# Preprocessing and statistical analysis
#=============================================================================#

# Load data
load("GSE123753/gxMatrix_raw.RData")
load("GSE123753/sampleInfo.RData")
all(sampleInfo$SampleID == colnames(gxMatrix))

# Make DGE list object
d <- DGEList(counts = gxMatrix,
             group = paste0(sampleInfo$Tissue,"_", sampleInfo$Genotype))

# remove lowly expressed miRNAs
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]

# Normalize the data
d <- calcNormFactors(d)
gxMatrix_norm <- log2(cpm(d) + 1)
save(gxMatrix_norm, file = "GSE123753/gxMatrix_norm.RData")

# Make design matrix
group <- factor(paste(sampleInfo$Tissue,sampleInfo$Genotype,sep="."))
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# Estimate dispersion
d <- estimateDisp(d, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(d, design)

# Compare RTT vs WT
contr <- makeContrasts(NPC.RTT - NPC.WT,
                       Neu.RTT - Neu.WT,
                       levels = colnames(coef(fit)))

# Test for significance
test <- glmQLFTest(fit, contrast = contr[,1])
top.table_NPC <- topTags(test, n = Inf)$table
head(top.table_NPC , 20)
save(top.table_NPC, file = "GSE123753/top.table_NPC.RData")

test <- glmQLFTest(fit, contrast = contr[,2])
top.table_Neuron <- topTags(test, n = Inf)$table
head(top.table_Neuron , 20)
save(top.table_Neuron, file = "GSE123753/top.table_Neuron.RData")


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
p <- ggplot(plotPCA, aes(x = PC1, y = PC2, color = Genotype, shape = Tissue)) +
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

# Save plot
ggsave(filename = "GSE123753/PCA_Scores.png", p, width = 7, height = 5)


################################################################################
#
#   GSE128380
#
################################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/RNAseq"))

#=============================================================================#
# Load raw data
#=============================================================================#
accessID <- "GSE128380"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

annotation <- fread("Human.GRCh38.p13.annot.tsv")
rownames(annotation) <- as.character(annotation$GeneID)

gene_names <- rep(NA, nrow(tbl))
for (i in 1:nrow(tbl)){
  gene_names[i] <- annotation[which(rownames(annotation) %in% rownames(tbl)[i])[1], "Symbol"]
}
rownames(tbl) <- gene_names

# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE128380_series_matrix.txt.gz@phenoData@data[,c(1:2,8,11)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Age")
sampleInfo$Age <- str_remove(sampleInfo$Age, "age: ")
sampleInfo$Age <- as.numeric(str_remove(sampleInfo$Age, " years"))

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))

#=============================================================================#
# Preprocessing and statistical analysis
#=============================================================================#

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
save(gxMatrix_norm, file = "GSE128380/gxMatrix_norm.RData")

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
p <- ggplot(plotPCA, aes(x = PC1, y = PC2, color = Genotype, shape = Tissue)) +
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

# Save plot
ggsave(filename = "GSE128380/PCA_Scores.png", p, width = 7, height = 5)

