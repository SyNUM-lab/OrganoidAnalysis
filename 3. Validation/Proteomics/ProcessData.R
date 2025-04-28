# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(DEP)
library(biomaRt)
library(readxl)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/Proteomics"))

# Read proteomics data
pxData_raw <- read_excel("Data/Proteomics_validation.xlsx")
pxMatrix_raw <- as.matrix(pxData_raw[,-c(1,2)])
rownames(pxMatrix_raw) <- as.vector(pxData_raw[,1])[[1]]

# Make sample information dataframe
sampleInfo <- data.frame(
  sampleID = colnames(pxMatrix_raw),
  Genotype = unlist(lapply(str_split(colnames(pxMatrix_raw), "_"),function(x)x[1])),
  Time = unlist(lapply(str_split(colnames(pxMatrix_raw), "_"),function(x)x[2])),
  Rep = unlist(lapply(str_split(colnames(pxMatrix_raw), "_"),function(x)x[3]))
)

save(sampleInfo, file = "Data/sampleInfo.RData")

# Exlude X and Y chromosomal proteins:

# Select biomaRt dataset
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# Get chromosome annotations for each protein
annotations <- getBM(attributes=c("uniprot_gn_id",
                                  "hgnc_symbol",
                                  "chromosome_name"), 
                     filters = 'uniprot_gn_id',
                     values = str_remove(rownames(pxMatrix_raw), "-.*"),
                     mart = ensembl)

save(annotations, file = "Data/ProteinAnnotations.RData")

# Get X/Y chromosomal proteins
XYproteins <- annotations$uniprot_gn_id[(annotations$chromosome_name == "X")|
                                          (annotations$chromosome_name == "Y")]
XYproteins1 <- NULL
for (i in 1:length(XYproteins)){
  XYproteins1 <- c(XYproteins1,
                   rownames(pxMatrix_raw)[str_detect(rownames(pxMatrix_raw),XYproteins[i])])
}

# Remove X/Y chromosomal proteins from the raw data
pxMatrix_raw <- pxMatrix_raw[!(rownames(pxMatrix_raw) %in% XYproteins1),]

#*****************************************************************************#
#   Filtering + Normalization
#*****************************************************************************#

# Make SummarizedExperiment object
pxData <- as.data.frame(pxMatrix_raw)
pxData$ID <- rownames(pxMatrix_raw)
pxData <- make_unique(pxData,ids="ID", names = "ID")
pxData <- make_se(pxData, 1:24, data.frame(label = sampleInfo$sampleID,
                                           condition = paste(sampleInfo$Genotype,
                                                             sampleInfo$Time,
                                                             sep = "_"),
                                           replicate = sampleInfo$Rep
))

plot_frequency(pxData)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(pxData, thr = 1)
plot_numbers(data_filt)
plot_coverage(data_filt)

# Normalize data
pxData_norm <- normalize_vsn(data_filt)
meanSdPlot(pxData_norm)
plot_normalization(data_filt, pxData_norm)

#*****************************************************************************#
#   Imputation
#*****************************************************************************#

# Extract protein names with missing values 
# in all replicates of at least one condition
# These are missing not at random
proteins_MNAR <- get_df_long(pxData_norm) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(pxData_norm) %in% proteins_MNAR

pxData_imp <- DEP::impute(
  pxData_norm, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "QRILC") 

plot_imputation(pxData_norm, pxData_imp)

# Save imputed data
save(pxData_imp, file = "Data/pxData_imp.RData")


#*****************************************************************************#
#   DE Analysis
#*****************************************************************************#

# Load expression data
load("Data/pxData_imp.RData")

# Differential enrichment analysis  based on linear models and empherical Bayes statistics
comparisons <- c("RTT_D3_vs_IC_D3",
                 "RTT_D9_vs_IC_D9",
                 "RTT_D15_vs_IC_D15",
                 "RTT_D22_vs_IC_D22"
)
data_diff <- test_diff(pxData_imp, type = "manual",test = comparisons)


# Denote significant proteins based on user defined cutoffs (adjusted p-value and log fold change)
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(2))

DEresults_px <- get_results(dep)

# Save DE results
save(DEresults_px, file = "Data/DEresults_px.RData")


#*****************************************************************************#
# PCA
#*****************************************************************************#

# Load expression data
load("Data/pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]

# Run PCA
pca <- prcomp(t(pxMatrix_imp), 
              retx = TRUE, # Give rotated variables (PCA loadings)
              center = TRUE,
              scale = TRUE) # Variables are scaled to unit variance

# Get expained variances
expl_var <- round((pca$sdev^2)/sum(pca$sdev^2),3)

# Combine sample info in PCA scores
plotPCA <- as.data.frame(pca$x)

# Add sample information to dataframe
all(rownames(plotPCA) == sampleInfo$sampleID)
plotPCA <- cbind(plotPCA, sampleInfo)

# Make plot
plotPCA$Colour <- paste0(plotPCA$Genotype, ": ", plotPCA$Time)

# Set colors
colors <- c("#6BAED6","#4292C6","#2171B5","#084594",
            "#FB6A4A","#EF3B2C","#CB181D","#99000D")

# Make plot
p <- ggplot(plotPCA, aes(x = PC1, y = PC2, color = Colour)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", expl_var[1]*100,"%)")) +
  ylab(paste0("PC2 (", expl_var[2]*100,"%)")) +
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
  guides(color = guide_legend("Group: Time"), shape = guide_legend("Region"))

ggsave(filename = "QCplots/PCAscores.png", p, width = 7, height = 5)

