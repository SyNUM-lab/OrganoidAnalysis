
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(DEP)
library(biomaRt)

# Set working directory
setwd("D:/RTTproject/CellAnalysis/Proteins/Preprocessing")
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")

# Read proteomics data
pxData_raw <- read.delim("D:/RTTproject/OriginalData/Proteomics/DEP analysis/Result files/20210310_Cell__total_scaling_raw.txt",
                     as.is= TRUE)

# Filter for contaminant proteins and low confidence hits
pxData_raw <- pxData_raw[!pxData_raw$Contaminant,]
pxData_raw <- pxData_raw[pxData_raw$Score.Sequest.HT>10,]
rownames(pxData_raw) <- pxData_raw[,1]

# Remove irrelevant columns
pxData_raw <- pxData_raw[,-c(1:9, 52:59)]

# Remove empty rows
pxData_raw <- pxData_raw[!apply(pxData_raw == "", 1, all),]

# Make dataframe numeric
pxData_raw <- mutate_all(pxData_raw, na_if,"")
for (i in 1:ncol(pxData_raw)){
  pxData_raw[,i] <- str_replace_all(pxData_raw[,i],",", ".")
  pxData_raw[,i] <- as.numeric(pxData_raw[,i])
}

# Chance column names
samples <- str_remove_all(colnames(pxData_raw), "\\.")
samples1 <- rep(0, length(samples))
for (i in 1:length(samples)){
  
  #IC or RTT
  if (str_detect(samples[i],"IC")){
    a <- "IC_"
  } else{
    a <- ""
  }
  
  a <- paste0(a, "MeCP2_R255X_")
  
  #Day
  if (str_detect(samples[i],"D0")){
    a <- paste0(a, "D0")
  } else if (str_detect(samples[i], "D13")){
    a <- paste0(a, "D13")
  } else if (str_detect(samples[i], "D40")){
    a <- paste0(a, "D40")
  } else if (str_detect(samples[i], "D75")){
    a <- paste0(a, "D75")
  }
  
  #Tissue
  if (str_detect(samples[i],"Dorsal")){
    a <- paste0(a, "_Dorsal")
  } else if (str_detect(samples[i], "Ventral")){
    a <- paste0(a, "_Ventral")
  }
  
  #Replicate
  if (str_detect(samples[i],"n1")){
    a <- paste0(a, "_1")
  } else if (str_detect(samples[i], "n2")){
    a <- paste0(a, "_2")
  } else if (str_detect(samples[i], "n3")){
    a <- paste0(a, "_3")
  }
  samples1[i] <- a
}

all(samples1 %in% sampleInfo$SampleID)


colnames(pxData_raw) <- samples1
pxData_raw <- pxData_raw[,sampleInfo$SampleID]
all(sampleInfo$SampleID == colnames(pxData_raw))

# Exlude X and Y chromosomal proteins
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

#biomaRt::listAttributes(ensembl)
annotations <- getBM(attributes=c("uniprot_gn_id",
                                  "hgnc_symbol",
                                  "chromosome_name"), 
                     filters = 'uniprot_gn_id',
                     values = str_remove(rownames(pxData_raw), "-.*"),
                     mart = ensembl)

XYproteins <- annotations$uniprot_gn_id[(annotations$chromosome_name == "X")|
                                          (annotations$chromosome_name == "Y")]
XYproteins1 <- NULL
for (i in 1:length(XYproteins)){
  XYproteins1 <- c(XYproteins1,
                  rownames(pxData_raw)[str_detect(rownames(pxData_raw),XYproteins[i])])
}


pxData_raw <- pxData_raw[!(rownames(pxData_raw) %in% XYproteins1),]

#*****************************************************************************#
#   Filtering + Normalization
#*****************************************************************************#

# Make SummarizedExperiment object
pxData <- pxData_raw
pxData$ID <- rownames(pxData_raw)
pxData <- make_unique(pxData,ids="ID", names = "ID")
pxData <- make_se(pxData, 1:42, data.frame(label = sampleInfo$SampleID,
                                                     condition = paste(sampleInfo$Group,
                                                                       sampleInfo$Tissue,
                                                                       sampleInfo$Time,
                                                                       sep = "_"),
                                                     replicate = sampleInfo$Replicate
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

save(pxData_imp, file = "Proteins/pxData_imp.RData")

#*****************************************************************************#
#   DE Analysis
#*****************************************************************************#

load("Proteins/pxData_imp.RData")

# Differential enrichment analysis  based on linear models and empherical Bayes statistics
comparisons <- c("RTT_Cell_D0_vs_IC_Cell_D0",
                 "RTT_Dorsal_D13_vs_IC_Dorsal_D13",
                 "RTT_Dorsal_D40_vs_IC_Dorsal_D40",
                 "RTT_Dorsal_D75_vs_IC_Dorsal_D75",
                 "RTT_Ventral_D13_vs_IC_Ventral_D13",
                 "RTT_Ventral_D40_vs_IC_Ventral_D40",
                 "RTT_Ventral_D75_vs_IC_Ventral_D75"
                 )
data_diff <- test_diff(pxData_imp, type = "manual",test = comparisons)


# Denote significant proteins based on user defined cutoffs (adjusted p-value and log fold change)
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(2))

DEresults_px <- get_results(dep)
save(DEresults_px, file = "Proteins/DEresults_px.RData")



#*****************************************************************************#
# PCA
#*****************************************************************************#
load("pxData_imp.RData")

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
if (all(rownames(plotPCA) == paste(sampleInfo$Group,
                                   sampleInfo$Tissue,
                                   sampleInfo$Time,
                                   sampleInfo$Replicate,
                                   sep = "_"))){
  plotPCA <- cbind(plotPCA, sampleInfo)
}


# Make plot
plotPCA$Colour <- paste0(plotPCA$Group, ": ", plotPCA$Time)
plotPCA$Tissue[plotPCA$Tissue == "Cell"] <- "iPSC"
plotPCA$Tissue <- factor(plotPCA$Tissue, levels = c("iPSC", "Dorsal", "Ventral"))

colors <- c("#6BAED6","#4292C6","#2171B5","#084594",
            "#FB6A4A","#EF3B2C","#CB181D","#99000D")
plot <- ggplot(plotPCA, aes(x = PC1, y = PC2, color = Colour, shape = Tissue)) +
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

ggsave(filename = "PCAscores_Proteins.png", plot, width = 7, height = 5)


#*****************************************************************************#
#   Boxplots
#*****************************************************************************#
load("pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]

# Format data for plotting
exprPlot <- gather(as.data.frame(pxMatrix_imp))
sampleInfo$SampleID1 <- paste(sampleInfo$Group,
                              sampleInfo$Tissue,
                              sampleInfo$Time,
                              sampleInfo$Replicate,
                              sep = "_")
exprPlot <- inner_join(exprPlot, sampleInfo, by = c("key" = "SampleID1"))

exprPlot$Colour <- paste0(exprPlot$Group, ": ", exprPlot$Time)
exprPlot$Tissue[exprPlot$Tissue == "Cell"] <- "iPSC"
exprPlot$Tissue <- factor(exprPlot$Tissue, levels = c("Dorsal","iPSC",  "Ventral"))

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
  ylab(expression(log[2]~"intensity")) +
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
ggsave(finalPlot, file = "Boxplots_Proteins.png", width = 8.5, height = 6)


# Make density plot
exprBoxplot <- ggplot() +
  geom_density(data = exprPlot, aes(x = value, color = key)) +
  ylab(expression(log[2]~"intensity")) +
  xlab("")+
  theme_bw() +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic")) +
  scale_color_viridis_d()

# Save plot
ggsave(exprBoxplot, file = "Density_Proteins.png", width = 8, height = 5)


# Make protein annotation data
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

#biomaRt::listAttributes(ensembl)
annotations1 <- getBM(attributes=c("uniprot_gn_id",
                                  "hgnc_symbol",
                                  "chromosome_name",
                                  "entrezgene_id"), 
                     filters = 'uniprot_gn_id',
                     values = str_remove(rownames(pxMatrix_imp), "-.*"),
                     mart = ensembl)

annotations <- annotations1[!duplicated(annotations1[,c(1,2,4)]),]

rowAnn <- data.frame(ID = rownames(pxMatrix_imp),
                     Name = str_remove(rownames(pxMatrix_imp), "-.*"))

annotations <- right_join(annotations, rowAnn, by = c("uniprot_gn_id" = "Name"))
all(rownames(pxMatrix_imp) %in% annotations$ID)

save(annotations, file = "Preprocessing/Proteins/proteinAnnotation.RData")

