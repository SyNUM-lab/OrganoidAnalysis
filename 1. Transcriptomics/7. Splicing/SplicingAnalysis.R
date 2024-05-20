# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(data.table)
library(biomaRt)
library(edgeR)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/7. Splicing"))

# Load data
load("Data/IsoPctMatrix.RData")
load("Data/txMatrix_raw.RData")
load("Data/txAnnotation.RData")
load(paste0(homeDir,"/sampleInfo.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_raw1.RData"))

#==============================================================================#
# PCA
#==============================================================================#

# Remove transcripts with zero variance
IsoPctMatrix_fil <- IsoPctMatrix[apply(IsoPctMatrix,1,var)>0,]

# Run PCA
pca <- prcomp(t(IsoPctMatrix_fil), 
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


# Save plot
ggsave(filename = "PCA_Scores_IsoPct.png", p, width = 7, height = 5)


#==============================================================================#
# differential splicing analysis
#==============================================================================#

# Filter lowly expressed transcripts:

# 1) transcript with a count of 10 or more in at least one of the experimental group
rownames(sampleInfo) <- sampleInfo$SampleID
RTT_samples <- sampleInfo$SampleID[sampleInfo$Group == "RTT"]
IC_samples <- sampleInfo$SampleID[sampleInfo$Group == "IC"]
keep_RTT <- rownames(txMatrix_raw)[rowSums(txMatrix_raw[,RTT_samples]<10)==0]
keep_IC <- rownames(txMatrix_raw)[rowSums(txMatrix_raw[,IC_samples]<10)==0]
IsoPctMatrix_fil <- IsoPctMatrix[union(keep_RTT,keep_IC),]

# 2) Transcripts from genes with a count of at least 50
keep_genes <- rownames(gxMatrix_raw)[rowSums(gxMatrix_raw < 50)==0]
keep_transcripts <- txAnnotation$transcript_id[txAnnotation$gene_id %in% keep_genes]
IsoPctMatrix_fil <- IsoPctMatrix_fil[rownames(IsoPctMatrix_fil) %in% keep_transcripts,]


# Remove transcripts with zero variance
IsoPctMatrix_fil <- IsoPctMatrix_fil[apply(IsoPctMatrix_fil,1,var)>0,]

# Remove transcripts that are the only transcript of a gene
sel_genes <- names(table(txAnnotation$gene_id))[table(txAnnotation$gene_id)>1]
sel_transcripts <- txAnnotation$transcript_id[txAnnotation$gene_id %in% sel_genes]
IsoPctMatrix_fil <- IsoPctMatrix_fil[rownames(IsoPctMatrix_fil) %in% sel_transcripts,]

# Compare alternative splicing between groups
comps <- unique(sampleInfo[,c("Time", "Tissue")])
pvalue <- matrix(NA, nrow = nrow(IsoPctMatrix_fil), ncol = nrow(comps))
diff <- matrix(NA, nrow = nrow(IsoPctMatrix_fil), ncol = nrow(comps))

for (j in 1:nrow(comps)){
  time <- comps$Time[j]
  tissue <- comps$Tissue[j]
  
  for (i in 1:nrow(IsoPctMatrix_fil)){
    
    # Get RTT samples
    samplesRTT <- sampleInfo$SampleID[(sampleInfo$Time == time 
                                       & sampleInfo$Tissue == tissue &
                                         sampleInfo$Group == "RTT")]
    
    # Get IC samples
    samplesIC <- sampleInfo$SampleID[(sampleInfo$Time == time 
                                      & sampleInfo$Tissue == tissue &
                                        sampleInfo$Group == "IC")]
    
    # Collect RTT isoform percentage values
    X_RTT <- IsoPctMatrix_fil[i,samplesRTT]
    
    # Collect IC isoform percentage values
    X_IC <- IsoPctMatrix_fil[i,samplesIC]
    
    # Independent 2-group Mann-Whitney U Test:
    # Is the isoform percentage different between the experimental groups?
    pvalue[i,j] <- wilcox.test(X_RTT, X_IC)[["p.value"]]
    
    # mean difference
    diff[i,j] <- mean(X_RTT) - mean(X_IC)
    
  } 
}
colnames(diff) <- paste0(comps$Time, "_", comps$Tissue)
rownames(diff) <- rownames(IsoPctMatrix_fil)

colnames(pvalue) <- paste0(comps$Time, "_", comps$Tissue)
rownames(pvalue) <- rownames(IsoPctMatrix_fil)

# Save results
save(pvalue,diff,file = "Data/DASResults.RData")

# Get genes that are differentially alternative spliced:
# genes with transcripts with a P <= 0.1 in more than three time points/regions
load("Data/DASResults.RData")
sel <- rownames(pvalue)[which(rowSums(pvalue <= 0.1) > 3)]
t <- table(txAnnotation$gene_id[txAnnotation$transcript_id %in% sel]) > 1
DASgenes <- names(t)[t]

# Save DAS genes
save(DASgenes, file = "Data/DASgenes.RData")

#==============================================================================#
# Make plots
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(data.table)
library(biomaRt)
library(edgeR)
library(patchwork)
library(scales)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/7. Splicing"))

# Load data
load("Data/IsoPctMatrix.RData")
load("Data/txMatrix_raw.RData")
load("Data/txAnnotation.RData")
load("Data/DASgenes.RData")
load(paste0(homeDir,"/sampleInfo.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_norm.RData"))
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/geneAnnotation.RData"))
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData"))

# Put samples in correct order
samples <- sampleInfo$SampleID
gxMatrix_norm <- gxMatrix_norm[,samples]
IsoPctMatrix <- IsoPctMatrix[,samples]
sampleInfo$Tissue[sampleInfo$Tissue == "Cell"] <- "iPSC"

# Select gene
gene <- "ENSG00000092964_DPYSL2" 
transcripts <- txAnnotation$transcript_id[txAnnotation$gene_id == gene]
geneExpr <- 2^gxMatrix_norm[gene,]-1
IsoPct <- IsoPctMatrix[transcripts,]/100

# Collect isoform percentage and expression values into data frame
IsoPctValue <- NULL
IsoExprValue <- NULL
tid <- NULL
sampleid <- NULL
for (t in 1:length(transcripts)){
  temp <- IsoPct[t,]
  IsoPctValue <- c(IsoPctValue, temp)
  temp <- geneExpr*IsoPct[t,]
  IsoExprValue <- c(IsoExprValue, temp)
  tid <- c(tid,rep(transcripts[t],length(temp)))
  sampleid <- c(sampleid,names(geneExpr))
}

plotDF <- data.frame(IsoPctValue = IsoPctValue,
                     IsoExprValue = IsoExprValue,
                     tid = tid,
                     sampleid = sampleid)

# Combine with sample information
plotDF <- inner_join(plotDF,sampleInfo, by = c("sampleid" = "SampleID"))
plotDF$TimeRegion <- paste0(plotDF$Tissue, "_",plotDF$Time)
plotDF$TimeRegion <- factor(plotDF$TimeRegion,levels = c("Dorsal_D75",
                                                         "Dorsal_D40",
                                                         "Dorsal_D13",
                                                         "iPSC_D0",
                                                         "Ventral_D13",
                                                         "Ventral_D40",
                                                         "Ventral_D75"))
plotDF <- arrange(plotDF, by = TimeRegion)
plotDF$sampleid <- factor(plotDF$sampleid, levels = unique(plotDF$sampleid))

# Get transcript ID
plotDF$tid <- unlist(lapply(str_split(plotDF$tid,"_"),function(x)x[[2]]))

# Get total gene expression value
plotDF1 <- plotDF %>%
  group_by(sampleid) %>%
  mutate(sumExpr = sum(IsoExprValue))
plotDF1 <- unique(plotDF1[,4:11])

# Plot of gene expression values
p1 <- ggplot(plotDF1) +
  geom_bar(aes(x = sampleid,y=sumExpr),
           stat = "identity") +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
  ylab("CPM") +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y  = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 13),
        strip.text = element_text(size = 12),
        panel.grid.minor = element_blank())


# Plot of isoform percentages
p2 <- ggplot(plotDF) +
  geom_bar(aes(x = sampleid,y=IsoPctValue,fill = tid), 
           color = "black", width = 0.95,
           stat = "identity") +
  geom_hline(yintercept = 1, color = "black", linewidth = 1) +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
  ylab("Relative\nisoform expression") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y  = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 16),
        strip.text = element_blank(),
        legend.title = element_blank()
        #,legend.position = "none"
        ) 

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotDF$sampleid,
  time = plotDF$Time,
  region = plotDF$Tissue,
  Group = plotDF$Group))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor[str_detect(colSideColor$sample, "_2"),],
            aes(x = sample, y = "label", label = time)) +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
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
  geom_text(data = colSideColor[(str_detect(colSideColor$sample, "_2")) & 
                                   (str_detect(colSideColor$sample, "D40") |
                                   str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = region)) +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# Combine all subplots
p_all <- p1 + p2 + colSideColorPlot_time + colSideColorPlot_tissue +
  patchwork::plot_layout(ncol = 1, nrow = 4,
                         heights = c(3,10,0.5,0.5))
p_all

# Save plot
ggsave(p_all, file = "DPYSL2_splicing.png", width = 8, height = 7)
