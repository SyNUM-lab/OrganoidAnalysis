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
setwd("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing")

# Load data
load("IsoPctMatrix.RData")
load("txMatrix_raw.RData")
load("txAnnotation.RData")
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_raw1.RData")

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


ggsave(filename = "PCA_Scores_IsoPct.png", p, width = 7, height = 5)


#==============================================================================#
# differential splicing analysis
#==============================================================================#

# Filter lowly expressed transcripts
rownames(sampleInfo) <- sampleInfo$SampleID
RTT_samples <- sampleInfo$SampleID[sampleInfo$Group == "RTT"]
IC_samples <- sampleInfo$SampleID[sampleInfo$Group == "IC"]

keep_RTT <- rownames(txMatrix_raw)[rowSums(txMatrix_raw[,RTT_samples]<10)==0]
keep_IC <- rownames(txMatrix_raw)[rowSums(txMatrix_raw[,IC_samples]<10)==0]
IsoPctMatrix_fil <- IsoPctMatrix[union(keep_RTT,keep_IC),]

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
    
    # Collect RTT expression values
    X_RTT <- IsoPctMatrix_fil[i,samplesRTT]
    
    # Collect IC expession values
    X_IC <- IsoPctMatrix_fil[i,samplesIC]
    
    # independent 2-group Mann-Whitney U Test
    pvalue[i,j] <- wilcox.test(X_RTT, X_IC)[["p.value"]]
    
    # mean difference
    diff[i,j] <- mean(X_RTT) - mean(X_IC)
    
  } 
}
colnames(diff) <- paste0(comps$Time, "_", comps$Tissue)
rownames(diff) <- rownames(IsoPctMatrix_fil)

colnames(pvalue) <- paste0(comps$Time, "_", comps$Tissue)
rownames(pvalue) <- rownames(IsoPctMatrix_fil)

save(pvalue,diff,file = "DASResults.RData")

load("DASResults.RData")
sel <- rownames(pvalue)[which(rowSums(pvalue <= 0.1) > 3)]

t <- table(txAnnotation$gene_id[txAnnotation$transcript_id %in% sel]) > 1
DASgenes <- names(t)[t]
save(DASgenes, file = "DASgenes.RData")

sel <- rownames(pvalue)[which(rowSums(pvalue <= 0.1) > 6)]
t <- table(txAnnotation$gene_id[txAnnotation$transcript_id %in% sel]) > 1
DASgenes <- names(t)[t]

#==============================================================================#
# Length isoforms
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
setwd("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing")

# Load data
load("IsoPctMatrix.RData")
load("txMatrix_raw.RData")
load("txAnnotation.RData")
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_norm.RData")
load("DASgenes.RData")
load("DASResults.RData")

load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/geneAnnotation.RData")
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOgenes_ENSEMBL.RData"))
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOannotation.RData"))

translation_genes <- geneAnnotation$gene_id[geneAnnotation$EnsemblID %in% GOgenes[["GO:0000375"]]]

plotDF <- NULL
for (time in colnames(diff)){
  sel <- rownames(pvalue)[which((pvalue[,time] <= 0.1))]
  sel_up <- rownames(pvalue)[which((pvalue[,time] <= 0.1) & (diff[,time] > 0))]
  sel_down <- rownames(pvalue)[which((pvalue[,time] <= 0.1) & (diff[,time] < 0))]
  
  # Select genes with at least one up and one downregulated isform
  genes_up <- unique(txAnnotation$gene_id[txAnnotation$transcript_id %in% sel_up])
  genes_down <- unique(txAnnotation$gene_id[txAnnotation$transcript_id %in% sel_down])
  genes_sel <- unique(intersect(genes_up,genes_down))
  genes_sel <- intersect(genes_sel,translation_genes)
  
  for (g in genes_sel){
    t_up <- intersect(txAnnotation$transcript_id[txAnnotation$gene_id == g], sel_up)
    t_down <- intersect(txAnnotation$transcript_id[txAnnotation$gene_id == g], sel_down)
    for (t in t_up){
      temp <- data.frame(Up = rep(t,length(t_down)),
                         Down = t_down,
                         Up_value = rep(txAnnotation$length[txAnnotation$transcript_id %in% t],length(t_down)),
                         Down_value = txAnnotation$length[txAnnotation$transcript_id %in% t_down],
                         TimeRegion = rep(time,length(t_down))
      )
      plotDF <- rbind.data.frame(plotDF, temp)
    }
  }
  
}

#[1] "D0_Cell"     "D13_Dorsal"  "D13_Ventral" "D40_Dorsal"  "D40_Ventral" "D75_Dorsal" 
# [7] "D75_Ventral"
mean(plotDF$Up_value[plotDF$TimeRegion =="D0_Cell"]/plotDF$Down_value[plotDF$TimeRegion =="D0_Cell"])
median(plotDF$Up_value[plotDF$TimeRegion =="D0_Cell"]/plotDF$Down_value[plotDF$TimeRegion =="D0_Cell"])


mean(log(plotDF$Up_value[plotDF$TimeRegion =="D13_Dorsal"]/plotDF$Down_value[plotDF$TimeRegion =="D13_Dorsal"]))
median(plotDF$Up_value[plotDF$TimeRegion =="D13_Dorsal"]/plotDF$Down_value[plotDF$TimeRegion =="D13_Dorsal"])


ggplot(plotDF) +
  geom_point(aes(y = log(Down_value), x = log(Up_value))) +
  geom_abline(intercept = 0, slope = 1, color = "red")



#==============================================================================#
# Investigate cytoplasmic translation
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
setwd("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing")

# Load data
load("IsoPctMatrix.RData")
load("txMatrix_raw.RData")
load("txAnnotation.RData")
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_norm.RData")
load("DASgenes.RData")
load("DASResults.RData")


load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/geneAnnotation.RData")
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOgenes_ENSEMBL.RData"))
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOannotation.RData"))

translation_genes <- geneAnnotation$gene_id[geneAnnotation$EnsemblID %in% GOgenes[["GO:0002181"]]]
translation_genes <- intersect(DASgenes,translation_genes)



sel <- rownames(pvalue)[which(rowSums(pvalue <= 0.1) > 3)]
#txAnnotation_sel <- txAnnotation[txAnnotation$transcript_id %in% sel,]
txAnnotation_sel <- txAnnotation[txAnnotation$gene_id %in% translation_genes,]
txAnnotation_sel <- txAnnotation_sel[txAnnotation_sel$gene_id %in% rownames(gxMatrix_norm),]

txAnnotation_sel$Value <- rep(NA, nrow(txAnnotation_sel))
all(colnames(gxMatrix_norm) == colnames(IsoPctMatrix))

for (i in 1:nrow(txAnnotation_sel)){
  IsoPct <- (IsoPctMatrix[txAnnotation_sel$transcript_id[i],]/100)
  gxExpr <- 2^(gxMatrix_norm[txAnnotation_sel$gene_id[i],])-1
  
  txAnnotation_sel$Value[i] <- cor(IsoPct,gxExpr, method = "spearman")
  
}

ggplot(txAnnotation_sel) +
  geom_point(aes(x = transcript_id, y = Value)) +
  facet_grid(cols = vars(gene_id), scale = "free", space = "free")


#==============================================================================#
# Make plots
#==============================================================================#

# [1] "ENSG00000100744_GSKIP"    "ENSG00000106484_MEST"     "ENSG00000114796_KLHL24"  
# [4] "ENSG00000120071_KANSL1"   "ENSG00000127022_CANX"     "ENSG00000129195_PIMREG"  
# [7] "ENSG00000129484_PARP2"    "ENSG00000143727_ACP1"     "ENSG00000147274_RBMX"    
# [10] "ENSG00000158373_H2BC5"    "ENSG00000163682_RPL9"     "ENSG00000168653_NDUFS5"  
# [13] "ENSG00000170291_ELP5"     "ENSG00000198453_ZNF568"   "ENSG00000198825_INPP5F"  
# [16] "ENSG00000221983_UBA52"    "ENSG00000224078_SNHG14"   "ENSG00000250337_PURPL"   
# [19] "ENSG00000274512_TBC1D3L"  "ENSG00000285756"          "ENSG00000289047"         
# [22] "ENSG00000290376_HERC2P3"  "ENSG00000291003"          "ENSG00000291093_SVIL-AS1"

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
setwd("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing")

# Load data
load("IsoPctMatrix.RData")
load("txMatrix_raw.RData")
load("txAnnotation.RData")
load("D:/RTTproject/CellAnalysis/Data/sampleInfo.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_norm.RData")
load("DASgenes.RData")


load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/geneAnnotation.RData")
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOgenes_ENSEMBL.RData"))
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOannotation.RData"))

translation_genes <- geneAnnotation$gene_id[geneAnnotation$EnsemblID %in% GOgenes[["GO:0002181"]]]
intersect(DASgenes,translation_genes)

# Put samples in correct order
samples <- sampleInfo$SampleID
gxMatrix_norm <- gxMatrix_norm[,samples]
IsoPctMatrix <- IsoPctMatrix[,samples]
sampleInfo$Tissue[sampleInfo$Tissue == "Cell"] <- "iPSC"

#gene <- DASgenes[2]
gene <- "ENSG00000092964_DPYSL2" 
transcripts <- txAnnotation$transcript_id[txAnnotation$gene_id == gene]
geneExpr <- 2^gxMatrix_norm[gene,]-1
IsoPct <- IsoPctMatrix[transcripts,]/100

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

plotDF$tid <- unlist(lapply(str_split(plotDF$tid,"_"),function(x)x[[2]]))

plotDF1 <- plotDF %>%
  group_by(sampleid) %>%
  mutate(sumExpr = sum(IsoExprValue))
plotDF1 <- unique(plotDF1[,4:11])

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


p_all <- p1 + p2 + colSideColorPlot_time + colSideColorPlot_tissue +
  patchwork::plot_layout(ncol = 1, nrow = 4,
                         heights = c(3,10,0.5,0.5))
p_all

ggsave(p_all, file = "DPYSL2_splicing.png", width = 8, height = 7)
