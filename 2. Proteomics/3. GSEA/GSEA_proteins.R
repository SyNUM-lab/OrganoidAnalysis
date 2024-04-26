

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

# Set working directory
setwd("D:/RTTproject/CellAnalysis")

# Load data
load("Proteins/Preprocessing/pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load("Proteins/Preprocessing/DEresults_px.RData")
load("Proteins/Preprocessing/proteinAnnotation.RData")
all(annotations$uniprot_gn_id %in% str_remove(rownames(pxMatrix_imp), "-.*"))
load("Data/sampleInfo.RData")

# Get p-values
pvalues <- DEresults_px[,3:9]
rownames(pvalues) <- DEresults_px$name

# get adjusted p-values
adjpvalues <- pvalues
for (i in 1:ncol(pvalues)){
  adjpvalues[,i] <- p.adjust(pvalues[,i], "fdr")
}

# Get logFCs
logFCs <- DEresults_px[,str_detect(colnames(DEresults_px), "ratio")]
rownames(logFCs) <- DEresults_px$name

###############################################################################

# GSEA

###############################################################################


for (i in 1:ncol(adjpvalues)){
  gsea_input <- -log(pvalues[,i])*sign(logFCs[,i])
  names(gsea_input) <- str_remove(rownames(pvalues), "-.*")
  keep <- sort(abs(gsea_input), decreasing = TRUE)
  keep <- names(keep[!duplicated(names(keep))])
  gsea_input <- gsea_input[keep]
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  
  set.seed(123)
  GOtest <- gseGO(
    geneList = gsea_input,
    keyType = "UNIPROT",
    OrgDb = org.Hs.eg.db,
    ont =  "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500)
  
  if (i == 1){
    GOresults <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                          GOtest@result$ID, ")"),
                            pvalue = GOtest@result$pvalue)
    GOresults_NES <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                          GOtest@result$ID, ")"),
                            NES = GOtest@result$NES)
  }
  if (i > 1){
    temp <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                     GOtest@result$ID, ")"),
                       pvalue = GOtest@result$pvalue)
    GOresults <- full_join(GOresults,temp, by = c("Name" = "Name"))
    
    temp_NES <- data.frame(Name = paste0(GOtest@result$Description, " (",
                                     GOtest@result$ID, ")"),
                       NES = GOtest@result$NES)
    GOresults_NES <- full_join(GOresults_NES,temp_NES, by = c("Name" = "Name"))
  }
  
}

save(GOresults, file = "Proteins/RTTvsIC/GSEA/GOresults_GSEA.RData")
save(GOresults_NES, file = "Proteins/RTTvsIC/GSEA/GOresults_NES_GSEA.RData")


# Save genes and GO terms
GOannotation <- data.frame(Name =  paste0(GOtest@result$Description, " (",
                                          GOtest@result$ID, ")"),
                           ID = GOtest@result$ID,
                           Description = GOtest@result$Description)
save(GOannotation, file = "D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")

GOgenes <- GOtest@geneSets
save(GOgenes, file = "D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_UNIPROT.RData")


###############################################################################

# Make figures

###############################################################################


# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

# Set working directory
setwd("D:/RTTproject/CellAnalysis")

# Load GSEA results
load("Proteins/RTTvsIC/GSEA/GOresults_GSEA.RData")
load("Proteins/RTTvsIC/GSEA/GOresults_NES_GSEA.RData")

# Set column names of GO results
colnames(GOresults) <- c("Name", "Cell_D0", "Dorsal_D13", "Dorsal_D40", 
                         "Dorsal_D75", "Ventral_D13", "Ventral_D40", "Ventral_D75")
colnames(GOresults_NES) <- c("Name", "Cell_D0", "Dorsal_D13", "Dorsal_D40", 
                             "Dorsal_D75", "Ventral_D13", "Ventral_D40", "Ventral_D75")


# Get adjusted p-values
GOresults_adj <- apply(GOresults[,2:8],2,function(x){p.adjust(x, method = "fdr")})
rownames(GOresults_adj) <- GOresults$Name

# Load GO annotation data
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOgenes_UNIPROT.RData"))
load(paste0("D:/RTTproject/CellAnalysis/GO_annotation/","GOannotation.RData"))

# Get GO terms
rownames(GOannotation) <- GOannotation$Name
GOterms_all <- GOannotation[GOresults$Name, "ID"]

# Make similarity matrix of GO terms
simMatrix <- calculateSimMatrix(
  GOterms_all,
  orgdb="org.Hs.eg.db",
  ont = "BP",
  method = "Resnik"
)

# Select GO term with highest average -log10 p-value
scores <- setNames(rowSums(-log10(GOresults[,2:8])), GOterms_all)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.8,
                                orgdb="org.Hs.eg.db")

save(reducedTerms, file = "Proteins/RTTvsIC/GSEA/reducedTerms_RTTvsIC_BP.RData")

load("Proteins/RTTvsIC/GSEA/reducedTerms_RTTvsIC_BP.RData")

# Remove GO:0045229, because it includes the same genes as GO:0030198
#reducedTerms <- reducedTerms[reducedTerms$parent != "GO:0045229",]

# Reduce number of GO terms based on previous selection
GOannotation_fil <- GOannotation[GOannotation$ID %in% unique(reducedTerms$parent),]
GOresults_fil <- GOresults[GOresults$Name %in% GOannotation_fil$Name,]
GOresults_NES_fil <- GOresults_NES[GOresults_NES$Name %in% GOannotation_fil$Name,]
GOresults_adj_fil <- GOresults_adj[rownames(GOresults_adj) %in% GOannotation_fil$Name,]


# Get terms that are significant (FDR < 0.05) in at least 1 time point
sigterms <- rownames(GOresults_adj_fil)[rowSums(GOresults_adj_fil < 0.05) > 2]


# Clustering
values <- -log10(GOresults[,2:8]) * sign(GOresults_NES[,2:8])
rownames(values) <- GOresults$Name
values_fil <- values[sigterms,]
GOresults_adj_fil <- GOresults_adj[sigterms,]

model <- hclust(dist(values_fil), "ward.D2")
terms_ordered <- model$labels[model$order]
save(terms_ordered, file = "Proteins/RTTvsIC/GSEA/terms_ordered.RData")

# Prepare data for plotting
plotData <- gather(values_fil)
plotData$Name <- rep(rownames(values_fil),7)

# combine with GO annotaiton
plotData <- inner_join(plotData, GOannotation, by = c('Name' = 'Name'))

# Combine with adjusted pvalues
GOresults_adj_fil <- as.data.frame(GOresults_adj_fil)
GOresults_adj_fil <- gather(GOresults_adj_fil[,1:7])
plotData <- cbind.data.frame(plotData, GOresults_adj_fil$value)
plotData$Sig <- ifelse(plotData$`GOresults_adj_fil$value` < 0.05, "Yes", "No")

# Add sample info
plotData$Region <- NA
plotData$Region[str_detect(plotData$key, "Cell")] <- "iPSC"
plotData$Region[str_detect(plotData$key, "Dorsal")] <- "Dorsal"
plotData$Region[str_detect(plotData$key, "Ventral")] <- "Ventral"

plotData$Time <- "D0"
plotData$Time[str_detect(plotData$key, "D13")] <- "D13"
plotData$Time[str_detect(plotData$key, "D40")] <- "D40"
plotData$Time[str_detect(plotData$key, "D75")] <- "D75"

plotData$key <- factor(plotData$key,
                       levels = c("Dorsal_D75",
                                  "Dorsal_D40",
                                  "Dorsal_D13",
                                  "Cell_D0",
                                  "Ventral_D13",
                                  "Ventral_D40",
                                  "Ventral_D75"))

plotData$Description <- factor(firstup(plotData$Description), 
                               levels = firstup(GOannotation[terms_ordered, "Description"]))
plotData$Name <- factor(plotData$Name, levels = terms_ordered)


# Make main plot: Heatmap
mainPlot <- ggplot(data = plotData, aes(x = key, y = Description, fill = value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(.~Region, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "Signed\n-log10 p-value") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_discrete(breaks = levels(plotData$Description),
                   labels = str_replace(levels(plotData$Description),
                                        "transesterification reactions with bulged adenosine as nucleophile",
                                        "...*"))

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotData$key,
  time = plotData$Time,
  tissue = plotData$Region))

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


# Combine plots into a single figure
finalPlot <- ggarrange(colSideColorPlot_tissue,
                       mainPlot,
                       colSideColorPlot_time,
                       heights = c(0.5,9,0.5), widths = c(2,8),nrow = 3,ncol = 1,
                       align = "hv",
                       common.legend = FALSE)

# Save plot
ggsave(finalPlot, file = "Proteins/RTTvsIC/GSEA/GSEA_proteins.png", width = 8, height = 8)


# Get legends:
legendPlot <- ggplot(data = plotData, aes(x = key, y = Description, fill = value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(.~Region, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "") +
  guides(color = "none") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))


legend <- cowplot::get_legend(legendPlot)

ggsave(legend, file = "Proteins/RTTvsIC/GSEA/GSEA_legend_proteins.png", width = 8, height = 8)
