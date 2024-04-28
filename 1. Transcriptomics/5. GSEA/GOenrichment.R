# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(ggpubr)
library(biomaRt)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/5. GSEA")

# Load data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/SampleInfo.RData")

genes_all <- rownames(topList[[1]])

# Get p-values
pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues) <- unlist(lapply(str_split(genes_all, "_"), function(x){x[1]}))
colnames(pvalues) <- names(topList)

# get adjusted p-values
adj_pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"FDR"]})), ncol = length(topList))
rownames(adj_pvalues) <- unlist(lapply(str_split(genes_all, "_"), function(x){x[1]}))
colnames(adj_pvalues) <- names(topList)

# Get logFCs
logFCs <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
rownames(logFCs) <- unlist(lapply(str_split(genes_all, "_"), function(x){x[1]}))
colnames(logFCs) <- names(topList)


###############################################################################

# GSEA

###############################################################################

# Perform GSEA at each time point/region
for (i in 1:ncol(adj_pvalues)){
  gsea_input <- -log(pvalues[,i])*sign(logFCs[,i])
  names(gsea_input) <- str_remove(rownames(pvalues), "-.*")
  keep <- sort(abs(gsea_input), decreasing = TRUE)
  keep <- names(keep[!duplicated(names(keep))])
  gsea_input <- gsea_input[keep]
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  
  set.seed(123)
  GOtest <- gseGO(
    geneList = gsea_input,
    keyType = "ENSEMBL",
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

colnames(GOresults) <- c("Name", colnames(pvalues))
colnames(GOresults_NES) <- c("Name", colnames(pvalues))

# Save results
save(GOresults, file ="Data/GOresults_GSEA_gx1.RData")
save(GOresults_NES, file = "Data/GOresults_NES_GSEA_gx1.RData")

###############################################################################

# Analyse results

###############################################################################

# Load GSEA results
load("Data/GOresults_GSEA_gx1.RData")
load("Data/GOresults_NES_GSEA_gx1.RData")

# Get adjusted p-values
GOresults_adj <- apply(GOresults[,2:8],2,function(x){p.adjust(x, method = "fdr")})
rownames(GOresults_adj) <- GOresults$Name

# Load GO annotation data
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData")
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOannotation.RData")

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
save(simMatrix, file = "Data/simMatrix_RTTvsIC_BP1.RData")


# Select GO term with highest average -log10 p-value
scores <- setNames(rowSums(-log10(GOresults[,2:8])), GOterms_all)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.8,
                                orgdb="org.Hs.eg.db")

# Remove GO:0045229, because it includes the same genes as GO:0030198
#reducedTerms <- reducedTerms[reducedTerms$parent != "GO:0045229",]

# Save reduced terms
save(reducedTerms, file = "Data/reducedTerms_RTTvsIC_BP1.RData")

# Load data
load("Data/reducedTerms_RTTvsIC_BP1.RData")

# Reduce number of GO terms based on previous selection
GOannotation_fil <- GOannotation[GOannotation$ID %in% unique(reducedTerms$parent),]
GOresults_fil <- GOresults[GOresults$Name %in% GOannotation_fil$Name,]
GOresults_NES_fil <- GOresults_NES[GOresults_NES$Name %in% GOannotation_fil$Name,]
GOresults_adj_fil <- GOresults_adj[rownames(GOresults_adj) %in% GOannotation_fil$Name,]

# Get terms that are significant (FDR < 0.05) in at least 4 time point
sigterms <- rownames(GOresults_adj_fil)[rowSums(GOresults_adj_fil < 0.05) > 3]

# Clustering
values <- -log10(GOresults[,2:8]) * sign(GOresults_NES[,2:8])
rownames(values) <- GOresults$Name
values_fil <- values[sigterms,]
GOresults_adj_fil <- GOresults_adj[sigterms,]

values_fil1 <- values_fil
values_fil1[values_fil1 > 3] <- 3
values_fil1[values_fil1 < -3] <- -3
model <- hclust(dist(values_fil1), "ward.D2")
terms_ordered <- model$labels[model$order]
save(terms_ordered, file = "Data/terms_ordered1.RData")

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
                       levels = c("Dorsal_D75_RTTvsIC",
                                  "Dorsal_D40_RTTvsIC",
                                  "Dorsal_D13_RTTvsIC",
                                  "Cell_D0_RTTvsIC",
                                  "Ventral_D13_RTTvsIC",
                                  "Ventral_D40_RTTvsIC",
                                  "Ventral_D75_RTTvsIC"))

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
        strip.text.x = element_blank())

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
ggsave(finalPlot, file = "Plots/GOEnrichment_RTTvsIC1.png", width = 8, height = 8)

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

ggsave(legend, file = "Plots/legend_GOEnrichment1.png", width = 8, height = 8)


###############################################################################

# Plot Individual GO term

###############################################################################

# Load GSEA results
load("Data/GOresults_GSEA_gx1.RData")
load("Data/GOresults_NES_GSEA_gx1.RData")

# Format data
plotExpr <- gather(as.data.frame(gxMatrix_norm))
plotExpr$GeneID <- rep(rownames(gxMatrix_norm),ncol(gxMatrix_norm))
plotExpr <- inner_join(plotExpr, sampleInfo, by = c("key" = 'SampleID'))
plotExpr$SampleGene <- paste(plotExpr$Time, 
                                plotExpr$Tissue, 
                                plotExpr$Group, 
                                plotExpr$GeneID, sep = "_")

plotExpr <- plotExpr %>% 
  group_by(SampleGene) %>%
  mutate(meanExpr = mean(value))

plotExpr <- plotExpr[,-c(1,2,7)]
plotExpr <- plotExpr[!duplicated(plotExpr),]

# Standardize
plotExpr <- plotExpr %>%
  group_by(GeneID) %>%
  mutate(stExpr = (meanExpr - mean(meanExpr))/sd(meanExpr))

# Load GO annotation data
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData")
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOannotation.RData")

# Select GO term
termID <- "GO:0097192"
# termName <- "regulation of postsynaptic membrane potential"
#termID <- GOannotation$ID[GOannotation$Description == termName]

# Get name of GO term
termName <- GOannotation$Description[GOannotation$ID == termID]

# Get genes of GO term
genes <- GOgenes[[termID]]
geneIDs <- geneAnnotation$gene_id[which(geneAnnotation$EnsemblID%in% genes)]

# Plot expression over time
plotDF <- plotExpr[plotExpr$GeneID %in% geneIDs,]

plotDF$Tissue[plotDF$Tissue == "Cell"] <- "iPSC"
plotDF$TimeRegion <- paste0(plotDF$Tissue, "_",plotDF$Time)
plotDF$TimeRegion <- factor(plotDF$TimeRegion,levels = c("Dorsal_D75",
                                                         "Dorsal_D40",
                                                         "Dorsal_D13",
                                                         "iPSC_D0",
                                                         "Ventral_D13",
                                                         "Ventral_D40",
                                                         "Ventral_D75"))

medianExpr <- plotDF %>%
  group_by(TimeRegion, Group) %>%
  summarise(test = median(stExpr))

# Main plot: expression over time
mainPlot <- ggplot(plotDF) +
  geom_point(aes(x = TimeRegion, y = stExpr, color = Group), size = 2, shape = "-") +
  geom_line(aes(x = TimeRegion, y = stExpr, color = Group, 
                group = paste0(Group, GeneID)), alpha = 0.3) +
  geom_line(data = medianExpr, aes(x = TimeRegion, y = test, color = Group,
                                   group = Group), linewidth = 1) +
  ggtitle(termName) +
  ylab(expression("Standardized"~log[2]~"CPM")) +
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
  sample = plotDF$TimeRegion,
  time = plotDF$Time,
  region = plotDF$Tissue))

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
ggsave(finalPlot, file = paste0("Plots/",termName,"_expr.png"), width = 10, height = 8)




###############################################################################

# Plot overlap of GO terms

###############################################################################


# Load data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/SampleInfo.RData")
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData")
load("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/GOannotation.RData")

# Get selected GO terms
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

# Prepare data
finalDF <- inner_join(TERM2GENE, geneAnnotation, by = c("GENE" = "EnsemblID"))
finalDF <- inner_join(finalDF, GOannotation, by = c("TERM" = "ID"))
finalDF$Name <- factor(firstup(finalDF$Name), levels = firstup(terms_ordered))

# Cluster genes
values <- matrix(0, nrow = length(unique(finalDF$TERM)), ncol = length(unique(finalDF$GENE)))
rownames(values) <- unique(finalDF$TERM)
colnames(values) <- unique(finalDF$GENE)
for (i in 1:nrow(finalDF)){
  values[finalDF$TERM[i], finalDF$GENE[i]] <- 1
}

model <- hclust(dist(t(values)), "ward.D2")
genes_ordered <- model$labels[model$order]
finalDF$GENE <- factor(finalDF$GENE, levels = genes_ordered)

# Make plot
p <- ggplot(finalDF) +
  geom_tile(aes(x = GENE, y = Name)) +
  ylab(NULL) +
  xlab("Genes") +
  theme_bw() +
  theme(axis.text.x = element_blank())

# Save plot
ggsave(p, file = paste0("Plots/common_genes.png"), width = 12, height = 8)
