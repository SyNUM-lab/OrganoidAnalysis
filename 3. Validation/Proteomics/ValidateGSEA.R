# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(stringr)
library(DEP)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/Proteomics"))

# Load data
load("Data/DEresults_px.RData")
load("Data/ProteinAnnotations.RData")
load("Data/sampleInfo.RData")
load("Data/pxData_imp.RData")
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_UNIPROT_Hs.RData"))
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
rownames(GOannotation) <- GOannotation$Name
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

# Get p-values
pvalues <- DEresults_px[,3:6]
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

# Prepare TERM2GENE dataframe
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(term = substr(names(unlist(GOgenes_fil)),1,10),
                        gene = unlist(GOgenes_fil))

# Perform GSEA for each time point/region
for (i in 1:ncol(adjpvalues)){
  gsea_input <- -log(pvalues[,i])*sign(logFCs[,i])
  names(gsea_input) <- str_remove(rownames(pvalues), "-.*")
  keep <- sort(abs(gsea_input), decreasing = TRUE)
  keep <- names(keep[!duplicated(names(keep))])
  gsea_input <- gsea_input[keep]
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  
  set.seed(123)
  GOtest <- GSEA(
    geneList = gsea_input,
    pAdjustMethod = "fdr",
    pvalueCutoff = 1,
    minGSSize = 1,
    maxGSSize = 1000,
    TERM2GENE = TERM2GENE)
  
  if (i == 1){
    GOresults <- data.frame(Name = GOtest@result$ID,
                            pvalue = GOtest@result$pvalue)
    GOresults_NES <- data.frame(Name =GOtest@result$ID,
                                NES = GOtest@result$NES)
  }
  if (i > 1){
    temp <- data.frame(Name = GOtest@result$ID,
                       pvalue = GOtest@result$pvalue)
    GOresults <- full_join(GOresults,temp, by = c("Name" = "Name"))
    
    temp_NES <- data.frame(Name = GOtest@result$ID,
                           NES = GOtest@result$NES)
    GOresults_NES <- full_join(GOresults_NES,temp_NES, by = c("Name" = "Name"))
  }
  
}

# Set column names
colnames(GOresults) <- c("GOname", "D15", "D22", "D3", "D9")
colnames(GOresults_NES) <- c("GOname","D15", "D22", "D3", "D9")

# Combine with GO annotation data
GOresults_NES <- inner_join(GOannotation, GOresults_NES,
                            by = c("ID" = "GOname"))
GOresults <- inner_join(GOannotation, GOresults,
                            by = c("ID" = "GOname"))

# Save GSEA results
save(GOresults, file = "Data/GOresults_GSEA.RData")
save(GOresults_NES, file = "Data/GOresults_NES_GSEA.RData")


###############################################################################

# Make validation heatmap

###############################################################################

# Function to capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Load data
load("Data/GOresults_GSEA.RData")
load("Data/GOresults_NES_GSEA.RData")

# Calculate signed -log10 P value
values <- -log10(GOresults[,4:7]) * sign(GOresults_NES[,4:7])
rownames(values) <- GOresults$Name
values_fil <- values[terms_ordered,]

# Prepare data for plotting
plotData <- gather(values_fil)
plotData$Name <- rep(rownames(values_fil),4)

# Combine with GO annotation
plotData <- inner_join(plotData, GOannotation, by = c('Name' = 'Name'))
plotData$Sig <- ifelse(abs(plotData$value) > -log10(0.05), "Yes", "No")


# Order rows
plotData$key <- factor(plotData$key,
                       levels = c("D3", "D9", "D15", "D22"))

# Order columns
plotData$Description <- factor(firstup(plotData$Description), 
                               levels = firstup(GOannotation[terms_ordered, "Description"]))
plotData$Name <- factor(plotData$Name, levels = terms_ordered)

# Make main plot: Heatmap
mainPlot <- ggplot(data = plotData, aes(x = key, y = Description, fill = value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", 
                       midpoint = 0, trans = "pseudo_log", lim = c(-6,6),
                       oob = scales::squish) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "Signed\n-log10 p-value") +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())


# Time
colSideColorPlot_time <- ggplot(data = data.frame(Time = factor(c("D3", "D9", "D15", "D22"),
                                                                levels = c("D3", "D9", "D15", "D22")))) +
  geom_tile(aes(x = Time, y = "label", fill = Time)) +
  #geom_text(aes(x = Time, y = "label", label = Time)) +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  #scale_fill_brewer(palette = "Oranges") +
  scale_fill_manual(values = rep("#D95F02",4))

# Make final plot
finalPlot <- ggarrange(colSideColorPlot_time, mainPlot,
                       heights = c(0.5,9), nrow = 2,ncol = 1,
                       align = "v",
                       common.legend = FALSE)

# Save final plot
ggsave(finalPlot, file = "GOEnrichment_ProteomicsValidation.png", width = 6, height = 8)



################################################################################

# Permutation analysis

################################################################################

# Prepare TERM2GENE dataframe
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(term = substr(names(unlist(GOgenes_fil)),1,10),
                        gene = unlist(GOgenes_fil))


# Perform 1000-permutation
results <- NULL
set.seed(345)
nPerm <- 1000
for (i in 1:nPerm){
  
  gsea_input <- sample(1:nrow(pvalues),nrow(pvalues))
  names(gsea_input) <- str_remove(rownames(pvalues), "-.*")
  keep <- sort(abs(gsea_input), decreasing = TRUE)
  keep <- names(keep[!duplicated(names(keep))])
  gsea_input <- gsea_input[keep]
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  
  suppressWarnings({
    GOtest <- GSEA(
      geneList = gsea_input,
      pAdjustMethod = "fdr",
      pvalueCutoff = 1,
      minGSSize = 1,
      maxGSSize = 1000,
      TERM2GENE = TERM2GENE)
  })
  
  
  # Get results
  temp <- data.frame(pvalue = GOtest@result$pvalue,
                     term = GOtest@result$ID,
                     Perm = i)
  
  results <- rbind.data.frame(results, temp)
  
}

# Save permutation results
save(results, file = "Data/permResults_validation1.RData")

# Load data (if needed)
load("permResults_validation1.RData")

# Make histogram to check uniform distribution
hist(results$pvalue)

# How many permutations have more than four significant GO terms?
results_all <- NULL
for (j in unique(results$term)){
  results_fil <- results[results$term == j,]
  results_fil$permN <- 1:nrow(results_fil)
  
  results_all <- rbind.data.frame(results_all, results_fil)
}

# Get significance
sum(table(results$Perm[results$pvalue<0.05])>5)
