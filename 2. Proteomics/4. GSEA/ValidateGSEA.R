
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

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/2. Proteomics/4. GSEA"))

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Load data
load("GOresults_GSEA.RData")
load("GOresults_NES_GSEA.RData")
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_UNIPROT_Hs.RData"))
rownames(GOannotation) <- GOannotation$Name
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

# Set column names of GO results
colnames(GOresults) <- c("Name", "Cell_D0", "Dorsal_D13", "Dorsal_D40", 
                         "Dorsal_D75", "Ventral_D13", "Ventral_D40", "Ventral_D75")
colnames(GOresults_NES) <- c("Name", "Cell_D0", "Dorsal_D13", "Dorsal_D40", 
                         "Dorsal_D75", "Ventral_D13", "Ventral_D40", "Ventral_D75")


values <- -log10(GOresults[,2:8]) * sign(GOresults_NES[,2:8])
rownames(values) <- GOresults$Name
values_fil <- values[terms_ordered,]

# Prepare data for plotting
plotData <- gather(values_fil)
plotData$Name <- rep(rownames(values_fil),7)

# combine with GO annotation
plotData <- inner_join(plotData, GOannotation, by = c('Name' = 'Name'))

# Combine with adjusted P values
pvalues <- GOresults
rownames(pvalues) <- GOresults$Name
pvalues <- gather(pvalues[terms_ordered,2:8])
plotData <- cbind.data.frame(plotData, pvalues$value)
plotData$Sig <- ifelse(plotData$`pvalues$value` < 0.05, "Yes", "No")

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
ggsave(finalPlot, file = "GOEnrichment_ProteomicsValidation.png", width = 8, height = 8)


# Make legend
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
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20))

legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "legend_GOEnrichment_ProteomicsValidation.png", width = 8, height = 8)


################################################################################

# Permutation analysis

################################################################################

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

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
      exponent = 1,
      minGSSize = 1,
      maxGSSize = 10000,
      eps = 1e-10,
      pvalueCutoff = Inf,
      pAdjustMethod = "fdr",
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      verbose = FALSE
    )
  })
  

  # Get results
  temp <- data.frame(pvalue = GOtest@result$pvalue,
                     term = GOtest@result$Description)
  
  results <- rbind.data.frame(results, temp)
  
}

# Save permutation results
save(results, file = "permResults_validation1.RData")

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

sum(table(results_all$permN[results_all$pvalue<0.05])>4)
