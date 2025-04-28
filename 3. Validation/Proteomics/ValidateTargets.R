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
library(patchwork)

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
load(paste0(homeDir,"/2. Proteomics/3. GetTargets/proteinSelection.RData"))

# Get p-values
pvalues <- DEresults_px[,3:6]
rownames(pvalues) <- DEresults_px$name

# Get logFCs
logFCs <- DEresults_px[,str_detect(colnames(DEresults_px), "ratio")]
rownames(logFCs) <- DEresults_px$name

# Put p values in data frame
test_p <- pvalues
test_p$ID <- rownames(pvalues)
test_p$Uniprot <- str_remove(rownames(pvalues), "-.*")
rownames(test_p) <- test_p$ID
colnames(test_p) <- c("D15", "D22", "D3", "D9")

# Put logFCs in data frame
test_logFC <- logFCs
test_logFC$ID <- rownames(logFCs)
test_logFC$Uniprot <- str_remove(rownames(logFCs), "-.*")
rownames(test_logFC) <- test_logFC$ID
colnames(test_logFC) <- c("D15", "D22", "D3", "D9")

# Prepare data for plotting
plotDF <- gather(test_p[selection$ID,1:4])
plotDF$ProteinID <- rep(selection$ID,4)
plotDF$logFC <- gather(test_logFC[selection$ID,1:4])$value
plotDF <- inner_join(plotDF, selection, by = c("ProteinID" = "ID"))

# Set order of proteins
protein_order <- rev(c("DCXR", "FLNC", "NEO1", "MAN2B1", "PTGR1", "TPD52", "LDAH",
                   "LDHA", "NAPRT", "CBS", "HAT1", "NXF1", "PCDH1", "MT1X",
                   "S100A11", "NEFM", "NEFL", "EPB41L2", "VGF", "SNCG",
                   "GBE1", "NTN1", "LY6H", "PSMD5", "GSTM1", "KLC2", "ACP1", "CRYZ"))
  
plotDF$hgnc_symbol <- factor(plotDF$hgnc_symbol, levels = protein_order)
plotDF <- arrange(plotDF, by = hgnc_symbol)
plotDF$Name <- factor(paste0(plotDF$hgnc_symbol, " (", plotDF$ProteinID, ")"),
                      levels = unique(paste0(plotDF$hgnc_symbol, " (", plotDF$ProteinID, ")")))

# Add significance (yes or no) to data frame
plotDF$Sig <- ifelse(plotDF$value < 0.05, "Yes", "No")

# Set order of time points
plotDF$key <- factor(plotDF$key, levels =  c("D3", "D9", "D15", "D22"))

# Main plot
mainPlot <- ggplot() +
  geom_bar(data = plotDF, aes(x = abs(logFC), y = Name, fill = logFC,  color = Sig),
           linewidth = 0.5, stat = "identity") +
  facet_grid(cols = vars(key), rows = vars(ProteinGroup),scales = "free", space = "free_y") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-3,3),trans = "pseudo_log", oob = scales::squish) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "logFC") +
  ylab("Protein Name\n") +
  theme_void() +
  theme(axis.text.y = element_text(size = 10, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())

# Colors of columns (time points)
colSideColorPlot_time <- ggplot(data = data.frame(Time = factor(c("D3", "D9", "D15", "D22"),
                                                                levels = c("D3", "D9", "D15", "D22")))) +
  geom_tile(aes(x = Time, y = "label", fill = Time)) +
  geom_text(aes(x = Time, y = "label", label = Time)) +
  facet_grid(cols = vars(Time),scales = "free", space = "free_y") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Oranges")

# Combine all plots
finalPlot <- mainPlot + colSideColorPlot_time + patchwork::plot_layout(heights = c(9,0.5), nrow = 2,ncol = 1)
fisher.test(matrix(c(6,2,5572,28),nrow = 2))

# Save final plot
ggsave(finalPlot, file = "Target_ProteomicsValidation.png", width = 5, height = 8)

# Get legend of plot
legendPlot <- ggplot() +
  geom_bar(data = plotDF, aes(x = abs(logFC), y = Name, fill = logFC,  color = Sig),
           linewidth = 0.5, stat = "identity") +
  facet_grid(cols = vars(key), rows = vars(ProteinGroup),scales = "free", space = "free_y") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-3,3),trans = "pseudo_log", oob = scales::squish) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "logFC") +
  ylab("Protein Name\n") +
  theme_void() +
  theme(axis.text.y = element_text(size = 10, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())

legend <- cowplot::get_legend(legendPlot)

# Save legend
ggsave(legend, file = "legend_ProteinTargets.png", width = 8, height = 8)
