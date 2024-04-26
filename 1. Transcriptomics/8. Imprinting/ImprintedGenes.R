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
library(readxl)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

setwd("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/")

# Load data
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_norm.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/DEresults_RTTvsIC_gx.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/geneAnnotation.RData")
load("D:/RTTproject/CellAnalysis/Data/SampleInfo.RData")

genes_all <- rownames(topList[[1]])

# Get p-values
pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues) <- genes_all
colnames(pvalues) <- names(topList)


# get adjusted p-values
adj_pvalues <- matrix(unlist(lapply(topList, function(x){x[genes_all,"FDR"]})), ncol = length(topList))
rownames(adj_pvalues) <- genes_all
colnames(adj_pvalues) <- names(topList)

# Get logFCs
logFCs <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
rownames(logFCs) <- genes_all
colnames(logFCs) <- names(topList)

all_genes <- unique(geneAnnotation$GeneName[geneAnnotation$gene_id %in% rownames(pvalues)])
all_genes <- all_genes[!is.na(all_genes)]

# Get imprinted genes: https://www.geneimprint.com/site/genes-by-species
imprintDF <- read_xlsx("D:/RTTproject/CellAnalysis/ImprintedGenes.xlsx")
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,3)],]
imprinted_genes <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[3])]
imprinted_genes <- imprinted_genes[!is.na(imprinted_genes)]
imprinted_genes <- intersect(imprinted_genes, all_genes)

# get non-imprinted genes
nimprinted_genes <- setdiff(all_genes, imprinted_genes)
nimprinted_genes <- nimprinted_genes[!is.na(nimprinted_genes)]

results_up <- matrix(NA,nrow = ncol(adj_pvalues), ncol = 4)
for (i in 1:ncol(adj_pvalues)){
  # get upregulated genes
  overall_up <- c(rownames(adj_pvalues)[(logFCs[,i] > 1) & (adj_pvalues[,i] < 0.05)])
  up_genes <- unique(geneAnnotation$GeneName[geneAnnotation$gene_id %in% overall_up])
  up_genes <- up_genes[!is.na(up_genes)]
  
  # get non-upregulated genes
  nup_genes <- setdiff(all_genes, up_genes)
  
  
  up_imp <- intersect(up_genes,imprinted_genes)
  nup_imp <- intersect(nup_genes,imprinted_genes)
  up_nimp <- intersect(up_genes, nimprinted_genes)
  nup_nimp <- intersect(nup_genes, nimprinted_genes)
  
  m <- matrix(c(length(up_imp), 
                length(nup_imp), 
                length(up_nimp),
                length(nup_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_up[i,] <- c(output$p.value,output$estimate,output$conf.int)
}





results_down <- matrix(NA,nrow = ncol(adj_pvalues), ncol = 4)
for (i in 1:ncol(adj_pvalues)){
  # get upregulated genes
  overall_down <- c(rownames(adj_pvalues)[(logFCs[,i] < -1) & (adj_pvalues[,i] < 0.05)])
  down_genes <- unique(geneAnnotation$GeneName[geneAnnotation$gene_id %in% overall_down])
  down_genes <- down_genes[!is.na(down_genes)]
  
  # get non-upregulated genes
  ndown_genes <- setdiff(all_genes, down_genes)
  
  
  down_imp <- intersect(down_genes,imprinted_genes)
  ndown_imp <- intersect(ndown_genes,imprinted_genes)
  down_nimp <- intersect(down_genes, nimprinted_genes)
  ndown_nimp <- intersect(ndown_genes, nimprinted_genes)
  
  m <- matrix(c(length(down_imp), 
                length(ndown_imp), 
                length(down_nimp),
                length(ndown_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_down[i,] <- c(output$p.value,output$estimate,output$conf.int)
}



results_both <- matrix(NA,nrow = ncol(adj_pvalues), ncol = 4)
for (i in 1:ncol(adj_pvalues)){
  # get upregulated genes
  overall_both <- c(rownames(adj_pvalues)[(abs(logFCs[,i]) > 1) & (adj_pvalues[,i] < 0.05)])
  both_genes <- unique(geneAnnotation$GeneName[geneAnnotation$gene_id %in% overall_both])
  both_genes <- both_genes[!is.na(both_genes)]
  
  # get non-upregulated genes
  nboth_genes <- setdiff(all_genes, both_genes)
  
  
  both_imp <- intersect(both_genes,imprinted_genes)
  nboth_imp <- intersect(nboth_genes,imprinted_genes)
  both_nimp <- intersect(both_genes, nimprinted_genes)
  nboth_nimp <- intersect(nboth_genes, nimprinted_genes)
  
  m <- matrix(c(length(both_imp), 
                length(nboth_imp), 
                length(both_nimp),
                length(nboth_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_both[i,] <- c(output$p.value,output$estimate,output$conf.int)
}





plotDF <- rbind.data.frame(results_up,
                           results_down,
                           results_both)
colnames(plotDF) <- c("Pvalue", "OR", "Lower", "Upper")
plotDF$TimeRegion <- factor(rep(str_remove(colnames(pvalues), "_RTTvsIC"),3),
                            levels = c("Dorsal_D75", "Dorsal_D40", "Dorsal_D13", "Cell_D0",
                                       "Ventral_D13" , "Ventral_D40", "Ventral_D75"))
plotDF$Set <- factor(c(rep("Upregulated",7), rep("Downregulated",7), rep("Both",7)),
                     levels = c("Downregulated", "Upregulated", "Both"))

plotDF$Time <- str_remove(plotDF$TimeRegion,".*_")
plotDF$Region <- str_remove(plotDF$TimeRegion, "_.*")
plotDF$Region[plotDF$Region == "Cell"] <- "iPSC"

mainPlot <- ggplot() +
  geom_hline(yintercept = 1, color = "black") +
  geom_point(data = plotDF, aes(x = Set, y = OR, color = Set), size = 3) +
  geom_segment(data = plotDF, aes(x = Set, xend = Set,y = Lower, yend = Upper, color = Set),
               linewidth = 1) +
  facet_grid(cols = vars(TimeRegion)) +
  labs(color = NULL) +
  scale_color_manual(values = c("#2171B5","#CB181D","#525252"))+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(hjust = 1, vjust = 0.5,
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjust = 1,
                                    margin = margin(t = 0, r = 5, b = 0, l = 0),
                                    size = 12))

sideColorPlot_time <- ggplot(data = plotDF) +
  geom_tile(aes(x = Set, y = "label", fill = Time)) +
  geom_text(data = plotDF[plotDF$Set == "Upregulated",],
            aes(x = Set, y = "label", label = Time)) +
  facet_grid(.~TimeRegion, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# region
sideColorPlot_region <- ggplot(data = plotDF) +
  geom_tile(aes(x = Set, y = "label", fill = Region)) +
  geom_text(data = plotDF[(plotDF$Set == "Upregulated") & 
                            (plotDF$Time == "D40" | plotDF$Time == "D0"),],
            aes(x = Set, y = "label", label = Region)) +
  facet_grid(.~TimeRegion, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots into a single figure
finalPlot <- ggarrange(mainPlot,
                       sideColorPlot_time,
                       sideColorPlot_region,
                       heights = c(8,0.5,0.5),nrow = 3,ncol = 1,
                       align = "v",
                       common.legend = TRUE,
                       legend = "right")

finalPlot


ggsave(finalPlot, 
       file = "D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Imprinting/OR_imprinting.png", 
       width = 8, height = 0.86*4.5)





# Permuation
nPerm <- 10000
results_perm <- matrix(NA,nrow = nPerm, ncol = 4)
for (i in 1:nPerm){
  # get upregulated genes
  overall_both <- rownames(adj_pvalues)[sample(1:nrow(adj_pvalues),2000)]
  both_genes <- unique(geneAnnotation$GeneName[geneAnnotation$gene_id %in% overall_both])
  both_genes <- both_genes[!is.na(both_genes)]
  
  # get non-upregulated genes
  nboth_genes <- setdiff(all_genes, both_genes)
  
  
  both_imp <- intersect(both_genes,imprinted_genes)
  nboth_imp <- intersect(nboth_genes,imprinted_genes)
  both_nimp <- intersect(both_genes, nimprinted_genes)
  nboth_nimp <- intersect(nboth_genes, nimprinted_genes)
  
  m <- matrix(c(length(both_imp), 
                length(nboth_imp), 
                length(both_nimp),
                length(nboth_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_perm[i,] <- c(output$p.value,output$estimate,output$conf.int)
}




















overall_down <- c(rownames(pvalues)[rowSums((logFCs < 0) & (pvalues < 0.05)) >= 6])

overall_down <- c(rownames(adj_pvalues)[(logFCs[,5] < -1) & (adj_pvalues[,5] < 0.05)])

overall_down <- rownames(adj_pvalues)[sample(1:nrow(adj_pvalues),1000)]

# get upgulated genes
down_genes <- unique(geneAnnotation$GeneName[geneAnnotation$EnsemblID %in% overall_down])
down_genes <- down_genes[!is.na(down_genes)]

# get non-upregulated genes
ndown_genes <- unique(geneAnnotation$GeneName[!(geneAnnotation$EnsemblID %in% overall_down)])
ndown_genes <- nup_genes[!is.na(ndown_genes)]


down_imp <- intersect(down_genes,imprinted_genes)
ndown_imp <- intersect(ndown_genes,imprinted_genes)
down_nimp <- intersect(down_genes, nimprinted_genes)
ndown_nimp <- intersect(ndown_genes, nimprinted_genes)

m <- matrix(c(length(down_imp), 
              length(ndown_imp), 
              length(down_nimp),
              length(ndown_nimp)),nrow = 2)


output <- chisq.test(m)
output$expected
output$observed
output

# odds for an imprinted gene to be upregulated is 9 times the odds for a non-imprinted gene
output <- fisher.test(m)
output
