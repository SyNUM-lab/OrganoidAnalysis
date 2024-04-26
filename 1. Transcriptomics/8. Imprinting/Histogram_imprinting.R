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

imprinted_genes_id <- unique(geneAnnotation$gene_id[geneAnnotation$GeneName %in% imprinted_genes])
imprinted_genes_id <- intersect(imprinted_genes_id, rownames(pvalues))
nimprinted_genes_id <- unique(geneAnnotation$gene_id[geneAnnotation$GeneName %in% nimprinted_genes])
nimprinted_genes_id <- intersect(nimprinted_genes_id, rownames(pvalues))

nSig <- rowSums((abs(logFCs) > 1)&(adj_pvalues<0.05))

nSig_imprinted <- nSig[names(nSig) %in% imprinted_genes_id]

nSig_nimprinted <- nSig[names(nSig) %in% nimprinted_genes_id]

names(nSig_imprinted)[nSig_imprinted == 6]
names(nSig_imprinted)[nSig_imprinted == 7]

table(nSig_imprinted)/length(nSig_imprinted)
table(nSig_nimprinted)/length(nSig_nimprinted)

plotDF <- data.frame(value_rel = c(table(nSig_imprinted)/length(nSig_imprinted),
                                   table(nSig_nimprinted)/length(nSig_nimprinted)),
                     value_abs = c(table(nSig_imprinted),
                                   table(nSig_nimprinted)),
                     nSig_num = as.numeric(c(names(table(nSig_imprinted)), 
                              names(table(nSig_nimprinted)))),
                     nSig_char = c(names(table(nSig_imprinted)), 
                                         names(table(nSig_nimprinted))),
                     Group = factor(c(rep("Imprinted",8),
                               rep("Non-imprinted",8)), 
                               levels = c("Non-imprinted", "Imprinted"))
                     )

plotDF$text_pos <- plotDF$nSig_num + ifelse(plotDF$Group == "Imprinted",0.2,-0.2)

colors <- rev(c("#EF4040", "#FFA732"))

p <- ggplot(plotDF) +
  geom_bar(aes(y = value_rel, x = nSig_num, fill = Group),
           stat = "identity", position = position_dodge2(padding = 0.1)) +
  geom_text(aes(y = value_rel + 0.01, x = text_pos, label = value_abs), 
            angle = 90, vjust = 0.5, hjust = 0) +
  xlab("# Time points/brain regions") +
  ylab("# DEGs") +
  ylim(0,0.8) +
  labs(fill = NULL) +
  scale_x_continuous(breaks = 0:7) +
  scale_fill_manual(values = colors) +
  theme_void() +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(linewidth = 0.5, lineend = "butt"),
        axis.text.x = element_text(),
        axis.title.x = element_text(size = 12,
                                    margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(linewidth = 0.5, lineend = "butt"),
        axis.text.y = element_text(hjust = 1, vjust = 0.5,
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjust = 1,
                                    margin = margin(t = 0, r = 5, b = 0, l = 0),
                                    size = 12),
        axis.ticks.length  = unit(0.1, "cm"),
        plot.margin = margin(1,1,1,1, "cm"),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 12))

ggsave(p, 
       file = "D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Imprinting/histogram_imprinting.png",
       width = 6, height = 4.5)

plotDF_pie <- data.frame(value = c(length(unique(nimprinted_genes_id)),
                                   length(unique(imprinted_genes_id))),
                         Group = factor(c("Non-imprinted", "Imprinted"),
                                        levels = c("Non-imprinted", "Imprinted")))

pie <- ggplot(plotDF_pie, aes(x="", y=value, fill=Group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = colors) +
  theme_void() +
  theme(legend.position = "none")
  
ggsave(pie, 
       file = "D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Imprinting/piechart_imprinting.png",
       width = 6, height = 4.5)