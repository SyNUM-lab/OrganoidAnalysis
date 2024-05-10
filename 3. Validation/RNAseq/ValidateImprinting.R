
# Clear workspace and console
rm(list = ls())
cat("\014")
gc()

# Load packages
library(tidyverse)
library(ggpubr)
library(readxl)
library(meta)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/RNAseq"))

#==============================================================================#
# Gather validation sets
#==============================================================================#

# Threshold for significance
pvalue_thres <- 0.05
logFC_thres <- 1
p_type <- "PValue"

# Make empty lists
up_genes <- list()
down_genes <- list()
both_genes <- list()
all_genes <- list()

# Load statistics
accessID <- "GSE123753"
load(file = paste0(accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

up_genes[[1]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[1]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[1]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[1]] <- rownames(output)

# Load statistics
accessID <- "GSE123753"
load(file = paste0(accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

up_genes[[2]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[2]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[2]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[2]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_iPSC.RData"))
output <- top.table_iPSC

up_genes[[3]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[3]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[3]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[3]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

up_genes[[4]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[4]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[4]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[4]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0(accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

up_genes[[5]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[5]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[5]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[5]] <- rownames(output)

# Load statistics
accessID <- "GSE117511"
load(file = paste0(accessID,"/top.table.RData"))
output <- top.table

up_genes[[6]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[6]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[6]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[6]] <- rownames(output)

# Load statistics
accessID <- "GSE128380"
load(file = paste0(accessID,"/top.table_CC.RData"))
output <- top.table_CC

up_genes[[7]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[7]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[7]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[7]] <- rownames(output)

# Load statistics
accessID <- "GSE128380"
load(file = paste0(accessID,"/top.table_TC.RData"))
output <- top.table_TC

up_genes[[8]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC > logFC_thres)]
down_genes[[8]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (output$logFC < -1*logFC_thres)]
both_genes[[8]] <- rownames(output)[(output[,p_type] < pvalue_thres) & (abs(output$logFC) > logFC_thres)]
all_genes[[8]] <- rownames(output)

#==============================================================================#
# Evaluate overrepresentation of imprinted genes
#==============================================================================#

# Get imprinted genes: https://www.geneimprint.com/site/genes-by-species
imprintDF <- read_xlsx(paste0(homeDir,"/1. Transcriptomics/8. Imprinting/ImprintedGenes.xlsx"))
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,3)],]

results_up <- matrix(NA,nrow = length(both_genes), ncol = 4)
results_down <- matrix(NA,nrow = length(both_genes), ncol = 4)
results_both <- matrix(NA,nrow = length(both_genes), ncol = 4)
for (d in 1:length(both_genes)){
  
  # get all genes
  all_genes1 <- all_genes[[d]]
  
  # Get imprinted genes
  imprinted_genes <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[3])]
  imprinted_genes <- imprinted_genes[!is.na(imprinted_genes)]
  imprinted_genes <- intersect(imprinted_genes, all_genes1)
  
  # get non-imprinted genes
  nimprinted_genes <- setdiff(all_genes1, imprinted_genes)
  nimprinted_genes <- nimprinted_genes[!is.na(nimprinted_genes)]
  
  
  # Evaluate upregulated genes:
  
  # get upregulated genes
  up_genes1 <- up_genes[[d]]
  
  # get non-upregulated genes
  nup_genes1 <- setdiff(all_genes1, up_genes1)
  
  
  up_imp <- intersect(up_genes1,imprinted_genes)
  nup_imp <- intersect(nup_genes1,imprinted_genes)
  up_nimp <- intersect(up_genes1, nimprinted_genes)
  nup_nimp <- intersect(nup_genes1, nimprinted_genes)
  
  m <- matrix(c(length(up_imp), 
                length(nup_imp), 
                length(up_nimp),
                length(nup_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_up[d,] <- c(output$p.value,output$estimate,output$conf.int)
  
  # Evaluate downregulated genes:
  
  # get downpregulated genes
  down_genes1 <- down_genes[[d]]
  
  # get non-downregulated genes
  ndown_genes1 <- setdiff(all_genes1, down_genes1)
  
  
  down_imp <- intersect(down_genes1,imprinted_genes)
  ndown_imp <- intersect(ndown_genes1,imprinted_genes)
  down_nimp <- intersect(down_genes1, nimprinted_genes)
  ndown_nimp <- intersect(ndown_genes1, nimprinted_genes)
  
  m <- matrix(c(length(down_imp), 
                length(ndown_imp), 
                length(down_nimp),
                length(ndown_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_down[d,] <- c(output$p.value,output$estimate,output$conf.int)
  
  
  # Evaluate up- and downregulated genes:
  
  # get DEd genes
  both_genes1 <- both_genes[[d]]
  
  # get non-DEd genes
  nboth_genes1 <- setdiff(all_genes1, both_genes1)
  
  
  both_imp <- intersect(both_genes1,imprinted_genes)
  nboth_imp <- intersect(nboth_genes1,imprinted_genes)
  both_nimp <- intersect(both_genes1, nimprinted_genes)
  nboth_nimp <- intersect(nboth_genes1, nimprinted_genes)
  
  m <- matrix(c(length(both_imp), 
                length(nboth_imp), 
                length(both_nimp),
                length(nboth_nimp)),nrow = 2)
  
  output <- fisher.test(m)
  results_both[d,] <- c(output$p.value,output$estimate,output$conf.int)
}

#==============================================================================#
# Make plot
#==============================================================================#

# Combine results into data frame
plotDF <- rbind.data.frame(results_up,
                           results_down,
                           results_both)
colnames(plotDF) <- c("Pvalue", "OR", "Lower", "Upper")
plotDF$Dataset <- factor(c(rep("GSE123753",2),
                           rep("GSE107399",3),
                           "GSE117511",
                           rep("GSE128380",2)),
                         levels = c("GSE107399","GSE123753","GSE117511","GSE128380"))

plotDF$Tissue <- factor(c("Neuron", "NPC",
                          "iPSC", "NPC","Neuron",
                          "Inter\nneuron",
                          "Cingulate\ncortex", "Temporal\ncortex"),
                        levels = c("iPSC", "NPC", "Neuron", 
                                   "Inter\nneuron",
                                   "Temporal\ncortex", "Cingulate\ncortex"))

plotDF$DatasetTissue <- factor(paste0(plotDF$Tissue, "_",
                                      plotDF$Dataset),
                               levels = paste0(c("iPSC", "NPC", "Neuron", 
                                                 "NPC", "Neuron",
                                                 "Inter\nneuron",
                                                 "Temporal\ncortex", "Cingulate\ncortex"), 
                                               "_",
                                               c(rep("GSE107399",3),
                                                 rep("GSE123753",2),
                                                 "GSE117511",
                                                 rep("GSE128380",2))))


plotDF$Set <- factor(c(rep("Upregulated",8), rep("Downregulated",8), rep("Both",8)),
                     levels = c("Downregulated", "Upregulated", "Both"))

# make plot
mainPlot <- ggplot() +
  geom_hline(yintercept = 1, color = "black") +
  geom_point(data = plotDF, aes(x = Set, y = OR, color = Set), size = 3) +
  geom_segment(data = plotDF, aes(x = Set, xend = Set,y = Lower, yend = Upper, color = Set),
               linewidth = 1) +
  facet_grid(cols = vars(DatasetTissue)) +
  labs(color = NULL) +
  scale_color_manual(values = c("#2171B5","#CB181D","#525252"))+
  theme_bw() +
  ylim(0,9.5) +
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

colors <- c("#E7298A","#66A61E","#E6AB02","#A6761D","#666666")
sideColorPlot_dataset <- ggplot(data = plotDF) +
  geom_tile(aes(x = Set, y = "label", fill = Dataset)) +
  facet_grid(.~DatasetTissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = colors)


sideColorPlot_tissue <- ggplot(data = plotDF) +
  geom_tile(aes(x = Set, y = "label"), fill = "white") +
  geom_text(data = plotDF[plotDF$Set == "Upregulated",], 
            aes(x = Set, y = "label", label = Tissue), angle = 90,hjust = 0.5) +
  facet_grid(.~DatasetTissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())


# Combine plots into a single figure
finalPlot <- ggarrange(mainPlot,
                       sideColorPlot_dataset,
                       sideColorPlot_tissue,
                       heights = c(8,0.5,2),nrow = 3,ncol = 1,
                       align = "v",
                       common.legend = TRUE,
                       legend = "right")

finalPlot


ggsave(finalPlot, 
       file = "OR_imprinting_validation.png", 
       width = 8, height = 4.5)

#==============================================================================#
# Perform meta analysis
#==============================================================================#

# Meta analysis for downregulated genes
m.gen <- metagen(TE = log(results_down[,2]),
                 lower = log(results_down[,3]),
                 upper = log(results_down[,4]),
                 studlab = paste0("test",1:8),
                 sm = "OR",
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 hakn = TRUE)

summary(m.gen)

# Meta analysis for upregulated genes
m.gen <- metagen(TE = log(results_up[,2]),
                 lower = log(results_up[,3]),
                 upper = log(results_up[,4]),
                 studlab = paste0("test",1:8),
                 sm = "OR",
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 hakn = TRUE)

summary(m.gen)

# Meta analysis for all differentially expressed genes
m.gen <- metagen(TE = log(results_both[,2]),
                 lower = log(results_both[,3]),
                 upper = log(results_both[,4]),
                 studlab = paste0("test",1:8),
                 sm = "OR",
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 hakn = TRUE)

summary(m.gen)

# Combine results from meta analysis in data frame
plotDF <- rbind.data.frame(c(0.0003, 2.0587, 1.5959, 2.6556),
                           c(0.0046, 2.3016, 1.4236, 3.7213),
                           c(0.0070, 1.9864,1.2915, 3.0552))
colnames(plotDF) <- c("Pvalue", "OR", "Lower", "Upper")
plotDF$Set <- factor(c("Upregulated", "Downregulated", "Both"),
                     levels = c("Downregulated", "Upregulated", "Both"))

# Make plot
mainPlot <- ggplot() +
  geom_hline(yintercept = 1, color = "black") +
  geom_point(data = plotDF, aes(x = Set, y = OR, color = Set), size = 3) +
  geom_segment(data = plotDF, aes(x = Set, xend = Set,y = Lower, yend = Upper, color = Set),
               linewidth = 1) +
  labs(color = NULL) +
  scale_color_manual(values = c("#2171B5","#CB181D","#525252"))+
  ylim(0,9.5) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = "#F0F0F0"),
        legend.text = element_text(size = 12),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"),
        axis.text.y = element_text(hjust = 1, vjust = 0.5,
                                   margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(angle = 90, vjust = 1,
                                    margin = margin(t = 0, r = 5, b = 0, l = 0),
                                    size = 12))

# Save plot
ggsave(mainPlot, 
       file = "OR_imprinting_validation_pooled.png", 
       width = 3.2, height = 4.5)
