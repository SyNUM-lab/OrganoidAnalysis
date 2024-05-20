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
library(ggrepel)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/7. Splicing"))

# Load sample information
load(paste0(homeDir,"/sampleInfo.RData"))
sampleInfo$SampleID2 <- paste0(sampleInfo$Group, "_",
                               sampleInfo$Tissue, "_",
                               sampleInfo$Time, "_",
                               sampleInfo$Replicate)

# Load proteomics data
load(paste0(homeDir,"/2. Proteomics/1. Preprocessing/pxData_imp.RData"))
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]

# Load transcriptomics data
preprocessing_dir <- paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/")
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))

# Load protein-gene associations
load(paste0(homeDir,"/2. Proteomics/1. Preprocessing/px_gx_annotations1.RData"))

# Load GO data
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_UNIPROT_Hs.RData"))
rownames(GOannotation) <- GOannotation$Name
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

# Format proteomics data
pxDF <- gather(as.data.frame(pxMatrix_imp))
pxDF$ProteinID <- rep(rownames(pxMatrix_imp), ncol(pxMatrix_imp))
pxDF <- inner_join(pxDF, sampleInfo, by  = c("key" = "SampleID2"))
pxDF <- inner_join(pxDF, px_gx_annotations, by = c("ProteinID" = "ID"))
pxDF$ProteinID_Sample <- paste0(pxDF$ProteinID, "_", pxDF$key)
pxDF <- unique(pxDF[,c("ProteinID_Sample", "value")])

# Format transcriptomics data
gxDF <- gather(as.data.frame(gxMatrix_norm))
gxDF$GeneID <- rep(rownames(gxMatrix_norm), ncol(gxMatrix_norm))
gxDF <- inner_join(gxDF, sampleInfo, by = c("key" = "SampleID"))
gxDF <- inner_join(gxDF, px_gx_annotations, by = c("GeneID"= "gene_id"))
gxDF$ProteinID_Sample <- paste0(gxDF$ID, "_", gxDF$SampleID2)
gxDF <- gxDF[!is.na(gxDF$ID),]
gxDF <- unique(gxDF[,c("ProteinID_Sample", "value", "uniprot_gn_id", "ID",
                       "Group", "Time", "Tissue", "Replicate")])

# Combine proteomics and transcriptomics
completeDF <- inner_join(pxDF, gxDF, by = c("ProteinID_Sample" = "ProteinID_Sample"))


resultDF <- NULL
for (t in 1:length(terms_ordered)){
  
  # Select genes from GO term
  GOterm <- GOannotation$ID[GOannotation$Name == terms_ordered[t]]
  selGenes <- GOgenes[[GOterm]]
  completeDF_fil <- completeDF[completeDF$uniprot_gn_id %in% selGenes,]
  
  # Fit linear model and collect stats
  p <- rep(NA, length(unique(completeDF_fil$ID)))
  for (i in 1:length(unique(completeDF_fil$ID))){
    testDF <- completeDF_fil[completeDF_fil$ID== unique(completeDF_fil$ID)[i],]
    p[i] <- summary(lm(value.x ~ value.y + Group, data = testDF))$coefficients[3,4]
  }
  tempDF <- data.frame(Term = terms_ordered[t],
                       P = p,
                       Protein = unique(completeDF_fil$ID))
  
  resultDF <- rbind.data.frame(resultDF, tempDF)
}


# Combine stats with gene/protein annotations
resultDF <- left_join(resultDF, unique(px_gx_annotations[,c(4,6)]), by = c("Protein" = "ID"))

# Combine with GO annotation
resultDF <- inner_join(resultDF, GOannotation, by = c("Term" = "Name"))
resultDF$Term<- factor(resultDF$Term, levels = terms_ordered)
resultDF <- arrange(resultDF, by = Term)
resultDF$Description <- factor(firstup(resultDF$Description),
                                    levels = firstup(unique(resultDF$Description)))
resultDF$Group <- NA
resultDF$Group[resultDF$Term %in% terms_ordered[1:6]] <- "  "
resultDF$Group[resultDF$Term %in% terms_ordered[7:18]] <- " "
resultDF$Group[resultDF$Term %in% terms_ordered[19:26]] <- ""

# Make plot
set.seed(123)
p <- ggplot(data = unique(resultDF[,-4])) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_boxplot(aes(x = Description, y = -log10(P))) +
  geom_point(data = resultDF[resultDF$Protein %in% c("P16152", "Q16555-1", "P62854", "P32969"),], 
             aes(x = Description, y = -log10(P)),
             color = "red") +
  ggrepel::geom_text_repel(data = resultDF[resultDF$Protein %in% c("P16152", "Q16555-1", "P62854", "P32969"),],
                           aes(x = Description, y = -log10(P), label = GeneName)) +
  ylab(expression(-log[10]~"P value")) +
  xlab(NULL) +
  coord_flip() +
  facet_grid(rows = vars(Group), scale = "free", space = "free") +
  theme_bw()


panel_colors <- c("#FEC44F","#74C476","#4292C6")

# convert to grob
gp <- ggplotGrob(p) # where p is the original ggplot object

# assign the first 4 right-side facet strips with blue fill
for(i in 1:3){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}
grid::grid.draw(gp)

# Save plot
ggsave(gp, file = "DifferentialProteinGeneAssociations.png",
       width = 6, height = 6)


#==============================================================================#
# Plot protein-gene associations
#==============================================================================#

# "P16152" CBR1 ENSG00000159228_CBR1
#"Q16555-1" DPYSL2 ENSG00000092964_DPYSL2
#"P62854" RPS26 ENSG00000197728_RPS26
#"P32969" RPL9 ENSG00000163682_RPL9
i = which(px_gx_annotations$uniprot_gn_id== "P32969")[1]
modelDF <- data.frame(gxExpr = gxMatrix_norm[px_gx_annotations$gene_id[i],],
                      pxExpr = pxMatrix_imp[px_gx_annotations$ID[i],],
                      Group = c(rep("IC",21), rep("RTT",21)))

lmRes <- lm(pxExpr ~ gxExpr + Group, data = modelDF)

lmSum <- summary(lmRes)
pg <- lmSum[["coefficients"]][3,4]
ei <- lmSum[["coefficients"]][1,1]
ee <- lmSum[["coefficients"]][2,1]
eg <- lmSum[["coefficients"]][3,1]

p <- ggplot(modelDF) +
  geom_abline(intercept = ei, 
              slope = ee,
              color ="#9ECAE1", linewidth = 1) +
  geom_abline(intercept = ei + eg, 
              slope = ee,
              color = "#FC9272", linewidth = 1) +
  geom_point(aes(x = gxExpr, y = pxExpr, color = Group)) +
  ggtitle(str_split(px_gx_annotations[i,3], "_")[[1]][2]) +
  xlab(expression("Gene expression ("~log[2]~"CPM )")) +
  ylab(expression("Protein expression ("~log[2]~"intensity )")) +
  scale_color_manual(values = c("#377EB8","#E41A1C")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        legend.position = "none",
        legend.title = element_blank())

p
ggsave(p, width = 5, height = 3.5,
       file = "RPL9.png")




#==============================================================================#
# Make splicing plots
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
# "P16152" CBR1 ENSG00000159228_CBR1
#"Q16555-1" DPYSL2 ENSG00000092964_DPYSL2
#"P62854" RPS26 ENSG00000197728_RPS26
#"P32969" RPL9 ENSG00000163682_RPL9
gene <- "ENSG00000163682_RPL9" 
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
ggsave(p_all, file = "RPL9_splicing.png", width = 8, height = 7)
