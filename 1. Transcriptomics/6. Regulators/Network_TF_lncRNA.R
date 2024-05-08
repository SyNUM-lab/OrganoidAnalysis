
# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(GENIE3)
library(tidyverse)
library(circlize)

# FUNCTION: Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
setwd("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/6. Regulators")


#==============================================================================#
# Perform network construction using GENIE3
#==============================================================================#

#..............................................................................#
# Transcription factors (proteins)
#..............................................................................#

# Load protein expression data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/2. Proteomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"pxData_imp.RData"))
pxMatrix_imp <- pxData_imp@assays@data@listData[[1]]
load(paste0(preprocessing_dir,"proteinAnnotation.RData"))
proteinAnnotation <- annotations

# Get transcription factors: http://humantfs.ccbr.utoronto.ca/download.php
TFs <- data.table::fread("Data/TFs.csv")
TFs <- TFs[TFs$`Is TF?`=="Yes"]

# Select expression values of TFs only
pxMatrix_fil <- pxMatrix_imp[rownames(pxMatrix_imp) %in% proteinAnnotation$ID[proteinAnnotation$hgnc_symbol %in% TFs$`HGNC symbol`],]

#..............................................................................#
# lncRNAs-regulators
#..............................................................................#

# Load lncRNA expression data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"gxMatrix_norm.RData"))
load(paste0(preprocessing_dir,"geneAnnotation.RData"))
load(paste0(preprocessing_dir,"Gene_biotype.RData"))

# Filter for lncRNA genes only
gxMatrix_lncRNA <- gxMatrix_norm[str_remove(rownames(gxMatrix_norm), "_.*") %in% annotations$ensembl_gene_id[annotations$gene_biotype == "lncRNA"],]

# Load lncRNA interaction data: Retrieved from NPInter V5
lncRNAint <- read.delim("Data/lncRNA_interaction.txt", header = FALSE)

# Select lncRNA-protein interactions
lncRNAint_fil <- lncRNAint[(lncRNAint$V7 == "protein") &
                             (lncRNAint$V11 == "Homo sapiens") &
                             (lncRNAint$V15 == "RNA-Protein"),
]

# Select lncRNAs that interact with TFs only
lncRNAint_fil <- lncRNAint_fil[lncRNAint_fil$V5 %in% TFs$`HGNC symbol`,]
lncRNAint_fil <- unique(lncRNAint_fil[,c(2,5)])
colnames(lncRNAint_fil) <- c("lncRNA", "Protein")

lncRNA_TF <- unique(lncRNAint_fil$lncRNA)

# Select lncRNA-DNA interactions
lncRNAint_fil <- lncRNAint[(lncRNAint$V11 == "Homo sapiens") &
                             (lncRNAint$V15 == "RNA-DNA"),
]

# Select lncRNAs that interact with TFs only
lncRNA_DNA <- unique(lncRNAint_fil$V2)

# get expression of these lncRNA regulators only
gxMatrix_lncRNA <- gxMatrix_lncRNA[rownames(gxMatrix_lncRNA) %in% 
                                     geneAnnotation$gene_id[geneAnnotation$GeneName 
                                                            %in% unique(c(lncRNA_DNA,lncRNA_TF))],]



#..............................................................................#
# Perform network interference
#..............................................................................#

# Eigengenes of GO terms
load(paste0("Data/eigengenes_all.RData"))
exprMatrix <- rbind(pxMatrix_fil,gxMatrix_lncRNA)
exprMatrix <- rbind(exprMatrix,eigengenes)

# Proteins are used as candidate regulators
regulators <- 1:(nrow(pxMatrix_fil)+nrow(gxMatrix_lncRNA))

# GO terms are used as candidate targets
targets <- ((nrow(pxMatrix_fil)+nrow(gxMatrix_lncRNA))+1):nrow(exprMatrix)

# Perform network construction
set.seed(123)
weightMat1 <- GENIE3(exprMatrix, regulators=regulators, targets = targets,
                     nTrees = 10000)

# Save results
save(weightMat1, file = "Data/weightMat_lncRNA_TF.RData")


#==============================================================================#
# Make plot
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(GENIE3)
library(tidyverse)
library(circlize)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Weight matrix
load("Data/weightMat_lncRNA_TF.RData")

# Protein annotation data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/2. Proteomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"proteinAnnotation.RData"))
proteinAnnotation <- annotations

# Gene annotation data
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"geneAnnotation.RData"))

# GO annotation data
load(paste0("D:/RTTproject/CellAnalysis/OrganoidAnalysis/GO_annotation/",
            "GOannotation.RData"))

# Selected GO terms
load(paste0("D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/5. GSEA/Data/",
            "terms_ordered1.RData"))


# Get regulatory interactions
linkList <- getLinkList(weightMat1)
linkList$regulatoryClass <- ifelse(str_detect(linkList$regulatoryGene, "ENSG"),
                                   "lncRNA", "TF")

# Top TF interactions
linkList_TF <- linkList[linkList$regulatoryClass == "TF",]
linkList_TF <- linkList_TF[linkList_TF$weight > quantile(linkList_TF$weight, 0.95),]
linkList_TF <- left_join(linkList_TF, unique(proteinAnnotation[,c(1,2,5)]),
                          by = c("regulatoryGene" = "ID"))

# Top lncRNA interactions
linkList_lncRNA <- linkList[linkList$regulatoryClass == "lncRNA",]
linkList_lncRNA <- linkList_lncRNA[linkList_lncRNA$weight > quantile(linkList_lncRNA$weight, 0.995),]
linkList_lncRNA <- left_join(linkList_lncRNA, unique(geneAnnotation[,c(1,2,3)]),
                         by = c("regulatoryGene" = "gene_id"))

# Combine TF and lncRNA interactions
linkList_fil <- data.frame(regulatoryGene = c(linkList_TF$hgnc_symbol,
                                              linkList_lncRNA$GeneName),
                           regulatoryID = c(linkList_TF$regulatoryGene,
                                            linkList_lncRNA$regulatoryGene),
                           targetTerm = c(linkList_TF$targetGene,
                                          linkList_lncRNA$targetGene),
                           Group = c(rep("TF", nrow(linkList_TF)),
                                     rep("lncRNA", nrow(linkList_lncRNA))))

linkList_fil <- left_join(linkList_fil, GOannotation,
                          by = c("targetTerm" = "Name"))
linkList_fil$Description <- firstup(linkList_fil$Description)
sel_IDs <- names(table(linkList_fil$regulatoryID))[table(linkList_fil$regulatoryID) > 2]
linkList_fil <- linkList_fil[linkList_fil$regulatoryID %in% sel_IDs,]

# Add ID to HGNC for proteins with multiple isoforms
isoGene <- unique(linkList_fil$regulatoryGene[which(duplicated(linkList_fil$regulatoryGene) & !duplicated(linkList_fil$regulatoryID))])

linkList_fil$regulatoryGene[linkList_fil$regulatoryGene == isoGene] <- paste0(linkList_fil$regulatoryGene[linkList_fil$regulatoryGene == isoGene], "\n(",
                                                                             linkList_fil$regulatoryID[linkList_fil$regulatoryGene == isoGene],")")

# Prepare data for plotting
plotData1 <- data.frame(source = as.character(linkList_fil$Description),
                        target = as.character(linkList_fil$regulatoryGene),
                        value = 1)


# Set order
rownames(GOannotation) <- GOannotation$Name
terms_ordered_descr <- firstup(GOannotation[terms_ordered,"Description"])
terms_ordered_descr[nchar(terms_ordered_descr)>20] <- paste0(substring(terms_ordered_descr[nchar(terms_ordered_descr)>20],1,17),"...")
plotData1$source[nchar(plotData1$source)>20] <- paste0(substring(plotData1$source[nchar(plotData1$source)>20],1,17),"...")
plotData1$order <- factor(plotData1$source,
                          levels = terms_ordered_descr)
plotData1 <- arrange(plotData1, by = order)
order <- factor(c(unique(plotData1$source), unique(plotData1$target)),
                levels = c(unique(plotData1$source), unique(plotData1$target)))

# make groups
groups <- c(rep("GO term", length(unique(plotData1$source))),
            ifelse(unique(plotData1$target) %in% unique(linkList_fil$regulatoryGene[linkList_fil$Group == "TF"]), "TF", "lncRNA")
            )
names(groups) <- order

# Set colors of groups
yellows <- grDevices::colorRampPalette(c("white","#DBA800"))
greens <- grDevices::colorRampPalette(c("white","darkgreen"))

colors <- c(RColorBrewer::brewer.pal(n = 6, "Blues"),
            greens(13)[2:13],
            yellows(9)[2:9],
            rep("darkgrey", length(unique(plotData1$target)))
)
names(colors) <- order

# Load proteomics statistics
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/2. Proteomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"DEresults_px.RData"))

# Get p-values
pvalues_px <- DEresults_px[,3:9]
rownames(pvalues_px) <- DEresults_px$name

# get adjusted p-values
adjpvalues_px <- pvalues_px
for (i in 1:ncol(pvalues_px)){
  adjpvalues_px[,i] <- p.adjust(pvalues_px[,i], "fdr")
}

# Get logFCs
logFCs_px <- DEresults_px[,str_detect(colnames(DEresults_px), "ratio")]
rownames(logFCs_px) <- DEresults_px$name

values_px <- -log10(pvalues_px)*sign(logFCs_px)

# Load transcriptomics statistics
preprocessing_dir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis/1. Transcriptomics/1. Preprocessing/"
load(paste0(preprocessing_dir,"DEresults_RTTvsIC_gx.RData"))
genes_all <- rownames(topList[[1]])

# Get p-values
pvalues_gx <- matrix(unlist(lapply(topList, function(x){x[genes_all,"PValue"]})), ncol = length(topList))
rownames(pvalues_gx) <- genes_all
colnames(pvalues_gx) <- names(topList)

# get adjusted p-values
adj_pvalues_gx <- matrix(unlist(lapply(topList, function(x){x[genes_all,"FDR"]})), ncol = length(topList))
rownames(adj_pvalues_gx) <- genes_all
colnames(adj_pvalues_gx) <- names(topList)

# Get logFCs
logFCs_gx <- matrix(unlist(lapply(topList, function(x){x[genes_all,"logFC"]})), ncol = length(topList))
rownames(logFCs_gx) <- genes_all
colnames(logFCs_gx) <- names(topList)


values_gx <- -log10(pvalues_gx)*sign(logFCs_gx)
colnames(values_gx) <- colnames(values_px)
values_all <- rbind(values_px, values_gx)
colnames(adj_pvalues_gx) <- colnames(adjpvalues_px)
adjpvalues_all <- rbind(adjpvalues_px, adj_pvalues_gx)

# Get regulatory interactions
values_fil <- values_all[unique(as.character(linkList_fil$regulatoryID)),]

# Prepare data for plotting
plotData <- gather(values_fil)
plotData$Name <- rep(rownames(values_fil),7)

# Combine with adjusted pvalues
adjpvalues_fil <- gather(adjpvalues_all[unique(as.character(linkList_fil$regulatoryID)),])
plotData <- cbind.data.frame(plotData, adjpvalues_fil$value)
plotData$Sig <- ifelse(plotData$`adjpvalues_fil$value` < 0.05, "black", "white")

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
                       levels = c("RTT_Dorsal_D75_vs_IC_Dorsal_D75_p.val",
                                  "RTT_Dorsal_D40_vs_IC_Dorsal_D40_p.val",
                                  "RTT_Dorsal_D13_vs_IC_Dorsal_D13_p.val",
                                  "RTT_Cell_D0_vs_IC_Cell_D0_p.val",
                                  "RTT_Ventral_D13_vs_IC_Ventral_D13_p.val",
                                  "RTT_Ventral_D40_vs_IC_Ventral_D40_p.val",
                                  "RTT_Ventral_D75_vs_IC_Ventral_D75_p.val"))

all_annotations <- data.frame(ID = c(proteinAnnotation$ID, geneAnnotation$gene_id),
                              hgnc_symbol = c(proteinAnnotation$hgnc_symbol, geneAnnotation$GeneName))


all_annotations <- data.frame(ID = linkList_fil$regulatoryID,
                              hgnc_symbol = linkList_fil$regulatoryGene)

plotData <- inner_join(plotData, unique(all_annotations),
                       by = c("Name" = "ID"))


# Set colors for statistics
colfunc <- grDevices::colorRampPalette(c("#000072", "white","red"))
vec_col <- colfunc(1000)
plotData$value1 <- log(abs(plotData$value))*sign(plotData$value)
test <- (round((plotData$value + max(abs(plotData$value))) * (999/max((plotData$value + max(abs(plotData$value))))))+1)
plotData$color <- vec_col[test]

tr_df <- data.frame(
  time = c("D75", "D40", "D13", "D0", "D13", "D40", "D75"),
  region = c(rep("Ventral",3), "iPSC", rep("Dorsal",3))
)


# Make plot
jpeg("TF_lncRNA_regulatoryNetwork.png", width = 8000, height = 8000, quality = 100)
circos.clear()
circos.par(start.degree = 257, 
           gap.degree = 4, 
           track.margin = c(-0.1, 0.1), 
           points.overflow.warning = FALSE,
           canvas.xlim = c(-1, 1), 
           canvas.ylim = c(-1, 1))
par(mar = c(0,0,0,0))

#grid.col = structure(c(rep(2, 5), rep(3, 5), rep(4, 5)), names = c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5)))

# Base plot
chordDiagram(
  x = plotData1, 
  order = order,
  group = groups,
  grid.col = colors,
  transparency = 0.25,
  directional = 0,
  direction.type = "diffHeight", 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05,0.05,0.05), #c(0.05, 0.1),
  #link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE,
  big.gap = 20,
  preAllocateTracks = list(
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(50)),
         track.margin = c(mm_h(3), 0)),
    list(track.height = c(mm_h(70)),
         track.margin = c(mm_h(3), 0))
  )
)


# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
    aa = c(1, 0.5)
    if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    #1.7
    circos.text(x=mean(xlim), y=-11.2, labels=sector.index, facing = dd, cex=7,  adj = aa)
    
  }
)
for (tr in 1:nrow(tr_df)){
  for (i in 1:length(unique(plotData1$target))){
    highlight.sector(unique(plotData1$target)[i], track.index = tr+1, 
                     col = plotData$color[(plotData$hgnc_symbol == unique(plotData1$target)[i]) &
                                            (plotData$Time == tr_df$time[tr]) &
                                            (plotData$Region == tr_df$region[tr])], 
                     border = plotData$Sig[(plotData$hgnc_symbol == unique(plotData1$target)[i]) &
                                           (plotData$Time == tr_df$time[tr]) &
                                           (plotData$Region == tr_df$region[tr])], #text
                     lwd = 3,
                     cex = 10, text.col = "black", niceFacing = TRUE)
  }
}

highlight.sector(names(groups)[groups == "lncRNA"], track.index = 1, 
                 col =  "#D9D9D9", border = "black",lwd = 2,
                 text = "lncRNA",
                 cex = 10, text.col = "black", niceFacing = TRUE)

highlight.sector(names(groups)[groups == "TF"], track.index = 1, 
                 col = "#D9D9D9", border = "black", lwd = 2,
                 text = "Protein",
                 cex = 10, text.col = "black", niceFacing = TRUE)

dev.off()
