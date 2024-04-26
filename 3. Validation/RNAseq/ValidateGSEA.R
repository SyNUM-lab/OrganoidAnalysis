# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(biomaRt)

# Set working directory
setwd("D:/RTTproject/CellAnalysis")

#==============================================================================#
# GSE123753: Neuron
#==============================================================================#

# Load statistics
accessID <- "GSE123753"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/gxMatrix_raw.RData"))
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

# Get ENSEMBL IDs
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
annotations <- getBM(attributes=c("entrezgene_id",
                                  "ensembl_gene_id",
                                  "hgnc_symbol"), 
                     filters = 'hgnc_symbol',
                     values = rownames(gxMatrix),
                     mart = ensembl)
save(annotations, file = "Genes/Validation/HGNC2ENSEMBL.RData")

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))


# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE123753_Neuron.RData")


#==============================================================================#
# GSE123753: NPC
#==============================================================================#

# Load statistics
accessID <- "GSE123753"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))

# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE123753_NPC.RData")


#==============================================================================#
# GSE107399: iPSC
#==============================================================================#

# Load statistics
accessID <- "GSE107399"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_iPSC.RData"))
output <- top.table_iPSC

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))


# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE107399_iPSC.RData")

#==============================================================================#
# GSE107399: NPC
#==============================================================================#

# Load statistics
accessID <- "GSE107399"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))

# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE107399_NPC.RData")

#==============================================================================#
# GSE107399: Neuron
#==============================================================================#

# Load statistics
accessID <- "GSE107399"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))


# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE107399_Neuron.RData")

#==============================================================================#
# GSE117511: Organoids at multiple time points
#==============================================================================#

# Load statistics
accessID <- "GSE117511"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table.RData"))
output <- top.table

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))


# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE117511.RData")

#==============================================================================#
# GSE128380: Cingulate cortex (CC)
#==============================================================================#

# Load statistics
accessID <- "GSE128380"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_CC.RData"))
output <- top.table_CC

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))


# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE128380_CC.RData")

#==============================================================================#
# GSE128380: Temporal cortex (TC)
#==============================================================================#

# Load statistics
accessID <- "GSE128380"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_TC.RData"))
output <- top.table_TC

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])
sum(duplicated(annotations$ensembl_gene_id))


# Prepare data for GSEA
gsea_input <- data.frame(value = -log(output$PValue)*sign(output$logFC),
                         key = rownames(output))
gsea_input <- inner_join(gsea_input, annotations, by = c("key" = "hgnc_symbol"))

# remove duplicated ensembl ids (remove with lowest significance)
gsea_input <- arrange(gsea_input, by = desc(abs(value)))
gsea_input <- gsea_input[!duplicated(gsea_input$ensembl_gene_id),]

# Make named vector for input
gsea_vector <- setNames(gsea_input$value, gsea_input$ensembl_gene_id)
gsea_vector <- sort(gsea_vector, decreasing = TRUE)

# Load GO annotation
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform GSEA
set.seed(123)
GOtest <- GSEA(
  geneList = gsea_vector,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 10000,
  eps = 1e-10,
  pvalueCutoff = 1.5,
  pAdjustMethod = "fdr",
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
)

# Get results
results <- GOtest@result

# Save output
save(GOtest,file =  "Genes/Validation/GOtest_GSE128380_TC.RData")

###############################################################################

# Make plot

###############################################################################
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggpubr)


# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
setwd("D:/RTTproject/CellAnalysis")

# Load GO terms
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")

rownames(GOannotation) <- GOannotation$Name
terms_ordered_descr <- GOannotation[terms_ordered, "Description"]

# Load results
statList <- list()

load("Genes/Validation/GOtest_GSE107399_NPC.RData")
statList[["GSE107399_NPC"]] <- GOtest@result

load("Genes/Validation/GOtest_GSE107399_iPSC.RData")
statList[["GSE107399_iPSC"]] <- GOtest@result

load("Genes/Validation/GOtest_GSE107399_Neuron.RData")
statList[["GSE107399_Neuron"]] <- GOtest@result

load("Genes/Validation/GOtest_GSE123753_NPC.RData")
statList[["GSE123753_NPC"]] <- GOtest@result

load("Genes/Validation/GOtest_GSE123753_Neuron.RData")
statList[["GSE123753_Neuron"]] <- GOtest@result

load("Genes/Validation/GOtest_GSE117511.RData")
statList[["GSE117511"]] <- GOtest@result

load("Genes/Validation/GOtest_GSE128380_TC.RData")
statList[["GSE128380_TC"]] <- GOtest@result

load("Genes/Validation/GOtest_GSE128380_CC.RData")
statList[["GSE128380_CC"]] <- GOtest@result


plotDF <- NULL
for (i in 1:length(statList)){
  temp <- data.frame(Description = statList[[i]]$Description,
                     Value = -log10(statList[[i]]$pvalue)*sign(statList[[i]]$NES),
                     Sig = ifelse(statList[[i]]$pvalue < 0.05,"Yes", "No"),
                     Set = names(statList)[i]
  )
  
  
  plotDF <- rbind.data.frame(plotDF, temp)
}

plotDF$Description <- factor(firstup(plotDF$Description), 
                             levels = firstup(terms_ordered_descr))

setInfo <- data.frame(Set = c("GSE107399_NPC", "GSE107399_iPSC", "GSE107399_Neuron",
                              "GSE123753_NPC", "GSE123753_Neuron", "GSE117511", 
                              "GSE128380_TC", "GSE128380_CC"),
                      Study = c("GSE107399", "GSE107399", "GSE107399",
                                "GSE123753", "GSE123753", "GSE117511", 
                                "GSE128380", "GSE128380"),
                      Source = c("NPC", "iPSC", "Neuron",
                                 "NPC", "Neuron", "Inter-\nneuron", 
                                 "Temporal\ncortex", 
                                 "Cingulate\ncortex")
                      )

plotDF <- inner_join(plotDF, setInfo, by = c("Set" = "Set"))
plotDF$Source <- factor(plotDF$Source,
                        levels = c("iPSC", "NPC", "Neuron",
                                   "Inter-\nneuron", 
                                   "Temporal\ncortex", 
                                   "Cingulate\ncortex"))

setInfo$Source <- factor(setInfo$Source,
                        levels = c("iPSC", "NPC", "Neuron",
                                   "Inter-\nneuron", 
                                   "Temporal\ncortex", 
                                   "Cingulate\ncortex"))

plotDF$Study <- factor(plotDF$Study,
                       levels = c("GSE107399", "GSE123753", "GSE117511", 
                                  "GSE128380"))

setInfo$Study<- factor(setInfo$Study,
                       levels = c("GSE107399", "GSE123753", "GSE117511", 
                                  "GSE128380"))

plotDF <- plotDF[plotDF$Description %in% firstup(terms_ordered_descr),]

colors <- c("#E7298A","#66A61E","#E6AB02","#A6761D","#666666")
mainPlot <- ggplot(data = plotDF, aes(x = Source, y = Description, fill = Value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(cols = vars(Study), scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "Signed\n-log10 p-value") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Time
colSideColorPlot <- ggplot(data = setInfo) +
  geom_tile(aes(x = Source, y = "label", fill = Study)) +
  facet_grid(.~Study, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = colors)


# Combine plots into a single figure
finalPlot <- ggarrange(colSideColorPlot,
                       mainPlot,
                       heights = c(0.5,9),nrow = 2,ncol = 1,
                       align = "v",
                       common.legend = FALSE)

ggsave(finalPlot, file = "Genes/Validation/GOEnrichment_RTTvsIC1.png", 
       width = 8, height = 8)




# Get legends:
legendPlot <- ggplot(data = setInfo) +
  geom_tile(aes(x = Source, y = "label", fill = Study)) +
  facet_grid(.~Study, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20)) +
  scale_fill_manual(values = colors)

legend <- cowplot::get_legend(legendPlot)

ggsave(legend, file = "Genes/Validation/legend1.png", width = 8, height = 8)

# Get legends:
legendPlot <- ggplot(data = plotDF, aes(x = Source, y = Description, fill = Value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(cols = vars(Study), scale = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "Signed\n-log10 p-value") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=20),
        legend.title = element_blank())

legend <- cowplot::get_legend(legendPlot)

ggsave(legend, file = "Genes/Validation/legend2.png", width = 8, height = 8)

#==============================================================================#
# Perform permutation analysis to validate findings
#==============================================================================#

all_genes <- list()
# Load statistics
accessID <- "GSE123753"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

all_genes[[1]] <- rownames(output)

# Load statistics
accessID <- "GSE123753"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

all_genes[[2]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_iPSC.RData"))
output <- top.table_iPSC

all_genes[[3]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_NPC.RData"))
output <- top.table_NPC

all_genes[[4]] <- rownames(output)

# Load statistics
accessID <- "GSE107399"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_Neuron.RData"))
output <- top.table_Neuron

all_genes[[5]] <- rownames(output)

# Load statistics
accessID <- "GSE117511"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table.RData"))
output <- top.table

all_genes[[6]] <- rownames(output)

# Load statistics
accessID <- "GSE128380"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_CC.RData"))
output <- top.table_CC

all_genes[[7]] <- rownames(output)

# Load statistics
accessID <- "GSE128380"
load(file = paste0("D:/RTTproject/PublicDatasets/", accessID,"/top.table_TC.RData"))
output <- top.table_TC

all_genes[[8]] <- rownames(output)

nGenes <- round(mean(unlist(lapply(all_genes, length))))
unionGenes <- unique(unlist(all_genes))

# Load annotation
load("Genes/Validation/HGNC2ENSEMBL.RData")
annotations <- unique(annotations[!is.na(annotations$ensembl_gene_id),2:3])

gene_vector  <- annotations$ensembl_gene_id[annotations$hgnc_symbol %in% unique(unlist(all_genes))]
gene_vector <- unique(gene_vector[!is.na(gene_vector)])


# Load GO annotation
setwd("D:/RTTproject/CellAnalysis")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOannotation.RData")
load("D:/RTTproject/CellAnalysis/GO_annotation/GOgenes_ENSEMBL.RData")
load("Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

results <- NULL
set.seed(345)
nPerm <- 1000
for (i in 1:nPerm){
  
  # select random genes from the union of the genes that pass QC
  # The length of the selected genes should be equal to the mean number of
  # genes that pass QC in the different datasets
  nSel <- round(nGenes*(length(gene_vector)/length(unionGenes)))
  
  # Put the selected genes in a random order
  genes <- gene_vector[sample(1:length(gene_vector),nSel)]
  gsea_vector <- setNames(rev(1:length(genes)),genes)
  
  # Perform GSEA
  GOtest <- GSEA(
    geneList = gsea_vector,
    exponent = 1,
    minGSSize = 0,
    maxGSSize = 10000,
    eps = 1e-10,
    pvalueCutoff = 1.5,
    pAdjustMethod = "fdr",
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    verbose = FALSE
  )
  
  # Get results
  temp <- data.frame(pvalue = GOtest@result$pvalue,
                     term = GOtest@result$Description)
  
  results <- rbind.data.frame(results, temp)
  
}

save(results, file = "Genes/Validation/permResults_validation1.RData")
hist(results$pvalue)

results_fil <- results[results$term != "external encapsulating structure organization",]

results$permN <- rep(1:1000,each = 26)

sum(table(results$permN[results$pvalue<0.05])>5)
