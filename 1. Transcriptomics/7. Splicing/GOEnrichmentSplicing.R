# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/7. Splicing"))

# Load DAS genes
load("Data/DASgenes.RData")
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_raw1.RData"))

# Perform GO ORA
gene <- unlist(lapply(str_split(DASgenes, "_"), function(x)x[[1]]))
universe <- unlist(lapply(str_split(rownames(gxMatrix_raw), "_"), function(x)x[[1]]))

GOtest <- enrichGO(
  gene = gene,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pvalueCutoff = Inf,
  pAdjustMethod = "fdr",
  universe = universe,
  qvalueCutoff = Inf,
  minGSSize = 0,
  maxGSSize = Inf,
  readable = FALSE,
  pool = FALSE
)

GOresults <- GOtest@result
GOresults$name <- paste0(GOresults$Description, " (", GOresults$ID, ")")

# Save results
save(GOresults, file = "Data/GOresults_splicing.RData")


#==============================================================================#
# Make plot
#==============================================================================#

# Load data
load("Data/GOresults_splicing.RData")
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

# Prepare data for plotting
GOresults_fil <- GOresults[GOresults$name %in% terms_ordered,]
GOresults_fil$name <- factor(GOresults_fil$name, levels = terms_ordered)
GOresults_fil <- arrange(GOresults_fil, by = name)
GOresults_fil$Description <- factor(firstup(GOresults_fil$Description),
                                    levels = firstup(GOresults_fil$Description))
GOresults_fil$Group <- NA
GOresults_fil$Group[1:6] <- "  "
GOresults_fil$Group[7:18] <- " "
GOresults_fil$Group[19:26] <- ""

# How many GO terms reach significance?
sum(GOresults_fil$pvalue < 0.05)

# Make plot
p <- ggplot(GOresults_fil) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_bar(aes(x = Description, y = -log10(pvalue)),
           stat = "identity", position = position_dodge()) +
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
ggsave(gp, file = "GOenrichment_splicing.png",
       width = 6, height = 6)

#==============================================================================#
# Perform permutation
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/7. Splicing"))

# Load DAS genes
load("Data/DASgenes.RData")
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/gxMatrix_raw1.RData"))

# Load GO annotation
load(paste0(homeDir,"/GO_annotation/GOannotation.RData"))
load(paste0(homeDir,"/GO_annotation/GOgenes_BP_ENSEMBL_Hs.RData"))
load(paste0(homeDir,"/1. Transcriptomics/5. GSEA/Data/terms_ordered1.RData"))

# Prepare TERM2GENE and TERM2NAME objects
GOgenes_fil <- GOgenes[GOannotation$ID[GOannotation$Name %in% terms_ordered]]
TERM2GENE <- data.frame(TERM = unlist(lapply(seq_along(GOgenes_fil), function(x){rep(names(GOgenes_fil)[x],length(GOgenes_fil[[x]]))})),
                        GENE = unlist(GOgenes_fil))

TERM2NAME <- GOannotation[GOannotation$Name %in% terms_ordered, c(2,3)]
colnames(TERM2NAME) <- c("TERM", "NAME")

# Perform permutation
nPerm <- 1000
all_genes <- rownames(gxMatrix_raw)
nGenes <- length(DASgenes)
permResults <- rep(NA, nPerm)
set.seed(123)
for (p in 1:nPerm){
  
  permGenes <- all_genes[sample(1:length(all_genes),nGenes)]
  gene <- unlist(lapply(str_split(permGenes, "_"), function(x)x[[1]]))
  universe <- unlist(lapply(str_split(all_genes, "_"), function(x)x[[1]]))
  
  GOtest <- enricher(
    gene = gene,
    pvalueCutoff = Inf,
    pAdjustMethod = "fdr",
    universe = universe,
    qvalueCutoff = Inf,
    minGSSize = 0,
    maxGSSize = Inf,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
  )
  
  GOresults <- GOtest@result
  permResults[p] <- sum(GOresults$pvalue < 0.05)
}

hist(permResults)
sum(permResults >= 15)


