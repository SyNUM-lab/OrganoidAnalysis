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

# Load DAS genes
load("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing/DASgenes.RData")
load("D:/RTTproject/CellAnalysis/Genes/Preprocessing/gxMatrix_raw1.RData")

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
save(GOresults, file = "D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing/GOresults_splicing.RData")

load("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing/GOresults_splicing.RData")
load("D:/RTTproject/CellAnalysis/Genes/RTTvsIC/GOEnrichment/Data/terms_ordered1.RData")
GOresults_fil <- GOresults[GOresults$name %in% terms_ordered,]
GOresults_fil$name <- factor(GOresults_fil$name, levels = terms_ordered)
GOresults_fil <- arrange(GOresults_fil, by = name)
GOresults_fil$Description <- factor(firstup(GOresults_fil$Description),
                                    levels = firstup(GOresults_fil$Description))

GOresults_fil$Group <- NA
GOresults_fil$Group[1:6] <- "  "
GOresults_fil$Group[7:18] <- " "
GOresults_fil$Group[19:26] <- ""

sum(GOresults_fil$pvalue < 0.05)

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

ggsave(gp, file = "D:/RTTproject/CellAnalysis/Genes/RTTvsIC/Splicing/GOenrichment_splicing.png",
       width = 6, height = 6)