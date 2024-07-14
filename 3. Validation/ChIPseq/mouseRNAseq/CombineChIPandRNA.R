# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(GEOquery)
library(rtracklayer)
library(tidyverse)
library(biomaRt)
library(Gviz)
library(readxl)
library(chipseq)
library(BSgenome.Mmusculus.UCSC.mm9)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/ChIPseq"))

# Binning function
BinChIPseq = function(reads, bins){
  mcols(bins)$score = countOverlaps(bins, reads) 
  return(bins) 
}

# Get sequence info
genome = BSgenome.Mmusculus.UCSC.mm9
si = seqinfo(genome)
si = si[ paste0('chr', c(1:19))]

# Import BW files (GSE71126, GEO)
MeCP2_rep1 <- import.bw("Data/GSE71126_RAW/GSM1827604_MeCP2_ChIP_WT_rep1.extend200.coverage.mm9.bw")
MeCP2_rep2 <- import.bw("Data/GSE71126_RAW/GSM1827605_MeCP2_ChIP_WT_rep2.extend200.coverage.mm9.bw")
Input_rep1 <- import.bw("Data/GSE71126_RAW/GSM1827606_Input_WT_rep1.extend200.coverage.mm9.bw")
Input_rep2 <- import.bw("Data/GSE71126_RAW/GSM1827607_Input_WT_rep2.extend200.coverage.mm9.bw")

# Select chr1 - chr19 (tremove sex and mitochrondrial genes)
MeCP2_rep1_fil <- MeCP2_rep1[!(seqnames(MeCP2_rep1) %in% c("chrM", "chrY", "chrX"))]
rm(MeCP2_rep1)
MeCP2_rep2_fil <- MeCP2_rep2[!(seqnames(MeCP2_rep2) %in% c("chrM", "chrY", "chrX"))]
rm(MeCP2_rep2)
Input_rep1_fil <- Input_rep1[!(seqnames(Input_rep1) %in% c("chrM", "chrY", "chrX"))]
rm(Input_rep1)
Input_rep2_fil <- Input_rep2[!(seqnames(Input_rep2) %in% c("chrM", "chrY", "chrX"))]
rm(Input_rep2)


#******************************************************************************#
# Get imprinted genes
#******************************************************************************#

# load data (if needed)
load("Data/egs.RData")

# All genes in data
all_genes <- egs$external_gene_id

# Imprinted genes in data (https://www.geneimprint.com/site/genes-by-species.Mus+musculus)
imprintDF <- read_xlsx("Data/imprinted_genes_mouse.xlsx")
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,4)],]
imprinted_genes <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[1])]
imprinted_genes <- imprinted_genes[!is.na(imprinted_genes)]
imprinted_genes <- intersect(imprinted_genes, all_genes)

# Imprinted vs non-imprinted genes
egs_imprinted <- egs[egs$external_gene_id %in% imprinted_genes,]
egs_nonimprinted <- egs[!(egs$external_gene_id %in% imprinted_genes),]

# Settings:
range <- 3000 # Range around TSS
lo <- 60      # length out (bin size)

#==============================================================================#
# 1) Get bin counts around TSS of imprinted genes
#==============================================================================#

# Set tiles around TSS of imprinted genes
tiles_imprinted = sapply(1:nrow(egs_imprinted), function(i)
  if(egs_imprinted$strand[i] == "1" )
    egs_imprinted$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_imprinted$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_imprinted = GRanges(tilename = paste(rep(egs_imprinted$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                          seqnames = Rle(rep(paste0('chr', egs_imprinted$chromosome_name), each=lo)), 
                          ranges = IRanges(start = as.vector(tiles_imprinted),
                                           width = 100),
                          strand = Rle(rep("*", length(as.vector(tiles_imprinted)))),
                          seqinfo=si)

# Get bin counts
imprintedLoci = countOverlaps(tiles_imprinted, MeCP2_rep1_fil) + 
  countOverlaps(tiles_imprinted, MeCP2_rep2_fil)

imprintedLoci_matrix = matrix(imprintedLoci, nrow=nrow(egs_imprinted), 
                              ncol=lo, byrow=TRUE)


# Set tiles around TSS of non-imprinted genes
tiles_nonimprinted = sapply(1:nrow(egs_nonimprinted), function(i)
  if(egs_nonimprinted$strand[i] == "1" )
    egs_nonimprinted$TSS[i] + seq(-1*range, range-100, length.out=lo)
  else
    egs_nonimprinted$TSS[i] + seq(range-100, -1*range, length.out=lo))

tiles_nonimprinted = GRanges(tilename = paste(rep(egs_nonimprinted$ensembl_gene_id, each=lo), 1:lo, sep="_"),
                             seqnames = Rle(rep(paste0('chr', egs_nonimprinted$chromosome_name), each=lo)), 
                             ranges = IRanges(start = as.vector(tiles_nonimprinted),
                                              width = 100),
                             strand = Rle(rep("*", length(as.vector(tiles_nonimprinted)))),
                             seqinfo=si)

# Get bin counts
nonimprintedLoci = countOverlaps(tiles_nonimprinted, MeCP2_rep1_fil) + 
  countOverlaps(tiles_nonimprinted, MeCP2_rep2_fil)

nonimprintedLoci_matrix = matrix(nonimprintedLoci, nrow=nrow(egs_nonimprinted), 
                                 ncol=lo, byrow=TRUE)

statistics <- t.test(rowMeans(nonimprintedLoci_matrix), rowMeans(imprintedLoci_matrix))
statistics$p.value
statistics$estimate[2] - statistics$estimate[1]

################################################################################

# Combine with RNA-seq

################################################################################

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/ChIPseq/mouseRNAseq"))

# Load statistics
load("topTable.RData")
load("geneAnnotation.RData")
topTable <- inner_join(topTable, geneAnnotation, by = c("ID" = "gene_id"))

enrichmentDF1 <- data.frame(Gene = c(egs_imprinted$external_gene_id,
                                    egs_nonimprinted$external_gene_id),
                           Enrichment =  c(rowMeans(imprintedLoci_matrix),
                                           rowMeans(nonimprintedLoci_matrix)),
                           Imprinted = c(rep("Yes", nrow(imprintedLoci_matrix)),
                                         rep("No", nrow(nonimprintedLoci_matrix))),
                           Distance = "TSS +/- 3000 bp"
                           )

enrichmentDF2 <- data.frame(Gene = c(egs_imprinted$external_gene_id,
                                    egs_nonimprinted$external_gene_id),
                           Enrichment =  c(rowMeans(imprintedLoci_matrix[,11:50]),
                                           rowMeans(nonimprintedLoci_matrix[,11:50])),
                           Imprinted = c(rep("Yes", nrow(imprintedLoci_matrix)),
                                         rep("No", nrow(nonimprintedLoci_matrix))),
                           Distance = "TSS +/- 2000 bp"
)

enrichmentDF3 <- data.frame(Gene = c(egs_imprinted$external_gene_id,
                                     egs_nonimprinted$external_gene_id),
                            Enrichment =  c(rowMeans(imprintedLoci_matrix[,21:40]),
                                            rowMeans(nonimprintedLoci_matrix[,21:40])),
                            Imprinted = c(rep("Yes", nrow(imprintedLoci_matrix)),
                                          rep("No", nrow(nonimprintedLoci_matrix))),
                            Distance = "TSS +/- 1000 bp"
)

enrichmentDF <- rbind.data.frame(enrichmentDF1, enrichmentDF2, enrichmentDF3)

plotDF <- inner_join(topTable, enrichmentDF, by = c("GeneName" = "Gene"))
plotDF <- unique(plotDF[,c("logFC", "PValue", "GeneName", "Enrichment", "Imprinted", "Distance")])

testDF <- plotDF[(plotDF$Distance == "TSS +/- 3000 bp") &
                   (plotDF$Imprinted == "Yes"),]
cor.test(abs(testDF$logFC), testDF$Enrichment, method = "spearman")

testDF <- plotDF[(plotDF$Distance == "TSS +/- 2000 bp") &
                   (plotDF$Imprinted == "Yes"),]
cor.test(abs(testDF$logFC), testDF$Enrichment, method = "spearman")

testDF <- plotDF[(plotDF$Distance == "TSS +/- 1000 bp") &
                   (plotDF$Imprinted == "Yes"),]
cor.test(abs(testDF$logFC), testDF$Enrichment, method = "spearman")


p <- ggplot() +
  geom_point(data = plotDF[plotDF$Imprinted == "Yes",], 
             aes(x = Enrichment, y = abs(logFC), color = logFC)) +
  facet_grid(rows = vars(Distance), space = "free", scale = "free") +
  xlab("Mean enrichment") +
  ylab(expression("|log"[2]~"FC|")) +
  scale_color_gradient2(low = "#2171B5", 
                       mid = "grey", 
                       high = "#CB181D", 
                       midpoint = 0,
                       limits = c(-2,2),
                       oob = scales::squish)+
  theme_bw()


ggsave(p, file = "logFCvsEnrichment.png", width = 7, height = 5)