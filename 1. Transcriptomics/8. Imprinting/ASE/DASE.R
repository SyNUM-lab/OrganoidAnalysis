# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readxl)
library(biomaRt)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/8. Imprinting/ASE"))

#==============================================================================#
# Get output files
#==============================================================================#

# List all files
files <- list.files("ASE_results", full.names = TRUE)

# Collect chromosome, position, reference, alternative allele, and ASE score
chr <- matrix(NA, nrow = 240621, ncol = 42)
pos <- matrix(NA, nrow = 240621, ncol = 42)
ref <- matrix(NA, nrow = 240621, ncol = 42)
alt <- matrix(NA, nrow = 240621, ncol = 42)
total <- matrix(NA, nrow = 240621, ncol = 42)
score <- matrix(NA, nrow = 240621, ncol = 42)
for (f in 1:length(files)){
  temp <- readLines(files[f])
  temp_fil <- temp[!str_detect(temp, "@")]
  
  
  ASE_matrix <- matrix(NA, nrow = length(temp_fil), ncol  = 6)
  for (i in 1:length(temp_fil)){
    ASE_matrix[i,] <- str_split(temp_fil[i], "\t")[[1]]
  }
  
  colnames(ASE_matrix) <- ASE_matrix[1,]
  ASE_matrix <- ASE_matrix[-1,]
  
  chr[,f] <- ASE_matrix[,"CONTIG"]
  pos[,f] <- as.numeric(ASE_matrix[, "POSITION"])
  ref[,f] <- as.numeric(ASE_matrix[, "REF_COUNT"])
  alt[,f] <- as.numeric(ASE_matrix[,"ALT_COUNT"])
  total[,f] <- ref[,f] + alt[,f]
  score[,f] <- (abs((ref[,f]/(total[,f]))-0.5)+0.5)
}

# Make dataframe with id, chromosome, and position for each SNP
pos_info <- data.frame(
  id = paste0(chr[,1], ":", pos[,1]),
  chr = chr[,1],
  pos = as.numeric(pos[,1])
)

rownames(ref) <- pos_info$id
rownames(alt) <- pos_info$id
rownames(total) <- pos_info$id
rownames(score) <- pos_info$id

# Add sample names to objects
sample_names <- str_remove(files, "D:/RTTproject/GeneticFiles/ASE/ASE_results/")
sample_names <- str_remove(sample_names, ".sorted.outputTable.tsv")
colnames(ref) <- sample_names
colnames(alt) <- sample_names
colnames(total) <- sample_names
colnames(score) <- sample_names

#==============================================================================#
# Get gene names for each position
#==============================================================================#

# Get HG38 annotations
granges <- TxDb.Hsapiens.UCSC.hg38.knownGene
granges <- transcriptsBy(granges, by = "gene")

# Make genomic ranges from SNP dataframe
gr_pos <- makeGRangesFromDataFrame(
  data.frame(
    chr=pos_info$chr,
    start=pos_info$pos,
    end=pos_info$pos+1),
  keep.extra.columns=TRUE)

# Overlap between SNPs and genes
overlaps <- findOverlaps(granges, gr_pos, type = 'any')

# Add gene names to SNP data frame
idx_query <- queryHits(overlaps)
genes1 <- names(granges)[idx_query] 

idx_subject <- subjectHits(overlaps)
pos_info1 <- pos_info[idx_subject,]

pos_info1$EntrezID <- genes1

#==============================================================================#
# Get imprinted genes
#==============================================================================#

# Get imprinted genes
imprintDF <- read_xlsx(paste0(homeDir,"/1. Transcriptomics/8. Imprinting/ImprintedGenes.xlsx"))
imprintDF <- imprintDF[imprintDF$Status %in% unique(imprintDF$Status)[c(1,2,3)],]
imprinted_genes <- imprintDF$Gene[!(imprintDF$Status == unique(imprintDF$Status)[3])]
imprinted_genes <- imprinted_genes[!is.na(imprinted_genes)]

# Get entrez gene ID for the imprinted genes
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
annotations <- getBM(attributes=c("hgnc_symbol", "entrezgene_id"), 
                     filters = "hgnc_symbol",
                     values = imprinted_genes,
                     mart = ensembl)

# Select SNP info for imprinted genes only
imprinted_genes_entrezid <- as.character(unique(annotations$entrezgene_id[!is.na(annotations$entrezgene_id)]))
pos_info_imprinted <- pos_info1[pos_info1$EntrezID %in% imprinted_genes_entrezid,]

# Get HGNC symbol for the imprinted genes
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
annotations <- getBM(attributes=c("hgnc_symbol", "entrezgene_id"), 
                     filters = "entrezgene_id",
                     values = pos_info_imprinted$EntrezID,
                     mart = ensembl)

annotations$entrezgene_id <- as.character(annotations$entrezgene_id)
pos_info_imprinted <- left_join(pos_info_imprinted,
                                annotations,
                                by = c("EntrezID" = "entrezgene_id"))

# Save data
save(ref, alt, total, score, pos_info1, pos_info_imprinted, 
     file = "outputASE.RData")


#==============================================================================#
# Identify genes with differential allele specific expression
#==============================================================================#

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readxl)
library(ggpubr)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/1. Transcriptomics/8. Imprinting/ASE"))

# Load data
load("outputASE.RData")
load(paste0(homeDir,"/1. Transcriptomics/1. Preprocessing/geneAnnotation.RData"))

# Filter for imprinted genes
ref_imp <- ref[rownames(ref) %in% pos_info_imprinted$id,]
alt_imp <- alt[rownames(alt) %in% pos_info_imprinted$id,]
total_imp <- total[rownames(total) %in% pos_info_imprinted$id,]
score_imp <- score[rownames(score) %in% pos_info_imprinted$id,]

# Make sample information data frame
sampleInfo <- data.frame(
  SampleID = colnames(alt),
  Region = "iPSC",
  Time = "D0",
  Group = "IC",
  Replicate = "1"
)

sampleInfo$Region[str_detect(sampleInfo$SampleID, "Ventral")] <- "Ventral"
sampleInfo$Region[str_detect(sampleInfo$SampleID, "Dorsal")] <- "Dorsal"

sampleInfo$Time[str_detect(sampleInfo$SampleID, "D13")] <- "D13"
sampleInfo$Time[str_detect(sampleInfo$SampleID, "D40")] <- "D40"
sampleInfo$Time[str_detect(sampleInfo$SampleID, "D75")] <- "D75"

sampleInfo$Group[str_detect(sampleInfo$SampleID, "RTT")] <- "RTT"

sampleInfo$TimeRegion <- paste0(sampleInfo$Time, "_", sampleInfo$Region)
sampleInfo$TimeRegionGroup <- paste0(sampleInfo$Time, "_", 
                                     sampleInfo$Region, "_",
                                     sampleInfo$Group)

sampleInfo$Replicate[str_detect(sampleInfo$SampleID, "_2")] <- "2"
sampleInfo$Replicate[str_detect(sampleInfo$SampleID, "_3")] <- "3"

#------------------------------------------------------------------------------#
# 1) Get bi-allelically expressed imprinted genes (in IC)
#------------------------------------------------------------------------------#

sel_ids <- list()
for (tr in unique(sampleInfo$TimeRegion)){
  
  # total reads per SNP
  total_TimeRegion <- total_imp[, colnames(total_imp) %in% sampleInfo$SampleID[sampleInfo$TimeRegion == tr]]
  
  # Reference and alternative allele count for IC samples
  ref_IC <- ref_imp[,(colnames(ref_imp) %in% colnames(total_TimeRegion))&
                      (str_detect(colnames(ref_imp), "IC"))]
  alt_IC <- alt_imp[,(colnames(alt_imp) %in% colnames(total_TimeRegion))&
                      (str_detect(colnames(alt_imp), "IC"))]
  score_IC <-  score_imp[,(colnames(score_imp) %in% colnames(total_TimeRegion))&
                         (str_detect(colnames(score_imp), "IC"))]
  
  # Total alt and ref count both larger than 5 and ASE score < 0.85
  score_IC <- score_imp[(rowSums((ref_IC > 5) & (alt_IC > 5) & (score_IC < 0.85))==3),
                        colnames(ref_IC)] 
  
  # ASE score < 0.85
  sel_ids[[tr]] <- rownames(score_IC)
}

ids_ic <- names(table(unlist(sel_ids))[table(unlist(sel_ids))>3])
pos_info_imprinted[pos_info_imprinted$id %in% ids_ic,]


#------------------------------------------------------------------------------#
# Get imprinted genes that are have differential ASE
#------------------------------------------------------------------------------#

sel_ids <- list()
for (tr in unique(sampleInfo$TimeRegion)){
  total_TimeRegion <- total_imp[, colnames(total_imp) %in% sampleInfo$SampleID[sampleInfo$TimeRegion == tr]]
  
  # Scores for RTT
  ref_RTT <- ref_imp[,(colnames(ref_imp) %in% colnames(total_TimeRegion))&
                       (str_detect(colnames(ref_imp), "RTT"))]
  alt_RTT <- alt_imp[,(colnames(alt_imp) %in% colnames(total_TimeRegion))&
                       (str_detect(colnames(alt_imp), "RTT"))]
  score_RTT <- score_imp[,(colnames(ref_imp) %in% colnames(total_TimeRegion))&
                           (str_detect(colnames(ref_imp), "RTT"))] 
  
  # Scores for IC
  ref_IC <- ref_imp[,(colnames(ref_imp) %in% colnames(total_TimeRegion))&
                      (str_detect(colnames(ref_imp), "IC"))]
  alt_IC <- alt_imp[,(colnames(alt_imp) %in% colnames(total_TimeRegion))&
                      (str_detect(colnames(alt_imp), "IC"))]
  score_IC <- score_imp[,(colnames(ref_imp) %in% colnames(total_TimeRegion))&
                           (str_detect(colnames(ref_imp), "IC"))] 
  
  
  # Get SNPs that are heterozygous in either RTT or IC
  het_ids <- ((rowSums((ref_RTT > 5) & (alt_RTT > 5) & (score_RTT < 0.85))==3)|
                (rowSums((ref_IC > 5) & (alt_IC > 5) & (score_IC < 0.85))==3)) 
  
  score_RTT <- score_RTT[het_ids,] 
  score_IC <- score_IC[het_ids,] 
  
  # Get P value
  pvalue <- rep(NA, nrow(score_RTT))
  FC <- rep(NA, nrow(score_RTT))
  for (i in 1:nrow(score_RTT)){
    pvalue[i] <- wilcox.test(score_RTT[i,], score_IC[i,])$p.value
    FC[i] <- mean(score_RTT[i,])/mean(score_IC[i,])
  }
  sel_ids[[tr]] <- rownames(score_RTT)[pvalue <= 0.1]
}

# differential ASE at more than 3 time points
ids <- names(table(unlist(sel_ids))[table(unlist(sel_ids))>3])
pos_info_imprinted[pos_info_imprinted$id %in% ids,]


#------------------------------------------------------------------------------#
# Make plot
#------------------------------------------------------------------------------#

# Prepare data for plotting
plotDF <- gather(as.data.frame(score_imp[ids,]))
plotDF$ID <- rep(ids,ncol(score_imp))
plotDF <- inner_join(plotDF,pos_info_imprinted, by = c("ID" = "id"))
plotDF <- inner_join(plotDF, sampleInfo, by = c("key" = "SampleID"))
plotDF$TimeRegion <- factor(plotDF$TimeRegion,levels = c("D75_Dorsal",
                                                         "D40_Dorsal",
                                                         "D13_Dorsal",
                                                         "D0_iPSC",
                                                         "D13_Ventral",
                                                         "D40_Ventral",
                                                         "D75_Ventral"))
# Set order of IDs
plotDF$ID <- factor(plotDF$ID, levels = c("chr1:228321726",
                                          "chr19:31276811",
                                          "chr19:58571238",
                                          "chr5:95755497",   
                                          "chr7:78018005",
                                          "chr14:100734523"))
empty_names <- setNames(rep("",6), levels(plotDF$ID))

# Main plot
mainPlot <- ggplot(plotDF) +
  geom_point(aes(y = value, x = TimeRegion, color = Group)) +
  geom_line(aes(y = value, x = TimeRegion, color = Group, group = paste0(Replicate, Group))) +
  facet_grid(cols = vars(ID)) +
  coord_flip() +
  xlab(NULL) +
  ylab("ASE score") +
  labs(color = NULL) +
  scale_y_continuous(breaks = c(0.5,0.6,0.7,0.8,0.9,1), 
                     labels = c("","0.6", "","","0.9","")) +
  scale_color_manual(values = c("#377EB8","#E41A1C")) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 7))

panel_colors <- c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC")

# convert to grob (to set colors of panels)
gp <- ggplotGrob(mainPlot) # where p is the original ggplot object
for(i in 1:6){
  grob.i <- grep("strip-t", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- panel_colors[i]
}

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotDF$TimeRegion,
  time = plotDF$Time,
  region = plotDF$Region))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time), angle = 90) +
  coord_flip()+
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
  geom_text(data = colSideColor[(str_detect(colSideColor$sample, "D40") |
                                   str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = region), angle = 90) +
  coord_flip() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# combine plots
p <- colSideColorPlot_tissue +
  colSideColorPlot_time +
  gp +
  patchwork::plot_layout(widths = c(0.5,0.5,6))

# Save plot
ggsave(p, file = "ASE_plot.png", width = 8, height = 4)

