# Clear workspace and console
rm(list = ls())
cat("\014") 

# Set working directory
setwd("E:/RTTproject/PublicDatasets")

# Load packages
library(tidyverse)
library(GEOquery)
library(data.table)


#*****************************************************************************#
#   GSE117511
#*****************************************************************************#
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE117511", "file=GSE117511_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

annotation <- fread("Human.GRCh38.p13.annot.tsv")
rownames(annotation) <- as.character(annotation$GeneID)

gene_names <- rep(NA, nrow(tbl))
for (i in 1:nrow(tbl)){
  gene_names[i] <- annotation[which(rownames(annotation) %in% rownames(tbl)[i])[1], "Symbol"]
}
rownames(tbl) <- gene_names


test <- getGEO("GSE117511")
sampleInfo <- test$GSE117511_series_matrix.txt.gz@phenoData@data[,c(1:2,11:14)]
colnames(sampleInfo) <- c("Title", "SampleID", "Genotype", "Stimulation", "Treatment", "Time")
sampleInfo$Genotype <- str_remove(sampleInfo$Genotype, "genotype/variation: ")
sampleInfo$Stimulation <- str_remove(sampleInfo$Stimulation, "stimulation: ")
sampleInfo$Treatment <- str_remove(sampleInfo$Treatment, "treatment: ")
sampleInfo$Time <- str_remove(sampleInfo$Time, "time: ")

# Filter for no treatment and stimulation
sampleInfo <- sampleInfo[(sampleInfo$Stimulation == "None") & (sampleInfo$Treatment == "None"),]
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = "GSE117511/gxMatrix_raw.RData")
save(sampleInfo, file = "GSE117511/sampleInfo.RData")


#*****************************************************************************#
#   GSE107399
#*****************************************************************************#
accessID <- "GSE107399"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

annotation <- fread("Human.GRCh38.p13.annot.tsv")
rownames(annotation) <- as.character(annotation$GeneID)

gene_names <- rep(NA, nrow(tbl))
for (i in 1:nrow(tbl)){
  gene_names[i] <- annotation[which(rownames(annotation) %in% rownames(tbl)[i])[1], "Symbol"]
}
rownames(tbl) <- gene_names


# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE107399_series_matrix.txt.gz@phenoData@data[5:30,c(1:2,8,12)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Genotype")
sampleInfo$Genotype <- str_remove(sampleInfo$Genotype, "phenotype: ")

samples <- intersect(colnames(tbl),sampleInfo$SampleID)
sampleInfo <- sampleInfo[samples,]

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))

#*****************************************************************************#
#   GSE123753
#*****************************************************************************#
accessID <- "GSE123753"

# Read data
tbl <- fread(paste0(accessID,"/GSE123753_counts_gene_replicates.csv"))

# Filter data
tbl_fil <- data.frame(tbl)[,str_detect(colnames(tbl), "Input")]
tbl_fil <- data.frame(tbl_fil)[,str_detect(colnames(tbl_fil), "NumReads")]

# Set column names
colnames(tbl_fil) <- str_remove(colnames(tbl_fil), "NumReads.Input_")

# Set row names
tbl_fil <- tbl_fil[!is.na(tbl$gene),]
rownames(tbl_fil) <- tbl$gene[!is.na(tbl$gene)]

# Make sample info
sampleInfo <- data.frame(SampleID = colnames(tbl_fil),
                         Tissue = unlist(lapply(str_split(colnames(tbl_fil), "_"), function(x) x[[1]])),
                         Genotype = unlist(lapply(str_split(colnames(tbl_fil), "_"), function(x) x[[2]])),
                         Replicate = unlist(lapply(str_split(colnames(tbl_fil), "_"), function(x) x[[3]]))
                         )

# Make final gene expression matrix
gxMatrix <- tbl_fil

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))


# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE123753_series_matrix.txt.gz@phenoData@data[,c(1:2,8,11)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Genotype")
sampleInfo$Genotype <- str_remove(sampleInfo$Genotype, "genotype: ")
sampleInfo <- sampleInfo[!str_detect(sampleInfo$Title, "TRAP"),]
sampleInfo <- sampleInfo[sampleInfo$Tissue != "Induced pluripotent stem cells",]

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))


#*****************************************************************************#
#   GSE165577
#*****************************************************************************#
accessID <- "GSE165577"

test <- fread(paste0(accessID, "/", "GSM5044261_D56_iCtrl_Raw_counts.csv"))

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)


# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE128380_series_matrix.txt.gz@phenoData@data[,c(1:2,8,11)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Genotype")
sampleInfo$Genotype <- str_remove(sampleInfo$Genotype, "genotype: ")
sampleInfo <- sampleInfo[!str_detect(sampleInfo$Title, "TRAP"),]
sampleInfo <- sampleInfo[sampleInfo$Tissue != "Induced pluripotent stem cells",]

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))

#*****************************************************************************#
#   GSE128380
#*****************************************************************************#
accessID <- "GSE128380"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

annotation <- fread("Human.GRCh38.p13.annot.tsv")
rownames(annotation) <- as.character(annotation$GeneID)

gene_names <- rep(NA, nrow(tbl))
for (i in 1:nrow(tbl)){
  gene_names[i] <- annotation[which(rownames(annotation) %in% rownames(tbl)[i])[1], "Symbol"]
}
rownames(tbl) <- gene_names

# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE128380_series_matrix.txt.gz@phenoData@data[,c(1:2,8,11)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Age")
sampleInfo$Age <- str_remove(sampleInfo$Age, "age: ")
sampleInfo$Age <- as.numeric(str_remove(sampleInfo$Age, " years"))

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))

#*****************************************************************************#
#   GSE128380
#*****************************************************************************#
accessID <- "GSE128380"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste0(urld, "&acc=",accessID, "&file=",accessID,"_raw_counts_GRCh38.p13_NCBI.tsv.gz");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Get meta data from GEO
test <- getGEO(accessID)
sampleInfo <- test$GSE128380_series_matrix.txt.gz@phenoData@data[,c(1:2,8,11)]
colnames(sampleInfo) <- c("Title", "SampleID", "Tissue", "Age")
sampleInfo$Age <- str_remove(sampleInfo$Age, "age: ")
sampleInfo$Age <- as.numeric(str_remove(sampleInfo$Age, " years"))

# filter matrix
gxMatrix <- tbl[,sampleInfo$SampleID]

# Save data
save(gxMatrix, file = paste0(accessID,"/gxMatrix_raw.RData"))
save(sampleInfo, file = paste0(accessID,"/sampleInfo.RData"))

