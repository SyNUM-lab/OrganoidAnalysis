# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(data.table)
library(biomaRt)

# Directory with RSEM output
dataDir <- "C:/Users/jarno/Desktop/Validation/FASTQ/Genes"

# Working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/ChIPseq/mouseRNAseq"))

# Get raw gene expression values:

  
  # Get gene expression files
  files <- list.files(path = dataDir, pattern = ".genes.results")
  
  # Get sample names
  sampleNames <- str_remove(files, ".genes.results")
  
  for (j in 1:length(files)){
    
    # Read gene expression file
    GeneExpr <- fread(paste0(dataDir,"/",files[j]))
    print(paste0(j,": ", nrow(GeneExpr)))
    
    # Retrieve expression value from file
    GeneExpr <- GeneExpr[,c("gene_id", "expected_count")]
    
    # Set column name to sample name
    colnames(GeneExpr) <- c("gene_id",sampleNames[j])
    
    # Combine with the already-collected gene expression values
    if (j==1){
      GeneExpr_all <- GeneExpr
    } else{
      GeneExpr_all <- inner_join(GeneExpr_all,GeneExpr, by = c("gene_id" = "gene_id"))
    }
  }


# Format data
gxMatrix_raw <- as.matrix(GeneExpr_all[,-1])
rownames(gxMatrix_raw) <- GeneExpr_all$gene_id

# Save expression matrix
save(gxMatrix_raw, file = "gxMatrix_raw.RData")

# Make meta data
metaData <- data.frame(
  SampleID = c("SRR2119607","SRR2119608","SRR2119609",
               "SRR2119610","SRR2119611","SRR2119612"),
  Group = c("WT", "WT", "WT",
            "KO", "KO", "KO")
)

save(metaData, file = "metaData.RData")


# Make gene annotation file:

# Read gene expression file
GeneExpr <- fread(paste0(dataDir,"/",files[j]))

# Retrieve relevant columns
geneAnnotation <- GeneExpr[,c("gene_id", "transcript_id(s)", "length", "effective_length")]

# Add gene name and Ensembl ID to the file
geneAnnotation$EnsemblID <- unlist(lapply(str_split(geneAnnotation$gene_id, "_"), function(x){x[1]}))
geneAnnotation$GeneName <- unlist(lapply(str_split(geneAnnotation$gene_id, "_"), function(x){x[2]}))

# Change order of columns
geneAnnotation <- geneAnnotation[,c("gene_id", "EnsemblID", "GeneName", "transcript_id(s)", "length", "effective_length")]

# Add entrez gene id and chromosome name to the gene annotation
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "entrezgene_id",
                                  "chromosome_name"), 
                     filters = 'ensembl_gene_id',
                     values = geneAnnotation$EnsemblID,
                     mart = ensembl)

geneAnnotation <- left_join(geneAnnotation, annotations, by = c("EnsemblID" = "ensembl_gene_id"))

# Save file
save(geneAnnotation, file = "geneAnnotation.RData")



################################################################################

# Preprocessing

################################################################################


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(stringr)
library(edgeR)
library(patchwork)

# Set working directory
homeDir <- "D:/RTTproject/CellAnalysis/OrganoidAnalysis"
setwd(paste0(homeDir,"/3. Validation/ChIPseq/mouseRNAseq"))


#*****************************************************************************#
# Normalization
#*****************************************************************************#

# Load data
load("gxMatrix_raw.RData")                                            # Raw count data
load("geneAnnotation.RData")                                          # Gene annotation data
load("metaData.RData")                                               # Sample Information
all(metaData$SampleID== colnames(gxMatrix_raw))

# Remove sex chromosomal genes
gxMatrix_fil <- gxMatrix_raw[rownames(gxMatrix_raw) %in% geneAnnotation$gene_id[(geneAnnotation$chromosome_name != "X") & (geneAnnotation$chromosome_name != "Y")],]

# Make DGE list object
y <- DGEList(counts = gxMatrix_fil,
             group = metaData$Group)

# remove lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)
gxMatrix_norm <- log2(cpm(y) + 1)
save(gxMatrix_norm, file = "gxMatrix_norm.RData")

#*****************************************************************************#
# DE analysis: MECP2 KO vs WT
#*****************************************************************************#

group <- metaData$Group

# Make design matrix
design <- model.matrix(~ 0 + group)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(y, design)

# Make Contrasts: RTT vs IC at every time point
all_contrasts <- makeContrasts(
  KOvsWT = groupKO-groupWT,
  levels = design
)

test <- glmQLFTest(fit, contrast = all_contrasts)

#Get top table
topTable <- topTags(test, n = Inf)$table
topTable$ID <- rownames(topTable)

# Save differential expression analysis results
save(topTable, file = "topTable.RData")