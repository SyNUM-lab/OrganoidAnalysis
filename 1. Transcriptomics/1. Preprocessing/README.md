# Transcriptomics - Preprocessing
1. `runFASTQC.sh` - run FASTQC to perform quality control on the FASTQ files.
2. `SetUpRSEM.sh` - set up RSEM.
3. `runRSEM.sh` - run RSEM to get the gene and transcript expression levels.
4. `PrepareData.R` - prepare RSEM output for preprocessing and differential expression analysis.
5. `Preprocessing.R` - perform preprocessing and differential expression analysis.
6. `EvaluateApoptosis.R` - evaluate the expression of apoptotic markers over time.
7. `GetBiotype.R` - get gene biotype using biomaRt.