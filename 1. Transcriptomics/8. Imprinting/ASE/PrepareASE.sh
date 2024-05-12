# Step 1
./faToTwoBit GRCh38.fa GRCh38.2bit

# Step 2
java -jar tools/picard.jar CreateSequenceDictionary REFERENCE=ref/GRCh38.fa OUTPUT=ref/GRCh38.dict

# Step 3: download reference genome
wget "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Step 4: index reference genome
samtools faidx ref/GRCh38.fa

# Step 5: Download VCF
wget "https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
wget "https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"

# More information: https://github.com/macsbio/AlleleSpecificExpression/blob/main/Code/ASEReadCounter/ASEReadCounter.sh