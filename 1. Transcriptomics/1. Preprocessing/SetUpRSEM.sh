
# Download reference genome
wget "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
wget "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz"
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.110.chr.gtf.gz

# Prepare reference genome for RSEM
rsem-prepare-reference --gtf Homo_sapiens.GRCh38.110.chr.gtf --bowtie2 \
Homo_sapiens.GRCh38.dna.primary_assembly.fa human_ref
