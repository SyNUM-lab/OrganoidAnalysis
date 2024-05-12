#!/bin/bash

# Reference genome:
FASTA=GRCh38

# Start of loop:
for f in $PWD/*.bam ; do
  FILE=`basename ${f}`
  FILE=${FILE%.*}

#     Reorder BAM file to match exactly with FASTA reference file contigs.
java -jar tools/picard.jar ReorderSam INPUT=${FILE}.bam \
 OUTPUT=${FILE}.r.bam CREATE_INDEX=true SEQUENCE_DICTIONARY=$PWD/ref/${FASTA}.dict \
 ALLOW_INCOMPLETE_DICT_CONCORDANCE=true ALLOW_CONTIG_LENGTH_DISCORDANCE=true

#     Generate additional read group characteristics (required by GATK)
java -jar $PWD/tools/picard.jar AddOrReplaceReadGroups \
       I=${FILE}.r.bam \
       O=${FILE}.rg.bam \
       RGID=test \
       RGLB=test \
       RGPL=illumina \
       RGPU=test \
       RGSM=test \
       CREATE_INDEX=true
rm -f ${FILE}.r.ba*

#	ASE variant calling:
./tools/gatk/gatk CollectAllelicCounts \
          -I ${FILE}.rg.bam \
          -R ref/GRCh38.fa \
          -L ref/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
          -O $PWD/ASE_results/${FILE}.outputTable.tsv

done

echo $SECONDS

# END
