#!/usr/bin/env bash
DATA_PATH=$PWD/FASTQ
OUTPUT_PATH=$PWD/BAM
for file in $(ls ${DATA_PATH}/*fastq.gz)
do
    SAMPLE=`basename $file`
    SAMPLE=${SAMPLE%.*}
    SAMPLE=${SAMPLE%.*}
    bowtie2 -x /mnt/d/RTTproject/GeneticFiles/referenceGenome/human_ref --no-unal -U ${DATA_PATH}/${SAMPLE}.fastq.gz | samtools view -bS - > ${OUTPUT_PATH}/${SAMPLE}.bam
    samtools sort ${OUTPUT_PATH}/${SAMPLE}.bam -o ${OUTPUT_PATH}/${SAMPLE}.sorted.bam
    samtools index ${OUTPUT_PATH}/${SAMPLE}.sorted.bam
done
