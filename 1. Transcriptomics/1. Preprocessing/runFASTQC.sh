#!/usr/bin/env bash
DATA_PATH=$PWD/FASTQ
OUTPUT_PATH=$PWD/FASTQC
for file in $(ls ${DATA_PATH}/*fastq.gz)
do
    SAMPLE=`basename $file`
    SAMPLE=${SAMPLE%.*}
    SAMPLE=${SAMPLE%.*}
    /mnt/d/RTTproject/GeneticFiles/FastQC/fastqc -t 5 ${DATA_PATH}/${SAMPLE}.fastq.gz -o ${OUTPUT_PATH}
done
