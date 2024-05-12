#!/bin/bash

for i in $(ls *.fq.1.gz | sed -r 's/[.]fq[.]1[.]gz//' | uniq)
do
rsem-calculate-expression --paired-end --bowtie2 --append-names --no-bam-output -p 7\
  "${i}.fq.1.gz" \
  "${i}.fq.2.gz" \
  referenceGenome/human_ref \
  "rsemOutput/${i}"
done
