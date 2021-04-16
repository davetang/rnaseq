#!/usr/bin/env bash

set -euo pipefail

thread=8
star_index=../raw/chrX_data/indexes/star_index_74
rsem_ref=../raw/chrX_data/indexes/rsem_reference/chrX

for r1 in $(ls ../raw/chrX_data/samples/*1.fastq.gz); do

   base=$(basename ${r1} _1.fastq.gz)
   mkdir -p ${base}
   r2=$(echo ${r1} | sed 's/1.fastq.gz/2.fastq.gz/')

   ./map.sh -a ${r1} -b ${r2} -p ${thread} -i ${star_index} -r true

   ./quant.sh -b ${base}.Aligned.toTranscriptome.out.bam -p ${thread} -i ${rsem_ref} -r true
   mv *.bam *.out *.tab *._STAR* *.junction *.rsem.stat *.rsem.genes.results.gz *.rsem.isoforms.results.gz ${base}

done

