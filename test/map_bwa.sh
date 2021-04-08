#!/usr/bin/env bash

set -euo pipefail

samtools=$(command -v samtools)
num_threads=8
index=../raw/chrX_data/genome/chrX.fa
fasta=eg.fa

bwa mem \
   -t ${num_threads} \
   ${index} \
   ${fasta} |
   ${samtools} sort -@ ${num_threads} -O BAM > bwa_default.bam

bwa mem \
   -t ${num_threads} \
   -h 10 \
   ${index} \
   ${fasta} |
   ${samtools} sort -@ ${num_threads} -O BAM > bwa_h10.bam

exit 0

