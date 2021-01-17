#!/usr/bin/env bash

set -euo pipefail

samtools=$(command -v samtools)
hisat2=../src/hisat2-2.2.0/hisat2
num_threads=8
index=../raw/chrX_data/indexes/chrX_tran
fasta=eg.fa

${hisat2} --dta \
   -f \
   -p ${num_threads} \
   -x ${index} \
   -U ${fasta} |
   ${samtools} sort -O BAM > hisat2_default.bam

${hisat2} --dta \
   -f \
   -p ${num_threads} \
   -k 20 \
   -x ${index} \
   -U ${fasta} |
   ${samtools} sort -O BAM > hisat2_k20.bam

exit 0

