#!/usr/bin/env bash

set -euo pipefail

if [[ ! -e ../src/samtools/bin/samtools ]]; then
  >&2 echo Please setup samtools first
  exit 1
fi

if [[ ! -e ../src/hisat2-2.2.0/hisat2 ]]; then
  >&2 echo Please setup hisat2 first
  exit 1
fi

samtools=../src/samtools/bin/samtools
hisat2=../src/hisat2-2.2.0/hisat2
num_threads=8
index=../raw/chrX_data/indexes/hisat2_index/chrX_tran

map_fastq () {
   fq1=$1
   fq2=$(echo $fq1 | sed 's/_1.fastq.gz/_2.fastq.gz/')
   prefix=$(basename $fq1 _1.fastq.gz)
   mkdir ${prefix}
   >&2 echo Mapping ${prefix}
   ${hisat2} --dta \
             -p ${num_threads} \
             -x ${index} \
             -1 ${fq1} \
             -2 ${fq2} |\
             ${samtools} sort -O BAM |\
             tee ${prefix}/${prefix}.bam |\
             ${samtools} index - ${prefix}/${prefix}.bam.bai
}

for fastq in $(ls ../raw/chrX_data/samples/*_1.fastq.gz); do
   map_fastq ${fastq}
done

>&2 echo Done

exit 0

