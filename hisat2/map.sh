#!/usr/bin/env bash

if ! [ -x "$(command -v samtools)" ]; then
  >&2 echo 'Please install samtools'
  exit 1
fi

samtools=$(command -v samtools)
hisat2=../src/hisat2-2.2.0/hisat2
num_threads=4
index=../raw/chrX_data/indexes/chrX_tran

map_fastq () {
   fq1=$1
   fq2=$(echo $fq1 | sed 's/_1.fastq.gz/_2.fastq.gz/')
   prefix=$(basename $fq1 _1.fastq.gz)
   >&2 echo "Mapping ${prefix}"
   ${hisat2} --dta \
             -p ${num_threads} \
             -x ${index} \
             -1 ${fq1} \
             -2 ${fq2} |\
             ${samtools} sort -O BAM |\
             tee ${prefix}.bam |\
             ${samtools} index - ${prefix}.bam.bai
}

for fastq in $(ls ../raw/chrX_data/samples/*_1.fastq.gz); do
   map_fastq ${fastq}
done

exit 0

