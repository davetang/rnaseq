#!/usr/bin/env bash

set -euo pipefail

if [[ ! -e ../src/stringtie-2.1.4/stringtie ]]; then
  >&2 echo Please setup stringtie first
  exit 1
fi

stringtie=../src/stringtie-2.1.4/stringtie
gtf=../raw/chrX_data/genes/chrX.gtf
num_threads=8
index=../raw/chrX_data/indexes/hisat2_index/chrX_tran

quant () {
   bam=$1
   echo >&2 Quantifying ${bam}
   prefix=$(basename ${bam} .bam)
   ${stringtie} -e -B -p ${num_threads} \
      -G ${gtf} \
      -o ${prefix}.gtf \
      -A ${prefix}.abundance.tsv \
      ${bam}

   mv *.ctab ${prefix}.abundance.tsv ${prefix}.gtf ${prefix}
}

for bam in $(find . -name "*.bam" -print); do
   quant ${bam}
done

>&2 echo Done

exit 0

