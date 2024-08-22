#!/usr/bin/env bash

set -euo pipefail

if [[ ! -e ../src/kallisto/kallisto ]]; then
  >&2 echo Please setup kallisto first
  exit 1
fi

kallisto=../src/kallisto/kallisto

for r1 in $(ls ../raw/chrX_data/samples/*_1.fastq.gz); do
   prefix=$(basename ${r1} _1.fastq.gz)
   r2=$(echo ${r1} | sed 's/_1.fastq.gz/_2.fastq.gz/')
   >&2 echo Quantifying ${prefix}
   ${kallisto} quant \
      -i ../raw/chrX_data/indexes/kallisto_index/transcripts.idx \
      -o ${prefix} \
      -b 100 \
      ${r1} ${r2}
   rm ${prefix}/abundance.h5
done

>&2 echo Done

exit 0

