#!/usr/bin/env bash

set -euo pipefail

if [[ ! -e chrX_data/indexes/star_index_74 ]]; then
../src/STAR \
   --runMode genomeGenerate \
   --genomeDir chrX_data/indexes/star_index_74 \
   --genomeFastaFiles chrX_data/genome/chrX.fa \
   --genomeSAindexNbases 12 \
   --sjdbGTFfile chrX_data/genes/chrX.gtf \
   --sjdbOverhang 74 \
   --runThreadN 8
fi

if [[ ! -e chrX_data/indexes/rsem_reference ]]; then
mkdir -p chrX_data/indexes/rsem_reference
../src/RSEM-1.3.3/rsem-prepare-reference \
	chrX_data/genome/chrX.fa \
   chrX_data/indexes/rsem_reference/chrX \
   --gtf chrX_data/genes/chrX.gtf \
   --num-threads 8
fi

