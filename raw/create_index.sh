#!/usr/bin/env bash

set -euo pipefail

num_threads=8

if [[ ! -d chrX_data/indexes/star_index_74 ]]; then
   >&2 echo Building STAR index
   ../src/STAR \
      --runMode genomeGenerate \
      --genomeDir chrX_data/indexes/star_index_74 \
      --genomeFastaFiles chrX_data/genome/chrX.fa \
      --genomeSAindexNbases 12 \
      --sjdbGTFfile chrX_data/genes/chrX.gtf \
      --sjdbOverhang 74 \
      --runThreadN ${num_threads}
fi

if [[ ! -d chrX_data/indexes/rsem_reference ]]; then
   >&2 echo Building RSEM index
   mkdir -p chrX_data/indexes/rsem_reference
   ../src/RSEM-1.3.3/rsem-prepare-reference \
      chrX_data/genome/chrX.fa \
      chrX_data/indexes/rsem_reference/chrX \
      --gtf chrX_data/genes/chrX.gtf \
      --num-threads ${num_threads}
fi

if [[ ! -d chrX_data/indexes/hisat2_index ]]; then
   >&2 echo Building HISAT2 index
   mkdir chrX_data/indexes/hisat2_index
   ../src/hisat2-2.2.0/extract_splice_sites.py chrX_data/genes/chrX.gtf > chrX.ss
   ../src/hisat2-2.2.0/extract_exons.py chrX_data/genes/chrX.gtf > chrX.exon
   ../src/hisat2-2.2.0/hisat2-build -p ${num_threads} --ss chrX.ss --exon chrX.exon chrX_data/genome/chrX.fa chrX_data/indexes/hisat2_index/chrX_tran
   rm chrX.ss chrX.exon
fi

>&2 echo Done

