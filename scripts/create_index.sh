#!/usr/bin/env bash

set -euo pipefail

num_threads=8
gtf=chrX_data/genes/gencode.v37.annotation.chrx.gtf
overhang=75

if [[ ! -d chrX_data/indexes/star_index_${overhang} ]]; then
   >&2 echo Building STAR index
   ../src/STAR \
      --runMode genomeGenerate \
      --genomeDir chrX_data/indexes/star_index_${overhang} \
      --genomeFastaFiles chrX_data/genome/chrX.fa \
      --genomeSAindexNbases 12 \
      --sjdbGTFfile ${gtf} \
      --sjdbOverhang ${overhang} \
      --runThreadN ${num_threads}
fi

if [[ ! -d chrX_data/indexes/rsem_reference ]]; then
   >&2 echo Building RSEM index
   mkdir -p chrX_data/indexes/rsem_reference
   ../src/RSEM-1.3.3/rsem-prepare-reference \
      chrX_data/genome/chrX.fa \
      chrX_data/indexes/rsem_reference/chrX \
      --gtf ${gtf} \
      --num-threads ${num_threads}
fi

if [[ ! -d chrX_data/indexes/hisat2_index ]]; then
   >&2 echo Building HISAT2 index
   mkdir chrX_data/indexes/hisat2_index
   ../src/hisat2-2.2.0/extract_splice_sites.py ${gtf} > chrX.ss
   ../src/hisat2-2.2.0/extract_exons.py ${gtf} > chrX.exon
   ../src/hisat2-2.2.0/hisat2-build -p ${num_threads} --ss chrX.ss --exon chrX.exon chrX_data/genome/chrX.fa chrX_data/indexes/hisat2_index/chrX_tran
   rm chrX.ss chrX.exon
fi

if [[ ! -d chrX_data/indexes/kallisto_index ]]; then
   >&2 echo Building kallisto index
   mkdir chrX_data/indexes/kallisto_index
   ../src/kallisto/kallisto index -i chrX_data/indexes/kallisto_index/transcripts.idx chrX_data/genes/gencode.v37.transcripts.chrx.fa
fi

>&2 echo Done
exit 0

