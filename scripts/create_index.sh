#!/usr/bin/env bash

set -euo pipefail

DIR=$(dirname $(realpath $0))
DATA_DIR=${DIR}/../raw
BIN_DIR=${DIR}/../bin
CHRX_DIR=${DIR}/../raw/chrX_data

RSEM_VER=1.3.3
HISAT2_VER=2.2.1
GENCODE_VER=46

THREADS=4
GTF=${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.gtf
OVERHANG=75

if [[ ! -d ${CHRX_DIR}/indexes/star_index_${OVERHANG} ]]; then
   >&2 echo Building STAR index
   ${BIN_DIR}/STAR \
      --runMode genomeGenerate \
      --genomeDir ${CHRX_DIR}/indexes/star_index_${OVERHANG} \
      --genomeFastaFiles ${CHRX_DIR}/genome/chrX.fa \
      --genomeSAindexNbases 12 \
      --sjdbGTFfile ${GTF} \
      --sjdbOverhang ${OVERHANG} \
      --runThreadN ${THREADS}
else
   >&2 echo ${CHRX_DIR}/indexes/star_index_${OVERHANG} already exists
fi

if [[ ! -d ${CHRX_DIR}/indexes/rsem_reference ]]; then
   >&2 echo Building RSEM index
   mkdir -p ${CHRX_DIR}/indexes/rsem_reference
   ${BIN_DIR}/RSEM-${RSEM_VER}/rsem-prepare-reference \
      ${CHRX_DIR}/genome/chrX.fa \
      ${CHRX_DIR}/indexes/rsem_reference/chrX \
      --gtf ${GTF} \
      --num-threads ${THREADS}
else
   >&2 echo ${CHRX_DIR}/indexes/rsem_reference already exists
fi

if [[ ! -d ${CHRX_DIR}/indexes/kallisto_index ]]; then
   >&2 echo Building kallisto index
   mkdir ${CHRX_DIR}/indexes/kallisto_index
   ${BIN_DIR}/kallisto/kallisto index -i ${CHRX_DIR}/indexes/kallisto_index/transcripts.idx ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.chrx.fa
else
   >&2 echo ${CHRX_DIR}/indexes/kallisto_index already exists
fi

if [[ ! -d ${CHRX_DIR}/indexes/hisat2_index ]]; then
   >&2 echo Building HISAT2 index
   mkdir ${CHRX_DIR}/indexes/hisat2_index
   ${BIN_DIR}/hisat2-${HISAT2_VER}/extract_splice_sites.py ${GTF} > chrX.ss
   ${BIN_DIR}/hisat2-${HISAT2_VER}/extract_exons.py ${GTF} > chrX.exon
   ${BIN_DIR}/hisat2-${HISAT2_VER}/hisat2-build -p ${THREADS} --ss chrX.ss --exon chrX.exon ${CHRX_DIR}/genome/chrX.fa ${CHRX_DIR}/indexes/hisat2_index/chrX_tran
   rm chrX.ss chrX.exon
else
   >&2 echo ${CHRX_DIR}/indexes/hisat2_index already exists
fi

>&2 echo Done
exit 0
