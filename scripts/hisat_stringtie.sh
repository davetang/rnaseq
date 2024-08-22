#!/usr/bin/env bash

set -euo pipefail

DIR=$(dirname $(realpath $0))
RESULTS_DIR=${DIR}/../results
BIN_DIR=${DIR}/../bin
DATA_DIR=${DIR}/../raw
CHRX_DIR=${DIR}/../raw/chrX_data

SAMTOOLS=${BIN_DIR}/samtools
HISAT2=${BIN_DIR}/hisat2-2.2.1/hisat2
THREADS=4
INDEX=${CHRX_DIR}/indexes/hisat2_index/chrX_tran

GENCODE_VER=46
STRINGTIE=${BIN_DIR}/stringtie-2.2.3.Linux_x86_64/stringtie
GTF=${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.gtf

map_fastq () {
   FQ1=$1
   FQ2=$(echo ${FQ1} | sed 's/_1.fastq.gz/_2.fastq.gz/')
   PREFIX=$(basename ${FQ1} _1.fastq.gz)

   mkdir ${PREFIX}
   >&2 echo Mapping ${PREFIX}
   ${HISAT2} --dta \
      -p ${THREADS} \
      -x ${INDEX} \
      -1 ${FQ1} \
      -2 ${FQ2} |\
      ${SAMTOOLS} sort -O BAM |\
      tee ${PREFIX}/${PREFIX}.bam |\
      ${SAMTOOLS} index - ${PREFIX}/${PREFIX}.bam.bai
}

quant () {
   BAM=$1
   echo >&2 Quantifying ${BAM}
   PREFIX=$(basename ${BAM} .bam)
   ${STRINGTIE} -e -B -p ${THREADS} \
      -G ${GTF} \
      -o ${PREFIX}.gtf \
      -A ${PREFIX}.abundance.tsv \
      ${BAM}

   mv *.ctab ${PREFIX}.abundance.tsv ${PREFIX}.gtf ${PREFIX}
}

if [[ ! -d ${RESULTS_DIR}/hisat_stringtie ]]; then
   mkdir -p ${RESULTS_DIR}/hisat_stringtie
fi

cd ${RESULTS_DIR}/hisat_stringtie

for FASTQ in $(ls ${CHRX_DIR}/samples/*_1.fastq.gz); do
   map_fastq ${FASTQ}
done

for BAM in $(find . -name "*.bam" -print); do
   quant ${BAM}
done

>&2 echo Done
exit 0
