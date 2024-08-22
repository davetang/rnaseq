#!/usr/bin/env bash

set -euo pipefail

DIR=$(dirname $(realpath $0))
RESULTS_DIR=${DIR}/../results
BIN_DIR=${DIR}/../bin
CHRX_DIR=${DIR}/../raw/chrX_data

if [[ ! -d ${RESULTS_DIR}/kallisto ]]; then
   mkdir ${RESULTS_DIR}/kallisto
fi
cd ${RESULTS_DIR}/kallisto

for R1 in $(ls ${CHRX_DIR}/samples/*_1.fastq.gz); do
   PREFIX=$(basename ${R1} _1.fastq.gz)
   R2=$(echo ${R1} | sed 's/_1.fastq.gz/_2.fastq.gz/')

   >&2 echo Quantifying ${PREFIX}
   ${BIN_DIR}/kallisto/kallisto quant \
      -i ${CHRX_DIR}/indexes/kallisto_index/transcripts.idx \
      -o ${PREFIX} \
      -b 100 \
      ${R1} ${R2}
   rm ${PREFIX}/abundance.h5
done

>&2 echo Done
exit 0
