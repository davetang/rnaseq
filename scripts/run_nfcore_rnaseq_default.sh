#!/usr/bin/env bash

set -euo pipefail

DIR=$(realpath $(dirname $(realpath $0))/..)
RESULTS_DIR=${DIR}/results
CHRX_DIR=${DIR}/raw/chrX_data
GENCODE_VER=46
FASTA=$(realpath ${CHRX_DIR}/genome/chrX.fa)
GTF=$(realpath ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.gtf)
THREADS=6
MAXMEM=16GB
OUTDIR=$(realpath ${RESULTS_DIR})/nfcore_rnaseq_default
DATETIME=$(date +%Y_%m_%d_%H_%M_%S)

if [[ ! -d ${OUTDIR} ]]; then
   mkdir -p ${OUTDIR}
fi

export NXF_SINGULARITY_CACHEDIR=${HOME}/nf-core/sif

nextflow run ${HOME}/nf-core/rnaseq/3_14_0/main.nf \
    -resume \
    -with-report execution_report_${DATETIME}.html\
    -with-trace \
    -with-dag flowchart_${DATETIME}.html \
    --input ${DIR}/samplesheet.csv \
    --outdir ${OUTDIR} \
    --fasta ${FASTA} \
    --gtf ${GTF} \
    --aligner star_rsem \
    --save_reference true \
    -profile singularity \
    --max_cpus ${THREADS} \
    --max_memory ${MAXMEM}

>&2 echo Done
exit 0
