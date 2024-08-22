#!/usr/bin/env bash

set -euo pipefail

DIR=$(realpath $(dirname $(realpath $0))/..)
RESULTS_DIR=${DIR}/results
CHRX_DIR=${DIR}/raw/chrX_data
GENCODE_VER=46
FASTA=$(realpath ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.chrx.fa)
GTF=$(realpath ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.gtf)
THREADS=6
OUTDIR=$(realpath ${RESULTS_DIR}/nfcore_rnaseq)

if [[ ! -d ${OUTDIR} ]]; then
   mkdir ${OUTDIR}
fi

export NXF_SINGULARITY_CACHEDIR=${HOME}/nf-core/sif

nextflow run ${HOME}/nf-core/rnaseq/3_14_0/main.nf \
    -resume \
    -with-report execution_report.html \
    -with-trace \
    -with-dag flowchart.html \
    --input ${DIR}/samplesheet.csv \
    --outdir ${OUTDIR} \
    --fasta ${FASTA} \
    --gtf ${GTF} \
    --aligner star_rsem \
    --save_reference \
    --skip_markduplicates \
    --skip_dupradar \
    --skip_deseq2_qc \
    --skip_stringtie \
    -profile singularity \
    --max_cpus ${THREADS} \
    --max_memory 60GB

>&2 echo Done
exit 0
