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
OUTDIR=$(realpath ${RESULTS_DIR})/nfcore_rnaseq_ref

TMPDIR=$(mktemp -d)
head -2 ${DIR}/samplesheet.csv > ${TMPDIR}/samplesheet.csv

if [[ ! -d ${OUTDIR} ]]; then
   mkdir -p ${OUTDIR}
fi

export NXF_SINGULARITY_CACHEDIR=${HOME}/nf-core/sif

nextflow run ${HOME}/nf-core/rnaseq/3_14_0/main.nf \
    -with-trace \
    --input ${TMPDIR}/samplesheet.csv \
    --outdir ${OUTDIR} \
    --fasta ${FASTA} \
    --gtf ${GTF} \
    --aligner star_rsem \
    -profile singularity \
    --save-reference true \
    --skip_gtf_filter true \
    --skip_gtf_transcript_filter true \
    --skip_trimming true \
    --skip_markduplicates true \
    --skip_bigwig true \
    --skip_stringtie true \
    --skip_fastqc true \
    --skip_dupradar true \
    --skip_qualimap true \
    --skip_rseqc true \
    --skip_biotype_qc true \
    --skip_deseq2_qc true \
    --skip_multiqc true \
    --skip_qc true \
    --max_cpus ${THREADS} \
    --max_memory ${MAXMEM}

rm -rf ${TMPDIR}

>&2 echo Done
exit 0
