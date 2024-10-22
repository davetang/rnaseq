#!/usr/bin/env bash

set -euo pipefail

FASTA=/home/dtang/data/ensembl/release-112/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
GTF=/home/dtang/data/ensembl/release-112/Homo_sapiens.GRCh38.112.gtf.gz
THREADS=4
MAXMEM=50GB
DATETIME=$(date +%Y_%m_%d_%H_%M_%S)

export NXF_SINGULARITY_CACHEDIR=${HOME}/nf-core/sif

nextflow run ${HOME}/nf-core/rnaseq/3_14_0/main.nf \
    -resume \
    -with-report execution_report_${DATETIME}.html\
    -with-trace \
    -with-dag flowchart_${DATETIME}.html \
    --input samplesheet.csv \
    --outdir results \
    --fasta ${FASTA} \
    --gtf ${GTF} \
    --aligner star_rsem \
    --save_reference true \
    -profile singularity \
    --max_cpus ${THREADS} \
    --max_memory ${MAXMEM}

>&2 echo Done
exit 0
