#!/usr/bin/env bash

set -euo pipefail

DIR=$(dirname $(realpath $0))
RESULTS_DIR=${DIR}/../results
BIN_DIR=${DIR}/../bin
CHRX_DIR=${DIR}/../raw/chrX_data

THREADS=4
OVERHANG=75
STAR_INDEX=${CHRX_DIR}/indexes/star_index_${OVERHANG}
RSEM_VER=1.3.3
RSEM_REF=${CHRX_DIR}/indexes/rsem_reference/chrX

if [[ ! -d ${RESULTS_DIR}/star_rsem ]]; then
   mkdir ${RESULTS_DIR}/star_rsem
fi
cd ${RESULTS_DIR}/star_rsem

for R1 in $(ls ${CHRX_DIR}/samples/*1.fastq.gz); do

   BASE=$(basename ${R1} _1.fastq.gz)
   mkdir -p ${BASE}
   R2=$(echo ${R1} | sed 's/1.fastq.gz/2.fastq.gz/')

   ${BIN_DIR}/STAR --runMode alignReads \
      --genomeDir ${STAR_INDEX} \
      --readFilesCommand "gunzip -c" \
      --outFileNamePrefix ${BASE}. \
      --readFilesIn ${R1} ${R2} \
      --runThreadN ${THREADS} \
      --twopassMode Basic \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.1 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --outFilterType BySJout \
      --outFilterScoreMinOverLread 0.33 \
      --outFilterMatchNminOverLread 0.33 \
      --limitSjdbInsertNsj 1200000 \
      --outSAMstrandField intronMotif \
      --outFilterIntronMotifs None \
      --alignSoftClipAtReferenceEnds Yes \
      --quantMode TranscriptomeSAM GeneCounts \
      --outSAMtype BAM Unsorted \
      --outSAMunmapped Within \
      --genomeLoad NoSharedMemory \
      --chimSegmentMin 15 \
      --chimJunctionOverhangMin 15 \
      --chimOutType Junctions WithinBAM SoftClip \
      --chimMainSegmentMultNmax 1 \
      --chimOutJunctionFormat 0 \
      --outSAMattributes NH HI AS nM NM ch \
      --outSAMattrRGline ID:rg1 SM:sm1

   ${BIN_DIR}/RSEM-${RSEM_VER}/rsem-calculate-expression \
      --num-threads ${THREADS} \
      --fragment-length-max 1000 \
      --paired-end \
      --estimate-rspd \
      --no-bam-output \
      --bam \
      ${BASE}.Aligned.toTranscriptome.out.bam ${RSEM_REF} ${BASE}.rsem

   gzip *.results
   mv *.bam *.out *.tab *._STAR* *.junction *.rsem.stat *.rsem.genes.results.gz *.rsem.isoforms.results.gz ${BASE}
done

>&2 echo Done
exit 0
