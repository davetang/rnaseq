#!/usr/bin/env bash

set -euo pipefail

DIR=$(dirname $(realpath $0))
DATA_DIR=${DIR}/../raw
CHRX_DIR=${DIR}/../raw/chrX_data

if [[ ! -d ${CHRX_DIR} ]]; then
   wget -c https://davetang.org/file/chrX_data.tar.gz -O ${DATA_DIR}/chrX_data.tar.gz
   tar --directory ${DATA_DIR} -xzf ${DATA_DIR}/chrX_data.tar.gz
   rm ${DATA_DIR}/chrX_data.tar.gz
else
   >&2 echo ${CHRX_DIR} already exists
fi

GENCODE_VER=46

if [[ ! -e ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.gtf ]]; then
   # wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VER}/gencode.v${GENCODE_VER}.annotation.gtf.gz -O ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.gtf.gz
   wget -c https://davetang.org/file/gencode/Gencode_human/release_${GENCODE_VER}/gencode.v${GENCODE_VER}.annotation.gtf.gz -O ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.gtf.gz
   gunzip ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.gtf.gz

   cat ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.gtf \
      | grep "^chrX" \
      > ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.gtf

   cat ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.gtf \
      | perl -lane 'next unless $F[2] eq "transcript"; $tid = $1 if /transcript_id\s\"(.*?)\";/; print "$tid"' \
      > ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.tid.txt
else
   >&2 echo ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.gtf already exists
fi

if [[ ! -e ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.fa ]]; then
   # wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VER}/gencode.v${GENCODE_VER}.transcripts.fa.gz -O ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.fa.gz
   wget -c https://davetang.org/file/gencode/Gencode_human/release_${GENCODE_VER}/gencode.v${GENCODE_VER}.transcripts.fa.gz -O ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.fa.gz
   gunzip ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.fa.gz

   ${DIR}/extract_fasta.pl -i ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.annotation.chrx.tid.txt -f ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.fa > ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.chrx.fa
else
   >&2 echo ${CHRX_DIR}/genes/gencode.v${GENCODE_VER}.transcripts.fa already exists
fi

>&2 echo Done
exit 0
