#!/usr/bin/env bash

set -euo pipefail

if [[ ! -e chrX_data.tar.gz ]]; then
   wget -c https://davetang.org/file/chrX_data.tar.gz
   tar xvzf chrX_data.tar.gz
fi

ver=37

if [[ ! -e chrX_data/genes/gencode.v${ver}.annotation.gtf ]]; then
   wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${ver}/gencode.v${ver}.annotation.gtf.gz
   gunzip gencode.v${ver}.annotation.gtf.gz
   cat gencode.v${ver}.annotation.gtf | grep "^chrX" > gencode.v${ver}.annotation.chrx.gtf
   cat gencode.v37.annotation.chrx.gtf | perl -lane 'next unless $F[2] eq "transcript"; $tid = $1 if /transcript_id\s\"(.*?)\";/; print "$tid"' > gencode.v37.annotation.chrx.tid.txt
   mv gencode.v${ver}.annotation.gtf gencode.v${ver}.annotation.chrx.gtf gencode.v37.annotation.chrx.tid.txt chrX_data/genes/
fi

if [[ ! -e chrX_data/genes/gencode.v${ver}.transcripts.fa ]]; then
   wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${ver}/gencode.v${ver}.transcripts.fa.gz
   gunzip gencode.v${ver}.transcripts.fa.gz
   ../scripts/extract_fasta.pl -i chrX_data/genes/gencode.v37.annotation.chrx.tid.txt -f gencode.v${ver}.transcripts.fa > gencode.v${ver}.transcripts.chrx.fa
   mv gencode.v${ver}.transcripts.fa gencode.v${ver}.transcripts.chrx.fa chrX_data/genes/
fi

>&2 echo Done
exit 0
