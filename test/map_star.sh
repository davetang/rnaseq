#!/usr/bin/env bash

set -euo pipefail

star_index=../raw/chrX_data/indexes/star_index_74
read1=eg.fa

base=default
../src/STAR --runMode alignReads \
   --genomeDir ${star_index} \
   --readFilesIn ${read1} \
   --outFileNamePrefix ${base}.
rm ${base}.Log.final.out ${base}.Log.out ${base}.Log.progress.out ${base}.SJ.out.tab

base=multimap_20
../src/STAR --runMode alignReads \
   --genomeDir ${star_index} \
   --readFilesIn ${read1} \
   --outFileNamePrefix ${base}. \
   --outFilterMultimapNmax 20
rm ${base}.Log.final.out ${base}.Log.out ${base}.Log.progress.out ${base}.SJ.out.tab


base=multimap_20_unmapped
../src/STAR --runMode alignReads \
   --genomeDir ${star_index} \
   --readFilesIn ${read1} \
   --outFileNamePrefix ${base}. \
   --outFilterMultimapNmax 20 \
   --outSAMunmapped Within
rm ${base}.Log.final.out ${base}.Log.out ${base}.Log.progress.out ${base}.SJ.out.tab

