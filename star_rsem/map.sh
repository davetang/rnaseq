#!/usr/bin/env bash

set -euo pipefail

usage() {
   >&2 echo "Usage: $0 [ -a r1.fastq.gz ] [ -b r2.fastq.gz ] [ -p num_thread ] [ -i star_index ] [ -r true ]"
   exit 1
}

while getopts ":a:b:p:i:r:" options; do
  case "${options}" in
    a)
      read1=${OPTARG}
      ;;
    b)
      read2=${OPTARG}
      ;;
    i)
      star_index=${OPTARG}
      ;;
    r)
      remove_boolean=${OPTARG}
      ;;
    p)
      num_thread=${OPTARG}
      regex='^[1-9][0-9]*$'
      if [[ ! ${num_thread} =~ $regex ]]; then
        usage
      fi
      ;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      exit 1
      ;;
    *)
      usage ;;
  esac
done

if [ $OPTIND -ne 11 ]; then
   usage
fi

base=$(basename $read1 _1.fastq.gz)

../src/STAR --runMode alignReads \
     --genomeDir ${star_index} \
     --readFilesCommand "gunzip -c" \
     --outFileNamePrefix ${base}. \
     --readFilesIn ${read1} ${read2} \
     --runThreadN ${num_thread} \
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

if [[ ${remove_boolean} == true ]]; then
   rm ${base}.Aligned.out.bam
fi

exit 0

