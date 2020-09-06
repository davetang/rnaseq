#!/usr/bin/env bash

set -euo pipefail

usage() {
   >&2 echo "Usage: $0 [ -b transcriptome.bam ] [ -p num_thread ] [ -i rsem_index ] [ -r true ]"
   exit 1
}

while getopts ":b:p:i:r:" options; do
  case "${options}" in
    b)
      bam=${OPTARG}
      ;;
    i)
      rsem_index=${OPTARG}
      ;;
    r)
      no_bam_output=${OPTARG}
      ;;
    p)
      num_thread=${OPTARG}
      regex='^[1-9][0-9]*$'
      if [[ ! $num_thread =~ $regex ]]; then
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

if [ $OPTIND -ne 9 ]; then
   usage
fi

base=$(basename ${bam} .bam)

if [[ ${no_bam_output} == true ]]; then
   ../src/RSEM-1.3.3/rsem-calculate-expression \
      --num-threads ${num_thread} \
      --fragment-length-max 1000 \
      --paired-end \
      --estimate-rspd \
      --no-bam-output \
      --bam \
      ${bam} ${rsem_index} ${base}.rsem
else
   ../src/RSEM-1.3.3/rsem-calculate-expression \
      --num-threads ${num_thread} \
      --fragment-length-max 1000 \
      --paired-end \
      --estimate-rspd \
      --bam \
      ${bam} ${rsem_index} ${base}.rsem
fi

gzip *.results

exit 0

