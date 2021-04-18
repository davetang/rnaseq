#!/usr/bin/env bash

set -euo pipefail

if ! [ -x "$(command -v esearch)" ]; then
  >&2 echo Please install Entrez Direct first
  exit 1
fi

runs=(ERR188044 ERR188104 ERR188234 ERR188245 ERR188257 ERR188273 ERR188337 ERR188383 ERR188401 ERR188428 ERR188454 ERR204916)

for run in ${runs[@]}; do
   echo ${run}
   esearch -db sra -query ${run} | efetch -format runinfo > ${run}.runinfo.txt; sleep 3
done

exit 0

