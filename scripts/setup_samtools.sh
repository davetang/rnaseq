#!/usr/bin/env bash

set -euo pipefail

VER=1.20
TOOL=samtools
URL=https://github.com/samtools/samtools/releases/download/${VER}/${TOOL}-${VER}.tar.bz2
SCRIPT_DIR=$(dirname $(realpath $0))
OUTDIR=$(realpath ${SCRIPT_DIR}/..)

if [[ -d ${OUTDIR}/bin/samtools ]]; then
   >&2 echo ${OUTDIR}/samtools already exists
   exit
fi

TMPDIR=$(mktemp -d)
trap "rm -rf ${TMPDIR}" SIGINT SIGTERM

cd ${TMPDIR}
wget ${URL}
tar xjf ${TOOL}-${VER}.tar.bz2
cd ${TOOL}-${VER}
./configure --prefix=${OUTDIR}
make && make install
cd

rm -rf ${TMPDIR}

${OUTDIR}/bin/samtools --version

>&2 echo Done
exit 0
