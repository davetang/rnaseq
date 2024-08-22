#!/usr/bin/env bash

set -euo pipefail

tool=samtools
ver=1.12
url=https://github.com/samtools/samtools/releases/download/${ver}/${tool}-${ver}.tar.bz2
dir=$(pwd)/${tool}

if [[ -d samtools ]]; then
   >&2 echo samtools already exists
   exit 0
else
   mkdir ${tool}
fi

wget ${url}
tar xjf ${tool}-${ver}.tar.bz2
cd ${tool}-${ver}
./configure --prefix=${dir}
make && make install
cd ..

rm -rf ${tool}-${ver} ${tool}-${ver}.tar.bz2

${dir}/bin/samtools --version

>&2 echo Done

exit 0

