#!/usr/bin/env bash

set -euo pipefail

VER=1.3.3
URL=https://github.com/deweylab/RSEM/archive/v${VER}.tar.gz
SCRIPT_DIR=$(dirname $(realpath $0))
BIN_DIR=$(realpath ${SCRIPT_DIR}/../bin)

if [[ -e ${BIN_DIR}/RSEM-${VER}/rsem-calculate-expression ]]; then
   >&2 echo ${BIN_DIR}/rsem-calculate-expression already exists
   exit
fi

wget ${URL} -O ${BIN_DIR}/RSEM-v${VER}.tar.gz
tar --directory ${BIN_DIR} -xzf ${BIN_DIR}/RSEM-v${VER}.tar.gz
cd ${BIN_DIR}/RSEM-${VER} && make
./rsem-calculate-expression --version
cd
rm ${BIN_DIR}/RSEM-v${VER}.tar.gz

>&2 echo Done
exit 0
