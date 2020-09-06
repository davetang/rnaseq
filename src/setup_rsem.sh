#!/usr/bin/env bash

set -euo pipefail

if ! [ -x "$(command -v wget)" ]; then
  >&2 echo 'Please install wget'
  exit 1
fi

wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz -O RSEM-v1.3.3.tar.gz
tar xzf RSEM-v1.3.3.tar.gz
rm RSEM-v1.3.3.tar.gz
cd RSEM-1.3.3 && make
./rsem-calculate-expression --version

