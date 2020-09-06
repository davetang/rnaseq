#!/usr/bin/env bash

set -euo pipefail

if ! [ -x "$(command -v wget)" ]; then
  >&2 echo 'Please install wget'
  exit 1
fi

wget https://github.com/alexdobin/STAR/archive/2.7.5b.tar.gz -O STAR-2.7.5b.tar.gz
tar xzf STAR-2.7.5b.tar.gz
rm STAR-2.7.5b.tar.gz
cd STAR-2.7.5b/source && make STAR
./STAR --version

