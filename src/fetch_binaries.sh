#!/usr/bin/env bash

set -euo pipefail

if ! [ -x "$(command -v wget)" ]; then
  >&2 echo 'Please install wget'
  exit 1
fi

if ! [ -x "$(command -v unzip)" ]; then
  >&2 echo 'Please install unzip'
  exit 1
fi

case "$OSTYPE" in
   darwin*)
      >&2 echo "Preparing binaries for macOS"
      if [[ ! -e hisat2-2.2.0 ]]; then
         wget -q -c https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-OSX_x86_64/download -O hisat2-2.2.0-OSX_x86_64.zip
         unzip -q hisat2-2.2.0-OSX_x86_64.zip
         rm hisat2-2.2.0-OSX_x86_64.zip
      fi
      if [[ ! -e stringtie-2.1.4.OSX_x86_64 ]]; then
         wget -q -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.4.OSX_x86_64.tar.gz
         tar -xzf stringtie-2.1.4.OSX_x86_64.tar.gz
         rm stringtie-2.1.4.OSX_x86_64.tar.gz
      fi
      if [[ ! -e STAR ]]; then
         wget -q -c https://github.com/alexdobin/STAR/blob/master/bin/MacOSX_x86_64/STAR?raw=true -O STAR
         chmod 755 STAR
      fi
   ;; 
   linux*)
      >&2 echo "Preparing binaries for Linux"
      if [[ ! -e hisat2-2.2.0 ]]; then
         wget -q -c https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download -O hisat2-2.2.0-Linux_x86_64.zip
         unzip -q hisat2-2.2.0-Linux_x86_64.zip
         rm hisat2-2.2.0-Linux_x86_64.zip
      fi
      if [[ ! -e stringtie-2.1.4.Linux_x86_64 ]]; then
         wget -q -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.4.Linux_x86_64.tar.gz
         tar -xzf stringtie-2.1.4.Linux_x86_64.tar.gz
         rm stringtie-2.1.4.Linux_x86_64.tar.gz
      fi
      if [[ ! -e STAR ]]; then
         wget -q -c https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64/STAR?raw=true -O STAR
         chmod 755 STAR
      fi
   ;;
  *)
      echo "Not sure which exe to download for $OSTYPE"
   ;;
esac

>&2 echo Done

