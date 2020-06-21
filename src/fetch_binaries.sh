#!/usr/bin/env bash

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
      wget -c https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-OSX_x86_64/download -O hisat2-2.2.0-OSX_x86_64.zip
      wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.3b.Linux_x86_64.tar.gz
      unzip hisat2-2.2.0-OSX_x86_64.zip
      tar -xzvf stringtie-2.1.3b.Linux_x86_64.tar.gz
      rm hisat2-2.2.0-OSX_x86_64.zip stringtie-2.1.3b.Linux_x86_64.tar.gz
   ;; 
   linux*)
      >&2 echo "Preparing binaries for Linux"
      wget -c https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download -O hisat2-2.2.0-Linux_x86_64.zip
      wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.3b.Linux_x86_64.tar.gz
      unzip hisat2-2.2.0-Linux_x86_64.zip
      tar -xzvf stringtie-2.1.3b.Linux_x86_64.tar.gz
      rm hisat2-2.2.0-Linux_x86_64.zip stringtie-2.1.3b.Linux_x86_64.tar.gz
   ;;
  *)
      echo "Not sure which exe to download for $OSTYPE"
   ;;
esac

