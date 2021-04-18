#!/usr/bin/env bash

set -euo pipefail

check_depend (){
   tool=$1
   if [[ ! -x $(command -v ${tool}) ]]; then
     >&2 echo Could not find ${tool}
     exit 1
   fi
}

dependencies=(wget unzip)
for tool in ${dependencies[@]}; do
   check_depend ${tool}
done

case ${OSTYPE} in
   darwin*)
      >&2 echo Downloading binaries for macOS
      if [[ ! -d hisat2-2.2.0 ]]; then
         wget -q -c https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-OSX_x86_64/download -O hisat2-2.2.0-OSX_x86_64.zip
         unzip -q hisat2-2.2.0-OSX_x86_64.zip
         rm hisat2-2.2.0-OSX_x86_64.zip
         >&2 echo Finished downloading hisat2
      fi
      if [[ ! -d stringtie-2.1.4 ]]; then
         wget -q -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.4.OSX_x86_64.tar.gz
         tar -xzf stringtie-2.1.4.OSX_x86_64.tar.gz
         mv stringtie-2.1.4.OSX_x86_64 stringtie-2.1.4
         rm stringtie-2.1.4.OSX_x86_64.tar.gz
         >&2 echo Finished downloading stringtie
      fi
      if [[ ! -e STAR ]]; then
         wget -q -c https://github.com/alexdobin/STAR/blob/master/bin/MacOSX_x86_64/STAR?raw=true -O STAR
         chmod 755 STAR
         >&2 echo Finished downloading STAR
      fi
      if [[ ! -d kallisto ]]; then
         wget -q -c https://github.com/pachterlab/kallisto/releases/download/v0.46.2/kallisto_mac-v0.46.2.tar.gz
         tar -xzf kallisto_mac-v0.46.2.tar.gz
         rm kallisto_mac-v0.46.2.tar.gz
         >&2 echo Finished downloading kallisto
      fi
   ;; 
   linux*)
      >&2 echo Downloading binaries for Linux
      if [[ ! -d hisat2-2.2.0 ]]; then
         wget -q -c https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download -O hisat2-2.2.0-Linux_x86_64.zip
         unzip -q hisat2-2.2.0-Linux_x86_64.zip
         rm hisat2-2.2.0-Linux_x86_64.zip
         >&2 echo Finished downloading hisat2
      fi
      if [[ ! -d stringtie-2.1.4.Linux_x86_64 ]]; then
         wget -q -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.4.Linux_x86_64.tar.gz
         tar -xzf stringtie-2.1.4.Linux_x86_64.tar.gz
         mv stringtie-2.1.4.Linux_x86_64 stringtie-2.1.4
         rm stringtie-2.1.4.Linux_x86_64.tar.gz
         >&2 echo Finished downloading stringtie
      fi
      if [[ ! -e STAR ]]; then
         wget -q -c https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64/STAR?raw=true -O STAR
         chmod 755 STAR
         >&2 echo Finished downloading STAR
      fi
      if [[ ! -d kallisto ]]; then
         wget -q -c https://github.com/pachterlab/kallisto/releases/download/v0.46.2/kallisto_linux-v0.46.2.tar.gz
         tar -xzf kallisto_linux-v0.46.2.tar.gz
         rm kallisto_linux-v0.46.2.tar.gz
         >&2 echo Finished downloading kallisto
      fi
   ;;
  *)
      >&2 echo Not sure which binaries to download for ${OSTYPE}
   ;;
esac

>&2 echo Done

