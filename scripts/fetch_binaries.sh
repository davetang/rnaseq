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

DIR=$(dirname $(realpath $0))
BIN_DIR=${DIR}/../bin

if [[ ! -d ${BIN_DIR} ]]; then
   mkdir ${BIN_DIR}
fi

if [[ ! -d ${BIN_DIR}/hisat2-2.2.1 ]]; then
   wget -c https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O ${BIN_DIR}/hisat2-2.2.1-Linux_x86_64.zip
   unzip -d ${BIN_DIR} ${BIN_DIR}/hisat2-2.2.1-Linux_x86_64.zip
   rm ${BIN_DIR}/hisat2-2.2.1-Linux_x86_64.zip
else
   >&2 echo ${BIN_DIR}/hisat2-2.2.1 already exists
fi

if [[ ! -d ${BIN_DIR}/stringtie-2.2.3.Linux_x86_64 ]]; then
   wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.3.Linux_x86_64.tar.gz -O ${BIN_DIR}/stringtie-2.2.3.Linux_x86_64.tar.gz
   tar --directory ${BIN_DIR} -xzf ${BIN_DIR}/stringtie-2.2.3.Linux_x86_64.tar.gz
   rm ${BIN_DIR}/stringtie-2.2.3.Linux_x86_64.tar.gz
else
   >&2 echo ${BIN_DIR}/stringtie-2.2.3.Linux_x86_64 already exists
fi

if [[ ! -e ${BIN_DIR}/STAR ]]; then
   wget -c https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64_static/STAR?raw=true -O ${BIN_DIR}/STAR
   chmod 755 ${BIN_DIR}/STAR
else
   >&2 echo ${BIN_DIR}/STAR already exists
fi

KALLISTO_VER=0.50.1
if [[ ! -d ${BIN_DIR}/kallisto ]]; then
   wget -c https://github.com/pachterlab/kallisto/releases/download/v${KALLISTO_VER}/kallisto_linux-v${KALLISTO_VER}.tar.gz -O ${BIN_DIR}/kallisto_linux-v${KALLISTO_VER}.tar.gz
   tar --directory ${BIN_DIR} -xzf ${BIN_DIR}/kallisto_linux-v${KALLISTO_VER}.tar.gz
   rm ${BIN_DIR}/kallisto_linux-v${KALLISTO_VER}.tar.gz
else
   >&2 echo ${BIN_DIR}/kallisto already exists
fi

>&2 echo Done
exit
