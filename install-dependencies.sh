#!/usr/bin/env sh

set -evx
currdir="${PWD}"

mkdir -p "${currdir}/ext/"
cd "${currdir}/ext/"
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 --inet4-only
tar -fvx bcftools-1.9.tar.bz2
cd "${currdir}/ext/bcftools-1.9"
./configure
make
cp bcftools "${currdir}/bin/"
make install # this command may fail without root privilege, but it does not matter much as bcftools is in the PATH variable by exporting in uvcTN.sh

