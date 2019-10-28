#!/usr/bin/env sh

set -evx
currdir="${PWD}"

if [ $(echo "${1}" | grep skip-bcftools | wc -l) -eq 0 ]; then
    mkdir -p "${currdir}/ext/"
    cd "${currdir}/ext/"
    wget --inet4-only https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 
    tar -xvf bcftools-1.9.tar.bz2
    cd "${currdir}/ext/bcftools-1.9"
    ./configure
    make -j 4
    cp bcftools "${currdir}/bin/"
    #make install # this command may fail without root privilege, but it does not matter much as bcftools is in the PATH variable by exporting in uvcTN.sh
fi

if [ $(echo "${1}" | grep skip-parallel | wc -l) -eq 0 ]; then
    cd "${currdir}/ext/"
    wget --inet4-only http://ftp.gnu.org/gnu/parallel/parallel-20191022.tar.bz2
    tar -xvf parallel-20191022.tar.bz2
    cd "${currdir}/ext/parallel-20191022"
    ./configure
    make
    cp src/parallel "${currdir}/bin/"
fi

