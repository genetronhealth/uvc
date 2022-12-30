#!/usr/bin/env sh

# cmdline param 1: option for this script
# cmdline param 2: option for the configure script in bcftools and samtools

set -evx
currdir="${PWD}"

if [ $(echo "${1}" | grep skip-bcftools | wc -l) -eq 0 ]; then
    mkdir -p "${currdir}/ext/"
    cd "${currdir}/ext/"
    if [ $(echo "${1}" | grep skip-downloading-bcftools | wc -l) -eq 0 ]; then
        wget --inet4-only https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2
    fi
    tar -xvf bcftools-1.11.tar.bz2
    cd "${currdir}/ext/bcftools-1.11"
    ./configure ${2}
    make -j 4
    cp bcftools "${currdir}/bin/"
    
    # htslib-*-lowdep is used for compiling UVC
    cp -r "${currdir}/ext/bcftools-1.11/htslib-1.11" "${currdir}/ext/htslib-1.11-lowdep"
    cd "${currdir}/ext/htslib-1.11-lowdep"
    ./configure -disable-plugins --disable-libcurl --disable-s3 --disable-largefile --without-libdeflate ${2} # --disable-bz2 and --disable-lzma are both for disabling CRAM files
    make -j 4
    
    # make install # this command may fail without root privilege, but it does not matter much as bcftools is in the PATH variable by exporting in uvcTN.sh
fi

# samtools is not directly used by UVC but is often used for processing the input to UVC.
if [ $(echo "${1}" | grep skip-samtools | wc -l) -eq 0 ]; then
    mkdir -p "${currdir}/ext/"
    cd "${currdir}/ext/"
    if [ $(echo "${1}" | grep skip-downloading-samtools | wc -l) -eq 0 ]; then
        wget --inet4-only https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
    fi
    tar -xvf samtools-1.11.tar.bz2
    cd "${currdir}/ext/samtools-1.11"
    ./configure ${2}
    make -j 4
    cp samtools "${currdir}/bin/"
fi

if [ $(echo "${1}" | grep skip-parallel | wc -l) -eq 0 ]; then
    cd "${currdir}/ext/"
    wget --inet4-only http://ftp.gnu.org/gnu/parallel/parallel-20201122.tar.bz2
    tar -xvf parallel-20201122.tar.bz2
    cd "${currdir}/ext/parallel-20201122"
    ./configure
    make
    cp src/parallel "${currdir}/bin/"
fi

