#!/usr/bin/env sh
set -evx

ref="$1"
tbam="$2"
nbam="$3"
outdir="$4"
samplename="$5"

scriptdir="$(dirname "$(which "$0")")"

tvcfgz="${outdir}/${samplename}_T_uvc1.vcf.gz"
nvcfgz="${outdir}/${samplename}_N_uvc1.vcf.gz"
avcfgz="${outdir}/${samplename}_TandN_uvc1.vcf.gz"

tlog="${outdir}/${samplename}_T_uvc1.stderr"
nlog="${outdir}/${samplename}_N_uvc1.stderr"

mkdir -p "${outdir}"

date && time -p "${scriptdir}/uvc1" -f "${ref}" -s "${samplename}_T" "${tbam}" -o "${tvcfgz}" "${@:6}" 2> "${tlog}"
date && time -p "${scriptdir}/uvc1" -f "${ref}" -s "${samplename}_N" "${nbam}" -o "${nvcfgz}" "${@:6}" 2> "${nlog}"

date && time -p bcftools index -t "${tvcfgz}"
date && time -p bcftools index -t "${nvcfgz}"

date && time -p bcftools merge -m none -Ou "${tvcfgz}" "${nvcfgz}" | "${scriptdir}/callTN1" - "${avcfgz}"
date && time -p bcftools index -t "${avcfgz}"

