#!/usr/bin/env sh

scriptdir="$(dirname "$(which "$0")")"
if [ $# -lt 5 ]; then
    echo "Usage: $0 REF tumor-bam normal-bam output-directory tumor-sample-name[,normal-sample-name] [--tumor-parameters <tumor-params>] [--normal-parameters <normal-params>]"
    echo "     --tumor-parameters is optional and is followed by the parameters to ${scriptdir}/uvc1 for the  tumor-sample"
    echo "    --normal-parameters is optional and is followed by the parameters to ${scriptdir}/uvc1 for the normal-sample"
    exit 1
fi

tstate=1
nstate=1
tparams=()
nparams=()
for p in "${@:6}"; do
    cstate=0
    if [[ "${p}" = "--tumor-params" ]]; then
        tstate=1
        nstate=0
        cstate=1
    elif [[ "${p}" = "--normal-params" ]]; then
        tstate=0
        nstate=1
        cstate=1
    fi
    if [ $cstate -eq 0 ]; then
        if [ $tstate -eq 1 ]; then
            tparams+=($p)
        fi
        if [ $nstate -eq 1 ]; then
            nparams+=($p)
        fi
    fi
done

echo TUMOR-PARAMETERS "${tparams[@]}" 
echo NORMAL-PARAMETERS "${nparams[@]}" 

set -evx

ref="$1"
tbam="$2"
nbam="$3"
outdir="$4"
samplename="$5"

tvcfgz="${outdir}/${samplename}_T_uvc1.vcf.gz"
nvcfgz="${outdir}/${samplename}_N_uvc1.vcf.gz"
avcfgz="${outdir}/${samplename}_TandN_uvc1.vcf.gz"

tlog="${outdir}/${samplename}_T_uvc1.stderr"
nlog="${outdir}/${samplename}_N_uvc1.stderr"

mkdir -p "${outdir}"

if [ $(echo "${samplename}" | awk -F "," '{print NF}') -eq 2 ]; then
    tsample=$(echo "${samplename}" | awk -F "," '{print $1}')
    nsample=$(echo "${samplename}" | awk -F "," '{print $2}')
else
    tsample="${samplename}_T"
    nsample="${samplename}_N"
fi

date && time -p "${scriptdir}/uvc1" -f "${ref}" -s "${tsample}" "${tbam}" -o "${tvcfgz}" "${tparams[@]}" 2> "${tlog}"
date && time -p bcftools index -t "${tvcfgz}"

date && time -p "${scriptdir}/uvc1" -f "${ref}" -s "${nsample}" "${nbam}" -o "${nvcfgz}" "${nparams[@]}" --tumor-vcf "${tvcfgz}" 2> "${nlog}"
date && time -p bcftools index -t "${nvcfgz}"

#date && time -p bcftools merge -m none -Ou "${tvcfgz}" "${nvcfgz}" | "${scriptdir}/callTN1" - "${avcfgz}"
#date && time -p bcftools index -t "${avcfgz}"

