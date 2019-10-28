#!/usr/bin/env sh

scriptdir="$(dirname "$(which "$0")")"
if [ $# -lt 5 ]; then
    echo "Usage: $0 <REF> <tumor-bam> <normal-bam> <output-directory> <tumor-sample-name>[,<normal-sample-name>][/<nprocs>/<parallel|qsub>] [<all-params>] [--tumor-params <tumor-params>] [--normal-params <normal-params>]"
    echo "    The output bgzipped vcf file is <output-directory>/<normal-sample-name>_uvc1.vcf.gz if normal-sample-name is provided or <output-directory>/<tumor-sample-name>_N_uvc1.vcf.gz if normal-sample-name is not provided"
    echo "    <nprocs> is the number of processes corresponding to the number of chromosomes that are run concurrently in parallel, where 0 (zero by default) means no chromosome-level parallelization"
    echo "    <parallel|qsub> means the string parallel or qsub, where parallel requires GNU parallel to be installed and qsub requires the variable UVC_QSUB_CMD to be set. For example, UVC_QSUB_CMD can be \"qsub -V -S /bin/sh\""
    echo "    <all-params> is the set of parameters to ${scriptdir}/uvc1 for both tumor and normal samples"
    echo "     --tumor-parameters is optional and is followed by the parameters to ${scriptdir}/uvc1 for only the  tumor-sample"
    echo "    --normal-parameters is optional and is followed by the parameters to ${scriptdir}/uvc1 for only the normal-sample"
    echo "    For help on the usage of \"${scriptdir}/uvc1\", please enter \"${scriptdir}/uvc1\" -h "
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
samplename=$(echo "$5/0/parallel" | awk -F"/" '{print $1}')
nprocs=$(echo "$5/0/parallel" | awk -F"/" '{print $2}')
paratool=$(echo "$5/0/parallel" | awk -F"/" '{print $2}')

mkdir -p "${outdir}"

if [ $(echo "${samplename}" | awk -F "," '{print NF}') -eq 2 ]; then
    tsample=$(echo "${samplename}" | awk -F "," '{print $1}')
    nsample=$(echo "${samplename}" | awk -F "," '{print $2}')
else
    tsample="${samplename}_T"
    nsample="${samplename}_N"
fi

tvcfgz="${outdir}/${tsample}_uvc1.vcf.gz"
nvcfgz="${outdir}/${nsample}_uvc1.vcf.gz"

tlog="${outdir}/${tsample}_uvc1.stderr"
nlog="${outdir}/${nsample}_uvc1.stderr"

export PATH="${scriptdir}:${PATH}" # remove this line in the rare case that an important executable is shadowed by this command

if [ "${nprocs}" -gt 0 ]; then
    tnames=$(cat "${ref}.fai" | awk '{print $1}')
    if [ "${paratool}" -eq "parallel" ]; then
        for tname in ${tnames}; do
            echo "${ref}" "${tbam}" "${nbam}" "${outdir}/${tname}" "${samplename}" -t "${tname}" "${@:6}" 
        done > "${outdir}/run_parallel.sh"
        cat "${outdir}/run_parallel.sh" | parallel -j "${nprocs}"
        bcftools concat -n -Oz -o "${nvcfgz}" "${outdir}/"*"/${nsample}_uvc1.vcf.gz"
    elif [ "${paratool}" -eq "qsub" ]; then
        if [ -z "${UVC_QSUB_CMD}" ]; then
            echo "The variable UVC_QSUB_CMD must be set and exported in order to use qsub!"
            exit -2
        fi
        for tname in ${tnames}; do
            echo "${ref}" "${tbam}" "${nbam}" "${outdir}/${tname}" "${samplename}" -t "${tname}" "${@:6}" "|" "${UVC_QSUB_CMD}" -o "${outdir}" -e "${outdir}" -j "${tname}.job"
        done > "${outdir}/run_qsub.sh"
        sh "${outdir}/run_qsub.sh"
    else
        echo "The multiprocessing tool ${paratool} is neither parallel nor qsub."
        exit -1
    fi
else
    date 
    "${scriptdir}/uvc1" -f "${ref}" -s "${tsample}" "${tbam}" -o "${tvcfgz}" "${tparams[@]}" 2> "${tlog}"
    date 
    bcftools index -ft "${tvcfgz}" # or use tabix, requires htslib 1.6 or plus

    date
    "${scriptdir}/uvc1" -f "${ref}" -s "${nsample}" "${nbam}" -o "${nvcfgz}" "${nparams[@]}" --tumor-vcf "${tvcfgz}" 2> "${nlog}"
    date
    bcftools index -ft "${nvcfgz}" # or use tabix, requires htslib 1.6 or plus
fi

date

