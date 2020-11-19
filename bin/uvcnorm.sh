#!/usr/bin/env bash

DEFAULT_NUM_THREADS=4
DEFAULT_MIN_SNV_QUAL=58.5
DEFAULT_MIN_NON_SNV_QUAL=49.5
DEFAULT_MIN_NLODQ=-9999

scriptdir="$(dirname "$(which "$0")")"

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input-vcf> <output-vcf> <multiallelic-control> <num-threads> [<min-SNV-QUAL>] [<min-non-SNV-QUAL>] [<min-NLODQ>]"
    echo "    input-vcf: uvc-generated VCF with non-normalized variants. "
    echo "    output-vcf: normalized VCF. "
    echo "    multiallelic-control: snps|indels|both|any, same as the --multiallelics option in bcftools. "
    echo "    num-threads: the number of threads to use, same as the --threads option in bcftools [${DEFAULT_NUM_THREADS}]. "
    echo "    min-SNV-qual: the minimum QUAL at each position below which the variant is always not merged [${DEFAULT_MIN_SNV_QUAL}]. "
    echo "    min-non-SNV-qual: the minimum QUAL  at each position below which the variant is always not merged [${DEFAULT_MIN_NON_SNV_QUAL}]. "
    echo "    min-NLODQ: the minimum NLODQ at each position below which the variant is filtered out. "
    echo "        This option can remove tumor SNV/InDel with InDel/SNV in the matched normal at the same position, respectively. "
    echo "        Set to a very negative value to disable this filter [${DEFAULT_MIN_NLODQ}]. "
    exit 1
fi

if [ -z "${4}" ]; then
    numthreads="${DEFAULT_NUM_THREADS}"
else
    numthreads="${4}"
fi

if [ -z "${5}" ]; then
    minSNVqual="${DEFAULT_MIN_SNV_QUAL}"
else
    minSNVqual="${5}"
fi

if [ -z "${6}" ]; then
    minNonSNVqual="${DEFAULT_MIN_NON_SNV_QUAL}"
else
    minNonSNVqual="${6}"
fi

if [ -z "${7}" ]; then
    minNLODQ="${DEFAULT_MIN_NLODQ}"
else
    minNLODQ="${7}"
fi

export PATH="${scriptdir}:${PATH}" # remove this line in the rare case that an important executable is shadowed by this command

bcftools view --threads $numthreads -i "(vNLODQ[0:0] > ${minNLODQ} && vNLODQ[0:1] > ${minNLODQ}) && (TYPE == 'snps' && QUAL  < ${minSNVqual} || TYPE != 'snps' && QUAL  < ${minNonSNVqual})" "${1}" -Oz -o "${2}.nonorm.vcf.gz"
bcftools index --threads $numthreads -ft "${2}.nonorm.vcf.gz"
bcftools view --threads $numthreads -i "(vNLODQ[0:0] > ${minNLODQ} && vNLODQ[0:1] > ${minNLODQ}) && (TYPE == 'snps' && QUAL >= ${minSNVqual} || TYPE != 'snps' && QUAL >= ${minNonSNVqual})" "${1}" | bcftools norm -m+${3} -Oz -o "${2}.norm.vcf.gz"
bcftools index --threads $numthreads -ft "${2}.norm.vcf.gz"
bcftools concat --threads $numthreads -a "${2}.nonorm.vcf.gz" "${2}.norm.vcf.gz" -Oz -o "${2}"
rm "${2}.nonorm.vcf.gz" "${2}.nonorm.vcf.gz.tbi" "${2}.norm.vcf.gz" "${2}.norm.vcf.gz.tbi"

