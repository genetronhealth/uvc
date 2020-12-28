#!/usr/bin/env bash

# WARNING: this script is still in alpha stage and not production-ready !!!

# params:  outvcf invcf ref bam[,bed] InDel-length
# output: outvcf: vcf with surrogate alignments

set -evx

scriptdir="$(dirname "$(which "$0")")"

export PATH="${scriptdir}:${PATH}" # remove this line in the rare case that an important executable is shadowed by this command

if [ -z "${UVC_BIN_EXE_FULL_NAME}" ]; then
    UVC_BIN_EXE_FULL_NAME="${scriptdir}/uvc1"
else
    echo "WARNING: using UVC_BIN_EXE_FULL_NAME=${UVC_BIN_EXE_FULL_NAME} from environment variable."
    echo "Please enter the shell command (unset UVC_BIN_EXE_FULL_NAME) before running uvcTN.sh if the default uvc binary exe full path should be used."
fi

BWA_SURROGATE_PARAMS=" -A 3 -B 12 -O 18 -E 1 -L 18 "
GROUNDTRUTH_ALLELE_IDENTITY_FLAG="both"

outvcf="${1}"
invcf="${2}"
ref="${3}"
bam=$(echo "${4}" | awk -F "," '{print $1}')
bed=$(echo "${4}" | awk -F "," '{print $2}')
indelsize=14

outdir="${outvcf}.surrogate"
#fq="${outdir}/surrogate.fastq.gz"
fq0="${outdir}/surrogate.SE.fastq.gz"
fq1="${outdir}/surrogate.R1.fastq.gz"
fq2="${outdir}/surrogate.R2.fastq.gz"
#sfq="${outdir}/surrogate.singleton.fastq.gz"

ncpus=32
nbams=8

mkdir -p "${outdir}"
if [ -z "${bed}" ]; then
    bed="${outdir}/superactive.bed"
    echo "Will generate the bed file ${bed} which contains superactive regions used for surrogate alignments. "
    printf "track name=superactive description=\"Containing super-active regions (regions with noisy alignments) for surrogate alignments.\"\n" > "${bed}"
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' -i "ALT = \"<ADDITIONAL_INDEL_CANDIDATE>\"" "${invcf}" | bedtools slop -b 120 -g "${ref}.fai" -i - | bedtools merge -i - >> "${bed}"
fi

samtools view -@ ${nbams} -L "${bed}" "${bam}" -bhu | samtools sort -@ ${nbams} -n -u -o - | samtools fastq -@ ${nbams} - -s "${fq0}" -1 "${fq1}" -2 "${fq2}"
bwa mem ${BWA_SURROGATE_PARAMS} -t ${ncpus} "${ref}" "${fq1}" "${fq2}" | samtools sort -@ ${nbams} - -o "${outdir}/surrogate.bam"
samtools index -@ ${nbams} "${outdir}/surrogate.bam"

sample=$(bcftools view "${invcf}" --header-only | tail -n1 |awk '{print $NF}')
"${UVC_BIN_EXE_FULL_NAME}" -t ${ncpus} --outvar-flag 0xF -f "${ref}" -s "${sample}" -o "${outdir}/surrogate.vcf.gz" "${outdir}/surrogate.bam"
bcftools index --threads ${nbams} -ft "${outdir}/surrogate.vcf.gz"

#bcftools view --threads ${nbams} -i "TYPE != \"indel\" || (abs(strlen(ALT)-strlen(REF)) < ${indelsize})" "${invcf}" -Oz -o "${outdir}/original-filt.vcf.gz"

bcftools view --threads ${nbams} "${invcf}" -Oz -o "${outdir}/original-filt.vcf.gz"
bcftools index --threads ${nbams} "${outdir}/original-filt.vcf.gz"
bcftools view --threads ${nbams} -i "TYPE = \"indel\" && (abs(strlen(ALT)-strlen(REF)) > ${indelsize}) && GERMLINE=1 && GT != \"ref\"" "${outdir}/surrogate.vcf.gz" -Oz -o "${outdir}/surrogate-filt.vcf.gz"
bcftools index --threads ${nbams} "${outdir}/surrogate-filt.vcf.gz"

bcftools concat --threads ${nbams} -a -d ${GROUNDTRUTH_ALLELE_IDENTITY_FLAG} "${outdir}/surrogate-filt.vcf.gz" "${outdir}/original-filt.vcf.gz" -Oz -o "${outvcf}"
bcftools index  --threads ${nbams} -ft "${outvcf}"

