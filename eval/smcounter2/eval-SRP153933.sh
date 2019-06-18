#!/usr/bin/env bash
# need GNU parallel, java8, fastq-dump (with aspera), bwa, samtools, bcftools, and vcfeval (from RTG)

set -evx
shopt -s expand_aliases
alias fastq-dump=$SOFT/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump

SRP=SRP153933
SOFTNAME=smcounter2
PAT="$1"

if [ -z "${java8}" ] ; then
    java8=$(which java)
fi

CURRTIME=$(date +%Y%m%d-%H%M%S)
SCRIPTDIR=$(dirname $(which $0))
ROOTDIR=${SCRIPTDIR}/../../
UVC=${ROOTDIR}/uvc1
V=$($UVC --version | head -n1 | awk '{print $(NF)}')

DATADIR="${SCRIPTDIR}/data/${SRP}/"
OUTDIR="${DATADIR}/${V}/${CURRTIME}"

mkdir -p "${OUTDIR}"

SRADIR="${ROOTDIR}/eval/smcounter2/data/${SRP}/sra/" 
mkdir -p "${SRADIR}"
FQDIR="${ROOTDIR}/eval/smcounter2/data/${SRP}/fq/"
mkdir -p "${FQDIR}"
for sra in $(ls "${SRADIR}" | grep ".sra$" | grep -P "${PAT}") ; do
    date # && time -p fastq-dump --gzip --split-files -O "${FQDIR}" "${SRADIR}/${sra}"
done

for fq1 in $(ls "${FQDIR}" | grep "_1.fastq.gz$" | grep -P "${PAT}") ; do
    fq2=${fq1/%_1.fastq.gz/_2.fastq.gz}
    out1=${fq1/%_1.fastq.gz/_1_barcoded.fastq.gz}
    out2=${fq2/%_2.fastq.gz/_2_barcoded.fastq.gz}
    echo date
    if [ $(echo "${PAT}" | grep run-extract-barcode | wc -l) -gt 0 ]; then
        echo time -p "${ROOTDIR}/scripts/extract-barcodes.py" "${FQDIR}/${fq2}" "${FQDIR}/${fq1}" "${FQDIR}/${out2}" "${FQDIR}/${out1}" 0 11
    fi
done | parallel -j 8

BAMDIR="${ROOTDIR}/eval/${SOFTNAME}/data/${SRP}/bam/"
mkdir -p ${BAMDIR}
for fq1 in $(ls "${FQDIR}" | grep "_1_barcoded.fastq.gz$" | grep -P "${PAT}") ; do
    fq2=${fq1/%_1_barcoded.fastq.gz/_2_barcoded.fastq.gz}
    bam=${fq1/%_1_barcoded.fastq.gz/_12.bam}
    if [ $(echo "${PAT}" | grep run-bwa | wc -l) -gt 0 ]; then
        date && time -p bwa mem -t 16 "${HG19}" "${FQDIR}/${fq1}" "${FQDIR}/${fq2}" | samtools view -bh1 | samtools sort -@2 -o "${BAMDIR}/${bam}"
        samtools index "${BAMDIR}/${bam}"
    fi
done

MAGERI="${ROOTDIR}/eval/mageri/mageri.jar"
PRESET1="${ROOTDIR}/eval/mageri/mageri-paper/processing/hd734_and_donors/preset.xml"
PRESET2="${ROOTDIR}/eval/mageri/mageri-paper/processing/hd734_and_donors/preset_1_1_1.xml"
cat "${PRESET1}" | sed 's/1.1.1-SNAPSHOT/1.1.1/g' > "${PRESET2}"
META="${ROOTDIR}/eval/mageri/mageri-paper/processing/hd734_and_donors/meta/"
BASELINEDIR=i"${ROOTDIR}/eval/${SOFTNAME}/data/${SRP}/baseline/"
mkdir -p "${BASELINEDIR}"
for fq1 in $(ls "${FQDIR}" | grep "_1_barcoded.fastq.gz$" | grep -P "${PAT}") ; do
    fq2=${fq1/%_1_barcoded.fastq.gz/_2_barcoded.fastq.gz}
    resultdir=${fq1/%_1_barcoded.fastq.gz/.resultdir}
    # must run smcounter2 in docker
    echo '''
    date && echo time -p ${java8} -Xmx32G -jar ${MAGERI} \
        --import-preset ${PRESET2} \
        --references    ${META}/refs.fa \
        --contigs       ${META}/contigs.txt \
        -M2             ${META}/primers.txt \
        --bed           ${META}/refs.bed \
        -R1             "${FQDIR}/${fq1}" \
        -R2             "${FQDIR}/${fq2}" \
        --threads       8 \
        -O              "${BASELINEDIR}/${resultdir}" ;
    '''
done

HOMVCF_NOHDR="${ROOTDIR}/eval/smcounter2/smcounter-v2-paper/N13532/misc/na24385.plus.na12878.hom.all.vcf"
HOMVCF="${ROOTDIR}/eval/smcounter2/smcounter-v2-paper/N13532/misc/na24385.plus.na12878.hom.all.hasheader.vcf"
HETVCF="${ROOTDIR}/eval/smcounter2/smcounter-v2-paper/N13532/misc/na12878.uniq.all.het.vcf"

bcftools view "${HETVCF}.gz" | grep "^#" > "${HOMVCF}"
cat "${HOMVCF_NOHDR}" | sed 's/^/chr/g' >> "${HOMVCF}"
bcftools view -Oz -o "${HOMVCF}.gz" "${HOMVCF}"
bcftools index -f -t "${HOMVCF}.gz"

TRUTHVCF="${HETVCF}" #"${ROOTDIR}/eval/smcounter2/smcounter-v2-paper/N13532/misc/truth.vcf"
TRUTHVCFGZ="${TRUTHVCF}.gz"

bcftools index -f -t "${TRUTHVCFGZ}"

CALLVCFGZ="${ROOTDIR}/eval/smcounter2/data/${SRP}/baseline/N13532-SRR7526729/NEB_S2.smCounter.anno.vcf.gz" 
bcftools index -f -t "${CALLVCFGZ}"

TRUTHBED="${ROOTDIR}/eval/smcounter2/smcounter-v2-paper/N13532/misc/N13532.hc.all.na12878.na24385.v3.3.2.bed" 

IEXPR="(bDP >= 10000 && cDP >= 3000 & bFA >= 0.001 && cFA >= 0.001 && QUAL>=65 && MAX(bSB1) <= 250 && MAX(cSB1) <= 250 && MAX(aDB) <= 150 && MAX(bPBL) <= 250 && MAX(bPBR) <= 250 && MAX(cPBL) <= 250 && MAX(cPBR) <= 250 && MAX(cVQ3) >= 100 && MAX(bMMB) <= 250 && MAX(cMMB) <= 250)" 
IMBA=200
IEXPR="(((1.0 - cFA - cFR)/ cFA) <= 0.75 && (cFA/bFA) <= 1.75 && MAX(bSB1) <= $IMBA && MAX(cSB1) <= $IMBA && MAX(bPBL) <= $IMBA && MAX(bPBR) <= $IMBA && MAX(cPBL) <= $IMBA && MAX(cPBR) <= $IMBA && MAX(bMMB) <= $IMBA && MAX(cMMB) <= $IMBA)"

VCFDIR="${ROOTDIR}/eval/${SOFTNAME}/data/${SRP}/vcf/"
for bam in $(ls ${BAMDIR} | grep ".bam$" | grep -P "${PAT}"); do
    vcf=${bam/%.bam/_uvc1.vcf.gz}
    if [ $(echo "${PAT}" | grep run-uvc | wc -l) -gt 0 ]; then
        date && time -p "${UVC}" -t 16 -f "${HG19}" -o "${VCFDIR}/${vcf}" -s "${bam}" "${BAMDIR}/${bam}" 2> "${VCFDIR}/${vcf}.stderr"
        date && time -p bcftools index -f -t "${VCFDIR}/${vcf}"
    fi
    ## 
    isecdir="${VCFDIR}/${vcf/%.vcf.gz/.germline2call_isec.dir}"
    bcftools isec -Oz -p "${isecdir}" "${HOMVCF}.gz" "${VCFDIR}/${vcf}" -R <(cat ${TRUTHBED} | sed 's/^/chr/g')  
    bcftools index -f -t "${isecdir}/0001.vcf.gz" 
    #bcftools view "${isecdir}/0001.vcf.gz" -Oz -o "${isecdir}/0001-flt.vcf.gz" -i "${IEXPR}" 
    bcftools view -i "QUAL>=40" "${isecdir}/0001.vcf.gz" \
        | bcftools filter -Ou -m+ -s moreThanTwoAlleles  -e "1 / cFA - 1 - cFR / cFA > 0.5" \
        | bcftools filter -Ou -m+ -s lowDupEfficiency    -e "(cFA/bFA) > 2" \
        | bcftools filter -Ou -m+ -s dupedStrandBias     -e "(bSBR[0:0] * bDP1[0:0] + bSBR[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBA" \
        | bcftools filter -Ou -m+ -s dedupStrandBias     -e "(cSBR[0:0] * cDP1[0:0] + cSBR[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBA" \
        | bcftools filter -Ou -m+ -s dupedPositionBias   -e "((bPBL[0:0] * bDP1[0:0] + bPBL[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBA) || ((bPBR[0:0] * bDP1[0:0] + bPBR[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > $IMBA)" \
        | bcftools filter -Ou -m+ -s dedupPositionBias   -e "((cPBL[0:0] * cDP1[0:0] + cPBL[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBA) || ((cPBR[0:0] * cDP1[0:0] + cPBR[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > $IMBA)" \
        | bcftools filter -Ou -m+ -s dupedManyMismatches -e "(bMMB[0:0] * bDP1[0:0] + bMMB[0:1] * bDP1[0:1]) / (bDP1[0:0] + bDP1[0:1]) > ($IMBA+100)" \
        | bcftools filter -Ou -m+ -s dedupManyMismatches -e "(cMMB[0:0] * cDP1[0:0] + cMMB[0:1] * cDP1[0:1]) / (cDP1[0:0] + cDP1[0:1]) > ($IMBA+100)" \
        | bcftools view   -Oz -o "${isecdir}/0001-flt.vcf.gz" 
    bcftools index -f -t "${isecdir}/0001-flt.vcf.gz"
    evalout="${bam/%.bam/_uvc1-flt.vcfeval.outdir}"
    ## 
    rm -r "${VCFDIR}/${evalout}" || true
    date && time -p "${java8}" -jar "${ROOTDIR}/eval/tools/rtg-tools-3.10.1/RTG.jar" vcfeval --squash-ploidy -f QUAL \
        -b "${TRUTHVCFGZ}" \
        -c "${isecdir}/0001-flt.vcf.gz" \
        -e <(cat ${TRUTHBED} | sed 's/^/chr/g') \
        -t "${ROOTDIR}/eval/datafiles/hg19_UCSC.sdf" \
        -o "${VCFDIR}/${evalout}" || true

    evalout="${bam/%.bam/_uvc1.vcfeval.outdir}"
    ## 
    rm -r "${VCFDIR}/${evalout}" || true
    date && time -p "${java8}" -jar "${ROOTDIR}/eval/tools/rtg-tools-3.10.1/RTG.jar" vcfeval --squash-ploidy --all-records -f QUAL \
        -b "${TRUTHVCFGZ}" \
        -c "${isecdir}/0001-flt.vcf.gz" \
        -e <(cat ${TRUTHBED} | sed 's/^/chr/g') \
        -t "${ROOTDIR}/eval/datafiles/hg19_UCSC.sdf" \
        -o "${VCFDIR}/${evalout}" || true
done

for sample in N13532-SRR7526729; do
    evalout=${sample}_smcounter2.vcfeval.outdir #"${bam/%.bam/_smcounter2.vcfeval.outdir}"
    absvcf="${ROOTDIR}/eval/smcounter2/data/SRP153933/baseline/${sample}/NEB_S2.smCounter.anno.vcf.gz"
    isecdir="${absvcf/%.vcf.gz/.germline2call_isec.dir}"
    bcftools isec -Oz -p "${isecdir}" "${HOMVCF}.gz" "${absvcf}" -R <(cat ${TRUTHBED} | sed 's/^/chr/g') 
    bcftools index -f -t "${isecdir}/0001.vcf.gz"

    rm -r "${VCFDIR}/${evalout}" || true 
    date && time -p "${java8}" -jar "${ROOTDIR}/eval/tools/rtg-tools-3.10.1/RTG.jar" vcfeval --squash-ploidy -f QUAL \
        -b "${TRUTHVCFGZ}" \
        -c "${isecdir}/0001.vcf.gz" \
        -e <(cat ${TRUTHBED} | sed 's/^/chr/g') \
        -t "${ROOTDIR}/eval/datafiles/hg19_UCSC.sdf" \
        -o "${VCFDIR}/${evalout}" || true
    date && time -p $java8 -jar /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc//eval/tools/rtg-tools-3.10.1/RTG.jar rocplot \
        "${VCFDIR}/"*SRR7526729*.vcfeval.outdir/*"snp_roc.tsv.gz" --scores --zoom 200,300 \
        --svg "${VCFDIR}/${evalout}/${evalout/%_smcounter2.vcfeval.outdir/_all_methods_all_vars_rocplot.svg}" \
        --png "${VCFDIR}/${evalout}/${evalout/%_smcounter2.vcfeval.outdir/_all_methods_all_vars_rocplot.png}"
    date && time -p $java8 -jar /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc//eval/tools/rtg-tools-3.10.1/RTG.jar rocplot \
        "${VCFDIR}/"*SRR7526729*.vcfeval.outdir/*"snp_roc.tsv.gz" --scores --zoom 20,300 \
        --svg "${VCFDIR}/${evalout}/${evalout/%_smcounter2.vcfeval.outdir/_all_methods_all_vars_rocplot_zoom.svg}" \
        --png "${VCFDIR}/${evalout}/${evalout/%_smcounter2.vcfeval.outdir/_all_methods_all_vars_rocplot_zoom.png}"
done

echo UNUSED '''
for bam in $(ls ${BAMDIR} | grep ".bam$" | grep "${PAT}"); do
    date && time -p "${java8}" -jar "${ROOTDIR}/eval/tools/rtg-tools-3.10.1/RTG.jar" vcfeval --squash-ploidy -f QUAL \
    -b /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc/eval/smcounter2/smcounter-v2-paper/N13532/misc/na12878.uniq.all.het.vcf.gz \
    -c /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc/eval/smcounter2/data/baseline/N13532-SRR7526729/NEB_S2.smCounter.anno.vcf.gz \
    -e <(cat /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc/eval/smcounter2/smcounter-v2-paper/N13532/misc/N13532.hc.all.na12878.na24385.v3.3.2.bed | sed 's/^/chr/g') \
    -t /bionfsdate/ctDNA/experiment/zhaoxiaofei/database/genome/hg19_UCSC.sdf \
    -o /hongshan/tmp/vcfeval-test5 \
done
'''
