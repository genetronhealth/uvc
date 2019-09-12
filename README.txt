
# Example usage, where the output file is /hongshan/tmp/subset/22/uvcTN-results/NA12878-NA24385-mixture_N_uvc1.vcf.gz
# -R /clinical/BEDfiles/TISSUE/WES/wes_509.merged.bed means using this bed
# file for analyzing both tumor and normal bam data
#  --tumor-params -t 16 -q 5 means using 16 CPUs for analyzing  tumor bam data
#  and using variant quality (QUAL) threshold of 5 for filtering out variant
#  candidates.
# --normal-params -t 12      means using 12 CPUs for analyzing normal bam data

./uvcTN.sh /biocluster/data/bioexec/database/genome/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc/eval/tn/data/bam//subset/24385-12878-30-200_RX_001-22.bam /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc/eval/tn/data/bam//subset/24385-200_AH5G7WCCXX_S4_L004_RX_001-22.bam /hongshan/tmp/subset/22/uvcTN-results/ NA12878-NA24385-mixture \
-R /clinical/BEDfiles/TISSUE/WES/wes_509.merged.bed \
 --tumor-params -t 16 -q 5 \
--normal-params -t 12

# 变异检测结果在后缀是 _N_uvc1.vcf.gz 文件中（选VCF/INFO 是 PASS的结果即可）
# 如果只有 tumor 结果，用 uvc1 检测即可。
# 如果read名称是readname#UMI的格式，那么此软件会自动识别UMI（分子标签）
# 这个脚本有五个必须输的参数，剩下的参数可有可无，如果有的话会修改 uvc1 默认参数。
# 
# 如果 bcftools 报错，请使用 
# BCFTOOLS_DIR=/biocluster/data/bioexec/software/bcftools-1.3.1
# export PATH="${BCFTOOLS_DIR}:${PATH}"
# 
# 如果程序运行报错并且提示没有找到 libhts.so*, 请执行一下命令：
# HTSLIB_DIR=/biocluster/data/bioexec/software/htslib-1.6
# export LD_LIBRARY_PATH="${HTSLIB_DIR}:${LD_LIBRARY_PATH}"

# 总共有四个可执行文件：
uvc1：运行速度最快的release版本，使用多线程。因为对运行速度做了极大的优化，所以如果运行出错不会产生有用的错误信息。
uvc.mt.out：多线程debug版本。如果release版本出问题，请使用这个多线程debug版本。然后会产生报错信息。报错信息用于开发人员debug。
uvc.st.out：单线程debug版本。如果多线程debug版本出问题，请使用这个单线程debug版本。然后会产生报错信息。报错信息用于开发人员debug。
uvcTN.sh: 做tumor-normal配对分析。

