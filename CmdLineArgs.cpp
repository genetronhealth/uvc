#ifndef CmdLineArgs_INCLUDED
#define CmdLineArgs_INCLUDED

#include "CmdLineArgs.hpp"
#include "version.h"

#include "htslib/sam.h"

#include <fstream>

#include <float.h>

SequencingPlatform 
CommandLineArgs::selfUpdateByPlatform() {
    SequencingPlatform inferred_sequencing_platform = this->sequencing_platform;
    if (SEQUENCING_PLATFORM_AUTO == this->sequencing_platform || SEQUENCING_PLATFORM_OTHER == this->sequencing_platform) {
        samFile *sam_infile = sam_open(this->bam_input_fname.c_str(), "r"); // AlignmentFile(samfname, "rb")
        if (NULL == sam_infile) {
            fprintf(stderr, "Failed to open the file %s", this->bam_input_fname.c_str());
            exit(-32);
        }
        bam_hdr_t * samheader = sam_hdr_read(sam_infile);
        bam1_t *b = bam_init1();
        unsigned int countPE = 0;
        unsigned int countSE = 0;
        std::vector<unsigned int> qlens;
        qlens.reserve(500+1);
        qlens.push_back(150);
        while (sam_read1(sam_infile, samheader, b) >= 0 && (countPE + countSE) < 500) {
            if (b->core.flag & 0x1) {
                countPE++;
            } else {
                countSE++;
            }
            qlens.push_back(b->core.l_qseq);
        }
        std::sort(qlens.begin(), qlens.end());
        if (0 == this->central_readlen) { this->central_readlen = qlens.at(qlens.size()/2); }
        bam_destroy1(b);
        bam_hdr_destroy(samheader);
        sam_close(sam_infile);
        if (0 < countPE) {
            inferred_sequencing_platform = SEQUENCING_PLATFORM_ILLUMINA;
        } else {
            inferred_sequencing_platform = SEQUENCING_PLATFORM_IONTORRENT;
        }
    }
    if (SEQUENCING_PLATFORM_IONTORRENT == inferred_sequencing_platform && SEQUENCING_PLATFORM_OTHER != this->sequencing_platform) {
        bq_phred_added_indel += 0;
        bq_phred_added_misma += 4; // it was 6
        minABQ_pcr_snv += 0;
        minABQ_pcr_indel += 0;
        minABQ_cap_snv += 0;
        minABQ_cap_indel += 0;
        if (0 == highqual_thres_indel) { highqual_thres_indel = highqual_thres_snv - 4; }
    }
    if (SEQUENCING_PLATFORM_ILLUMINA == inferred_sequencing_platform && SEQUENCING_PLATFORM_OTHER != this->sequencing_platform) {
        bq_phred_added_indel += 10; // 0; // 6; //10;
        bq_phred_added_misma += 0;
        minABQ_pcr_snv += 25;
        minABQ_pcr_indel += 15; // 18;
        minABQ_cap_snv += 20;
        minABQ_cap_indel += 0; // 5; // 13;
        if (0 == highqual_thres_indel) { highqual_thres_indel = highqual_thres_snv + 6; }
    }
    return inferred_sequencing_platform;
}

void 
check_file_exist(const std::string & fname, const std::string ftype) {
    std::ifstream ifile(fname.c_str());
    if (!(bool)ifile) {
        std::cerr << "The file " << fname << " of type (" << ftype << ") does not exist." << std::endl;
        exit(-4);
    }
}

const std::string 
stringvec_to_descstring(const std::vector<std::string> & v) {
    std::string ret;
    for (unsigned int i = 0; i < v.size(); i++) {
        ret += std::to_string(i) + " : " + v[i] + ". ";
    }
    return ret;
};

int
CommandLineArgs::initFromArgCV(int & parsing_result_flag, SequencingPlatform & inferred_sequencing_platform, int argc, const char *const* argv) {
    parsing_result_flag = -1;
    auto version_cb = [](int count){
        std::cout << "uvc-" << VERSION << std::endl;
        exit(0);
    };
    CLI::App app{(std::string("Unified Variant Caller (UVC) version ") + VERSION_DETAIL)};
    app.add_flag_function("-v,--version", version_cb,   "Show the version of this program (打印此软件版本号).");    
    app.add_option("inputBAM",       bam_input_fname,   "Input coordinate-sorted BAM file that is supposed to contain raw reads (按照位置排序的有原始reads的BAM文件).")->required()->check(CLI::ExistingFile);
    app.add_option("--output-bpRES", vcf_output_fname,  "Output bgzipped VCF file in the format of base-pair resolution. "
            "This option is deprecated, please use -A instead. (每个位置都输出检测信息的VCF输出文件，文件有可能非常大。已经废除，建议使用-A参数)");
    app.add_option("-o,--output",    vcf_out_pass_fname,"Output bgzipped VCF file in the format of blocked gVCF. (区块gVCF输出文件).", true);
    app.add_flag("-A,--All-out",     should_let_all_pass,"All possible alleles including REF allele at each position is in <--output>. The output is base-pair resolution VCF.(在--output文件输出每个位置的所有等位基因形)");
    app.add_option("--tumor-vcf",    vcf_tumor_fname,   "Block-gzipped VCF file of the tumor sample as input to the BAM file of normal sample. If specified/unspecified, then the input BAM is considered to be from normal/tumor. (肿瘤样本的VCF文件，用于输入。如果提供/不提供则认为BAM文件是normal/tumor的).", true);
    app.add_option("-f,--fasta",     fasta_ref_fname,   "Input reference fasta file, where the special value NA means not available (FASTA参考基因组).");
    app.add_option("--targets",     tier1_target_region,"Input region string (for example, chr1:2-3 or chr1) which is optional and globally delimits the genomic region to call variants from, "
            "if not provided then call variants from all the sequenced genomic regions. This parameter does not affect memory consumption. "
            "(基因组区域字符串，可有可无，如果没提供则从reads覆盖信息自动计算出BED文件，此参数不影响内存消耗).");
    app.add_option("-R,--regions-file", bed_region_fname,  "Input BED region file which is optional and delimits the genomic regions to call variants from, "
            "if not provided then call variants from all the sequenced genomic regions. "
            "This overrides the -t parameter. Please ensure that each region here is small enough to prevent memory exhaustion. "
            "Typically, each region covers one exon. "
            "(BED文件，可有可无，如果没提供则从reads覆盖信息自动计算出BED文件，此参数导致-t参数失效。"
            "如果BED文件中任何区域太大可能会耗尽内存，每个区域应该大概对应一个外显子).")->check(CLI::ExistingFile);
    app.add_option("-s,--sample",    sample_name,       "Sample name which is optional (样本名称，可有可无).");
    // app.add_option("--primers",      tsv_primer_fname,  "primer files")->check(CLI::ExistingFile);
    app.add_option("--bed-outfname", bed_out_fname,     "The BED file to which genomic-region information will be written. Empty string means no output is written. This BED file can be generated by tumor and used by normal. (BED输出文件).", true);
    app.add_option("--bed-infname",  bed_in_fname,      "The BED file from which genomic-region information is read. Empty string means no input is read. This BED file can be generated by tumor and used by normal. This param overrides the -R param (BED输入文件).", true);
    
    app.add_option("-q,--vqual",     vqual,             "Every variant that satisfies this minimum variant quality is in the -o file. (如果变异质量低于此值，则不输出到-o文件).", true);
    app.add_option("--v-min-altdp",  vad,               "Every variant that satisfies this minimum allele depth is in the -o file. (如果变异质量低于此值，则不输出到-o文件).", true);

    app.add_option("-d,--min-depth", min_depth_thres,   "Minimum depth below which results are fitlered out and therefore not in the output VCF (如果低于此原始深度则在VCF不输入任何结果).", true);
    app.add_option("-D,--min-altdp", min_altdp_thres,   "Minimum depth of ALT below which results are filtered out (如果ALT深度低于此数则不输出结果).", true);
    app.add_option("-t,--threads",   max_cpu_num,       "Number of cpu cores or equivalently threads to use (使用CPU线程的数量).", true);
    app.add_option("--is-tumor-format-retrieved", is_tumor_format_retrieved, "Boolean (0: false, 1: true) indicating if the format from the tumor VCF should be retrieved in the tumor-normal comparison. This boolean is only useful if tumor VCF is provided. Notice that a true value for this boolean disables the generation of genomic block so that the output is no longer gvcf.", true);
    
    app.add_option("--alnlen",       min_aln_len,       "Minimum alignment length below which the alignment is filtered out (如果比对长度低于比值则过滤掉一行的比对结果).", true);
    app.add_option("--mapqual",      min_mapqual,       "Minimum mapping  quality below which the alignment is filtered out (如果比对质量低于此值则过滤掉一行的比对结果).", true);
    app.add_option("--fixedthresBQ", fixedthresBQ,      "Base quality cutoff. This parameter is only for generating statistics and therefore does not affect variant quality (碱基质量阈值，只用于统计，不影响变异质量).", true);
    app.add_option("--nogap-phred",  nogap_phred,       "Base quality for the symbol denoting non-InDel, can solve disconcordant alignment problem in the overlap between R1 and R2, This parameter is now obsolete because the current R1-R2 merge considers symbols to be merged so that InDels dominate over non-Indels. (合并 R1 和 R2 有可能遇到没有InDel和有InDel这种不一致情况，因此没有InDel的质量有这个上限，目前已废除).", true);
    app.add_option("--vc-stats-fname", vc_stats_fname  ,"Output TSV file containing variant-call statistics. The default is the standard error stream. (含有检测变异统计信息的TSV文件，默认为标准错误输出).", true);
     
    app.add_option("--uni-bias-thres", uni_bias_thres,  "Unified-bias threshold for generating the filter strings in FORMAT/FT. This parameter is only for generating statistics and therefore does not affect variant quality. Downstream hard filtering with FORMAT/FT is possible (统一偏好性的阈值，用于生成FORMAT/FT信息，只用于统计，不影响变异质量，FORMAT/FT可用于下游硬过滤). ", true);
    app.add_option("--uni-bias-r-max", uni_bias_r_max,  "Maximum unified-bias threshold used for reducing variant read support. This parameter does affect variant quality. (统一偏好性的最大值，用于减少变异支持，会影响变异质量). ", true);
    //app.add_option("--diffVAQfrac",    diffVAQfrac,     "Experimental real-numbered parameter that should be set to either zero or one (实验性的实数参数，理论值要么是零要么是一). ", true);
    
    app.add_option("--highqual-thres-snv",          highqual_thres_snv,
            "The SNV quality threshold above which the family quality is considered to be high", true);
    app.add_option("--highqual-thres-indel",        highqual_thres_indel,
            "The InDel quality threshold above which the family quality is considered to be high, "
            "zero means auto infer to highqual-thres-snv + 6 for Illumina/BGI and -4 for IonTorrent)", true);
    app.add_option("--highqual-min-ratio",          highqual_min_ratio,
            "The mininum ratio of the raw non-deduplicated read depth to the deduplicated read family depth to trigger tumor-normal comparison with high quality families only", true);
    //app.add_option("--highqual-min-vardep",         highqual_min_vardep,
    //        "the mininum number of families suporting the variant in the tumor sample to trigger tumor-normal comparison with families of duplicated reads only", true);
    //app.add_option("--highqual-min-totdep",         highqual_min_vardep,
    //        "the mininum number of familiess supporting any allele in the tumor sample to trigger tumor-normal comparison with families of duplicated reads only", true);
    
    app.add_option("--phred-frag-indel-ext",        phred_max_frag_indel_ext,
            "Maximum phred score fo the indel of one additional base (excluding the one base required for opening indel), capped at two additional bases", true);
    app.add_option("--phred-frag-indel-basemax",    phred_max_frag_indel_basemax,
            "Maximum phred score fo the opening of an indel (including the one base required for opening indel)", true);
    app.add_option("--phred-sscs-transition-CG-TA", phred_max_sscs_transition_CG_TA, 
            "Maximum phred score for single-strand consensus sequences (SSCSs) for C:G > T:A transition", true);
    app.add_option("--phred-sscs-transition-TA-CG", phred_max_sscs_transition_TA_CG, 
            "Maximum phred score for single-strand consensus sequences (SSCSs) for T:A > C:G transition", true);
    app.add_option("--phred-sscs-transversion-any", phred_max_sscs_transversion_any, 
            "Maximum phred score for single-strand consensus sequences (SSCSs) for any transversion", true);
    app.add_option("--phred-sscs-indel-open",       phred_max_sscs_indel_open, 
            "Maximum phred score for single-strand consensus sequences (SSCSs) for the opening of indel gap (the opening includes the insertion/deletion of one base)", true);
    app.add_option("--phred-sscs-indel-ext" ,       phred_max_sscs_indel_ext, 
            "Maximum phred score for single-strand consensus sequences (SSCSs) for the extension of indel gap (excluding the extension of one base) per base", true);
    app.add_option("--phred-dscs-all",              phred_max_dscs_all, 
            "Maximum phred score for double-strand consensus sequences (DSCSs) for all types of mutations", true);
    app.add_option("--phred-pow-sscs-origin",       phred_pow_sscs_origin, 
            "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment", true);
    app.add_option("--phred-pow-dscs-origin",       phred_pow_dscs_origin, 
            "The phred-score that is subtracted from phred-dscs to get the power-law quality adjustment", true);
    
    app.add_option("--ess-georatio-dedup-cap", ess_georatio_dedup_cap, 
                   "Geometric common ratio of the increase in the observed number of deduped reads as a function of the effective number of deduped reads (effective sample size) for capture-based assays.", true);
    app.add_option("--ess-georatio-dedup-pcr", ess_georatio_dedup_pcr, 
                   "Geometric common ratio of the increase in the observed number of deduped reads as a function of the effective number of deduped reads (effective sample size) for PCR-based assays.", true);
    app.add_option("--ess-georatio-duped-pcr", ess_georatio_duped_pcr, 
                   "Geometric common ratio of the increase in the observed number of duped reads as a function of the effective number of duped reads (effective sample size) for reads derived from one template molecule.", true);
 
    app.add_option("--minABQ-pcr-snv",   minABQ_pcr_snv,   "Minimum average base quality below which variant quality is capped to average base quality for PCR assay and SNVs, "
                   "recommend 25 for Illumina and 0 for IonTorrent "
                   " (如果位点平均碱基质量低于此值则变异质量不会超过平均碱基质量，建议对Illumina用25并且对IonTorrent用0).", true); 
    app.add_option("--minABQ-pcr-indel", minABQ_pcr_indel, "Minimum average base quality below which variant quality is capped to average base quality for PCR assay and InDels, "
                   "recommend 18 for Illumina and 0 for IonTorrent "
                   " (如果位点平均碱基质量低于此值则变异质量不会超过平均碱基质量，建议对Illumina用25并且对IonTorrent用0).", true); 
    app.add_option("--minABQ-cap-snv",   minABQ_cap_snv,   "Minimum average base quality below which variant quality is capped to average base quality for capture assay and SNVs."
                   "recommend 20 for Illumina and 0 for IonTorrent "
                   " (如果位点平均碱基质 量低于此值则变异质量不会超过平均碱基质量(捕获试验)，建议对Illumina用25并且对IonTorrent用0).", true); 
    app.add_option("--minABQ-cap-indel", minABQ_cap_indel, "Minimum average base quality below which variant quality is capped to average base quality for capture assay and InDels."
                   "recommend 13 for Illumina and 0 for IonTorrent "
                   " (如果位点平均碱基质量低于此值则变异质量不会超过平均碱基质量(捕获试验)，建议对Illumina用25并且对IonTorrent用0).", true); 
    
    app.add_option("--minMQ1", minMQ1, "Minimum root-mean-square (RMS) mapping quality (MQ) of non-dedupped raw reads below which variant quality is capped to this RMS MQ.", true);
    app.add_option("--maxMQ" , maxMQ , "Maximum mapping quality (MQ) of the aligned reads, highly recommended to be the 60 from BWA.", true);
    app.add_option("--min-edge-dist", min_edge_dist, "Minimum average number of bases to the left and right aligned positions below which variant quality is capped to 4 times this average number.", true);
    
    app.add_option("--central-readlen",      central_readlen,      "Central (median) value for read lengths, 0 means estimate from the data.");
    app.add_option("--bq-phred-added-misma", bq_phred_added_misma, "Additional base-quality phred score added to match and mismatch, recommend 4 for IonTorrent from Life Technologies.");
    app.add_option("--bq-phred-added-indel", bq_phred_added_indel, "Additional base-quality phred score added to indel and no-indel, recommend 10 for Illumina and BGI.");
    
    app.add_option("--phred-homref-gt",       phred_germline_polymorphism,"Phred-scaled prior probability of germline polymorphism event at a loci.", true);
    app.add_option("--phred-triallelic-indel",phred_triallelic_indel,     "Phred-scaled prior probability that two or more forms of InDel occur at a loci.", true);
 
    app.add_option("--any-mul-contam-frac",  any_mul_contam_frac,         "Multiplicative contamination rate for the fraction of reads generated by any source of contamination in any sample", true);
    app.add_option("--t2n-mul-contam-frac",  t2n_mul_contam_frac,         "Multiplicative contamination rate for the fraction of tumor reads in the normal", true); 
    app.add_option("--t2n-add-contam-frac",  t2n_add_contam_frac,         "Additive contamination rate for the fraction of tumor reads in the normal", true);
    app.add_option("--t2n-add-contam-transfrac",t2n_add_contam_transfrac, "Additive contamination rate for the fraction of tumor reads in the normal by read transfer", true); 
    
    app.add_option("--phred-frac-indel-error-before-barcode-labeling", phred_frac_indel_error_before_barcode_labeling, 
                                                                       "PHRED-scaled fraction of InDel errors that occurred before the attachment of UMI single-strand barcodes.", true); 
    app.add_option("--baq-per-aligned-base", baq_per_aligned_base,     "PHRED-scaled base alignment quality (BAQ) for each additional base between the variant and the read edge.", true); 
   
    app.add_option("--is-somatic-snv-filtered-by-any-nonref-germline-snv", is_somatic_snv_filtered_by_any_nonref_germline_snv,
                   "Set to 0 (zero, false) if reject any nonref germline and to 1 (one, true) if only reject the specific ALT germline for SNV candidate.", true); 
    app.add_option("--is-somatic-indel-filtered-by-any-nonref-germline-indel", is_somatic_indel_filtered_by_any_nonref_germline_indel,
                   "Set to 0 (zero, false) if reject any nonref germline and to 1 (one, true) if only reject the specific ALT germline for InDel candidate.", true); 
    app.add_option("--illumina-BQ-pow2-div-coef", illumina_BQ_pow2_div_coef,
                   "The square of Illumina root-mean-square base quality divided by this number is the maximum possible TLOD part of variant quality.", true); 
    app.add_option("--varqual-per-mapqual", varqual_per_mapqual,
                   "The root-mean square mapping quality multiplied by this number is the maximum possible TLOD part of variant quality.", true); 
    
    app.add_option("--powlaw-exponent",     powlaw_exponent,
                   "Exponent of the power-law relationship between error probability and allele-fraction deviation, strongly recommended to use the default value.", true);     
    app.add_option("--powlaw-anyvar-base", powlaw_anyvar_base,
                   "Error probability at allele-fraction of 1 or 100%, storngly recommended to use the default value", true); 
    app.add_option("--syserr-maxqual", syserr_maxqual,
                   "PHRED-scaled probability that a candidate of systematic error is actually non-systematic.", true); 
    app.add_option("--syserr-norm-devqual", syserr_norm_devqual,
                   "PHRED-scaled likelihood that the observed allele fraction additively deviates from the expected allele fraction by a multiplicative factor of two.", true); 
    
    app.add_flag("--Should-add-note",        should_add_note,             "Flag indicating if the program generates more detail (can be used for debugging) in the VCF result file.");
    
    app.add_option("--disable-dup-read-merge", disable_dup_read_merge,      "Boolean (0: false, 1: true) indicating if the program should disable the merge of duplicate reads.", true);
    app.add_option("--enable-dup-read-vqual",  enable_dup_read_vqual,       "Boolean (0: false, 1: true) indicating if the program should enable the use of raw non-dedupped reads in the calculation of variant quality.", true);
    unsigned int assay_type_uint = (unsigned int)assay_type;
    unsigned int molecule_tag_uint = (unsigned int)molecule_tag;
    unsigned int sequencing_platform_uint = (unsigned int)sequencing_platform;
    unsigned int pair_end_merge_uint = (unsigned int)pair_end_merge;
    app.add_option("--assay-type",           assay_type_uint,           "Assay type. " + stringvec_to_descstring(ASSAY_TYPE_TO_MSG), true);
    app.add_option("--molecule-tag",         molecule_tag_uint,         "Molecule tag. " + stringvec_to_descstring(MOLECULE_TAG_TO_MSG), true);
    app.add_option("--sequencing-platform",  sequencing_platform_uint,  "Sequencing platform. " + stringvec_to_descstring(SEQUENCING_PLATFORM_TO_MSG), true);
    app.add_option("--pair-end-merge",       pair_end_merge_uint,       "Mode for the merge of R1 and R2 in a read pair. " + stringvec_to_descstring(PAIR_END_MERGE_TO_MSG), true);
    
    app.add_flag("--Disable-duplex",         disable_duplex,            "Flag indicating if the program should disable the merge of two SSCSs (single-strand-consensus-sequences) into a DSCS (double strand consensus sequence). The UMI of the duplex tag should be in the form of <alpha>+<beta>.");
    
    app.add_option("--dedup-center-mult",    dedup_center_mult,         "Exponential decay per additional base of distance." , true);
    app.add_option("--dedup-amplicon-count-to-surrcount-frac", dedup_amplicon_count_to_surrcount_frac, 
                   "centroidCount/surroundingCount of reads ending at a position above which the assay is inferred to be amplicon. Assay type can be override on command-line", true);
    app.add_option("--dedup-yes-umi-2ends-peak-frac",    dedup_yes_umi_2ends_peak_frac, 
                   "If the ratio of the numbers of reads ending on the two ends with UMI is above this value, then UMI is assumed to be single-end." , true);
    app.add_option("--dedup-non-umi-2ends-peak-frac",    dedup_non_umi_2ends_peak_frac, 
                   "Same as above except for L8G" , true); 
    
    app.callback([&]() {
        assay_type = (AssayType)assay_type_uint;
        molecule_tag_uint = (MoleculeTag)molecule_tag_uint;
        sequencing_platform = (SequencingPlatform)sequencing_platform_uint;
        pair_end_merge = (PairEndMerge)pair_end_merge_uint;
        
        check_file_exist(bam_input_fname, "BAM");
        check_file_exist(bam_input_fname + ".bai", "BAM index");
        if (fasta_ref_fname.compare(std::string("NA")) != 0) {
            check_file_exist(fasta_ref_fname, "FASTA");
            check_file_exist(fasta_ref_fname + ".fai", "FASTA index");
        } else {
            fasta_ref_fname = "";
        }
        inferred_sequencing_platform = this->selfUpdateByPlatform();
        if (vcf_tumor_fname != NOT_PROVIDED) {
            // vqual -= (double)10; // maybe useful but not now
        }
        if (bed_in_fname != NOT_PROVIDED) {
            bed_region_fname = bed_in_fname;
        }
        parsing_result_flag = 0;
        if (t2n_add_contam_transfrac < DBL_EPSILON) { t2n_add_contam_transfrac = DBL_EPSILON; }
    });
    
    CLI11_PARSE(app, argc, argv);
    return 0;
}

#endif

