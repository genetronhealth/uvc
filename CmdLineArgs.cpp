#ifndef CmdLineArgs_INCLUDED
#define CmdLineArgs_INCLUDED

#include "CmdLineArgs.hpp"
#include "version.h"

#include "htslib/sam.h"

#include <fstream>

#include <float.h>

#define ADD_OPTDEF(app, k, v, msg) (app.add_option(k, v, msg, true))

bool
is_bitflag_checked(uint32_t bitflag_InDel_penal_t_UMI_n_UMI, bool is_InDel, bool is_penal, bool is_t_UMI, bool is_n_UMI) {
    uint32_t bitflag_index = 0;
    if (!is_InDel) {
        bitflag_index += 8;
    }
    if (!is_penal) {
        bitflag_index += 4;
    } 
    if (!is_t_UMI) {
        bitflag_index += 2;
    }
    if (!is_n_UMI) {
        bitflag_index += 1;
    }
    return ((0x1L << bitflag_index) & bitflag_InDel_penal_t_UMI_n_UMI);
}

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
        unsigned int q30_n_fail_bases = 0;
        unsigned int q30_n_pass_bases = 0;
        while (sam_read1(sam_infile, samheader, b) >= 0 && (countPE + countSE) < 500) {
            if (b->core.flag & 0x1) {
                countPE++;
            } else {
                countSE++;
            }
            qlens.push_back(b->core.l_qseq);
            for (int qpos = 0; qpos < b->core.l_qseq; qpos++) {
                unsigned int bq = (bam_get_qual((b))[(qpos)]);
                if (bq < 30) {
                    q30_n_fail_bases++;
                } else {
                    q30_n_pass_bases++;
                }
            }
        }
        std::sort(qlens.begin(), qlens.end());
        if (0 == this->central_readlen) { this->central_readlen = qlens.at(qlens.size()/2); }
        bam_destroy1(b);
        bam_hdr_destroy(samheader);
        sam_close(sam_infile);
        if ((0 < countPE) || (q30_n_fail_bases * 3 < q30_n_pass_bases)) {
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
        bq_phred_added_indel += 17; // 10; // 0; // 6; //10;
        bq_phred_added_misma += 0;
        minABQ_pcr_snv += 175; // 19; // 25;
        minABQ_pcr_indel += 90; // 10; // 15; // 18;
        minABQ_cap_snv += 175; // 19; // 23;
        minABQ_cap_indel += 90; // 10; // 5; // 13;
        if (0 == highqual_thres_indel) { highqual_thres_indel = highqual_thres_snv + 6; }
    }
    return inferred_sequencing_platform;
}

void 
check_file_exist(const std::string & fname, const std::string ftype) {
    std::ifstream ifile(fname.c_str());
    if (!(bool)ifile) {
        std::cerr << "The file " << fname << " of type (" << ftype << ") does not exist. " << std::endl;
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
    CLI::App app{(std::string("UVC version ") + VERSION_DETAIL)};
    
    // more frequently used options
    app.add_flag_function("-v,--version", version_cb, "Show the version of this program. ");
    ADD_OPTDEF(app, 
        "inputBAM", 
           bam_input_fname, 
        "Input coordinate-sorted BAM file that is supposed to contain raw reads. ")->required()->check(CLI::ExistingFile);
    ADD_OPTDEF(app, 
        "-f,--fasta", 
           fasta_ref_fname, 
        "Input reference fasta file, where the special value NA means not available. ");
    
    //ADD_OPTDEF(app, "--output-bpRES", vcf_output_fname,  "Output bgzipped VCF file in the format of base-pair resolution. "
    //        "This option is deprecated, please use -A instead. ");
    ADD_OPTDEF(app, 
        "-o,--output", 
           vcf_out_pass_fname, 
        "Output bgzipped VCF file in the format of blocked gVCF. ");
    ADD_OPTDEF(app, 
        "-R,--regions-file", 
           bed_region_fname, 
        "Input BED region file which is optional and delimits the genomic regions to call variants from, "
        "if not provided then call variants from all the sequenced genomic regions. "
        "This overrides the -t parameter. Please ensure that each region here is small enough to prevent memory exhaustion. "
        "Typically, each region covers one exon. ")->check(CLI::ExistingFile);
    ADD_OPTDEF(app, 
        "-s,--sample", 
           sample_name, 
        "Sample name which is optional. ");
    ADD_OPTDEF(app, 
        "--outvar-flag", 
           outvar_flag, 
        std::string("Output-variant flag in bits specifying which type of variants are in the VCF output. ") +
        "The " + std::to_string(OUTVAR_GERMLINE) + " bit indicates germline variant. " +
        // "The " + std::to_string(OUTVAR_HOMREF) + " bit indicates homozygous reference allele (<NO_SNV> and/or <NO_INDEL>). " +
        "The " + std::to_string(OUTVAR_SOMATIC) + " bit indicates somatic variant. " +
        "The " + std::to_string(OUTVAR_ANY) + " bit indicates variant of any origin. " +
        "The " + std::to_string(OUTVAR_NONREF) + " bit indicates non-ref symbolic region indicating the absence of variants. ");

    ADD_OPTDEF(app, 
        "-q,--vqual", 
           vqual, 
        "Every variant that satisfies this minimum variant quality is in the <--output> file. ");
    app.add_flag(
        "-A,--all-out", 
           should_output_all, 
        "All possible alleles including REF allele at each position is in <--output>. The output is base-pair resolution VCF. ");
    app.add_flag(
        "--all-germline-out", 
           should_output_all_germline,
        "All possible alleles including REF allele at each position is in <--output> for germline variants. ");
    ADD_OPTDEF(app, 
        "-d,--min-depth", 
           min_depth_thres, 
        "Minimum depth below which results are fitlered out and therefore not in the output VCF. ");
    ADD_OPTDEF(app, 
        "-a,--min-altdp", 
           min_altdp_thres, 
        "Minimum depth of ALT below which results are filtered out. ");
    ADD_OPTDEF(app, 
        "-t,--threads", 
           max_cpu_num, 
        "Number of cpu cores or equivalently threads to use. ");
    
    ADD_OPTDEF(app,
        "--primerlen", 
           primerlen, 
        "Number of bases from each end of each insert that are part of the primers (for only amplicon assay). ");

    ADD_OPTDEF(app, 
        "--tumor-vcf", 
           vcf_tumor_fname, 
        "Block-gzipped VCF file of the tumor sample as input to the BAM file of normal sample. "
        "If specified/unspecified, then the input BAM is considered to be from normal/tumor. ");
    ADD_OPTDEF(app, 
        "--targets", 
           tier1_target_region, 
        "Optional input region string (for example, chr1:2-3 or chr1) "
        "which globally delimits the genomic region to call variants from, "
        "if not provided then call variants from all the sequenced genomic regions. This parameter does not affect memory consumption. ");
    
    ADD_OPTDEF(app, 
        "--somaticGT", 
           somaticGT, 
        "Boolean (0: false, 1: true) indicating if the GT in the output VCF file refers to the genotype of the tumor cells "
        "in a heterogeneous mixture of tumor cells and normal cells. "
        "Please set to zero to call germline variants in a homogeneous mixture of cells. ");
    ADD_OPTDEF(app, 
        "--is-tumor-format-retrieved", 
           is_tumor_format_retrieved, 
        "Boolean (0: false, 1: true) indicating if the format from the tumor VCF should be retrieved in the tumor-normal comparison. "
        "This boolean is only useful if tumor VCF is provided. "
        "Notice that a true value for this boolean disables the generation of genomic block so that the output is no longer gvcf. ");
    ADD_OPTDEF(app, 
        "--alnlen", 
           min_aln_len,
        "Minimum alignment length below which the alignment is filtered out. ");
    ADD_OPTDEF(app, 
        "--mapqual", 
           min_mapqual,
        "Minimum mapping quality below which the alignment is filtered out. ");
    ADD_OPTDEF(app, 
        "--v-min-alldp", 
           vdp,
        "Every variant at each locus with at least this total depth is in the -o file. ");
    ADD_OPTDEF(app, 
        "--v-min-altdp", 
           vad,
        "Every variant that satisfies this minimum allele depth is in the -o file. ");
    
    ADD_OPTDEF(app, 
        "--dedup-flag", 
           dedup_flag, 
        "Flag determinating what information is used for dedupping. "
        "(0x1: read begin position. 0x2: read end position. 0x4: read QNAME. 0x8: UMI molecular-barcode). "
        "If this flag is zero, then infer such information from sequencing data. "
        "For example, if set to 0x3, then reads having the same begin and end positions are considered to be duplicates of each other. ");
    
    // less frequently used options
    ADD_OPTDEF(app, 
        "--bed-outfname", 
           bed_out_fname, 
        "The BED file to which genomic-region information will be written. "
        "Empty string means no output is written. This BED file can be generated by tumor and used by normal. ");
    ADD_OPTDEF(app,
        "--bed-infname",
           bed_in_fname,
        "The BED file from which genomic-region information is read. Empty string means no input is read. "
        "This BED file can be generated by tumor and used by normal. This param overrides the -R param. ");
    
    ADD_OPTDEF(app,
        "--fixedthresBQ", 
           fixedthresBQ, 
        "Base quality cutoff. This parameter is only for generating statistics and does not affect variant quality. ");
    ADD_OPTDEF(app,
        "--vc-stats-fname", 
           vc_stats_fname, 
        "Output TSV file containing variant-call statistics. The default is the standard error stream. ");
    ADD_OPTDEF(app,
        "--uni-bias-thres", 
           uni_bias_thres, 
        "Unified-bias threshold for generating the filter strings in FORMAT/FTS. "
        "This parameter is only for generating statistics and therefore does not affect variant quality. "
        "Downstream hard filtering with FORMAT/FTS is possible. ");
ADD_OPTDEF(app,
        "--uni-bias-r-max", 
           uni_bias_r_max, 
        "Maximum unified-bias threshold used for reducing variant read support. This parameter does affect variant quality. ");
    ADD_OPTDEF(app,
        "--bias-flag-amp-snv", 
           bias_flag_amp_snv, 
        "Flag where each bit enables evidence reduction by each bias for SNVs and amplicon assay. "
        + stringvec_to_descstring(BIAS_TYPE_TO_MSG));
    ADD_OPTDEF(app,
        "--bias-flag-amp-indel", 
           bias_flag_amp_indel, 
        "Flag where each bit enables evidence reduction by each bias for InDels and amplicon assay. ");
    ADD_OPTDEF(app,
        "--bias-flag-cap-snv", 
           bias_flag_cap_snv, 
        "Flag where each bit enables evidence reduction by each bias for SNVs and capture assay. ");
    ADD_OPTDEF(app,
        "--bias-flag-cap-indel", 
           bias_flag_cap_indel, 
        "Flag where each bit enables evidence reduction by each bias for InDels and capture assay. ");
    
    ADD_OPTDEF(app,
        "--highqual-thres-snv", 
           highqual_thres_snv,
         "The SNV quality threshold above which the family quality is considered to be high");
    ADD_OPTDEF(app,
        "--highqual-thres-indel", 
           highqual_thres_indel,
        "The InDel quality threshold above which the family quality is considered to be high, "
        "zero means auto infer to highqual-thres-snv + 6 for Illumina/BGI and -4 for IonTorrent)");
    ADD_OPTDEF(app,
        "--highqual-min-ratio", 
           highqual_min_ratio,
        "The mininum ratio of the raw non-deduplicated read depth to the deduplicated read family depth "
        "to trigger tumor-normal comparison with high quality families only");
     
    //ADD_OPTDEF(app, "--phred-frag-indel-ext",        phred_max_frag_indel_ext,
    //        "Maximum phred score for the indel of one additional base (excluding the one base required for opening indel), capped at two additional bases");
    ADD_OPTDEF(app, 
        "--phred-frag-indel-basemax", 
           phred_max_frag_indel_basemax,
        "Maximum phred score for the opening of an indel (including the one base required for opening indel). ");
    ADD_OPTDEF(app, 
        "--phred-sscs-transition-CG-TA", 
           phred_max_sscs_transition_CG_TA, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for C:G > T:A transition");
    ADD_OPTDEF(app, 
        "--phred-sscs-transition-TA-CG", 
           phred_max_sscs_transition_TA_CG, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for T:A > C:G transition");
    ADD_OPTDEF(app, 
        "--phred-sscs-transversion-any", 
           phred_max_sscs_transversion_any, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for any transversion");
    ADD_OPTDEF(app, 
        "--phred-sscs-indel-open", 
           phred_max_sscs_indel_open, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for the opening of indel gap "
        "(the opening includes the insertion/deletion of one base)");
    ADD_OPTDEF(app, 
        "--phred-sscs-indel-ext", 
           phred_max_sscs_indel_ext, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for the extension of indel gap "
        "(excluding the extension of one base) per base");
    ADD_OPTDEF(app, 
        "--phred-dscs-all", 
           phred_max_dscs_all, 
        "Maximum phred score for double-strand consensus sequences (DSCSs) for all types of mutations");
    ADD_OPTDEF(app, 
        "--phred-pow-sscs-origin", 
           phred_pow_sscs_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for SNVs");
    ADD_OPTDEF(app, 
        "--phred-pow-sscs-indel-origin", 
           phred_pow_sscs_indel_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for InDels");
    ADD_OPTDEF(app, 
        "--phred-pow-dscs-origin", 
           phred_pow_dscs_origin, 
        "The phred-score that is subtracted from phred-dscs to get the power-law quality adjustment");
    
    ADD_OPTDEF(app, 
        "--phred-snv-to-indel-ratio", 
           phred_snv_to_indel_ratio,
        "The phred-score for the ratio of SNVs to InDels. ");
    ADD_OPTDEF(app, 
        "--phred-umi-dimret-mult-snv", 
           phred_umi_dimret_mult_snv,
        "The diminishing-return factor for the PHRED-scaled raw variant quality for amplicon SNVs");
    ADD_OPTDEF(app, 
        "--phred-umi-dimret-mult-indel", 
           phred_umi_dimret_mult_indel,
        "The diminishing-return factor for the PHRED-scaled raw variant quality for amplicon InDels");
    
    ADD_OPTDEF(app, 
        "--ess-georatio-dedup-cap", 
           ess_georatio_dedup_cap,
        "Geometric common ratio of the increase in the observed number of deduped reads as a function of the effective number of deduped reads "
        "(effective sample size) for capture assay. ");
    ADD_OPTDEF(app, 
        "--ess-georatio-dedup-amp", 
           ess_georatio_dedup_amp,
        "Geometric common ratio of the increase in the observed number of deduped reads as a function of the effective number of deduped reads "
        "(effective sample size) for amplicon assay. ");
    ADD_OPTDEF(app, 
        "--ess-georatio-duped-pcr", 
           ess_georatio_duped_pcr, 
        "Geometric common ratio of the increase in the observed number of duped reads as a function of the effective number of duped reads "
        "(effective sample size) for reads derived from one template molecule. ");
 
    ADD_OPTDEF(app, 
        "--minABQ-pcr-snv", 
           minABQ_pcr_snv, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for SNVs and PCR amplicon, "
        "recommend 25 for Illumina and 0 for IonTorrent. "); 
    ADD_OPTDEF(app, 
        "--minABQ-pcr-indel", 
           minABQ_pcr_indel, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for InDels and PCR amplicon, "
        "recommend 18 for Illumina and 0 for IonTorrent. ");
    ADD_OPTDEF(app, 
        "--minABQ-cap-snv", 
           minABQ_cap_snv, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for SNVs and hybrid selection, "
        "recommend 20 for Illumina and 0 for IonTorrent. ");
    ADD_OPTDEF(app, 
        "--minABQ-cap-indel", 
           minABQ_cap_indel, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for InDels and hybrid selection, "
        "recommend 13 for Illumina and 0 for IonTorrent. ");
    
    ADD_OPTDEF(app, 
        "--minMQ1", 
           minMQ1, 
        "Minimum root-mean-square (RMS) mapping quality (MQ) of non-dedupped raw reads below which the variant quality is capped to this RMS MQ. ");
    ADD_OPTDEF(app, 
        "--maxMQ",
           maxMQ,
        "Maximum mapping quality (MQ) of the aligned reads, highly recommended to be the 60 from BWA. ");
    ADD_OPTDEF(app, 
        "--min-edge-dist",
           min_edge_dist,
        "Minimum average number of bases to the left and right aligned positions "
        "below which variant quality is capped to 4 times this average number. ");
    
    ADD_OPTDEF(app, 
        "--central-readlen", 
           central_readlen, 
        "Central (median) value for read lengths, 0 means estimate from the data. ");
    ADD_OPTDEF(app, 
        "--bq-phred-added-misma", 
           bq_phred_added_misma,
        std::string("Additional base-quality phred score added to match and mismatch, recommend 4 for IonTorrent from Life Technologies. "
        "This parameter value is automatically inferred if the sequencing platform is not ") + SEQUENCING_PLATFORM_TO_DESC[SEQUENCING_PLATFORM_OTHER]);
    ADD_OPTDEF(app, 
        "--bq-phred-added-indel", 
           bq_phred_added_indel, 
        std::string("Additional base-quality phred score added to indel and no-indel, recommend 17 for Illumina and BGI. "
        "This parameter value is automatically inferred if the sequencing platform is not ") + SEQUENCING_PLATFORM_TO_DESC[SEQUENCING_PLATFORM_OTHER]);
    
    ADD_OPTDEF(app, 
        "--phred-hetero-gt-snp", 
           phred_germline_polymorphism,
        "Phred-scaled prior probability of germline polymorphism at a genomic position for SNPs. ");
    ADD_OPTDEF(app, 
        "--phred-hetero-gt-indel", 
           phred_germline_indel,
           "hred-scaled prior probability of germline polymorphism at a genomic position for InDels."),
    ADD_OPTDEF(app, 
        "--phred-triallelic-snp",
           phred_triallelic_indel,     
        "Phred-scaled prior probability that two or more forms of germline SNPs occur at a genomic position. ");
    ADD_OPTDEF(app, 
        "--phred-triallelic-indel", 
           phred_triallelic_indel,     
        "Phred-scaled prior probability that two or more forms of germline InDels occur at a genomic position. ");
    ADD_OPTDEF(app, 
        "--any-mul-contam-frac", 
           any_mul_contam_frac,
        "Multiplicative contamination rate for the fraction of reads generated by any source of contamination in any sample");
    ADD_OPTDEF(app, 
        "--t2n-mul-contam-frac", 
        t2n_mul_contam_frac,         
        "Multiplicative contamination rate for the fraction of tumor reads in the normal"); 
    ADD_OPTDEF(app, 
        "--t2n-add-contam-frac", 
           t2n_add_contam_frac,         
        "Additive contamination rate for the fraction of tumor reads in the normal");
    ADD_OPTDEF(app, 
        "--t2n-add-contam-transfrac", 
           t2n_add_contam_transfrac, 
        "Additive contamination rate for the fraction of tumor reads in the normal by read transfer");
    
    ADD_OPTDEF(app, 
        "--phred-frac-indel-error-before-barcode-labeling", 
           phred_frac_indel_error_before_barcode_labeling, 
        "PHRED-scaled fraction of InDel errors that occurred before the attachment of UMI single-strand barcodes. ");
    ADD_OPTDEF(app, 
        "--baq-per-aligned-base", 
           baq_per_aligned_base,
        "PHRED-scaled base alignment quality (BAQ) for each additional base between the variant and the read edge. "); 
   
    ADD_OPTDEF(app, 
        "--is-somatic-snv-filtered-by-any-nonref-germline-snv", 
           is_somatic_snv_filtered_by_any_nonref_germline_snv,
        "Set to 0 (zero, false) if reject any nonref germline and to 1 (one, true) if only reject the specific ALT germline for SNV candidate. "); 
    ADD_OPTDEF(app, 
        "--is-somatic-indel-filtered-by-any-nonref-germline-indel", 
           is_somatic_indel_filtered_by_any_nonref_germline_indel,
        "Set to 0 (zero, false) if reject any nonref germline and to 1 (one, true) if only reject the specific ALT germline for InDel candidate. "); 
    ADD_OPTDEF(app,
        "--amp-BQ-sqr-coef",
           amp_BQ_sqr_coef,
        "The square of Illumina root-mean-square base quality multiplied by this number is the max possible TLOD of variant quality for PCR-amplicon assays. "); 
    ADD_OPTDEF(app,
        "--cap-BQ-sqr-coef",
           cap_BQ_sqr_coef,
        "The square of Illumina root-mean-square base quality multiplied by this number is the max possible TLOD of variant quality for hybrid-capture assays. "); 
    
    ADD_OPTDEF(app, 
        "--phred-varcall-err-per-map-err-per-base", 
           phred_varcall_err_per_map_err_per_base,
        "The root-mean square mapping quality added by this number is the maximum possible TLOD part of variant quality. "); 

    ADD_OPTDEF(app, 
        "--powlaw-exponent", 
           powlaw_exponent,
        "Exponent of the power-law relationship between error probability and allele-fraction deviation, "
        "strongly recommended to use the default value. ");
    ADD_OPTDEF(app, 
        "--powlaw-anyvar-base", 
           powlaw_anyvar_base,
        "Error probability at allele-fraction of 1 or 100%, storngly recommended to use the default value");
    ADD_OPTDEF(app, 
        "--syserr-maxqual", 
           syserr_maxqual,
            "PHRED-scaled probability that a candidate of systematic error is actually non-systematic. "); 
    ADD_OPTDEF(app, 
        "--syserr-norm-devqual", 
           syserr_norm_devqual,
        "PHRED-scaled likelihood that the observed allele fraction additively deviates from the expected allele fraction "
        "by a multiplicative factor of two. "); 
    
    app.add_flag(
        "--Should-add-note", 
           should_add_note, 
        "Flag indicating if the program generates more detail (can be used for debugging) in the VCF result file. ");
    
    ADD_OPTDEF(app, 
        "--disable-dup-read-merge", 
           disable_dup_read_merge, 
        "Boolean (0: false, 1: true) indicating if the program should disable the merge of duplicate reads. ");
    ADD_OPTDEF(app, 
        "--enable-dup-read-vqual", 
           enable_dup_read_vqual, 
        "Boolean (0: false, 1: true) indicating if the program should enable the use of raw non-dedupped reads in the calculation of variant quality. ");
    
    unsigned int assay_type_uint = (unsigned int)assay_type;
    unsigned int molecule_tag_uint = (unsigned int)molecule_tag;
    unsigned int sequencing_platform_uint = (unsigned int)sequencing_platform;
    unsigned int pair_end_merge_uint = (unsigned int)pair_end_merge;
    ADD_OPTDEF(app, 
        "--assay-type", 
           assay_type_uint, 
        "Assay type. " + stringvec_to_descstring(ASSAY_TYPE_TO_MSG));
    ADD_OPTDEF(app, 
        "--molecule-tag", 
           molecule_tag_uint, 
        "Molecule tag. " + stringvec_to_descstring(MOLECULE_TAG_TO_MSG));
    ADD_OPTDEF(app, 
        "--sequencing-platform", 
           sequencing_platform_uint, 
        "Sequencing platform. " + stringvec_to_descstring(SEQUENCING_PLATFORM_TO_MSG));
    ADD_OPTDEF(app, 
        "--pair-end-merge", 
           pair_end_merge_uint, 
        "Mode for the merge of R1 and R2 in a read pair. " + stringvec_to_descstring(PAIR_END_MERGE_TO_MSG));   
    app.add_flag(
        "--Disable-duplex", 
           disable_duplex, 
        "Flag indicating if the program should disable the merge of "
        "two SSCSs (single-strand-consensus-sequences) into one DSCS (double-strand-consensus-sequence). "
        "The UMI of the duplex tag should be in the form of <alpha>+<beta>. ");
    app.add_flag(
        "--always-log", 
           always_log, 
        "Flag indicating if the program should generate detailed log results to stderr");
    
    ADD_OPTDEF(app, 
        "--dedup-center-mult", 
           dedup_center_mult, 
        "Exponential decay per additional base of distance. ");
    ADD_OPTDEF(app, 
        "--dedup-amplicon-count-to-surrcount-ratio", 
           dedup_amplicon_count_to_surrcount_ratio, 
        "centroidCount/surroundingCount of reads ending at a position above which the assay is inferred to be amplicon. "
        "Assay type can be override on command-line");
    ADD_OPTDEF(app, 
        "--dedup-amplicon-end2end-ratio", 
           dedup_amplicon_end2end_ratio,
        "oneEndSegmentCount/otherEndSegmentCount above which the assay is inferred to be amplicon UMI");
    
    ADD_OPTDEF(app, 
        "--bitflag-InDel-penal-t-UMI-n-UMI", 
           bitflag_InDel_penal_t_UMI_n_UMI, 
        "Advanced flag for comparing tumor with normal"); 
    ADD_OPTDEF(app, 
        "--ref-bias-awareness", 
           ref_bias_awareness, 
        "Boolean bit-vector indicating if germline calls are aware of reference bias. "
        "The 0x1 bit is for amplicon and the 0x2 bit is for non-amplicon. ");
    
    ADD_OPTDEF(app, 
        "--haplo-in-diplo-allele-perc", 
           haplo_in_diplo_allele_perc, 
        "Percent allele of a haplotype observed in a diploid sample");
    ADD_OPTDEF(app, 
        "--diplo-oneside-posbias-perc", 
           diplo_oneside_posbias_perc, 
        "Percent ratio of sequencing-segments not biased to either left or right side "
        "below which position bias is added to the germline filter for homozygous variants");
    ADD_OPTDEF(app, 
        "--diplo-twoside-posbias-perc", 
           diplo_twoside_posbias_perc, 
        "Percent ratio of sequencing-segments not biased to both left and right sides "
        "below which position bias is added to the germline filter for homozygous variants");
    ADD_OPTDEF(app, 
        "--haplo-oneside-posbias-perc", 
           haplo_oneside_posbias_perc, 
        "Percent ratio of sequencing-segments not biased to either left or right side "
        "below which position bias is added to the germline filter for heterozygous variants");

    ADD_OPTDEF(app, 
        "--haplo-twoside-posbias-perc", 
           haplo_twoside_posbias_perc, 
        "Percent ratio of sequencing-segments not biased to both left and right sides "
        "below which position bias is added to the germline filter for heterozygous variants");
    ADD_OPTDEF(app, 
        "--regside-nbases", 
           regside_nbases, 
        "A variant is in the side region (left, right, or both) if and only if "
        "the number of bases to the sequencing-segment end is at most this number. ");
    
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
        if (t2n_add_contam_transfrac < (double)FLT_MIN) { t2n_add_contam_transfrac = (double)FLT_MIN; }
    });
    
    CLI11_PARSE(app, argc, argv);
    return 0;
}

#endif

