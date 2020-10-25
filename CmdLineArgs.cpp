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
        samFile *sam_infile = sam_open(this->bam_input_fname.c_str(), "r");
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
        bq_phred_added_misma += 4;
        syserr_minABQ_pcr_snv += 0;
        syserr_minABQ_pcr_indel += 0;
        syserr_minABQ_cap_snv += 0;
        syserr_minABQ_cap_indel += 0;
    }
    if (SEQUENCING_PLATFORM_ILLUMINA == inferred_sequencing_platform && SEQUENCING_PLATFORM_OTHER != this->sequencing_platform) {
        bq_phred_added_indel += 0;
        bq_phred_added_misma += 0;
        syserr_minABQ_pcr_snv += 19 * 10;
        syserr_minABQ_pcr_indel += syserr_minABQ_pcr_snv - 9 * 10;
        syserr_minABQ_cap_snv += 19 * 10;
        syserr_minABQ_cap_indel += syserr_minABQ_cap_snv - 9 * 10;
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
    
    unsigned int assay_type_uint = (unsigned int)assay_type;
    unsigned int molecule_tag_uint = (unsigned int)molecule_tag;
    unsigned int sequencing_platform_uint = (unsigned int)sequencing_platform;
    unsigned int pair_end_merge_uint = (unsigned int)pair_end_merge;
    
    app.add_flag_function("-v,--version", version_cb, "Show the version of this program. ");
    
// *** 00. frequently used parameters
    ADD_OPTDEF(app, 
        "inputBAM", 
           bam_input_fname, 
        "Input coordinate-sorted BAM file that is supposed to contain raw reads. ")->required()->check(CLI::ExistingFile);
    ADD_OPTDEF(app, 
        "-f,--fasta", 
           fasta_ref_fname, 
        "Input reference fasta file, where the special value NA means not available. ");
    ADD_OPTDEF(app, 
        "-o,--output", 
           vcf_out_pass_fname, 
        "Output bgzipped VCF file. ");
    
    ADD_OPTDEF(app, 
        "-R,--regions-file", 
           bed_region_fname, 
        "Input BED region file which is optional and delimits the genomic regions to call variants from. "
        "If not provided then call variants from all the sequenced genomic regions. "
        "This overrides the -t parameter. Please ensure that each region here is small enough to prevent memory exhaustion. "
        "Typically, each region covers one exon. ")->check(CLI::ExistingFile);
    ADD_OPTDEF(app, 
        "--targets", 
           tier1_target_region, 
        "Optional input region string (for example, chr1:2-3 or chr1) which globally delimits the genomic region to call variants from. "
        "If not provided then call variants from all the sequenced genomic regions. "
        "This parameter does not affect memory consumption. ");
    ADD_OPTDEF(app, 
        "-s,--sample", 
           sample_name, 
        "Sample name which is optional. ");
    
    ADD_OPTDEF(app, 
        "-t,--threads", 
           max_cpu_num, 
        "Number of cpu cores or equivalently threads to use. ");
    
    ADD_OPTDEF(app, 
        "--outvar-flag", 
           outvar_flag, 
        std::string("Output-variant flag in bits specifying which type of variants are in the VCF output. ") +
        "The " + std::to_string(OUTVAR_GERMLINE) + " bit indicates germline variant. " +
        "The " + std::to_string(OUTVAR_SOMATIC) + " bit indicates somatic variant. " +
        "The " + std::to_string(OUTVAR_ANY) + " bit indicates variant of any origin. " +
        "The " + std::to_string(OUTVAR_NONREF) + " bit indicates non-ref symbolic region indicating the absence of variants. ");
    app.add_flag(
        "-A,--all-out", 
           should_output_all, 
        "All possible alleles including REF allele at each position is in <--output>. "
        "The output is base-pair resolution VCF. ");
    app.add_flag(
        "--all-germline-out", 
           should_output_all_germline,
        "All possible alleles including REF allele at each position is in <--output> for germline variants. ");
    ADD_OPTDEF(app, 
        "-q,--vqual", 
           vqual, 
        "Every variant that satisfies this minimum variant quality is in the <--output> file. ");
    
    ADD_OPTDEF(app,
        "--assay-type", 
           assay_type_uint, 
        "Assay type. " + stringvec_to_descstring(ASSAY_TYPE_TO_MSG));
    
    ADD_OPTDEF(app,
        "--fam-thres-highBQ-snv", 
           fam_thres_highBQ_snv,
        "Threshold of base quality below which the base support is discarded in a barcode family for SNVs (PMC6477992, in code). ");
    ADD_OPTDEF(app,
        "--fam-thres-highBQ-indel", 
           fam_thres_highBQ_indel,
        "Threshold of base quality below which the base support is discarded in a barcode family for InDels for the Illumina/BGI platforms. ");
    ADD_OPTDEF(app,
        "--fam-thres-dup1add", 
           fam_thres_dup1add,
        "Tier-1 threshold of barcode-family size (10.1073/pnas.1105422108, supermutant). ");
    ADD_OPTDEF(app,
        "--fam-thres-dup1perc", 
           fam_thres_dup1perc,
        "Tier-1 threshold of barcode-family percent identity derived from allele consensus. ");
    ADD_OPTDEF(app,
        "--fam-thres-dup2add", 
           fam_thres_dup2add,
        "Tier-2 threshold of barcode-family size (PMC4271547). ");
    ADD_OPTDEF(app,
        "--fam-thres-dup2perc", 
           fam_thres_dup2perc,
        "Tier-2 threshold of barcode-family percent identity derived from allele consensus (PMC4271547). ");
    
// *** 01. parameters of the names of files, samples, regions, etc
    
    ADD_OPTDEF(app, 
        "--tumor-vcf",
           vcf_tumor_fname,
        "Block-gzipped VCF file of the tumor sample as input to the BAM file of normal sample. "
        "If specified/unspecified, then the input BAM is considered to be from normal/tumor. ");
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
   
// *** 02. parameters that control input, output, and logs (driven by computational requirements and resources)
    
    ADD_OPTDEF(app, 
        "--somaticGT", 
           somaticGT, 
        "Boolean (0: false, 1: true) indicating if the GT in the output VCF file refers to the genotype of the tumor cells "
            "in a heterogeneous mixture of tumor cells and normal cells. "
        "Please set to zero to call germline variants in a homogeneous mixture of cells (general format specification is mentioned in PMC3137218). ");
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
        "-d,--min-depth", 
           min_depth_thres, 
        "Minimum depth of all alleles below which results are fitlered out and therefore not in the output VCF. ");
    ADD_OPTDEF(app, 
        "-a,--min-altdp", 
           min_altdp_thres, 
        "Minimum depth of ALT below which results are filtered out. ");
    ADD_OPTDEF(app, 
        "--v-min-alldp", 
           vdp,
        "Every variant at each locus with at least this total depth is always in the -o file. ");
    ADD_OPTDEF(app, 
        "--v-min-altdp", 
           vad,
        "Every variant that satisfies this minimum allele depth is always in the -o file. ");
    
    app.add_flag(
        "--Should-add-note", 
           should_add_note, 
        "Flag indicating if the program generates more detail (can be used for debugging) in the VCF result file. ");
    app.add_flag(
        "--always-log", 
           always_log, 
        "Flag indicating if the program should generate detailed log results to stderr. ");
    
// *** 03. parameters that are driven by the properties of the assay
    
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
    ADD_OPTDEF(app,
        "--primerlen",
           primerlen,
        "Usual maximum number of bases from each end of each insert that are part of the primers for amplicon assay "
        "(22 from https://doi.org/10.1007/978-1-4020-6241-4_5 and 24 from https://doi.org/10.1101/gr.3.3.s30). ");
    
    ADD_OPTDEF(app, 
        "--central-readlen", 
           central_readlen, 
        "Central (median) value for read lengths, 0 means estimate from the data. ");
    ADD_OPTDEF(app, 
        "--bq-phred-added-misma", 
           bq_phred_added_misma,
        std::string("Additional base-quality phred score added to match and mismatch, recommend 4 for IonTorrent from Life Technologies. "
        "This parameter value is automatically inferred if the sequencing platform is not provided. ") + SEQUENCING_PLATFORM_TO_DESC[SEQUENCING_PLATFORM_OTHER]);
    ADD_OPTDEF(app,
        "--bq-phred-added-indel", 
           bq_phred_added_indel, 
        std::string("Additional base-quality phred score added to InDel and no-InDel, recommend 17 for Illumina and BGI. "
        "This parameter value is automatically inferred if the sequencing platform is not provided. ") + SEQUENCING_PLATFORM_TO_DESC[SEQUENCING_PLATFORM_OTHER]);
    
    ADD_OPTDEF(app, 
        "--powlaw-exponent", 
           powlaw_exponent,
        "Exponent of the power-law relationship between error probability and allele-fraction deviation. "
        "It is strongly recommended to use the default value. "); // novel
    ADD_OPTDEF(app, 
        "--powlaw-anyvar-base", 
           powlaw_anyvar_base,
        "Error probability at allele-fraction of 1 or 100%. "
        "It is strongly recommended to use the default value. "); // novel
    
    ADD_OPTDEF(app, 
        "--penal4lowdep",
           penal4lowdep,
        "Maximum penalty in variant quality for low read support. ");
        
// *** 04. parameters for dedupping reads
    
    ADD_OPTDEF(app, 
        "--dedup-center-mult", 
           dedup_center_mult, 
        "Exponential decay per additional base of distance. ");
    ADD_OPTDEF(app, 
        "--dedup-amplicon-count-to-surrcount-ratio", 
           dedup_amplicon_count_to_surrcount_ratio, 
        "centroidCount/surroundingCount of reads ending at a position above which the assay is inferred to be amplicon. "
        "Assay type can be override on command-line. ");
    ADD_OPTDEF(app, 
        "--dedup-amplicon-count-to-surrcount-ratio-twosided", 
           dedup_amplicon_count_to_surrcount_ratio_twosided, 
        "centroidCount/surroundingCount of reads ending at a position above which the assay is inferred to be amplicon. "
        "Assay type can be override on command-line. ");
    
    ADD_OPTDEF(app, 
        "--dedup-amplicon-end2end-ratio", 
           dedup_amplicon_end2end_ratio,
        "oneEndSegmentCount/otherEndSegmentCount above which the assay is inferred to be amplicon UMI. ");
    ADD_OPTDEF(app, 
        "--dedup-flag", 
           dedup_flag, 
        "Flag determinating what information is used for dedupping. "
        "(0x1: read begin position. 0x2: read end position. 0x4: read QNAME. 0x8: UMI molecular-barcode). "
        "If this flag is zero, then infer such information from sequencing data. "
        "For example, if set to 0x3, then reads having the same begin and end positions are considered to be duplicates of each other. ");
    
// *** 05. parameters related to bias thresholds
    
    ADD_OPTDEF(app,
        "--bias-thres-highBQ", 
           bias_thres_highBQ, 
        "Threshold of base quality (BQ) above which the base is considered to be of good BQ. ");
    ADD_OPTDEF(app,
        "--bias-thres-highBAQ", 
           bias_thres_highBAQ, 
        "Threshold of base alignment quality (BAQ) above which the base is considered to be of good BAQ. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-aLPxT-add",
           bias_thres_aLPxT_add,
        "The threshold for the number of bases to the left/right sequenced-segment end below which the segment is considered to be affected by position bias. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLPxT-perc",
           bias_thres_aLPxT_perc,
        "The threshold --bias_thres_aLPxT_add is increased by thisparam/100 times the average length of nearby InDels. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-PFXM1T-add",
           bias_thres_PFXM1T_add,
        "The tier1-threshold of 10x mismatch (XM) below which the estimated 100x number of passing-filter (PF) reads decreases according to the inverse-square law. ");
    ADD_OPTDEF(app,
        "--bias-thres-PFXM2T-add",
           bias_thres_PFXM1T_add,
        "The tier2-threshold of 10x mismatch (XM) below which the estimated 100x number of passing-filter (PF) reads decreases quadratically. ");
    ADD_OPTDEF(app,
        "--bias-thres-PFGO1T-add",
           bias_thres_PFGO1T_add,
        "The tier1-threshold of 10x gap-opening (GO) below which the estimated 100x number of PF reads decreases according to the inverse-square law. ");
    ADD_OPTDEF(app,
        "--bias-thres-PFGO2T-add",
           bias_thres_PFGO2T_add,
        "The tier2-threshold of 10x gap-opening (GO) below which the estimated 100x number of PF reads decreases according to the inverse-square law. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-PFXM1T-perc",
           bias_thres_PFXM1T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    ADD_OPTDEF(app,
        "--bias-thres-PFXM2T-perc",
           bias_thres_PFXM1T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    ADD_OPTDEF(app,
        "--bias-thres-PFGO1T-perc",
           bias_thres_PFGO1T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    ADD_OPTDEF(app,
        "--bias-thres-PFGO2T-perc",
           bias_thres_PFGO2T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-aLRP1t-minus",
           bias_thres_aLRP1t_minus,
        "The tier-1-threshold of the number of bases to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRP2t-minus",
           bias_thres_aLRP2t_minus,
        "The tier-2-threshold of the number of bases to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRB1t-minus",
           bias_thres_aLRB1t_minus,
        "The tier-1-threshold of BAQ (base alignment quality) to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRB2t-minus",
           bias_thres_aLRB2t_minus,
        "The tier-1-threshold of BAQ (base alignment quality) to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-aLRP1t-avgmul-perc",
           bias_thres_aLRP1t_avgmul_perc,
        "The percent of the tier-1 average number of bases to the left/right segment ends (LRP) below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRP2t-avgmul-perc",
           bias_thres_aLRP2t_avgmul_perc,
        "The percent of the tier-2 average number of bases to the left/right segment ends (LRP) below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRB1t-avgmul-perc",
           bias_thres_aLRB1t_avgmul_perc,
        "The percent of the tier-1 average BAQ to the left/right segment ends (LRB) below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRB2t-avgmul-perc",
           bias_thres_aLRB2t_avgmul_perc,
        "The percent of the tier-2 average BAQ to the left/right segment ends (LRB) below which the read support is not effective. ");

    ADD_OPTDEF(app,
        "--bias-thres-aLRI1T-perc",
           bias_thres_aLRI1T_perc,
        "Tier-1 percent coefficient of the average number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRI2T-perc",
           bias_thres_aLRI2T_perc,
        "Tier-2 percent coefficient of the average number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRI1t-perc",
           bias_thres_aLRI1t_perc,
        "Tier-1 percent coefficient of the average BAQ (base alignment quality) to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRI2t-perc",
           bias_thres_aLRI2t_perc,
        "Tier-2 percent coefficient of the average BAQ (base alignment quality) to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-aLRI1T-add",
           bias_thres_aLRI1T_add,
        "Tier-1 additive constant to the number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-aLRI2T-add",
           bias_thres_aLRI2T_add,
        "Tier-2 additive constant to the number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-PFBQ1",
           bias_thres_PFBQ1,
        "The tier-1 threshold of base quality (BQ) below which "
            "the estimated 100x number of PF (passing filter) reads decreases according to the inverse-square law. ");
    ADD_OPTDEF(app,
        "--bias-thres-PFBQ2",
           bias_thres_PFBQ2,
        "The tier-2 threshold of base quality (BQ) below which "
            "the estimated 100x number of PF (passing filter) reads decreases according to the inverse-square law. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-aXM1T-add",
           bias_thres_aXM1T_add,
        "The tier-1-threshold of 10x mismatch (XM) below which the read support is not effective. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-interfering-indel",
           bias_thres_interfering_indel,
        "If the reference allele is at most this number of bases away from the nearest low-quality InDel, then this reference allele is interfered by this InDel. ");
    ADD_OPTDEF(app,
        "--bias-thres-interfering-indel-BQ",
           bias_thres_interfering_indel_BQ,
        "The nearby InDel is of low-quality if the minimum BQ (base quality) of the two bases anchoring the InDel (plus all inserted bases if applicable) "
            "is lower than this value. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-BAQ1",
           bias_thres_BAQ1,
        "The tier-1-threshold of base alignment quality (BAQ) below which the read support is not effective. ");
    ADD_OPTDEF(app,
        "--bias-thres-BAQ2",
           bias_thres_BAQ2,
        "The tier-2-threshold of base alignment quality (BAQ) below which the read support is not effective. ");
    
    ADD_OPTDEF(app,
        "--bias-thres-FTS-FA",
           bias_thres_FTS_FA,
        "If the bias-reduced allele fraction multiplied by this parameter is less than the non bias-reduced allele fraction, then flag the variant for bias in FTS. ");

// *** 06. parameters related to the priors of bias
    
    ADD_OPTDEF(app,
        "--bias-prior-DPadd-perc",
           bias_prior_DPadd_perc,
        "The percentage of variant InDel read support that is not considered when checking if any InDel is nearby. ");
    
    ADD_OPTDEF(app,
        "--bias-priorfreq-pos",
           bias_priorfreq_pos,
        "The prior weight of the null hypothesis that there is no segment position bias given no other information. ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-indel-in-read-div",
           bias_priorfreq_indel_in_read_div,
        "Reduction in the prior weight of the null hypothesis that there is no segment position bias, given that the reads support both the ALT allele and some InDel allele(s). ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-indel-in-var-div2",
           bias_priorfreq_indel_in_var_div2,
        "Additional reduction in the prior weight of the null hypothesis of no segment position bias, given that some InDel allele and this ALT allle overlap. ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-indel-in-STR-div2",
           bias_priorfreq_indel_in_var_div2,
        "Additional reduction in the prior weight of the null hypothesis of no segment position bias, given that the variant is in an InDel region. ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-var-in-STR-div2",
           bias_priorfreq_var_in_STR_div2,
        "Additional reduction in the prior weight of the null hypothesis of no segment position bias, given that this variant is in a STR region. ");
    
    ADD_OPTDEF(app,
        "--bias-prior-var-DP-mul",
           bias_prior_var_DP_mul,
        "Minimum read depth odds ratio of nearby InDel allele to the variant allele above which the additional reduction in the prior weight is in effect. ");
    
    ADD_OPTDEF(app,
        "--bias-priorfreq-ipos-snv",
           bias_priorfreq_ipos_snv,
        "Prior weight of the null hypothesis of no insert position bias for SNVs. ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-ipos-indel",
           bias_priorfreq_ipos_indel,
        "Prior weight of the null hypothesis of no insert position bias for InDels. ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-strand-snv-base",
           bias_priorfreq_strand_snv_base,
        "Prior weight of the null hypothesis of no strand bias for SNVs at the average base quality of zero. ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-strand-indel",
           bias_priorfreq_strand_indel,
        "Prior weight of the null hypothesis of no strand bias for InDels. ");
    
    ADD_OPTDEF(app,
        "--bias-FA-pseudocount-indel-in-read",
          bias_FA_pseudocount_indel_in_read,
        "Pseudocount of the effective variant allele depth for segment position bias in read regions affected by InDels. ");
    
    ADD_OPTDEF(app,
        "--bias-priorfreq-orientation-snv-base",
           bias_priorfreq_orientation_snv_base,
        "The prior weight of the null hypothesis of no read-orientation bias at 100\% allele fraction for SNVs. ");
    ADD_OPTDEF(app,
        "--bias-priorfreq-orientation-indel-base",
           bias_priorfreq_orientation_indel_base,
        "The prior weight of the null hypothesis of no read-orientation bias at 100\% allele fraction for InDels. ");
    ADD_OPTDEF(app,
        "--bias-removal-orientation-avg-end-len",
           bias_orientation_counter_avg_end_len,
        "If the number of bases to segment end for the ALT allele is at least this param and the one for all alleles is at most this param, "
            "then read-orientation bias is removed. ");
    
    ADD_OPTDEF(app,
        "--bias-FA-powerlaw-noUMI-phred-inc-snv",
           bias_FA_powerlaw_noUMI_phred_inc_snv,
        "The Phred-scale decrease in false positive rates applied on top of universality for non-UMI SNVs. ");
    ADD_OPTDEF(app,
        "--bias-FA-powerlaw-noUMI-phred-inc-indel",
           bias_FA_powerlaw_noUMI_phred_inc_indel,
        "The Phred-scale decrease in false positive rates applied on top of universality for non-UMI InDels. ");
    ADD_OPTDEF(app,
        "--bias-FA-powerlaw-withUMI-phred-inc-snv",
           bias_FA_powerlaw_withUMI_phred_inc_snv,
        "The Phred-scale decrease in false positive rates applied on top of universality for UMI SNVs. ");
    ADD_OPTDEF(app,
        "--bias-FA-powerlaw-withUMI-phred-inc-indel",
           bias_FA_powerlaw_withUMI_phred_inc_indel,
        "The Phred-scale decrease in false positive rates applied on top of universality for UMI Indels. ");

    ADD_OPTDEF(app,
        "--nobias-flag",
           nobias_flag,
        "Advanced flag for reducing one bias by another bias. ");
    
    ADD_OPTDEF(app,
        "--nobias-pos-indel-lenfrac-thres",
           nobias_pos_indel_lenfrac_thres,
        "Threshold of (F*L) above which segment position bias is nullified, where F is allele fraction and L is InDel length (aka number of gap extensions). ");
    ADD_OPTDEF(app,
        "--nobias-pos-indel-STR-track-len",
           nobias_pos_indel_STR_track_len,
        "Threshold of short tandem repeat (STR) track length above which segment position bias is nullified if the variant is also the majority in the STR track. ");
    

// *** 07. parameters related to read families
    
    ADD_OPTDEF(app,
        "--fam-thres-emperr-all-flat-snv",
           fam_thres_emperr_all_flat_snv, 
        "Mininum number of read support of all alleles needed to compute the empirical basecall error for SNVs. ");
    ADD_OPTDEF(app,
        "--fam-thres-emperr-all-flat-indel",
           fam_thres_emperr_all_flat_indel, 
        "Mininum number of read support of all alleles needed to compute the empirical basecall-like error for InDels. ");
    ADD_OPTDEF(app,
        "--fam-thres-emperr-con-perc-snv",
           fam_thres_emperr_con_perc_snv,
        "Mininum percent of read support of consensus allele needed to compute the empirical basecall error for SNVs. ");
    ADD_OPTDEF(app,
        "--fam-thres-emperr-con-perc-indel",
           fam_thres_emperr_con_perc_indel, 
        "Mininum percent of read support of consensus alllele needed to compute the empirical basecall-like error for InDels. ");
    
    ADD_OPTDEF(app,
        "--fam-min-n-copies",
           fam_min_n_copies,
        "Mininum number of DNA copies (1 ng contains approximately 300 copies of human genome) to have zero penalty for UMI-labeled barcode-family (PMC6197739 and PMC5856404). ");
    ADD_OPTDEF(app,
        "--fam-min-overseq-perc",
           fam_min_overseq_perc,
        "Mininum percentage of over-sequencing (one plus average family size) to have zero penalty for UMI-labeled barcode-family. ");
    

    ADD_OPTDEF(app,
        "--fam-phred-indel-inc-before-barcode-labeling",
           fam_phred_indel_inc_before_barcode_labeling,
        "PHRED-scaled fraction of InDel errors that occurred before the attachment of UMI single-strand barcodes. ");
    ADD_OPTDEF(app, 
        "--fam-phred-sscs-transition-CG-TA", 
           fam_phred_sscs_transition_CG_TA, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for C > T transition which is mainly caused by cytosine deanimation (PMC3437896). ");
    ADD_OPTDEF(app, 
        "--fam-phred-sscs-transition-AT-GC", 
           fam_phred_sscs_transition_AT_GC, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for A:T > G:C transition (PMC3437896). ");
    ADD_OPTDEF(app, 
        "--fam-phred-sscs-transversion-CG-AT",
           fam_phred_sscs_transversion_CG_AT, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for the G > T transversion which is mainly caused by 8-oxo-guanine (PMC3437896). ");

    ADD_OPTDEF(app, 
        "--fam-phred-sscs-transversion-other", 
           fam_phred_sscs_transversion_other, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for any other transversion (PMC3437896). ");
    ADD_OPTDEF(app, 
        "--fam-phred-sscs-indel-open", 
           fam_phred_sscs_indel_open, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for the opening of InDel gap "
        "(the opening includes the insertion/deletion of one base). ");
    ADD_OPTDEF(app, 
        "--fam-phred-sscs-indel-ext", 
           fam_phred_sscs_indel_ext, 
        "Maximum phred score for single-strand consensus sequences (SSCSs) for the extension of InDel gap "
            "(excluding the extension of one base) per base. ");
    ADD_OPTDEF(app,
        "--fam-phred-dscs-all", 
           fam_phred_dscs_all, 
        "Maximum phred score for double-strand consensus sequences (DSCSs) for all types of mutations (PMC3437896). ");
    
    ADD_OPTDEF(app, 
        "--fam-phred-pow-sscs-transversion-AT-TA-origin",
           fam_phred_pow_sscs_transversion_AT_TA_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for the A:T > T:A type of SNV (PMC3973132). "
        "It should be lower if the mutations are known to be caused by aristolochic acid (PMC3973132). ");
    ADD_OPTDEF(app, 
        "--fam-phred-pow-sscs-snv-origin",
           fam_phred_pow_sscs_snv_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for all other types of SNVs (PMC3437896). ");
    ADD_OPTDEF(app, 
        "--fam-phred-pow-sscs-indel-origin", 
           fam_phred_pow_sscs_indel_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for InDels. ");
    ADD_OPTDEF(app, 
        "--fam-phred-pow-dscs-all-origin", 
           fam_phred_pow_dscs_all_origin,
        "The phred-score that is subtracted from phred-dscs to get the power-law quality adjustment (PMC3437896). ");
    
// *** 08. parameters related to systematic errors
    
    ADD_OPTDEF(app,
        "--syserr-BQ-prior", 
           syserr_BQ_prior, 
        "Phred-scale probability of the systematic error (defined at PMC3295828). "); 
    
    ADD_OPTDEF(app,
        "--syserr-BQ-sbratio-q-add",
           syserr_BQ_sbratio_q_add,
        "Squared relative deviation for the odds ratio of strand read support by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF(app,
        "--syserr-BQ-sbratio-q-max",
           syserr_BQ_sbratio_q_max,
        "The maximum increased quality for the odds ratio of strand read support by which the threshold of low-quality bases is increased. ");
    ADD_OPTDEF(app,
        "--syserr-BQ-xmratio-q-add",
           syserr_BQ_xmratio_q_add,
        "Squared relative deviation for the odds ratio of mismatch by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF(app,
        "--syserr-BQ-xmratio-q-max",
           syserr_BQ_xmratio_q_max,
        "The maximum increased quality for the odds ratio of mismatch by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF(app,
        "--syserr-BQ-bmratio-q-add",
           syserr_BQ_bmratio_q_add,
        "Squared relative deviation for the odds ratio of same-base mismatch by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF(app,
        "--syserr-BQ-bmratio-q-max",
           syserr_BQ_bmratio_q_max,
        "The maximum increased quality for the odds ratio of same-base mismatch by which the threshold of low-quality bases is increased. "); 
    
    ADD_OPTDEF(app,
        "--syserr-BQ-strand-favor-mul",
           syserr_BQ_strand_favor_mul,
        "The strand with higher average base quality (BQ) is given this much increase in weight to overcome the potential low average BQ on the other strand. "); 

    ADD_OPTDEF(app,
        "--syserr-minABQ-pcr-snv", 
           syserr_minABQ_pcr_snv, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for SNVs and PCR amplicon. "
        "The recommend value is 25 for Illumina and 0 for IonTorrent. "); 
    ADD_OPTDEF(app, 
        "--syserr-minABQ-pcr-indel", 
           syserr_minABQ_pcr_indel, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for InDels and PCR amplicon. "
        "The recommend value is 18 for Illumina and 0 for IonTorrent. ");
    ADD_OPTDEF(app, 
        "--syserr-minABQ-cap-snv", 
           syserr_minABQ_cap_snv, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for SNVs and hybrid selection. "
        "The recommended is 20 for Illumina and 0 for IonTorrent. ");
    ADD_OPTDEF(app,
        "--syserr-minABQ-cap-indel", 
           syserr_minABQ_cap_indel, 
        "Average base quality (ABQ) below which variant quality is capped to ABQ for InDels and hybrid selection. "
        "The ecommended value is 13 for Illumina and 0 for IonTorrent. ");
    
    ADD_OPTDEF(app,
        "--syserr-maxMQ",
           syserr_maxMQ,
        "Maximum mapping quality (MQ) of the aligned reads, highly recommended to be the 60 from BWA. ");
    ADD_OPTDEF(app, 
        "--syserr-phred-varcall-err-per-map-err-per-base", 
           syserr_phred_varcall_err_per_map_err_per_base,
        "The root-mean square mapping quality added by this number is the maximum possible TLOD part of variant quality. "); 

// *** 09. parameters related to germline vars 
    
    ADD_OPTDEF(app, 
        "--germ-hetero-FA", 
           germ_hetero_FA,
        "ALT allele fraction of heterozygous germline polymorphism (PMC6587442). ");
    
    ADD_OPTDEF(app,
        "--germ-phred-hetero-snp", 
           germ_phred_hetero_snp,
        "Phred-scaled prior probability of 0/1 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
    ADD_OPTDEF(app, 
        "--germ-phred-homalt-snp", 
           germ_phred_homalt_snp,
        "Phred-scaled prior probability of 1/1 homozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
    ADD_OPTDEF(app, 
        "--germ-phred-het3al-snp",
           germ_phred_het3al_snp,     
        "Phred-scaled prior probability of 1/2 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
    
    ADD_OPTDEF(app, 
        "--germ-phred-hetero-indel", 
           germ_phred_hetero_indel,
        "Phred-scaled prior probability of 0/1 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). "),
    ADD_OPTDEF(app, 
        "--germ-phred-homalt-indel", 
           germ_phred_homalt_indel,
        "Phred-scaled prior probability of 1/1 homozygous germline polymorphism at a genomic position for InDels (common empirical knowledge). "),
    ADD_OPTDEF(app, 
        "--germ-phred-het3al-indel", 
           germ_phred_het3al_indel,   
        "Phred-scaled prior probability of 1/2 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
        
// *** 10. parameters related to tumor-normal-pairs.
    ADD_OPTDEF(app,
        "--tn-q-inc-max",
           tn_q_inc_max,
        "Maximum Phred-scale increase in variant quality by comparing the tumor with its matched normal. Theoretically, it should be Phred-scaled 2 to the power of --powlaw-exponent. ");
    ADD_OPTDEF(app,
        "--tn-syserr-norm-devqual",
           tn_syserr_norm_devqual,
        "Phred-scale decrease in the variant quality subtracted from the tumor if the tumor FA (allele fraction) deviates by +100\%/-50\% from its matched normal FA. ");

// *** 11. parameters related to InDels.
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/
    ADD_OPTDEF(app,
        "--indel-BQ-max",
           indel_BQ_max,
        "For InDels, the maximum base quality of the InDel of one base (PMC4719071). ");
    ADD_OPTDEF(app,
        "--indel-str-repeatsize-max",
           indel_str_repeatsize_max,
        "For InDels, the maximum length of a repeating substring for the repeat to be considered as eligible for short tandem repeat (STR) (PMC5054066). ");
    ADD_OPTDEF(app,
        "--indel-polymerase-size",
           indel_polymerase_size,
        "For InDels, the number of bases that can be fit within the catalytic unit of the polymerase used for PCR (PMC149199). ");
    ADD_OPTDEF(app,
        "--indel-polymerase-slip-rate",
           indel_polymerase_slip_rate,
        "For InDels, the number of bases that can slip when processed by the catalytic unit of the polymerase used for PCR (PMC149199). ");
    ADD_OPTDEF(app,
        "--indel-del-to-ins-err-ratio",
           indel_del_to_ins_err_ratio,
        "For InDels, the number of times that deletion is more likely (PMC149199). ");
    ADD_OPTDEF(app,
        "--indel-adj-tracklen-dist",
           indel_adj_tracklen_dist,
        "InDels within this number of bases of each other are considered to be close to each other. ");
    ADD_OPTDEF(app,
        "--indel-adj-indellen-perc",
           indel_adj_indellen_perc,
        "For InDels, the multiplicative percent coefficient of the InDel length used to determine the base-distance threshold. "
        "Other InDels found within this threshold make the InDel at this position more likely to be false positive. ");
    
    ADD_OPTDEF(app,
        "--indel-multiallele-samepos-penal",
           indel_multiallele_samepos_penal,
        "For InDels, the Phred-scale penalty of having more than two alleles at the same VCF position. ");
    ADD_OPTDEF(app,
        "--indel-multiallele-diffpos-penal",
           indel_multiallele_diffpos_penal,
        "For InDels, the Phred-scale penalty of having more than two alelles that are overlapping with each other. ");
    ADD_OPTDEF(app,
        "--indel-multiallele-soma-penal-thres",
           indel_multiallele_soma_penal_thres,
        "For InDels, the Phred-scale penalty threshold above which penalty applies for somatic variants. ");
    ADD_OPTDEF(app,
        "--indel-tetraallele-germ-penal-value",
           indel_tetraallele_germline_penal_value,
        "For InDels, the Phres-scale penalty of having more than three alleles at the same VCF position for germline variants. ");
    ADD_OPTDEF(app,
        "--indel-tetraallele-germ-penal-thres",
           indel_tetraallele_germline_penal_thres,
        "For InDels, the Phred-scale penalty threshold above which penalty applies for germline variants. ");
    
    ADD_OPTDEF(app,
        "--indel-ins-penal-pseudocount",
           indel_ins_penal_pseudocount,
        "For InDels, the length of inserted sequence at which at most half of the inserted sequences are erroneous. ");
    
    ADD_OPTDEF(app,
        "--indel-nonSTR-phred-per-base",
           indel_nonSTR_phred_per_base,
        "The Phred quality per each additional base to the sequenced-segment left/right end in non-STR regions. ");
    ADD_OPTDEF(app,
        "--indel-STR-phred-per-region",
           indel_STR_phred_per_region,
        "The Phred quality per each additional STR track to the sequenced-segment left/right end. ");
        
// *** 12. parameters related to contamination
    
    ADD_OPTDEF(app, 
        "--contam-any-mul-frac", 
           contam_any_mul_frac,
        "Multiplicative contamination rate for the fraction of reads generated by any source of contamination in any sample (PMC3167057). ");
    ADD_OPTDEF(app, 
        "--contam-t2n-mul-frac", 
           contam_t2n_mul_frac,
        "Multiplicative contamination rate for the fraction of tumor reads in the normal (PMC6528031). "); 

// *** 13. parameters related to micro-adjustment (they do not have any clear theory support)
    
    ADD_OPTDEF(app, 
        "--microadjust_xm", 
           microadjust_xm,
        "Average number of mismatches over 150bp above which one-base insertion has quality capped at the base-quality of the inserted base. ");
    ADD_OPTDEF(app, 
        "--microadjust_delFAQmax", 
           microadjust_delFAQmax,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");

    ADD_OPTDEF(app, 
        "--microadjust-bias-pos-indel-fold", 
           microadjust_bias_pos_indel_fold,
        "If the number of reads supporting nearby InDels multiplied by this number is more than the number of reads supporting the variant, then the variant is considered to be affected by nearby InDels. ");
    ADD_OPTDEF(app, 
        "--microadjust-bias-pos-indel-misma-to-indel-ratio",
           microadjust_bias_pos_indel_misma_to_indel_ratio,
        "If the number of reads supporting nearby mismatches multiplied by this number is more than the number of reads supporting the InDel variant, then the InDel variant is considered to be affected by nearby mismatches. ");
    
    ADD_OPTDEF(app, 
        "--microadjust-nobias-pos-indel-misma-to-indel-ratio",
           microadjust_nobias_pos_indel_misma_to_indel_ratio,
        "If the number of reads supporting nearby mismatches multiplied by this number is less than the number of reads supporting the InDel variant, the STR track length is long enough, and the majority of nearby InDels are of this variant type, then this InDel variant is not affected by any position bias. ");
    ADD_OPTDEF(app, 
        "--microadjust-nobias-pos-indel-maxlen",
           microadjust_nobias_pos_indel_maxlen,
        "The maximum effective InDel length for <--nobias-pos-indel-lenfrac-thres>. ");
    ADD_OPTDEF(app, 
        "--microadjust-nobias-pos-indel-bMQ",
           microadjust_nobias_pos_indel_bMQ,
        "If the root-mean-square mapping quality (RMS-MQ) of the InDel variant is above this threshold and <--microadjust-nobias-pos-indel-perc> is passed, then this variant is not subject to insert bias. ");
    ADD_OPTDEF(app, 
        "--microadjust-nobias-pos-indel-perc",
           microadjust_nobias_pos_indel_perc,
        "The ratio of XM2 to (aDPff + aDPfr + aDPrf + aDPrr) in VCF tags. ");
    ADD_OPTDEF(app, 
        "--microadjust-nobias-strand-all-fold",
           microadjust_nobias_strand_all_fold,
        "The ratio of read count on the unbiased strand to the one on the biased strand, above which the strand bias and insert bias are removed. ");
    ADD_OPTDEF(app, 
        "--microadjust-refbias-indel-max",
           microadjust_refbias_indel_max,
        "Maximum reference bias for InDels. ");
    
    ADD_OPTDEF(app, 
        "--microadjust-fam-binom-qual-halving-thres",
           microadjust_fam_binom_qual_halving_thres,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");
    ADD_OPTDEF(app, 
        "--microadjust-fam-lowfreq-invFA",
           microadjust_fam_lowfreq_invFA,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");
    ADD_OPTDEF(app, 
        "--microadjust-ref-MQ-dec-max",
           microadjust_ref_MQ_dec_max,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");

/// *** end

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
        // if (t2n_add_contam_transfrac < (double)FLT_MIN) { t2n_add_contam_transfrac = (double)FLT_MIN; }
    });
    
    CLI11_PARSE(app, argc, argv);
    return 0;
}

#endif

