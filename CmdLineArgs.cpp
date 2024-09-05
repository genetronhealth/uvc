#ifndef CmdLineArgs_INCLUDED
#define CmdLineArgs_INCLUDED

#include "CmdLineArgs.hpp"
#include "version.h"

#include "CLI11-1.7.1/CLI11.hpp"
#include "htslib/sam.h"

#include <fstream>
#include <float.h>

#define BQ_PHRED_ADDED_MISMA_IONTORRENT 8
#define SYSERR_MINABQ_SNV_ILLUMINA 200
#define SYSERR_MINABQ_INDEL_ILLUMINA 100

#define REPLACE_UNDERSCORE_WITH_HYPHEN(x) (std::string("--") + replace_underscore_with_hyphen( #x ))
#define ADD_OPTDEF(app, k, v, msg) (app.add_option(k, v, msg, true))
#define ADD_OPTDEF2(app, k, msg) (app.add_option(REPLACE_UNDERSCORE_WITH_HYPHEN(k), k, msg, true))

#define UPDATE_NON_NEG_MINUS(a, b) { (a) = ((a) - MIN((a), (b))); }

const std::string 
replace_underscore_with_hyphen(const std::string astring) {
    std::string ret = "";
    for (auto letter : astring) {
        if (letter == '_') {
            ret.push_back('-');
        } else {
            ret.push_back(letter);
        }
    }
    return ret;
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
        uvc1_unsigned_int_t countPE = 0;
        uvc1_unsigned_int_t countSE = 0;
        std::vector<uvc1_unsigned_int_t> qlens;
        qlens.reserve(5000+1);
        qlens.push_back(150);
        uvc1_unsigned_int_t q20_n_fail_bases = 0;
        uvc1_unsigned_int_t q30_n_fail_bases = 0;
        uvc1_unsigned_int_t q30_n_pass_bases = 0;
        while (sam_read1(sam_infile, samheader, b) >= 0 && (countPE + countSE) < 5000) {
            this->inferred_maxMQ = MAX(this->inferred_maxMQ, b->core.qual);
            if (b->core.flag & 0x1) {
                countPE++;
            } else {
                countSE++;
            }
            qlens.push_back(b->core.l_qseq);
            /*
            int skip_l_qlen = 0;
            int skip_r_qlen = 0;
            if ((!(b->core.flag & 0x4)) && (b->core.n_cigar > 0)) {
                 const auto cigar = bam_get_cigar(b);
                 if (BAM_CSOFT_CLIP == bam_cigar_op(cigar[0])) {
                    skip_l_qlen = bam_cigar_oplen(cigar[0]);
                 }
                 if (BAM_CSOFT_CLIP == bam_cigar_op(cigar[b->core.n_cigar-1])) {
                    skip_r_qlen = bam_cigar_oplen(cigar[b->core.n_cigar-1]);
                 }
            }
            
            for (int qpos = skip_l_qlen ; qpos < ((int)b->core.l_qseq) - skip_r_qlen; qpos++) */
            for (int qpos = 0; qpos < ((int)b->core.l_qseq); qpos++) {
                uvc1_unsigned_int_t bq = (bam_get_qual((b))[(qpos)]);
                if (bq < 30) {
                    q30_n_fail_bases++;
                } else {
                    q30_n_pass_bases++;
                }
                if (bq < 20) {
                    q20_n_fail_bases++;
                }
            }
        }
        std::sort(qlens.begin(), qlens.end());
        if (0 == this->central_readlen) { this->central_readlen = qlens.at(qlens.size()/2); }
        bam_destroy1(b);
        bam_hdr_destroy(samheader);
        sam_close(sam_infile);
        const bool isPE = (0 < countPE);
        const bool is_2x_Q20toQ30_lessthan_Q30plus = (2 * (q30_n_fail_bases - q20_n_fail_bases) < q30_n_pass_bases);
        const bool is_4x_Q20toQ30_lessthan_Q30plus = (4 * (q30_n_fail_bases - q20_n_fail_bases) < q30_n_pass_bases);
        // const bool isQ30BQ = (q30_n_fail_bases * 3 < q30_n_pass_bases);
        // const bool isQ30BQsoft = (q30_n_fail_bases < q30_n_pass_bases);
        const bool isfixqlen = (qlens.at(qlens.size()/2) * 100 > qlens.at(qlens.size()-1) * 95);
        if (isPE || is_4x_Q20toQ30_lessthan_Q30plus || (is_2x_Q20toQ30_lessthan_Q30plus && isfixqlen)) {
            inferred_sequencing_platform = SEQUENCING_PLATFORM_ILLUMINA;
        } else {
            inferred_sequencing_platform = SEQUENCING_PLATFORM_IONTORRENT;
        }
        std::cerr << "Inferred_sequencing_platform=" <<  SEQUENCING_PLATFORM_TO_NAME.at(inferred_sequencing_platform)
                << "\n IsPairedEnd=" << isPE 
                << "\n is 2*BQ30<(BQ30-BQ20) passed = " << is_2x_Q20toQ30_lessthan_Q30plus
                << "\n is 4*BQ30<(BQ30-BQ20) passed = " << is_4x_Q20toQ30_lessthan_Q30plus
                << "\n isFixedReadQuerySeqLength = " << isfixqlen 
                << std::endl;

    }
    if (SEQUENCING_PLATFORM_IONTORRENT == inferred_sequencing_platform && SEQUENCING_PLATFORM_OTHER != this->sequencing_platform) {
        bq_phred_added_indel += 0;
        bq_phred_added_misma += BQ_PHRED_ADDED_MISMA_IONTORRENT;
        syserr_minABQ_pcr_snv += 0;
        syserr_minABQ_pcr_indel += 0;
        syserr_minABQ_cap_snv += 0;
        syserr_minABQ_cap_indel += 0;
        
        UPDATE_NON_NEG_MINUS(fam_thres_highBQ_snv, 30);
        UPDATE_NON_NEG_MINUS(fam_thres_highBQ_indel, 30);
        UPDATE_NON_NEG_MINUS(bias_thres_PFBQ1, 30);
        UPDATE_NON_NEG_MINUS(bias_thres_PFBQ2, 30);
        UPDATE_NON_NEG_MINUS(bias_thres_highBQ, 13);
    }
    if (SEQUENCING_PLATFORM_ILLUMINA == inferred_sequencing_platform && SEQUENCING_PLATFORM_OTHER != this->sequencing_platform) {
        bq_phred_added_indel += 0;
        bq_phred_added_misma += 0;
        syserr_minABQ_pcr_snv += SYSERR_MINABQ_SNV_ILLUMINA; // it was 19 * 10
        syserr_minABQ_pcr_indel += SYSERR_MINABQ_INDEL_ILLUMINA;  // modified accordingly
        syserr_minABQ_cap_snv += SYSERR_MINABQ_SNV_ILLUMINA; 
        syserr_minABQ_cap_indel += SYSERR_MINABQ_INDEL_ILLUMINA; 
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

template <class TContainer>
const std::string 
stringvec_to_phrase(const TContainer & v) {
    if (v.size() == 0) {
        return "";
    }
    if (v.size() == 1) {
        return v[0];
    }
    std::string ret;
    for (uvc1_unsigned_int_t i = 0; i < v.size(); i++) {
        if ((i + 1) == v.size()) { ret += " and "; }
        else if (i > 0 && i < (i + 1)) { ret += ", "; }
        ret += v[i];
    }
    return ret;
};

const std::string 
stringvec_to_descstring(const std::vector<std::string> & v) {
    std::string ret;
    for (uvc1_unsigned_int_t i = 0; i < v.size(); i++) {
        ret += std::to_string(i) + " : " + v[i] + ". ";
    }
    return ret;
};

int
CommandLineArgs::initFromArgCV(int & parsing_result_flag, int argc, const char *const* argv) {
    parsing_result_flag = -1;
    auto version_cb = [](int n_values IGNORE_UNUSED_PARAM){
        std::cout << "uvc-" << VERSION << std::endl;
        exit(0);
    };
    CLI::App app{(std::string("UVC version ") + VERSION_DETAIL)};
    
    uvc1_unsigned_int_t assay_type_uint = (uvc1_unsigned_int_t)assay_type;
    uvc1_unsigned_int_t molecule_tag_uint = (uvc1_unsigned_int_t)molecule_tag;
    uvc1_unsigned_int_t sequencing_platform_uint = (uvc1_unsigned_int_t)sequencing_platform;
    uvc1_unsigned_int_t pair_end_merge_uint = (uvc1_unsigned_int_t)pair_end_merge;
    
    app.add_flag_function("-v,--version", version_cb, "Show the version of this program. ");
    
// *** 00. frequently used parameters
    ADD_OPTDEF(app, 
        "inputBAM", 
           bam_input_fname, 
        ("The input coordinate-sorted and indexed BAM file that is supposed to contain raw reads. "
        "If set to " OPT_ONLY_PRINT_VCF_HEADER ", then only print the VCF header, that describes the output format and is not instantiated from the input files, and then exit with the exit code of zero. "
        "Important warnings about potential mis-use and mis-understanding are mentioned with the keyword CAVEAT in the VCF header. "))->required();
    ADD_OPTDEF(app, 
        "-f,--fasta", 
           fasta_ref_fname, 
        "The input reference FASTA file to which the inputBAM is aligned, where the special value NA means not available. ");
    ADD_OPTDEF(app, 
        "-o,--output", 
           vcf_out_pass_fname, 
        "The output bgzipped VCF file. "
        "If this parameter is set to either the empty string (\"\") or dot (\".\") and the --fam-consensus-out-fastq parameter is not set to either the empty string (\"\") or dot (\".\"), "
        "then do not generate any VCF. ");
    
    ADD_OPTDEF(app, 
        "-R,--regions-file", 
           bed_region_fname, 
        "Optional input BED region file which delimits the genomic regions to call variants from. "
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
    ADD_OPTDEF2(app, mem_per_thread,
        "Approximate amount of RAM, in terms of mega-bytes, used per thread. Please note that this is not the exact amount consumed at runtime. ");
    
    ADD_OPTDEF(app, 
        "--outvar-flag", 
           outvar_flag, 
        std::string("Output-variant flag in bits specifying which type of variants are in the VCF output. ")
        + "The " + std::to_string(OUTVAR_GERMLINE) + " bit indicates germline variant. "
        + "The " + std::to_string(OUTVAR_SOMATIC) + " bit indicates somatic variant. "
        + "The " + std::to_string(OUTVAR_ANY) + " bit indicates variant of any origin. "
        + "The " + std::to_string(OUTVAR_MGVCF) + " bit indicates symbolic MGVCF (multi-sample genomic VCF). "
                "Each line of MGVCF contains multiple regions of potentially very different depths. "
                "MGVCF is similar to GVCF but allows for easy comparison of sequencing depths of multiple samples at any arbitrary position). "
        + "The " + std::to_string(OUTVAR_ADDITIONAL_INDEL_CANDIDATE) + " bit indicates position with noisy alignments nearby, "
                "which can be a candidate for long InDel, CNV, SV, etc. "
        + "The " + std::to_string(OUTVAR_BASE_NN) + " bit indicates padded deletion at a nucleotide-base position between two adjacent gap positions. "
        + "The " + std::to_string(OUTVAR_LINK_NN) + " bit indicates padded deletion at a gap position between two adjacent nucleotide-base positions. "
        );
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
        "--tumor-vcf",
           vcf_tumor_fname,
        "Block-gzipped VCF file of the tumor sample as input to the BAM file of normal sample. "
        "If specified/unspecified, then the input BAM is considered to be from normal/tumor. ");
    
    ADD_OPTDEF(app, 
        "--fam-consensus-out-fastq",
           fam_consensus_out_fastq,
        "Filename prefix of the output FASTQ file generated by the tier-1 consensus of fragments in the same UMI families. "
        "The filename suffixes are (" + stringvec_to_phrase(FASTQ_LIKE_SUFFIXES) + ") "
        "for (PE-R1-FASTQ, PE-R2-FASTQ, and SE-FASTQ) files, respectively, "
        "where PE and SE mean paire-end and single-end, respectively. "
        "If this parameter is set to either the empty string (\"\") or dot (\".\"), then consensus by UMI-group is not performed so that no FASTQ file is generated. ");
    ADD_OPTDEF2(app, fam_consensus_out_fastq_thres_dup1add,
        "Tier-1 threshold of barcode-family size (PMC3111315, supermutant) for the generation of consensus FASTQ file. "); 

// *** 01. parameters of the names of files, samples, regions, etc
    
    ADD_OPTDEF2(app, bed_out_fname,
        "The BED file to which genomic-region information will be written. "
        "Empty string (\"\") and dot (\".\") mean no output is written. This BED file can be generated by tumor and used by normal. ");
    ADD_OPTDEF2(app, bed_in_fname,
        "The BED file from which genomic-region information is read. Empty string (\"\") and dot (\".\") mean no input is read. "
        "This BED file can be generated by tumor and used by normal. This param overrides the <--regions-file> parameter. ");
    ADD_OPTDEF2(app, bed_in_avg_sequencing_DP,
        "Average sequencing depth in the BED input file specified above. If set to -1, then infer from the input BAM file. ");
    ADD_OPTDEF2(app, bed_in_avg_sequencing_DP_n_from_t,
        "If set to 1, then the number of reads covered by each BED region is taken from the --bed-out-fname generated from the corresponding tumor BAM. ");

// *** 02. parameters that control input, output, and logs (driven by computational requirements and resources)
    
    ADD_OPTDEF2(app, is_tumor_format_retrieved, 
        "Boolean (0: false, 1: true) indicating if the format from the tumor VCF should be retrieved in the tumor-normal comparison. "
        "This boolean has no effect if <--tumor-vcf> is not provided. ");
    
    ADD_OPTDEF2(app, kept_aln_min_aln_len,
        "Minimum alignment length below which the alignment is filtered out. ");
    ADD_OPTDEF2(app, kept_aln_min_mapqual,
        "Minimum mapping quality below which the alignment is filtered out. ");
    ADD_OPTDEF2(app, kept_aln_min_isize,
        "Minimum insert size below which the alignment is filtered out (not applicable to zero insert size). ");
    ADD_OPTDEF2(app, kept_aln_max_isize,
        "Maximum insert size above which the alignment is filtered out (not applicable to zero insert size). "
        "It is recommended to set this param to 167*2-30 for ctDNA assay (https://link.springer.com/article/10.1007/s40291-021-00534-6). ");
    ADD_OPTDEF2(app, kept_aln_is_zero_isize_discarded,
        "The number 0 (false) or 1 (true) indicating if an alignment with zero insert size should be filtered out. ");
    
    ADD_OPTDEF2(app, min_altdp_thres, 
        "Minimum allele depth of fragments below which allele record is not in the <--output> VCF. "
        "The paramters <--all-out> and <--all-germline-out> take precedence over this parameter. ");
    
    ADD_OPTDEF2(app, vdp1,
        "Every variant allele with at least this total depth of highBQ segments is always in the <--output> VCF if the <--vad1> and <--vfa1> conditions are also satisfied. "
        "This is used to rescue potential variant candidates in regions with high sequencing depth. ");
    ADD_OPTDEF2(app, vad1,
        "Every variant allele with at least this allele depth of highBQ segments is always in the <--output> VCF if the <--vdp1> and <--vfa1> conditions are also satisfied. ");
    ADD_OPTDEF2(app, vfa1,
        "Every variant allele with at least this allele fraction of highBQ segments is always in the <--output> VCF if the <--vdp1> and <--vad1> conditions are also satisfied. ");
    ADD_OPTDEF2(app, vdp2,
        "Every variant allele with at least this total depth of fragments is always in the <--output> VCF if the <--vad2> and <--vfa2> conditions are also satisfied. "
        "This is used to rescue potential variant candidates in regions with high fragment depth. ");
    ADD_OPTDEF2(app, vad2,
        "Every variant allele with at least this allele depth of fragments is always in the <--output> VCF if the <--vdp2> and <--vfa2> conditions are also satisfied. ");
    ADD_OPTDEF2(app, vfa2,
        "Every variant allele with at least this allele fraction of fragments is always in the <--output> VCF if the <--vdp2> and <--vad2> conditions are also satisfied. ");
    
    ADD_OPTDEF2(app, min_r_ad,
        "Every reference allele with less than this allele depth of fragments is always not in the <--output> VCF. "
        "This paramter takes precedence over the parameters <--all-out> and <--all-germline-out>. ");
    ADD_OPTDEF2(app, min_a_ad,
        "Every variant allele with less than this allele depth of fragments is always not in the <--output> VCF. "
        "This paramter takes precedence over the parameters <--all-out> and <--all-germline-out>. ");
    
    ADD_OPTDEF2(app, should_add_note, 
        "Flag indicating if the program generates more detail (can be used for debugging) in the VCF result file. ");
    ADD_OPTDEF2(app, always_log, 
        "Flag indicating if the program should generate detailed log results to stderr. ");
    
// *** 03. parameters that are driven by the properties of the assay
    
    ADD_OPTDEF(app, 
        "--molecule-tag", 
           molecule_tag_uint, 
        std::string("Molecule tag. ") + stringvec_to_descstring(MOLECULE_TAG_TO_MSG)
        + "The UMI of a single-strand family shoud be in the form of #label "
        + "(for example, the UMI tag for the read with name SRR123.456#ATCC is ATCC). "
        + "The UMI of a duplex family should be in the form of #<alpha>+<beta> "
        + "(for example, the alpha and beta UMI tags for the read with name SRR123.456#ATCC+TGAA are ATCC and TGAA, respectively). ");
    
    ADD_OPTDEF(app, 
        "--sequencing-platform", 
           sequencing_platform_uint, 
        "Sequencing platform. " + stringvec_to_descstring(SEQUENCING_PLATFORM_TO_MSG));
    ADD_OPTDEF(app, 
        "--pair-end-merge", 
           pair_end_merge_uint, 
        "Mode for merging R1 and R2 in a read pair. " + stringvec_to_descstring(PAIR_END_MERGE_TO_MSG));   
    ADD_OPTDEF2(app, disable_duplex, 
        "Flag indicating if the program should not merge "
        "two SSCSs (single-strand-consensus-sequences) into one DSCS (double-strand-consensus-sequence). ");
    ADD_OPTDEF2(app, primerlen,
        "Usual maximum number of bases from each end of each insert that are part of the primers for amplicon assay "
        "(22 from https://doi.org/10.1007/978-1-4020-6241-4_5 and 24 from https://doi.org/10.1101/gr.3.3.s30). ");
    ADD_OPTDEF2(app, primerlen2,
        "Primer length that is used to cancel out bias. ");
    ADD_OPTDEF2(app, primer_flag,
        "Bitwise flag for PCR primers. "
        "The 0x1 bit means: use the matched normal sample to filter out false positive calls on primers to get more true positive calls. "
        "The 0x2 bit means: even if the value of the --primerlen parameter is greater than zero, still do not skip variants in the primer read regions, but instead check if the read is part of an amplicon to determine the presence/absence of primers (By default: if --primerlen is set to a strictly positive value, then always trim read support in read bases that are inferred to be of primer origin accordingly). "
        "The 0x4 bit means: even if the value of the --primerlen parameter is greater than zero, still do not assume that the variant is generated by PCR-amplicons, but instead check it, and then perform amplicon-specific variant-quality adjustment only if the variant is generated by PCR amplicons. ");
    
    ADD_OPTDEF2(app, central_readlen, 
        "Central (median) value for read lengths, 0 means estimate from the data. ");
    ADD_OPTDEF2(app, bq_phred_added_misma,
        "Additional base-quality Phred score added to match and mismatch. "
        "The recommended value is " STRINGIZE(BQ_PHRED_ADDED_MISMA_IONTORRENT) " for " PLAT_ION_LIKE " and zero for " PLAT_ILLUMINA_LIKE ". "
        "This parameter value is automatically inferred if the sequencing platform is not provided. ");
    ADD_OPTDEF2(app, bq_phred_added_indel, 
        std::string("Additional base-quality Phred score added to InDel and no-InDel. "));
    
    ADD_OPTDEF2(app, powlaw_exponent,
        "Exponent of the power-law relationship between error probability and allele-fraction deviation. "
        "It is strongly recommended to use the default value. "); // novel
    ADD_OPTDEF2(app, powlaw_anyvar_base,
        "Phred-scaled error probability at the maximum possible allele-fraction of 1 or 100%. "
        "It is strongly recommended to use the default value. "); // novel
    ADD_OPTDEF2(app, powlaw_amplicon_allele_fraction_coef,
        "Multiplicative coefficient for the effective allele fraction in PCR-amplicon sequencing compared with hybrid-capture sequencing, WES, or WGS. "
        "It is strongly recommended to use the default value. ");

    ADD_OPTDEF2(app, penal4lowdep,
        "Maximum penalty in variant quality for low read support. ");
    ADD_OPTDEF2(app, assay_sequencing_BQ_max,
        "Maximum basecall quality (BQ) determined by empirical error rate. ");
    ADD_OPTDEF2(app, assay_sequencing_BQ_inc,
        "Each base-call quality (BQ) in the input BAM file is added by this parameter. If this parameter is set to a negative integer, then reduce BQ. ");
    
    ADD_OPTDEF2(app, phasing_haplotype_max_count,
        "The maximum number of haplotypes covering each genomic position used to generate the bHap and cHap format fields in the output VCF. ");
    ADD_OPTDEF2(app, phasing_haplotype_min_ad,
        "This value plus the number of variants in a haplotype is the minimum haplotype allele depth for the haplotype (bHap and cHap) to be present in the output VCF. "
        "For example, if a haplotype links variants A and B, then this haplotype requires the linkage between A and B to have a depth of at least (<this-value> + 2), "
        "and this depth becomes (<this-value> + 3) if this haplotype links three variants instead of two, etc. ");
    ADD_OPTDEF2(app, phasing_haplotype_max_detail_cnt,
        "The maximum number of haplotypes that have detailed read counts. ");
   
// *** 04. parameters for dedupping reads
    
    ADD_OPTDEF2(app, dedup_center_mult, 
        "Exponential decay per additional base of distance. ");
    /*
    ADD_OPTDEF2(app, dedup_amplicon_count_to_surrcount_ratio, 
        "centroidCount/surroundingCount of reads ending at a position above which the assay is inferred to be amplicon. "
        "Assay type can be override on command-line. ");
    ADD_OPTDEF2(app, dedup_amplicon_count_to_surrcount_ratio_twosided, 
        "centroidCount/surroundingCount of reads ending at a position above which the assay is inferred to be amplicon. "
        "Assay type can be override on command-line. ");
    */
    ADD_OPTDEF2(app, dedup_amplicon_border_to_insert_cov_weak_avgDP_ratio,
        "centroidInsertBorderSequencingDepth/meanInsertSequencingDepth above which the assay is inferred to be amplicon (require both left and right borders). "
        "Assay type can be override on command-line. This condition is necessary but not sufficient. ");
    ADD_OPTDEF2(app, dedup_amplicon_border_to_insert_cov_strong_avgDP_ratio,
        "centroidInsertBorderSequencingDepth/meanInsertSequencingDepth above which the assay is inferred to be amplicon (require either left or right border). "
        "Assay type can be override on command-line. This condition is necessary but not sufficient. ");
    ADD_OPTDEF2(app, dedup_amplicon_border_to_insert_cov_weak_totDP_ratio,
        "centroidInsertBorderSequencingDepth/totalInsertSequencingDepth above which the assay is inferred to be amplicon (require both left and right borders). "
        "Assay type can be override on command-line. This condition is necessary but not sufficient. ");
    ADD_OPTDEF2(app, dedup_amplicon_border_to_insert_cov_strong_totDP_ratio,
        "centroidInsertBorderSequencingDepth/totalInsertSequencingDepth above which the assay is inferred to be amplicon (require either left or right border). "
        "Assay type can be override on command-line. This condition is necessary but not sufficient. ");
    ADD_OPTDEF2(app, dedup_amplicon_border_weak_minDP,
        "centroidInsertBorderSequencingDepth above which the assay is inferred to be amplicon (require both left and right borders). "
        "Assay type can be override on command-line. This condition is necessary but not sufficient. ");
    ADD_OPTDEF2(app, dedup_amplicon_border_strong_minDP,
        "centroidInsertBorderSequencingDepth above which the assay is inferred to be amplicon (require either left or right border). "
        "Assay type can be override on command-line. This condition is necessary but not sufficient. ");
    
    ADD_OPTDEF2(app, dedup_amplicon_end2end_ratio,
        "oneEndSegmentCount/otherEndSegmentCount above which the amplicon assay is inferred to be using anchored-PCR (a.k.a. one-sided PCR). ");
    ADD_OPTDEF2(app, dedup_flag, 
        "Flag determinating what information is used for dedupping. "
        "(0x1: read begin position. 0x2: read end position. 0x4: read QNAME. 0x8: UMI molecular-barcode). "
        "If this flag is zero, then infer such information from sequencing data. "
        "For example, if set to 0x3, then reads having the same begin and end positions are considered to be duplicates of each other. "
        "If the 0x4 bit is set, then only reads having the same read name can be assumed to originate from the same original DNA molecule, which essentially disables duplicate detection. ");
    
// *** 05. parameters related to bias thresholds
    
    ADD_OPTDEF2(app, bias_thres_highBQ, 
        "Threshold of base quality (BQ) above which the base is considered to be of good BQ. "
        "Please note that this value is automatically decreased by 13 for IonTorrent/LifeTech/ThermoFisher sequencing data. ");
    ADD_OPTDEF2(app, bias_thres_highBAQ, 
        "Threshold of base alignment quality (BAQ) above which the base is considered to be of good BAQ. ");
    
    ADD_OPTDEF2(app, bias_thres_aLPxT_add,
        "The threshold for the number of bases to the left/right sequenced-segment end below which the segment is considered to be affected by position bias. ");
    ADD_OPTDEF2(app, bias_thres_aLPxT_perc,
        "The threshold --bias_thres_aLPxT_add is increased by thisparam/100 times the average length of nearby InDels. ");

#if COMPILATION_ENABLE_XMGOT
    ADD_OPTDEF2(app, bias_thres_PFXM1T_add,
        "The tier1-threshold of 10x mismatch (XM) below which the estimated 100x number of passing-filter (PF) reads decreases according to the inverse-square law. ");
    ADD_OPTDEF2(app, bias_thres_PFXM2T_add,
        "The tier2-threshold of 10x mismatch (XM) below which the estimated 100x number of passing-filter (PF) reads decreases quadratically. ");
    ADD_OPTDEF2(app, bias_thres_PFGO1T_add,
        "The tier1-threshold of 10x gap-opening (GO) below which the estimated 100x number of PF reads decreases according to the inverse-square law. ");
    ADD_OPTDEF2(app, bias_thres_PFGO2T_add,
        "The tier2-threshold of 10x gap-opening (GO) below which the estimated 100x number of PF reads decreases according to the inverse-square law. ");
    
    ADD_OPTDEF2(app, bias_thres_PFXM1T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    ADD_OPTDEF2(app, bias_thres_PFXM2T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    ADD_OPTDEF2(app, bias_thres_PFGO1T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    ADD_OPTDEF2(app, bias_thres_PFGO2T_perc,
        "The percent increase in the corresponding threshold relative to the background noise level. ");
    
    ADD_OPTDEF2(app, bias_thres_PFXM1NT_perc,
        "The percent increase in the corresponding threshold relative to the background noise level for the matched normal sample. ");
    ADD_OPTDEF2(app, bias_thres_PFGO1NT_perc,
        "The percent increase in the corresponding threshold relative to the background noise level for the matched normal sample. ");
#endif
    
    ADD_OPTDEF2(app, bias_thres_aLRP1t_minus,
        "The tier-1-threshold of the number of bases to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    ADD_OPTDEF2(app, bias_thres_aLRP2t_minus,
        "The tier-2-threshold of the number of bases to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    ADD_OPTDEF2(app, bias_thres_aLRB1t_minus,
        "The tier-1-threshold of BAQ (base alignment quality) to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    ADD_OPTDEF2(app, bias_thres_aLRB2t_minus,
        "The tier-1-threshold of BAQ (base alignment quality) to the left/right segment ends (LRP) below which the read support is not effective, "
            "substracted from the average. ");
    
    ADD_OPTDEF2(app, bias_thres_aLRP1t_avgmul_perc,
        "The percent of the tier-1 average number of bases to the left/right segment ends (LRP) below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_aLRP2t_avgmul_perc,
        "The percent of the tier-2 average number of bases to the left/right segment ends (LRP) below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_aLRB1t_avgmul_perc,
        "The percent of the tier-1 average BAQ to the left/right segment ends (LRB) below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_aLRB2t_avgmul_perc,
        "The percent of the tier-2 average BAQ to the left/right segment ends (LRB) below which the read support is not effective. ");
    
    ADD_OPTDEF2(app, bias_thres_aLRP1Nt_avgmul_perc,
        "The percent of the tier-1 average number of bases to the left/right segment ends (LRP) below which the read support is not effective "
        "for the matched normal sample. ");
    ADD_OPTDEF2(app, bias_thres_aLRB1Nt_avgmul_perc,
        "The percent of the tier-1 average BAQ to the left/right segment ends (LRP) below which the read support is not effective "
        "for the matched normal sample. ");
    
    ADD_OPTDEF2(app, bias_thres_aLRI1T_perc,
        "Tier-1 percent coefficient of the average number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_aLRI2T_perc,
        "Tier-2 percent coefficient of the average number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_aLRI1t_perc,
        "Tier-1 percent coefficient of the average BAQ (base alignment quality) to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_aLRI2t_perc,
        "Tier-2 percent coefficient of the average BAQ (base alignment quality) to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    
    ADD_OPTDEF2(app, bias_thres_aLRI1NT_perc,
        "Tier-1 percent coefficient of the average number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective for the matched normal. ");
    ADD_OPTDEF2(app, bias_thres_aLRI1Nt_perc,
        "Tier-1 percent coefficient of the average BAQ (base alignment quality) to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective for the matched normal. ");

    ADD_OPTDEF2(app, bias_thres_aLRI1T_add,
        "Tier-1 additive constant to the number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_aLRI2T_add,
        "Tier-2 additive constant to the number of bases to the left/right insert ends (LRI), "
            "used for determinating the threshold below which the read support is not effective. ");
    
    ADD_OPTDEF2(app, bias_thres_PFBQ1,
        "The tier-1 threshold of base quality (BQ) below which "
            "the estimated 100x number of PF (passing filter) reads decreases according to the inverse-square law for " PLAT_ILLUMINA_LIKE ". "
            "For " PLAT_ION_LIKE ", this parameter will be subtracted by 30 before being used. ");
    ADD_OPTDEF2(app, bias_thres_PFBQ2,
        "The tier-2 threshold of base quality (BQ) below which "
            "the estimated 100x number of PF (passing filter) reads decreases according to the inverse-square law for " PLAT_ILLUMINA_LIKE ". "
            "For " PLAT_ION_LIKE ", this parameter will be subtracted by 30 before being used. ");
    
    ADD_OPTDEF2(app, bias_thres_aXM1T_add,
        "The tier-1-threshold of 10x mismatch (XM) below which the read support is not effective. ");
    
    ADD_OPTDEF2(app, bias_thres_interfering_indel,
        "If the reference allele is at most this number of bases away from the nearest low-quality InDel, then this reference allele is interfered by this InDel. ");
    ADD_OPTDEF2(app, bias_thres_interfering_indel_BQ,
        "The nearby InDel is of low-quality if the minimum BQ (base quality) of the two bases anchoring the InDel (plus all inserted bases if applicable) "
            "is lower than this value. ");
    
    ADD_OPTDEF2(app, bias_thres_BAQ1,
        "The tier-1-threshold of base alignment quality (BAQ) below which the read support is not effective. ");
    ADD_OPTDEF2(app, bias_thres_BAQ2,
        "The tier-2-threshold of base alignment quality (BAQ) below which the read support is not effective. ");
    
    ADD_OPTDEF2(app, bias_thres_strict_c2LRP0,
        "The strictly defined position bias without any prior of non-bias, "
            "which uses SSCS support not within X-bps of the left/right SSCS mapping coordinates to eliminate biased SSCS support, "
            "where X is the value of this command-line parameter. "
        "This bias is always assumed to be present (although its strength can be weak to have any effect at all) "
            "and hence not subject to the Bayes-factor with the null hypothesis of having no bias. ");
    
    ADD_OPTDEF2(app, bias_thres_FTS_FA,
        "If the non-bias-reduced allele fraction multiplied by this parameter is greater than the bias-reduced allele fraction, then flag the variant for bias in FTS. ");
    ADD_OPTDEF2(app, bias_is_orientation_artifact_mixed_with_sequencing_error,
        "Model false positive calls as a mixture of read-orientation-specific error and sequencing error. "
            "Improve specificity at the cost of sensitivity. "
            "This may be useful if both heavy FFPE artifact and heavy sequencing error are present. ");
    ADD_OPTDEF2(app, bias_orientation_min_effective_allelefrac,
        "The minimum allele fraction that is effective for orientation bias. "
            "This parameter should be set to the allele fraction of oxidation artifact. ");


// *** 06. parameters related to the priors of bias
    
    ADD_OPTDEF2(app, bias_prior_DPadd_perc,
        "The percentage of variant InDel read support that is not considered when checking if any InDel is nearby. ");
    
    ADD_OPTDEF2(app, bias_priorfreq_pos,
        "The prior weight of the null hypothesis that there is no segment position bias given no other information. ");
    ADD_OPTDEF2(app, bias_priorfreq_indel_in_read_div,
        "Reduction in the prior weight of the null hypothesis that there is no segment position bias, "
        "given that the reads support both the ALT allele and some InDel allele(s). ");
    ADD_OPTDEF2(app, bias_priorfreq_indel_in_var_div2,
        "Additional reduction in the prior weight of the null hypothesis of no segment position bias, given that some InDel allele and this ALT allle overlap. ");
    ADD_OPTDEF2(app, bias_priorfreq_indel_in_str_div2,
        "Additional reduction in the prior weight of the null hypothesis of no segment position bias, given that the variant is in an InDel region. ");
    ADD_OPTDEF2(app, bias_priorfreq_var_in_str_div2,
        "Additional reduction in the prior weight of the null hypothesis of no segment position bias, given that this variant is in a STR region. ");
    
    ADD_OPTDEF2(app, bias_prior_var_DP_mul,
        "Minimum read depth odds ratio of nearby InDel allele to the variant allele above which the additional reduction in the prior weight is in effect. ");
    
    ADD_OPTDEF2(app, bias_priorfreq_ipos_snv,
        "Prior weight of the null hypothesis of no insert position bias for SNVs. ");
    ADD_OPTDEF2(app, bias_priorfreq_ipos_indel,
        "Prior weight of the null hypothesis of no insert position bias for InDels. ");
    ADD_OPTDEF2(app, bias_priorfreq_strand_snv_base,
        "Prior weight of the null hypothesis of no strand bias for SNVs at the average base quality of zero. ");
    ADD_OPTDEF2(app, bias_priorfreq_strand_indel,
        "Prior weight of the null hypothesis of no strand bias for InDels. ");
    
    ADD_OPTDEF2(app, bias_FA_pseudocount_indel_in_read,
        "Pseudocount of the effective variant allele depth for segment position bias in read regions affected by InDels. ");
    
    ADD_OPTDEF2(app, bias_priorfreq_orientation_snv_base,
        "The prior weight of the null hypothesis of no read-orientation bias at 100\% allele fraction for SNVs. ");
    ADD_OPTDEF2(app, bias_priorfreq_orientation_indel_base,
        "The prior weight of the null hypothesis of no read-orientation bias at 100\% allele fraction for InDels. ");
    ADD_OPTDEF2(app, bias_orientation_counter_avg_end_len,
        "If the number of bases to segment end for the ALT allele is at least this param and the one for all alleles is at most this param, "
            "then read-orientation bias is removed. ");
    
    ADD_OPTDEF2(app, bias_FA_powerlaw_noUMI_phred_inc_snv,
        "The Phred-scaled decrease in false positive rates applied on top of universality for non-UMI SNVs. ");
    ADD_OPTDEF2(app, bias_FA_powerlaw_noUMI_phred_inc_indel,
        "The Phred-scaled decrease in false positive rates applied on top of universality for non-UMI InDels. ");
    ADD_OPTDEF2(app, bias_FA_powerlaw_withUMI_phred_inc_snv,
        "The Phred-scaled decrease in false positive rates applied on top of universality for UMI SNVs. ");
    ADD_OPTDEF2(app, bias_FA_powerlaw_withUMI_phred_inc_indel,
        "The Phred-scaled decrease in false positive rates applied on top of universality for UMI Indels. ");
    ADD_OPTDEF2(app, bias_reduction_by_high_sequencingDP_min_n_totDepth,
        "Minimum all-allele total sequencing depth at each locus above which the number of effective reads supporting each bias is reduced in inverse proportion "
        "if the condition below is also satisfied. ");
    ADD_OPTDEF2(app, bias_reduction_by_high_sequencingDP_min_n_altDepth,
        "Minimum allele-specific  sequencing depth at each locus above which the number of effective reads supporting each bias is reduced in inverse proportion "
        "if the above condition is also satisfied. ");
    
    ADD_OPTDEF2(app, nobias_flag,
        "Advanced flag for reducing one bias by another bias. ");
    
    ADD_OPTDEF2(app, nobias_pos_indel_lenfrac_thres,
        "Threshold of (F*L) above which segment position bias is nullified, where F is allele fraction and L is InDel length (aka number of gap extensions). ");
    ADD_OPTDEF2(app, nobias_pos_indel_str_track_len,
        "Threshold of short tandem repeat (STR) track length above which segment position bias is nullified if the variant is also the majority in the STR track. ");
    
// *** 07. parameters related to read families
    
    ADD_OPTDEF2(app, fam_thres_highBQ_snv,
        "Threshold of base quality below which the base support is discarded in a barcode family for " PLAT_ILLUMINA_LIKE  " (PMC6477992, in code) for SNVs. "
        "For " PLAT_ION_LIKE ", this paramter will be subtracted by 30 before being used. ");
    ADD_OPTDEF2(app, fam_thres_highBQ_indel,
        "Threshold of base quality below which the base support is discarded in a barcode family for " PLAT_ILLUMINA_LIKE  " for InDels. "
        "For " PLAT_ION_LIKE ", this paramter will be subtracted by 30 before being used. ");
    ADD_OPTDEF2(app, fam_thres_dup1add,
        "Tier-1 threshold of barcode-family size (PMC3111315, supermutant). ");
    ADD_OPTDEF2(app, fam_thres_dup1perc,
        "Tier-1 threshold of barcode-family percent identity derived from allele consensus. ");
    ADD_OPTDEF2(app, fam_thres_dup2add,
        "Tier-2 threshold of barcode-family size (PMC4271547). ");
    ADD_OPTDEF2(app, fam_thres_dup2perc,
        "Tier-2 threshold of barcode-family percent identity derived from allele consensus (PMC4271547). ");
    ADD_OPTDEF2(app, fam_thres_qseqlen,
        "The average length of sequenced segments in an SSCS below which the SSCS support is discarded. ");
    
    ADD_OPTDEF2(app, fam_thres_emperr_all_flat_snv, 
        "Mininum number of read support of all alleles needed to compute the empirical basecall error for SNVs. ");
    ADD_OPTDEF2(app, fam_thres_emperr_all_flat_indel, 
        "Mininum number of read support of all alleles needed to compute the empirical basecall-like error for InDels. ");
    ADD_OPTDEF2(app, fam_thres_emperr_con_perc_snv,
        "Mininum percent of read support of consensus allele needed to compute the empirical basecall error for SNVs. ");
    ADD_OPTDEF2(app, fam_thres_emperr_con_perc_indel, 
        "Mininum percent of read support of consensus alllele needed to compute the empirical basecall-like error for InDels. ");
    
    ADD_OPTDEF2(app, fam_min_n_copies,
        "Mininum number of DNA copies (1 ng contains approximately 300 copies of human genome) to have zero penalty for UMI-labeled barcode-family (PMC6197739 and PMC5856404). ");
    ADD_OPTDEF2(app, fam_min_n_copies_DPxAD,
        "Mininum UMI-consensus (DP*AD) above which there is zero penalty for UMI-labeled barcode-family even if fam-min-n-copies is not satisfied. ");
    ADD_OPTDEF2(app, fam_min_overseq_perc,
        "Mininum percentage of over-sequencing (one plus average family size) to have zero penalty for UMI-labeled barcode-family. ");
    ADD_OPTDEF2(app, fam_bias_overseq_perc,
        "Mininum percentage of over-sequencing (one plus average family size) to trigger the computation of biases for tier-2 single-strand consensus-sequence (SSCS) families. ");
    ADD_OPTDEF2(app, fam_tier3DP_bias_overseq_perc,
        "Mininum percentage of over-sequencing (one plus average family size) to enable the computation of biases of the discordance between tier-2 and tier-3 single-strand consensus-sequence (SSCS) families allele fractions. ");
    ADD_OPTDEF2(app, fam_indel_nonUMI_phred_dec_per_fold_overseq,
        std::string("Phred decrease in non-UMI variant quality for InDels per fold over-sequencing above (")
        + REPLACE_UNDERSCORE_WITH_HYPHEN(fam_thres_emperr_all_flat_indel) + " + 1). ");
    
    ADD_OPTDEF2(app, fam_phred_indel_inc_before_barcode_labeling,
        "Phred-scaled fraction of InDel errors that occurred before the attachment of UMI single-strand barcodes (PMC3111315). ");
    ADD_OPTDEF2(app, fam_phred_sscs_transition_CG_TA, 
        "Maximum Phred score for single-strand consensus sequences (SSCSs) for the C > T transition which is mainly caused by cytosine deanimation (PMC3437896). ");
    ADD_OPTDEF2(app, fam_phred_sscs_transition_AT_GC, 
        "Maximum Phred score for single-strand consensus sequences (SSCSs) for the A:T > G:C transition (PMC3437896). ");
    ADD_OPTDEF2(app, fam_phred_sscs_transversion_CG_AT, 
        "Maximum Phred score for single-strand consensus sequences (SSCSs) for the G > T transversion which is mainly caused by 8-oxo-guanine (PMC3437896). ");

    ADD_OPTDEF2(app, fam_phred_sscs_transversion_other, 
        "Maximum Phred score for single-strand consensus sequences (SSCSs) for any other transversion (PMC3437896). ");
    ADD_OPTDEF2(app, fam_phred_sscs_indel_open, 
        "Maximum Phred score for single-strand consensus sequences (SSCSs) for the opening of InDel gap "
        "(the opening includes the insertion/deletion of one base). ");
    ADD_OPTDEF2(app, fam_phred_sscs_indel_ext, 
        "Maximum Phred score for single-strand consensus sequences (SSCSs) for the extension of InDel gap "
            "(excluding the extension of one base) per base. ");
    ADD_OPTDEF2(app, fam_phred_dscs_all, 
        "Maximum Phred score for double-strand consensus sequences (DSCSs) for all types of mutations (PMC3437896). ");
    ADD_OPTDEF2(app, fam_phred_dscs_inc_max, 
        "Maximum increase in variant quality by merging SSCSs into DSCSs. "
        "This value defaults to the power-law exponent times the increase in variant quality per DSCS). ");
    
    ADD_OPTDEF2(app, fam_phred_pow_sscs_transversion_AT_TA_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for the A:T > T:A type of SNV (PMC3973132). "
        "It should be lower if the mutations are known to be caused by aristolochic acid (PMC3973132). ");
    ADD_OPTDEF2(app, fam_phred_pow_sscs_snv_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for all other types of SNVs (PMC3437896). ");
    ADD_OPTDEF2(app, fam_phred_pow_sscs_indel_origin, 
        "The phred-score that is subtracted from phred-sscs to get the power-law quality adjustment for InDels (PMC3287198). "
        "Note: somatic expansion/contraction naturally occurs in vivo with a probablity of approximately 0.001 per STR locus per cell generation "
        "(pubmed.ncbi.nlm.nih.gov/9111867 : 1e-2, PMC5054066 : 1e-5 to 1e-3, PMID 8401493: 1e-3). "
        "Hence, STR calls with low UMI allele-fractions that are not in germline high-confidence regions are not necessarily false positive. ");
    ADD_OPTDEF2(app, fam_phred_pow_dscs_all_origin,
        "The phred-score that is subtracted from phred-dscs to get the power-law quality adjustment (PMC3437896). ");
    ADD_OPTDEF2(app, fam_flag,
        "Additional flag. "
        "0x1 sets maximum values to base-call quality (BQ) with the --fam-phred-sscs-* parameters. "
        "0x2 disables the requirement of having UMI labels to consensus duplicated reads into SSCS families. ");

// *** 08. parameters related to systematic errors
    
    ADD_OPTDEF2(app, syserr_BQ_prior, 
        "Phred-scaled probability of the systematic error (defined at PMC3295828). "); 
    
    ADD_OPTDEF2(app, syserr_BQ_sbratio_q_add,
        "Squared relative deviation for the odds ratio of strand read support by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF2(app, syserr_BQ_sbratio_q_max,
        "The maximum increased quality for the odds ratio of strand read support by which the threshold of low-quality bases is increased. ");
    ADD_OPTDEF2(app, syserr_BQ_xmratio_q_add,
        "Squared relative deviation for the odds ratio of mismatch by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF2(app, syserr_BQ_xmratio_q_max,
        "The maximum increased quality for the odds ratio of mismatch by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF2(app, syserr_BQ_bmratio_q_add,
        "Squared relative deviation for the odds ratio of same-base mismatch by which the threshold of low-quality bases is increased. "); 
    ADD_OPTDEF2(app, syserr_BQ_bmratio_q_max,
        "The maximum increased quality for the odds ratio of same-base mismatch by which the threshold of low-quality bases is increased. "); 
    
    ADD_OPTDEF2(app, syserr_BQ_strand_favor_mul,
        "The strand with higher average base quality (BQ) is given this much increase in weight to overcome the potential low average BQ on the other strand. "); 

    ADD_OPTDEF2(app, syserr_minABQ_pcr_snv, 
        "Average base quality (ABQ) in deciPhred below which variant quality is capped to ABQ for SNVs and PCR amplicon. "
        "The recommended value is " STRINGIZE(SYSERR_MINABQ_SNV_ILLUMINA) " for " PLAT_ILLUMINA_LIKE " and 0 for " PLAT_ION_LIKE ". "); 
    ADD_OPTDEF2(app, syserr_minABQ_pcr_indel, 
        "Average base quality (ABQ) in deciPhred below which variant quality is capped to ABQ for InDels and PCR amplicon. "
        "The recommended value is " STRINGIZE(SYSERR_MINABQ_INDEL_ILLUMINA) " for " PLAT_ILLUMINA_LIKE " and 0 for " PLAT_ILLUMINA_LIKE ". ");
    ADD_OPTDEF2(app, syserr_minABQ_cap_snv, 
        "Average base quality (ABQ) in deciPhred below which variant quality is capped to ABQ for SNVs and hybrid selection. "
        "The recommended value is " STRINGIZE(SYSERR_MINABQ_SNV_ILLUMINA) " for " PLAT_ILLUMINA_LIKE " and 0 for " PLAT_ILLUMINA_LIKE ". ");
    ADD_OPTDEF2(app, syserr_minABQ_cap_indel, 
        "Average base quality (ABQ) in deciPhred below which variant quality is capped to ABQ for InDels and hybrid selection. "
        "The recommended value is " STRINGIZE(SYSERR_MINABQ_INDEL_ILLUMINA) " for " PLAT_ILLUMINA_LIKE " and 0 for " PLAT_ILLUMINA_LIKE ". ");
    
    ADD_OPTDEF2(app, syserr_mut_region_n_bases,
        "The mutations that are within this number of bases from each other are assumed to be caused by the same mutation event. ");
    
    ADD_OPTDEF2(app, syserr_MQ_max,
        "Maximum mapping quality (MAPQ or MQ) of the aligned reads, highly recommended to be the 60 from BWA. ");
    ADD_OPTDEF2(app, syserr_MQ_min,
        "Minimum mapping quality (MAPQ or MQ) of the aligned reads, which should be approximately the fraction of regions that are not unique in the entire genome. ");
    
    ADD_OPTDEF2(app, syserr_MQ_NMR_expfrac,
        "Expected fraction of XM-affected base positions in the sequencing fragments covering a locus. "
        "By default, a position is affected if a SNV or InDel is found within 11 base pairs (one turn of DNA helix). ");
    ADD_OPTDEF2(app, syserr_MQ_NMR_altfrac_coef,
        "Base and exponent multiplicative factor for the ALT allele for computing the likelihood of variant being true positive. ");
    ADD_OPTDEF2(app, syserr_MQ_NMR_nonaltfrac_coef,
        "Base and exponent multiplicative factor for the non-ALT alleles for computing the likelihood of variant being true positive. ");
    ADD_OPTDEF2(app, syserr_MQ_NMR_pl_exponent,
        "Exponent multiplicative factor for XM-induced systematic error associated with mapping qualities. ");
    ADD_OPTDEF2(app, syserr_MQ_nonref_base,
        "The variant quality induced by zero mapping quality. \n");
    
// *** 09. parameters related to germline vars 
    
    ADD_OPTDEF2(app, germ_hetero_FA,
        "ALT allele fraction of heterozygous germline polymorphism (PMC6587442). ");
    
    ADD_OPTDEF2(app, germ_phred_hetero_snp,
        "Phred-scaled prior probability of 0/1 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
    ADD_OPTDEF2(app, germ_phred_homalt_snp,
        "Phred-scaled prior probability of 1/1 homozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
    ADD_OPTDEF2(app, germ_phred_het3al_snp,     
        "Phred-scaled prior probability of 1/2 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
    
    ADD_OPTDEF2(app, germ_phred_hetero_indel,
        "Phred-scaled prior probability of 0/1 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). "),
    ADD_OPTDEF2(app, germ_phred_homalt_indel,
        "Phred-scaled prior probability of 1/1 homozygous germline polymorphism at a genomic position for InDels (common empirical knowledge). "),
    ADD_OPTDEF2(app, germ_phred_het3al_indel,   
        "Phred-scaled prior probability of 1/2 heterozygous germline polymorphism at a genomic position for SNPs (common empirical knowledge). ");
        
// *** 10. parameters related to tumor-normal-pairs.
    
    ADD_OPTDEF2(app, tn_q_inc_max,
        "Maximum Phred-scaled increase in variant quality by comparing the tumor with its matched normal. Theoretically, it should be Phred-scaled 2 to the power of --powlaw-exponent. ");
    ADD_OPTDEF2(app, tn_q_inc_max_sscs_CG_AT,
        "Maximum Phred-scaled increase in variant quality by comparing the tumor with its matched normal "
        "for UMI families of single strand consensus sequences (SSCSs) for the G > T transversion. ");
    ADD_OPTDEF2(app, tn_q_inc_max_sscs_other,
        "Maximum Phred-scaled increase in variant quality by comparing the tumor with its matched normal "
        "for UMI families of single strand consensus sequences (SSCSs) for other types of mutations. ");
    
    ADD_OPTDEF2(app, tn_syserr_norm_devqual,
        "Phred-scaled decrease in the variant quality subtracted from the tumor if the tumor FA (allele fraction) deviates by +100\%/-50\% from its matched normal FA. ");
    ADD_OPTDEF2(app, tn_is_paired,
        "The boolean with value 0 (false) or 1 (true) indicating if tumor-normal paired sequencing is used. ");
    
    //ADD_OPTDEF2(app, tn_flag,
    //    "Bitwise flag for tumor-normal comparison. "
    //    "The 0x1 bit means: use the matched normal sample to filter out false positive calls on primers to get more true positive calls. ");
    
// *** 11. parameters related to InDels.
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/
    
    ADD_OPTDEF2(app, indel_BQ_max,
        "For InDels, the maximum base quality of the InDel of one base (PMC4719071). ");
    ADD_OPTDEF2(app, indel_str_repeatsize_max,
        "For InDels, the maximum length of a repeating substring for the repeat to be considered as eligible for short tandem repeat (STR) (PMC5054066). ");
    ADD_OPTDEF2(app, indel_vntr_repeatsize_max,
        "For InDels, the maximum length of a repeating substring for the repeat to be considered as eligible for variable-number tandem repeat (VNTR). ");
    ADD_OPTDEF2(app, indel_polymerase_size,
        "For InDels, the number of bases that can be fit within the catalytic unit of the polymerase used for PCR (PMC149199). ");
    ADD_OPTDEF2(app, indel_polymerase_slip_rate,
        "For InDels, the number of bases that can slip when processed by the catalytic unit of the polymerase used for PCR (PMC149199). ");
    ADD_OPTDEF2(app, indel_del_to_ins_err_ratio,
        "For InDels, the number of times that deletion is more likely (PMC149199). ");
    ADD_OPTDEF2(app, indel_adj_tracklen_dist,
        "InDels within this number of bases of each other are considered to be close to each other. ");
    ADD_OPTDEF2(app, indel_adj_indellen_perc,
        "For InDels, the multiplicative percent coefficient of the InDel length used to determine the base-distance threshold. "
        "Other InDels found within this threshold make the InDel at this position more likely to be false positive. ");
    
    ADD_OPTDEF2(app, indel_multiallele_samepos_penal,
        "For InDels, the Phred-scaled penalty of having more than two alleles at the same VCF position. ");
    ADD_OPTDEF2(app, indel_multiallele_diffpos_penal,
        "For InDels, the Phred-scaled penalty of having more than two alelles that are overlapping with each other. ");
    ADD_OPTDEF2(app, indel_multiallele_soma_penal_thres,
        "For InDels, the Phred-scaled penalty threshold above which penalty applies for somatic variants. ");
    ADD_OPTDEF2(app, indel_tetraallele_germline_penal_value,
        "For InDels, the Phres-scale penalty of having more than three alleles at the same VCF position for germline variants. ");
    ADD_OPTDEF2(app, indel_tetraallele_germline_penal_thres,
        "For InDels, the Phred-scaled penalty threshold above which penalty applies for germline variants. ");
    
    ADD_OPTDEF2(app, indel_ins_penal_pseudocount,
        "For InDels, the length of inserted sequence at which at most half of the inserted sequences are erroneous. ");
    
    ADD_OPTDEF2(app, indel_nonSTR_phred_per_base,
        "The Phred quality per each additional base to the sequenced-segment left/right end in non-STR regions. ");
    ADD_OPTDEF2(app, indel_str_phred_per_region,
        "The Phred quality per each additional STR track to the sequenced-segment left/right end. ");
    
    ADD_OPTDEF2(app, indel_filter_edge_dist,
        "For InDels, the number of bases to the sequenced segment left and right ends. "
        "If the minimum of these two numbers is below this threshold, then the InDel support is filtered out. ");

// *** 12. parameters related to contamination
    
    ADD_OPTDEF2(app, contam_any_mul_frac,
        "Multiplicative contamination rate for the fraction of reads generated by any source of contamination in any sample (PMC3167057). ");
    ADD_OPTDEF2(app, contam_t2n_mul_frac,
        "Multiplicative contamination rate for the fraction of tumor reads in the normal (PMC6528031). "); 

// *** 13. parameters related to micro-adjustment (they do not have any clear theory support)
    
    ADD_OPTDEF2(app, microadjust_xm,
        "Average number of mismatches over 150bp above which one-base insertion has quality capped at the base-quality of the inserted base. ");
    ADD_OPTDEF2(app, microadjust_delFAQmax,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");

    ADD_OPTDEF2(app, microadjust_bias_pos_indel_fold,
        "If the number of reads supporting nearby InDels multiplied by this number is more than the number of reads supporting the variant, then the variant is considered to be affected by nearby InDels. ");
    ADD_OPTDEF2(app, microadjust_bias_pos_indel_misma_to_indel_ratio,
        "If the number of reads supporting nearby mismatches multiplied by this number is more than the number of reads supporting the InDel variant, then the InDel variant is considered to be affected by nearby mismatches. ");
    
    ADD_OPTDEF2(app, microadjust_nobias_pos_indel_misma_to_indel_ratio,
        "If the number of reads supporting nearby mismatches multiplied by this number is less than the number of reads supporting the InDel variant, the STR track length is long enough, and the majority of nearby InDels are of this variant type, then this InDel variant is not affected by any position bias. ");
    ADD_OPTDEF2(app, microadjust_nobias_pos_indel_maxlen,
        "The maximum effective InDel length for <--nobias-pos-indel-lenfrac-thres>. ");
    ADD_OPTDEF2(app, microadjust_nobias_pos_indel_bMQ,
        "If the root-mean-square mapping quality (RMS-MQ) of the InDel variant is above this threshold and <--microadjust-nobias-pos-indel-perc> is passed, then this variant is not subject to insert bias. ");
    ADD_OPTDEF2(app, microadjust_nobias_pos_indel_perc,
        "The ratio of XM2 to (aDPff + aDPfr + aDPrf + aDPrr) in VCF tags. ");
    ADD_OPTDEF2(app, microadjust_nobias_strand_all_fold,
        "The ratio of read count on the unbiased strand to the one on the biased strand, above which the strand bias and insert bias are removed. ");
    ADD_OPTDEF2(app, microadjust_refbias_indel_max,
        "Maximum reference bias for InDels. ");
    
    ADD_OPTDEF2(app, microadjust_counterbias_pos_odds_ratio,
        "The ratio of (allele odds computed using reads with no position bias) to (allele odds computed using reads with position bias), "
        "above which rescue by position counterbias can be activated beause the non-ALT alleles are biased. ");
    ADD_OPTDEF2(app, microadjust_counterbias_pos_fold_ratio,
        "The ratio of (number of reads with no position bias) to (number of reads with position bias), "
        "below which rescue by position counterbias can be activated beause sufficient part of sequencing depth is from reads with position bias. ");
    
    ADD_OPTDEF2(app, microadjust_fam_binom_qual_halving_thres,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");
    ADD_OPTDEF2(app, microadjust_fam_lowfreq_invFA,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");
    ADD_OPTDEF2(app, microadjust_ref_MQ_dec_max,
        "The base-alignment quality of a deletion cannot be below this value at 100\% allele fraction. ");
    
    ADD_OPTDEF2(app, microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual_min,
        "If the decrease in MQ due to being in XM-affected region is above this threshold, "
        "then the normal-sample variant quality is always subtracted from the tumor-sample variant quality in the T-vs-N comparison. ");
    ADD_OPTDEF2(app, microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual_max,
        "The maximum effective decrease in MQ in microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual. ");
    
    ADD_OPTDEF2(app, microadjust_near_clip_dist,
        "A clip (soft or hard) affects this number of bases around it by adding depth to reduce allele fraction for amplicon data. "
        "This prevent false positive calls around highly clipped positions. "); 
    
    ADD_OPTDEF2(app, microadjust_longfrag_sidelength_min,
        "Number of bases to fragment end above which the fragment side is considered to be long for increasing mapping quality. ");
    ADD_OPTDEF2(app, microadjust_longfrag_sidelength_max,
        "Fragment side length above which the fragment side is cappped for increasing mapping quality. ");
    ADD_OPTDEF2(app, microadjust_longfrag_sidelength_zeroMQpenalty,
        "The maximum increase in mapping quality per 100 bases on each side due to long fragment side length. ");
    
    ADD_OPTDEF2(app, microadjust_alignment_clip_min_len,
        "If an alignment position has a clip spanning at least this number of bases besides it, "
        "then this alignment is a candidate for non-small variants and is not counted in APDP[10]. ");
    ADD_OPTDEF2(app, microadjust_alignment_clip_min_frac,
        "If a position has at least this fraction of alignments not couted in APDP[10] and satisfies <--microadjust-alignment-clip-min-count>, "
        "then it generates a <ADDITIONAL_INDEL_CANDIDATE> record in the VCF. ");
    ADD_OPTDEF2(app, microadjust_alignment_clip_min_count,
        "If a position has at least this number of alignments not couted in APDP[10] and satisfies <--microadjust-alignment-clip-min-frac>, "
        "then it generates a <ADDITIONAL_INDEL_CANDIDATE> record in the VCF. ");
    ADD_OPTDEF2(app, microadjust_alignment_tracklen_min,
        "If a position has an STR track containing at least this number of bases, "
        "then it generates a <ADDITIONAL_INDEL_CANDIDATE> record in the VCF. ");
    
    ADD_OPTDEF2(app, microadjust_germline_mix_with_del_snv_penalty,
        "Phred penalty applied to a non-germline-level SNV if this SNV is near both germline-level SNV and germline-level InDel. ");
    ADD_OPTDEF2(app, microadjust_padded_deletion_flag,
        "Bitwise flag indicating treatment of padded deletion. "
        "Bits 0x1 and 0x2: padded deletion is ignored in the nucleotide-base consensus for the Illumina/BGI and ThermoFisher/LifeTech/IonTorrent sequencing platforms, respectively. ");
    
    ADD_OPTDEF2(app, microadjust_strand_orientation_absence_DP_fold,
        "The total-sequencing-depth ratio of one strand/read-orientation with high read support to the other strand/read-orientation with low read support, "
        "above which the strand/read-orientation with low read support is considered to be absent for adjusting the single-sample power-law. ");
    ADD_OPTDEF2(app, microadjust_orientation_absence_snv_penalty,
        "Phred penalty applied to the single-sample power-law quality if read orientations are only in one direction for all SNV alleles including the REF. "
        "This case usually occurs in amplicon-sequencing data. "
        "Please note that this case indicates a lack of information to detect read-orientation bias instead of the result of having the bias being detected. ");
    ADD_OPTDEF2(app, microadjust_strand_absence_snv_penalty,
        "Phred penalty applied to the single-sample power-law quality if sequenced strands are only in one direction for all SNV alleles including the REF. "
        "This case usually occurs in amplicon-sequencing data. "
        "Please note that this case indicates a lack of information to detect sequenced-strand bias instead of the result of having the bias being detected. ");
    ADD_OPTDEF2(app, microadjust_dedup_absence_indel_penalty,
        "Phred penalty applied to the single-sample power-law quality if no dedupping is performed for all InDel alleles including the REF. "
        "This case usually occurs in amplicon-sequencing data. ");
    ADD_OPTDEF2(app, microadjust_median_readlen_thres,
        "Mininum median read length below which microadjust-BAQ-per-base-x1024 is applied as an alternative way to compute base-alignment quality (BAQ). ");
    ADD_OPTDEF2(app, microadjust_BAQ_per_base_x1024,
        "Minimum base-alignment quality (BAQ) BAQ per base/1024 between the locus of interest and the end of segment "
        "to rescue potential false negative calls due to increased bias at short read length. ");
    
    ADD_OPTDEF2(app, lib_wgs_min_avg_fraglen,
        "Minimum average fragment length, below which assay type is inferred to be non-WGS (whole-genome sequencing). ");
    ADD_OPTDEF2(app, lib_nonwgs_ad_pseudocount,
        "For non-WGS: pseudocount for ALT allele-depth for deduplicated reads");
    ADD_OPTDEF2(app, lib_nonwgs_clip_penal_min_indelsize,
        "For non-WGS: if this minimum InDel size is reached, then the effective allele fraction of this InDel is reduced by excessive clipping.");
    ADD_OPTDEF2(app, lib_nonwgs_normal_full_self_rescue_fa,
        "For non-WGS: raw normal allele-fraction (FA) above which the corresponding effective FA is always not reduced.");
    ADD_OPTDEF2(app, lib_nonwgs_normal_min_self_rescue_fa_ratio,
        "For non-WGS: normal effective allele-fraction (FA) cannot be reduced below this times the corresponding raw FA.");
    ADD_OPTDEF2(app, lib_nonwgs_normal_add_mul_ad,
        "For non-WGS: normal cross-sample read count is increased by this multiplicative factor. ");
    ADD_OPTDEF2(app, lib_nonwgs_normal_max_rescued_MQ,
        "For non-WGS: additional maximum Phred-scaled increase in the mapping quality part of the variant quality by comparing the tumor with its matched normal. ");
    ADD_OPTDEF2(app, lib_wgs_normal_max_rescued_MQ,
        "For WGS: additional maximum Phred-scaled increase in the mapping quality part of the variant quality by comparing the tumor with its matched normal. ");
    
// *** 14 debugging
    
    ADD_OPTDEF2(app, debug_note_flag,
        "The flag used for advanced debugging. Please do not activate this option in normal production environments.");
    ADD_OPTDEF2(app, debug_warn_min_read_end_ins_cigar_oplen,
        "If an insertion cigar operation (op) is the first or last cigar op in the cigar string and the cigar-op length is above this threshold, "
        "then a warning is generated to stderr. Does not affect the VCF output result. "); 
    ADD_OPTDEF2(app, debug_tid,
        "This option is used for advanced debugging. This option set the tid of the locus for which detailed debug message will be printed to stdout. ");
    ADD_OPTDEF2(app, debug_pos,
        "This option is used for advanced debugging. This option set the pos of the locus for which detailed debug message will be printed to stdout. ");
    
/// *** end

    app.callback([&]() {
        assay_type = (AssayType)assay_type_uint;
        molecule_tag_uint = (MoleculeTag)molecule_tag_uint;
        sequencing_platform = (SequencingPlatform)sequencing_platform_uint;
        pair_end_merge = (PairEndMerge)pair_end_merge_uint;
        
        if (bam_input_fname.compare(OPT_ONLY_PRINT_VCF_HEADER) == 0) {
            return;
        }
        if (bam_input_fname.compare(OPT_ONLY_PRINT_DEBUG_DETAIL) == 0) {
            return;
        }
        check_file_exist(bam_input_fname, "BAM");
        check_file_exist(bam_input_fname + ".bai", "BAM index");
        if (fasta_ref_fname.compare(std::string("NA")) != 0) {
            check_file_exist(fasta_ref_fname, "FASTA");
            check_file_exist(fasta_ref_fname + ".fai", "FASTA index");
        } else {
            fasta_ref_fname = "";
        }
        this->inferred_sequencing_platform = this->selfUpdateByPlatform();
        if (IS_PROVIDED(vcf_tumor_fname)) {
            // vqual -= (double)10; // maybe useful but not now
        }
        if (IS_PROVIDED(bed_in_fname)) {
            bed_region_fname = bed_in_fname;
        }
        
        this->inferred_is_fastq_generated = (this->fam_consensus_out_fastq.size() > 0);
        this->inferred_is_vcf_generated = ((this->fam_consensus_out_fastq.size() == 0) || (vcf_out_pass_fname.size() > 0));
        
        parsing_result_flag = 0;
    });
    
    CLI11_PARSE(app, argc, argv);
    return 0;
}

#endif

