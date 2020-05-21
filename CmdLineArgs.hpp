#ifndef CmdLineArgs_hpp_INCLUDED
#define CmdLineArgs_hpp_INCLUDED

#include "common.hpp"
#include "CLI11-1.7.1/CLI11.hpp"
#include <string>
#include <limits.h>

struct CommandLineArgs {
    std::string bam_input_fname = NOT_PROVIDED; // not missing
    std::string vcf_output_fname = NOT_PROVIDED;
    std::string vcf_out_pass_fname = "-";
    std::string fasta_ref_fname = NOT_PROVIDED;
    std::string tier1_target_region = NOT_PROVIDED; // bcftools view -t
    std::string bed_region_fname = NOT_PROVIDED;    // bcftools view -R
    std::string vcf_tumor_fname = NOT_PROVIDED;
    std::string sample_name = "-";
    std::string vc_stats_fname = NOT_PROVIDED;
    std::string bed_out_fname = NOT_PROVIDED;
    std::string bed_in_fname = NOT_PROVIDED;
    std::string tsv_primer_fname = NOT_PROVIDED;
    
    bool somaticGT = true;
    bool is_tumor_format_retrieved = false;
    bool should_let_all_pass = false;
    bool disable_dup_read_merge = false;
    bool  enable_dup_read_vqual = true; // set to false can result in about 10% increase in false positive (FP) variants at around 1 FP/kilobase rate.
    bool disable_duplex = false;
    unsigned int nogap_phred = 8; // obsolete, not used anymore
    AssayType assay_type = ASSAY_TYPE_AUTO;
    MoleculeTag molecule_tag = MOLECULE_TAG_AUTO;
    SequencingPlatform sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    PairEndMerge pair_end_merge = PAIR_END_MERGE_YES;
    unsigned int fixedthresBQ = 20; // to count the number of bases with base quality of more than this value at each position
    unsigned int uni_bias_thres = 180; // sampple FT includes the filter string if the filter value is higher than this value
    unsigned int uni_bias_r_max = (unsigned int)(INT32_MAX); // is used as infinity here
    
    // it is used to decide whether UMI or non-UMI tumor-vs-normal quality should be used
    unsigned int highqual_thres_snv = 44;
    unsigned int highqual_thres_indel = 0; // 44+6;
    double highqual_min_ratio = 2.5;
    
    // https://www.biostars.org/p/110670/
    uint32_t    min_depth_thres = 4;
    uint32_t    min_altdp_thres = 2;
    uint32_t    min_aln_len = 0;
    uint32_t    min_mapqual = 0; // 40; // from GATK
    uint32_t    max_cpu_num = 8;
    uint32_t    primerlen = 0;
    uint32_t    phred_max_frag_indel_ext = 5;
    uint32_t    phred_max_frag_indel_basemax = 34; // 30; // 40-10 // 37; // 42; // 37; // 35; // 25; // 34;
    uint32_t    phred_max_sscs_transition_CG_TA = 44; // Cytosine deamination into Uracil, especially in FFPE samples, also by UV light radiation, more upstream
    uint32_t    phred_max_sscs_transition_TA_CG = 48; // https://en.wikipedia.org/wiki/DNA_oxidation, DNA synthesis error, more downstream
    uint32_t    phred_max_sscs_transversion_any = 52;
    uint32_t    phred_max_sscs_indel_open = 60; // 34;
    uint32_t    phred_max_sscs_indel_ext  = 0; // 5; // 0;  // 5;
    uint32_t    phred_max_dscs_all = 60;
    uint32_t    phred_pow_sscs_origin = 48-38; // 10*log(2.7e-3/1.5e-4)/log(10)*3 = 37.658175 from https://doi.org/10.1073/pnas.1208715109
    uint32_t    phred_pow_sscs_indel_origin = 0; 
    uint32_t    phred_pow_dscs_origin = 0;
    
    uint32_t    phred_umi_indel_dimret_qual = 40;
    double      phred_umi_indel_dimret_fold = 2;
    
    double      vqual = (double)15; // 10; set to 20 for less output
    uint32_t    vad = (uint32_t)INT_MAX;
    //std::string platform = "auto";
    uint32_t    minABQ_pcr_snv = 0;
    uint32_t    minABQ_pcr_indel = 0;
    uint32_t    minABQ_cap_snv = 0;
    uint32_t    minABQ_cap_indel = 0;
    
    double      ess_georatio_dedup_cap = 1.25;
    double      ess_georatio_dedup_pcr = 1.50; // increase to 1.65 does not help in matching empirical variant score
    double      ess_georatio_duped_pcr = 2.00;
    
    uint32_t    minMQ1 = 40; // from GATK
    uint32_t    maxMQ  = 60; // from bwa
    uint32_t    min_edge_dist = 15; // heuristic (may not work well in STR region)
    
    uint32_t    central_readlen = 0; // estimate from the data
    uint32_t    bq_phred_added_indel = 0;
    uint32_t    bq_phred_added_misma = 0;
    bool        should_add_note = false;
    uint32_t    phred_germline_polymorphism = 31; // +5; // 30+3; // https://www.biostars.org/p/6177/ probablity of hetero is 0.8e-3 for non-african, it should be 32 for african.
    uint32_t    phred_triallelic_indel = 30; // +5; // 30+3; // https://www.biostars.org/p/6177/ probablity of hetero is 0.8e-3 for non-african, it should be 32 for african.
    
    // PMC4271055: probablity of germline call error is between 1/100kb and 1/200kb
    
    double      any_mul_contam_frac = 0.02; // 1e-10; 
    double      t2n_mul_contam_frac = 0.02; // 1e-10; // 0.050; // 0.04 * 2.0; // ;
    double      t2n_add_contam_frac = 0.02;
    double      t2n_add_contam_transfrac = 0.0; // 0.02; // 1e-10; // 0.025; // 0.04; // 0.125*1.5;
    
    uint32_t    phred_frac_indel_error_before_barcode_labeling = 23;
    uint32_t    baq_per_aligned_base = 4; 
    
    bool        is_somatic_snv_filtered_by_any_nonref_germline_snv = true;
    bool        is_somatic_indel_filtered_by_any_nonref_germline_indel = true;
    double      illumina_BQ_pow2_div_coef= 12.0;
    double      varqual_per_mapqual = 4.0/3.0;
    
    double      powlaw_exponent = 3.0;
    double      powlaw_anyvar_base = (double)(60+25+5);
    double      syserr_maxqual = (double)(25.0); // PHRED-scaled probability that a candidate of systematic error is actually non-systematic
    double      syserr_norm_devqual = (double)(12.5); // PHRED-scaled likelihood that the observed allele fraction additively deviates from the expected allele fraction by a multiplicative factor of two

    uint32_t    dedup_center_mult = 5;
    uint32_t    dedup_amplicon_count_to_surrcount_frac = 16;
    uint32_t    dedup_yes_umi_2ends_peak_frac = 4;
    uint32_t    dedup_non_umi_2ends_peak_frac = 8;
    
    /*      Reward                                  Penalty         
            
            (B)     N/raw   N/UMI   N/DCS   (D-F)   N/raw   N/UMI   N/DCS
    SNV:    T/raw   1       no      no      T/raw   1       yes yes
            T/UMI   1       1       no      T/UMI   0.5     1   yes
            T/DCS   1       1       1       T/DCS   0       0   1
        
    InDel   (B)     N/raw   N/UMI   N/DCS   (5)     N/raw   N/UMI   N/DCS
            T/raw   1       no      no      T/raw   0       yes yes
            T/UMI   1       1       no      T/UMI   0       1   yes
            T/DCS   1       1       1       T/DCS   0       0   1
    */
    // (1011 11?1 1011 0101 or equivalently 11 13-15 11 5)
    uint32_t    bitflag_InDel_penal_t_UMI_n_UMI = 0xBFB5; // combine the four above hex letters in parentheses from left-to-right then top-to-bottom.
    
    int
    initFromArgCV(int & parsing_result_flag, SequencingPlatform & inferred_sequencing_platform, int argc, const char *const* argv);
    
    SequencingPlatform 
    selfUpdateByPlatform(void);
};

bool
is_bitflag_checked(uint32_t bitflag_InDel_penal_t_UMI_n_UMI, bool is_InDel, bool is_penal, bool is_t_UMI, bool is_n_UMI);

#endif
