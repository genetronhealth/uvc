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
    // std::string tsv_primer_fname = NOT_PROVIDED;
    
    uint32_t outvar_flag = OUTVAR_SOMATIC + OUTVAR_ANY + OUTVAR_NONREF; // 4 anyvar, 2 somatic, 1 germline
    bool somaticGT = true;
    bool is_tumor_format_retrieved = true;
    bool should_output_all = false;
    bool disable_dup_read_merge = false;
    bool  enable_dup_read_vqual = true; // set to false can result in about 10% increase in false positive (FP) variants at around 1 FP/kilobase rate.
    bool disable_duplex = false;
    bool always_log = false;
    unsigned int nogap_phred = 8; // obsolete, not used anymore
    AssayType assay_type = ASSAY_TYPE_AUTO;
    MoleculeTag molecule_tag = MOLECULE_TAG_AUTO;
    SequencingPlatform sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    PairEndMerge pair_end_merge = PAIR_END_MERGE_YES;
    unsigned int fixedthresBQ = 20; // to count the number of bases with base quality of more than this value at each position
    unsigned int uni_bias_thres = 180; // sampple FT includes the filter string if the filter value is higher than this value
    unsigned int uni_bias_r_max = (unsigned int)(INT32_MAX); // is used as infinity here
    // segment: strand, side, position
    //   duped: position, strand, mismatch, deduplication
    // deduped: position, strand, mismatch, deduplication
    // BIAS_SSEG_END may be useful for discarding special false positive variants
    uint32_t bias_flag_amp_snv = (BIAS_FRAG_DUP | BIAS_FRAG_POS | BIAS_FRAG_STR | BIAS_FRAG_MIS);
    uint32_t bias_flag_amp_indel = (BIAS_FRAG_DUP | BIAS_FRAG_POS | BIAS_FRAG_STR | BIAS_FRAG_MIS);
    uint32_t bias_flag_cap_snv = (BIAS_FRAG_DUP | BIAS_FRAG_POS | BIAS_FRAG_STR | BIAS_FRAG_MIS | BIAS_SSEG_POS | BIAS_SSEG_STR);
    uint32_t bias_flag_cap_indel = (BIAS_FRAG_DUP | BIAS_FRAG_POS | BIAS_FRAG_STR | BIAS_FRAG_MIS | BIAS_SSEG_POS | BIAS_SSEG_STR);
    
    // it is used to decide whether UMI or non-UMI tumor-vs-normal quality should be used
    uint8_t highqual_thres_snv = 44;
    uint8_t highqual_thres_indel = 0; // 44+6;
    double highqual_min_ratio = 2.5;
    
    // https://www.biostars.org/p/110670/
    uint32_t    min_depth_thres = 3; // 4; // 3 is used for germline
    uint32_t    min_altdp_thres = 2;
    uint32_t    min_aln_len = 0;
    uint8_t     min_mapqual = 0; // 40; // from GATK
    uint16_t    max_cpu_num = 8;
    // uint32_t    primerlen = 0;
    uint8_t     phred_max_frag_indel_ext = 5; // unused
    uint8_t     phred_max_frag_indel_basemax = 30; // 39; // 34; // 40-10 // 37; // 42; // 37; // 35; // 25;
    uint8_t     phred_max_sscs_transition_CG_TA = 44; // Cytosine deamination into Uracil, especially in FFPE samples, also by UV light radiation, more upstream
    uint8_t     phred_max_sscs_transition_TA_CG = 48; // https://en.wikipedia.org/wiki/DNA_oxidation, DNA synthesis error, more downstream
    uint8_t     phred_max_sscs_transversion_any = 52;
    uint8_t     phred_max_sscs_indel_open = 60; // 34;
    uint8_t     phred_max_sscs_indel_ext  = 0; // 5; // 0;  // 5;
    uint8_t     phred_max_dscs_all = 60;
    double      phred_pow_sscs_origin = 48 - 41; // 10*log((2.7e-3-3.5e-5)/(1.5e-4-3.5e-5))/log(10)*3 = 41 from https://doi.org/10.1073/pnas.1208715109
    double      phred_pow_sscs_indel_origin = 60 - 60; // 60 - 38; 
    double      phred_pow_dscs_origin = 0;
    
    double      phred_umi_dimret_mult_snv = 0.75; //0.4;
    double      phred_umi_dimret_mult_indel = 0.75; // 0.4;
    double      phred_snv_to_indel_ratio = 10.0;
    double      vqual = (double)15; // 10; set to 20 for less output
    uint32_t    vdp = (uint32_t)INT_MAX;
    uint32_t    vad = (uint32_t)INT_MAX;
    uint8_t     minABQ_pcr_snv = 0;
    uint8_t     minABQ_pcr_indel = 0;
    uint8_t     minABQ_cap_snv = 0;
    uint8_t     minABQ_cap_indel = 0;
    
    double      ess_georatio_dedup_cap = 1.25;
    double      ess_georatio_dedup_amp = 1.50; // 1.5 increase to 1.65 does not help in matching empirical variant score
    double      ess_georatio_duped_pcr = 2.00;
    
    uint8_t     minMQ1 = 40; // from GATK
    uint8_t     maxMQ  = 60; // from bwa
    uint8_t     min_edge_dist = 15; // heuristic (may not work well in STR region)
    
    uint8_t     central_readlen = 0; // estimate from the data
    uint8_t     bq_phred_added_indel = 0;
    uint8_t     bq_phred_added_misma = 0;
    bool        should_add_note = false;
    uint8_t     phred_germline_polymorphism = 31; // https://www.biostars.org/p/6177/ probablity of hetero is 0.8e-3 for non-african, it should be 32 for african.
    uint8_t     phred_germline_indel = 40;        // https://www.biostars.org/p/6177/ probablity of hetero is 0.8e-3 for non-african, it should be 32 for african.
    uint8_t     phred_triallelic_snp = 55; // 30;
    uint8_t     phred_triallelic_indel = 48; // 30;
    
    // PMC4271055: probablity of germline call error is between 1/100kb and 1/200kb
    
    double      any_mul_contam_frac = 0.02; // 1e-10; 
    double      t2n_mul_contam_frac = 0.02; // 1e-10; // 0.050; // 0.04 * 2.0; // ;
    double      t2n_add_contam_frac = 0.02;
    double      t2n_add_contam_transfrac = 0.0; // 0.02; // 1e-10; // 0.025; // 0.04; // 0.125*1.5;
    
    uint8_t     phred_frac_indel_error_before_barcode_labeling = 23; // 12, 18, 24 // 23;
    // uint32_t    baq_per_aligned_base = 3; // 4; // The BAQ of 3 seems to make more sense after reading the feedback from reviewer #2
    uint8_t     baq_per_aligned_base = 6; // 5; // + 2; // 7; // 6; // 5; // According to "A New Lossless DNA Compression Algorithm Based on A Single-Block Encoding Scheme" Table 7 Korea2009024, there is 2*577/800 bits of info per nucleotide for the human genome.
    
    bool        is_somatic_snv_filtered_by_any_nonref_germline_snv = true;
    bool        is_somatic_indel_filtered_by_any_nonref_germline_indel = true;
    double      amp_BQ_sqr_coef = 2.0; // 0.09; // 22.0/256.0; 
    double      cap_BQ_sqr_coef = 2.0; // 0.12; // .0/256.0;
    uint8_t     phred_varcall_err_per_map_err_per_base = 20; // HLA pairwise sequence identity is about 1% from https://www.ebi.ac.uk/ipd/mhc/alignment/hla/
    
    double      powlaw_exponent = 3.0;
    double      powlaw_anyvar_base = (double)(60+25+5);
    double      syserr_maxqual = (double)(25.0); // PHRED-scaled probability that a candidate of systematic error is actually non-systematic
    double      syserr_norm_devqual = (double)(12.5); // PHRED-scaled likelihood that the observed allele fraction additively deviates from the expected allele fraction by a multiplicative factor of two

    uint8_t     dedup_center_mult = 5; // PCR stutter noise at (di,tri,tetra,...)-nucleotide generates +-(2,3,4...) shift in read end position, so more accurate dedupping requires us to consider these cases. This constant is good enough for the general case.
    uint8_t     dedup_amplicon_count_to_surrcount_ratio = 16;
    double      dedup_amplicon_end2end_ratio = 1.5;
    
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
    uint32_t    ref_bias_awareness = 0x2; // The 0x1 bit is for amplicon and the 0x2 bit is for non-amplicon
    
    uint32_t haplo_in_diplo_allele_perc = 50;
    uint32_t diplo_oneside_posbias_perc = 20;
    uint32_t diplo_twoside_posbias_perc = 10;
    uint32_t haplo_oneside_posbias_perc = 25;
    uint32_t haplo_twoside_posbias_perc = 15;
    uint32_t regside_nbases = 10; // 30;
    uint32_t dedup_flag = 0x0;
    
    int
    initFromArgCV(int & parsing_result_flag, SequencingPlatform & inferred_sequencing_platform, int argc, const char *const* argv);
    
    SequencingPlatform 
    selfUpdateByPlatform(void);
};

bool
is_bitflag_checked(uint32_t bitflag_InDel_penal_t_UMI_n_UMI, bool is_InDel, bool is_penal, bool is_t_UMI, bool is_n_UMI);

#endif
