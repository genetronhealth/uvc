#ifndef CmdLineArgs_hpp_INCLUDED
#define CmdLineArgs_hpp_INCLUDED

#include "common.hpp"
#include "CLI11-1.7.1/CLI11.hpp"
#include <string>

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
     
    bool is_tumor_format_retrieved = false;
    bool should_let_all_pass = false;
    bool disable_dup_read_merge = false;
    bool  enable_dup_read_vqual = false;
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
    uint32_t    phred_pow_sscs_origin = 44-26; // 44-30;
    uint32_t    phred_pow_dscs_origin = 0;
    
    double      vqual = (double)15; // 10; set to 20 for less output
    uint32_t    vad = 20; 
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
    
    int 
    initFromArgCV(int & parsing_result_flag, SequencingPlatform & inferred_sequencing_platform, int argc, const char *const* argv);
    
    SequencingPlatform 
    selfUpdateByPlatform(void);
};
#endif
