#ifndef CmdLineArgs_hpp_INCLUDED
#define CmdLineArgs_hpp_INCLUDED

#include <string>
#include "CLI11-1.7.1/CLI11.hpp"
#include "common.h"

struct CommandLineArgs {
    std::string bam_input_fname = NOT_PROVIDED; // not missing
    std::string vcf_output_fname = NOT_PROVIDED;
    std::string vcf_out_pass_fname = "-";
    std::string fasta_ref_fname = NOT_PROVIDED;
    std::string tier1_target_region = NOT_PROVIDED; // bcftools view -t
    std::string bed_region_fname = NOT_PROVIDED;    // bcftools view -R
    std::string vcf_tumor_fname = NOT_PROVIDED;
    std::string sample_name = "-";
    std::string tsv_primer_fname = NOT_PROVIDED;
    
    bool is_tumor_format_retrieved = false;
    bool should_let_all_pass = false;
    bool disable_dup_read_merge = false;
    bool  enable_dup_read_vqual = true;
    bool disable_duplex = false;
    unsigned int nogap_phred = 8; // obsolete, not used anymore
    AssayType assay_type = ASSAY_TYPE_AUTO;
    MoleculeTag molecule_tag = MOLECULE_TAG_AUTO;
    SequencingPlatform sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    PairEndMerge pair_end_merge = PAIR_END_MERGE_YES;
    unsigned int fixedthresBQ = 20; // to count the number of bases with base quality of more than this value at each position
    unsigned int uni_bias_thres = 180; // sampple FT includes the filter string if the filter value is higher than this value
    unsigned int uni_bias_r_max = 1900900900; // is used as infinity here
    double diffVAQfrac = 0; // set to 1 to get the old behavior
    
    // it is used to decide whether UMI or non-UMI tumor-vs-normal quality should be used
    unsigned int highqual_thres_snv = 44;
    unsigned int highqual_thres_indel = 0; // 44+6;
    double highqual_min_ratio = 2.5;
    // removed two variables because too hard to use
    //unsigned int highqual_min_vardep = 3;
    //unsigned int highqual_min_totdep = 500;
    
    // https://www.biostars.org/p/110670/
    uint32_t    min_depth_thres = 4;
    uint32_t    min_altdp_thres = 2;
    uint32_t    min_aln_len = 0;
    uint32_t    min_mapqual = 0; // 40; // from GATK
    uint32_t    max_cpu_num = 8;
    uint32_t    primerlen = 0;
    uint32_t    phred_max_frag_indel_ext = 5;
    uint32_t    phred_max_frag_indel_basemax = 34;
    uint32_t    phred_max_sscs_transition_CG_TA = 44; // Cytosine deamination into Uracil, especially in FFPE samples, also by UV light radiation, more upstream
    uint32_t    phred_max_sscs_transition_TA_CG = 48; // https://en.wikipedia.org/wiki/DNA_oxidation, DNA synthesis error, more downstream
    uint32_t    phred_max_sscs_transversion_any = 52;
    uint32_t    phred_max_sscs_indel_open = 50; // 34;
    uint32_t    phred_max_sscs_indel_ext  = 0;  // 5;
    uint32_t    phred_dscs_minus_sscs = 10;
    double      vqual = (double)15; // 10; set to 20 for less output
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
    uint32_t    bq_phred_added_indel = 0;
    uint32_t    bq_phred_added_misma = 0;
    bool        should_add_note = false;
    uint32_t    phred_germline_polymorphism = 30; // probablity of germline polymorphism is between 1/500 and 1/1kb
    //uint32_t    phred_sys_bias = 0;
    uint32_t    phred_sys_artifact_snv   = phred_germline_polymorphism * 2; // or 55; // PMC4271055: probablity of germline call error is between 1/100kb and 1/200kb
    uint32_t    phred_sys_artifact_indel = phred_germline_polymorphism * 2; // 3 / 2;
    double      nonref_to_alt_frac_snv   = 0.50; // 0.50 for practically removing tri-allelic sites.
    double      nonref_to_alt_frac_indel = 0.20;
    double      tnq_mult_snv   = 1.0; // 2.5; // 0.05;
    double      tnq_mult_indel = 1.0; // 2.5; // 0.05; // * 1.5;
    double      tnq_mult_tADadd_snv   = 4.000; // not used anymore
    double      tnq_mult_tADadd_indel = 4.000; // not used anymore // * 1.5;
    
    double      ldi_tier_qual = 1e9; // always disabled // 0; // strongly enabled ; // 20;
    uint32_t    ldi_tier1cnt  = 150; //     disabled 300;  
    uint32_t    ldi_tier2cnt  = 100; //     disabled weakly enabled with add-one smoothing
    double      mai_tier_qual = 0; // 20;  // // 20;  // enabled // 40; // 40; // = 40; // PHRED-scale probability that an InDel locus is noisy
    uint32_t    mai_tier1abq  = 1.0; //     highly enabled // 40; // = 40; // approximately one extra indel is added as a pseudocount
    uint32_t    mai_tier2abq  = 0.0; //     always disabled // 1024*1024*1024; // disabled
    double      str_tier_qual = 1e9; // always disabled 50; // 50; // 45; // = 50; // is approximately the same as phred_sys_artifact_indel
    uint32_t    str_tier1len  = 15;  //     disabled // = 16; // critical STR region size at which polymerase slippage error reaches a plateau
    uint32_t    str_tier2len  = 15;  //     disabled // enabled 
    
    double      contam_est_qual_thres = 50.0; // used for bootstrapping
    double      contam_ratio = 0.0; // can be estimated by bootstrapping
    double      sys_to_nonsys_err_ratio = 1.0; // can be estimated by bootstrapping
    
    int 
    initFromArgCV(int & parsing_result_flag, SequencingPlatform & inferred_sequencing_platform, int argc, const char *const* argv);
    
    SequencingPlatform 
    selfUpdateByPlatform(void);
};
#endif
