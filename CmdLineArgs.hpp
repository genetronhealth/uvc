#ifndef CmdLineArgs_hpp_INCLUDED
#define CmdLineArgs_hpp_INCLUDED

#include "common.hpp"
#include "CLI11-1.7.1/CLI11.hpp"
#include <string>
#include <limits.h>

struct CommandLineArgs {

// *** 00. frequently used parameters
    
    std::string bam_input_fname = NOT_PROVIDED; // not missing
    std::string fasta_ref_fname = NOT_PROVIDED;
    std::string vcf_out_pass_fname = "-";
    
    std::string bed_region_fname = NOT_PROVIDED;    // bcftools view -R
    std::string tier1_target_region = NOT_PROVIDED; // bcftools view -t
    std::string sample_name = "-";
    
    uint16_t    max_cpu_num = 8;
    
    uint32_t    outvar_flag = OUTVAR_SOMATIC + OUTVAR_ANY + OUTVAR_NONREF; // 4 anyvar, 2 somatic, 1 germline
    bool        should_output_all = false;
    bool        should_output_all_germline = false;
    double      vqual = (double)15; // 10; set to 20 for less output
    
    AssayType assay_type = ASSAY_TYPE_AUTO;
    
    uint32_t fam_thres_highBQ_snv = 25;
    uint32_t fam_thres_highBQ_indel = 13;
    uint32_t fam_thres_dup1add = 2;
    uint32_t fam_thres_dup1perc = 80;
    uint32_t fam_thres_dup2add = 3;
    uint32_t fam_thres_dup2perc = 85;
    
// *** 01. parameters of the names of files, samples, regions, etc.
    
    std::string vcf_tumor_fname = NOT_PROVIDED;
    std::string bed_out_fname = NOT_PROVIDED;
    std::string bed_in_fname = NOT_PROVIDED;
    
// *** 02. parameters that control input, output, and logs (driven by computational requirements and resources)
    
    bool somaticGT = true;
    bool is_tumor_format_retrieved = true;
    
    // https://www.biostars.org/p/110670/
    
    uint32_t    min_aln_len = 0;
    uint16_t    min_mapqual = 0; // 40; // from GATK

    uint32_t    min_depth_thres = 3; // 4; // 3 is used for germline
    uint32_t    min_altdp_thres = 2;
    uint32_t    vdp = (uint32_t)INT_MAX;
    uint32_t    vad = (uint32_t)INT_MAX;

    bool        should_add_note = false;
    bool        always_log = false;
       
// *** 03. parameters that are driven by the properties of the assay
    
    MoleculeTag molecule_tag = MOLECULE_TAG_AUTO;
    SequencingPlatform sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    PairEndMerge pair_end_merge = PAIR_END_MERGE_YES;
    bool        disable_duplex = false;
    uint32_t    primerlen = 23;
    
    uint16_t    central_readlen = 0; // estimate from the data
    uint16_t    bq_phred_added_misma = 0; // estiamte from the data
    uint16_t    bq_phred_added_indel = 0; // estimate from the data
    
    double      powlaw_exponent = 3.0; // universality constant
    double      powlaw_anyvar_base = (double)(60+25+5); // universality constant
    
    int32_t     penal4lowdep = 37;

// *** 04. parameters for dedupping reads
    
    // PCR stutter noise at (di,tri,tetra,...)-nucleotide generates +-(2,3,4...) shift in read end position, 
    // so more accurate dedupping requires us to consider these cases. This constant is good enough for the general case.
    uint16_t    dedup_center_mult = 5; 
    uint16_t    dedup_amplicon_count_to_surrcount_ratio = 16;
    uint16_t    dedup_amplicon_count_to_surrcount_ratio_twosided = 4; // 6;
    double      dedup_amplicon_end2end_ratio = 1.5;
    
    uint32_t    dedup_flag = 0x0;
    
// *** 05. parameters related to bias thresholds
    
    uint32_t bias_thres_highBQ = 20;
    uint32_t bias_thres_highBAQ = 23;
    
    uint32_t bias_thres_aLPxT_add = 5;
    uint32_t bias_thres_aLPxT_perc = 160;
        
    uint32_t bias_thres_PFXM1T_add = 35; // set very high to disable mismatch bias
    uint32_t bias_thres_PFXM2T_add = 20;
    uint32_t bias_thres_PFGO1T_add = 25; // set very high to disable gap bias
    uint32_t bias_thres_PFGO2T_add = 15;
    
    uint32_t bias_thres_PFXM1T_perc = 50;
    uint32_t bias_thres_PFXM2T_perc = 70;
    uint32_t bias_thres_PFGO1T_perc = 50;
    uint32_t bias_thres_PFGO2T_perc = 70;
    
    int32_t  bias_thres_aLRP1t_minus = 10;
    int32_t  bias_thres_aLRP2t_minus = 5;
    int32_t  bias_thres_aLRB1t_minus = 50;
    int32_t  bias_thres_aLRB2t_minus = 25;
    
    uint32_t bias_thres_aLRP1t_avgmul_perc = 100;
    uint32_t bias_thres_aLRP2t_avgmul_perc = 100;
    uint32_t bias_thres_aLRB1t_avgmul_perc = 100;
    uint32_t bias_thres_aLRB2t_avgmul_perc = 100;
    
    uint32_t bias_thres_aLRI1T_perc = 200;
    uint32_t bias_thres_aLRI2T_perc = 150;
    uint32_t bias_thres_aLRI1t_perc = 50;
    uint32_t bias_thres_aLRI2t_perc = 67;
    
    uint32_t bias_thres_aLRI1T_add = 180;
    uint32_t bias_thres_aLRI2T_add = 150;
    
    uint32_t bias_thres_PFBQ1 = 25;
    uint32_t bias_thres_PFBQ2 = 30;
    
    uint32_t bias_thres_aXM1T_add = 30;
    
    uint32_t bias_thres_interfering_indel = 5;
    int32_t bias_thres_interfering_indel_BQ = 21;
    uint32_t bias_thres_BAQ1 = 23;
    uint32_t bias_thres_BAQ2 = 33;
    
// *** 06. parameters related to the priors of bias
    
    double   bias_prior_pseudocount = 2; // set very high to disable position, mismatch, and 
    uint32_t bias_prior_DPadd_perc = 50;
   
    double bias_prior_pos = 1024*9/2; // set very high to disable position bias, insert-end bias, strand bias, and orientation bias.
    double bias_prior_indel_in_read_div = 32+64;
    double bias_prior_indel_in_var_div2 = 24;
    double bias_prior_indel_in_STR_div2 = 8;
    double bias_prior_var_in_STR_div2 = 3;
    
    double   bias_prior_var_DP_mul = 1.25 + DBLFLT_EPS;
    
    uint32_t bias_prior_ipos_SNV = 1e5; // set very high to disable insert-end bias
    uint32_t bias_prior_ipos_InDel = 3e3;
    uint32_t bias_prior_strand_SNV_base = 10; // set very high to disable strand bias
    uint32_t bias_prior_strand_InDel = 3e3;
    
    double   bias_FA_pseudocount_indel_in_read = 0.1;
    
    double   nobias_pos_indel_lenfrac_thres = 2.0; // set very low to disable position bias for InDels
    uint32_t nobias_pos_indel_STR_track_len = 16;

    double   bias_prior_orientation_SNV_base = 1e5; // set very high to disable orientation bias
    double   bias_prior_orientation_InDel_base = 1e5; 
    double   bias_orientation_counter_avg_end_len = 20;
    
    int32_t  bias_FA_powerlaw_noUMI_phred_inc_snv = 5;
    int32_t  bias_FA_powerlaw_noUMI_phred_inc_indel = 7; // this is actually the intrinsic lower error rate of indel instead of the one after reduction by bias.
    int32_t  bias_FA_powerlaw_withUMI_phred_inc_snv = 5; // +3; // 7+7;
    int32_t  bias_FA_powerlaw_withUMI_phred_inc_indel = 7; // +3; // 7+7;

// *** 07. parameters related to read families
    
    uint32_t fam_thres_emperr_all_flat_snv = 4;
    uint32_t fam_thres_emperr_con_perc_snv = 67;
    uint32_t fam_thres_emperr_all_flat_indel = 4; // 5;
    uint32_t fam_thres_emperr_con_perc_indel = 67; // 75;
    
    int32_t fam_min_n_copies = 300 * 3; // 300 DNA copies per nanogram of DNA
    int32_t fam_min_overseq_perc = 250; // percent fold of over-sequencing
    
    // 10: error of 10 PCR cycles using low-fidelity polymerase, https://www.nature.com/articles/s41598-020-63102-8
    // 13: reduction in error by using high-fidelity polymerase for UMI assay, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3287198/ 
    // https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_7076.pdf
    // uint16_t fam_phred_indel_err_red_by_high_fidelity_pol = 10; // 10 + 13;
    // https://www.nature.com/articles/s41598-018-31064-7 : All libraries included PCR steps totaling 37 cycles. During Step 4, at cycles 21, 23, 25, 27,
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3111315/ : Following 25 additional cycles of PCR, There were 19 cycles of PCR 
    uint16_t fam_phred_indel_inc_before_barcode_labeling = 14; // 10 + 13;
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3616734/ : all major library-prep artifacts
    uint16_t fam_phred_sscs_transition_CG_TA = 44; // Cytosine deamination into Uracil, especially in FFPE samples, also by UV light radiation, more upstream
    uint16_t fam_phred_sscs_transition_TA_CG = 48; // https://en.wikipedia.org/wiki/DNA_oxidation, DNA synthesis error, more downstream
    uint16_t fam_phred_sscs_transversion_CG_AT = 52; // https://www.pnas.org/content/109/36/14508 : there can be C->A artifact
    uint16_t fam_phred_sscs_transversion_any = 52;
    uint16_t fam_phred_sscs_indel_open = 60;
    uint16_t fam_phred_sscs_indel_ext  = 0;
    uint16_t fam_phred_dscs_all = 60;
    
    double   fam_phred_pow_sscs_SNV_origin = 48 - 41; // 10*log((2.7e-3-3.5e-5)/(1.5e-4-3.5e-5))/log(10)*3 = 41 from https://doi.org/10.1073/pnas.1208715109
    double   fam_phred_pow_sscs_indel_origin = fam_phred_sscs_indel_open - (fam_phred_indel_inc_before_barcode_labeling * 3);
    double   fam_phred_pow_dscs_all_origin = 0;
    
// *** 08. parameters related to systematic errors
    
    uint32_t syserr_BQ_prior = 30; // https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451
    
    uint32_t syserr_BQ_smratio_q_add = 5;
    uint32_t syserr_BQ_smratio_q_max = 40;
    uint32_t syserr_BQ_xmratio_q_add = 5;
    uint32_t syserr_BQ_xmratio_q_max = 30;
    uint32_t syserr_BQ_bmratio_q_add = 5;
    uint32_t syserr_BQ_bmratio_q_max = 40;
    
    int16_t  syserr_minABQ_pcr_snv = 0;   // will be inferred from sequencing platform
    int16_t  syserr_minABQ_pcr_indel = 0; // will be inferred from sequencing platform
    int16_t  syserr_minABQ_cap_snv = 0;
    int16_t  syserr_minABQ_cap_indel = 0;
    
    uint16_t syserr_maxMQ = 60; // from bwa
    uint16_t syserr_phred_varcall_err_per_map_err_per_base = 10; // this is the max phred probability of varcall error per base per mapping error
    
    // Make sure that, by default, all variants (which usually include hotspot variants) are found in the vcf output regardless of mapping quality.
    uint16_t syserr_minMQ = ((vqual > syserr_phred_varcall_err_per_map_err_per_base) ? (vqual - syserr_phred_varcall_err_per_map_err_per_base) : 0); 
 
// *** 09. parameters related to germline vars // PMC4271055: probablity of germline call error is between 1/100kb and 1/200kb

    double      germ_hetero_FA = 0.47;
    
    // https://www.biostars.org/p/6177/ probablity of hetero is 0.8e-3 for non-african, it should be 32 for african.
    uint16_t    germ_phred_hetero_snp = 31; 
    uint16_t    germ_phred_hetero_indel = 41-1;
    uint16_t    germ_phred_homalt_snp = 31+2;
    uint16_t    germ_phred_homalt_indel = 41-1+2;
    uint16_t    germ_phred_het3al_snp = 54+5;
    uint16_t    germ_phred_het3al_indel = 41-1+9;
    
// *** 10. parameters related to tumor-normal-pairs.
    
    uint32_t    tn_q_inc_max = 9;
    // PHRED-scaled likelihood that the observed allele fraction additively deviates from the expected allele fraction by a multiplicative factor of 2
    double      tn_syserr_norm_devqual = 15.0; // (double)(12.5); 
    

// *** 11. parameters related to InDels.
    
    uint32_t    indel_BQ_max = 43-1;
    uint32_t    indel_str_repeatsize_max = 6;
    double      indel_polymerase_size = 8.0;
    double      indel_polymerase_slip_rate = 8.0;
    double      indel_del_to_ins_err_ratio = 5.0; // 4.0; // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/ Table 1 homopolymer error
    uint32_t    indel_adj_tracklen_dist = 6;
    uint32_t    indel_adj_indellen_perc = 160; 

    double      indel_multiallele_samepos_penal = 11.0;
    double      indel_multiallele_diffpos_penal = 8.0;
    double      indel_multiallele_soma_penal_thres = 11.0;
    double      indel_tetraallele_germline_penal_value = 8.0 * 2;
    double      indel_tetraallele_germline_penal_thres = 22.0;
    
    uint32_t    indel_ins_penal_pseudocount = 16;
    
    // According to "A New Lossless DNA Compression Algorithm Based on A Single-Block Encoding Scheme" Table 7 Korea2009024, 
    // there is 2*577/800 bits of info per nucleotide for the human genome.
    uint32_t    indel_nonSTR_phred_per_base = 5;
    uint32_t    indel_STR_phred_per_region = 5*2; // 15, set to 10 to allow some correlation, https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1505-2
    
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2734402/#bib41 : powlaw exponent of 1.5-1.6 for mut rate vs indel len.
    // https://pubmed.ncbi.nlm.nih.gov/18641631/ : SNV mutation rate near (up to a few hundred bp) heterozygous InDels are higher than expected.
    
// *** 12. parameters related to contamination
    
    double      contam_any_mul_frac = 0.02; // from the ContEst paper at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3167057/
    double      contam_t2n_mul_frac = 0.05; // from the DeTiN paper at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6528031/ 

// *** extra useful info
    // https://www.biostars.org/p/254467/#254868 : Question: Are these false somatic variants? Visual inspection with IGV
    // How to tell the difference between HDR and kataegis?
// *** end 
    
    int
    initFromArgCV(int & parsing_result_flag, SequencingPlatform & inferred_sequencing_platform, int argc, const char *const* argv);
    
    SequencingPlatform 
    selfUpdateByPlatform(void);
};

bool
is_bitflag_checked(uint32_t bitflag_InDel_penal_t_UMI_n_UMI, bool is_InDel, bool is_penal, bool is_t_UMI, bool is_n_UMI);

#endif
