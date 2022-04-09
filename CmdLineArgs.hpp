#ifndef CmdLineArgs_hpp_INCLUDED
#define CmdLineArgs_hpp_INCLUDED

#include "common.hpp"
#include "CLI11-1.7.1/CLI11.hpp"
#include <string>
#include <limits.h>
#include <math.h>

#define DEBUG_NOTE_FLAG_BITMASK_BAQ_OFFSETARR 0x1

struct CommandLineArgs {

// *** 00. frequently used parameters
    
    std::string bam_input_fname = NOT_PROVIDED; // not missing
    std::string fasta_ref_fname = NOT_PROVIDED;
    std::string vcf_out_pass_fname = "-";
    
    std::string bed_region_fname = NOT_PROVIDED;    // bcftools view -R
    std::string tier1_target_region = NOT_PROVIDED; // bcftools view -t
    std::string sample_name = "-";
    
    size_t      max_cpu_num = 8;
    size_t      mem_per_thread = (1024 * 3 / 2); // MegaBytes, for samtools sort 1.11 the value is 768. Variant calling is more complex and mem-consuming than sorting.
    
    uvc1_flag_t outvar_flag = OUTVAR_SOMATIC + OUTVAR_ANY + OUTVAR_MGVCF + OUTVAR_BASE_NN + OUTVAR_ADDITIONAL_INDEL_CANDIDATE;
    bool        should_output_all = false;
    bool        should_output_all_germline = false;
    double      vqual = (double)15; // set to 20 for less output
    
    AssayType assay_type = ASSAY_TYPE_AUTO;
    
    uvc1_qual_t        fam_thres_highBQ_snv = 25;
    uvc1_qual_t        fam_thres_highBQ_indel = 13;
    uvc1_readnum_t     fam_thres_dup1add = 2;
    uvc1_readnum100x_t fam_thres_dup1perc = 80;
    uvc1_readnum_t     fam_thres_dup2add = 3;
    uvc1_readnum100x_t fam_thres_dup2perc = 70; // 85 for more consensus specificity
    std::string        fam_consensus_out_fastq = "";
    
// *** 01. parameters of the names of files, samples, regions, etc.
    
    std::string vcf_tumor_fname = NOT_PROVIDED;
    std::string bed_out_fname = NOT_PROVIDED;
    std::string bed_in_fname = NOT_PROVIDED;
    uvc1_readnum_t bed_in_avg_sequencing_DP = -1; // infer from input BAM data
    
// *** 02. parameters that control input, output, and logs (driven by computational requirements and resources)
    
    bool is_tumor_format_retrieved = true;
    
    // https://www.biostars.org/p/110670/
    
    uvc1_readpos_t    kept_aln_min_aln_len = 0;
    uvc1_qual_t       kept_aln_min_mapqual = 0; // 40; // from GATK
    uvc1_readpos_t    kept_aln_min_isize = 0;
    uvc1_readpos_t    kept_aln_max_isize = INT32_MAX;
    bool              kept_aln_is_zero_isize_discarded = false;
    
    uvc1_readnum_t    min_altdp_thres = 2;
    
    uvc1_readnum_t     vdp1 = 1000;
    uvc1_readnum_t     vad1 = 4;
    double             vfa1 = 0.002;
    uvc1_readnum_t     vdp2 = 10000;
    uvc1_readnum_t     vad2 = 8;
    double             vfa2 = 0.0002;
    
    uvc1_readnum_t     min_r_ad = 0;
    uvc1_readnum_t     min_a_ad = 0;
    
    bool        should_add_note = false;
    bool        always_log = false;
   // *** 03. parameters that are driven by the properties of the assay
    
    MoleculeTag molecule_tag = MOLECULE_TAG_AUTO;
    SequencingPlatform sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    
    // NOTE: these two inferred parameters are not shown (and cannot be passed in as parameters) on the command-line.
    SequencingPlatform inferred_sequencing_platform = sequencing_platform;
    uvc1_qual_t        inferred_maxMQ = 0;
    
    PairEndMerge pair_end_merge = PAIR_END_MERGE_YES;
    bool              disable_duplex = false;
    uvc1_readpos_t    primerlen = 0; // 23; // https://link.springer.com/chapter/10.1007/978-1-4020-6241-4_5 : 18 - 22 bps
    uvc1_readpos_t    primerlen2 = 23; // https://genome.cshlp.org/content/3/3/S30.full.pdf : 18 - 24 bps
    uvc1_flag_t       primer_flag = 0x0; 
    uvc1_readpos_t    central_readlen = 0; // estimate from the data
    uvc1_qual_t       bq_phred_added_misma = 0; // estimate from the data
    uvc1_qual_t       bq_phred_added_indel = 0; // estimate from the data
    
    // http://snap.stanford.edu/class/cs224w-2015/slides/04-powerlaws.pdf
    // https://cs.brynmawr.edu/Courses/cs380/spring2013/section02/slides/10_ScaleFreeNetworks.pdf
    double           powlaw_exponent = 3.0; // universality constant
    double           powlaw_anyvar_base = (double)(60+25+5); // universality constant
    double           powlaw_amplicon_allele_fraction_coef = (5.0/8.0); // universality constant
    
    uvc1_qual_t      penal4lowdep = 37;
    uvc1_qual_t      assay_sequencing_BQ_max = 37;
    uvc1_qual_t      assay_sequencing_BQ_inc = 0;

    uvc1_readnum_t   phasing_haplotype_max_count = 8;
    uvc1_readnum_t   phasing_haplotype_min_ad = 1;
    uvc1_readnum_t   phasing_haplotype_max_detail_cnt = 3;
    
// *** 04. parameters for dedupping reads
    
    // PCR stutter noise at (di,tri,tetra,...)-nucleotide generates +-(2,3,4...) shift in read end position, 
    // so more accurate dedupping requires us to consider these cases. This constant is good enough for the general case.
    double           dedup_center_mult = 5; 
    //double           dedup_amplicon_count_to_surrcount_ratio = 16; // can be 20
    //double           dedup_amplicon_count_to_surrcount_ratio_twosided = 4; // can be 5 or 6
    double           dedup_amplicon_end2end_ratio = 1.5;
    double           dedup_amplicon_border_to_insert_cov_weak_avgDP_ratio = 4; // *1.5;
    double           dedup_amplicon_border_to_insert_cov_strong_avgDP_ratio = 16; // *1.5;
    
    uvc1_flag_t      dedup_flag = 0x0;
    
// *** 05. parameters related to bias thresholds
    
    uvc1_qual_t      bias_thres_highBQ = 20;
    uvc1_qual_t      bias_thres_highBAQ = 20; // is 20+3 for SNVs
    
    uvc1_readpos_t   bias_thres_aLPxT_add = 5;
    uvc1_readpos_t   bias_thres_aLPxT_perc = 160;
    
#if COMPILATION_ENABLE_XMGOT
    uvc1_base1500x_t   bias_thres_PFXM1T_add = 130; // 35; // set very high to disable mismatch bias
    uvc1_base1500x_t   bias_thres_PFXM2T_add = 20;
    uvc1_base1500x_t   bias_thres_PFGO1T_add = 125; // set very high to disable gap bias
    uvc1_base1500x_t   bias_thres_PFGO2T_add = 15;
    
    uvc1_readnum100x_t bias_thres_PFXM1T_perc = 50;
    uvc1_readnum100x_t bias_thres_PFXM2T_perc = 70;
    uvc1_readnum100x_t bias_thres_PFGO1T_perc = 50;
    uvc1_readnum100x_t bias_thres_PFGO2T_perc = 70;
    
    uvc1_readnum100x_t bias_thres_PFXM1NT_perc = 70; // for normal
    uvc1_readnum100x_t bias_thres_PFGO1NT_perc = 70; // for normal
#endif
    
    uvc1_readpos_t     bias_thres_aLRP1t_minus = 10;
    uvc1_readpos_t     bias_thres_aLRP2t_minus = 5;
    uvc1_readpos_t     bias_thres_aLRB1t_minus = 50;
    uvc1_readpos_t     bias_thres_aLRB2t_minus = 25;
    
    uvc1_readnum100x_t bias_thres_aLRP1t_avgmul_perc = 100;
    uvc1_readnum100x_t bias_thres_aLRP2t_avgmul_perc = 100;
    uvc1_readnum100x_t bias_thres_aLRB1t_avgmul_perc = 100;
    uvc1_readnum100x_t bias_thres_aLRB2t_avgmul_perc = 100;
    
    uvc1_readnum100x_t bias_thres_aLRP1Nt_avgmul_perc = 80; // for normal
    uvc1_readnum100x_t bias_thres_aLRB1Nt_avgmul_perc = 80; // for normal

    uvc1_readnum100x_t bias_thres_aLRI1T_perc = 200;
    uvc1_readnum100x_t bias_thres_aLRI2T_perc = 150;
    uvc1_readnum100x_t bias_thres_aLRI1t_perc = 50;
    uvc1_readnum100x_t bias_thres_aLRI2t_perc = 67;

    uvc1_readnum100x_t bias_thres_aLRI1NT_perc = 250; // for normal
    uvc1_readnum100x_t bias_thres_aLRI1Nt_perc = 40;  // for normal
    
    uvc1_readpos_t     bias_thres_aLRI1T_add = 180; // 200;
    uvc1_readpos_t     bias_thres_aLRI2T_add = 150;
    
    uvc1_qual_t        bias_thres_PFBQ1 = 25;
    uvc1_qual_t        bias_thres_PFBQ2 = 30;
    
    uvc1_base1500x_t   bias_thres_aXM1T_add = 30;
    
    uvc1_readpos_t     bias_thres_interfering_indel = 5;
    uvc1_qual_t        bias_thres_interfering_indel_BQ = 21;
    uvc1_qual_t        bias_thres_BAQ1 = 23;
    uvc1_qual_t        bias_thres_BAQ2 = 33;
    
    double             bias_thres_FTS_FA = 2.0;
    bool               bias_is_orientation_artifact_mixed_with_sequencing_error = false;
// *** 06. parameters related to the priors of bias
    
    uvc1_readnum100x_t bias_prior_DPadd_perc = 50;
   
    double         bias_priorfreq_pos = 40; // set very high to disable position bias, insert-end bias, strand bias, and orientation bias.
    double         bias_priorfreq_indel_in_read_div = 20;
    double         bias_priorfreq_indel_in_var_div2 = 15;
    double         bias_priorfreq_indel_in_str_div2 = 10;
    double         bias_priorfreq_var_in_str_div2 = 5;
    
    double         bias_prior_var_DP_mul = 1.25 + DBLFLT_EPS;
    
    uvc1_qual_t    bias_priorfreq_ipos_snv = 60-15; // set very high to disable insert-end bias
    uvc1_qual_t    bias_priorfreq_ipos_indel = 60-15;
    uvc1_qual_t    bias_priorfreq_strand_snv_base = 10; // set very high to disable strand bias
    uvc1_qual_t    bias_priorfreq_strand_indel = 60-15;
    
    double         bias_FA_pseudocount_indel_in_read = 0.5/10.0;
    
    double         bias_priorfreq_orientation_snv_base = 60-15; // set very high to disable orientation bias
    double         bias_priorfreq_orientation_indel_base = 60-15;
    double         bias_orientation_counter_avg_end_len = 20;
    
    uvc1_qual_t    bias_FA_powerlaw_noUMI_phred_inc_snv = 5;
    uvc1_qual_t    bias_FA_powerlaw_noUMI_phred_inc_indel = 7; // this is actually the intrinsic lower error rate of indel instead of the one after reduction by bias.
    uvc1_qual_t    bias_FA_powerlaw_withUMI_phred_inc_snv = 5+3;
    uvc1_qual_t    bias_FA_powerlaw_withUMI_phred_inc_indel = 7;
    
    uvc1_flag_t    nobias_flag = 0x2;
    double         nobias_pos_indel_lenfrac_thres = 2.0; // set very low to disable position bias for InDels
    uvc1_readpos_t nobias_pos_indel_str_track_len = 16;
    
// *** 07. parameters related to read families
    
    uvc1_readnum_t      fam_thres_emperr_all_flat_snv = 4;
    uvc1_readnum100x_t  fam_thres_emperr_con_perc_snv = 67;
    uvc1_readnum_t      fam_thres_emperr_all_flat_indel = 4; // can be 5
    uvc1_readnum100x_t  fam_thres_emperr_con_perc_indel = 67; // can be 75
    
    uvc1_readnum_t      fam_min_n_copies = 300 * 3; // 300 DNA copies per nanogram of DNA
    uvc1_readnum100x_t  fam_min_overseq_perc = 250; // percent fold of over-sequencing
    uvc1_readnum100x_t  fam_bias_overseq_perc = 150; // percent fold of over-sequencing
    
    uvc1_readnum100x_t  fam_indel_nonUMI_phred_dec_per_fold_overseq = 9;

    // 10: error of 10 PCR cycles using low-fidelity polymerase, https://www.nature.com/articles/s41598-020-63102-8
    // 13: reduction in error by using high-fidelity polymerase for UMI assay, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3287198/ 
    // https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_7076.pdf
    // uint16_t fam_phred_indel_err_red_by_high_fidelity_pol is 10; // 10 + 13;
    // https://www.nature.com/articles/s41598-018-31064-7 : All libraries included PCR steps totaling 37 cycles. During Step 4, at cycles 21, 23, 25, 27,
    // 14: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3111315/ : Following 25 additional cycles of PCR, There were 19 cycles of PCR 
    uvc1_qual_t         fam_phred_indel_inc_before_barcode_labeling = 14; // can be 13, 13 + 14, or 10 + 13
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3616734/ : all major library-prep artifacts
    uvc1_qual_t         fam_phred_sscs_transition_CG_TA = 40; // Cytosine deamination into Uracil, especially in FFPE samples, also by UV light radiation, more upstream
    uvc1_qual_t         fam_phred_sscs_transition_AT_GC = 44; // https://en.wikipedia.org/wiki/DNA_oxidation, DNA synthesis error, more downstream
    uvc1_qual_t         fam_phred_sscs_transversion_CG_AT = 48; // https://www.pnas.org/content/109/36/14508 : there can be C->A artifact
    uvc1_qual_t         fam_phred_sscs_transversion_other = 48;
    uvc1_qual_t         fam_phred_sscs_indel_open = 58;
    uvc1_qual_t         fam_phred_sscs_indel_ext  = 0;
    uvc1_qual_t         fam_phred_dscs_all = 58;
    uvc1_qual_t         fam_phred_dscs_max = 68; // theoretical max
    uvc1_qual_t         fam_phred_dscs_inc_max = (68-48);
    
    // 10*log((2.7e-3-3.5e-5)/(1.5e-4-3.5e-5))/log(10)*3 is 41 from https://doi.org/10.1073/pnas.1208715109 PMC3437896
    // The -6 is to accommodate for the fact that most BQs are strictly above 30.
    uvc1_qual_t         fam_phred_pow_sscs_transversion_AT_TA_origin = 44-(41-2-6); // A:T > T:A somatic mutations are uncommon
    double              fam_phred_pow_sscs_snv_origin = 44 - (41-6);
    double              fam_phred_pow_sscs_indel_origin = fam_phred_sscs_indel_open - 9 * 3;
    double              fam_phred_pow_dscs_all_origin = 0;
    uvc1_flag_t         fam_flag = 0x0; // 0x1: set cap to BQ. 0x2: disable requirement of having UMI to consensus.

// *** 08. parameters related to systematic errors
    
    uvc1_qual_t         syserr_BQ_prior = 30; // https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451
    
    uvc1_deciphred_t    syserr_BQ_sbratio_q_add = 5; // note: this is 10*phred, or equivalently in deciphred
    uvc1_deciphred_t    syserr_BQ_sbratio_q_max = 40;
    uvc1_deciphred_t    syserr_BQ_xmratio_q_add = 5;
    uvc1_deciphred_t    syserr_BQ_xmratio_q_max = 40;
    uvc1_deciphred_t    syserr_BQ_bmratio_q_add = 5;
    uvc1_deciphred_t    syserr_BQ_bmratio_q_max = 40;
    
    int                 syserr_BQ_strand_favor_mul = 3;
    
    uvc1_deciphred_t    syserr_minABQ_pcr_snv = 0;   // will be inferred from sequencing platform
    uvc1_deciphred_t    syserr_minABQ_pcr_indel = 0; // will be inferred from sequencing platform
    uvc1_deciphred_t    syserr_minABQ_cap_snv = 0;
    uvc1_deciphred_t    syserr_minABQ_cap_indel = 0;
    
    uvc1_refgpos_t      syserr_mut_region_n_bases = 11;

    uvc1_qual_t         syserr_MQ_min = 0; 
    uvc1_qual_t         syserr_MQ_max = 60; // from bwa
    
    double              syserr_MQ_NMR_expfrac = 0.03; // 23/750
    double              syserr_MQ_NMR_altfrac_coef = 2.0; // base and exponent multiplicative factor for the ALT allele
    // SRR7890876_SRR7890881_fp_chr7_100955016_T_C in the MUC3A gene can be a true positive variant
    double              syserr_MQ_NMR_nonaltfrac_coef = 2.0; // base and exponent multiplicative factor for the non-ALT alleles
    double              syserr_MQ_NMR_pl_exponent = 3.0; // power-law exponent for penalty to the the region of high-basecall-quality XM regions.
    double              syserr_MQ_nonref_base = 40; // power-law exponent for penalty to the the region of high-basecall-quality XM regions.
    
    // Make sure that, by default, all variants (which usually include hotspot variants) are found in the vcf output regardless of mapping quality.
     
// *** 09. parameters related to germline vars // PMC4271055: probablity of germline call error is between 1/100kb and 1/200kb

    double      germ_hetero_FA = 0.47;
    
    // https://www.biostars.org/p/6177/ probablity of hetero is 0.8e-3 for non-african, it should be 32 for african.
    uvc1_qual_t germ_phred_hetero_snp = 31; 
    uvc1_qual_t germ_phred_hetero_indel = 41-1;
    uvc1_qual_t germ_phred_homalt_snp = 31+2;
    uvc1_qual_t germ_phred_homalt_indel = 41-1+2;
    uvc1_qual_t germ_phred_het3al_snp = 54+5;
    uvc1_qual_t germ_phred_het3al_indel = 41-1+9;
    
// *** 10. parameters related to tumor-normal-pairs.
    
    uvc1_qual_t tn_q_inc_max = 9;
    uvc1_qual_t tn_q_inc_max_sscs_CG_AT = 0;
    uvc1_qual_t tn_q_inc_max_sscs_other = 5;
    
    // PHRED-scaled likelihood that the observed allele fraction additively deviates from the expected allele fraction by a multiplicative factor of 2
    double      tn_syserr_norm_devqual = 15.0; // can be (double)(12.5);
    uvc1_flag_t tn_is_paired = false;
    //uvc1_flag_t tn_flag = 0x0;

// *** 11. parameters related to InDels.
    
    uvc1_qual_t         indel_BQ_max = 43-1;
    uvc1_readpos_t      indel_str_repeatsize_max = 6;
    uvc1_readpos_t      indel_vntr_repeatsize_max = 35;
    double              indel_polymerase_size = 8.0;
    double              indel_polymerase_slip_rate = 8.0;
    double              indel_del_to_ins_err_ratio = 5.0; // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/ Table 1 homopolymer error
    uvc1_readpos_t      indel_adj_tracklen_dist = 6;
    uvc1_readnum100x_t  indel_adj_indellen_perc = 160;
    
    double      indel_multiallele_samepos_penal = 11.0;
    double      indel_multiallele_diffpos_penal = 8.0;
    double      indel_multiallele_soma_penal_thres = 11.0;
    double      indel_tetraallele_germline_penal_value = 8.0 * 2;
    double      indel_tetraallele_germline_penal_thres = 22.0;
    
    uvc1_readpos_t      indel_ins_penal_pseudocount = 16;
    
    // According to "A New Lossless DNA Compression Algorithm Based on A Single-Block Encoding Scheme" Table 7 Korea2009024, 
    // there is 2*577/800 bits of info per nucleotide for the human genome.
    uvc1_qual_t         indel_nonSTR_phred_per_base = 5;
    // https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1505-2
    uvc1_qual_t         indel_str_phred_per_region = 5*2; // should be 15 but set to 10 to allow some correlation
    uvc1_readpos_t      indel_filter_edge_dist = 5;

    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2734402/#bib41 : powlaw exponent of 1.5-1.6 for mut rate vs indel len.
    // https://pubmed.ncbi.nlm.nih.gov/18641631/ : SNV mutation rate near (up to a few hundred bp) heterozygous InDels are higher than expected.
    
// *** 12. parameters related to contamination
    
    double              contam_any_mul_frac = 0.02; // from the ContEst paper at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3167057/
    double              contam_t2n_mul_frac = 0.05; // from the DeTiN paper at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6528031/ 

// *** 13. parameters related to micro-adjustment (they do not have any clear theory support)
    int32_t             microadjust_xm = 7;
    uvc1_readpos_t      microadjust_cliplen = 5;
    uvc1_qual_t         microadjust_delFAQmax = 10+9+30; // to override 
    
    double              microadjust_bias_pos_indel_fold = 2;
    double              microadjust_bias_pos_indel_misma_to_indel_ratio = 4 * (1.0 - DBL_EPSILON);
    
    double              microadjust_nobias_pos_indel_misma_to_indel_ratio = 4 * (1.0 - DBL_EPSILON);
    uvc1_readpos_t      microadjust_nobias_pos_indel_maxlen = 16;
    uvc1_qual_t         microadjust_nobias_pos_indel_bMQ = 50;
    uvc1_readnum100x_t  microadjust_nobias_pos_indel_perc = 50;
    double              microadjust_nobias_strand_all_fold = 5; // it was reduced from 20
    double              microadjust_refbias_indel_max = 2.0;
    
    double              microadjust_counterbias_pos_odds_ratio = 3.5;
    double              microadjust_counterbias_pos_fold_ratio = 5.0;
    
    uvc1_qual_t         microadjust_fam_binom_qual_halving_thres = 70; // 22; // =x where x means halved effect of read support at one-FP/(10^(x/10)) base-pairs
    int32_t             microadjust_fam_lowfreq_invFA = 1000;
    uvc1_qual_t         microadjust_ref_MQ_dec_max = 15;
    
    uvc1_qual_t         microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual_min = 30;
    uvc1_qual_t         microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual_max = 30+12;
    uvc1_readpos_t      microadjust_near_clip_dist = 2;
    
    uvc1_readpos_t      microadjust_longfrag_sidelength_min = 300; // both sides span at least one exon
    uvc1_readpos_t      microadjust_longfrag_sidelength_max = 600;
    double              microadjust_longfrag_sidelength_zeroMQpenalty = 300;    

    uvc1_readpos_t      microadjust_alignment_clip_min_len = 12;
    double              microadjust_alignment_clip_min_frac = 0.05;
    uvc1_readpos_t      microadjust_alignment_clip_min_count = 2;
    uvc1_readpos_t      microadjust_alignment_tracklen_min = 25;
    
    uvc1_qual_t         microadjust_germline_mix_with_del_snv_penalty = 9;
    uvc1_flag_t         microadjust_padded_deletion_flag = 0x2;
    
    uvc1_readnum_t      microadjust_strand_orientation_absence_DP_fold = 5;
    uvc1_qual_t         microadjust_orientation_absence_snv_penalty = 4;
    uvc1_qual_t         microadjust_strand_absence_snv_penalty = 4;
    uvc1_qual_t         microadjust_dedup_absence_indel_penalty = 8;
    
    uvc1_readpos_t      lib_wgs_min_avg_fraglen = 300;
    double              lib_nonwgs_ad_pseudocount = 0.1;
    uvc1_readpos_t      lib_nonwgs_clip_penal_min_indelsize = 8;
    double              lib_nonwgs_normal_full_self_rescue_fa = 0.1;
    double              lib_nonwgs_normal_min_self_rescue_fa_ratio = 0.2;
    double              lib_nonwgs_normal_add_mul_ad = 1.0;
        
// *** 14. parameters related to debugging in vcf
    uvc1_flag_t         debug_note_flag = 0x0;
    uvc1_readpos_t      debug_warn_min_read_end_ins_cigar_oplen = 16;

// *** extra useful info
    // https://www.biostars.org/p/254467/#254868 : Question: Are these false somatic variants? Visual inspection with IGV
    // How to tell the difference between HDR and kataegis?
// *** end 
    
    int
    initFromArgCV(int & parsing_result_flag, int argc, const char *const* argv);
    
    SequencingPlatform 
    selfUpdateByPlatform(void);
};

#endif
