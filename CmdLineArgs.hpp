#ifndef CmdLineArgs_hpp_INCLUDED
#define CmdLineArgs_hpp_INCLUDED

#include <string>
#include "CLI11-1.7.1/CLI11.hpp"
#include "common.h"

struct CommandLineArgs {
    std::string bam_input_fname = ""; // not missing
    std::string vcf_output_fname = "";
    std::string vcf_out_pass_fname = "-";
    std::string fasta_ref_fname = "";
    std::string bed_region_fname = "";
    std::string vcf_tumor_fname = "";
    std::string sample_name = "-";
    std::string tsv_primer_fname = "";
    
    bool disable_dup_read_merge = false;
    AssayType assay_type = ASSAY_TYPE_AUTO;
    MoleculeTag molecule_tag = MOLECULE_TAG_AUTO;
    SequencingPlatform sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    PairEndMerge pair_end_merge = PAIR_END_MERGE_YES;
    // https://www.biostars.org/p/110670/
    uint32_t    min_depth_thres = 4;
    uint32_t    min_altdp_thres = 2;
    uint32_t    min_aln_len = 0;
    uint32_t    min_mapqual = 0; // 40; // from GATK
    uint32_t    max_cpu_num = 8;
    uint32_t    primerlen = 0;
    uint32_t    phred_max_frag_indel_ext = 5;
    uint32_t    phred_max_sscs_transition_CG_TA = 44; // Cytosine deamination into Uracil, especially in FFPE samples, also by UV light radiation, more upstream
    uint32_t    phred_max_sscs_transition_TA_CG = 48; // https://en.wikipedia.org/wiki/DNA_oxidation, DNA synthesis error, more downstream
    uint32_t    phred_max_sscs_transversion_any = 52;
    uint32_t    phred_max_sscs_indel_open = 34;
    uint32_t    phred_max_sscs_indel_ext  = 5;
    uint32_t    phred_dscs_minus_sscs = 10;
    double      vqual = 10;
    //std::string platform = "auto";
    uint32_t    minABQ_pcr_snv = 0;
    uint32_t    minABQ_pcr_indel = 0;
    uint32_t    minABQ_cap_snv = 0;
    uint32_t    minABQ_cap_indel = 0;
    uint32_t    minMQ1 = 40; // from GATK
    uint32_t    bq_phred_added_indel = 0;
    uint32_t    bq_phred_added_misma = 0;
    bool        should_add_note = false;
    uint32_t    phred_germline_polymorphism = 30;
    double      nonref_to_alt_frac_snv = 0.001;
    double      nonref_to_alt_frac_indel = 0.25;
    double      tnq_mult_snv = 1.0;
    double      tnq_mult_indel = 2.0;
    int initFromArgCV(int & parsing_result_flag, SequencingPlatform & inferred_sequencing_platform, int argc, const char *const* argv);
    SequencingPlatform selfUpdateByPlatform(void);
};
#endif
