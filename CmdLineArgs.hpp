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
    std::string sample_name = "-";
    std::string tsv_primer_fname = "";
    
    bool is_dup_aware = true;
    AssayType assay_type = ASSAY_TYPE_AUTO;
    MoleculeTag molecule_tag = MOLECULE_TAG_AUTO;
    SequencingPlatform sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    PairEndMerge pair_end_merge = PAIR_END_MERGE_YES;
    // https://www.biostars.org/p/110670/
    uint32_t    min_depth_thres = 4;
    uint32_t    min_altdp_thres = 2;
    uint32_t    seq_data_type = 0;
    uint32_t    min_aln_len = 0;
    uint32_t    min_mapqual = 40; // from GATK
    uint32_t    max_cpu_num = 8;
    uint32_t    primerlen = 0;
    uint32_t    phred_max_sscs = 50;
    uint32_t    phred_max_dscs = 60;
    double      vqual = 10;
    std::string platform = "auto";
    uint32_t    minABQ = 0;
    uint32_t    bq_phred_added_indel = 0;
    uint32_t    bq_phred_added_misma = 0;
    bool        should_add_note = false;
    int initFromArgCV(int & parsing_result_flag, int argc, const char *const* argv);
    std::string selfUpdateByPlatform(void); 
};
#endif
