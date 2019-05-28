#ifndef CmdLineArgs_INCLUDED
#define CmdLineArgs_INCLUDED

#include <fstream>
#include "htslib/sam.h"

#include "common.h"
#include "version.h"

#include "CmdLineArgs.hpp"

SequencingPlatform CommandLineArgs::selfUpdateByPlatform() {
    SequencingPlatform inferred_sequencing_platform = this->sequencing_platform;
    if (SEQUENCING_PLATFORM_AUTO == this->sequencing_platform || SEQUENCING_PLATFORM_OTHER == this->sequencing_platform) {
        // unsigned int bqsum = 0;
        unsigned int bqcnt = 0;

        samFile *sam_infile = sam_open(this->bam_input_fname.c_str(), "r"); // AlignmentFile(samfname, "rb")
        if (NULL == sam_infile) {
            fprintf(stderr, "Failed to open the file %s", this->bam_input_fname.c_str());
            exit(-32);
        }
        std::array<unsigned int, 96> phred2countArr = {0};
        bam_hdr_t * samheader = sam_hdr_read(sam_infile);
        bam1_t *b = bam_init1();
        unsigned int indel_len = 0;
        unsigned int countPE = 0;
        unsigned int countSE = 0;
        while (sam_read1(sam_infile, samheader, b) >= 0 && bqcnt <= 100*1000) {
            if (b->core.flag & 0x1) {
                countPE++;
            } else {
                countSE++;
            }
            for (int i = 0; i < b->core.l_qseq; i++) {
                unsigned int bq = bam_get_qual(b)[i];
                phred2countArr[bq]++;
            }
            bqcnt += b->core.l_qseq;
            const uint32_t n_cigar = b->core.n_cigar;
            const uint32_t *cigar =  bam_get_cigar(b);
            for (unsigned int i = 0; i < n_cigar; i++) {
                uint32_t c = cigar[i];
                unsigned int cigar_op = bam_cigar_op(c);
                if (BAM_CINS == cigar_op || BAM_CDEL ==  cigar_op) {
                    unsigned int cigar_oplen = bam_cigar_oplen(c);
                    indel_len += cigar_oplen;
                }
            }
        }
        bam_destroy1(b);
        bam_hdr_destroy(samheader);
        sam_close(sam_infile);
        
        // The 2/3 of max base quality (BQ) as threshold has no theoretical foundation, so commented out
        /*
        unsigned int phredmax = 0;
        for (int i = 0; i < 96; i++) { 
            if (0 < phred2countArr[i]) { 
                phredmax = i; 
            }
        }
        unsigned int phredcut = phredmax * 2 / 3;
        unsigned int phredpass = 0;
        unsigned int phredfail = 0;
        for (int i = 0; i < 96; i++) { 
            if (i >= phredcut) { 
                phredpass += phred2countArr[i]; 
            } else {
                phredfail += phred2countArr[i];
            }
        }
        bool has_enough_pass = phredpass >= phredfail * 3;
        */
        if (0 < countPE) {
            inferred_sequencing_platform = SEQUENCING_PLATFORM_ILLUMINA;
        } else {
            inferred_sequencing_platform = SEQUENCING_PLATFORM_IONTORRENT;
        }
    }
    if (SEQUENCING_PLATFORM_IONTORRENT == inferred_sequencing_platform && SEQUENCING_PLATFORM_OTHER != this->sequencing_platform) {
        bq_phred_added_indel += 0;
        bq_phred_added_misma += 6;
        minABQ += 25;
    }
    if (SEQUENCING_PLATFORM_ILLUMINA == inferred_sequencing_platform && SEQUENCING_PLATFORM_OTHER != this->sequencing_platform) {
        bq_phred_added_indel += 6;
        bq_phred_added_misma += 0;
        minABQ += 0;
    }
    return inferred_sequencing_platform;
}

void check_file_exist(const std::string & fname, const std::string ftype) {
    std::ifstream ifile(fname.c_str());
    if (!(bool)ifile) {
        std::cerr << "The file " << fname << " of type (" << ftype << ") does not exist." << std::endl;
        exit(-4);
    }
}

const std::string stringvec_to_descstring(const std::vector<std::string> & v) {
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
    CLI::App app{(std::string("Universal Variant Caller (UVC) version ") + VERSION_DETAIL)};
    app.add_flag_function("-v,--version", version_cb,   "Show the version of this program (打印此软件版本号).");    
    app.add_option("inputBam",       bam_input_fname,   "Input coordinate-sorted BAM file that is supposed to contain raw reads (按照位置排序的有原始reads的BAM文件).")->required()->check(CLI::ExistingFile);
    app.add_option("--output-bpRES", vcf_output_fname,  "Output bgzipped vcf file in the format of base-pair resolution (每个位置都输出检测信息的VCF输出文件，文件有可能非常大).");
    app.add_option("-o,--output",    vcf_out_pass_fname,"Output bgzipped vcf file in the format of blocked gVCF. (区块gVCF输出文件).", true);
    app.add_option("--tumor-vcf",    vcf_tumor_fname,   "Bgzipped vcf file of tumor sample as input to the bam file of normal sample (肿瘤样本的VCF文件，用于输入).", true);
    
    app.add_option("-f,--fasta",     fasta_ref_fname,   "Input reference fasta file, where the special value NA means not available (FASTA参考基因组).");
    app.add_option("-R,--region",    bed_region_fname,  "Input BED region file which is optional and delimits the genomic regions to call variants from, "
            "if not provided then call variants from all the sequenced genomic regions. "
            "(BED文件，可有可无，如果没提供则从reads覆盖信息自动计算出BED文件).")->check(CLI::ExistingFile);
    app.add_option("-s,--sample",    sample_name,       "Sample name which is optional (样本名称，可有可无).");
    // app.add_option("--primers",      tsv_primer_fname,  "primer files")->check(CLI::ExistingFile);
    
    app.add_option("-q,--vqual",     vqual,             "Minimum variant quality to be present in --out-pass (如果变异质量低于此值，则不输出到--out-pass文件).", true);
    app.add_option("-d,--min-depth", min_depth_thres,   "Minimum depth below which results are fitlered out and therefore not in the output vcf (如果低于此原始深度则在VCF不输入任何结果).", true);
    app.add_option("-D,--min-altdp", min_altdp_thres,   "Minimum depth of ALT below which results are filtered out (如果ALT深度低于此数则不输出结果).", true);
    app.add_option("-t,--threads",   max_cpu_num,       "Number of cpu cores or equivalently threads to use (使用CPU线程的数量).", true);
    app.add_option("--alnlen",       min_aln_len,       "Minimum alignment length below which the alignment is filtered out (如果比对长度低于比值则过滤掉一行的比对结果).", true);
    app.add_option("--mapqual",      min_mapqual,       "Minimum mapping  quality below which the alignment is filtered out (如果比对质量低于此值则过滤掉一行的比对结果).", true);
    app.add_option("--phred-sscs",   phred_max_sscs,    "Maximum phred score for single-strand consensus sequences (SSCSs)", true);
    app.add_option("--phred-dscs",   phred_max_dscs,    "Maximum phred score for double-strand consensus sequences (DSCSs)", true);
    //app.add_option("--platform",     platform,          "Platform or the sequencer that generated the data, which is either illumina or iontorrent."); 
    app.add_option("--minABQ",       minABQ,            "Minimum average base quality below which variant quality is capped to average base quality, "
                   "recommend 25 for Illumina and 0 for IonTorrent "
                   " (如果位点平均碱基质量低于此值则变异质量不会超过平均碱基质量，建议对Illumina用25并且对IonTorrent用0).", true); 
    app.add_option("--bq-phred-added-misma", bq_phred_added_misma, "Additional base-quality phred score added to match and mismatch, recommend 6 for Illumina and BGI.");
    app.add_option("--bq-phred-added-indel", bq_phred_added_indel, "Additional base-quality phred score added to indel and no-indel, recommend 6 for IonTorrent.");
    app.add_option("--should-add-note",      should_add_note,      "Boolean indicating if the program generates more detail in the vcf result file.");

    app.add_option("--is-dup-aware",         is_dup_aware,             "Is aware of the duplicated reads in BAM data. " + stringvec_to_descstring(ASSAY_TYPE_TO_MSG), true);
    unsigned int assay_type_uint = (unsigned int)assay_type;
    unsigned int molecule_tag_uint = (unsigned int)molecule_tag;
    unsigned int sequencing_platform_uint = (unsigned int)sequencing_platform;
    unsigned int pair_end_merge_uint = (unsigned int)pair_end_merge;
    app.add_option("--assy-type",            assay_type_uint,           "Assay type. " + stringvec_to_descstring(ASSAY_TYPE_TO_MSG), true);
    app.add_option("--molecule-tag",         molecule_tag_uint,         "Molecule tag. " + stringvec_to_descstring(MOLECULE_TAG_TO_MSG), true);
    app.add_option("--sequencing-platform",  sequencing_platform_uint,  "Sequencing platform. " + stringvec_to_descstring(SEQUENCING_PLATFORM_TO_MSG), true);
    app.add_option("--pair-end-merge",       pair_end_merge_uint,       "Merge stat. " + stringvec_to_descstring(ASSAY_TYPE_TO_MSG), true);
    
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
        parsing_result_flag = 0;
    });
    
    CLI11_PARSE(app, argc, argv);
    return 0;
}

#endif

