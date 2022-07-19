#include <algorithm>
#include <array>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <assert.h>
#include <float.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#include <math.h>
#include <unistd.h>

#define VERSION3 "0.1.5"

#define BCF_NUM_A (-1)
#define BCF_NUM_R (-2)
#define BCF_NUM_G (-3)
#define BCF_NUM_D (-4)

#define BCF_TYPE_FLAG 0
#define BCF_TYPE_INTEGER 1
#define BCF_TYPE_FLOAT 2
#define BCF_TYPE_STRING 3

struct BcfInfo {
    const std::string ID;
    const int number;
    const int type;
    const std::string description;
    
    BcfInfo(std::string aID, int anumber, int atype, std::string adescription) : ID(aID), number(anumber), type(atype), description(adescription) {};
    
    const std::string to_header_string() const {
        std::string string_number = (
                   (BCF_NUM_A == number) ? "A" 
                : ((BCF_NUM_R == number) ? "R" 
                : ((BCF_NUM_G == number) ? "G" 
                : ((BCF_NUM_D == number) ? "D" 
                : (std::to_string(number))))));
        std::string string_type = (
                   (BCF_TYPE_FLAG == type) ? "Flag" 
                : ((BCF_TYPE_INTEGER == type) ? "Integer" 
                : ((BCF_TYPE_FLOAT == type) ? "Float" 
                : ((BCF_TYPE_STRING == type) ? "String" : ""))));
        return std::string("##INFO=<ID=") + ID + ",Number=" + string_number  + ",Type=" + string_type + ",Description=\"" + description + "\">"; 
    }
};

enum BcfInfoTag {
    delinsHap,
    diPRA,
    diADA,
    diDPm,
    diDPM,
    diADm,
    diADM,
    diRDm,
    diRDM,
    diAD2F,
    diCVQ,
    diHVQ
};

const BcfInfo BCF_INFO_LIST[] = {
        
        [delinsHap] = BcfInfo(
        "delinsHap", BCF_NUM_A, BCF_TYPE_STRING,
            "For each delins ALT allele, the value of this tag is the key of the other tag "
            "used by the upstream SNV/InDel caller to represent linked variants that are phased into one single haplotype candidate. "
            "In [x]Hap, [x] can be b, c, or c2 if UVC was used upstream. "),
    
        [diPRA] = BcfInfo(
        "diPRA",     BCF_NUM_A, BCF_TYPE_STRING, 
            "Tumor position_REF_ALT of each delins ALT allele, with the three VCF fields separated by underscore. "),
    
        [diADA] = BcfInfo(
        "diADA",     BCF_NUM_A, BCF_TYPE_INTEGER, 
            "Depth of each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the depth exactly supporting its corresponding delins haplotype. "
            "For example, at read level, only reads supporting the exact delins haplotype are counted. "),

        [diDPm] = BcfInfo(
        "diDPm",     BCF_NUM_A, BCF_TYPE_INTEGER, 
            "Total depth of all alleles located at each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the minimum depth among the SNV and/or InDel alleles that the delins variant is composed of. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. "),
                
        [diDPM] = BcfInfo(
        "diDPM",     BCF_NUM_A, BCF_TYPE_INTEGER,
            "Total depth of all alleles located at each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the maximum depth among the SNV and/or InDel alleles that the delins variant is composed of. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. "),
    
        [diADm] = BcfInfo(
        "diADm",     BCF_NUM_A, BCF_TYPE_INTEGER,
            "Depth of each delins (MVV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the minimum delins-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. "),
    
        [diADM] = BcfInfo(
        "diADM",     BCF_NUM_A, BCF_TYPE_INTEGER,
            "Depth of each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the maximum delins-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. "), // this should rarely be used
                
        [diRDm] = BcfInfo(
        "diRDm",     BCF_NUM_A, BCF_TYPE_INTEGER,
            "Depth of the reference allele corresponding to each delins (MVV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the minimum REF-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. "),
        
        [diRDM] = BcfInfo(
        "diRDM",     BCF_NUM_A, BCF_TYPE_INTEGER,
            "Depth of the reference allele corresponding to each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the maximum REF-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. "), // this should rarely be used
    
        [diAD2F] = BcfInfo(
        "diAD2F",    BCF_NUM_A, BCF_TYPE_INTEGER, 
            "Percentage (between 0 and 100) of reads that support the delins variant "
            "among the individual SNV and/or InDel variants that constitute the delins variant. "),
    
        [diCVQ] = BcfInfo(
        "diCVQ",     BCF_NUM_A, BCF_TYPE_INTEGER, 
            "Combined (by combining SNVs and InDels) variant Quality of each delins ALT variant. "),
    
        [diHVQ] = BcfInfo(
        "diHVQ",     BCF_NUM_A, BCF_TYPE_INTEGER, 
            "Haplotyped (non-SNV and non-InDel small) variant Quality of each delins ALT variant. "),
};


class VariantInfo {
public:
    float qual;
    int tbDP;
    int tDP;
    std::array<int, 2> tADR;
    VariantInfo(float q, int dp1, int dp2, std::array<int, 2> tADR1) {
        qual = q;
        tbDP = dp1;
        tDP = dp2;
        for (int i = 0; i < 2; i++) {
            tADR[i] = tADR1[i];
        }
    };
    bool operator < (const VariantInfo & vi2) const {
        return 0;
    }
};

const auto MIN(const auto a, const auto b) { return ((a) < (b) ? (a) : (b)); }
const auto MAX(const auto a, const auto b) { return ((a) > (b) ? (a) : (b)); }
const auto UPDATE_MIN(auto & a, const auto b) { a = MIN(a, b); }
const auto UPDATE_MAX(auto & a, const auto b) { a = MAX(a, b); }

template <class T>
std::string
other_join(const T & container, std::string sep = std::string(",")) {
    std::string ret = "";
    for (const auto & e : container) {
        ret += std::to_string(e) + sep;
    }
    if (ret.size() > 0) { ret.pop_back(); }
    return ret;
}

int 
cHapSubstr_to_totDP(const std::string & cHapSubstr) {
    int ret = 0;
    int tot_n_amperstands = 0;
    for (const auto ch : cHapSubstr) {
        if ('&' == ch) {
            tot_n_amperstands++;
        }
    }
    std::string tmp;
    std::stringstream ss(cHapSubstr);
    int ntokens = 0;
    while(getline(ss, tmp, '&')) {
        if ((tot_n_amperstands - 1 == ntokens) || (tot_n_amperstands == ntokens)) {
            ret += atoi(tmp.c_str());
        }
        ntokens++;
    }
    return ret;
}

inline int varlen2reflen(int varlen) {
    if (varlen > 0) { return varlen * 2; }
    else { return -varlen; }
}

std::vector<std::vector<std::tuple<int, std::string, std::string, VariantInfo>>> 
vecof_pos_ref_alt_tup_split(const std::vector<std::tuple<int, std::string, std::string, int, int, VariantInfo>> & vecof_pos_ref_alt_begpos_endpos_tup,
        int defaultCB, int defaultCO, int defaultCE, bool enable_short_tandem_repeat_adjust) {
    std::vector<std::vector<std::tuple<int, std::string, std::string, VariantInfo>>> vecof_vecof_pos_ref_alt_tup;
    int prev_pos = INT32_MIN;
    int prev_varlen = 0;
    int nextof_prev = 0;
    int prevof_curr = 0;
    int prev_link_n_bases = 0;
    int delim_pos = 0;
    for (const auto & pos_ref_alt_begpos_endpos_tuple : vecof_pos_ref_alt_begpos_endpos_tup) {
        int varlen = ((int)(std::get<1>(pos_ref_alt_begpos_endpos_tuple).size()) - (int)(std::get<2>(pos_ref_alt_begpos_endpos_tuple).size()));
        const VariantInfo & variantinfo = std::get<5>(pos_ref_alt_begpos_endpos_tuple);
        
        int link_n_bases = 0;
        if (0 == varlen && 0 == prev_varlen) {
            link_n_bases = defaultCB;
        } else {
            link_n_bases = defaultCO + (MAX(varlen2reflen(varlen), abs(prev_varlen)) * defaultCE);
        }
        UPDATE_MAX(delim_pos, prev_pos + MAX(link_n_bases, prev_link_n_bases));
        int curr_pos = std::get<0>(pos_ref_alt_begpos_endpos_tuple);
        if (curr_pos >= delim_pos && ((!enable_short_tandem_repeat_adjust) || nextof_prev <= prevof_curr)) {
            vecof_vecof_pos_ref_alt_tup.push_back(std::vector<std::tuple<int, std::string, std::string, VariantInfo>>());
        }
        prev_link_n_bases = link_n_bases;
        prev_pos = std::get<0>(pos_ref_alt_begpos_endpos_tuple) + (int)MAX(std::get<1>(pos_ref_alt_begpos_endpos_tuple).size(), std::get<2>(pos_ref_alt_begpos_endpos_tuple).size());
        prev_varlen = varlen;
        nextof_prev = std::get<4>(pos_ref_alt_begpos_endpos_tuple);
        vecof_vecof_pos_ref_alt_tup.back().push_back(
            std::make_tuple(
                std::get<0>(pos_ref_alt_begpos_endpos_tuple),
                std::get<1>(pos_ref_alt_begpos_endpos_tuple),
                std::get<2>(pos_ref_alt_begpos_endpos_tuple),
                variantinfo
            )
        );
    }
    return vecof_vecof_pos_ref_alt_tup;
}

const std::vector<std::string> 
cHapString_to_cHapSubstrs(const std::string & cHapString) {
    std::vector<std::string> cHap_substrs;
    std::string cHap_substr = "";
    int parenlevel = 0;
    for (const auto cHap_char : cHapString) {
        if ('(' == cHap_char) {
            parenlevel++; 
        }
        if (')' == cHap_char) {
            parenlevel--; 
        }
        cHap_substr.push_back(cHap_char);
        if (0 == parenlevel) {
            cHap_substrs.push_back(cHap_substr);
            cHap_substr = "";
        }
    }
    return cHap_substrs;
}

std::map<std::string, int> build_tname2tid_from_faidx(const faidx_t *faidx) {
    std::map<std::string, int> ret;
    for (int i = 0; i < faidx_nseq(faidx); i++) {
        ret.insert(std::make_pair(faidx_iseq(faidx, i), i));
    }
    return ret;
}

const double DEFAULT_C = 0.75 + 1e-6;
const double DEFAULT_C2 = 0.50 + 1e-6;
const double DEFAULT_C3 = 0.75 + 1e-6;
const int DEFAULT_D = 3;
const double DEFAULT_F = 0.1 + 1e-6;
const char *DEFAULT_H = "cHap"; // bHap, cHap, c2Hap
const char *DEFAULT_A = "AD";   // bAD,  AD,   c2AD
const char *DEFAULT_M = "w";
const int DEFAULT_CB = 4; // SNV to SNV
const int DEFAULT_CO = 6; // (SNV to InDel gap-open) and (InDel to InDel gap-open)
const int DEFAULT_CE = 1; // (SNV to InDel gap-ext ) and (InDel to InDel gap-ext )
// const int DEFAULT_CR = 10;   // (SNV to InDel repeat-decrease) and (InDel to InDel repeat-decrease)
const double  POWLAW_EXPONENT = 3.0;

void help(int argc, char **argv) {
    fprintf(stderr, "Program %s version %s.%s ( %s )\n", argv[0], VERSION3, COMMIT_VERSION, COMMIT_DIFF_SH);
    fprintf(stderr, "  This program combines simple variants into delins variants and prints the result VCF to stdout. \n");
    
    fprintf(stderr, "Usage: %s <REFERENCE-FASTA> <UVC-VCF-GZ> \n", argv[0]);
    fprintf(stderr, "Optional parameters:\n");
    
    fprintf(stderr, " -2 the fraction of simple-variant depth used to construct the delins variant with the highest depth, above which the simple variant is discarded in -D if -3 is also satisfied [default to %f].\n", DEFAULT_C2);
    fprintf(stderr, " -3 min ratio of min-depth to max-depth of the simple variants, above which the simple variants are discarded in -D if -2 is also satisfied [default to %f].\n", DEFAULT_C3);

    fprintf(stderr, " -c the fraction of simple-variant depth used to construct the delins variant with the highest depth, above which the simple variant is discarded in -D [default to %f].\n", DEFAULT_C);
    fprintf(stderr, " -d minimum allele depth of the linked variants [default to %d].\n", DEFAULT_D);
    fprintf(stderr, " -f minimum fraction of the linked variants [default to %f].\n", DEFAULT_F);
    fprintf(stderr, " -p the power-law exponent for computing tumor haplotype variant qualities [default to %f].\n", POWLAW_EXPONENT);    
    
    fprintf(stderr, " -A FORMAT tag in the UVC-VCF-GZ file indicating allele depths. "
            "The 2 values bAD/AD/c2AD at -A correspond to bHap/cHap/c2Hap at -H. [default to %s].\n", DEFAULT_A);
    fprintf(stderr, " -B from BWA: maximum number of bases between SNV and SNV to be considered as linked [default to %d].\n", DEFAULT_CB);
    fprintf(stderr, " -C the VCF file containing simple variants discarded by constructing delins variants [default to NULL, generating no output].\n");
    fprintf(stderr, " -D the VCF file containing simple variants kept after constructing delins variants [default to NULL, generating no output].\n");
    fprintf(stderr, " -E from BWA: gap extension for the maximum number of bases between InDel and SNV/InDel to be considered as linked [default to %d].\n", DEFAULT_CE);
    fprintf(stderr, " -H FORMAT tag in the UVC-VCF-GZ file used to contain the haplotype information [default to %s].\n", DEFAULT_H);
    fprintf(stderr, " -M mode for the output VCF file containing simple variants that are entirely parts of some delins variant [default to %s].\n", DEFAULT_M);
    fprintf(stderr, " -O from BWA: gap opening for the maximum number of bases between InDel and SNV/InDel to be considered as linked [default to %d].\n", DEFAULT_CO);
    fprintf(stderr, " -T the bed file that overrides the -d, -f, -B -O, and -E parameters in the defined regions [default to None].\n");
    
    fprintf(stderr, " -I boolean flag indicating if the generation of delins variants is not required for the elimination of simple variants. [default to false].\n");
    fprintf(stderr, " -L boolean flag indicating if left-trimming of bases occurring in both REF and ALT should be disabled [default to false].\n");
    fprintf(stderr, " -R boolean flag indicating if right-trimming of bases occurring in both REF and ALT should be disabled [default to false].\n");
    fprintf(stderr, " -S boolean flag indicating if short-tandem-repeats (STRs) should be considered in the merging of simple variants [default to false].\n");
    
    fprintf(stderr, " -h print this usage help then exit with zero.\n");
    fprintf(stderr, " -v print version number then exit with zero.\n");
    
    fprintf(stderr, "\nNote: the output VCF (which is printed to stdout) has the following INFO tags\n");
    for (const auto bcf_info : BCF_INFO_LIST) {
        fprintf(stderr, "%s\n", bcf_info.to_header_string().c_str());
    }
}

int main(int argc, char **argv) {
    
    std::clock_t c_beg = std::clock(); 
    auto t_beg = std::chrono::high_resolution_clock::now();
    
    time_t rawtime;
    time(&rawtime);
    char timestring[80];
    strftime(timestring, 80, "%F %T", localtime(&rawtime));

    char *fastaref = NULL;
    char *uvcvcf = NULL;
    
    char *bedfile = NULL;
    double powlaw_exponent = POWLAW_EXPONENT;
    double defaultC = DEFAULT_C;
    double defaultC2 = DEFAULT_C2;
    double defaultC3 = DEFAULT_C3;
    int linkdepth1 = DEFAULT_D;
    double linkfrac1 = DEFAULT_F;
    int defaultCB1 = DEFAULT_CB;
    int defaultCO1 = DEFAULT_CO;
    int defaultCE1 = DEFAULT_CE;
    bool disable_is_part_of_delinsvar_3tups_check = false;
    bool enable_short_tandem_repeat_adjust = false;
    bool disable_left_trim = false;
    bool disable_right_trim = false;
    const char *defaultH1 = DEFAULT_H;
    const char *defaultAD = DEFAULT_A; // bAD, AD, c2AD for bHap, cHap, and c2Hap
    const char *simple_outvcfname = NULL;
    const char *non_delins_outvcfname = NULL;
    const char *defaultMode = DEFAULT_M;
    int opt = -1;
    while ((opt = getopt(argc, argv, "2:3:c:d:f:p:A:B:C:D:E:M:H:O:T:ILRShv")) != -1) {
        switch (opt) {
            case 'c': defaultC = atof(optarg); break; // delins2simple_var_frac_above_which_discard_simple 
            case '2': defaultC2 = atof(optarg); break;
            case '3': defaultC3 = atof(optarg); break;
            
            case 'd': linkdepth1 = atoi(optarg); break;
            case 'f': linkfrac1 = atof(optarg); break;
            case 'p': powlaw_exponent = atof(optarg); break;
            
            case 'A': defaultAD = optarg; ; break;            
            case 'B': defaultCB1 = atoi(optarg); break;
            case 'C': simple_outvcfname = optarg; break;
            case 'D': non_delins_outvcfname = optarg; break;
            case 'E': defaultCE1 = atoi(optarg); break;
            case 'M': defaultMode = optarg; break;
            case 'H': defaultH1 = optarg; break;
            case 'O': defaultCO1 = atoi(optarg); break;
            case 'T': bedfile = optarg; break;
            
            case 'I': disable_is_part_of_delinsvar_3tups_check = true; break;
            case 'L': disable_left_trim = true; break;
            case 'R': disable_right_trim = true; break;
            case 'S': enable_short_tandem_repeat_adjust = true; break;
            
            case 'h': help(argc, argv); exit(0);
            case 'v': printf("UVC-delins %s.%s ( %s )\n", VERSION3, COMMIT_VERSION, COMMIT_DIFF_SH); exit(0);
            default: help(argc, argv); exit(1);
        }
    }
    for (int posidx = 0; optind < argc; optind++, posidx++) {
        if (0 == posidx) { fastaref = argv[optind]; }
        else if (1 == posidx) { uvcvcf = argv[optind]; }
    }
    if (NULL == fastaref || NULL == uvcvcf) {
        help(argc, argv); exit(2);
    }
    
    faidx_t *faidx = fai_load(fastaref);
    htsFile *fp = vcf_open(uvcvcf, "r");
    bcf_hdr_t *bcf_hdr = vcf_hdr_read(fp);
    
    for (const auto bcf_info : BCF_INFO_LIST) {
        bcf_hdr_append(bcf_hdr, bcf_info.to_header_string().c_str());
    }
    /*
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=delinsHap,Number=A,Type=String,Description=\"For each delins ALT allele, the value of this tag is the key of the other tag "
            "used by the upstream SNV/InDel caller to represent linked variants that are phased into one single haplotype candidate. "
            "In [x]Hap, [x] can be b, c, or c2 if UVC was used upstream. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diPRA,Number=A,Type=String,Description=\"Tumor position_REF_ALT of each delins ALT allele, with the three VCF fields separated by underscore. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diADA,Number=A,BCF_TYPE_INTEGER,Description=\"Depth of each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the depth exactly supporting its corresponding delins haplotype. "
            "For example, at read level, only reads supporting the exact delins haplotype are counted. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diDPm,Number=A,BCF_TYPE_INTEGER,Description=\"Total depth of all alleles located at each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the minimum depth among the SNV and/or InDel alleles that the delins variant is composed of. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diDPM,Number=A,BCF_TYPE_INTEGER,Description=\"Total depth of all alleles located at each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the maximum depth among the SNV and/or InDel alleles that the delins variant is composed of. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diADm,BCF_NUM_R,BCF_TYPE_INTEGER,Description=\"Depth of each delins (MVV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the minimum delins-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diADM,BCF_NUM_R,BCF_TYPE_INTEGER,Description=\"Depth of each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the maximum delins-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. \">"); // this should rarely be used
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diRDm,BCF_NUM_R,BCF_TYPE_INTEGER,Description=\"Depth of the reference allele corresponding to each delins (MVV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the minimum REF-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diRDM,BCF_NUM_R,BCF_TYPE_INTEGER,Description=\"Depth of the reference allele corresponding to each delins (MNV and/or complex InDel) ALT allele. "
            "Here, each depth is defined as the maximum REF-allele depth among the individual SNVs and/or InDels that constitute the delins variant. "
            "For example, at read level, reads supporting an individual SNV and/or InDel but not supporting the delins haplotype are still counted. \">"); // this should rarely be used
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diAD2F,Number=A,BCF_TYPE_INTEGER,Description=\"Percentage (between 0 and 100) of reads that support the delins variant "
            "among the individual SNV and/or InDel variants that constitute the delins variant. \">");
    
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diCVQ,Number=A,BCF_TYPE_INTEGER,Description=\"Combined (by combining SNVs and InDels) variant Quality of each delins ALT variant. \">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=diHVQ,Number=A,BCF_TYPE_INTEGER,Description=\"Haplotyped (non-SNV and non-InDel small) variant Quality of each delins ALT variant. \">");
    */
    
    bcf_hdr_append(bcf_hdr, (std::string("##delinsVariantDate=") + timestring).c_str());
    bcf_hdr_append(bcf_hdr, "##delinsVariantCallerVersion=" COMMIT_VERSION "(" COMMIT_DIFF_SH ")");
    std::string vcfcmd = "##delinsVariantCallerCommand=";
    for (int i = 0; i < argc; i++) {
        vcfcmd += std::string("") + std::string(argv[i]) + "  ";
    }
    bcf_hdr_append(bcf_hdr, vcfcmd.c_str());

    bcf_hdr_t *bcf_hdr2 = bcf_hdr_dup(bcf_hdr);
    // int set_samples_ret = bcf_hdr_set_samples(bcf_hdr, bcf_hdr->samples[bcf_hdr->nsamples_ori - 1], false);
    int set_samples_ret1 = bcf_hdr_set_samples(bcf_hdr2, NULL, false);
    assert(0 == set_samples_ret1);
    kstring_t hdr_str = {0, 0, NULL};
    bcf_hdr_format(bcf_hdr2, true, &hdr_str);
    std::cout << hdr_str.s;
    bcf_hdr_destroy(bcf_hdr2);
    // int set_samples_ret2 = bcf_hdr_set_samples(bcf_hdr, "-", false);
    // assert(0 == set_samples_ret2);
    htsFile *simple_outvcf = ((NULL == simple_outvcfname) ? NULL : vcf_open(simple_outvcfname, (defaultMode)));
    htsFile *non_delins_outvcf = ((NULL == non_delins_outvcfname) ? NULL : vcf_open(non_delins_outvcfname, (defaultMode)));
    
    if (NULL != simple_outvcf) {
        int vcf_hdr_write_ret = vcf_hdr_write(simple_outvcf, bcf_hdr);
        assert (vcf_hdr_write_ret >= 0);
    }
    if (NULL != non_delins_outvcf) {
        int vcf_hdr_write_ret = vcf_hdr_write(non_delins_outvcf, bcf_hdr);
        assert (vcf_hdr_write_ret >= 0);
    }
    int vcf_nseqs = -1;
    const char **seqnames = bcf_hdr_seqnames(bcf_hdr, &vcf_nseqs);
    
    if (faidx_nseq(faidx) != vcf_nseqs) {
        fprintf(stderr, "FATAL_ERROR : the FASTA file %s contains %d reference sequences whereas the VCF file %s contains %d reference sequences. "
                "Please ensure that both FASTA and VCF files contain the same number of reference sequences. \n", fastaref, faidx_nseq(faidx), uvcvcf, vcf_nseqs);
        exit (1);
    }
    
    for (int i = 0; i < faidx_nseq(faidx); i++) {
        const char *seqname = faidx_iseq(faidx, i);
        if (strcmp(seqname, seqnames[i])) {
            fprintf(stderr, "FATAL_ERROR : the FASTA file %s contains the reference sequence with ID %s "
                    "whereas the VCF file %s contains the reference sequence with ID %s when comparing the %d-th reference sequence. "
                    "Please ensure that both FASTA and VCF files contain exactly the same list of reference sequences. \n", fastaref, seqname, uvcvcf, seqnames[i], i);
            exit (2);
        }
    }
    
    std::ifstream bedstream;
    std::map<std::string, int> tname2tid;
    if (NULL != bedfile) { 
        bedstream.open(bedfile); 
        tname2tid = build_tname2tid_from_faidx(faidx);
    }
    int bedtid = -1;
    int bedbeg = -1;
    int bedend = -1;

    for (int i = 0; i < vcf_nseqs; i++) {
        const char *tname = faidx_iseq(faidx, i);
        int regionlen;
        char *fetchedseq = fai_fetch(faidx, tname, &regionlen);
        const std::string refstring = fetchedseq;
        free(fetchedseq);
        
        bcf_srs_t *const sr = bcf_sr_init();
        if (NULL == sr) {
            std::cerr << "Failed to initialize bcf sr\n";
            exit(-6);
        }

        bcf_sr_set_regions(sr, tname, false);
        int sr_set_opt_retval = bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
        if (sr_set_opt_retval < 0) {
            exit(-8);
        }
        int sr_add_reader_retval = bcf_sr_add_reader(sr, uvcvcf);
        if (sr_add_reader_retval != 1) {
            exit(-9);
        }

        int valsize = 0;
        int ndst_val = 0;
        int string_valsize = 0;
        int string_ndst_val = 0;
        char **bcfstring = NULL;
        int32_t *bcfints = NULL;
        
        std::map<std::string, std::vector<std::tuple<int, std::string, std::string, int, int, VariantInfo>>> map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup;
        
        int nsamples = bcf_hdr_nsamples(bcf_hdr); 
        int sampleidx = nsamples - 1; // last sample
        
        while (bcf_sr_next_line(sr)) {
            bcf1_t *line = bcf_sr_get_line(sr, 0);
            bcf_unpack(line, BCF_UN_ALL);
            ndst_val = 0;
            
            fprintf(stderr, "Processing the line: tid = %d pos = %ld\n", line->rid, line->pos);
            valsize = bcf_get_info_int32(bcf_hdr, line, "R3X2", &bcfints, &ndst_val);
            if (valsize < 6) { continue; }
            const int posleft = bcfints[0];
            const int posright = bcfints[3] + (bcfints[4] * bcfints[5]);
            
            if (bcfstring) {
                free(bcfstring[0]);
                free(bcfstring);
                bcfstring = NULL;
                string_ndst_val = 0;
            }
            string_valsize = bcf_get_format_string(bcf_hdr, line, defaultH1, &bcfstring, &string_ndst_val);
            assert (1 <= string_valsize || !fprintf(stderr, "The size of cHap is %d instead of at least 1!\n", string_valsize));
            valsize = bcf_get_info_int32(bcf_hdr, line, "tbDP", &bcfints, &ndst_val);
            int tbDP = bcfints[0];
            valsize = bcf_get_info_int32(bcf_hdr, line, "tDP", &bcfints, &ndst_val);
            int tDP = bcfints[0];
            valsize = bcf_get_info_int32(bcf_hdr, line, "tADR", &bcfints, &ndst_val);
            std::array<int, 2> tADR = std::array<int, 2>({bcfints[0], bcfints[1]});
            
            const auto pos_ref_alt_begpos_endpos_tup = std::make_tuple(line->pos, std::string(line->d.allele[0]), std::string(line->d.allele[1]), posleft, posright,
                    VariantInfo(line->qual, tbDP, tDP, tADR));
            // This assertion may fail (for example, SRR12100766 - hs37d5 21645705 . T C) due to corner-case variant-call and is therefore disabled. 
            // assert (tDP > 0 || !fprintf(stderr, "%d > 0 failed for rid - %d pos - %ld ref - %s alt - %s!\n", tDP, line->rid, line->pos, line->d.allele[0], line->d.allele[1]));
            for (int j = sampleidx; j < nsamples; j++) {
                std::vector<std::string> cHap_substrs = cHapString_to_cHapSubstrs(bcfstring[j]);
                for (const std::string & cHap_string : cHap_substrs) {
                    map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup.insert(std::make_pair(cHap_string, std::vector<std::tuple<int, std::string, std::string, int, int, VariantInfo>>()));
                    map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup[cHap_string].push_back(pos_ref_alt_begpos_endpos_tup);
                }
            }
        }
        
        bcf_sr_seek(sr, tname, 0);
        std::cerr << "Will finish processing tname " << tname << "\n";
        std::set<std::tuple<int, std::string, std::string>> delinsvar_3tups;
        while (bcf_sr_next_line(sr)) {
           
            int linkdepth = linkdepth1;
            int linkfrac = linkfrac1;
            int defaultCB = defaultCB1;
            int defaultCO = defaultCO1;
            int defaultCE = defaultCE1;
            double delins2simple_var_frac_above_which_discard_simple = defaultC;
            
            bcf1_t *line = bcf_sr_get_line(sr, 0);
            bcf_unpack(line, BCF_UN_ALL); 
            
            if (bedfile != NULL && bedstream.good()) {
                int vcftid = line->rid;
                int vcfbeg = line->pos;
                while (bedtid < vcftid || (bedtid == vcftid && bedbeg < vcfbeg)) {
                    std::string line;
                    getline(bedstream, line);
                    if (line.empty()) { break; }
                    std::istringstream linestream(line);
                    std::string bedchrom;
                    linestream >> bedchrom;
                    linestream >> bedbeg;
                    linestream >> bedend;
                    if (tname2tid.find(bedchrom) == tname2tid.end()) {
                        fprintf(stderr, "The bed file %s with line %s has invalid tname %s so this line is skipped!\n", bedfile, line.c_str(), bedchrom.c_str());
                        bedtid = -1;
                    } else {
                        bedtid = tname2tid[bedchrom];
                        std::string token;
                        while (linestream.good()) {
                            linestream >> token;
                            if (!token.compare("-c")) {
                                linestream >> token;
                                delins2simple_var_frac_above_which_discard_simple = atof(token.c_str());
                            }
                            if (!token.compare("-d")) {
                                linestream >> token;
                                linkdepth = atoi(token.c_str());
                            }
                            if (!token.compare("-f")) {
                                linestream >> token;
                                linkfrac = atof(token.c_str());
                            }
                            if (!token.compare("-B")) {
                                linestream >> token;
                                defaultCB = atoi(token.c_str());
                            }
                            if (!token.compare("-O")) {
                                linestream >> token;
                                defaultCO = atoi(token.c_str());
                            }
                            if (!token.compare("-E")) {
                                linestream >> token;
                                defaultCE = atoi(token.c_str());
                            }
                        }
                    }
                }
            }
            
            ndst_val = 0;
            valsize = bcf_get_format_string(bcf_hdr, line, defaultH1, &bcfstring, &ndst_val);
            if (valsize <= 0) { continue; }
            ndst_val = 0;
            valsize = bcf_get_format_int32(bcf_hdr, line, defaultAD, &bcfints, &ndst_val);
            if (valsize <= 0) { continue; } 
            int vcflineAD = bcfints[ndst_val - 1];
            const auto pos_ref_alt_tup_from_vcfline = std::make_tuple(line->pos, std::string(line->d.allele[0]), std::string(line->d.allele[1]));
            
            for (int j = sampleidx; j < nsamples; j++) {
                bool is_part_of_delinsvar_3tups = false;
                int max_delinsvarDP = 0;
                std::vector<std::string> cHap_substrs = cHapString_to_cHapSubstrs(bcfstring[j]);
                int tADRmin_overall = 0;
                int tADRmax_overall = 0;
                for (std::string cHap_string : cHap_substrs) {
                    int delinsvarDP = cHapSubstr_to_totDP(cHap_string);
                    UPDATE_MAX(max_delinsvarDP, delinsvarDP);
                    if (delinsvarDP < linkdepth || delinsvarDP < linkfrac * vcflineAD) {
                        // fprintf(stderr, "Skipping the variant %s %d because it has low DP\n", tname, line->pos);
                        continue; // this variant has low DP
                    }
                    auto vecof_pos_ref_alt_begpos_endpos_tup = map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup.find(cHap_string)->second;
                    if (vecof_pos_ref_alt_begpos_endpos_tup.size() == 0) {
                        // fprintf(stderr, "Skipping the variant %s %d because it has no variants!\n", tname, line->pos);
                        continue; // no variant is found
                    }
                    
                    std::sort(vecof_pos_ref_alt_begpos_endpos_tup.begin(), vecof_pos_ref_alt_begpos_endpos_tup.end());
                    
                    const auto vecof_vecof_pos_ref_alt_tup = vecof_pos_ref_alt_tup_split(vecof_pos_ref_alt_begpos_endpos_tup,
                        defaultCB, defaultCO, defaultCE, enable_short_tandem_repeat_adjust);
                    
                        for (const auto & vecof_pos_ref_alt_tup1 : vecof_vecof_pos_ref_alt_tup) {
                        
                        if (vecof_pos_ref_alt_tup1.size() <= 1) { 
                            // fprintf(stderr, "Skipping the variant %s %d because it is not delins\n", tname, line->pos);
                            continue; // this variant is not delins
                        }
                        double errmodel_delinsvarfrac_phred = 0.0;
                        int delinsvar_begpos = INT32_MAX;
                        int delinsvar_endpos = 0;
                        float cv_qual = FLT_MAX;
                        int tDPmin = INT_MAX;
                        int tDPmax = 0;
                        std::array<int, 2> tADRmin = {INT_MAX, INT_MAX};
                        std::array<int, 2> tADRmax = {0};
                        for (auto pos_ref_alt_tup : vecof_pos_ref_alt_tup1) {
                            int pos = std::get<0>(pos_ref_alt_tup);
                            std::string ref = std::get<1>(pos_ref_alt_tup);
                            std::string alt = std::get<2>(pos_ref_alt_tup);
                            int endpos = pos + (int)MAX(alt.size(), ref.size());
                            UPDATE_MIN(delinsvar_begpos, pos);
                            UPDATE_MAX(delinsvar_endpos, endpos);
                            float qual = std::get<3>(pos_ref_alt_tup).qual;
                            int tDP = std::get<3>(pos_ref_alt_tup).tDP;
                            const auto &tADR = std::get<3>(pos_ref_alt_tup).tADR;
                            
                            UPDATE_MIN(tDPmin, tDP);
                            UPDATE_MAX(tDPmax, tDP);
                            for (int j = 0; j < 2; j++) {
                                UPDATE_MIN(tADRmin[j], tADR[j]);
                                UPDATE_MAX(tADRmax[j], tADR[j]);
                                assert (tADR[1] <= tDP);
                                assert (delinsvarDP <= tADR[1]);
                            }
                            errmodel_delinsvarfrac_phred += 10.0 / log(10.0) * log((double)(tADR[1] - delinsvarDP + 0.5) / (double)(tDP + 1.0));
                            cv_qual = MIN(qual, cv_qual);
                        }
                        double delinsvarfrac = (double)(delinsvarDP + 0.5) / (double)(tDPmax + 1.0);
                        double delinsvarfrac_ratio_phred = 10.0 / log(10.0) * log(delinsvarfrac) - errmodel_delinsvarfrac_phred;
                        
                        std::string delins_ref = refstring.substr(delinsvar_begpos, delinsvar_endpos - delinsvar_begpos);
                        std::vector<std::string> delins_alt_;
                        for (const auto base : delins_ref) {
                            delins_alt_.push_back(std::string(&base, 1));
                        }
                        for (auto pos_ref_alt_tup : vecof_pos_ref_alt_tup1) {
                            int pos = std::get<0>(pos_ref_alt_tup);
                            std::string ref = std::get<1>(pos_ref_alt_tup);
                            std::string alt = std::get<2>(pos_ref_alt_tup);
                            for (int k = pos; k < pos + (int)ref.size(); k++) {
                                delins_alt_[k - delinsvar_begpos] = "";
                            }
                            delins_alt_[pos - delinsvar_begpos] = alt;
                        }
                        std::string delins_alt = "";
                        for (const auto subalt : delins_alt_) {
                            delins_alt.append(subalt);
                        }
                        int vcfline_pos = (std::get<0>(pos_ref_alt_tup_from_vcfline));
                        if (delinsvar_begpos > vcfline_pos || vcfline_pos >= (delinsvar_begpos + (int)MAX(delins_ref.size(), delins_alt.size()))) { 
                            /*
                            fprintf(stderr, "Skipping the variant %s %d %s %s because it is outside the ROI\n", 
                                tname, vcfline_pos, 
                                std::get<1>(pos_ref_alt_tup_from_vcfline).c_str(), 
                                std::get<2>(pos_ref_alt_tup_from_vcfline).c_str());
                            */
                            continue; // this var is outside the delins var region
                        }
                        const auto delinsvar_3tup = std::make_tuple(delinsvar_begpos, delins_ref, delins_alt);
                        
                        if (delinsvar_3tups.find(delinsvar_3tup) != delinsvar_3tups.end()) { 
                            fprintf(stderr, "Skipping the variant %s %d %s %s because it is already visited\n", 
                                tname, vcfline_pos, 
                                std::get<1>(pos_ref_alt_tup_from_vcfline).c_str(), 
                                std::get<2>(pos_ref_alt_tup_from_vcfline).c_str());
                            is_part_of_delinsvar_3tups = true;
                            tADRmin_overall += tADRmin[1];
                            tADRmax_overall += tADRmax[1];
                            continue; // this var has already been printed
                        }
                        int cv_begpos = delinsvar_begpos;
                        std::string cv_ref = delins_ref;
                        std::string cv_alt = delins_alt;
                        if (!disable_left_trim) {
                            size_t begpos_inc = 0;
                            while ((begpos_inc + 1 < cv_ref.size())
                                    && (begpos_inc + 1 < cv_alt.size())
                                    && (cv_ref[begpos_inc] == (cv_alt[begpos_inc]))) {
                                begpos_inc++;
                            }
                            cv_begpos += begpos_inc;
                            cv_ref = cv_ref.substr(begpos_inc);
                            cv_alt = cv_alt.substr(begpos_inc);
                        }
                        if (!disable_right_trim) {
                            size_t endpos_dec = 0;
                            while (cv_ref.size() > (1 + endpos_dec) && cv_alt.size() > (1 + endpos_dec)
                                    && cv_ref[cv_ref.size() - (1 + endpos_dec)] == cv_alt[cv_alt.size() - (1 + endpos_dec)]) {
                                endpos_dec++;
                            }
                            cv_ref = cv_ref.substr(0, cv_ref.size() - endpos_dec);
                            cv_alt = cv_alt.substr(0, cv_alt.size() - endpos_dec);
                        }
                        double hv_qual = (powlaw_exponent * delinsvarfrac_ratio_phred); // c and h : combined and haplotyped combined
                        std::cout 
                            << tname 
                            << "\t" << (cv_begpos + 1) 
                            << "\t.\t" << cv_ref 
                            << "\t" << cv_alt 
                            << "\t" << std::to_string(MIN(cv_qual, hv_qual)) 
                            << "\t.\t" // FILTER
                            <<        BCF_INFO_LIST[delinsHap].ID << "=" << cHap_string 
                            << ";" << BCF_INFO_LIST[diPRA].ID     << "="
                                << (std::get<0>(pos_ref_alt_tup_from_vcfline)) << "_"
                                << (std::get<1>(pos_ref_alt_tup_from_vcfline)) << "_"
                                << (std::get<2>(pos_ref_alt_tup_from_vcfline))
                            << ";" << BCF_INFO_LIST[diADA].ID    << "=" << delinsvarDP
                            << ";" << BCF_INFO_LIST[diDPm].ID    << "=" << tDPmin
                            << ";" << BCF_INFO_LIST[diDPM].ID    << "=" << tDPmax
                            << ";" << BCF_INFO_LIST[diADm].ID    << "=" << tADRmin[1]
                            << ";" << BCF_INFO_LIST[diADM].ID    << "=" << tADRmax[1] 
                            << ";" << BCF_INFO_LIST[diRDm].ID    << "=" << tADRmin[0] 
                            << ";" << BCF_INFO_LIST[diRDM].ID    << "=" << tADRmax[0]
                            << ";" << BCF_INFO_LIST[diAD2F].ID   << "=" << (tADRmin[1] * 100 / MAX(1, tADRmax[1]))
                            << ";" << BCF_INFO_LIST[diHVQ].ID    << "=" << hv_qual
                            << ";" << BCF_INFO_LIST[diCVQ].ID    << "=" << cv_qual
                            << "\n";
                        for (auto it = delinsvar_3tups.begin(); it != delinsvar_3tups.end(); ) {
                            int endpos = std::get<0>(*it) + (int)MAX(std::get<1>(*it).size(), std::get<2>(*it).size());
                            if (endpos < vcfline_pos) {
                                it = delinsvar_3tups.erase(it);
                            } else {
                                it++;
                            }
                        }
                        tADRmin_overall += tADRmin[1];
                        tADRmax_overall += tADRmax[1];
                        delinsvar_3tups.insert(delinsvar_3tup);
                        if (delinsvar_3tups.find(delinsvar_3tup) != delinsvar_3tups.end()) {
                            is_part_of_delinsvar_3tups = true;
                        }
                    }
                }
                const bool are_all_vars_delins1 = (max_delinsvarDP > vcflineAD * delins2simple_var_frac_above_which_discard_simple);
                const bool are_all_vars_delins2 = ((max_delinsvarDP > vcflineAD * defaultC2) && (tADRmin_overall > tADRmax_overall * defaultC3));
                const bool are_all_vars_delins = (are_all_vars_delins1 || are_all_vars_delins2);
                if ((simple_outvcf != NULL)
                        && (is_part_of_delinsvar_3tups || disable_is_part_of_delinsvar_3tups_check)
                        && are_all_vars_delins) { 
                    // is mostly part of a delins variant
                    int vcf_write_ret = vcf_write(simple_outvcf, bcf_hdr, line);
                    assert(vcf_write_ret >= 0);
                }
                if ((non_delins_outvcf != NULL) && (
                        (!is_part_of_delinsvar_3tups)
                        || (!are_all_vars_delins))) {
                    // is not mostly part of a delins variant
                    int vcf_write_ret = vcf_write(non_delins_outvcf, bcf_hdr, line);
                    assert(vcf_write_ret >= 0 || !fprintf(stderr, "The VCF record at tid-pos-ref-alt %d %ld %s %s\n", line->rid, line->pos, line->d.allele[0], line->d.allele[1]));
                }
            }
        }
        bcf_sr_destroy(sr);
    }
    if (NULL != bedfile) { 
        bedstream.close(); 
    }
    if (NULL != simple_outvcf) {
        int vcf_close_ret = vcf_close(simple_outvcf);
        assert (vcf_close_ret == 0);
    }
    if (NULL != non_delins_outvcf) {
        int vcf_close_ret = vcf_close(non_delins_outvcf);
        assert (vcf_close_ret == 0);
    }
    bcf_hdr_destroy(bcf_hdr);
    vcf_close(fp);
    fai_destroy(faidx);

    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
    
    std::cerr << std::fixed << std::setprecision(2) << "CPU time used: "
              << 1.0 * (c_end - c_beg) / CLOCKS_PER_SEC << " seconds\n"
              << "Wall clock time passed: "
              << std::chrono::duration<double>(t_end - t_beg).count()
              << " seconds\n";
}

