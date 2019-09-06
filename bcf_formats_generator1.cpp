// This is the framework for generating VCF/BCF header and its associated data on each line.
// Manual update of VCF/BCF header and its associated data on each line is extremely error-prone.
// Therefore, this framework is used.
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#define BCF_NUM_A (-1)
#define BCF_NUM_R (-2)
#define BCF_NUM_G (-3)
#define BCF_NUM_D (-4)

enum BCF_DATA_TYPE {
    BCF_STRING,
    BCF_INTEGER,
    BCF_FLOAT,
    BCF_SEP,
    NUM_BCF_DATA_TYPES
};

const char * CPP_DATA_STRING[] = {
    "std::string",
    "uint32_t",
    "float",
    "bool",
    NULL
};

const char * CPP_DATA_VALUES[] = {
    "\"\"",
    "0",
    "0",
    "false",
    NULL,
};

const char * BCF_DATA_STRING[] = {
    "String",
    "Integer",
    "Float",
    "String",
    NULL
};

const std::string vcf_number_to_header_str(int num) {
    if (0 < num) { return std::to_string(num); }
    else if (BCF_NUM_A == num) { return "A"; }
    else if (BCF_NUM_R == num) { return "R"; }
    else if (BCF_NUM_G == num) { return "G"; }
    else if (BCF_NUM_D == num) { return "."; } 
}

const std::vector<std::pair<std::string, std::string>> FILTER_VEC = {
    std::make_pair("noVar",         "Not a variant (for example, when REF and ALT are the same), but still included to get all statistics."),
    std::make_pair("upstreamDel",   "Deletion extended from another upstream deletion"),
    std::make_pair("cad3",          "Less than 3 clean deduppd reads"),
    std::make_pair("caf3",          "Less than 3/10000 allele fraction base on clean deduped reads"),
    std::make_pair("s50",           "Less than 50\% of samples have data"),
    std::make_pair("Q10",           "Quality below 10 and no other filters"),
    std::make_pair("Q20",           "Quality below 20 and no other filters"),
    std::make_pair("Q30",           "Quality below 30 and no other filters"),
    std::make_pair("Q40",           "Quality below 40 and no other filters"),
    std::make_pair("Q50",           "Quality below 50 and no other filters"),
    std::make_pair("Q60",           "Quality below 60 and no other filters"),
};

struct BcfFormatStruct {
    const char *id;
    int number;
    BCF_DATA_TYPE type;
    const char * description;
    BcfFormatStruct(const char *i, unsigned int n, BCF_DATA_TYPE t, const char *desc) {
        id = i;
        number = n;
        type = t;
        description = desc;
    }
};

const std::vector<BcfFormatStruct> FORMAT_VEC = {
    BcfFormatStruct("GT"       , 1, BCF_STRING,  "Genotype"),
    BcfFormatStruct("GQ"       , 1, BCF_INTEGER, "Genotype Quality"),
    BcfFormatStruct("HQ"       , 2, BCF_INTEGER, "Haplotype Quality"),
    BcfFormatStruct("DP"       , 1, BCF_INTEGER, "Fragment depth supporting any allele [consensus family, deduped]"),
    BcfFormatStruct("FA"       , 1, BCF_FLOAT,   "Frequency of the ALT allele [consensus family, deduped]"),
    BcfFormatStruct("FR"       , 1, BCF_FLOAT,   "Frequency of the REF allele [consensus family, deduped]"),
    BcfFormatStruct("DPHQ"     , 1, BCF_INTEGER, "Fragment depth supporting any allele if low-quality bases are ignored which means only high quality (HQ) bases are used [base read, duped]"),
    BcfFormatStruct("ADHQ"     , 1, BCF_INTEGER, "Fragment depth supporting the ALT allele if low-quality bases are ignored which means only high quality (HQ) bases are used [base read, duped]"),
    BcfFormatStruct("MQ"       , 1, BCF_FLOAT,   "Root mean square (RMS) mapping quality [base read, duped]"), 
    BcfFormatStruct("a_a"      , 1, BCF_SEP,     "Depth and REF/ALT allele frequency for base read and consensus family"),  
    BcfFormatStruct("bDP"      , 1, BCF_INTEGER, "Fragment depth supporting any allele [base read, duped]"),
    BcfFormatStruct("bFA"      , 1, BCF_FLOAT,   "Frequency of the ALT allele [base read, duped]"),
    BcfFormatStruct("bFR"      , 1, BCF_FLOAT,   "Frequency of the REF allele [base read, duped]"),
    BcfFormatStruct("bFO"      , 1, BCF_FLOAT,   "Frequency of all the other ALT alleles [base read, duped]"),
    BcfFormatStruct("cDP"      , 1, BCF_INTEGER, "Fragment depth supporting any allele [consensus family, deduped]"),
    BcfFormatStruct("cFA"      , 1, BCF_FLOAT,   "Frequency of the ALT allele [consensus family, deduped]"),
    BcfFormatStruct("cFR"      , 1, BCF_FLOAT,   "Frequency of the REF allele [consensus family, deduped]"),
    BcfFormatStruct("cFO"      , 1, BCF_FLOAT,   "Frequency of all the other ALT alleles [consensus family, deduped]"),
    BcfFormatStruct("a_b"      , 1, BCF_SEP,     "Consensus/variant allele type/quality"),  
    BcfFormatStruct("CType"    , 1, BCF_STRING,  "Type of consensus allele"),
    BcfFormatStruct("CAQ"      , 1, BCF_FLOAT,   "Consensus Allele Quality"),
    BcfFormatStruct("VType"    , 1, BCF_STRING,  "Type of variant allele"),
    BcfFormatStruct("VAQ"      , 1, BCF_FLOAT,   "Variant Allele Quality"),
    BcfFormatStruct("VAQ2"     , 1, BCF_FLOAT,   "Variant Allele Quality of the specific form(s) of InDel in ALT assuming other forms of InDels are noise"),
    BcfFormatStruct("a_c"      , 1, BCF_SEP,     "Average base quality for ALL/ALT (bAllBQ/bAltBQ) and duplication bias (aDB). All biases are percentages where 100% means no bias."), // global 
    BcfFormatStruct("bAllBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bAltBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bRefBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAllBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAltBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRefBQ"   , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("aDB"      , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("b_a"      , 1, BCF_SEP,     "Forward&reverse max-bias base distances to left/right end positions (T1PTL/T1PTR) and position bias (T1PBL/T1PBR) [base read, duped]"), 
    BcfFormatStruct("bPTL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bPTR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bPBL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bPBR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bMMT"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bMMB"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("b_b"      , 1, BCF_SEP,     "Forward&reverse sum of : base distances to left/right end positions (T1SDL/T1SDR) and non-normalized/normalized strand bias (T1SB1/T1SBR) [base read, duped]"),
    BcfFormatStruct("bSDL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bSDR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bSB1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bSBR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("b_c"      , 1, BCF_SEP,     "Forward&reverse depth (DP1), allele depth (AD1), allele depth after weak filter (AD2), quality threshold for weak filter (QT2) [base read, duped]"),
    BcfFormatStruct("bDP1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bAD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bRD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bAD2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bQT2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("b_d"      , 1, BCF_SEP,     "Forward&reverse allele depth after strong filter (AD3), quality threshold for strong filter (QT3), and variant quality after strong filter (VQ3) [duped]"),
    BcfFormatStruct("bAD3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bADB"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bQT3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bVQ3"     , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("b_e"      , 1, BCF_SEP,     "Forward&reverse allele root-mean-squre MQ (bMQ1) and (sum of squared MQ) / (sum of MQ) (bMQ2), where MQ is mapping quality, and the same for base quality BQ [duped]"),
    BcfFormatStruct("bMQ1"     , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("bMQ2"     , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("bBQ1"     , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("bBQ2"     , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("b_f"      , 1, BCF_SEP,     "Forward&reverse total depth (bDPLQ), ALT depth (bADLQ), and REF depth (bRDLQ) supported by low-quality (LQ) bases [duped]"),
    BcfFormatStruct("bDPLQ"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bADLQ"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bRDLQ"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("c_a"      , 1, BCF_SEP,     "Same as b_a but for single-strand families instead of reads [consensus family, deduped]"),
    BcfFormatStruct("cPTL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cPTR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cPBL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cPBR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cMMT"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cMMB"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("c_b"      , 1, BCF_SEP,     "Same as b_b but for single-strand families instead of reads [consensus family, deduped]"), 
    BcfFormatStruct("cSDL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cSDR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cSB1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cSBR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("c_c"      , 1, BCF_SEP,     "Same as b_c but for single-strand families instead of reads [consensus family, deduped]"), 
    BcfFormatStruct("cDP1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAD2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cQT2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("c_d"      , 1, BCF_SEP,     "Same as b_d but for single-strand families instead of reads [consensus family, deduped]"),
    BcfFormatStruct("cAD3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cADB"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cQT3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cVQ3"     , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("c_e"      , 1, BCF_SEP,     "Forward&reverse number of read supports generated by signal (cMajor) and noise (cMinor) inferred from heterogeneity of families, "
                                                 "and strand-specific VAQ (non-adjusted: cVAQ1; adjusted with BQ, MQ, and strand balance: cVAQ2) [consensus family, deduped]"),
    BcfFormatStruct("cMajor"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cMinor"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cVAQ1"    , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("cVAQ2"    , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("c_f"      , 1, BCF_SEP,     "Forward&reverse number of families without any filtering (total: cDPTT, allele: cADTT, reference: cRDTT, "
                                                 "allele singleton: cADT1, allele failing 80\% consensus : cADTN). "),
    BcfFormatStruct("cDPTT"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cADTT"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRDTT"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cADT1"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cADTN"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("d_a"      , 1, BCF_SEP,     "Same as b_a and c_a but for duplex families [duplex family, DSCS1-deduped]"),
    BcfFormatStruct("dDP1"     , 1, BCF_INTEGER, "see above"),
    BcfFormatStruct("dAD1"     , 1, BCF_INTEGER, "see above"),
    BcfFormatStruct("dVQ3"     , 1, BCF_INTEGER, "see above"), // AD2=QT2=NotApplicable QT3=60 (all passed)
    BcfFormatStruct("dAD3"     , 1, BCF_INTEGER, "see above"),
    // BcfFormatStruct("RVAQ"     , 1, BCF_FLOAT,   "Raw Variant Allele Quality"),
    BcfFormatStruct("e_gap"    , 1, BCF_SEP,     "Number of indels in forward&reverse strand (gapNum), max cAD diff (gapcADD), total cAD depth (gapcADT), indel sequences (gapSeq), duped read count of each gapSeq, dedup family count of each gapSeq"), 
    BcfFormatStruct("gapNum"   , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    BcfFormatStruct("gapcADD"  , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    BcfFormatStruct("gapcADT"  , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    BcfFormatStruct("gapSeq"   , BCF_NUM_D, BCF_STRING,  "see above"),
    BcfFormatStruct("gapbAD1"  , BCF_NUM_D, BCF_INTEGER, "see above"),
    BcfFormatStruct("gapcAD1"  , BCF_NUM_D, BCF_INTEGER, "see above"),
    BcfFormatStruct("bHap"     , 1, BCF_STRING,  "Duped forward&reverse linkage in the format of ((position&variantType)...depth)... "
                                                 "where ()... means more elements following the format in the preceding parenthesis. "),
    BcfFormatStruct("cHap"     , 1, BCF_STRING,  "Dedup forward&reverse linkage in the format of ((position&variantType)...depth)... "
                                                 "where ()... means more elements following the format in the preceding parenthesis. "),    
    BcfFormatStruct("note"     , 1, BCF_STRING,  "Additional note as comment for the given variant")
};

/*
 * for FORMAT, we have to keep in sync : 
 * C++ class variables per program, 
 * C++ class variable to stream methods per program, 
 * header lines per vcf file, 
 * tag-names (can be multiple sets of) per vcf line record
 **/
int main(int argc, char **argv) {
    unsigned int itnum = 0;
    std::cout << "// This file is automatically generated by " << argv[0] << ". All changes to this file will be lost after recompilation!!!\n"; 
    // std::cout << "namespace bcfhrec {\n";
    
    std::cout << "#ifndef bcf_formats_step1_INCLUDED\n";
    std::cout << "#define bcf_formats_step1_INCLUDED\n"; 
    std::cout << "#include<array>\n";
    std::cout << "#include<ostream>\n";
    std::cout << "#include<string>\n";
    std::cout << "#include<vector>\n";
    std::cout << "namespace bcfrec {\n";

    std::cout << "static const unsigned int FILTER_NUM = " << FILTER_VEC.size() << ";\n";
    
    std::cout << "enum FILTER_ENUM {\n";
    for (auto filter : FILTER_VEC) {
        std::string format_key = filter.first;
        std::cout << "    " << format_key << ",\n";
    }
    std::cout << "};\n";
    
    std::cout << "const char *const FILTER_IDS[] = {\n";
    for (auto filter : FILTER_VEC) {
        std::string filter_key = filter.first;
        std::cout << "    " << "\"" << filter_key << "\"" << ",\n";
    }
    std::cout << "};\n";
    
    std::cout << "const char *const FILTER_LINES[] = {\n";
    for (auto filter : FILTER_VEC) {
        std::string filter_key = filter.first;
        std::cout << "    "
                  << "\"" << "##FILTER=" 
                          << "<" << "ID=" << filter.first 
                                 << ",Description=" << "\\\"" << filter.second << "\\\"" 
                          << ">" 
                  << "\"" << ",\n";
    }
    std::cout << "};\n";
    
#ifdef OUTPUT_EXTRA_VCF_INFO 
    std::cout << "const char *const FILTER_DESCS[] = {\n";
    for (auto filter : FILTER_VEC) {
        std::string filter_desc = filter.second;
        std::cout << "    " << "\"" << filter_desc << "\"" << ",\n";
    }
    std::cout << "};\n";
#endif
    
    std::cout << "static const unsigned int FORMAT_NUM = " << FORMAT_VEC.size() << ";\n";
    
    std::cout << "enum FORMAT_ENUM {\n";
    for (auto fmt : FORMAT_VEC) {
        std::cout << "    " << fmt.id << ",\n";
    }
    std::cout << "};\n";
   
    std::cout << "const char *const FORMAT_STR_PER_REC = \"";  
    itnum = 0;
    for (auto fmt : FORMAT_VEC) {
        if (itnum) {
            std::cout << ":";
        }
        std::cout << fmt.id;
        itnum++;
    }
    std::cout << "\";\n";
     
    std::cout << "const char *const FORMAT_IDS[] = {\n";
    for (auto fmt : FORMAT_VEC) {
        std::cout << "    " << "\"" << fmt.id << "\"" << ",\n";
    }
    std::cout << "};\n";

    std::cout << "struct BcfFormat {\n";
    for (auto fmt : FORMAT_VEC) {
        if (0 == fmt.number || 1 == fmt.number) {
            if (0 == fmt.number) { assert (CPP_DATA_STRING[fmt.type] == std::string("bool")); }
            std::cout << "    " << CPP_DATA_STRING[fmt.type] << " " << fmt.id << " = " << CPP_DATA_VALUES[fmt.type] << ";" << "\n";
        } else if (1 < fmt.number) {
            std::cout << "    std::array <" << CPP_DATA_STRING[fmt.type] << ", " << fmt.number << ">" << fmt.id << ";" << "\n";
        } else {
            std::cout << "    std::vector<" << CPP_DATA_STRING[fmt.type] << ">" << fmt.id << ";" << "\n";
        }
    }
    std::cout << "};\n";
    
    std::cout << "static int streamAppendBcfFormat(std::string & outstring, const BcfFormat & fmt) {\n";
    itnum = 0;
    for (auto fmt : FORMAT_VEC) {
        if (itnum) { 
            std::cout << "    outstring += \":\";" << ";\n";
        }
        if (BCF_SEP == fmt.type) {
            std::cout << "    outstring += std::string(FORMAT_IDS[" << itnum << "]) + \"_\";\n";
        } else if (0 == fmt.number || 1 == fmt.number) {
            std::cout << "    outstring += " << (fmt.type == BCF_STRING ? "" : "std::to_string") << "(fmt." << fmt.id << ");\n";
        } else {
            std::cout << "    for (unsigned int i = 0; i < fmt." << fmt.id << ".size(); i++) {\n";
            
            //std::string customcode = ((BCF_STRING == fmt.type) ? 
            //        (std::string(" fmt.") + fmt.id + "[i].size() == 0 ? \".\" : fmt." + fmt.id + "[i]") 
            //        : (std::string("fmt.") + std::string(fmt.id) + "[i]"));
            std::cout << "        if (0 != i) { outstring += \",\"; }; outstring += " << (fmt.type == BCF_STRING ? "" : "std::to_string") << "(" << 
                    // customcode 
                    "fmt." << fmt.id << "[i]"
                    << ");\n";
            std::cout << "    };\n";
            if (BCF_STRING == fmt.type) {
                std::cout << "        if (fmt. " << fmt.id << ".size() == 1 && fmt." << fmt.id << "[0].size() == 0) { outstring += \",\"; };";
            }
            if (BCF_NUM_D == fmt.type) {
                std::cout << "        if (fmt. " << fmt.id << "[0].size() == 0) { outstring += \"0\"; };";
            }
        }
        itnum++;
    }
    std::cout << "};\n";
    
    std::cout << "const char *const FORMAT_LINES[] = {\n";
    for (auto fmt : FORMAT_VEC) {
        std::cout << "    " 
                  << "\"" << "##FORMAT=" 
                          << "<" << "ID=" << fmt.id
                                 << ",Number=" << vcf_number_to_header_str(fmt.number)
                                 << ",Type=" << BCF_DATA_STRING[fmt.type]
                                 << ",Description=" << "\\\"" << fmt.description << "\\\"" 
                          << ">" 
                  << "\"" << ",\n";
    }
    std::cout << "};\n";
    
#ifdef OUTPUT_EXTRA_VCF_FORMAT
    std::cout << "const int FORMAT_NUMBERS[] = {\n";
    for (auto filter : FORMAT_KEY_TO_CONTENT_MAP) {
        auto format_content = filter.second;
        std::cout << "    " << format_content.number << ",\n";
    }
    std::cout << "};\n";
    
    std::cout << "const char *const FORMAT_DATA_TYPES[] = {\n";
    for (auto filter : FORMAT_KEY_TO_CONTENT_MAP) {
        auto format_content = filter.second;
        std::cout << "    " << "\"" << BCF_DATA_STRING[format_content.type] << "\"" << ",\n";
    }
    std::cout << "};\n";
    
    std::cout << "  const char *const FORMAT_DESCRIPTIONS[] = {\n";
    for (auto filter : FORMAT_KEY_TO_CONTENT_MAP) {
        auto format_content = filter.second;
        std::cout << "    " << "\"" << format_content.description << "\"" << ",\n";
    }
    std::cout << "  };\n";
#endif
    
    std::cout << "};\n";
    std::cout << "#endif\n";
}

