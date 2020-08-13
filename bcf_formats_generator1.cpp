// This is the framework for generating VCF/BCF header and its associated data on each line.
// Manual update of VCF/BCF header and its associated data on each line is extremely error-prone.
// Therefore, this framework is used.
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "common.hpp"

#define BCF_NUM_A (-1)
#define BCF_NUM_R (-2)
#define BCF_NUM_G (-3)
#define BCF_NUM_D (-4)

enum BCF_DATA_TYPE {
    BCF_STRING,
    BCF_INTEGER,
    BCF_SIG_INT,
    BCF_FLOAT,
    BCF_SEP,
    NUM_BCF_DATA_TYPES
};

const char * CPP_DATA_STRING[] = {
    "std::string",
    "uint32_t",
    "int32_t",
    "float",
    "bool",
    NULL
};

const char * CPP_DATA_VALUES[] = {
    "\"\"",
    "0",
    "0",
    "0",
    "false",
    NULL,
};

const char * BCF_DATA_STRING[] = {
    "String",
    "Integer",
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
    fprintf(stderr, "%d is not a valid BCF type!\n", num);
    abort();
}

const std::vector<std::pair<std::string, std::string>> FILTER_VEC = {
    std::make_pair("noVar",         "Not a variant (for example, when REF and ALT are the same), but still included to get all statistics."),
    std::make_pair("upstreamDel",   "Deletion extended from another upstream deletion"),
    std::make_pair("cad3",          "Less than 3 clean deduped reads"),
    std::make_pair("caf3",          "Less than 3/10000 allele fraction base on clean deduped reads"),
    std::make_pair("s50",           "Less than 50\% of samples have data"),
    std::make_pair("Q10",           "Quality below 10 and no other filters"),
    std::make_pair("Q20",           "Quality below 20 and no other filters"),
    std::make_pair("Q30",           "Quality below 30 and no other filters"),
    std::make_pair("Q40",           "Quality below 40 and no other filters"),
    std::make_pair("Q50",           "Quality below 50 and no other filters"),
    std::make_pair("Q60",           "Quality below 60 and no other filters"),
    std::make_pair("GPBL1",         "For FORMAT/FT: Haplotype position bias on the left  mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBR1",         "For FORMAT/FT: Haplotype position bias on the right mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBLR1",        "For FORMAT/FT: Haplotype position bias on left and right mapping coordinates for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBL2",         "For FORMAT/FT: Diplotype position bias on the left  mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBR2",         "For FORMAT/FT: Diplotype position bias on the right mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBLR2",        "For FORMAT/FT: Diplotype position bias on left and right mapping coordinates for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("DB1",           "For FORMAT/FTS: Deduplication bias for duped base reads, meaning a high portion of variant evidence is from singleton families"),
    std::make_pair("DB2",           "For FORMAT/FTS: Deduplication bias for deduped consensus families, meaning a high portion of variant evidence is from singleton families"),
    std::make_pair("MB1",           "For FORMAT/FTS: Mismatch bias for duped base reads, meaning reads supporting the ALT allele have many more mismatches than reads supporting other alleles."),
    std::make_pair("MB2",           "For FORMAT/FTS: Mismatch bias for deduped consensus families, meaning families supporting the ALT allele have many more mismatches than families supporting other alleles."),
    std::make_pair("PB1L",          "For FORMAT/FTS: Position bias on the left  mappping coordinate for duped base reads"),
    std::make_pair("PB1R",          "For FORMAT/FTS: Position bias on the right mappping coordinate for duped base reads"),
    std::make_pair("PB2L",          "For FORMAT/FTS: Position bias on the left  mappping coordinate for deduped consensus families"),
    std::make_pair("PB2R",          "For FORMAT/FTS: Position bias on the right mappping coordinate for deduped consensus families"),
    std::make_pair("SB1",           "For FORMAT/FTS: Strand bias for duped base reads"),
    std::make_pair("SB2",           "For FORMAT/FTS: Strand bias for deduped consensus families"),
    std::make_pair("QTD1",          "For FORMAT/FTS: Quality-threshold difference for duped base reads, meaning variant evidence is from outlier reads"),
    std::make_pair("QTD2",          "For FORMAT/FTS: Quality-threshold difference for deduped consensus families, meaning variant evidence is from outlier families"), 
    std::make_pair("DBthis",        "For FORMAT/FTS: Allele-fraction-based deduplication bias for this ALT allele, meaning variant evidence is from abnormally small families"),
    std::make_pair("DBrest",        "For FORMAT/FTS: Allele-fraction-based deduplication bias for the non-ALT alelles, meaning variant evidence is from abnormally large families")
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

// philosophy : record signal instead of noise if possible.
const std::vector<BcfFormatStruct> FORMAT_VEC = {
    BcfFormatStruct("GT"    , 1,         BCF_STRING,  "Genotype (this is a guess for the GT of tumor cells if --somaticGT was set to true)"),
    BcfFormatStruct("GQ"    , 1,         BCF_INTEGER, "Genotype Quality"),
    BcfFormatStruct("HQ"    , 2,         BCF_INTEGER, "Haplotype Quality"),
    BcfFormatStruct("DP"    , 1,         BCF_INTEGER, "Fragment depth supporting any allele [consensus family, deduped]"),
    BcfFormatStruct("FA"    , 1,         BCF_FLOAT,   "Fraction of the ALT allele [consensus family, deduped]"),
    // BcfFormatStruct("FR"       , 1, BCF_FLOAT,   "Fraction of the REF allele [consensus family, deduped]"),
    // BcfFormatStruct("BQ"       , 1, BCF_INTEGER, "Root mean square (RMS) base quality of the ALT [base read, duped]"), 
    //BcfFormatStruct("MQ"       , 1, BCF_INTEGER, "Root mean square (RMS) mapping quality of the ALT [base read, duped]"), 
    //BcfFormatStruct("FT"       , 1, BCF_STRING,  "Sample genotype filter indicating if this genotype was 'called' (similar in concept to the FILTER field). "
    //                                             "Again, use PASS to indicate that all filters have been passed, a semi-colon separated list of codes for filters "
    //                                             "that fail, or ‘.’ to indicate that filters have not been applied. "
    //                                             "These values should be described in the meta-information in the same way as FILTERs "
    //                                             "(String, no white-space or semi-colons permitted)"),
    
    //BcfFormatStruct("__A1"     , 1, BCF_SEP,     "Dummy header (separator) for genotype-related FORMAT fields. The description of each associated FORMAT is shown below. "),
    //BcfFormatStruct("FTS"      , 1, BCF_STRING,  "Sample variant filter indicating if this variant was 'called' (similar in concept to the FILTER field). "
    //                                             "Again, use PASS to indicate that all filters have been passed, "
    //                                             "an amperstand-separated list of codes for filters that fail, "
    //                                             "or '.' to indicate that filters have not been applied. "
    //                                             "These values should be described in the meta-information in the same way as FILTERs. "
    //                                             "No white-space, semi-colons, or amperstand permitted."),
    //BcfFormatStruct("FTSV"     , BCF_NUM_D, BCF_INTEGER, "Percent bias values for the FTS strings"),
     
    // BcfFormatStruct("FA1"   , BCF_NUM_A, BCF_FLOAT,   "Normalized fraction of the ALT allele after throwing out outlier signals."),
    
    //BcfFormatStruct("__A2"     , 1, BCF_SEP,     "Dummy header (separator) for genotype-related FORMAT fields. The description of each associated FORMAT is shown below. "),
    
    //BcfFormatStruct("OType"    , 1, BCF_STRING,  "The non-reference allele type with the most evidence other than the ALT allele type of this record"),
    //BcfFormatStruct("ORAQs"    , 2, BCF_FLOAT,   "Qualities of OType and reference allele type"),
        
    BcfFormatStruct("__Aa"  , 1,         BCF_SEP,     "Preparation statistics for segment biases at this position."),
    BcfFormatStruct("APDP"  , 1,         BCF_INTEGER, "Total segment depth"),
    BcfFormatStruct("APGap" , 4,         BCF_INTEGER, "Estimated average numbers of insertion to the left and right sides (LRS) and of deletion to the LRS."),
    BcfFormatStruct("APXM"  , 1,         BCF_INTEGER, "Average number of mismatches."),
    BcfFormatStruct("APLRI" , 4,         BCF_INTEGER, "Summed distance to left insert end and the number of such inserts, and similarly for right insert end."),

    BcfFormatStruct("__Ab"  , 1,         BCF_SEP,     "Threshold for each type of bias (tier-1 means weak bias and tier-2 means strong bias)."),

    BcfFormatStruct("AEPT"  , 2,         BCF_INTEGER, "Number of bases to read-segment end for tier-1 and tier-2 edge position bias."),
    BcfFormatStruct("AXMT"  , 2,         BCF_INTEGER, "Number of mismatches on read-segment for tier-1 and tier-2 mismatch bias."),
    BcfFormatStruct("ALIT"  , 4,         BCF_INTEGER, "Number of bases to left insert end for tier-1 and tier-2 left insert bias."),
    BcfFormatStruct("ARIT"  , 4,         BCF_INTEGER, "Number of bases to right insert end for tier-1 and tier-2 right insert bias."),
    
    BcfFormatStruct("__A1"  , 1,         BCF_SEP,     "Depths of the raw sequencing segments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("ADPff" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R1-forward orientation and strand."),
    BcfFormatStruct("ADPfr" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R2-reverse orientation and strand."),
    BcfFormatStruct("ADPrf" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R2-forward orientation and strand."),
    BcfFormatStruct("ADPrr" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R1-reverse orientation and strand."),
    BcfFormatStruct("__A2"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("ABQ1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 base quality bias."),
    BcfFormatStruct("ABQ2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 base quality bias."),
    BcfFormatStruct("AEP1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 edge position bias."),
    BcfFormatStruct("AEP2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 edge position bias."),
    BcfFormatStruct("AXM1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 mismatch bias."),
    BcfFormatStruct("AXM2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 mismatch bias."),
    BcfFormatStruct("__A3"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("ALI1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("ALI2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the right side."),
    BcfFormatStruct("ARI1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("ARI2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the right side."),
    
    BcfFormatStruct("__a1"  , 1,         BCF_SEP,     "Depths of the raw sequencing segments for the REF and each ALT allele."),
    BcfFormatStruct("aDPff" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R1-forward orientation and strand."),
    BcfFormatStruct("aDPfr" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R2-reverse orientation and strand."),
    BcfFormatStruct("aDPrf" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R2-forward orientation and strand."),
    BcfFormatStruct("aDPrr" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R1-reverse orientation and strand."),
    BcfFormatStruct("__a2"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("aBQ1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 base quality bias."),
    BcfFormatStruct("aBQ2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 base quality bias."),
    BcfFormatStruct("aEP1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 edge position bias."),
    BcfFormatStruct("aEP2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 edge position bias."),
    BcfFormatStruct("aXM1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 mismatch bias."),
    BcfFormatStruct("aXM2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 mismatch bias."),
    BcfFormatStruct("__a3"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("aLI1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("aLI2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the right side."),
    BcfFormatStruct("aRI1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("aRI2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the right side."),
    
    BcfFormatStruct("__Bb"  , 1,         BCF_SEP,     "Non-deduped depths of the fragments."),
    BcfFormatStruct("BDPf"  , 2,         BCF_INTEGER, "Non-deduped fragment depth on the forward strand for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("BDPr"  , 2,         BCF_INTEGER, "Non-deduped fragment depth on the reverse strand for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("bDPf"  , BCF_NUM_R, BCF_INTEGER, "Non-deduped fragment depth on the forward strand for the REF allele and each ALT allele."),
    BcfFormatStruct("bDPr"  , BCF_NUM_R, BCF_INTEGER, "Non-deduped fragment dpeth on the reverse strand for the REF allele and each ALT allele."),
    
    // (10 duped reads and 0.8 id) (2 duped reads and 0.8 id)
    BcfFormatStruct("__C1"  , 1,         BCF_SEP,     "Deduped depths of the unique molecular fragments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("CDP1f" , 2,         BCF_INTEGER, "Total deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("CDP2f" , 2,         BCF_INTEGER, "Very weak consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("CDP3f" , 2,         BCF_INTEGER, "Weak consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("C1DPf" , 2,         BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("CDPMf" , 2,         BCF_INTEGER, "Major duped fragment depth on the forward read orientation."),
    BcfFormatStruct("CDPmf" , 2,         BCF_INTEGER, "Minor duped fragment depth on the forward read orientation."),
    
    BcfFormatStruct("__C2"  , 1,         BCF_SEP,     "Deduped depths of the unique molecular fragments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("CDP1r" , 2,         BCF_INTEGER, "Total deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("CDP2r" , 2,         BCF_INTEGER, "Very weak consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("CDP3r" , 2,         BCF_INTEGER, "Weak consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("C1DPr" , 2,         BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("CDPMr" , 2,         BCF_INTEGER, "Major duped fragment depth on the reverse read orientation."),
    BcfFormatStruct("CDPmr" , 2,         BCF_INTEGER, "Minor duped fragment depth on the reverse read orientation."),
    
    BcfFormatStruct("__c1"  , 1,         BCF_SEP,     "Deduped depths of the unique molecular fragments for the REF allele and each ALT allele."),
    BcfFormatStruct("cDP1f" , BCF_NUM_R, BCF_INTEGER, "Total deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("cDP2f" , BCF_NUM_R, BCF_INTEGER, "Very weak consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("cDP3f" , BCF_NUM_R, BCF_INTEGER, "Weak consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("c1DPf" , BCF_NUM_R, BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("cDPMf" , BCF_NUM_R, BCF_INTEGER, "Major duped fragment depth on the forward read orientation."),
    BcfFormatStruct("cDPmf" , BCF_NUM_R, BCF_INTEGER, "Minor duped fragment depth on the forward read orientation."),
    
    BcfFormatStruct("__c2"  , 1,         BCF_SEP,     "Deduped depths of the unique molecular fragments for the REF allele and each ALT allele."),
    BcfFormatStruct("cDP1r" , BCF_NUM_R, BCF_INTEGER, "Total deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("cDP2r" , BCF_NUM_R, BCF_INTEGER, "Very weak consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("cDP3r" , BCF_NUM_R, BCF_INTEGER, "Weak consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("c1DPr" , BCF_NUM_R, BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("cDPMr" , BCF_NUM_R, BCF_INTEGER, "Major duped fragment depth on the reverse read orientation."),
    BcfFormatStruct("cDPmr" , BCF_NUM_R, BCF_INTEGER, "Minor duped fragment depth on the reverse read orientation."),
    
    BcfFormatStruct("__Dd"  , 1,         BCF_SEP,     "Duplex depths of the original double-strand molecular fragments."),
    BcfFormatStruct("DDP1"  , 2,         BCF_INTEGER, "Duplex depth with allele diaagreement on the two strands for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("DDP2"  , 2,         BCF_INTEGER, "Duplex depth with allele agreement on the two strands for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("dDP1"  , BCF_NUM_R, BCF_INTEGER, "Duplex depth with allele disagreement on the two strands for the REF allele and each ALT allele."),
    BcfFormatStruct("dDP2"  , BCF_NUM_R, BCF_INTEGER, "Duplex depth with allele agreement on the two strands for the REF allele and each ALT allele."),
    
    BcfFormatStruct("__e1"   , 1,         BCF_SEP,     "Quality-related variables."),
    BcfFormatStruct("aSBQf" , BCF_NUM_R, BCF_INTEGER, "Summed sequencing-segment base quality on the forward strand."),
    BcfFormatStruct("aSBQr" , BCF_NUM_R, BCF_INTEGER, "Summed sequencing-segment base quality on the reverse strand."),     
    BcfFormatStruct("aBQQ"  , BCF_NUM_R, BCF_SIG_INT, "Variant quality capped by raw base qualities."),
    BcfFormatStruct("bMQ"   , BCF_NUM_R, BCF_SIG_INT, "Root-mean-square mapping quality."),
    // BcfFormatStruct("bMQQ"  , BCF_NUM_R, BCF_INTEGER, "Duplex depth with allele disagreement on the two strands for (all alleles) and (the padded deletion allele)."),
   
    BcfFormatStruct("__e2"  , 1,         BCF_SEP,     "Quality-related variables assuming read supports are IID (IID: independent and identically distributed) for duped reads."),
    BcfFormatStruct("bIAQb" , BCF_NUM_R, BCF_SIG_INT, "IID allele quality maximized with IAD and IDQ."),
    BcfFormatStruct("bIADb" , BCF_NUM_R, BCF_INTEGER, "IID allele depth (number of reads supporting each ALT)"),
    BcfFormatStruct("bIDQb" , BCF_NUM_R, BCF_SIG_INT, "IDD quality threshold per read support."),
    
    BcfFormatStruct("__e3"  , 1,         BCF_SEP,     "Quality-related variables assuming read supports are IID for deduped reads."),
    BcfFormatStruct("cIAQf" , BCF_NUM_R, BCF_SIG_INT, "IID allele quality maximized with IAD and IDQ on the forward read orientation."),
    BcfFormatStruct("cIADf" , BCF_NUM_R, BCF_INTEGER, "IID allele depth (number of reads supporting each ALT) on the forward read orientation."),
    BcfFormatStruct("cIDQf" , BCF_NUM_R, BCF_SIG_INT, "IDD quality threshold per read support on the forward read orientation."),
    
    BcfFormatStruct("__e4"  , 1,         BCF_SEP,     "Quality-related variables assuming read supports are IID for deduped reads."),
    BcfFormatStruct("cIAQr" , BCF_NUM_R, BCF_SIG_INT, "IID allele quality maximized with IAD and IDQ on the reverse read orientation."),
    BcfFormatStruct("cIADr" , BCF_NUM_R, BCF_INTEGER, "IID allele depth (number of reads supporting each ALT) on the reverse read orientation."),
    BcfFormatStruct("cIDQr" , BCF_NUM_R, BCF_SIG_INT, "IDD quality threshold per read support on the reverse read orientation."),
    
    BcfFormatStruct("__e5"  , 1,         BCF_SEP,     "Power-law variant qualities."),
    //BcfFormatStruct("aPLQ"  , BCF_NUM_A, BCF_INTEGER, "The sequencing-segment power-law variant quality."),
    BcfFormatStruct("bIAQ"  , BCF_NUM_R, BCF_SIG_INT, "The duped fragment binomial variant quality by assuming statistical independence."),
    BcfFormatStruct("cIAQ"  , BCF_NUM_R, BCF_SIG_INT, "The deduplicated fragment binomial variant quality by assuming statistical independence."),
    
    BcfFormatStruct("__e6"  , 1,         BCF_SEP,     "Power-law variant quality statistics for deduped read fragments."),
    BcfFormatStruct("cPLQ1" , BCF_NUM_R, BCF_SIG_INT, "The deduplicatd read fragment power-law variant allele quality."),
    BcfFormatStruct("cVQ1"  , BCF_NUM_R, BCF_SIG_INT, "The final variant quality computed with deduplicated read fragments."),
    BcfFormatStruct("cDP1v" , BCF_NUM_R, BCF_INTEGER, "The effective number of deduplicated read fragments supporting each allele multiplied by 100."),
    BcfFormatStruct("CDP1v" , 2        , BCF_INTEGER, "The effective number of deduplicated read fragments supporting all alleles multiplied by 100."),
    
    BcfFormatStruct("__e7"  , 1,         BCF_SEP,     "Power-law variant quality statistics for consensus barcode families."),
    BcfFormatStruct("cPLQ2" , BCF_NUM_R, BCF_SIG_INT, "The single-strand-consensus-sequence (SSCS) UMI-barcoded power-law variant allele quality."),
    BcfFormatStruct("cVQ2"  , BCF_NUM_R, BCF_SIG_INT, "The final variant allele quality computed with SSCS UMI-barcoded families"),
    BcfFormatStruct("cDP2v" , BCF_NUM_R, BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting each allele multiplied by 100."),
    BcfFormatStruct("CDP2v" , 2        , BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting all alleles multiplied by 100."), // TODO: implement?
    // BcfFormatStruct("dPLQ"  , BCF_NUM_A, BCF_INTEGER, "The double-strand-consensus-sequence (DSCS) power-law variant quality."),
    BcfFormatStruct("CONTQ" , BCF_NUM_R, BCF_SIG_INT, "Likelihood of the variant signal if the variant signal is contaminated."),
    
    BcfFormatStruct("__gap"  , 1,        BCF_SEP,     "InDel-related information."), 
    BcfFormatStruct("gapNf"  ,BCF_NUM_D, BCF_INTEGER, "Number of InDel sequences on the forward read orientation."),
    BcfFormatStruct("gapNr"  ,BCF_NUM_D, BCF_INTEGER, "Number of InDel sequences on the reverse read orientation."),
    BcfFormatStruct("gapSeq" ,BCF_NUM_D, BCF_STRING,  "InDel sequences"),
    BcfFormatStruct("gapbAD1",BCF_NUM_D, BCF_INTEGER, "Duped read count of each gapSeq"),
    BcfFormatStruct("gapcAD1",BCF_NUM_D, BCF_INTEGER, "Deduped read count of each gapSeq"),
    
    BcfFormatStruct("bDPa"   ,BCF_NUM_R, BCF_INTEGER, "Number of non-deduplicated fragments supporting each ALT allele which is more precise for InDels."),
    BcfFormatStruct("cDP1a"  ,BCF_NUM_R, BCF_INTEGER, "Number of deduplicated fragments supporting each ALT allele which is more precise for InDels"),
    BcfFormatStruct("gapSa"  ,BCF_NUM_R, BCF_STRING,  "InDel string of each allele"),
    BcfFormatStruct("VTI"    ,BCF_NUM_R, BCF_INTEGER, "Variant-type ID of each allele."),
    BcfFormatStruct("VTD"    ,BCF_NUM_R, BCF_STRING,  "Variant-type description of each allele"),
    
    BcfFormatStruct("bHap"  , 1,         BCF_STRING,  "Duped forward&reverse linkage in the format of ((position&variantType)...depth)... "
                                                      "where ()... means more elements following the format in the preceding parenthesis. "),
    BcfFormatStruct("cHap"  , 1,         BCF_STRING,  "Dedup forward&reverse linkage in the format of ((position&variantType)...depth)... "
                                                      "where ()... means more elements following the format in the preceding parenthesis. "),    
    
    // BcfFormatStruct("__ea"     , 1, BCF_SEP,     "Number of indels in forward&reverse strand (gapNum), max cAD diff (gapcADD), total cAD depth (gapcADT), indel sequences (gapSeq), duped read count of each gapSeq, dedup family count of each gapSeq, and duped/deduped sub/all allele read counts (gapDP4)"), 
    // BcfFormatStruct("gapNum"   , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    //BcfFormatStruct("gapcADD"  , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    //BcfFormatStruct("gapcADT"  , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    // BcfFormatStruct("gapSeq"   , BCF_NUM_D, BCF_STRING,  "see above"),
    // BcfFormatStruct("gapbAD1"  , BCF_NUM_D, BCF_INTEGER, "see above"),
    // BcfFormatStruct("gapcAD1"  , BCF_NUM_D, BCF_INTEGER, "see above"),
    // BcfFormatStruct("gapDP4"   , 4, BCF_INTEGER, "see above"),
    //BcfFormatStruct("gapbNRD"  , 4, BCF_INTEGER, "Forward&reverse (02&13) strand-specific duped read count for any insertion&deletion (01&23)"),
    BcfFormatStruct("note"  , 1,         BCF_STRING,  "Additional note as comment for the given variant")
};

/*
 * for FORMAT, we have to keep in sync : 
 * C++ class variables per program, 
 * C++ class variable to stream methods per program, 
 * header lines per vcf file, 
 * tag-names (can be multiple sets of) per vcf line record
 **/
int 
main(int argc, char **argv) {
    unsigned int itnum = 0;
    std::cout << "// This file is automatically generated by " << argv[0] << ". All changes to this file will be lost after recompilation!!!\n"; 
    
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
            std::cout << "    std::array <" << CPP_DATA_STRING[fmt.type] << ", " << fmt.number << ">" << fmt.id << " = {{" << CPP_DATA_VALUES[fmt.type] << "}};" << "\n";
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
            std::cout << "    outstring += std::string(FORMAT_IDS[" << itnum << "]) + \"\";\n";
        } else if (0 == fmt.number || 1 == fmt.number) {
            std::cout << "    outstring += " << (fmt.type == BCF_STRING ? "" : "std::to_string") << "(fmt." << fmt.id << ");\n";
        } else {
            std::cout << "    for (unsigned int i = 0; i < fmt." << fmt.id << ".size(); i++) {\n";
            std::cout << "        if (0 != i) { outstring += \",\"; }; outstring += " << (fmt.type == BCF_STRING ? "" : "std::to_string") 
                    << "(" << "fmt." << fmt.id << "[i]" << ");\n";
            std::cout << "    };\n";
            /*
            if (BCF_STRING == fmt.type) {
                std::cout << "        if (fmt. " << fmt.id << ".size() == 1 && fmt." << fmt.id << "[0].size() == 0) { outstring += \",\"; };";
            }
            if (BCF_NUM_D == fmt.type) {
                std::cout << "        if (fmt. " << fmt.id << "[0].size() == 0) { outstring += \"0\"; };";
            }
            */
        }
        itnum++;
    }
    std::cout << "\n    return 0;};\n";
    
    std::cout << "#include <assert.h>\n";
    std::cout << "static int streamFrontPushBcfFormatR(BcfFormat & dst, const BcfFormat & src) {\n";
    
    for (auto fmt : FORMAT_VEC) {
        if (BCF_NUM_R == fmt.number) {
            std::cout << "    assert(dst." << fmt.id << ".size() == 1 || !fprintf(stderr, \"\%d == 1 failed for " << fmt.id 
                    << "\", dst." << fmt.id << ".size() ) );\n";
            std::cout << "    assert(src." << fmt.id << ".size() == 1 || !fprintf(stderr, \"\%d == 1 failed for " << fmt.id 
                    << "\", src." << fmt.id << ".size() ) );\n";
            std::cout << "    auto " << fmt.id << "_tmp = dst." << fmt.id << "[0];\n";
            std::cout << "    dst." << fmt.id << "[0] = src." << fmt.id << "[0];\n";
            std::cout << "    dst." << fmt.id << ".push_back(" << fmt.id << "_tmp);\n";
        }
    }
    std::cout << "\n    return 0;};\n";

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
    return 0;
}

