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
    "int32_t", // "uint32_t",
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
    
    /*
    std::make_pair("GPBL1",         "For FORMAT/FT: Haplotype position bias on the left  mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBR1",         "For FORMAT/FT: Haplotype position bias on the right mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBLR1",        "For FORMAT/FT: Haplotype position bias on left and right mapping coordinates for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBL2",         "For FORMAT/FT: Diplotype position bias on the left  mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBR2",         "For FORMAT/FT: Diplotype position bias on the right mapping coordinate for raw reads with unmerged R1 and R2 ends"),
    std::make_pair("GPBLR2",        "For FORMAT/FT: Diplotype position bias on left and right mapping coordinates for raw reads with unmerged R1 and R2 ends"),
    */
    std::make_pair("ASI",           "For FORMAT/FTS: Stranded insert bias, meaning the most-supported strand has abnormal insert size."),
    std::make_pair("AXMB",          "For FORMAT/FTS: Absolute mismatch bias, meaning the variant is suppported by reads with a high number of mismatches"),
    std::make_pair("ABQB",          "For FORMAT/FTS: Absolute base-quality (BQ) bias, meaning the variant is suppported by reads with low base qualities at the variant site"),

    std::make_pair("PFB",           "For FORMAT/FTS: Passing-filter bias, meaning the variant allele is supported by reads with low base qualities at the variant site and/or with high number of mismatches relative to all alleles."),
    
    std::make_pair("DB1",           "For FORMAT/FTS: Deduplication bias for the under-amplification of variant reads, meaning the variant is under-amplified by PCR relative to all alleles"),
    std::make_pair("DB2",           "For FORMAT/FTS: Deduplication bias for the over-amplification of variant reads, meaning the variant is over-amplified by PCR relative to all alleles"),
    
    std::make_pair("BB1L",          "For FORMAT/FTS: Alignment bias on the left  mappping coordinate of the sequenced segment relative to all alleles."),
    std::make_pair("BB1R",          "For FORMAT/FTS: Alignment bias on the right mappping coordinate of the sequenced segment relative to all alleles."),
    std::make_pair("PB1L",          "For FORMAT/FTS: Position bias on the left  mappping coordinate of the sequenced segment relative to all alleles."),
    std::make_pair("PB1R",          "For FORMAT/FTS: Position bias on the right mappping coordinate of the sequenced segment relative to all alleles."),
    std::make_pair("PB2L",          "For FORMAT/FTS: Position bias on the left mapping coordinate of the insert relative to all alleles."),
    std::make_pair("PB2R",          "For FORMAT/FTS: Position bias on the right mapping coordinate of the insert relative to all alleles"),
    
    std::make_pair("SB1",           "For FORMAT/FTS: Strand bias relative to all alleles"),
    // std::make_pair("ROB0",          "For FORMAT/FTS: Read-orientation bias using all deduplicated reads relative to all alleles"),
    std::make_pair("ROB1",          "For FORMAT/FTS: Read-orientation bias using deduplicated reads families passing the base-quality threshold for generating barcode families relative to all alleles"),
    std::make_pair("ROB2",          "For FORMAT/FTS: Read-orientation bias using tier-2 barcode families relative to all alleles"),
    // std::make_pair("ROB3",          "For FORMAT/FTS: Read-orientation bias using tier-3 barcode families relative to all alleles"),
        
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
    BcfFormatStruct("FT"       , 1, BCF_STRING,  "Sample genotype filter indicating if this genotype was 'called' (similar in concept to the FILTER field). "
                                                 "Again, use PASS to indicate that all filters have been passed, a semi-colon separated list of codes for filters "
                                                 "that fail, or ‘.’ to indicate that filters have not been applied. "
                                                 "These values should be described in the meta-information in the same way as FILTERs "
                                                 "(String, no white-space or semi-colons permitted)"),
    BcfFormatStruct("FTS"      , 1, BCF_STRING,  "Sample variant filter indicating if this variant was 'called' (similar in concept to the FILTER field). "
                                                 "Again, use PASS to indicate that all filters have been passed, "
                                                 "an amperstand-separated list of codes for filters that fail, "
                                                 "or '.' to indicate that filters have not been applied. "
                                                 "These values should be described in the meta-information in the same way as FILTERs. "
                                                 "No white-space, semi-colons, or amperstand permitted."),
    
    BcfFormatStruct("__A_"  , 1,         BCF_SEP,     "Summary statistics."),
    BcfFormatStruct("DP"    , 1,         BCF_INTEGER, "Fragment depth of coverage with duplicates removed. "),
    BcfFormatStruct("AD"    , BCF_NUM_R, BCF_INTEGER, "Fragment depth supporting the ALT allele with duplicates removed. "),
    BcfFormatStruct("bDP"   , 1,         BCF_INTEGER, "Fragment depth of coverage with duplicates kept. "),
    BcfFormatStruct("bAD"   , BCF_NUM_R, BCF_INTEGER, "Fragment depth supporting the ALT allele with duplicates kept. "),
    BcfFormatStruct("c2DP"  , 1,         BCF_INTEGER, "Consensus barcode family depth of coverage using tier-2 thresholds for grouping fragments into families. "),
    BcfFormatStruct("c2AD"  , BCF_NUM_R, BCF_INTEGER, "Consensus barcode family depth of coverage supporting the ALT allele using tier-2 thresholds for grouping fragments into families. "),
    
    // BcfFormatStruct("BQ"       , 1, BCF_INTEGER, "Root mean square (RMS) base quality of the ALT [base read, duped]"), 
    // BcfFormatStruct("MQ"       , 1, BCF_INTEGER, "Root mean square (RMS) mapping quality of the ALT [base read, duped]"), 

    BcfFormatStruct("__Aa"  , 1,         BCF_SEP,     "Preparation statistics for segment biases at this position."),
    BcfFormatStruct("APDP"  , 1+4+4+1,   BCF_INTEGER, "Total segment depth, "
                              "segment depths within the indel length of insertion/deletion, segment depths within the tandem-repeat track length of insertion/deletion, "
                              "segment depth of high quality, "
                              "sum of squares of insertion lengths, sum of squares of deletion lengths, sum of the inverses of insertion lengths, sum of the inverses of deletion lengths, "
                              "and total segment depth of PCR amplicons."),
    BcfFormatStruct("APXM"  , 4+1+3,     BCF_INTEGER, "Total number of mismatches and total number of gap openings."),
    
    BcfFormatStruct("__Ab"  , 1,         BCF_SEP,     "Preparation statistics for segment biases at this position."),
    BcfFormatStruct("APLRID", 4,         BCF_INTEGER, "Total number of mismatches and total number of gap openings."),
    BcfFormatStruct("APPB"  , 4+4+2,     BCF_INTEGER, "Preparation statistics for position bias."),
    BcfFormatStruct("APLRI" , 4,         BCF_INTEGER, "Summed distance to left insert end and the number of such inserts, and similarly for right insert end."),
    BcfFormatStruct("APLRP" , 4,         BCF_INTEGER, "Summed distance to left and right ends, summed insertion length, and summed deletion length."),
    
    BcfFormatStruct("__Ac"  , 1,         BCF_SEP,     "Threshold for each type of bias (tier-1 means weak bias and tier-2 means strong bias)."),

    BcfFormatStruct("AXMT"  , 2,         BCF_INTEGER, "Number of mismatches on read-segment above which there is tier-1 and tier-2 mismatch bias."),
    BcfFormatStruct("ALRIT" , 4,         BCF_INTEGER, "Number of bases to left (01) and right (23) insert ends above which there is tier-1 and tier-2 insert bias."),
    BcfFormatStruct("ALRIt" , 4,         BCF_INTEGER, "Number of bases to left (01) and right (23) insert ends below which there is tier-1 and tier-2 insert bias."),
    BcfFormatStruct("ALRPt" , 4,         BCF_INTEGER, "Number of bases to left (01) and right (23) read-segment ends below which there is tier-1 and tier-2 position bias."),
    BcfFormatStruct("ALRBt" , 4,         BCF_INTEGER, "Base alignment quality (BAQ) to left (01) and right (23) read-segment ends below which there is tier-1 and tier-2 position bias."),

    BcfFormatStruct("__AQ"  , 1,         BCF_SEP,     "Statistics of the raw sequencing segments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("AMQs"  , 2,         BCF_INTEGER, "Raw sequencing segment sum of mapping qualities."),
    //  BcfFormatStruct("AXMp1" , 2,         BCF_INTEGER, "Raw sequencing segment sum of 100x depths normalized by the arithmetic inverse of mismatches in a 150-bp region window "
   //                           "(zero or one mismatch counts as 100, two mismatches count as 100/2=50, threee mismatches count as 100/3=33, etc.)."),
    
    BcfFormatStruct("A1BQf" , 2,         BCF_INTEGER, "Summed sequencing-segment base quality on the forward strand."),
    BcfFormatStruct("A1BQr" , 2,         BCF_INTEGER, "Summed sequencing-segment base quality on the reverse strand."),     
    
    BcfFormatStruct("__A1"  , 1,         BCF_SEP,     "Depths of the raw sequencing segments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("ADPff" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R1-forward orientation and strand."),
    BcfFormatStruct("ADPfr" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R2-reverse orientation and strand."),
    BcfFormatStruct("ADPrf" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R2-forward orientation and strand."),
    BcfFormatStruct("ADPrr" , 2,         BCF_INTEGER, "Raw sequencing segment depth with the R1-reverse orientation and strand."),
    
    BcfFormatStruct("__A2"  , 1,         BCF_SEP,     "Depths of the raw sequencing segments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("ALP1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 left-side position bias."),
    BcfFormatStruct("ALP2"  , 2,         BCF_INTEGER, "RSEaw sequencing segment depth unaffected by tier-2 left-side position bias."),
    BcfFormatStruct("ALPL"  , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the left-side sequencing-segment end."),
    BcfFormatStruct("ARP1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 right-side position bias."),
    BcfFormatStruct("ARP2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 right-side position bias."),
    BcfFormatStruct("ARPL"  , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the right-side sequencing-segment end."),
    
    BcfFormatStruct("__A3"  , 1,         BCF_SEP,     "Depths of the raw sequencing segments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("ALB1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 left-side position bias."),
    BcfFormatStruct("ALB2"  , 2,         BCF_INTEGER, "RSEaw sequencing segment depth unaffected by tier-2 left-side position bias."),
    BcfFormatStruct("ALBL"  , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the left-side sequencing-segment end."),
    BcfFormatStruct("ARB1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 right-side position bias."),
    BcfFormatStruct("ARB2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 right-side position bias."),
    BcfFormatStruct("ARBL"  , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the right-side sequencing-segment end."),
    
    BcfFormatStruct("__A4"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("ABQ2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 base quality bias."),
    BcfFormatStruct("APF1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 mismatch bias."),
    BcfFormatStruct("APF2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 mismatch bias."),
    BcfFormatStruct("AP1"   , 2,         BCF_INTEGER, "Raw sequencing segment depth of reads passing left and right number of bases threshold of distance."),
    BcfFormatStruct("AP2"   , 2,         BCF_INTEGER, "Raw sequencing segment depth of reads that are either labeled with UMIs or not coming from PCR amplicons."),


    //BcfFormatStruct("AXM1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 relative mismatch bias "
    //                "(for capping quality from tumor-normal comparison)."),
    //BcfFormatStruct("AXM2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 relative mismatch bias "
    //                "(for capping quality from tumor-normal comparison)."),
    
    BcfFormatStruct("__A5"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("ALI1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("ALI2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the right side."),
    BcfFormatStruct("ALIr"  , 2,         BCF_INTEGER, "Raw sequencing segment depth eligible for  left-side reverse-strand bias computation."),
    
    // BcfFormatStruct("ALILf" , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the left-side insert end on the forward read orientation."),
    // BcfFormatStruct("ALILr" , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the left-side insert end on the reverse read orientation."),
    BcfFormatStruct("ARI1"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("ARI2"  , 2,         BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the right side."),
    BcfFormatStruct("ARIf"  , 2,         BCF_INTEGER, "Raw sequencing segment depth eligible for right-side forward-strand bias computation."),
    
    // BcfFormatStruct("ARILf" , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the right-side insert end on the forward read orientation."),
    //BcfFormatStruct("ARILr" , 2,         BCF_INTEGER, "Raw summed distance (number of bases) to the right-side insert end on the reverse read orientation."),
    
    BcfFormatStruct("__aQ"  , 1,         BCF_SEP,     "Statistics of the raw sequencing segments for the REF and each ALT allele."),
    BcfFormatStruct("aMQs"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment sum of mapping qualities."),
    // BcfFormatStruct("aXMp1" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment sum of 100x depths normalized by the arithmetic inverse of mismatches in a 150-bp region window "
    //                          "(zero or one mismatch counts as 100, two mismatches count as 100/2=50, threee mismatches count as 100/3=33, etc.)."),
    BcfFormatStruct("a1BQf" , BCF_NUM_R, BCF_INTEGER, "Summed sequencing-segment base quality on the forward strand."),
    BcfFormatStruct("a1BQr" , BCF_NUM_R, BCF_INTEGER, "Summed sequencing-segment base quality on the reverse strand."),     
    
    BcfFormatStruct("__a1"  , 1,         BCF_SEP,     "Depths of the raw sequencing segments for the REF and each ALT allele."),
    BcfFormatStruct("aDPff" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R1-forward orientation and strand."),
    BcfFormatStruct("aDPfr" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R2-reverse orientation and strand."),
    BcfFormatStruct("aDPrf" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R2-forward orientation and strand."),
    BcfFormatStruct("aDPrr" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth with the R1-reverse orientation and strand."),
    
    BcfFormatStruct("__a2"  , 1,         BCF_SEP,     "As before."),
    
    BcfFormatStruct("aLP1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 left-side position bias."),
    BcfFormatStruct("aLP2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 left-side position bias."),
    BcfFormatStruct("aLPL"  , BCF_NUM_R, BCF_INTEGER, "Raw summed distance (number of bases) to the left-side sequencing-segment end."),
    BcfFormatStruct("aRP1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 right-side position bias."),
    BcfFormatStruct("aRP2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 right-side position bias."),
    BcfFormatStruct("aRPL"  , BCF_NUM_R, BCF_INTEGER, "Raw summed distance (number of bases) to the right-side sequencing-segment end."),
    
    BcfFormatStruct("__a3"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("aLB1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 left-side position bias."),
    BcfFormatStruct("aLB2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 left-side position bias."),
    BcfFormatStruct("aLBL"  , BCF_NUM_R, BCF_INTEGER, "Raw summed distance (number of bases) to the left-side sequencing-segment end."),
    BcfFormatStruct("aRB1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 right-side position bias."),
    BcfFormatStruct("aRB2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 right-side position bias."),
    BcfFormatStruct("aRBL"  , BCF_NUM_R, BCF_INTEGER, "Raw summed distance (number of bases) to the right-side sequencing-segment end."),
    
    BcfFormatStruct("__a4"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("aBQ2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 base quality bias."),
    BcfFormatStruct("aPF1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 relative mismatch bias and base quality bias."),
    BcfFormatStruct("aPF2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 relative mismatch bias and base quality bias."),
    BcfFormatStruct("aP1"   , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth of reads passing left and right number of bases threshold of distance."),
    BcfFormatStruct("aP2"   , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth of reads that are either labeled with UMIs or not coming from PCR amplicons."),
    
    // BcfFormatStruct("aXM1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 relative mismatch bias "
    //                "(for capping quality from tumor-normal comparison)."),
    // BcfFormatStruct("aXM2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 relative mismatch bias "
    //                "(for capping quality from tumor-normal comparison)."),
    
    BcfFormatStruct("__a5"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("aLI1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("aLI2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the right side."),
    BcfFormatStruct("aRI1"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-1 insert bias on the left side."),
    BcfFormatStruct("aRI2"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by tier-2 insert bias on the left side."),
    BcfFormatStruct("aLIr"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth eligible for  left-side reverse-strand bias computation."),
    BcfFormatStruct("aRIf"  , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth eligible for right-side forward-strand bias computation."),
    
    BcfFormatStruct("__B1"  , 1,         BCF_SEP,     "Non-deduped fragment depth for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("BDPf"  , 2,         BCF_INTEGER, "Non-deduped fragment depth on the forward strand for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("BTAf"  , 2,         BCF_INTEGER, "Non-deduped sum of sequenced fragment positions on the forward strand for "
                    "(all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("BTBf"  , 2,         BCF_INTEGER, "Non-deduped sum of sequenced fragment positions near mutations on the forward strand for "
                    "(all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("BDPr"  , 2,         BCF_INTEGER, "Non-deduped fragment depth on the reverse strand for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("BTAr"  , 2,         BCF_INTEGER, "Non-deduped sum of sequenced fragment positions on the reverse strand for "
                    "(all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("BTBr"  , 2,         BCF_INTEGER, "Non-deduped sum of sequenced fragment positions near mutations on the reverse strand for "
                    "(all alleles) and (the padded deletion allele)."),
    
    BcfFormatStruct("__b1"  , 1,         BCF_SEP,     "Non-deduped depths of the fragments for the REF allele and each ATL allele."),
    
    BcfFormatStruct("bDPf"  , BCF_NUM_R, BCF_INTEGER, "See BDPf."),
    BcfFormatStruct("bTAf"  , BCF_NUM_R, BCF_INTEGER, "See BTAf."),
    BcfFormatStruct("bTBf"  , BCF_NUM_R, BCF_INTEGER, "See BTBf."),
    BcfFormatStruct("bDPr"  , BCF_NUM_R, BCF_INTEGER, "See BDPr."),
    BcfFormatStruct("bTAr"  , BCF_NUM_R, BCF_INTEGER, "See BTAr."),
    BcfFormatStruct("bTBr"  , BCF_NUM_R, BCF_INTEGER, "See BTBr."),
    
    BcfFormatStruct("__C1"  , 1,         BCF_SEP,     "Deduped depths of the unique molecular fragments for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("CDP1f" , 2,         BCF_INTEGER, "Nonfiltered      deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("CDP12f", 2,         BCF_INTEGER, "BQ-filtered      deduped depth of the unique molecular fragments on the forward read orientation ."),
    BcfFormatStruct("CDP2f" , 2,         BCF_INTEGER, "Weak consensus   deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("CDP3f" , 2,         BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("C1DPf" , 2,         BCF_INTEGER, "Singleton        deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("CDPMf" , 2,         BCF_INTEGER, "Major duped fragment depth on the forward read orientation."),
    BcfFormatStruct("CDPmf" , 2,         BCF_INTEGER, "Minor duped fragment depth on the forward read orientation."),
    
    BcfFormatStruct("__C2"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("CDP1r" , 2,         BCF_INTEGER, "Nonfiltered      deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("CDP12r", 2,         BCF_INTEGER, "BQ-Filtered      deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("CDP2r" , 2,         BCF_INTEGER, "Weak consensus   deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("CDP3r" , 2,         BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("C1DPr" , 2,         BCF_INTEGER, "Singleton        deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("CDPMr" , 2,         BCF_INTEGER, "Major duped fragment depth on the reverse read orientation."),
    BcfFormatStruct("CDPmr" , 2,         BCF_INTEGER, "Minor duped fragment depth on the reverse read orientation."),
    
    BcfFormatStruct("__c1"  , 1,         BCF_SEP,     "Deduped depths of the unique molecular fragments for the REF allele and each ALT allele."),
    BcfFormatStruct("cDP1f" , BCF_NUM_R, BCF_INTEGER, "Nonfiltered      deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("cDP12f", BCF_NUM_R, BCF_INTEGER, "BQ-filtered      deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("cDP2f" , BCF_NUM_R, BCF_INTEGER, "Weak   consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("cDP3f" , BCF_NUM_R, BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("c1DPf" , BCF_NUM_R, BCF_INTEGER, "Singleton        deduped depth of the unique molecular fragments on the forward read orientation."),
    BcfFormatStruct("cDPMf" , BCF_NUM_R, BCF_INTEGER, "Major duped fragment depth on the forward read orientation."),
    BcfFormatStruct("cDPmf" , BCF_NUM_R, BCF_INTEGER, "Minor duped fragment depth on the forward read orientation."),
    
    BcfFormatStruct("__c2"  , 1,         BCF_SEP,     "As before."),
    BcfFormatStruct("cDP1r" , BCF_NUM_R, BCF_INTEGER, "Nonfiltered      deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("cDP12r", BCF_NUM_R, BCF_INTEGER, "BQ-filtered      deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("cDP2r" , BCF_NUM_R, BCF_INTEGER, "Weak   consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("cDP3r" , BCF_NUM_R, BCF_INTEGER, "Strong consensus deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("c1DPr" , BCF_NUM_R, BCF_INTEGER, "Singleton        deduped depth of the unique molecular fragments on the reverse read orientation."),
    BcfFormatStruct("cDPMr" , BCF_NUM_R, BCF_INTEGER, "Major duped fragment depth on the reverse read orientation."),
    BcfFormatStruct("cDPmr" , BCF_NUM_R, BCF_INTEGER, "Minor duped fragment depth on the reverse read orientation."),
    
    BcfFormatStruct("__Dd"  , 1,         BCF_SEP,     "Duplex depths of the original double-strand molecular fragments."),
    BcfFormatStruct("DDP1"  , 2,         BCF_INTEGER, "Duplex depth with allele diaagreement on the two strands for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("DDP2"  , 2,         BCF_INTEGER, "Duplex depth with allele agreement on the two strands for (all alleles) and (the padded deletion allele)."),
    BcfFormatStruct("dDP1"  , BCF_NUM_R, BCF_INTEGER, "Duplex depth with allele disagreement on the two strands for the REF allele and each ALT allele."),
    BcfFormatStruct("dDP2"  , BCF_NUM_R, BCF_INTEGER, "Duplex depth with allele agreement on the two strands for the REF allele and each ALT allele."),
    
    BcfFormatStruct("__e1"  , 1,         BCF_SEP,     "Error variables inferred from systematically low basecall qualities (BQs)."),
    
    // BcfFormatStruct("a1XM"  , BCF_NUM_R, BCF_INTEGER, "Total number of mismatches normalized with read length."),
    BcfFormatStruct("aBQ"   , BCF_NUM_R, BCF_SIG_INT, "Root-mean-square base quality for sequencing segments."),
    BcfFormatStruct("a2BQf" , BCF_NUM_R, BCF_INTEGER, "Summed squared/32 sequencing-segment base quality on the forward strand."),
    BcfFormatStruct("a2BQr" , BCF_NUM_R, BCF_INTEGER, "Summed squared/32 sequencing-segment base quality on the reverse strand."),     
    BcfFormatStruct("a2XM2" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by absolute mismatch bias."),
    BcfFormatStruct("a2BM2" , BCF_NUM_R, BCF_INTEGER, "Raw sequencing segment depth unaffected by absolute base-specific mismatch bias."),
    BcfFormatStruct("aBQQ"  , BCF_NUM_R, BCF_SIG_INT, "Variant quality capped by systematic error inferred from base qualities."),

    BcfFormatStruct("__e2"  , 1,         BCF_SEP,     "Error variables inferred from systematically low mapping qualities (MAPQs or MQs)."),
    
    BcfFormatStruct("bMQ"   , BCF_NUM_R, BCF_SIG_INT, "Root-mean-square mapping quality for read fragments (max of R1 and R2 MQ is taken)."),
    BcfFormatStruct("aAaMQ" , BCF_NUM_R, BCF_SIG_INT, "Read segment difference of average mapping quality between the ALT allele and the non-ALT alleles."),
    BcfFormatStruct("bNMQ"  , BCF_NUM_R, BCF_INTEGER, "Phred penalty inferred from high-BQ mismatches to the mapping quality (MAPQ). "
                    "This penalty increases the MAPQ-related systematic error."),
    BcfFormatStruct("bNMa"  , BCF_NUM_R, BCF_INTEGER, "Percent of the number of fragment positions that are affected by nearby high-BQ mutation "
                    "(by default, nearby means within one turn of DNA helix) for the ALT alleles."),
    BcfFormatStruct("bNMb"  , BCF_NUM_R, BCF_INTEGER, "Percent of the number of fragment positions that are affected by nearby high-BQ mutation "
                    "(by default, nearby means within one turn of DNA helix) for all non-ALT alleles."),
    BcfFormatStruct("bMQQ"  , BCF_NUM_R, BCF_SIG_INT, "Variant quality capped by systematic error inferred from mapping qualities."),

    BcfFormatStruct("__e3"  , 1,         BCF_SEP,     "Quality-related variables assuming IID read support (IID: independent and identically distributed) for duped reads."),
    BcfFormatStruct("bIAQb" , BCF_NUM_R, BCF_SIG_INT, "IID allele quality maximized with IAD and IDQ."),
    BcfFormatStruct("bIADb" , BCF_NUM_R, BCF_INTEGER, "IID allele depth (number of reads supporting each ALT)"),
    BcfFormatStruct("bIDQb" , BCF_NUM_R, BCF_SIG_INT, "IDD quality threshold per read support."),
    
    BcfFormatStruct("__e4"  , 1,         BCF_SEP,     "Quality-related variables assuming read supports are IID for deduped reads."),
    BcfFormatStruct("cIAQf" , BCF_NUM_R, BCF_SIG_INT, "IID allele quality maximized with IAD and IDQ on the forward read orientation."),
    BcfFormatStruct("cIADf" , BCF_NUM_R, BCF_INTEGER, "IID allele depth (number of reads supporting each ALT) on the forward read orientation."),
    BcfFormatStruct("cIDQf" , BCF_NUM_R, BCF_SIG_INT, "IDD quality threshold per read support on the forward read orientation."),
    
    BcfFormatStruct("__e5"  , 1,         BCF_SEP,     "Quality-related variables assuming read supports are IID for deduped reads."),
    BcfFormatStruct("cIAQr" , BCF_NUM_R, BCF_SIG_INT, "IID allele quality maximized with IAD and IDQ on the reverse read orientation."),
    BcfFormatStruct("cIADr" , BCF_NUM_R, BCF_INTEGER, "IID allele depth (number of reads supporting each ALT) on the reverse read orientation."),
    BcfFormatStruct("cIDQr" , BCF_NUM_R, BCF_SIG_INT, "IDD quality threshold per read support on the reverse read orientation."),
    
    BcfFormatStruct("__e6"  , 1,         BCF_SEP,     "Binomial variant qualities."),
    BcfFormatStruct("bIAQ"  , BCF_NUM_R, BCF_SIG_INT, "The duped fragment binomial variant quality by assuming statistical independence."),
    BcfFormatStruct("cIAQ"  , BCF_NUM_R, BCF_SIG_INT, "The deduplicated fragment binomial variant quality by assuming statistical independence."),
    BcfFormatStruct("bTINQ" , BCF_NUM_R, BCF_SIG_INT, "The fragment binomial tumor-in-normal quality by assuming statistical independence."),
    BcfFormatStruct("cTINQ" , BCF_NUM_R, BCF_SIG_INT, "The single-strand-consensus-sequence binomial tumor-in-normal quality by assuming statistical independence."),
 
    BcfFormatStruct("__e7"  , 1,         BCF_SEP,     "Power-law variant quality statistics for deduped read fragments."),
    BcfFormatStruct("cPCQ1" , BCF_NUM_R, BCF_SIG_INT, "The deduplicatd read fragment power-law variant allele quality cap that is not lowered by using matched normal."),
    BcfFormatStruct("cPLQ1" , BCF_NUM_R, BCF_SIG_INT, "The deduplicatd read fragment power-law variant allele quality."),
    BcfFormatStruct("cVQ1"  , BCF_NUM_R, BCF_SIG_INT, "The final variant quality computed with deduplicated read fragments."),
    BcfFormatStruct("gVQ1"  , BCF_NUM_R, BCF_SIG_INT, "The final variant quality computed with deduplicated read fragments used for germline variant calls."),

    BcfFormatStruct("__e8"  , 1,         BCF_SEP,     "Power-law variant quality statistics for SSCS UMI-barcoded families."),
    BcfFormatStruct("cPCQ2" , BCF_NUM_R, BCF_SIG_INT, "The SSCS power-law variant allele quality cap that is not lowered by using matched normal."),
    BcfFormatStruct("cPLQ2" , BCF_NUM_R, BCF_SIG_INT, "The single-strand-consensus-sequence (SSCS) UMI-barcoded power-law variant allele quality."),
    BcfFormatStruct("cVQ2"  , BCF_NUM_R, BCF_SIG_INT, "The final variant allele quality computed with SSCS UMI-barcoded families."),
    BcfFormatStruct("cMmQ"  , BCF_NUM_R, BCF_INTEGER, "The empirical base quality estimated with high-quality UMI barcode families. "
                                                      "This field is defined to be the Phred-scale fraction of minor read support "
                                                      "to the sum of major and major read support, "
                                                       "where major means agreement with UMI consensus, and minor means disagreement with UMI consensus. "),
    BcfFormatStruct("dVQinc", BCF_NUM_R, BCF_SIG_INT, "The increase in cVQ2 (excluding systematic error) contributed by "
                                                      "double-strand-consensus-sequences (DSCSs) of duplex barcode famillies. "
                                                      "Negative value means no increase."),

    BcfFormatStruct("__Ef1" , 1,         BCF_SEP,     "Effective read support for deduplicated read fragments for all alleles and the padded-deletion allele."),
    BcfFormatStruct("CDP1v" , 2        , BCF_INTEGER, "The effective number of deduplicated read fragments supporting all alleles multiplied by 100 "
                                                      "for within-sample comparison."),
    BcfFormatStruct("CDP1w" , 2        , BCF_INTEGER, "The effective number of deduplicated read fragments supporting all alleles multiplied by 100 "
                                                      "for sample-specific variant-quality cap."), 
    BcfFormatStruct("CDP1x" , 2        , BCF_INTEGER, "The effective number of deduplicated read fragments supporting all alleles multiplied by 100 "),

    BcfFormatStruct("__ef1" , 1,         BCF_SEP,     "Effective read support for deduplicated read fragments for the REF allele and each ALT allele."),
    BcfFormatStruct("cDP1v" , BCF_NUM_R, BCF_INTEGER, "The effective number of deduplicated read fragments supporting each allele multiplied by 100 "
                                                      "for within-sample comparison."),
    BcfFormatStruct("cDP1w" , BCF_NUM_R, BCF_INTEGER, "The effective number of deduplicated read fragments supporting each allele multiplied by 100 "
                                                      "for sample-specific variant-quality cap."),
    BcfFormatStruct("cDP1x" , BCF_NUM_R, BCF_INTEGER, "The effective number of deduplicated read fragments supporting each allele multiplied by 100 "
                                                      "for between-sample comparison."),
    
    BcfFormatStruct("__Ef2" , 1,         BCF_SEP,     "Effective read support for SSCS UMI-barcoded families."),
    BcfFormatStruct("CDP2v" , 2        , BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting all alleles multiplied by 100 "
                                                      "for within-sample comparison."),
    BcfFormatStruct("CDP2w" , 2        , BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting all alleles multiplied by 100 "
                                                      "for sample-specific variant-quality cap."), 
    BcfFormatStruct("CDP2x" , 2        , BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting all alleles multiplied by 100 "
                                                      "for between-sample comparison."),

    BcfFormatStruct("__ef2" , 1,         BCF_SEP,     "Effective read support for SSCS UMI-barcoded families."),
    BcfFormatStruct("cDP2v" , BCF_NUM_R, BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting each allele multiplied by 100 "
                                                      "for within-sample comparison."),
    BcfFormatStruct("cDP2w" , BCF_NUM_R, BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting each allele multiplied by 100 "
                                                      "for sample-specific variant-quality cap."),
    BcfFormatStruct("cDP2x" , BCF_NUM_R, BCF_INTEGER, "The effective number of SSCS UMI-barcoded families supporting each allele multiplied by 100 "
                                                      "for between-sample comparison."),
        
    BcfFormatStruct("__f1"  , 1,         BCF_SEP,     "Filter-related informationn including counter-filters that rescue variants."),
    BcfFormatStruct("CONTQ" , BCF_NUM_R, BCF_SIG_INT, "Likelihood of the variant signal if the variant signal is contaminated. "
                                                      "This value rescues variants in the matched normal control sample."),
    BcfFormatStruct("nPF"   , BCF_NUM_D, BCF_SIG_INT, "Phred prior bias probability for base position and BAQ."),
    BcfFormatStruct("nNFA"  , BCF_NUM_D, BCF_SIG_INT, "DeciPhred allele fractions computed using nullified bias (meaning that bias in ALT is countered by bias in REF). "),
    BcfFormatStruct("nAFA"  , BCF_NUM_D, BCF_SIG_INT, "DeciPhred allele fractions computed with sequencing-segment depths, "
            "reduced by none, left base position, base right position, left BAQ, right BAQ, left insert positoin, right insert position, strand, and pasing-filter biases (8 biases), respectively. "),
    BcfFormatStruct("nBCFA",  BCF_NUM_D, BCF_SIG_INT, "DeciPhred allele fractions computed with duped and deduped depths, "
            "using duped depth, deduped depth, tier-2 consensus family depth, tier-3 consensus family depth, tier-1 read-orientation depth, and tier-2 read-orientation depth, respectively. "),
    
    BcfFormatStruct("__g1"  , 1,         BCF_SEP,     "General variant-related information."),
    BcfFormatStruct("VTI"   , BCF_NUM_R, BCF_INTEGER, "Variant-type ID of each allele."),
    BcfFormatStruct("VTD"   , BCF_NUM_R, BCF_STRING,  "Variant-type description of each allele"),
    BcfFormatStruct("cVQ1M" , 2,         BCF_SIG_INT, "Consensus allele qualities for deduped fragments and UMI families"),
    BcfFormatStruct("cVQ2M" , 2,         BCF_SIG_INT, "Consensus allele qualities for deduped fragments and UMI families"),
    BcfFormatStruct("cVQAM" , 2,         BCF_STRING,  "Consensus allele symbolic descriptions for deduped fragments and UMI families"),
    BcfFormatStruct("cVQSM" , 2,         BCF_STRING,  "Consensus allele InDel strings for deduped fragments and UMI families"),
    
    BcfFormatStruct("__g2"   , 1,        BCF_SEP,     "Gap-related information for InDels."), 
    BcfFormatStruct("gapNf"  ,BCF_NUM_D, BCF_INTEGER, "Number of InDel sequences on the forward read orientation."),
    BcfFormatStruct("gapNr"  ,BCF_NUM_D, BCF_INTEGER, "Number of InDel sequences on the reverse read orientation."),
    BcfFormatStruct("gapSeq" ,BCF_NUM_D, BCF_STRING,  "InDel sequences"),
    BcfFormatStruct("gapbAD1",BCF_NUM_D, BCF_INTEGER, "Duped read count of each gapSeq"),
    BcfFormatStruct("gapcAD1",BCF_NUM_D, BCF_INTEGER, "Deduped read count of each gapSeq"),
    
    BcfFormatStruct("bDPa"   ,BCF_NUM_R, BCF_INTEGER, "Number of non-deduplicated fragments supporting each ALT allele which is more precise for InDels."),
    BcfFormatStruct("cDP0a"  ,BCF_NUM_R, BCF_INTEGER, "Number of deduplicated fragments supporting each ALT allele which is more precise for InDels"),
    BcfFormatStruct("gapSa"  ,BCF_NUM_R, BCF_STRING,  "InDel string of each allele"),
    
    BcfFormatStruct("__h1"   , 1,        BCF_SEP,     "Haplotype-related information."), 

    BcfFormatStruct("bHap"  , 1,         BCF_STRING,  "Duped forward&reverse linkage in the format of ((position&variantType)...depth)... "
                                                      "where ()... means more elements following the format in the preceding parenthesis. "),
    BcfFormatStruct("cHap"  , 1,         BCF_STRING,  "Dedup forward&reverse linkage in the format of ((position&variantType)...depth)... "
                                                      "where ()... means more elements following the format in the preceding parenthesis. "),    
    BcfFormatStruct("vAC"   , 2,         BCF_INTEGER,  "Number of SNVs and InDels that passed their variant quality thresholds at this position. This field can be used to filter out multiallelic variants or to merge alleles at the same position."),
    BcfFormatStruct("vNLODQ", 2,         BCF_SIG_INT,  "Number of SNVs and InDels that passed their variant quality thresholds at this position. This field can be used to filter out SNV-in-tumor with InDel-in-normal or InDel-in-tumor with SNV-in-normal at the same position."),
    

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
    
    std::cout << "static int resetBcfFormatD(BcfFormat & fmt) {\n";
    
    for (auto fmt : FORMAT_VEC) {
        if (BCF_NUM_D == fmt.number) {
            std::cout << "    fmt." << fmt.id << ".clear();\n";
        }
    }
    std::cout << "\n    return 0;};\n";
    
    std::cout << "static int streamFrontPushBcfFormatR(BcfFormat & dst, const BcfFormat & src) {\n";
    
    for (auto fmt : FORMAT_VEC) {
        if (BCF_NUM_R == fmt.number) {
            std::cout << "    assert(dst." << fmt.id << ".size() == 1 || !fprintf(stderr, \"\%lu == 1 failed for " << fmt.id 
                    << "\", dst." << fmt.id << ".size() ) );\n";
            std::cout << "    assert(src." << fmt.id << ".size() == 1 || !fprintf(stderr, \"\%lu == 1 failed for " << fmt.id 
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

