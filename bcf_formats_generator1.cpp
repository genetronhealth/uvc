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
    std::make_pair("cad3",          "Less than 3 clean deduppd reads"),
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

const std::vector<BcfFormatStruct> FORMAT_VEC = {
    BcfFormatStruct("GT"       , 1, BCF_STRING,  "Genotype (this is a guess for the GT of tumor cells if --somaticGT was set to true)"),
    BcfFormatStruct("GQ"       , 1, BCF_INTEGER, "Genotype Quality"),
    BcfFormatStruct("HQ"       , 2, BCF_INTEGER, "Haplotype Quality"),
    BcfFormatStruct("DP"       , 1, BCF_INTEGER, "Fragment depth supporting any allele [consensus family, deduped]"),
    BcfFormatStruct("FA"       , 1, BCF_FLOAT,   "Frequency of the ALT allele [consensus family, deduped]"),
    BcfFormatStruct("FR"       , 1, BCF_FLOAT,   "Frequency of the REF allele [consensus family, deduped]"),
    BcfFormatStruct("BQ"       , 1, BCF_INTEGER, "Root mean square (RMS) base quality of the ALT [base read, duped]"), 
    BcfFormatStruct("MQ"       , 1, BCF_INTEGER, "Root mean square (RMS) mapping quality of the ALT [base read, duped]"), 
    BcfFormatStruct("FT"       , 1, BCF_STRING,  "Sample genotype filter indicating if this genotype was 'called' (similar in concept to the FILTER field). "
                                                 "Again, use PASS to indicate that all filters have been passed, a semi-colon separated list of codes for filters "
                                                 "that fail, or ‘.’ to indicate that filters have not been applied. "
                                                 "These values should be described in the meta-information in the same way as FILTERs "
                                                 "(String, no white-space or semi-colons permitted)"),
    
    BcfFormatStruct("__A1"     , 1, BCF_SEP,     "Dummy header (separator) for genotype-related FORMAT fields. The description of each associated FORMAT is shown below. "),
    BcfFormatStruct("FTS"      , 1, BCF_STRING,  "Sample variant filter indicating if this variant was 'called' (similar in concept to the FILTER field). "
                                                 "Again, use PASS to indicate that all filters have been passed, "
                                                 "an amperstand-separated list of codes for filters that fail, "
                                                 "or '.' to indicate that filters have not been applied. "
                                                 "These values should be described in the meta-information in the same way as FILTERs. "
                                                 "No white-space, semi-colons, or amperstand permitted."),
    BcfFormatStruct("FTSV"     , BCF_NUM_D, BCF_INTEGER, "Percent bias values for the FTS strings"),
    BcfFormatStruct("DPHQ"     , 1, BCF_INTEGER, "Fragment depth supporting any allele if low-quality bases are ignored which means only high quality (HQ) bases are used [base read, duped]"),
    BcfFormatStruct("ADHQ"     , 1, BCF_INTEGER, "Fragment depth supporting the ALT allele if low-quality bases are ignored which means only high quality (HQ) bases are used [base read, duped]"),
    BcfFormatStruct("ALODQ"    , 1, BCF_INTEGER, "Artifact log-likelihood of data quality in PHRED-scale."),
    BcfFormatStruct("BLODQ"    , 1, BCF_INTEGER, "Bias-of-allele log-likelihood of data quality in PHRED-scale."),
    BcfFormatStruct("EROR"     , 5, BCF_INTEGER, "Deduplication bias, position bias, strand bias, and mismatch bias of the evidence reduction odds ratio for deduplicated reads, and the maximum of all biases in FTS"), 
 
    BcfFormatStruct("__A2"     , 1, BCF_SEP,     "Dummy header (separator) for genotype-related FORMAT fields. The description of each associated FORMAT is shown below. "),
    BcfFormatStruct("RefBias"  , 1, BCF_INTEGER, "Reference bias (on average, the read supporting the ALT is shorter than the read supporting the REF by this many bases)"),
    BcfFormatStruct("GTa"      , 1, BCF_STRING,  "Genotype of this alelle"),
    BcfFormatStruct("GQa"      , 1, BCF_INTEGER, "Genotype Quality of this allele"),
    BcfFormatStruct("GLa"      , 3, BCF_SIG_INT, "Genotype Likelihood of this allele"),
    BcfFormatStruct("GSTa"     ,10, BCF_SIG_INT, "Genotype Statistics of this allele"),
    BcfFormatStruct("GTb"      , 1, BCF_STRING,  "Genotype of other alelles"),
    BcfFormatStruct("GQb"      , 1, BCF_INTEGER, "Genotype Quality of other alleles"),
    BcfFormatStruct("GLb"      , 3, BCF_SIG_INT, "Genotype Likelihood of other alleles"),
    BcfFormatStruct("GSTb"     ,10, BCF_SIG_INT, "Genotype Statistics of other alleles"),
    
    BcfFormatStruct("__A3"     , 1, BCF_SEP,     "Depth and REF/ALT allele frequency for base read and consensus family"),  
    BcfFormatStruct("bDP"      , 1, BCF_INTEGER, "Fragment depth supporting any allele [base read, duped]"),
    BcfFormatStruct("bADR"     , BCF_NUM_R, BCF_INTEGER, "Fragment depth supporting each REF and ALT allele [base read, duped]"),
    BcfFormatStruct("cDP"      , 1, BCF_INTEGER, "Fragment depth supporting any allele [consensus family, deduped]"),
    BcfFormatStruct("cADR"     , BCF_NUM_R, BCF_INTEGER, "Fragment depth supporting each REF and ALT allele [consensus family, deduped]"),
    
    BcfFormatStruct("__A4"     , 1, BCF_SEP,     "Consensus/variant allele type/quality"),  
    BcfFormatStruct("OType"    , 1, BCF_STRING,  "The non-reference allele type with the most evidence other than the ALT allele type of this record"),
    BcfFormatStruct("ORAQs"    , 2, BCF_FLOAT,   "Qualities of OType and reference allele type"),
    BcfFormatStruct("VType"    , 1, BCF_STRING,  "Variant type for the ALT allele"),
    // BcfFormatStruct("VAQs"     , 2, BCF_FLOAT,   "Raw variant allele quality (VAQ) and VAQ of the specific form(s) of InDel in ALT assuming other forms of InDels are noise"),
    // BcfFormatStruct("VAQAB"    , 1, BCF_FLOAT,   "Variant Allele Quality adjusted with bias"),
    BcfFormatStruct("VQ1"      , BCF_NUM_A, BCF_INTEGER, "Variant allele quality capped by base alignment quality"),
    BcfFormatStruct("VQ2"      , BCF_NUM_A, BCF_INTEGER, "Variant allele quality capped by base quality"),
    BcfFormatStruct("VQ3"      , BCF_NUM_A, BCF_INTEGER, "Variant allele quality capped by base quality"),
    BcfFormatStruct("VAQ"      , BCF_NUM_A, BCF_INTEGER, "Variant allele quality of the call"),
    
    BcfFormatStruct("__A5"     , 1, BCF_SEP,     "Sum of base qualities (bAllBQ/bAltBQ/bRefBQ) and consensus-family qualities (cAllBQ/cAltBQ/cRefBQ) for ALL/ALT/REF alleles"), // global 
    BcfFormatStruct("bAllBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bAltBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bRefBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAllBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAltBQ"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRefBQ"   , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__A6"     , 1, BCF_SEP,     "Sum of square/10 consensus-family qualities and number of high-confidence families, and and duplication bias (aDB). All biases are percentages where 100% means no bias."), // global
    //BcfFormatStruct("cAllBQ2"  , 2, BCF_INTEGER, "see above"),
    //BcfFormatStruct("cAltBQ2"  , 2, BCF_INTEGER, "see above"),
    //BcfFormatStruct("cRefBQ2"  , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAllHD"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAltHD"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRefHD"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("aDB"      , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__aa"     , 1, BCF_SEP,     "Sequencing-segment statistics. In aDP and aAD, the four-directions are (fw-fw, fw-rv, rv-fw, rv-rv) strand orientations corresponding to (R1+, R2-, R2+, R1-)"),
    BcfFormatStruct("aDP"      , 4, BCF_INTEGER, "Four-direction raw depths total and summed for ALL alleles [unmerged]"),
    BcfFormatStruct("aRD"      , 4, BCF_INTEGER, "Four-direction raw depths specific to only the ALT allele  [unmerged]"),
    BcfFormatStruct("aAD"      , 4, BCF_INTEGER, "Four-direction raw depths specific to only the ALT allele  [unmerged]"),
    
    BcfFormatStruct("__ab"     , 1, BCF_SEP,     "More sequencing-segment statistics."),
    BcfFormatStruct("aNMRD"    , 4, BCF_INTEGER, "Four-direction total numbers of mismatches (the NM sam tag) specific to only the REF allele  [unmerged"),
    BcfFormatStruct("aNMAD"    , 4, BCF_INTEGER, "Four-direction total numbers of mismatches (the NM sam tag) specific to only the ALT allele  [unmerged"),
    BcfFormatStruct("aBQAD"    , 4, BCF_INTEGER, "Four-direction total sum of base qualities of the ALT allele  [unmerged"),
    
    BcfFormatStruct("__ac"     , 1, BCF_SEP,     "More sequencing-segment statistics."),
    BcfFormatStruct("aPBDP"    , 3, BCF_INTEGER, "Raw depths within the left&right, left, and right regions for all the alleles [unmerged]"),
    BcfFormatStruct("aPBAD"    , 3, BCF_INTEGER, "Raw depths within the left&right, left, and right regions for the ALT allele  [unmerged]"),
    BcfFormatStruct("aB"       , 4, BCF_INTEGER, "For sequencing-segments: strand, left&right position, left position, and right position biases [unmerged]"),
    BcfFormatStruct("aBAQDP"   , 1, BCF_INTEGER, "Total summed base alignment quality (BAQ) of all alleles [unmerged]"),
    BcfFormatStruct("aBAQADR"  , BCF_NUM_R, BCF_INTEGER, "Total summed base alignment quality (BAQ) of each REF and ALT allele [unmerged]"),

    BcfFormatStruct("__ba"     , 1, BCF_SEP,     "Forward&reverse max-bias base distances to left/right end positions (T1PTL/T1PTR) and position bias (T1PBL/T1PBR) [base read, duped]"), 
    BcfFormatStruct("bPTL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bPTR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bPBL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bPBR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bMMT"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bMMB"     , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__bb"     , 1, BCF_SEP,     "Forward&reverse sum of : base distances to left/right end positions (T1SDL/T1SDR) and non-normalized/normalized strand bias (T1SB1/T1SBR) [base read, duped]"),
    BcfFormatStruct("bSDL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bSDR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bSB1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bSBR"     , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__bc"     , 1, BCF_SEP,     "Forward&reverse depth (DP1), allele depth (AD1), allele depth after weak filter (AD2), quality threshold for weak filter (QT2) [base read, duped]"),
    BcfFormatStruct("bDP1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bAD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bRD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bAD2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bQT2"     , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__bd"     , 1, BCF_SEP,     "Forward&reverse allele depth after strong filter (AD3), quality threshold for strong filter (QT3), and variant quality after strong filter (VQ3) [duped]"),
    BcfFormatStruct("bAD3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bADB"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bQT3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bVQ3"     , 2, BCF_FLOAT,   "see above"),
    
    BcfFormatStruct("__be"     , 1, BCF_SEP,     "Forward&reverse  allele root-mean-squre MQ (bMQ1) and (sum of squared MQ) / (sum of MQ) (bMQ2), where MQ is mapping quality, and the same for base quality BQ [duped]"),
    BcfFormatStruct("bMQ1"     , 2, BCF_INTEGER, "see above"),
    //BcfFormatStruct("bMQ2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bBQ1"     , 2, BCF_INTEGER, "see above"),
    //BcfFormatStruct("bBQ2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bNSB"     , 2, BCF_INTEGER, "Total number of sequenced bases on all reads covering this position [base read, duped]"), 
    // BcfFormatStruct("bEDRD"    , 2, BCF_INTEGER, "Summed Phred-scaled edit distance of reads supporting the REF where an InDel corresponds to 2 edits [base read, duped]"),
    // BcfFormatStruct("bEDAD"    , 2, BCF_INTEGER, "Summed Phred-scaled edit distance of reads supporting the ALT where an InDel corresponds to 2 edits [base read, duped]"),

    BcfFormatStruct("__bf"     , 1, BCF_SEP,     "Forward&reverse total depth (bDPLQ), ALT depth (bADLQ), and REF depth (bRDLQ) supported by low-quality (LQ) bases [duped]"),
    BcfFormatStruct("bDPLQ"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bADLQ"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("bRDLQ"    , 2, BCF_INTEGER, "see above"),
    
      
    BcfFormatStruct("__ca"     , 1, BCF_SEP,     "Same as __ba but for single-strand families instead of reads [consensus family, deduped]"),
    BcfFormatStruct("cPTL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cPTR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cPBL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cPBR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cMMT"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cMMB"     , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__cb"     , 1, BCF_SEP,     "Same as __bb but for single-strand families instead of reads [consensus family, deduped]"), 
    BcfFormatStruct("cSDL"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cSDR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cSB1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cSBR"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cNSB"     , 2, BCF_INTEGER, "Total number of sequenced bases on all reads covering this position [consensus family, deduped]"),

    BcfFormatStruct("__cc"     , 1, BCF_SEP,     "Same as __bc but for single-strand families instead of reads [consensus family, deduped]"), 
    BcfFormatStruct("cDP1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRD1"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cAD2"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cQT2"     , 2, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__cd"     , 1, BCF_SEP,     "Same as __bd but for single-strand families instead of reads [consensus family, deduped]"),
    BcfFormatStruct("cAD3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cADB"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cQT3"     , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cVQ3"     , 2, BCF_FLOAT,   "see above"),
    
    BcfFormatStruct("__ce"     , 1, BCF_SEP,     "number of read supports generated by signal (cMajor) and noise (cMinor) inferred from heterogeneity of families, "
                                                 "and strand-specific VAQ (non-adjusted: cVAQ1; adjusted with BQ, MQ, and strand balance: cVAQ2) [consensus family, deduped]"),
    BcfFormatStruct("cMajor"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cMinor"   , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cVAQ1"    , 2, BCF_FLOAT,   "see above"),
    BcfFormatStruct("cVAQ2"    , 2, BCF_FLOAT,   "see above"),
    
    BcfFormatStruct("__cf"     , 1, BCF_SEP,     "Forward&reverse number of families without any filtering (REF/ALT/ALL total: cRDTT/cADTT/cDPTT, "
                                                 "REF/ALT/ALL allele singleton: cRDT1/cADT1/cDPT1, REF/ALT/ALL allele failing 80\% consensus : cRDTN/cADTN/cDPTN). "),
    BcfFormatStruct("cADTT"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cADT1"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cADTC"    , 2, BCF_INTEGER, "(TT - T1 - TN) for ALT"),
    BcfFormatStruct("cRDTT"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRDT1"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cRDTC"    , 2, BCF_INTEGER, "(TT - T1 - TN) for REF"),
    BcfFormatStruct("cDPTT"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cDPT1"    , 2, BCF_INTEGER, "see above"),
    BcfFormatStruct("cDPTC"    , 2, BCF_INTEGER, "(TT - T1 - TN) for ALL"),
    
    BcfFormatStruct("__da"     , 1, BCF_SEP,     "Same as __ba and __ca but for duplex families [duplex family, DSCS1-deduped]"),
    BcfFormatStruct("dDP1"     , 1, BCF_INTEGER, "see above"),
    BcfFormatStruct("dAD1"     , 1, BCF_INTEGER, "see above"),
    BcfFormatStruct("dVQ3"     , 1, BCF_INTEGER, "see above"),
    BcfFormatStruct("dAD3"     , 1, BCF_INTEGER, "see above"),
    
    BcfFormatStruct("__ea"    , 1, BCF_SEP,     "Number of indels in forward&reverse strand (gapNum), max cAD diff (gapcADD), total cAD depth (gapcADT), indel sequences (gapSeq), duped read count of each gapSeq, dedup family count of each gapSeq, and duped/deduped sub/all allele read counts (gapDP4)"), 
    BcfFormatStruct("gapNum"   , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    BcfFormatStruct("gapcADD"  , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    BcfFormatStruct("gapcADT"  , 2, BCF_INTEGER, "see above"), // 2 * number-of-alts
    BcfFormatStruct("gapSeq"   , BCF_NUM_D, BCF_STRING,  "see above"),
    BcfFormatStruct("gapbAD1"  , BCF_NUM_D, BCF_INTEGER, "see above"),
    BcfFormatStruct("gapcAD1"  , BCF_NUM_D, BCF_INTEGER, "see above"),
    BcfFormatStruct("gapDP4"   , 4, BCF_INTEGER, "see above"),
    BcfFormatStruct("gapbNRD"  , 4, BCF_INTEGER, "Forward&reverse (02&13) strand-specific duped read count for any insertion&deletion (01&23)"),
    // BcfFormatStruct("gapbNNRD" , 2, BCF_INTEGER, "Highest duped read count for any insertion&deletion nearby"),
    BcfFormatStruct("RCC"      , RCC_NFS*RCC_NUM, BCF_INTEGER, 
                                "STR-unit position of the mode, two indel counts of -2 and -1 STR units, mode count, two ins counts of +1 and +2 STR units with respect to the mode"),
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
            std::cout << "    outstring += std::string(FORMAT_IDS[" << itnum << "]) + \"\";\n";
        } else if (0 == fmt.number || 1 == fmt.number) {
            std::cout << "    outstring += " << (fmt.type == BCF_STRING ? "" : "std::to_string") << "(fmt." << fmt.id << ");\n";
        } else {
            std::cout << "    for (unsigned int i = 0; i < fmt." << fmt.id << ".size(); i++) {\n";
            std::cout << "        if (0 != i) { outstring += \",\"; }; outstring += " << (fmt.type == BCF_STRING ? "" : "std::to_string") 
                    << "(" << "fmt." << fmt.id << "[i]" << ");\n";
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

