#ifndef common_hpp_INCLUDED
#define common_hpp_INCLUDED

#define COMPILATION_ENABLE_XMGOT 0
#define COMPILATION_TRY_HIGH_DEPTH_POS_BIAS 0

// #include "precompiled/precompiled_main.hpp"

#include <array>
#include <string>
#include <vector>

#include <float.h>
#include <limits.h>

#ifdef __GNUC__
#if __GNUC__ > 3
#define IGNORE_UNUSED_PARAM __attribute__((unused))
#else
#define IGNORE_UNUSED_PARAM // let warning pop up but still compile fine
#endif
#else
#define IGNORE_UNUSED_PARAM 
#endif

#define SIGN2UNSIGN(x) ((x)) // disabled
#define UNSIGN2SIGN(x) ((const int64_t)(x))

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// conversion between Phred, nat, bit, frac, and states
#define  phred2nat(x) ((log(10.0)/10.0) * (x))
#define  nat2phred(x) ((10.0/log(10.0)) * (x))
#define frac2phred(x) (-(10.0/log(10.0)) * log(x))
#define phred2frac(x) (pow(10.0, (-x)/10.0))
#define numstates2phred(x) ((10.0/log(10.0)) * log(x))
#define phred2numstates(x) (pow(10.0, (x)/10.0))

#define numstates2deciphred(x) ((uvc1_qual_t)round((100.0/log(10.0)) * log(x)))

#define OUTVAR_GERMLINE 0x1
#define OUTVAR_SOMATIC 0x2
#define OUTVAR_ANY 0x4
#define OUTVAR_MGVCF 0x8
#define OUTVAR_BASE_NN 0x10
#define OUTVAR_LINK_NN 0x20

#define OPT_ONLY_PRINT_VCF_HEADER "/only-print-vcf-header/"
#define PLAT_ILLUMINA_LIKE "Illumina/BGI"
#define PLAT_ION_LIKE "IonTorrent/LifeTechnologies/ThermoFisher"

#define MAX_INSERT_SIZE 2000 // (1024*2) // https://doi.org/10.2147/AGG.S162531
#define NORM_INSERT_SIZE(b) { if (abs((b)->core.isize) >= MAX_INSERT_SIZE) { (b)->core.isize = 0; } }
#define DBLFLT_EPS ((double)FLT_EPSILON)

#define INT64MUL(x, y) ((int64_t)(x) * (int64_t)(y))

typedef uint64_t uvc1_unsigned_int_t; // It seems that a bug in g++ causes compiling error if this type is defined as (unsigned int)

typedef int32_t uvc1_qual_t;    // quality (usually Phred-scaled)
typedef int32_t uvc1_deciphred_t; // 10x Phred

typedef int32_t uvc1_readnum_t; // depth of segment, fragment, family, etc. // max 2 billion reads
typedef int32_t uvc1_readnum100x_t; // 100x depth of segment, fragment, etc. // max 20 million reads
typedef int32_t uvc1_readpos_t; // position with respect to the read, fragment, or insert
typedef int32_t uvc1_refgpos_t; // position with respect to the reference genome
typedef int32_t uvc1_rp_diff_t;

typedef int32_t uvc1_base_t;
typedef int32_t uvc1_base1500x_t;

typedef int64_t uvc1_readnum_big_t;
typedef int64_t uvc1_readpos_big_t;
typedef int64_t uvc1_qual_big_t; // big qual

typedef uint32_t uvc1_flag_t;
typedef uint64_t uvc1_hash_t;

enum BiasType {
    BIAS_FRAG_DUP = 1,
    BIAS_FRAG_POS = 2,
    BIAS_FRAG_STR = 4,
    BIAS_FRAG_MIS = 8,
    BIAS_SSEG_POS = 32,
    BIAS_SSEG_STR = 64,
    BIAS_SSEG_END = 128,
};

extern const std::vector<std::string> BIAS_TYPE_TO_MSG;

// const char *const NOT_PROVIDED = "";
#define NOT_PROVIDED ("")

enum AssayType {
    ASSAY_TYPE_AUTO,
    ASSAY_TYPE_CAPTURE,
    ASSAY_TYPE_AMPLICON,
};

extern const std::vector<std::string> ASSAY_TYPE_TO_MSG;

enum MoleculeTag {
    MOLECULE_TAG_AUTO,
    MOLECULE_TAG_NONE,
    MOLECULE_TAG_BARCODING,
    MOLECULE_TAG_DUPLEX,
};

extern const std::vector<std::string> MOLECULE_TAG_TO_MSG;

enum SequencingPlatform {
    SEQUENCING_PLATFORM_AUTO,
    SEQUENCING_PLATFORM_ILLUMINA,
    SEQUENCING_PLATFORM_IONTORRENT,
    SEQUENCING_PLATFORM_OTHER,
};

extern const std::vector<std::string> SEQUENCING_PLATFORM_TO_MSG;

extern const std::vector<std::string> SEQUENCING_PLATFORM_TO_NAME;

enum PairEndMerge {
    PAIR_END_MERGE_YES,
    PAIR_END_MERGE_NO,
};

extern const std::vector<std::string> PAIR_END_MERGE_TO_MSG;

extern const std::array<std::string, 2> GT_HETERO;
extern const std::array<std::string, 2> GT_HOMREF;
extern const std::array<std::string, 2> GT_HOMALT;

extern const std::array<std::string, 2> TT_HETERO;
extern const std::array<std::string, 2> TT_HOMREF;
extern const std::array<std::string, 2> TT_HOMALT;

/*
struct TnDP4 {
    unsigned int nvars= 0;
    unsigned int tuAD = 0;
    unsigned int noAD = 0;
    unsigned int tuDP = 0;
    unsigned int noDP = 0;
    TnDP4() {}
};

struct VcStats {
    std::array<TnDP4, VCFQUAL_NUM_BINS> vcfqual_to_count = {{TnDP4()}};
    
    VcStats() {}
    
    void 
    update(const VcStats & other) {
        for (unsigned int i = 0; i < VCFQUAL_NUM_BINS; i++) { 
            this->vcfqual_to_count[i].nvars+= other.vcfqual_to_count[i].nvars;
            this->vcfqual_to_count[i].tuAD += other.vcfqual_to_count[i].tuAD;
            this->vcfqual_to_count[i].noAD += other.vcfqual_to_count[i].noAD;
            this->vcfqual_to_count[i].tuDP += other.vcfqual_to_count[i].tuDP;
            this->vcfqual_to_count[i].noDP += other.vcfqual_to_count[i].noDP;
        }
    };
    
    int
    write_tsv(auto & ostream) {
        ostream <<"##Mass-function-aka-density-histogram\n";
        ostream <<"#Variant-quality\tnumber-of-calls\ttumor-AD\tnormal-AD\ttumor-DP\tnormal-DP\n";
        for (unsigned int i = 0; i < VCFQUAL_NUM_BINS; i++) { 
            ostream << i << "\t"; 
            ostream << this->vcfqual_to_count[i].nvars<< "\t";
            ostream << this->vcfqual_to_count[i].tuAD << "\t";
            ostream << this->vcfqual_to_count[i].noAD << "\t";
            ostream << this->vcfqual_to_count[i].tuDP << "\t";
            ostream << this->vcfqual_to_count[i].noDP << "\n";
        }
        unsigned int tot_nvars = 0;
        unsigned int tot_tuAD = 0;
        unsigned int tot_noAD = 0;
        unsigned int tot_tuDP = 0;
        unsigned int tot_noDP = 0;
        bool vcfqual_is_too_low = false;
        ostream <<"##Complementary-cumulative-function\n";
        ostream <<"#variant-quality\tnumber-of-calls\ttumor-AD\tnormal-AD\ttumor-DP\tnormal-DP"
                << "\testimated-additive-contamination-fraction\testimated-multiplicative-contamination-fraction"
                << "\tTUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION_STATUS_FOR_VCF_QUAL\n";
        for (int i = VCFQUAL_NUM_BINS - 1; i >= 0; i--) { 
            tot_nvars += this->vcfqual_to_count[i].nvars;
            tot_tuAD += this->vcfqual_to_count[i].tuAD;
            tot_noAD += this->vcfqual_to_count[i].noAD;
            tot_tuDP += this->vcfqual_to_count[i].tuDP;
            tot_noDP += this->vcfqual_to_count[i].noDP;
            
            ostream << i << "\t"; 
            ostream << tot_nvars<< "\t";
            ostream << tot_tuAD << "\t";
            ostream << tot_noAD << "\t";
            ostream << tot_tuDP << "\t";
            ostream << tot_noDP << "\t";

            ostream << (0.03 + (double)tot_noAD)                  / (1.0 + (double)tot_tuAD                 ) << "\t"; 
            ostream << (0.03 + (double)tot_noAD*(double)tot_tuDP) / (1.0 + (double)tot_tuAD*(double)tot_noDP) << "\t"; 
            
            if ((!vcfqual_is_too_low) && tot_nvars >= 25 && tot_noDP >= 25*100 && tot_tuDP >= 25*100) {
                ostream << "VCF_QUAL_EQUALS_THE_THRESHOLD_FOR_TUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION" << "\n";  
                vcfqual_is_too_low = true;
            } else if (vcfqual_is_too_low) {
                ostream << "VCF_QUAL_IS_BELOW_THE_THRESHOLD_FOR_TUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION" << "\n"; 
            } else {
                ostream << "VCF_QUAL_IS_ABOVE_THE_THRESHOLD_FOR_TUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION" << "\n";  
            }
        }
        return 0;
    }
};
*/
// region_pos32_unitlen8_repeatnum16_qual8_vec
struct RegionalTandemRepeat {
    uvc1_refgpos_t begpos = 0;
    uvc1_readpos_t tracklen = 0;
    uvc1_readpos_t unitlen = 0;

    uvc1_qual_t indelphred = 40 + 3;

    uvc1_refgpos_t anyTR_begpos = 0;
    uvc1_readpos_t anyTR_tracklen = 0;
    uvc1_readpos_t anyTR_unitlen = 0;
    
    // uint8_t edgeBAQ;
};

#define rtr_endpos(r) ((r).begpos + (r).tracklen)

// I tried to use the following correction types in the past. 
// However, correction types were very confusing to the end user and always missed some edge cases.
// Hence, I used assay types, sequencing platforms, UMI-awareness, etc. instead of correction types.
// If anyone thinks that correction types can be reworked into something intuitive to the end user, please let me know.
/*
enum ErrorCorrectionType {
    CORRECTION_AUTO,
    CORRECTION_NONE,
    CORRECTION_BASEQUAL,
    CORRECTION_SINGLETON,
    CORRECTION_DUPLICATE,
    CORRECTION_BARCODE,
    CORRECTION_DUPLEX,
    END_ERROR_CORRECTION_TYPES
};
*/
/*
bool ispowerof2(auto num) {
    return (num & (num-1)) == 0;
}
*/
#define ispowerof2(num) (((num) & ((num)-1)) == 0)

#endif

