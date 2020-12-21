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

// constants

#define MGVCF_REGION_MAX_SIZE 1000

#define OUTVAR_GERMLINE 0x1
#define OUTVAR_SOMATIC 0x2
#define OUTVAR_ANY 0x4
#define OUTVAR_MGVCF 0x8
#define OUTVAR_LONG_CLIP 0x10
#define OUTVAR_BASE_NN 0x20
#define OUTVAR_LINK_NN 0x40

#define NOT_PROVIDED ("")
#define OPT_ONLY_PRINT_VCF_HEADER "/only-print-vcf-header/"
#define PLAT_ILLUMINA_LIKE "Illumina/BGI"
#define PLAT_ION_LIKE "IonTorrent/LifeTechnologies/ThermoFisher"

#define MAX_INSERT_SIZE 2000 // (1024*2) // https://doi.org/10.2147/AGG.S162531
#define DBLFLT_EPS ((double)FLT_EPSILON)

// substitutions

#define SIGN2UNSIGN(x) ((x)) // disabled
#define UNSIGN2SIGN(x) ((const int64_t)(x))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define INT64MUL(x, y) ((int64_t)(x) * (int64_t)(y))

#define NORM_INSERT_SIZE(b) { if (abs((b)->core.isize) >= MAX_INSERT_SIZE) { (b)->core.isize = 0; } }

#define rtr_endpos(r) ((r).begpos + (r).tracklen)
#define ispowerof2(num) (((num) & ((num)-1)) == 0)

// conversion between Phred, nat, bit, frac, and states

#define  phred2nat(x) ((log(10.0)/10.0) * (x))
#define  nat2phred(x) ((10.0/log(10.0)) * (x))
#define frac2phred(x) (-(10.0/log(10.0)) * log(x))
#define phred2frac(x) (pow(10.0, (-x)/10.0))
#define numstates2phred(x) ((10.0/log(10.0)) * log(x))
#define phred2numstates(x) (pow(10.0, (x)/10.0))
#define numstates2deciphred(x) ((uvc1_qual_t)round((100.0/log(10.0)) * log(x)))

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

struct RegionalTandemRepeat {
    uvc1_refgpos_t begpos = 0;
    uvc1_readpos_t tracklen = 0;
    uvc1_readpos_t unitlen = 0;

    uvc1_qual_t indelphred = 40 + 3;

    uvc1_refgpos_t anyTR_begpos = 0;
    uvc1_readpos_t anyTR_tracklen = 0;
    uvc1_readpos_t anyTR_unitlen = 0;
};

#endif

