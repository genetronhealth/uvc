#ifndef conversion_hpp_INCLUDED
#define conversion_hpp_INCLUDED

#include "common.hpp"

#include "htslib/sam.h"
#include "htslib/vcf.h"

#include <algorithm>
#include <array>
#include <iostream>

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define SQR_QUAL_DIV 32

// math functions

template <class T>
auto MEDIAN(const T & v) {
    assertUVC(v.size() > 0);
    return (v[(v.size() - 1) / 2] + v[(v.size()) / 2]) / 2;
}

template <class T>
auto FIRST(const T & v) {
    assertUVC(v.size() > 0);
    return v[0];
}

template <class T>
auto LAST(const T & v) {
    assertUVC(v.size() > 0);
    return v[(v.size()-1)];
}

template <class T1, class T2>
int64_t
int64mul(const T1 a, const T2 b) {
    return ((int64_t)a) * ((int64_t)b);
}

template <class T1, class T2>
uint64_t
uint64mul(const T1 a, const T2 b) {
    return ((uint64_t)a) * ((uint64_t)b);
}

template <class T1, class T2>
auto
CENTER(T1 a, T2 b, int center = 0) {
    return ((abs(a - center) < abs(b - center)) ? a : b);
}

template <class T1, class T2, class T3>
auto
MIN3(T1 a, T2 b, T3 c) {
    return MIN(a, MIN(b, c));
}

template <class T1, class T2, class T3, class T4>
auto
MIN4(T1 a, T2 b, T3 c, T4 d) {
    return MIN(a, MIN3(b, c, d));
}

template <class T1, class T2, class T3, class T4, class T5>
auto
MIN5(T1 a, T2 b, T3 c, T4 d, T5 e) {
    return MIN(a, MIN4(b, c, d, e));
}

template <class T>
auto
MINVEC(T v) {
    assertUVC(v.size() > 0);
    auto ret = v[0];
    for (auto e : v) {
        ret = MIN(ret, e);
    }
    return ret;
}

template <class T1, class T2, class T3>
auto
MAX3(T1 a, T2 b, T3 c) {
    return MAX(a, MAX(b, c));
}

template <class T1, class T2, class T3, class T4>
auto
MAX4(T1 a, T2 b, T3 c, T4 d) {
    return MAX(a, MAX3(b, c, d));
}

template <class T1, class T2, class T3, class T4, class T5>
auto
MAX5(T1 a, T2 b, T3 c, T4 d, T5 e) {
    return MAX(a, MAX4(b, c, d, e));
}

template <class T1, class T2, class T3, class T4, class T5, class T6>
auto
MAX6(T1 a, T2 b, T3 c, T4 d, T5 e, T6 f) {
    return MAX(a, MAX5(b, c, d, e, f));
}

template <class T>
auto
MAXVEC(T v) {
    assertUVC(v.size() > 0);
    auto ret = v[0];
    for (auto e : v) {
        ret = MAX(ret, e);
    }
    return ret;
}

template <class T1, class T2, class T3>
auto
BETWEEN(T1 v, T2 a, T3 b) {
    return MIN(MAX(a, v), b);
}

template <class T1, class T2>
void 
UPDATE_MIN(T1 & a, const T2 & b) {
    a = MIN(a, b);
}

template <class T1, class T2>
void 
UPDATE_MAX(T1 & a, const T2 & b) {
    a = MAX(a, b);
}

template <class T>
auto
SUMVEC(const T & vec) {
    auto r = 0;
    for (size_t i = 0; i < vec.size(); i++) {
        r += vec[i];
    }
    return r;
}

template <class T>
T
calc_non_negative(const T v, T base = pow(10.0, 0.1), T thres = 10.0) {
    if (v < thres) {
        return log1p(pow(base, v)) / log(base);
    } else {
        return v;
    }
}

template <class T>
T
calc_score_with_penal_at_low_val(const T v, T penal_mult, T thres = 60.0) {
    return v * penal_mult;
    // enable the following if it makes sense.
    if (v <= thres) {
        return v * penal_mult;
    } else {
        return thres * penal_mult + (v - thres);
    }
}

template <class T>
auto
mathcube(T x) {
    return x * x * x;
}

constexpr
double
prob2odds(double p) {
    assertUVC((0.0 < p && p < 1.0) || !fprintf(stderr, "%f is not between 0 and 1!", p));
    return p /  (1.0 - p);
}

constexpr
double
odds2prob(double odds) {
    assertUVC((0.0 < odds) || !fprintf(stderr, "%f is not greater than zero!", odds));
    return odds / (odds + 1.0);
}

static_assert(prob2odds(odds2prob(1)) > 0.99);
static_assert(prob2odds(odds2prob(1)) < 1.01);

static_assert(odds2prob(prob2odds(0.66)) > 0.65);
static_assert(odds2prob(prob2odds(0.66)) < 0.67);

double
logit(double p) {
    return log(prob2odds(p));
}

double
logit2(double a, double b) {
    return logit((a + DBL_EPSILON)/(a+b + 2.0*DBL_EPSILON));
}

// always at least zero
template <bool TIsBiDirectional = false, bool TSetMaxProbToOne = false>
constexpr double
calc_binom_10log10_likeratio(double prob, double a, double b) {
    if (TSetMaxProbToOne) { prob = MIN(1.0, prob); }
    prob = (prob + DBL_EPSILON) / (1.0 + (2.0 * DBL_EPSILON));
    assertUVC((prob > 0 && prob < 1) || !fprintf(stderr, "The assertUVCion 0 < %f < 1 failed!\n", prob));
    a += DBL_EPSILON;
    b += DBL_EPSILON;
    double A = (      prob) * (a + b);
    double B = (1.0 - prob) * (a + b);
    if (TIsBiDirectional || a > A) {
        return 10.0 / log(10.0) * (a * log(a / A) + b * log(b / B));
    } else {
        return 0.0;
    }
}

#ifdef TEST_calc_binom_10log10_likeratio
int 
main(int argc, char **argv) {
    double prob = atof(argv[1]);
    double a = atof(argv[2]);
    double b = atof(argv[3]);
    double ret1 = calc_binom_10log10_likeratio<true>(prob, a, b);
    double ret2 = calc_binom_10log10_likeratio<false>(prob, a, b);
    printf("calc_binom_10log10_likeratio<false, true>(%f, %f, %f) = <%f, %f>\n", prob, a, b, ret1, ret2);
}
#endif

STATIC_ASSERT_WITH_DEFAULT_MSG(abs(calc_binom_10log10_likeratio(0.1, 10, 90)) < 1e-4);
STATIC_ASSERT_WITH_DEFAULT_MSG(calc_binom_10log10_likeratio(0.1, 90, 10) > 763); // 10/log(10) * (90*log(9)+10*log(1/9))
STATIC_ASSERT_WITH_DEFAULT_MSG(calc_binom_10log10_likeratio(0.1, 90, 10) < 764); // 10/log(10) * (90*log(9)+10*log(1/9))
STATIC_ASSERT_WITH_DEFAULT_MSG(abs(calc_binom_10log10_likeratio(0.1, 1, 99)) < 1e-4); // 10/log(10) * (90*log(9)+10*log(1/9))

template <class T=int, class T1>
T
collectget(const T1 & collection, size_t idx, T defaultval = 0) {
    return (idx < collection.size() ? collection[idx] : defaultval);
}

template <class T1, class T2>
void
clear_push(T1 & collection, T2 v, size_t idx = 0) {
    if (0 == idx) {
        collection.clear();
    }
    collection.push_back(v);
}

template <class T> 
void 
autoswap ( T& a, T& b ) {
    T c(a); a=b; b=c;
}

template <class T>
std::string
string_join(const T & container, std::string sep = std::string(",")) {
    std::string ret = "";
    for (auto e : container) {
        ret += e + sep;
    }
    ret.pop_back();
    return ret;
}

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

template <class T>
std::string 
int32t_join(const T & container, std::string sep = std::string(",")) {
    std::string ret = "";
    for (const auto & e : container) {
        if (e == INT32_MIN) {
            ret += ".,";
        } else {
            ret += std::to_string(e) + sep;
        }
    }
    if (ret.size() > 0) { ret.pop_back(); }
    return ret;
}

// variant-call data structures and functions

enum AlignmentSymbol {
    BASE_A,  //   = 0,
    BASE_C,  //   = 1,
    BASE_G,  //   = 2,
    BASE_T,  //   = 3,
    BASE_N,  //   = 4, // ambigous in the original sequencing data
    BASE_NN, //   = 5, // NA not available
    LINK_M,  //   = 6, // absence of any gap
    LINK_D3P,// = 7, // deletion of length 3 or plus
    LINK_D2, //  = 8,  // deletion of length 2
    LINK_D1, //  = 9,
    LINK_I3P,//  = 10, // insertion of length 1 // where the inserted sequence is not a repeat
    LINK_I2, //  = 11, 
    LINK_I1, // = 12, 
    LINK_NN, //  = 13, // padded InDel-noInDel symbol in deleted sequence
    END_ALIGNMENT_SYMBOLS,
    MGVCF_SYMBOL,
    ADDITIONAL_INDEL_CANDIDATE_SYMBOL,
};

const char* SYMBOL_TO_DESC_ARR[] = {
    [BASE_A] = "A", [BASE_C] = "C", [BASE_G] = "G", [BASE_T] = "T", [BASE_N] = "N",
    [BASE_NN] = "*",
    [LINK_M] = "<LR>",
    [LINK_D3P] = "<LD3P>", [LINK_D2] = "<LD2>", [LINK_D1] = "<LD1>",
    [LINK_I3P] = "<LI3P>", [LINK_I2] = "<LI2>", [LINK_I1] = "<LI1>",
    [LINK_NN] = "*",
    [END_ALIGNMENT_SYMBOLS] = "<NONE>",
    [MGVCF_SYMBOL] = "<NON_REF>",
    [ADDITIONAL_INDEL_CANDIDATE_SYMBOL] = "<ADDITIONAL_INDEL_CANDIDATE>",
};

#define NUM_ALIGNMENT_SYMBOLS 14
STATIC_ASSERT_WITH_DEFAULT_MSG(NUM_ALIGNMENT_SYMBOLS == END_ALIGNMENT_SYMBOLS);

// Please note that there are left-to-right clips and right-to-left clips, but I cannot know in advance if the direction matters in consensus. 
// My intuition tells me that it matters but only for some extremely rare situations, so I did not divide soft-clips according to their strands. 
// #define NUM_CLIP_SYMBOLS 2
// const std::array<ClipSymbol, NUM_CLIP_SYMBOLS> CLIP_SYMBOLS = {{CLIP_LEFT_TO_RIGHT, CLIP_RIGHT_TO_LEFT}};

#define NUM_INS_SYMBOLS 3
const std::array<AlignmentSymbol, NUM_INS_SYMBOLS> INS_SYMBOLS = {{LINK_I1, LINK_I2, LINK_I3P}};

#define NUM_DEL_SYMBOLS 3
const std::array<AlignmentSymbol, NUM_DEL_SYMBOLS> DEL_SYMBOLS = {{LINK_D1, LINK_D2, LINK_D3P}};

const std::array<AlignmentSymbol, (NUM_INS_SYMBOLS+NUM_DEL_SYMBOLS)> INDEL_SYMBOLS = {{LINK_I1, LINK_I2, LINK_I3P, LINK_D1, LINK_D2, LINK_D3P}};

constexpr bool
areSymbolsMutated(AlignmentSymbol ref, AlignmentSymbol alt) {
    if (alt <= BASE_NN) {
        return ref != alt && ref < BASE_N && alt < BASE_N;
    } else {
        return alt != LINK_M && alt != LINK_NN;
    }
};



enum SymbolType {
    BASE_SYMBOL,
    LINK_SYMBOL,
    NUM_SYMBOL_TYPES,
};

enum LinkType {
    MAT_LINK,
    INS_LINK,
    DEL_LINK,
    NUM_LINK_TYPES,
};

const AlignmentSymbol SYMBOL_TYPE_TO_INCLU_BEG[NUM_SYMBOL_TYPES] = {
    [BASE_SYMBOL] = BASE_A,
    [LINK_SYMBOL] = LINK_M,
};

const std::array<SymbolType, 2> SYMBOL_TYPE_ARR = {
    BASE_SYMBOL, LINK_SYMBOL
};

const std::array<std::vector<AlignmentSymbol>, NUM_SYMBOL_TYPES> SYMBOL_TYPE_TO_SYMBOLS = {{
    [BASE_SYMBOL] = std::vector<AlignmentSymbol>{{BASE_A,   BASE_C,   BASE_G,  BASE_T,   BASE_N,   BASE_NN}},
    [LINK_SYMBOL] = std::vector<AlignmentSymbol>{{LINK_M,   LINK_I1,  LINK_I2, LINK_I3P, LINK_D1,  LINK_D2,  LINK_D3P, LINK_NN}}
}};

const std::array<std::vector<AlignmentSymbol>, NUM_SYMBOL_TYPES> SYMBOL_TYPE_TO_NON_NN_SYMBOLS = {{
    [BASE_SYMBOL] = std::vector<AlignmentSymbol>{{ BASE_A, BASE_C, BASE_G,  BASE_T,   BASE_N }},
    [LINK_SYMBOL] = std::vector<AlignmentSymbol>{{ LINK_M, LINK_I1,LINK_I2, LINK_I3P, LINK_D1, LINK_D2, LINK_D3P }}
}};

const AlignmentSymbol SYMBOL_TYPE_TO_INCLU_END[NUM_SYMBOL_TYPES] = {
    [BASE_SYMBOL] = BASE_NN,
    [LINK_SYMBOL] = LINK_NN,
};

const AlignmentSymbol SYMBOL_TYPE_TO_AMBIG[NUM_SYMBOL_TYPES] = {
    [BASE_SYMBOL] = BASE_NN,
    [LINK_SYMBOL] = LINK_NN,
};

constexpr bool
isSymbolIns(const AlignmentSymbol symbol) {
    if (LINK_I3P == symbol || LINK_I2 == symbol || LINK_I1 == symbol) {
        return true;
    } else {
        return false;
    }
}

constexpr bool
isSymbolDel(const AlignmentSymbol symbol) {
     if (LINK_D3P == symbol || LINK_D2 == symbol || LINK_D1 == symbol) {
        return true;
    } else {
        return false;
    }
}

constexpr AlignmentSymbol
insLenToSymbol(uvc1_readpos_t len, const bam1_t *b) {
    assertUVC(len >= 0 || !fprintf(stderr, "Error: the bam record with qname %s at tid %d pos %ld has insertion of length %d !\n",
            bam_get_qname(b), b->core.tid, b->core.pos, len));
    return (1 == len ? LINK_I1 : ((2 == len) ? LINK_I2 : LINK_I3P));
}

constexpr AlignmentSymbol
delLenToSymbol(uvc1_readpos_t len, const bam1_t *b) {
    assertUVC(len >= 0 || !fprintf(stderr, "Error: the bam record with qname %s at tid %d pos %ld has deletion of length %d !\n",
            bam_get_qname(b), b->core.tid, b->core.pos, len));
    return (1 == len ? LINK_D1 : ((2 == len) ? LINK_D2 : LINK_D3P));
}



int
insSymbolToInsIdx(AlignmentSymbol s) {
    return (LINK_I1 == s ? 0 : ((LINK_I2 == s) ? 1: 2));
}

int
delSymbolToDelIdx(AlignmentSymbol s) {
    return (LINK_D1 == s ? 0 : ((LINK_D2 == s) ? 1: 2));
}

const std::array<SymbolType, 2> SYMBOL_TYPES_IN_VCF_ORDER = {{LINK_SYMBOL, BASE_SYMBOL}};

bool
isSymbolSubstitution(AlignmentSymbol symbol) {
    return (SYMBOL_TYPE_TO_INCLU_BEG[BASE_SYMBOL] <= symbol && symbol <= SYMBOL_TYPE_TO_INCLU_END[BASE_SYMBOL]);
}






struct _CharToSymbol {
    std::array<AlignmentSymbol, 128> data;
    _CharToSymbol() {
        for (size_t i = 0; i < 128; i++) {
            data[i] = BASE_N;
        }
        data['A'] = data['a'] = BASE_A;
        data['C'] = data['c'] = BASE_C;
        data['G'] = data['g'] = BASE_G;
        data['T'] = data['t'] = BASE_T;
        data['I'] = data['i'] = LINK_M;
        data['-'] = data['_'] = LINK_D1;
    }
};

const _CharToSymbol CHAR_TO_SYMBOL;

struct TumorKeyInfo {
    std::string ref_alt;
    int32_t VTI = -1;
    bool enable_tier2_consensus_format_tags = false;
    uvc1_refgpos_t pos = 0;
    
    uvc1_readnum_t BDP = 0;
    uvc1_readnum_t bDP = 0;
    
    uvc1_readnum100x_t CDP1x = 0;
    uvc1_readnum100x_t cDP1x = 0;
    uvc1_qual_t cVQ1 = 0;
    uvc1_qual_t cPCQ1 = 0;
    
    uvc1_readnum100x_t CDP2x = 0;
    uvc1_readnum100x_t cDP2x = 0;
    uvc1_qual_t cVQ2 = 0;
    uvc1_qual_t cPCQ2 = 0;
    
    uvc1_qual_t bNMQ = 0;
    uvc1_qual_t vHGQ = 0;
    bcf1_t *bcf1_record = NULL;
    
    uvc1_readnum_t tDP = 0;
    std::array<uvc1_readnum_t, 2> tADR = {{ 0 }};
    uvc1_readnum_t nDP = 0;
    std::array<uvc1_readnum_t, 2> nADR = {{ 0 }};
    uvc1_readnum_t tDPC = 0;
    std::array<uvc1_readnum_t, 2> tADCR = {{ 0 }};
    std::array<uvc1_readnum_t, 2> nADCR = {{ 0 }}; 
    
    /*
    ~TumorKeyInfo() {
        if (bcf1_record != NULL) {
            // this line must be elsewhere apparently due to the subtle differences between C and C++.
            // bcf_destroy(bcf1_record);
        }
    }
    */
};

std::basic_string<AlignmentSymbol>
string2symbolseq(const std::string & instring) {
    std::basic_string<AlignmentSymbol> ret;
    ret.reserve(instring.size());
    for (size_t i = 0; i < instring.size(); i++) {
        ret.push_back(CHAR_TO_SYMBOL.data[instring[i]]);
    }
    return ret;
};

struct SegFormatPrepSet {
    uvc1_readnum_t segprep_a_dp;
    uvc1_readnum_t segprep_a_near_ins_dp;
    uvc1_readnum_t segprep_a_near_del_dp;
    uvc1_readnum_t segprep_a_near_RTR_ins_dp;
    uvc1_readnum_t segprep_a_near_RTR_del_dp;
    
    // std::array<uvc1_readnum_t, 2> segprep_a_pcr_dps = {{ 0 }}; // depth of PCR amplicons on the R1R2 and R2R1 orientations
    uvc1_readnum_t segprep_a_pcr_dp;       // depth of PCR-amplicon sequencing segments
    uvc1_readnum_t segprep_a_umi_dp;       // depth of UMI-enabled sequencing segments
    uvc1_readnum_t segprep_a_snv_dp;
    uvc1_readnum_t segprep_a_dnv_dp;
    uvc1_readnum_t segprep_a_highBQ_dp;    // depth of high-BQ bases
    
    uvc1_readnum_t segprep_a_near_pcr_clip_dp;
    uvc1_readnum_t segprep_a_near_long_clip_dp;

    uvc1_readnum_t segprep_a_at_ins_dp;
    uvc1_readnum_t segprep_a_at_del_dp;

    // not really 100x but at the same order of magnitude
    uvc1_readnum100x_t segprep_a_XM1500;   // number of mismatches per 1500 bases
    uvc1_readnum100x_t segprep_a_GO1500;   // number of gap openings per 1500 bases
    uvc1_readnum100x_t segprep_a_GAPLEN;   
    uvc1_readnum100x_t segprep_a_qlen;
    
    uvc1_readpos_big_t segprep_a_near_ins_pow2len;
    uvc1_readpos_big_t segprep_a_near_del_pow2len;
    uvc1_readnum100x_t segprep_a_near_ins_inv100len;
    uvc1_readnum100x_t segprep_a_near_del_inv100len;

    uvc1_readpos_big_t segprep_a_near_ins_l_pow2len;    
    uvc1_readpos_big_t segprep_a_near_ins_r_pow2len;    
    uvc1_readpos_big_t segprep_a_near_del_l_pow2len;    
    uvc1_readpos_big_t segprep_a_near_del_r_pow2len;    
    
    uvc1_readpos_big_t segprep_a_LI;
    uvc1_readnum_t segprep_a_LIDP;
    uvc1_readpos_big_t segprep_a_RI;
    uvc1_readnum_t segprep_a_RIDP;
    
    uvc1_readpos_t segprep_a_l_dist_sum;
    uvc1_readpos_t segprep_a_r_dist_sum;
    uvc1_readpos_t segprep_a_inslen_sum;
    uvc1_readpos_t segprep_a_dellen_sum;
    
    uvc1_qual_big_t segprep_a_l_BAQ_sum;
    uvc1_qual_big_t segprep_a_r_BAQ_sum;
    uvc1_qual_big_t segprep_a_insBAQ_sum;
    uvc1_qual_big_t segprep_a_delBAQ_sum;
    
#if COMPILATION_TRY_HIGH_DEPTH_POS_BIAS
    // data-driven border for position bias
    uvc1_readnum_t segprep_aa_l_ins_dist_x_wei;
    uvc1_readnum_t segprep_aa_l_ins_weight;
    uvc1_readnum_t segprep_aa_r_ins_dist_x_wei;
    uvc1_readnum_t segprep_aa_r_ins_weight;
    
    uvc1_readnum_t segprep_aa_l_del_dist_x_wei;
    uvc1_readnum_t segprep_aa_l_del_weight;
    uvc1_readnum_t segprep_aa_r_del_dist_x_wei;
    uvc1_readnum_t segprep_aa_r_del_weight;
#endif

};
#define NUM_SEG_FORMAT_PREP_SETS ((size_t)SEG_FORMAT_PREP_SET_END)

template <class T1, class T2>
uvc1_readnum_big_t
calc_indel_weight(const T1 indelsize, const T2 borderlen) {
    return (1024L * 1024L) * mathcube(indelsize) / mathcube(MAX(borderlen, 8));
}

struct SegFormatThresSet {
    uvc1_readpos_t segthres_aLPxT;
    uvc1_readpos_t segthres_aRPxT;
    
#if COMPILATION_ENABLE_XMGOT
    uvc1_base1500x_t segthres_aXM1T; // mismatch, higher means more bias
    uvc1_base1500x_t segthres_aXM2T; 
    uvc1_base1500x_t segthres_aGO1T; // gap-open, higher means more bias
    uvc1_base1500x_t segthres_aGO2T; 
#endif
    
    uvc1_readpos_t segthres_aLI1T; // distance to left insert end, higher means more bias
    uvc1_readpos_t segthres_aLI2T; 
    uvc1_readpos_t segthres_aRI1T; // distance to right insert end
    uvc1_readpos_t segthres_aRI2T; 
    uvc1_readpos_t segthres_aLI1t; // distance to left insert end, lower means more bias
    uvc1_readpos_t segthres_aLI2t; 
    uvc1_readpos_t segthres_aRI1t; // distance to right insert end, lower means more bias
    uvc1_readpos_t segthres_aRI2t;  
    
    uvc1_readpos_t segthres_aLP1t;
    uvc1_readpos_t segthres_aLP2t;
    uvc1_readpos_t segthres_aRP1t;
    uvc1_readpos_t segthres_aRP2t;
    
    uvc1_qual_t segthres_aLB1t;
    uvc1_qual_t segthres_aLB2t;
    uvc1_qual_t segthres_aRB1t;
    uvc1_qual_t segthres_aRB2t;
};

struct SegFormatInfoSet {
    // allele-specific
    uvc1_readnum100x_t seginfo_a2XM2;
    uvc1_readnum100x_t seginfo_a2BM2;    
   
    uvc1_readnum100x_t seginfo_aPF1; // BQ without mismatches
    uvc1_readnum100x_t seginfo_aPF2; 
    
    uvc1_readnum_t seginfo_aBQ2; // number passing BQ filter
    
    uvc1_qual_t seginfo_aMQs;
    uvc1_readnum_t seginfo_aP1;
    uvc1_readnum_t seginfo_aP2;
    uvc1_readnum_t seginfo_aP3;
    uvc1_readnum_t seginfo_aNC;

    uvc1_readnum_t seginfo_aDPff;
    uvc1_readnum_t seginfo_aDPfr;
    uvc1_readnum_t seginfo_aDPrf;
    uvc1_readnum_t seginfo_aDPrr;
    
    uvc1_readnum_t seginfo_aLP1; // left seg pos
    uvc1_readnum_t seginfo_aLP2;
    uvc1_readpos_t seginfo_aLPL;
    uvc1_readnum_t seginfo_aRP1; // right seg pos
    uvc1_readnum_t seginfo_aRP2;
    uvc1_readpos_t seginfo_aRPL;
    
    uvc1_readnum_t seginfo_aLB1; // left seg pos
    uvc1_readnum_t seginfo_aLB2;
    uvc1_readpos_big_t seginfo_aLBL;
    uvc1_readnum_t seginfo_aRB1; // right seg pos
    uvc1_readnum_t seginfo_aRB2;
    uvc1_readpos_big_t seginfo_aRBL;
    
    uvc1_readnum_t seginfo_aLI1; // left insert
    uvc1_readnum_t seginfo_aLI2;
    
    uvc1_readnum_t seginfo_aRI1; // right insert
    uvc1_readnum_t seginfo_aRI2;
    
    uvc1_readnum_t seginfo_aRIf;
    uvc1_readnum_t seginfo_aLIr;

    uvc1_readpos_big_t seginfo_aLIT;
    uvc1_readpos_big_t seginfo_aRIT;
};

enum FragFormatDepthSet {
    FRAG_bDP, // raw
    FRAG_bTA, // FRAG_b10xSeqTlen,
    FRAG_bTB, // FRAG_b10xSeqTNevents,
    FRAG_FORMAT_DEPTH_SET_END
};
#define NUM_FRAG_FORMAT_DEPTH_SETS ((size_t)FRAG_FORMAT_DEPTH_SET_END)

struct FamFormatInfoSet {
    uvc1_readnum_t faminfo_c2LP1; // left tier-2-consensus family pos
    uvc1_readnum_t faminfo_c2LP2;
    uvc1_readpos_t faminfo_c2LPL;
    uvc1_readnum_t faminfo_c2RP1; // right tier-2-consensus family pos
    uvc1_readnum_t faminfo_c2RP2;
    uvc1_readpos_t faminfo_c2RPL;
    
    uvc1_readnum_t faminfo_c2LP0;
    uvc1_readnum_t faminfo_c2RP0;
    
    uvc1_readnum_t faminfo_c2LB1; // left tier-2-consensus family pos
    uvc1_readnum_t faminfo_c2LB2;
    uvc1_readpos_big_t faminfo_c2LBL;
    uvc1_readnum_t faminfo_c2RB1; // right tier-2-consensus family pos
    uvc1_readnum_t faminfo_c2RB2;
    uvc1_readpos_big_t faminfo_c2RBL;

    uvc1_readnum_t faminfo_c2BQ2;
};

enum FamFormatDepthSet {
    FAM_cDP1, // raw
    FAM_cDP12,// filtered
    FAM_cDP2, //  2, 0.8, family-consensus
    FAM_cDP3, // 10, 0.8, family-consensus
    FAM_cDPM, // duped match
    FAM_cDPm, // duped mismatch
    FAM_cDP21, // singleton
    
    FAM_FORMAT_DEPTH_SET_END
};
#define NUM_FAM_FORMAT_DEPTH_SETS ((size_t)FAM_FORMAT_DEPTH_SET_END)

enum DuplexFormatDepthSet {
    DUPLEX_dDP1,  // raw
    DUPLEX_dDP2, // double-strand consensus
    DUPLEX_FORMAT_TAG_SET_END
};
#define NUM_DUPLEX_FORMAT_DEPTH_SETS ((size_t)DUPLEX_FORMAT_TAG_SET_END)

enum VQFormatTagSet {
    
    VQ_a1BQf,
    VQ_a1BQr,
    VQ_a2BQf,
    VQ_a2BQr,
        
    VQ_bMQ,
    
    VQ_bIAQb, // prefinal
    VQ_bIADb, 
    VQ_bIDQb,
    
    VQ_cIAQf,
    VQ_cIADf,
    VQ_cIDQf,
    
    VQ_cIAQr,
    VQ_cIADr,
    VQ_cIDQr,
    
    // later computed
    VQ_aBQQ, // prefinal
    VQ_bIAQ, // prefinal
    VQ_cIAQ, // prefinal
    
    VQ_aPLQ, // preprefinal
    VQ_c1PLQ, // prefinal, deduped
    VQ_c2PLQ, // prefinal, consensus-applied
    VQ_dPLQ, // prefinal less priority
    
    VQ_C1DPv, // 100 times higher than expected
    VQ_c1DPv, // 100 times higher than expected
    VQ_c1VQ,  // final VarQual
    VQ_C2DPv, // 100 times higher than expected
    VQ_c2DPv, // ...
    VQ_c2VQ,  // ...
    
    VQ_FORMAT_TAG_SET_END
};
#define NUM_VQ_FORMAT_TAG_SETS ((size_t)VQ_FORMAT_TAG_SET_END)

template <class T>
uvc1_readnum_t
seg_format_get_ad(const T & s) {
    return s.seginfo_aDPff + s.seginfo_aDPfr + s.seginfo_aDPrf + s.seginfo_aDPrr;
};

template <class T1, class T2, class T3>
uvc1_qual_t
get_avgBQ(const T1 & bg_seg_bqsum_conslogo, const T2 & symbol_to_seg_format_depth_sets, const uvc1_refgpos_t epos, const T3 s) {
    const auto denom = seg_format_get_ad(symbol_to_seg_format_depth_sets.getByPos(epos)[s]);
    return bg_seg_bqsum_conslogo.getByPos(epos).getSymbolCount(s) / MAX(1, denom);
}

template
<bool TBidirectional = true, bool TIsOverseqFracDisabled = false>
std::array<double, 2> 
dp4_to_pcFA(double overseq_frac, double aADpass, double aADfail, double aDPpass, double aDPfail, 
        double pl_exponent = 3.0, double n_nats = log(500+1),
        double aADavgKeyVal = -1, double aDPavgKeyVal = -1, double priorAD = 0.5, double priorDP = 1.0) {
    assertUVC(aADpass >= 0.0 || !fprintf(stderr, "%f >= %f failed for pass!\n", aADpass, 0.0));
    assertUVC(aADfail >= 0.0 || !fprintf(stderr, "%f >= %f failed for fail!\n", aADfail, 0.0));
    assertUVC(aADpass <= aDPpass || !fprintf(stderr, "%f <= %f failed for pass!\n", aADpass, aDPpass));
    assertUVC(aADfail <= aDPfail || !fprintf(stderr, "%f <= %f failed for fail!\n", aADfail, aDPfail));
    if (!TIsOverseqFracDisabled) {
        aDPfail *= overseq_frac;
        aDPpass *= overseq_frac;
        aADfail *= overseq_frac;
        aADpass *= overseq_frac;
    }
    aDPfail += priorDP;
    aDPpass += priorDP;
    aADfail += priorAD;
    aADpass += priorAD;
    const double nobiasFA = (aADfail + aADpass) / (aDPfail + aDPpass); 
    if ((aADpass / aDPpass) >= (aADfail / aDPfail)) {
        if (TBidirectional) {
            autoswap(aDPfail, aDPpass);
            autoswap(aADfail, aADpass);
        } else {
            return std::array<double, 2> {{ (aADpass / aDPpass), nobiasFA }};
        }
    }
    auto aBDfail = aDPfail * 2 - aADfail * 1;
    auto aBDpass = aDPpass * 2 - aADpass * 1;
    assertUVC (aBDfail > 0);
    assertUVC (aBDpass > 0);
    double aADpassfrac = aADpass / (aADpass + aADfail);
    double aBDpassfrac = aBDpass / (aBDpass + aBDfail);
    if ((!TBidirectional) && (aADavgKeyVal >= 0) && (aDPavgKeyVal >= 0)) {
        aADpassfrac = aADavgKeyVal / (aADavgKeyVal + aDPavgKeyVal * 0.9); // interpolate
        aBDpassfrac = 1.0 - aADpassfrac;
    }
    double infogain = aADfail * log((1.0 - aADpassfrac) / (1.0 - aBDpassfrac));
    if (TBidirectional) { 
        infogain += aADpass * log(aADpassfrac / aBDpassfrac); 
    }
#ifdef TEST_dp4_to_pcFA
    printf("infogain = %f\n", infogain);
#endif
    if (infogain <= n_nats) {
        return std::array<double, 2> {{ aADfail / aDPfail, nobiasFA }};
    } else {
        return std::array<double, 2> {{ MAX(aADpass / aDPpass, (aADfail / aDPfail) * exp((n_nats - infogain) / pl_exponent)), nobiasFA }};
    }
}

#ifdef TEST_dp4_to_pcFA
int 
main(int argc, char **argv) {
    double adpass = atof(argv[1]);
    double adfail = atof(argv[2]);
    double dppass = atof(argv[3]);
    double dpfail = atof(argv[4]);
    double entropy = atof(argv[5]);
    double entrmax = atof(argv[6]);
    double ldist = atof(argv[7]);
    double rdist = atof(argv[8]);
    double pca = atof(argv[9]);
    double pcb = atof(argv[10]);
    
    const auto ret1 = dp4_to_pcFA<false>(-1, adpass, adfail, dppass, dpfail, entropy, entrmax, ldist, rdist, pca, pcb);
    const auto ret2 = dp4_to_pcFA<true>(-1, adpass, adfail, dppass, dpfail, entropy, entrmax, ldist, rdist, pca, pcb);
    printf("dp4_to_pcFA<(false AND true)>(%f, %f, %f, %f, %f, %f, %f, %f, %f, %f) = ({%f, %f} AND {%f, %f})\n", 
            adpass, adfail, dppass, dpfail, entropy, entrmax, ldist, rdist, pca, pcb, ret1[0], ret1[1], ret2[0], ret2[1]);
}

#endif

// conversion between different defintions in bioinformatics

uvc1_qual_t
char2phred(const char charvalue) {
    return charvalue - 33;
}

char 
phred2char(const uvc1_qual_t phredvalue) {
    return phredvalue + 33;
}

double 
phred2prob(const uvc1_qual_t phredvalue) {
    return pow(10, -((float)phredvalue) / 10);
}

uvc1_qual_t 
prob2phred(const double probvalue) {
    return floor(-10 * log(probvalue) / log(10));
}

double 
prob2realphred(const double probvalue) {
    return -10 * log(probvalue) / log(10);
}

template <class T1, class T2, class T3, class T4>
void
process_cigar(T1 & qpos, T2 & rpos, T3 cigar_op, T4 cigar_oplen) {
    if (cigar_op == BAM_CREF_SKIP) {
        rpos += cigar_oplen;
    } else if (cigar_op == BAM_CSOFT_CLIP) {
        qpos += cigar_oplen;
    } else if (cigar_op == BAM_CHARD_CLIP) {
        // pass
    } else if (cigar_op == BAM_CPAD) {
        // pass
    } else if (cigar_op == BAM_CBACK) {
        throw -1;
    } else {
        throw -2;
    }
}

// variant-call math functions

#define NUM_BUCKETS 16

template <class T>
uvc1_qual_t 
proton_cigarlen2phred(T cigarlen) {
    uvc1_qual_t oplen2phred[12+1] = {
        0, // generated by the python code: for i in range(1,12+1):  print('(uvc_qual_t)round({}), //{}'.format(10/log(10)*log(i**3), i))
        (uvc1_qual_t)round(0.0), //1
        (uvc1_qual_t)round(9.030899869919434), //2
        (uvc1_qual_t)round(14.313637641589871), //3
        (uvc1_qual_t)round(18.061799739838868), //4
        (uvc1_qual_t)round(20.969100130080562), //5
        (uvc1_qual_t)round(23.344537511509305), //6
        (uvc1_qual_t)round(25.352941200427697), //7
        (uvc1_qual_t)round(27.092699609758302), //8
        (uvc1_qual_t)round(28.627275283179742), //9
        (uvc1_qual_t)round(29.999999999999993), //10
        (uvc1_qual_t)round(31.241780554746747), //11
        (uvc1_qual_t)round(32.37543738142874), //12
    };
    return oplen2phred[MIN(cigarlen, 12)]; 
}

template <class T>
int 
infer_max_qual_assuming_independence(
        uvc1_qual_t & maxvqual,
        uvc1_readnum_t & argmaxAD,
        uvc1_qual_t & argmaxBQ,
        const uvc1_qual_t max_qual, 
        const uvc1_qual_t dec_qual,
        const std::array<T, NUM_BUCKETS> & qual_distr, 
        const uvc1_readnum_t totDP,
        const uvc1_hash_t specialflag IGNORE_UNUSED_PARAM) {
    
    uvc1_qual_t currvqual = 0;
    uvc1_readnum_t currAD = 0;
    maxvqual = 0; 
    argmaxAD = 0;
    argmaxBQ = 0;
    for (uvc1_qual_t idx = 0; idx < MIN(NUM_BUCKETS, max_qual / dec_qual); idx++) {
        const auto currQD = qual_distr[idx];
        if (0 == currQD) { continue; }
        currAD += currQD;
        auto currBQ = max_qual - (dec_qual * idx);
        double expBQ = 10.0 / log(10.0) * log(((double)totDP / (double)currAD) + DBL_EPSILON);
        currvqual = (uvc1_qual_t)(currAD * (currBQ - expBQ));
        if (currvqual > maxvqual) {
            argmaxAD = currAD;
            argmaxBQ = currBQ;
            maxvqual = currvqual;
        }
    }
    return 0;
}

#endif

