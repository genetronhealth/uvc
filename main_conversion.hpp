#ifndef conversion_hpp_INCLUDED
#define conversion_hpp_INCLUDED

#include "common.hpp"

#include "htslib/sam.h"
#include "htslib/vcf.h"

#include <algorithm>
#include <iostream>

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define SQR_QUAL_DIV 32

// math functions

auto MEDIAN(const auto & v) {
    assert(v.size() > 0);
    return (v[(v.size() - 1) / 2] + v[(v.size()) / 2]) / 2;
}

auto FIRST(const auto & v) {
    assert(v.size() > 0);
    return v[0];
}

const auto LAST(const auto & v) {
    assert(v.size() > 0);
    return v[(v.size()-1)];
}

int64_t
int64mul(const auto a, const auto b) {
    return ((int64_t)a) * ((int64_t)b);
}

uint64_t
uint64mul(const auto a, const auto b) {
    return ((uint64_t)a) * ((uint64_t)b);
}

auto 
non_neg_minus(auto a, auto b) {
    return (a > b ? (a - b) : 0);
}

auto
CENTER(auto a, auto b, int center = 0) {
    return ((abs(a - center) < abs(b - center)) ? a : b);
}

/*
auto 
MIN(auto a, auto b) {
    return (a < b ? a : b);
}
*/

auto
MIN3(auto a, auto b, auto c) {
    return MIN(a, MIN(b, c));
}

auto
MIN4(auto a, auto b, auto c, auto d) {
    return MIN(a, MIN3(b, c, d));
}

auto
MIN5(auto a, auto b, auto c, auto d, auto e) {
    return MIN(a, MIN4(b, c, d, e));
}

auto
MINVEC(auto v) {
    assert(v.size() > 0);
    auto ret = v[0];
    for (auto e : v) {
        ret = MIN(ret, e);
    }
    return ret;
}

/*
auto 
MAX(auto a, auto b) {
    return (a > b ? a :b);
}
*/

auto 
MAX3(auto a, auto b, auto c) {
    return MAX(a, MAX(b, c));
}

auto 
MAX4(auto a, auto b, auto c, auto d) {
    return MAX(a, MAX3(b, c, d));
}

auto 
MAX5(auto a, auto b, auto c, auto d, auto e) {
    return MAX(a, MAX4(b, c, d, e));
}

auto 
MAX6(auto a, auto b, auto c, auto d, auto e, auto f) {
    return MAX(a, MAX5(b, c, d, e, f));
}

auto
MAXVEC(auto v) {
    assert(v.size() > 0);
    auto ret = v[0];
    for (auto e : v) {
        ret = MAX(ret, e);
    }
    return ret;
}

auto
BETWEEN(auto v, auto a, auto b) {
    return MIN(MAX(a, v), b);
}

void 
UPDATE_MIN(auto & a, const auto & b) {
    a = MIN(a, b);
}

void 
UPDATE_MAX(auto & a, const auto & b) {
    a = MAX(a, b);
}

auto
SUMVEC(const auto & vec) {
    auto r = 0;
    for (unsigned int i = 0; i < vec.size(); i++) {
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
calc_score_with_dimret(const T v, T penal_mult, T thres) {
    if (v > thres) {
        return (v - thres) * penal_mult + thres;
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

auto 
safediv0(auto a, auto b) {
    return (b != 0 ? a / b : 0);
}

auto 
mathsquare(auto x) { 
    return x * x; 
}

double
prob2odds(double p) {
    assert((0.0 < p && p < 1.0) || !fprintf(stderr, "%f is not between 0 and 1!", p));
    return p /  (1.0 - p);
}

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
    assert((prob > 0 && prob < 1) || !fprintf(stderr, "The assertion 0 < %f < 1 failed!\n", prob));
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

static_assert(abs(calc_binom_10log10_likeratio(0.1, 10, 90)) < 1e-4);
static_assert(calc_binom_10log10_likeratio(0.1, 90, 10) > 763); // 10/log(10) * (90*log(9)+10*log(1/9))
static_assert(calc_binom_10log10_likeratio(0.1, 90, 10) < 764); // 10/log(10) * (90*log(9)+10*log(1/9))
static_assert(abs(calc_binom_10log10_likeratio(0.1, 1, 99)) < 1e-4); // 10/log(10) * (90*log(9)+10*log(1/9))

//template <class T, class V>
template <class T=int>
T
collectget(const auto & collection, size_t idx, T defaultval = 0) {
    return (idx < collection.size() ? collection[idx] : defaultval);
}

void
clear_push(auto & collection, auto v, size_t idx = 0) {
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

std::string
string_join(const auto & container, std::string sep = std::string(",")) {
    std::string ret = "";
    for (auto e : container) {
        ret += e + sep;
    }
    ret.pop_back();
    return ret;
}

std::string 
other_join(const auto & container, std::string sep = std::string(",")) {
    std::string ret = "";
    for (const auto & e : container) {
        ret += std::to_string(e) + sep;
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
    BASE_NN, //   = 5, // padded match-mismatch symbol in deleted sequence
    LINK_M,  //   = 6, // absence of any gap
    LINK_D3P,// = 7, // deletion of length 3 or plus
    LINK_D2, //  = 8,  // deletion of length 2
    LINK_D1, //  = 9,
    LINK_I3P,//  = 10, // insertion of length 1 // where the inserted sequence is not a repeat
    LINK_I2, //  = 11, 
    LINK_I1, // = 12, 
    LINK_NN, //  = 13, // padded InDel-noInDel symbol in deleted sequence
    END_ALIGNMENT_SYMBOLS,
    GVCF_SYMBOL,
};

#define NUM_ALIGNMENT_SYMBOLS 14
struct _CharToSymbol {
    std::array<AlignmentSymbol, 128> data;
    _CharToSymbol() {
        for (int i = 0; i < 128; i++) {
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
    // std::string FTS;
    uvc1_refgpos_t pos = 0;
    
    uvc1_readnum_t BDP = 0; // 
    uvc1_readnum_t bDP = 0; // 
    
    uvc1_readnum_t CDP12 = 0;
    uvc1_readnum_t cDP12 = 0;
    
    uvc1_readnum100x_t CDP1x = 0;
    uvc1_readnum100x_t cDP1x = 0;
    uvc1_qual_t cVQ1 = 0;
    uvc1_qual_t cPCQ1 = 0;
    
    uvc1_readnum100x_t CDP2x = 0;
    uvc1_readnum100x_t cDP2x = 0;
    uvc1_qual_t cVQ2 = 0;
    uvc1_qual_t cPCQ2 = 0;
    
    bcf1_t *bcf1_record = NULL;
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
    uvc1_readnum_t segprep_a_pcr_dp; // SEG_a_PCR_DP,    // depth of PCR amplicons
    uvc1_readnum_t segprep_a_highBQ_dp;      // SEG_a_highBQ_DP, // depth of high-BQ bases
    uvc1_readnum100x_t segprep_a_XM1500;         // SEG_a_XM,        // number of mismatches per 1500 bases
    uvc1_readnum100x_t segprep_a_GO1500;         // SEG_a_GO,        // number of gap openings per 1500 bases
    uvc1_readnum100x_t segprep_a_XM100inv;       // SEG_a_XM100INV,  // number of inverse 
    uvc1_readpos_t segprep_a_LI;    // SEG_a_LI,
    uvc1_readnum_t segprep_a_LIDP;  // SEG_a_LIDP,
    uvc1_readpos_t segprep_a_RI;    // SEG_a_RI,
    uvc1_readnum_t segprep_a_RIDP;  // SEG_a_RIDP,
    
    uvc1_readnum_t segprep_a_snv_dp;
    uvc1_readnum_t segprep_a_dnv_dp;

    uvc1_readpos_t segprep_a_l_dist_sum; // SEG_a_L_DIST_SUM,
    uvc1_readpos_t segprep_a_r_dist_sum; // SEG_a_R_DIST_SUM,
    uvc1_readpos_t segprep_a_inslen_sum; // SEG_a_INSLEN_SUM,
    uvc1_readpos_t segprep_a_dellen_sum; // SEG_a_DELLEN_SUM,
    
    uvc1_qual_t segprep_a_l_BAQ_sum;  //  SEG_a_L_BAQ_SUM,
    uvc1_qual_t segprep_a_r_BAQ_sum;  //  SEG_a_R_BAQ_SUM,
    uvc1_qual_t segprep_a_insBAQ_sum; //  SEG_a_INSBAQ_SUM,
    uvc1_qual_t segprep_a_delBAQ_sum; //  SEG_a_DELBAQ_SUM,
    
    uvc1_readpos_t segprep_a_near_ins_pow2len; //  SEG_a_NEAR_INS_LEN,
    uvc1_readpos_t segprep_a_near_del_pow2len; //  SEG_a_NEAR_DEL_LEN,
    uvc1_readnum100x_t segprep_a_near_ins_inv100len; //  SEG_a_NEAR_INS_INV100LEN,
    uvc1_readnum100x_t segprep_a_near_del_inv100len; //  SEG_a_NEAR_DEL_INV100LEN,

    uvc1_readnum_t segprep_a_at_ins_dp; //  SEG_a_AT_INS_DP,
    uvc1_readnum_t segprep_a_at_del_dp; //  SEG_a_AT_DEL_DP,
    uvc1_readnum_t segprep_a_near_ins_dp; // SEG_a_NEAR_INS_DP,
    uvc1_readnum_t segprep_a_near_del_dp; // SEG_a_NEAR_DEL_DP,
    uvc1_readnum_t segprep_a_near_RTR_ins_dp; // SEG_a_NEAR_RTR_INS_DP,
    uvc1_readnum_t segprep_a_near_RTR_del_dp; // SEG_a_NEAR_RTR_DEL_DP,
    
    // SEG_FORMAT_PREP_SET_END,
};
#define NUM_SEG_FORMAT_PREP_SETS ((size_t)SEG_FORMAT_PREP_SET_END)

/*
enum SegFormatThresSet {
    // SEG_aEP1t, // edge position, closer means more bias
    // SEG_aEP2t, 
    SEG_aLPxT,
    SEG_aRPxT,
    
    SEG_aXM1T, // mismatch, higher means more bias
    SEG_aXM2T, 
    SEG_aGO1T, // gap-open, higher means more bias
    SEG_aGO2T, 
    SEG_aLI1T, // distance to left insert end, higher means more bias
    SEG_aLI2T, 
    SEG_aRI1T, // distance to right insert end
    SEG_aRI2T, 
    SEG_aLI1t, // distance to left insert end, lower means more bias
    SEG_aLI2t, 
    SEG_aRI1t, // distance to right insert end, lower means more bias
    SEG_aRI2t,  
    
    SEG_aLP1t,
    SEG_aLP2t,
    SEG_aRP1t,
    SEG_aRP2t,
    
    SEG_aLB1t,
    SEG_aLB2t,
    SEG_aRB1t,
    SEG_aRB2t,
    
    SEG_FORMAT_THRES_SET_END
};
#define NUM_SEG_FORMAT_THRES_SETS ((size_t)SEG_FORMAT_THRES_SET_END)
*/
struct SegFormatThresSet {
    // SEG_aEP1t, // edge position, closer means more bias
    // SEG_aEP2t, 
    uvc1_readpos_t segthres_aLPxT;
    uvc1_readpos_t segthres_aRPxT;
    
    uvc1_base1500x_t segthres_aXM1T; // mismatch, higher means more bias
    uvc1_base1500x_t segthres_aXM2T; 
    uvc1_base1500x_t segthres_aGO1T; // gap-open, higher means more bias
    uvc1_base1500x_t segthres_aGO2T; 
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
// #define NUM_SEG_FORMAT_THRES_SETS ((size_t)SEG_FORMAT_THRES_SET_END)

/*
enum SegFormatDepthSet {
    SEG_aMQs,
    SEG_aXMp1,
    SEG_aP1,

    SEG_aDPff,
    SEG_aDPfr,
    SEG_aDPrf,
    SEG_aDPrr,
    
    SEG_aXM1,
    SEG_aXM2,
    SEG_aBM2,
    
    SEG_aBQ1, // base-quality bias
    SEG_aBQ2,
    
    SEG_aPF1, // mismatch
    SEG_aPF2, 
    
    SEG_aLP1, // left seg pos
    SEG_aLP2,
    SEG_aLPL,
    SEG_aRP1, // right seg pos
    SEG_aRP2,
    SEG_aRPL,
    
    SEG_aLB1, // left seg pos
    SEG_aLB2,
    SEG_aLBL,
    SEG_aRB1, // right seg pos
    SEG_aRB2,
    SEG_aRBL,
    
    SEG_aLI1, // left insert
    SEG_aLI2,
    SEG_aLILf,
    SEG_aLILr,
    
    SEG_aRI1, // right insert
    SEG_aRI2,
    SEG_aRILf,
    SEG_aRILr,
    
    SEG_FORMAT_DEPTH_SET_END
};
#define NUM_SEG_FORMAT_DEPTH_SETS ((size_t)SEG_FORMAT_DEPTH_SET_END)
*/

struct SegFormatInfoSet {
    uvc1_qual_t seginfo_aMQs;
    uvc1_readnum_t seginfo_aXMp1;
    uvc1_readnum_t seginfo_aP1;
    
    uvc1_readnum_t seginfo_aDPff;
    uvc1_readnum_t seginfo_aDPfr;
    uvc1_readnum_t seginfo_aDPrf;
    uvc1_readnum_t seginfo_aDPrr;
    
    uvc1_readnum100x_t seginfo_aXM1;
    uvc1_readnum100x_t seginfo_aXM2;
    uvc1_readnum100x_t seginfo_aBM2;
    
    uvc1_qual_t seginfo_aBQ1; // base-quality bias
    uvc1_qual_t seginfo_aBQ2;
    
    uvc1_readnum100x_t seginfo_aPF1; // mismatch
    uvc1_readnum100x_t seginfo_aPF2; 
    
    uvc1_readnum_t seginfo_aLP1; // left seg pos
    uvc1_readnum_t seginfo_aLP2;
    uvc1_readpos_t seginfo_aLPL;
    uvc1_readnum_t seginfo_aRP1; // right seg pos
    uvc1_readnum_t seginfo_aRP2;
    uvc1_readpos_t seginfo_aRPL;
    
    uvc1_readnum_t seginfo_aLB1; // left seg pos
    uvc1_readnum_t seginfo_aLB2;
    uvc1_readpos_t seginfo_aLBL;
    uvc1_readnum_t seginfo_aRB1; // right seg pos
    uvc1_readnum_t seginfo_aRB2;
    uvc1_readpos_t seginfo_aRBL;
    
    uvc1_readnum_t seginfo_aLI1; // left insert
    uvc1_readnum_t seginfo_aLI2;
    uvc1_readpos_t seginfo_aLILf;
    uvc1_readpos_t seginfo_aLILr;
    
    uvc1_readnum_t seginfo_aRI1; // right insert
    uvc1_readnum_t seginfo_aRI2;
    uvc1_readpos_t seginfo_aRILf;
    uvc1_readpos_t seginfo_aRILr;
};

enum FragFormatDepthSet {
    FRAG_bDP, // raw
    FRAG_FORMAT_DEPTH_SET_END
};
#define NUM_FRAG_FORMAT_DEPTH_SETS ((size_t)FRAG_FORMAT_DEPTH_SET_END)

enum FamFormatDepthSet {
    FAM_cDP1, // raw
    FAM_cDP12,// filtered
    FAM_cDP2, //  2, 0.8, family-consensus
    FAM_cDP3, // 10, 0.8, family-consensus
    FAM_cDPM, // duped match
    FAM_cDPm, // duped mismatch
    FAM_c1DP, //  singleton
    
    FAM_FORMAT_DEPTH_SET_END
};
#define NUM_FAM_FORMAT_DEPTH_SETS ((size_t)FAM_FORMAT_DEPTH_SET_END)

enum DuplexFormatDepthSet {
    DUPLEX_dDP1,  // raw
    DUPLEX_dDP2, // double-strand consensus
    DUPLEX_FORMAT_TAG_SET_END
};
#define NUM_DUPLEX_FORMAT_DEPTH_SETS ((size_t)DUPLEX_FORMAT_TAG_SET_END)

// MIN(A-BQ-syserr-QUAL-bidir (onlyIllumina),
//     B-MQ-QUAL,
//     // precursor A-FA-QUAL,
//     MAX(MIN(B-FA-QUAL (base0 ), B-DPQ-iiderr-QUAL), 
//         MIN(C-FA-QUAL (base30), C-DPQ-iiderr-QUAL-bidir))
//         MIN(D-FA-QUAL (base60), D-DPQ-iiderr-QUAL))
enum VQFormatTagSet {
    
    VQ_a1XM,
    
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

uvc1_readnum_t
seg_format_get_ad(const auto & s) {
    return s.seginfo_aDPff + s.seginfo_aDPfr + s.seginfo_aDPrf + s.seginfo_aDPrr;
};

uvc1_qual_t
get_avgBQ(const auto & bg_seg_bqsum_conslogo, const auto & symbol_to_seg_format_depth_sets, const uvc1_refgpos_t epos, const auto s) {
    const auto denom = seg_format_get_ad(symbol_to_seg_format_depth_sets.getByPos(epos)[s]);
    return bg_seg_bqsum_conslogo.getByPos(epos).getSymbolCount(s) / MAX(1, denom);
}

template
<bool TBidirectional = true>
std::array<double, 2> 
dp4_to_pcFA(double aADpass, double aADfail, double aDPpass, double aDPfail, 
        double pl_exponent = 3.0, double n_nats = log(500+1),
        double aADavgKeyVal = -1, double aDPavgKeyVal = -1, double priorAD = 0.5, double priorDP = 1.0) {
    assert(aADpass <= aDPpass || !fprintf(stderr, "%f <= %f failed for pass!\n", aADpass, aDPpass));
    assert(aADfail <= aDPfail || !fprintf(stderr, "%f <= %f failed for fail!\n", aADfail, aDPfail));
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
    assert (aBDfail > 0);
    assert (aBDpass > 0);
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
    
    const auto ret1 = dp4_to_pcFA<false>(adpass, adfail, dppass, dpfail, entropy, entrmax, ldist, rdist, pca, pcb);
    const auto ret2 = dp4_to_pcFA<true>(adpass, adfail, dppass, dpfail, entropy, entrmax, ldist, rdist, pca, pcb);
    printf("dp4_to_pcFA<false, true>(%f, %f, %f, %f, %f, %f, %f, %f, %f, %f) = <{%f, %f}, {%f, %f}>\n", adpass, adfail, dppass, dpfail, entropy, entrmax, ldist, rdist, pca, pcb, ret1[0], ret1[1], ret2[0], ret2[1]);
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

void
process_cigar(auto & qpos, auto & rpos, auto cigar_op, auto cigar_oplen) {
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

uvc1_qual_t 
proton_cigarlen2phred(auto cigarlen) {
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
    for (int idx = 0; idx < MIN(NUM_BUCKETS, max_qual / dec_qual); idx++) {
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

