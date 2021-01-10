#ifndef main_hpp_INCLUDED
#define main_hpp_INCLUDED

#include "bcf_formats.step1.hpp" // auto-generated

#include "CmdLineArgs.hpp"
#include "common.hpp"
#include "logging.hpp"
#include "main_conversion.hpp"

#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/synced_bcf_reader.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// T=uvc1_readnum_t for deletion and T=string for insertion
template <class T>
uvc1_readnum_t
posToIndelToData_get(const std::map<uvc1_readnum_t, std::map<T, uvc1_readnum_t>> & pos2indel2data, const uvc1_readnum_t refpos, const T & indel, uvc1_readnum_t default1 = 0, uvc1_readnum_t default2 = 0) {
    if (pos2indel2data.find(refpos) == pos2indel2data.end()) { return default1; }
    if (pos2indel2data.at(refpos).find(indel) == pos2indel2data.at(refpos).end()) { return default2; }
    return pos2indel2data.at(refpos).at(indel);
}

template <class T>
void
posToIndelToCount_inc(std::map<uvc1_readnum_t, std::map<T, uvc1_readnum_t>> & pos2indel2count, uvc1_readpos_t pos, const T indel, uvc1_readnum_t incvalue = 1) {
    assert(incvalue > 0);
    auto pos2indel2count4it = pos2indel2count.insert(std::make_pair(pos, std::map<T, uvc1_readnum_t>()));
    auto indel2count4it = pos2indel2count4it.first->second.insert(std::make_pair(indel, 0));
    indel2count4it.first->second += incvalue;
    assert (posToIndelToData_get(pos2indel2count, pos, indel) > 0);
}

template <class T>
T
posToIndelToCount_updateByConsensus(std::map<uvc1_readnum_t, std::map<T, uvc1_readnum_t>> & dst, const std::map<uvc1_readpos_t, std::map<T, uvc1_readnum_t>> & src, uvc1_readnum_t epos, uvc1_readnum_t incvalue = 1) {
    const auto pos2indel2count4it = src.find(epos);
    assert(pos2indel2count4it != src.end());
    auto indel2count = pos2indel2count4it->second;
    const T & src_indel = (indel2count.size() > 1 ? T() : indel2count.begin()->first);
    assert (indel2count.begin()->second > 0);
    posToIndelToCount_inc<T>(dst, epos, src_indel, incvalue);
    return src_indel;
}

template <bool TIsIncVariable, class T>
uvc1_readnum_t
posToIndelToCount_updateByRepresentative(std::map<uvc1_readnum_t, std::map<T, uvc1_readnum_t>> & dst, const std::map<uvc1_readnum_t, std::map<T, uvc1_readnum_t>> & src, uvc1_readnum_t epos, uvc1_readnum_t incvalue = 1) {
    const auto pos2indel2count4it = src.find(epos);
    assert(pos2indel2count4it != src.end());

    T max_indel;
    uvc1_readnum_t max_count = 0;
    uvc1_readnum_t sum_count = 0;
    for (auto indel2count : pos2indel2count4it->second) {
        auto indel = indel2count.first;
        auto count = indel2count.second;
        if (count > max_count) {
            max_indel = indel;
            max_count = count;
        }
        if (TIsIncVariable) {sum_count += count; }
    }
    assert(0 < max_count);
    if (TIsIncVariable) {
        posToIndelToCount_inc<T>(dst, epos, max_indel, sum_count);
    } else {
        posToIndelToCount_inc<T>(dst, epos, max_indel, incvalue);
    }
    return max_count;
}

template <class T>
void
posToIndelToCount_updateBySummation(std::map<uvc1_readnum_t, std::map<T, uvc1_readnum_t>> & dst, const std::map<uvc1_readnum_t, std::map<T, uvc1_readnum_t>> & src) {
    for (auto src_pos2indel2count4it : src) {
        auto src_pos = src_pos2indel2count4it.first;
        auto src_indel2count = src_pos2indel2count4it.second;
        for (auto src_indel2count4it : src_indel2count) {
            assert(src_indel2count4it.second > 0 || !(
                std::cerr << src_indel2count4it.second << " > 0 failed for the key " << src_indel2count4it.first << " at position " << src_pos << std::endl
            ));
            posToIndelToCount_inc<T>(dst, src_pos, src_indel2count4it.first, src_indel2count4it.second);
        }
    }
}

void
posToIndelToCount_DlenToDseq(std::map<uvc1_readnum_t, std::map<std::string, uvc1_readnum_t>> & dst, const std::map<uvc1_readnum_t, std::map<uvc1_readnum_t, uvc1_readnum_t>> & src,
        const std::string & refchars, uvc1_readnum_t incluBegPos) {
    for (auto src_pos2indel2count4it : src) {
        uvc1_readnum_t src_pos = src_pos2indel2count4it.first;
        auto src_del2count = src_pos2indel2count4it.second;
        for (auto src_del2count4it : src_del2count) {
            std::string dseq = refchars.substr(src_pos - incluBegPos, src_del2count4it.first);
            posToIndelToCount_inc(dst, src_pos, dseq, src_del2count4it.second);
        }
    }
}

#define INS_N_ANCHOR_BASES 1 //  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/
#define TVN_MICRO_VQ_DELTA 3
#define TIN_CONTAM_MICRO_VQ_DELTA 0
#define BAM_PHREDI(b, i) (bam_get_qual((b))[(i)])

const bcfrec::BcfFormat FORMAT_UNCOV = bcfrec::BcfFormat();

enum ValueType {
    SYMBOL_COUNT_SUM,
    BASE_QUALITY_MAX, // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6302405/
    VALUE_TYPE_END,
};

#define NUM_ALIGNMENT_SYMBOLS 14
static_assert(NUM_ALIGNMENT_SYMBOLS == END_ALIGNMENT_SYMBOLS);

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
    assert(len >= 0 || !fprintf(stderr, "Error: the bam record with qname %s at tid %d pos %d has insertion of length %d !\n", 
            bam_get_qname(b), b->core.tid, b->core.pos, len));
    return (1 == len ? LINK_I1 : ((2 == len) ? LINK_I2 : LINK_I3P));
}

constexpr AlignmentSymbol 
delLenToSymbol(uvc1_readpos_t len, const bam1_t *b) {
    assert(len >= 0 || !fprintf(stderr, "Error: the bam record with qname %s at tid %d pos %d has deletion of length %d !\n", 
            bam_get_qname(b), b->core.tid, b->core.pos, len));
    return (1 == len ? LINK_D1 : ((2 == len) ? LINK_D2 : LINK_D3P));
}

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

const std::array<std::vector<AlignmentSymbol>, NUM_SYMBOL_TYPES> SYMBOL_TYPE_TO_SYMBOLS= {{
    [BASE_SYMBOL] = std::vector<AlignmentSymbol>{{BASE_A,   BASE_C,   BASE_G,  BASE_T,   BASE_N,   BASE_NN}},
    [LINK_SYMBOL] = std::vector<AlignmentSymbol>{{LINK_M,   LINK_I1,  LINK_I2, LINK_I3P, LINK_D1,  LINK_D2,  LINK_D3P, LINK_NN}}
}};

const AlignmentSymbol SYMBOL_TYPE_TO_INCLU_END[NUM_SYMBOL_TYPES] = {
    [BASE_SYMBOL] = BASE_NN,
    [LINK_SYMBOL] = LINK_NN,
};

const AlignmentSymbol SYMBOL_TYPE_TO_AMBIG[NUM_SYMBOL_TYPES] = {
    [BASE_SYMBOL] = BASE_NN,
    [LINK_SYMBOL] = LINK_NN,
};

const std::array<SymbolType, 2> SYMBOL_TYPES_IN_VCF_ORDER = {{LINK_SYMBOL, BASE_SYMBOL}};

bool 
isSymbolSubstitution(AlignmentSymbol symbol) {
    return (SYMBOL_TYPE_TO_INCLU_BEG[BASE_SYMBOL] <= symbol && symbol <= SYMBOL_TYPE_TO_INCLU_END[BASE_SYMBOL]);
}

// is not WGS (i.e., is hybrid-capture-WES or amplicon-PCR)
bool
does_fmt_imply_short_frag(const auto & fmt) {
    return (fmt.APLRI[0] + fmt.APLRI[2]) < int64mul(fmt.APLRI[1] + fmt.APLRI[3], 300);
}

struct PhredMutationTable {
    const uvc1_qual_t transition_CG_TA;
    const uvc1_qual_t transition_AT_CG;
    const uvc1_qual_t transversion_CG_AT;
    const uvc1_qual_t transversion_other;
    const uvc1_qual_t indel_open;
    const uvc1_qual_t indel_ext;
    PhredMutationTable(
            const uvc1_qual_t transition_CG_TA,
            const uvc1_qual_t transition_AT_CG,
            const uvc1_qual_t transversion_CG_AT,
            const uvc1_qual_t transversion_other,
            const uvc1_qual_t indel_open,
            const uvc1_qual_t indel_ext)
            : 
            transition_CG_TA(transition_CG_TA),
            transition_AT_CG(transition_AT_CG),
            transversion_CG_AT(transversion_CG_AT),
            transversion_other(transversion_other),
            indel_open(indel_open),
            indel_ext(indel_ext)
            {
    }
    uvc1_qual_t toPhredErrRate(const AlignmentSymbol con_symbol, const AlignmentSymbol alt_symbol) const {
        if (isSymbolIns(con_symbol) || isSymbolDel(con_symbol)) {
            return indel_open;
        } else if (con_symbol == LINK_M) {
            if (LINK_D1 == alt_symbol || LINK_I1 == alt_symbol) {
                return indel_open;
            } else if (LINK_D2== alt_symbol || LINK_I2 == alt_symbol) {
                return indel_open + indel_ext * 1;
            } else {
                return indel_open + indel_ext * 2;
            }
        } else if ((con_symbol == BASE_C && alt_symbol == BASE_T) || (con_symbol == BASE_G && alt_symbol == BASE_A)) {
            return transition_CG_TA;
        } else if ((con_symbol == BASE_A && alt_symbol == BASE_G) || (con_symbol == BASE_T && alt_symbol == BASE_C)) {
            return transition_AT_CG;
        } else if ((con_symbol == BASE_C && alt_symbol == BASE_A) || (con_symbol == BASE_G && alt_symbol == BASE_T)) {
            return transversion_CG_AT;
        } else {
            return transversion_other;
        }
    }
};

constexpr bool 
is_mut_transition(const AlignmentSymbol con_symbol, const AlignmentSymbol alt_symbol) {
    return ((con_symbol == BASE_C && alt_symbol == BASE_T) || (con_symbol == BASE_G && alt_symbol == BASE_A)
        || (con_symbol == BASE_T && alt_symbol == BASE_C) || (con_symbol == BASE_A && alt_symbol == BASE_G)
    );
}

const uvc1_refgpos_t SYMBOL_TO_INDEL_N_UNITS[] = {
    [BASE_A] = 0, [BASE_C] = 0, [BASE_G] = 0, [BASE_T] = 0, [BASE_N] = 0,
    [BASE_NN] = 0, 
    [LINK_M] = 0, 
    [LINK_D3P] = -3, [LINK_D2] = -2, [LINK_D1] = -1,
    [LINK_I3P] = 3, [LINK_I2] = 2, [LINK_I1] = 1,
    [LINK_NN] = 0,
    [END_ALIGNMENT_SYMBOLS] = 0,
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

template <class T>
class TDistribution {
protected:
    std::array<T, NUM_ALIGNMENT_SYMBOLS> symbol2data;
};

typedef uvc1_readnum_t molcount_t;

template <class TB2C>
class GenericSymbol2Bucket2Count : TDistribution<TB2C> {
public:    
    molcount_t
    getSymbolBucketCount(AlignmentSymbol symbol, size_t bucket) const {
        assert(bucket < NUM_BUCKETS);
        return this->symbol2data[symbol][bucket];
    };
    
    const TB2C &
    getSymbolCounts(AlignmentSymbol symbol) const {
        return this->symbol2data[symbol];
    };

    void
    incSymbolBucketCount(AlignmentSymbol symbol, size_t bucket, uvc1_qual_t increment) {
        assert(bucket < NUM_BUCKETS || !fprintf(stderr, "%lu < %d failed!", bucket, NUM_BUCKETS));
        this->symbol2data[symbol].at(bucket) += increment;
    };
    
    const TB2C
    vectorsumBySymbolType(const SymbolType symboltype) const {
        TB2C ret = {0};
        for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symboltype]) {
            for (size_t i = 0; i < this->symbol2data[symbol].size(); i++) {
                ret[i] += this->getSymbolBucketCount(symbol, i);
            }
        }
        return ret;
    };

    void
    clearSymbolBucketCount() {
        assert(sizeof(TB2C) * NUM_ALIGNMENT_SYMBOLS == sizeof(this->symbol2data) || !fprintf(stderr, "%lu * %u != %lu\n", sizeof(TB2C), NUM_ALIGNMENT_SYMBOLS, sizeof(this->symbol2data)));
        for (size_t i = 0; i < NUM_ALIGNMENT_SYMBOLS; i++) {
            for (size_t j = 0; j < this->symbol2data[0].size(); j++) {
                this->symbol2data[i][j] = 0;
            }
        }
    };
};

template <class TInteger>
class GenericSymbol2Count : TDistribution<TInteger> {
public:
    TInteger
    getSymbolCount(AlignmentSymbol symbol) const {
        return this->symbol2data[symbol];
    };
    
    template <ValueType TUpdateType = SYMBOL_COUNT_SUM>
    int // update_max_inc : high GC : 3, even distribution of nucleotides : 6, conservative : 0
    incSymbolCount(const AlignmentSymbol symbol, const TInteger increment, const TInteger update_max_inc = 0) {
        static_assert(BASE_QUALITY_MAX == TUpdateType || SYMBOL_COUNT_SUM == TUpdateType);
        if (SYMBOL_COUNT_SUM == TUpdateType) {
            this->symbol2data[symbol] += increment;
        } else {
            this->symbol2data[symbol] = MAX(this->symbol2data[symbol], increment) + (this->symbol2data[symbol] ? update_max_inc : 0);
        }
        return 0;
    };
    
    const TInteger
    _sumBySymbolType(AlignmentSymbol beg, AlignmentSymbol end) const {
        assert (beg <= end);
        
        TInteger alpha_sum = 0;
        for (AlignmentSymbol symb = beg; symb <= end; symb = AlignmentSymbol(((uvc1_unsigned_int_t)symb) + 1)) {
            alpha_sum += this->symbol2data[symb];
        }
        return alpha_sum;
    };
    
    const TInteger
    sumBySymbolType(const SymbolType symboltype) const {
        if (symboltype == BASE_SYMBOL) {
            return this->_sumBySymbolType(BASE_A, BASE_NN);
        } else if (symboltype == LINK_SYMBOL) {
            return this->_sumBySymbolType(LINK_M, LINK_NN);
        } else {
            abort();
            return -1;
        }
    };
    
    template <bool TIndelIsMajor>
    int
    _fillConsensusCounts(
            AlignmentSymbol & count_argmax, uvc1_qual_t & count_max, uvc1_qual_t & count_sum,
            AlignmentSymbol incluBeg, AlignmentSymbol incluEnd) const {
        assert (incluBeg <= incluEnd);
        
        count_argmax = incluEnd; // assign with the value AlignmentSymbol(NUM_ALIGNMENT_SYMBOLS) to flag for error
        count_max = 0;
        count_sum = 0;
        for (AlignmentSymbol symb = incluBeg; symb <= incluEnd; symb = AlignmentSymbol(((uvc1_unsigned_int_t)symb) + 1)) {
            if (TIndelIsMajor) {
                if (count_max < this->symbol2data[symb] || (LINK_M == count_argmax && (0 < this->symbol2data[symb]))) {
                    count_argmax = symb;
                    count_max = this->symbol2data[symb];
                    count_sum = count_max;
                }
            } else {
                if (count_max < this->symbol2data[symb]) {
                    count_argmax = symb;
                    count_max = this->symbol2data[symb];
                }
                count_sum += this->symbol2data[symb];
            }
        }
        
        assert((incluBeg <= count_argmax && count_argmax <= incluEnd) || !fprintf(stderr, "The value %u is not between %u and %u", count_argmax, incluBeg, incluEnd));
        return 0;
    };
    
    template <bool TIndelIsMajor = false>
    int
    fillConsensusCounts(
            AlignmentSymbol & count_argmax, uvc1_qual_t & count_max, uvc1_qual_t & count_sum,
            const SymbolType symboltype) const {
        if (symboltype == BASE_SYMBOL) {
            return this->template _fillConsensusCounts<false        >(count_argmax, count_max, count_sum, BASE_A, BASE_NN);
        } else if (symboltype == LINK_SYMBOL) {
            return this->template _fillConsensusCounts<TIndelIsMajor>(count_argmax, count_max, count_sum, LINK_M, LINK_NN);
        } else {
            abort();
            return -1;
        }
    };
    
    template<bool TIndelIsMajor>
    AlignmentSymbol
    _updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count,
            const SymbolType symboltype, const AlignmentSymbol ambig_pos, const uvc1_qual_t incvalue) {
        AlignmentSymbol argmax_count = END_ALIGNMENT_SYMBOLS; // AlignmentSymbol(0) is not fully correct
        uvc1_qual_t max_count = 0;
        uvc1_qual_t sum_count = 0;
        thatSymbol2Count.template fillConsensusCounts<TIndelIsMajor>(argmax_count, max_count, sum_count, symboltype);
        
        if (max_count > 0) {
            if ((sum_count - max_count) == 0) {
                this->symbol2data[argmax_count] += incvalue;
                return argmax_count;
            } else {
                this->symbol2data[ambig_pos] += incvalue;
                return ambig_pos;
            }
        } else {
            return END_ALIGNMENT_SYMBOLS;
        }
    };

    template<bool TIndelIsMajor = false>
    std::array<AlignmentSymbol, 2>
    updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count, uvc1_qual_t incvalue = 1) {
        AlignmentSymbol baseSymb = this->template _updateByConsensus<false        >(thatSymbol2Count, BASE_SYMBOL, BASE_NN, incvalue);
        AlignmentSymbol linkSymb = this->template _updateByConsensus<TIndelIsMajor>(thatSymbol2Count, LINK_SYMBOL, LINK_NN, incvalue);
        return {baseSymb, linkSymb};
    };
    
    template <bool TIsIncVariable = true>
    AlignmentSymbol
    updateByRepresentative(const GenericSymbol2Count<TInteger> & other, uvc1_qual_t incvalue = 1) {
        AlignmentSymbol consalpha; 
        uvc1_qual_t countalpha, totalalpha;
        for (SymbolType symboltype : SYMBOL_TYPE_ARR) { 
            other.fillConsensusCounts(consalpha, countalpha, totalalpha, symboltype);
            if (countalpha > 0) {
                this->symbol2data[consalpha] += (TIsIncVariable ? totalalpha : incvalue);
            }
        }
        return consalpha;
    };
    
    int
    updateByFiltering(
            std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> & con_symbols, 
            const GenericSymbol2Count<TInteger> & other, 
            const std::array<uvc1_qual_t, NUM_SYMBOL_TYPES> thres = {{0}},
            const uvc1_qual_t incvalue = 1) {
        int ret = 0;
        AlignmentSymbol consalpha;
        uvc1_qual_t countalpha, totalalpha;
        for (SymbolType symboltype : SYMBOL_TYPE_ARR) {
            if (LINK_SYMBOL == symboltype) {
                other.template fillConsensusCounts<true >(consalpha, countalpha, totalalpha, symboltype);
            } else {
                other.template fillConsensusCounts<false>(consalpha, countalpha, totalalpha, symboltype);
            }
            auto adjcount = MAX(countalpha * 2, totalalpha) - totalalpha;
            if (adjcount >= thres[symboltype] && adjcount > 0) {
                this->symbol2data[consalpha] += incvalue;
                ret++;
            }
            con_symbols[symboltype] = consalpha;
        }
        return ret;
    };
    
    int
    updateByMajorMinusMinor(
            std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> & con_symbols,
            std::array<uvc1_qual_t, NUM_SYMBOL_TYPES> & con_counts,
            const GenericSymbol2Count<TInteger> & other) {
        int ret = 0;
        AlignmentSymbol consalpha;
        uvc1_qual_t countalpha, totalalpha;
        for (SymbolType symboltype : SYMBOL_TYPE_ARR) {
            if (LINK_SYMBOL == symboltype) {
                other.template fillConsensusCounts<true >(consalpha, countalpha, totalalpha, symboltype);
            } else {
                other.template fillConsensusCounts<false>(consalpha, countalpha, totalalpha, symboltype);
            }
            auto adjcount = MAX(countalpha * 2, totalalpha) - totalalpha;
            if (adjcount > 0) {
                this->symbol2data[consalpha] += adjcount;
                ret++;
            }
            con_symbols[symboltype] = consalpha;
            con_counts[symboltype] = adjcount;
        }
        return ret;
    };
};

template<class T>
class CoveredRegion {
         
protected:
    
    std::vector<T> idx2symbol2data;
    std::array<std::map<uvc1_refgpos_t, std::map<uvc1_readpos_t     , uvc1_readnum_t>>, NUM_INS_SYMBOLS> pos2dlen2data;
    std::array<std::map<uvc1_refgpos_t, std::map<std::string, uvc1_readnum_t>>, NUM_DEL_SYMBOLS> pos2iseq2data;
public:
    
    const uvc1_refgpos_t tid;
    const uvc1_refgpos_t incluBegPosition; // where end_pos = incluBegPosition + idx2symbol2data.size()
    
    CoveredRegion() {};
    CoveredRegion(uvc1_refgpos_t tid, uvc1_refgpos_t beg, uvc1_refgpos_t end): tid(tid), incluBegPosition(beg)  {
        assert (beg < end || !fprintf(stderr, "assertion %d < %d failed!\n", beg, end));
        this->idx2symbol2data = std::vector<T>(end-beg); // end should be real end position plus one
    };
    
    T &
    getRefByPos(const uvc1_refgpos_t pos, const bam1_t *bam = NULL) {
        assert(pos >= this->incluBegPosition || !fprintf(stderr, "%d >= %d failed for qname %s !\n", pos, this->incluBegPosition, (NULL != bam ? bam_get_qname(bam) : "?")));
        uvc1_refgpos_t pos2 = pos - this->incluBegPosition;
        assert(pos2 < UNSIGN2SIGN(idx2symbol2data.size()) || !fprintf(stderr, "%d  < %d failed for qname %s !\n", pos, this->incluBegPosition + (uvc1_refgpos_t)(idx2symbol2data.size()), (NULL != bam ? bam_get_qname(bam) : "?")));
        return this->idx2symbol2data[pos2];
    };
    
    const T &
    getByPos(const uvc1_refgpos_t pos, const bam1_t *bam = NULL) const {
        assert(pos >= this->incluBegPosition || !fprintf(stderr, "%d >= %d failed for qname %s !\n", pos, this->incluBegPosition, (NULL != bam ? bam_get_qname(bam) : "?")));
        uvc1_refgpos_t pos2 = pos - this->incluBegPosition;
        assert((pos2 < UNSIGN2SIGN(idx2symbol2data.size())) || !fprintf(stderr, "%d  < %d failed for qname %s !\n", pos, this->incluBegPosition + (uvc1_refgpos_t)(idx2symbol2data.size()), (NULL != bam ? bam_get_qname(bam) : "?")));
        return this->idx2symbol2data[pos2];
    };
    
    uvc1_refgpos_t // size_t 
    getIncluBegPosition() const {
        return this->incluBegPosition;
    };
    uvc1_refgpos_t // size_t 
    getExcluEndPosition() const {
        return this->incluBegPosition + UNSIGN2SIGN(idx2symbol2data.size());
    };
    
    const std::map<uvc1_refgpos_t, std::map<uvc1_refgpos_t    , uvc1_readnum_t>> & 
    getPosToDlenToData(const AlignmentSymbol s) const { 
        int idx = (LINK_D1 == s ? 0 : ((LINK_D2 == s) ? 1: 2));
        return pos2dlen2data[idx];
    };
    const std::map<uvc1_refgpos_t, std::map<std::string, uvc1_readnum_t>> & 
    getPosToIseqToData(const AlignmentSymbol s) const {
        int idx = (LINK_I1 == s ? 0 : ((LINK_I2 == s) ? 1: 2));
        return pos2iseq2data[idx];
    };
    
    std::map<uvc1_refgpos_t, std::map<uvc1_readpos_t   , uvc1_readnum_t>> & 
    getRefPosToDlenToData(const AlignmentSymbol s) {
        int idx = (LINK_D1 == s ? 0 : ((LINK_D2 == s) ? 1: 2));
        return pos2dlen2data[idx];
    };
    std::map<uvc1_refgpos_t, std::map<std::string, uvc1_refgpos_t>> & 
    getRefPosToIseqToData(const AlignmentSymbol s) {
        int idx = (LINK_I1 == s ? 0 : ((LINK_I2 == s) ? 1: 2));
        return pos2iseq2data[idx];
    };
};

template <class TSymbol2Bucket2Count>
class GenericSymbol2Bucket2CountCoverage : public CoveredRegion<TSymbol2Bucket2Count> {
    public:
    GenericSymbol2Bucket2CountCoverage() : CoveredRegion<TSymbol2Bucket2Count>(0, 0, 1) { }
    GenericSymbol2Bucket2CountCoverage(auto tid, auto beg, auto end) : CoveredRegion<TSymbol2Bucket2Count>(tid, beg, end) {}
};

typedef std::array<molcount_t, NUM_BUCKETS> Bucket2Count;

typedef GenericSymbol2Bucket2Count<Bucket2Count> Symbol2Bucket2Count;

typedef GenericSymbol2Count<uvc1_readnum_t> Symbol2Count;

typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2Count> Symbol2Bucket2CountCoverage;

typedef CoveredRegion<SegFormatPrepSet> SegFormatPrepSets;
typedef CoveredRegion<SegFormatThresSet> SegFormatThresSets;
typedef CoveredRegion<std::array<SegFormatInfoSet, NUM_ALIGNMENT_SYMBOLS>> Symbol2SegFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_FRAG_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2FragFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_FAM_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2FamFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_DUPLEX_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2DuplexFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_VQ_FORMAT_TAG_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2VQFormatTagSets;

int
formatSumBySymbolType(
        const auto & xFormatYSets,
        const auto symboltype,
        const auto formatSet) {
    int ret = 0;
    for (const auto symbol : SYMBOL_TYPE_TO_SYMBOLS[symboltype]) {
        ret += xFormatYSets[symbol][formatSet];
    }
    return ret;
};

void
initTidBegEnd(uvc1_refgpos_t & tid, uvc1_refgpos_t & inc_beg, uvc1_refgpos_t & exc_end) {
   tid = -1;
   inc_beg = INT32_MAX;
   exc_end = 0;
};

int 
fillTidBegEndFromAlns1(uvc1_refgpos_t & tid, uvc1_refgpos_t & inc_beg, uvc1_refgpos_t & exc_end, const std::vector<bam1_t *> & alns1, bool initialized=false) {
    assert(alns1.size() > 0);
    if (!initialized) {
        initTidBegEnd(tid, inc_beg, exc_end);
    }
    for (const bam1_t *aln : alns1) {
        assert (tid == -1 || SIGN2UNSIGN(aln->core.tid) == tid);
        tid = aln->core.tid;
        inc_beg = MIN(inc_beg, SIGN2UNSIGN(aln->core.pos));
        exc_end = MAX(exc_end, SIGN2UNSIGN(bam_endpos(aln))) + 1; // accounts for insertion at the end 
    }
    assert (tid != -1);
    assert (inc_beg < exc_end);
    return 0;
};

int 
fillTidBegEndFromAlns2(uvc1_refgpos_t  & tid, uvc1_refgpos_t & inc_beg, uvc1_refgpos_t & exc_end, const std::vector<std::vector<bam1_t *>> &alns2, bool initialized=false) {
    assert(alns2.size() > 0);
    if (!initialized) {
        initTidBegEnd(tid, inc_beg, exc_end);
    }
    for (auto alns1 : alns2) {
        fillTidBegEndFromAlns1(tid, inc_beg, exc_end, alns1, true);
    }
    return 0;
};

int 
fillTidBegEndFromAlns3(uvc1_refgpos_t & tid, uvc1_refgpos_t & inc_beg, uvc1_refgpos_t & exc_end, const std::vector<std::vector<std::vector<bam1_t *>>> & alns3, bool initialized=false) {
    assert(alns3.size() > 0);
    if (!initialized) {
        initTidBegEnd(tid, inc_beg, exc_end);
    }
    for (auto alns2 : alns3) {
        fillTidBegEndFromAlns2(tid, inc_beg, exc_end, alns2, true);
    }
    return 0;
};

bool 
is_indel_context_more_STR(uvc1_refgpos_t rulen1, uvc1_refgpos_t rc1, uvc1_refgpos_t rulen2, uvc1_refgpos_t rc2, auto indel_str_repeatsize_max) {
    if (rulen2 * rc2 == 0) {
        return true;
    }
    if (rulen1 > indel_str_repeatsize_max || rulen2 > indel_str_repeatsize_max) {
        return ((rulen1 < rulen2 || (rulen1 == rulen2 && rc1 > rc2))  ? true : false);
    }
    int rank1 = (rc1 <= 1 ? (-(int)rc1 * rulen1) : ((int)(rc1 - 1) * rulen1));
    int rank2 = (rc2 <= 1 ? (-(int)rc2 * rulen1) : ((int)(rc2 - 1) * rulen2));
    if (0 == rc1 || 0 == rulen1) {
        rank1 = -100;
    }
    if (0 == rc2 || 0 == rulen2) {
        rank2 = -100;
    }
    if (rank1 > rank2) {
        return true;
    } else {
        return false;
    }
}

uvc1_refgpos_t 
indelpos_repeatsize_to_repeatnum(const auto & refstring, const auto refpos, const auto repeatsize){
    uvc1_refgpos_t qidx = refpos;
    while ((qidx + repeatsize < UNSIGN2SIGN(refstring.size())) && refstring[qidx] == refstring[qidx+repeatsize]) {
        qidx++;
    }
    return (qidx - refpos) / repeatsize + 1;
}

int 
indelpos_to_context(
        std::string & repeatunit, 
        uvc1_refgpos_t & max_repeatnum,
        const std::string & refstring, 
        uvc1_refgpos_t refpos, 
        uvc1_refgpos_t indel_str_repeatsize_max) {
    max_repeatnum = 0;
    if (refpos >= UNSIGN2SIGN(refstring.size())) {
        repeatunit = "";
        return -1;
    }
    uvc1_refgpos_t repeatsize_at_max_repeatnum = 0;
    for (uvc1_refgpos_t repeatsize = 1; repeatsize <= indel_str_repeatsize_max; repeatsize++) {
        uvc1_refgpos_t repeatnum = indelpos_repeatsize_to_repeatnum(refstring, refpos, repeatsize);
        if (is_indel_context_more_STR(repeatsize, repeatnum, repeatsize_at_max_repeatnum, max_repeatnum, indel_str_repeatsize_max)) {
            max_repeatnum = repeatnum;
            repeatsize_at_max_repeatnum = repeatsize;
        }
    }
    repeatunit = refstring.substr(refpos, repeatsize_at_max_repeatnum);
    return 0;
}

uvc1_qual_t
indel_len_rusize_phred(uvc1_refgpos_t indel_len, uvc1_refgpos_t repeatunit_size) {
    assert (indel_len > 0 && repeatunit_size > 0);
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2734402/#bib41 suggests power law with exponent of 1.5-1.6
    // TODO: investigate this
    // derived from the python code: for i in range(1, 20): print("{}  {}".format(  int( round(10.0/log(10.0)*log(i)) ) , i  )) 
    const std::array<uvc1_qual_t, 18+1> n_units_to_phred = {{
        0 , 
        0 , // 1
        3 , // 2
        5 , // 3
        6 , // 4
        7 , // 5
        8 , // 6
        8 , // 7
        9 , // 8
        10 , // 9
        10 , // 10
        10 , // 11
        11 , // 12
        11 , // 13
        11 , // 14
        12 , // 15
        12 , // 16
        12 , // 17
        13 , // 18
    }};
    if (0 == (indel_len % repeatunit_size)) {
        auto n_units = indel_len / repeatunit_size;
        return n_units_to_phred[MIN(n_units  , UNSIGN2SIGN(n_units_to_phred.size()-1))];
    } else {
        return n_units_to_phred[MIN(indel_len, UNSIGN2SIGN(n_units_to_phred.size()-1))];
    }
}

// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/ for artifactual errors
// https://www.cell.com/trends/genetics/fulltext/S0168-9525(13)00070-X, ref12, for germline variants
uvc1_qual_t
indel_phred(double ampfact, uvc1_refgpos_t repeatsize_at_max_repeatnum, uvc1_refgpos_t max_repeatnum) {
    uvc1_refgpos_t region_size = repeatsize_at_max_repeatnum * max_repeatnum;
    double num_slips = (region_size > 64 ? (double)(region_size - 8) : log1p(exp((double)region_size - (double)8))) 
            * ampfact / ((double)(repeatsize_at_max_repeatnum * repeatsize_at_max_repeatnum));
    return prob2phred((1.0 - DBL_EPSILON) / (num_slips + 1.0));
    // e.g. for the short tandem repeat ACACAC : repeatsize_at_max_repeatnum = 2, indel_n_units = 3, repeatunit = AC
}

std::vector<RegionalTandemRepeat>
refstring2repeatvec(
        const auto & refstring,
        const uvc1_refgpos_t indel_str_repeatsize_max,
        const uvc1_refgpos_t indel_minisattelite_repeatsize_max,
        const uvc1_qual_t indel_BQ_max,
        const double indel_polymerase_slip_rate,
        const double indel_del_to_ins_err_ratio,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    std::vector<RegionalTandemRepeat> region_repeatvec(refstring.size());
    for (auto & rtr : region_repeatvec) {
        rtr.indelphred = indel_BQ_max;
    }
    for (uvc1_refgpos_t refpos = 0; refpos < UNSIGN2SIGN(refstring.size());) {
        uvc1_refgpos_t repeatsize_at_max_repeatnum = 0;
        uvc1_refgpos_t max_repeatnum = 0;
        uvc1_refgpos_t repeat_endpos = refpos;
        
        uvc1_refgpos_t anyTR_repeatsize_at_max_repeatnum = 0;
        uvc1_refgpos_t anyTR_max_repeatnum = 0;
        uvc1_refgpos_t anyTR_repeat_endpos = refpos;
        
        for (uvc1_refgpos_t repeatsize = 1; repeatsize <= indel_minisattelite_repeatsize_max; repeatsize++) {
            
            uvc1_refgpos_t qidx = refpos;
            while (qidx + repeatsize < UNSIGN2SIGN(refstring.size()) && refstring[qidx] == refstring[qidx+repeatsize]) {
                qidx++;
            }
            
            uvc1_refgpos_t repeatnum = (qidx - refpos) / repeatsize + 1;
            if (repeatsize <= indel_str_repeatsize_max 
                    && is_indel_context_more_STR(repeatsize, repeatnum, repeatsize_at_max_repeatnum, max_repeatnum, indel_str_repeatsize_max)) {
                repeatsize_at_max_repeatnum = repeatsize;
                max_repeatnum = repeatnum;
                repeat_endpos = qidx + repeatsize;
            }
            if (is_indel_context_more_STR(repeatsize, repeatnum, anyTR_repeatsize_at_max_repeatnum, anyTR_max_repeatnum, indel_minisattelite_repeatsize_max)) {
                anyTR_repeatsize_at_max_repeatnum = repeatsize;
                anyTR_max_repeatnum = repeatnum;
                anyTR_repeat_endpos = qidx + repeatsize;
            }
        }
        assert(repeat_endpos > refpos);
        {
            uvc1_refgpos_t tl = MIN(repeat_endpos, UNSIGN2SIGN(refstring.size())) - refpos;
            const uvc1_qual_t decphred = indel_phred(indel_polymerase_slip_rate * (indel_del_to_ins_err_ratio), 
                    repeatsize_at_max_repeatnum, tl / repeatsize_at_max_repeatnum);
            for (uvc1_refgpos_t i = refpos; i != MIN(repeat_endpos, UNSIGN2SIGN(refstring.size())); i++) {
                if (tl > region_repeatvec[i].tracklen) {
                    region_repeatvec[i].begpos = refpos;
                    region_repeatvec[i].tracklen = tl;
                    region_repeatvec[i].unitlen = repeatsize_at_max_repeatnum;
                    region_repeatvec[i].indelphred = indel_BQ_max - MIN(indel_BQ_max - 1, decphred);
                }
            }
        }
        {
            uvc1_refgpos_t anyTR_tl = MIN(anyTR_repeat_endpos, UNSIGN2SIGN(refstring.size())) - refpos;
            for (uvc1_refgpos_t i = refpos; i != MIN(anyTR_repeat_endpos, UNSIGN2SIGN(refstring.size())); i++) {
                if (anyTR_tl > region_repeatvec[i].anyTR_tracklen) {
                    region_repeatvec[i].anyTR_begpos = refpos;
                    region_repeatvec[i].anyTR_tracklen = anyTR_tl;
                    region_repeatvec[i].anyTR_unitlen = anyTR_repeatsize_at_max_repeatnum;
                }
            }
        }
        const auto nbases_to_next = indel_str_repeatsize_max + repeatsize_at_max_repeatnum;
        refpos += MAX(repeatsize_at_max_repeatnum * max_repeatnum, nbases_to_next + 1) - (nbases_to_next);
    }
    region_repeatvec.push_back(LAST(region_repeatvec));
    return region_repeatvec;
}

uvc1_qual_t
ref_to_phredvalue(uvc1_refgpos_t & n_units, 
        const auto & refstring, 
        const uvc1_refgpos_t refpos, 
        const uvc1_qual_t max_phred, 
        double ampfact, 
        const uvc1_refgpos_t cigar_oplen, 
        const auto cigar_op, 
        const uvc1_refgpos_t indel_str_repeatsize_max,
        const double indel_del_to_ins_err_ratio,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {

    uvc1_refgpos_t max_repeatnum = 0;
    uvc1_refgpos_t repeatsize_at_max_repeatnum = 0;
    for (uvc1_refgpos_t repeatsize = 1; repeatsize <= indel_str_repeatsize_max; repeatsize++) {
        uvc1_refgpos_t qidx = refpos;
        while (qidx + repeatsize < UNSIGN2SIGN(refstring.size()) && refstring[qidx] == refstring[qidx+repeatsize]) {
            qidx++;
        }
        uvc1_refgpos_t repeatnum = (qidx - refpos) / repeatsize + 1;
        if (is_indel_context_more_STR(repeatsize, repeatnum, repeatsize_at_max_repeatnum, max_repeatnum, indel_str_repeatsize_max)) {
            max_repeatnum = repeatnum;
            repeatsize_at_max_repeatnum = repeatsize;
        }
    }
    if (cigar_oplen == repeatsize_at_max_repeatnum && cigar_op == BAM_CDEL) {
        ampfact *= indel_del_to_ins_err_ratio;
    }
    // Because of a higher number of PCR cycles, it is higher than the one set in Fig. 3 at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/
    uvc1_qual_t decphred = indel_phred(ampfact, repeatsize_at_max_repeatnum, max_repeatnum);
    if (0x1 == (specialflag & 0x1)) { 
        LOG(logINFO) << "indel_phred(" 
                << ampfact << ", " 
                << repeatsize_at_max_repeatnum << ", "
                << max_repeatnum << ") = " << decphred; 
    }
    // The number of units should minimize InDel collision.
    if (repeatsize_at_max_repeatnum * (max_repeatnum - 1) >= 6 - 1) {
        n_units = ((0 == cigar_oplen % repeatsize_at_max_repeatnum) ? (cigar_oplen / repeatsize_at_max_repeatnum) : ((1 == cigar_oplen) ? 1 : 0)); 
    } else {
        n_units = 1 + (cigar_oplen / 6);
    }
    return max_phred - MIN(max_phred, decphred) + indel_len_rusize_phred(cigar_oplen, repeatsize_at_max_repeatnum); 
}

int
update_seg_format_prep_sets_by_aln(
        SegFormatPrepSets & seg_format_prep_sets,
        const bam1_t *aln,
        const std::vector<RegionalTandemRepeat> & rtr_vec,
        const CoveredRegion<uvc1_qual_big_t> & baq_offsetarr,
        const uvc1_refgpos_t region_offset,
        const auto dflag,
        const std::basic_string<AlignmentSymbol> & region_symbolvec,
        
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    const uvc1_refgpos_t rend = bam_endpos(aln);
    const auto cigar = bam_get_cigar(aln);
    uvc1_refgpos_t nge_cnt = 0; // number of gap extensions.
    uvc1_refgpos_t ngo_cnt = 0; // number of gap openings.
    
    uvc1_qual_t insbaq_sum = 0;
    uvc1_qual_t delbaq_sum = 0;
    uvc1_readpos_t inslen_sum = 0;
    uvc1_readpos_t dellen_sum = 0;
    
    uvc1_refgpos_t qpos = 0;
    uvc1_refgpos_t rpos = aln->core.pos;
    
    uvc1_readpos_t max_clip_len = 0;
    uvc1_readpos_t max_indel_len = 0;
    for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
        const auto c = cigar[i];
        const auto cigar_op = bam_cigar_op(c);
        const auto cigar_oplen = bam_cigar_oplen(c);
        if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) { 
            nge_cnt += bam_cigar_oplen(c);
            ngo_cnt++;
            if (BAM_CINS == cigar_op) {
                insbaq_sum += baq_offsetarr.getByPos(MIN(rpos + UNSIGN2SIGN(cigar_oplen), baq_offsetarr.getExcluEndPosition() - 1)) - baq_offsetarr.getByPos(rpos);
                inslen_sum += bam_cigar_oplen(c);
                qpos += cigar_oplen;
            } else if (BAM_CDEL == cigar_op) {
                delbaq_sum += baq_offsetarr.getByPos(MIN(rpos + UNSIGN2SIGN(cigar_oplen), baq_offsetarr.getExcluEndPosition() - 1)) - baq_offsetarr.getByPos(rpos);
                dellen_sum += bam_cigar_oplen(c);
                rpos += cigar_oplen;
            }
        } else if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            UPDATE_MAX(max_indel_len, UNSIGN2SIGN(cigar_oplen));
            qpos += cigar_oplen;
            rpos += cigar_oplen;
        } else {
            if (cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
                UPDATE_MAX(max_clip_len, UNSIGN2SIGN(cigar_oplen));
            }
            process_cigar(qpos, rpos, cigar_op, cigar_oplen);
        }
    }
    
    const auto *bam_aux_data = bam_aux_get(aln, "NM");
    const uvc1_refgpos_t nm_cnt = ((bam_aux_data != NULL) ? bam_aux2i(bam_aux_data) : nge_cnt);
    assert (nm_cnt >= nge_cnt);
    const uvc1_base1500x_t xm_cnt = nm_cnt - nge_cnt;
    const uvc1_base1500x_t xm1500 = xm_cnt * 1500 / (rend - aln->core.pos);
    const uvc1_base1500x_t go1500 = ngo_cnt * 1500 / (rend - aln->core.pos);
    const uvc1_base1500x_t avg_gaplen = nge_cnt / MAX(1, ngo_cnt);
    
    const uvc1_refgpos_t frag_pos_L = MIN(aln->core.pos, aln->core.mpos); 
    const uvc1_refgpos_t frag_pos_R = (frag_pos_L + abs(aln->core.isize));
    const bool isrc = ((aln->core.flag & 0x10) == 0x10); 
    
    const auto pcr_dp_inc = ((dflag & 0x4) ? 1 : 0);
    qpos = 0;
    rpos = aln->core.pos;
    
    for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
        const auto c = cigar[i];
        const auto cigar_op = bam_cigar_op(c);
        const auto cigar_oplen = bam_cigar_oplen(c);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            for (uint32_t j = 0; j < cigar_oplen; j++) {
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_pcr_dp += pcr_dp_inc;
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_dp += 1;
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_qlen += (rend - aln->core.pos);
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_XM1500 += xm1500;
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_GO1500 += go1500;
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_GAPLEN += avg_gaplen;

                if (aln->core.isize != 0) {
                    if (isrc) { 
                        seg_format_prep_sets.getRefByPos(rpos).segprep_a_LI += MIN(rpos - frag_pos_L + 1, MAX_INSERT_SIZE);
                        seg_format_prep_sets.getRefByPos(rpos).segprep_a_LIDP += 1;
                    } else { 
                        seg_format_prep_sets.getRefByPos(rpos).segprep_a_RI += MIN(frag_pos_R - rpos    , MAX_INSERT_SIZE);
                        seg_format_prep_sets.getRefByPos(rpos).segprep_a_RIDP += 1;
                    }
                }
                auto refsymbol = BASE_NN;
                auto readsymbol = END_ALIGNMENT_SYMBOLS;
                auto next_qpos = qpos;
                auto next_rpos = rpos;
                while (refsymbol != readsymbol && next_qpos < aln->core.l_qseq && next_rpos < rend) {
                    const auto base4bit = bam_seqi(bam_get_seq(aln), next_qpos);
                    const auto base3bit = seq_nt16_int[base4bit];
                    refsymbol = region_symbolvec[next_rpos - region_offset];
                    readsymbol = AlignmentSymbol(base3bit);
                    next_qpos++;
                    next_rpos++;
                } 
                if (next_rpos == rpos + 2) {
                    for (auto r = MAX(aln->core.pos, rpos - 1); r < MIN(next_rpos, rend); r++) {
                        seg_format_prep_sets.getRefByPos(r).segprep_a_snv_dp += 1;
                    }
                }
                if (next_rpos > rpos + 2) {
                    for (auto r = MAX(aln->core.pos, rpos - 1); r < MIN(next_rpos, rend); r++) {
                        seg_format_prep_sets.getRefByPos(r).segprep_a_dnv_dp += 1;
                    }
                }
                if (BAM_PHREDI(aln, qpos) >= paramset.bias_thres_highBQ) {
                    uvc1_refgpos_t ldist = rpos - aln->core.pos + 1;
                    uvc1_refgpos_t rdist = rend - rpos;
                    
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_l_dist_sum += ldist;
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_r_dist_sum += rdist;
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_inslen_sum += inslen_sum;
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_dellen_sum += dellen_sum;
                    
                    uvc1_qual_t lbaq = baq_offsetarr.getByPos(rpos) - baq_offsetarr.getByPos(aln->core.pos) + 1;
                    uvc1_qual_t rbaq = baq_offsetarr.getByPos(rend-1) - baq_offsetarr.getByPos(rpos) + 1;
                    
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_l_BAQ_sum += lbaq;
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_r_BAQ_sum += rbaq;
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_insBAQ_sum += insbaq_sum;
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_delBAQ_sum += delbaq_sum;
                    
                    seg_format_prep_sets.getRefByPos(rpos).segprep_a_highBQ_dp += 1;
                }
                qpos++;
                rpos++;
            }
        } else if (cigar_op == BAM_CINS) {
            const auto & rtr1 = rtr_vec[MAX(paramset.indel_adj_tracklen_dist, rpos - region_offset) - paramset.indel_adj_tracklen_dist];
            const auto & rtr2 = rtr_vec[MIN(rpos - region_offset + paramset.indel_adj_tracklen_dist, UNSIGN2SIGN(rtr_vec.size()) - 1)];
            
            // this code has the potential to be useful at extremely high sequencing depth, which usually does not occur in practice
#if COMPILATION_TRY_HIGH_DEPTH_POS_BIAS
            const auto ins_rbeg = MAX(aln->core.pos, rpos - 100);
            const auto ins_rend = MIN(rend + UNSIGN2SIGN(cigar_oplen), rpos + 100);
            for (auto rpos2 = ins_rbeg + 1; rpos2 < ins_rend; rpos2++) {
                auto ldist = rpos2 - ins_rbeg;
                auto rdist = ins_rend - rpos2;
                auto rdist2 = MAX(1, ins_rend - (rpos2 + rtr2.tracklen));
                auto lweight = calc_indel_weight(UNSIGN2SIGN(cigar_oplen), ldist);
                auto rweight = calc_indel_weight(UNSIGN2SIGN(cigar_oplen), rdist2);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_l_ins_dist_x_wei += ldist * lweight;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_l_ins_weight += lweight;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_r_ins_dist_x_wei += rdist * rweight;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_r_ins_weight += rweight;
            }
#endif
            const auto unitlen2 = MAX(1, (rtr1.tracklen > rtr2.tracklen) ? rtr1.unitlen : rtr2.unitlen);
            const uvc1_refgpos_t nbases = (cigar_oplen * paramset.indel_adj_indellen_perc / 100);
            for (uvc1_refgpos_t rpos2 = MAX(rpos - nbases, aln->core.pos); rpos2 < MIN(rpos + nbases, rend); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_ins_dp += 1;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_ins_pow2len += mathsquare(cigar_oplen);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_ins_l_pow2len += mathsquare(rpos2 + 1 - (rpos - nbases));
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_ins_r_pow2len += mathsquare((rpos + nbases) - rpos2);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_ins_inv100len += 100 / ((0 == cigar_oplen % unitlen2) ? (cigar_oplen / unitlen2) : 4);
            }
            
            for (uvc1_refgpos_t rpos2 = MAX((region_offset + rtr1.begpos) - paramset.indel_adj_tracklen_dist, aln->core.pos); 
                    rpos2 < MIN((region_offset + rtr2.begpos + rtr2.tracklen) + paramset.indel_adj_tracklen_dist, rend); 
                    rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_RTR_ins_dp += 1;
            }
            seg_format_prep_sets.getRefByPos(rpos).segprep_a_at_ins_dp += 1;
            qpos += cigar_oplen;
        } else if (cigar_op == BAM_CDEL) {
            const auto & rtr1 = rtr_vec[MAX(paramset.indel_adj_tracklen_dist, rpos - region_offset) - paramset.indel_adj_tracklen_dist];
            const auto & rtr2 = rtr_vec[MIN(rpos - region_offset + paramset.indel_adj_tracklen_dist, UNSIGN2SIGN(rtr_vec.size()) - 1)];
            assert (rtr1.begpos <= rtr2.begpos || !fprintf(stderr, "AssertionError: %d <= %d failed for rtr1 and rtr2!\n", rtr1.begpos, rtr2.begpos));
            
            // this code has the potential to be useful at extremely high sequencing depth, which usually does not occur in practice
#if COMPILATION_TRY_HIGH_DEPTH_POS_BIAS
            const auto del_rbeg = MAX(aln->core.pos + UNSIGN2SIGN(cigar_oplen), rpos - 100);
            const auto del_rend = MIN(rend - UNSIGN2SIGN(cigar_oplen), rpos + 100);
            for (auto rpos2 = del_rbeg + 1; rpos2 < del_rend; rpos2++) {
                auto ldist = rpos2 - del_rbeg;
                auto rdist = del_rend - rpos2; 
                auto rdist2 = MAX(1, del_rend - (rpos2 + rtr2.tracklen));
                auto lweight = calc_indel_weight(UNSIGN2SIGN(cigar_oplen), ldist);
                auto rweight = calc_indel_weight(UNSIGN2SIGN(cigar_oplen), rdist2);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_l_del_dist_x_wei += ldist * lweight;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_l_del_weight += lweight;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_r_del_dist_x_wei += rdist * rweight;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_aa_r_del_weight += rweight;
            }
#endif
            for (uvc1_refgpos_t rpos2 = rpos; rpos2 < rpos + UNSIGN2SIGN(cigar_oplen); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_pcr_dp += pcr_dp_inc;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_dp += 1;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_qlen += (rend - aln->core.pos);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_highBQ_dp += 1;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_XM1500 += xm1500;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_GO1500 += go1500;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_GAPLEN += avg_gaplen;
                
                if (aln->core.isize != 0) {
                    if (isrc) { 
                        seg_format_prep_sets.getRefByPos(rpos2).segprep_a_LI += MIN(rpos - frag_pos_L + 1, MAX_INSERT_SIZE); 
                        seg_format_prep_sets.getRefByPos(rpos2).segprep_a_LIDP += 1;
                    } else { 
                        seg_format_prep_sets.getRefByPos(rpos2).segprep_a_RI += MIN(frag_pos_R - rpos    , MAX_INSERT_SIZE); 
                        seg_format_prep_sets.getRefByPos(rpos2).segprep_a_RIDP += 1;
                    }
                }
                uvc1_refgpos_t ldist = rpos - aln->core.pos + 1;
                uvc1_refgpos_t rdist = rend - rpos;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_l_dist_sum += ldist;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_r_dist_sum += rdist;
                
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_inslen_sum += inslen_sum;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_dellen_sum += dellen_sum;
                
                uvc1_refgpos_t lbaq = baq_offsetarr.getByPos(rpos) - baq_offsetarr.getByPos(aln->core.pos) + 1;
                uvc1_refgpos_t rbaq = baq_offsetarr.getByPos(rend - 1) - baq_offsetarr.getByPos(rpos) + 1;
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_l_BAQ_sum += lbaq;
                seg_format_prep_sets.getRefByPos(rpos).segprep_a_r_BAQ_sum += rbaq;
                
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_insBAQ_sum += insbaq_sum;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_delBAQ_sum += delbaq_sum;
            }
            
            const auto unitlen2 = MAX(1, (rtr1.tracklen > rtr2.tracklen) ? rtr1.unitlen : rtr2.unitlen);
            const uvc1_refgpos_t nbases_l = (cigar_oplen * (paramset.indel_adj_indellen_perc - 100) / 100);
            const uvc1_refgpos_t nbases_r = (cigar_oplen * paramset.indel_adj_indellen_perc / 100);
            const auto indelregion_lpos = MAX(rpos - nbases_l, aln->core.pos);
            const auto indelregion_rpos = MIN(rpos + nbases_r, rend) - 1;
            for (uvc1_refgpos_t rpos2 = indelregion_lpos; rpos2 <= indelregion_rpos; rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_del_dp += 1;
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_del_pow2len += mathsquare(cigar_oplen);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_del_l_pow2len += mathsquare(rpos2 - indelregion_lpos + 1);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_del_r_pow2len += mathsquare(indelregion_rpos - rpos2 + 1);
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_del_inv100len += 100 / ((0 == cigar_oplen % unitlen2) ? (cigar_oplen / unitlen2) : 4);
            }
            
            for (uvc1_refgpos_t rpos2 = MAX((region_offset + rtr1.begpos) - paramset.indel_adj_tracklen_dist, aln->core.pos);
                    rpos2 < MIN((region_offset + rtr2.begpos + rtr2.tracklen) + paramset.indel_adj_tracklen_dist, rend);
                    rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_RTR_del_dp += 1;
            }
            seg_format_prep_sets.getRefByPos(rpos).segprep_a_at_del_dp += 1;
            rpos += cigar_oplen;
        } else {
            const uvc1_refgpos_t rpos_delta = ((0 == i) ? (0) : (-1));
            if ((BAM_CSOFT_CLIP == cigar_op || BAM_CHARD_CLIP == cigar_op) && pcr_dp_inc) {
                for (auto rpos2 = rpos + rpos_delta - paramset.microadjust_near_clip_dist; rpos2 <= rpos + rpos_delta + paramset.microadjust_near_clip_dist; rpos2++) {
                    if (seg_format_prep_sets.getIncluBegPosition() <= rpos2 && rpos2 < seg_format_prep_sets.getExcluEndPosition()) {
                        seg_format_prep_sets.getRefByPos(rpos2).segprep_a_near_pcr_clip_dp += pcr_dp_inc;
                    }
                }
            }
            // if a clipped sequence is beyond the primer, then the clipped sequence should not be of biological origin anyway
            if ((BAM_CSOFT_CLIP == cigar_op || BAM_CHARD_CLIP == cigar_op) 
                    && (0 == pcr_dp_inc) 
                    && (UNSIGN2SIGN(cigar_oplen) >= paramset.microadjust_alignment_clip_min_len)) {
                seg_format_prep_sets.getRefByPos(rpos + rpos_delta).segprep_a_near_long_clip_dp += 1;
            }
            process_cigar(qpos, rpos, cigar_op, cigar_oplen);
        }
    }
    assert(bam_endpos(aln) == rpos || !fprintf(stderr, "%u == %u failed for bam %s at tid %d position %d", 
            bam_endpos(aln), rpos, bam_get_qname(aln), aln->core.tid, aln->core.pos)); 
    return 0;
}

int
update_seg_format_thres_from_prep_sets(
        auto & region_repeatvec,
        SegFormatThresSets & seg_format_thres_sets,
        const SegFormatPrepSets & seg_format_prep_sets,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    assert(seg_format_thres_sets.getIncluBegPosition() + UNSIGN2SIGN(region_repeatvec.size()) == seg_format_thres_sets.getExcluEndPosition() 
            || !fprintf(stderr, "%d + %lu == %d failed !\n", 
           seg_format_thres_sets.getIncluBegPosition(), region_repeatvec.size(), seg_format_thres_sets.getExcluEndPosition()));
    assert(seg_format_thres_sets.getIncluBegPosition() == seg_format_prep_sets.getIncluBegPosition());
    assert(seg_format_thres_sets.getExcluEndPosition() == seg_format_prep_sets.getExcluEndPosition());
    for (uvc1_refgpos_t epos = seg_format_prep_sets.getIncluBegPosition(); epos != seg_format_prep_sets.getExcluEndPosition(); epos++) {
#if COMPILATION_ENABLE_XMGOT
        const uvc1_readnum_t segprep_a_dp = MAX(seg_format_prep_sets.getByPos(epos).segprep_a_dp, 1);
#endif
        const auto & p = seg_format_prep_sets.getByPos(epos);
        auto & t = seg_format_thres_sets.getRefByPos(epos);
        
        auto segLIDP = MAX(seg_format_prep_sets.getByPos(epos).segprep_a_LIDP, 1);
        auto segRIDP = MAX(seg_format_prep_sets.getByPos(epos).segprep_a_RIDP, 1);
        const auto ins_border_l_len = ceil(sqrt((p.segprep_a_near_ins_l_pow2len) / MAX(p.segprep_a_near_ins_dp, 1)));
        const auto del_border_l_len = ceil(sqrt((p.segprep_a_near_del_l_pow2len) / MAX(p.segprep_a_near_del_dp, 1)));
        const auto ins_border_r_len = ceil(sqrt((p.segprep_a_near_ins_r_pow2len) / MAX(p.segprep_a_near_ins_dp, 1)));
        const auto del_border_r_len = ceil(sqrt((p.segprep_a_near_del_r_pow2len) / MAX(p.segprep_a_near_del_dp, 1)));
        
        const auto dnv_border_len = ((SEQUENCING_PLATFORM_IONTORRENT == paramset.inferred_sequencing_platform 
                && (p.segprep_a_dnv_dp * 2 > p.segprep_a_snv_dp)) ? (10) : 0);
        
        const auto max_border_l_len = MAX3(ins_border_l_len, del_border_l_len, dnv_border_len);
        const auto max_border_r_len = MAX3(ins_border_r_len, del_border_r_len, dnv_border_len);
        t.segthres_aLPxT = max_border_l_len + paramset.bias_thres_aLPxT_add;
        t.segthres_aRPxT = max_border_r_len + paramset.bias_thres_aLPxT_add;
        
        auto & rtr = region_repeatvec[epos - seg_format_prep_sets.getIncluBegPosition()];
        if (p.segprep_a_near_ins_dp * paramset.indel_del_to_ins_err_ratio < p.segprep_a_near_del_dp) {
            rtr.indelphred += (uvc1_qual_t)round(numstates2phred(paramset.indel_del_to_ins_err_ratio)) / 2;
        }
        if (p.segprep_a_near_del_dp * paramset.indel_del_to_ins_err_ratio < p.segprep_a_near_ins_dp) {
            rtr.indelphred -= (uvc1_qual_t)round(numstates2phred(paramset.indel_del_to_ins_err_ratio)) / 2;
        }
        const auto pc_inc1 = (uvc1_qual_t)(3 * 100 * MAX(1, p.segprep_a_near_ins_dp + p.segprep_a_near_del_dp)
                / (MAX(1, p.segprep_a_near_ins_inv100len + p.segprep_a_near_del_inv100len))) - 3;
        assert (epos >= seg_format_prep_sets.getIncluBegPosition() || !fprintf(stderr, "%d >= %d failed!", epos, seg_format_prep_sets.getIncluBegPosition()));
        assert (epos - seg_format_prep_sets.getIncluBegPosition() < UNSIGN2SIGN(region_repeatvec.size())
                || !fprintf(stderr, "%d - %d < %lu failed!", epos, seg_format_prep_sets.getIncluBegPosition(), region_repeatvec.size()));
        rtr.indelphred += BETWEEN(pc_inc1, 0, 6);
        UPDATE_MAX(rtr.indelphred, 0);
        
        const bool is_normal = (NOT_PROVIDED != paramset.vcf_tumor_fname);
        
#if COMPILATION_ENABLE_XMGOT
        const auto bias_thres_PFXM1T_perc = (is_normal ? paramset.bias_thres_PFXM1NT_perc : paramset.bias_thres_PFXM1T_perc);
        const auto bias_thres_PFGO1T_perc = (is_normal ? paramset.bias_thres_PFGO1NT_perc : paramset.bias_thres_PFGO1T_perc);
        
        t.segthres_aXM1T = int64mul(p.segprep_a_XM1500,          bias_thres_PFXM1T_perc) / (segprep_a_dp * 100) + paramset.bias_thres_PFXM1T_add; // easier to pass
        t.segthres_aXM2T = int64mul(p.segprep_a_XM1500, paramset.bias_thres_PFXM2T_perc) / (segprep_a_dp * 100) + paramset.bias_thres_PFXM2T_add; // the higher the stronger the bias
        t.segthres_aGO1T = int64mul(p.segprep_a_GO1500,          bias_thres_PFGO1T_perc) / (segprep_a_dp * 100) + paramset.bias_thres_PFGO1T_add; // easier to pass
        t.segthres_aGO2T = int64mul(p.segprep_a_GO1500, paramset.bias_thres_PFGO2T_perc) / (segprep_a_dp * 100) + paramset.bias_thres_PFGO2T_add; // the higher the stronger the bias
#endif
        
        const auto bias_thres_aLRI1T_perc = (is_normal ? paramset.bias_thres_aLRI1NT_perc : paramset.bias_thres_aLRI1T_perc);
        const auto bias_thres_aLRI1t_perc = (is_normal ? paramset.bias_thres_aLRI1Nt_perc : paramset.bias_thres_aLRI1t_perc);
        
        t.segthres_aLI1T = int64mul(p.segprep_a_LI,          bias_thres_aLRI1T_perc) / (segLIDP * 100) + paramset.bias_thres_aLRI1T_add;
        t.segthres_aLI2T = int64mul(p.segprep_a_LI, paramset.bias_thres_aLRI2T_perc) / (segLIDP * 100) + paramset.bias_thres_aLRI2T_add; // higher > stronger bias
        t.segthres_aLI1t = int64mul(p.segprep_a_LI,          bias_thres_aLRI1t_perc) / (segLIDP * 100);
        t.segthres_aLI2t = int64mul(p.segprep_a_LI, paramset.bias_thres_aLRI2t_perc) / (segLIDP * 100); // lower > stronger bias

        t.segthres_aRI1T = int64mul(p.segprep_a_RI,          bias_thres_aLRI1T_perc) / (segRIDP * 100) + paramset.bias_thres_aLRI1T_add;
        t.segthres_aRI2T = int64mul(p.segprep_a_RI, paramset.bias_thres_aLRI2T_perc) / (segRIDP * 100) + paramset.bias_thres_aLRI2T_add; // higher > stronger bias
        t.segthres_aRI1t = int64mul(p.segprep_a_RI,          bias_thres_aLRI1t_perc) / (segRIDP * 100);
        t.segthres_aRI2t = int64mul(p.segprep_a_RI, paramset.bias_thres_aLRI2t_perc) / (segRIDP * 100); // lower > stronger bias
        
        const auto aLRP1t_avgmul_perc = (is_normal ? paramset.bias_thres_aLRP1Nt_avgmul_perc : paramset.bias_thres_aLRP1t_avgmul_perc);
        const auto aLRP2t_avgmul_perc = paramset.bias_thres_aLRP2t_avgmul_perc;
        const auto aLRB1t_avgmul_perc = (is_normal ? paramset.bias_thres_aLRB1Nt_avgmul_perc : paramset.bias_thres_aLRB1t_avgmul_perc);
        const auto aLRB2t_avgmul_perc = paramset.bias_thres_aLRB2t_avgmul_perc;
        
        t.segthres_aLP1t = non_neg_minus(int64mul(p.segprep_a_l_dist_sum, aLRP1t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRP1t_minus);
        t.segthres_aLP2t = non_neg_minus(int64mul(p.segprep_a_l_dist_sum, aLRP2t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRP2t_minus);
        t.segthres_aRP1t = non_neg_minus(int64mul(p.segprep_a_r_dist_sum, aLRP1t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRP1t_minus);
        t.segthres_aRP2t = non_neg_minus(int64mul(p.segprep_a_r_dist_sum, aLRP2t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRP2t_minus);
        
        const auto pdel = p.segprep_a_delBAQ_sum / MAX(1, p.segprep_a_highBQ_dp); 
        t.segthres_aLB1t = non_neg_minus(int64mul(p.segprep_a_l_BAQ_sum, aLRB1t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRB1t_minus + pdel);
        t.segthres_aLB2t = non_neg_minus(int64mul(p.segprep_a_l_BAQ_sum, aLRB2t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRB2t_minus);
        t.segthres_aRB1t = non_neg_minus(int64mul(p.segprep_a_r_BAQ_sum, aLRB1t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRB1t_minus + pdel);
        t.segthres_aRB2t = non_neg_minus(int64mul(p.segprep_a_r_BAQ_sum, aLRB2t_avgmul_perc) / MAX(1, p.segprep_a_highBQ_dp * 100), paramset.bias_thres_aLRB2t_minus);
    }
    return 0;
}

template <bool isGap>
inline
int
dealwith_segbias(
        const auto bq,
        const uvc1_refgpos_t rpos,
        auto & symbol_to_seg_format_depth_set,
        auto & symbol_to_VQ_format_tag_set,
        const auto & seg_format_thres_set,
        const bam1_t *aln,
        const auto xm1500,
        const auto go1500 IGNORE_UNUSED_PARAM,
        const auto bm1500,
        const CoveredRegion<uvc1_qual_big_t> & baq_offsetarr,
        const CoveredRegion<uvc1_qual_big_t> & baq_offsetarr2,
        
        const auto cigar_op,
        const auto indel_len_arg,
        const uvc1_refgpos_t dist_to_interfering_indel,
        const uvc1_flag_t dflag,

        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    const bool is_assay_amplicon = (dflag & 0x4);
    const bool is_normal_used_to_filter_vars_on_primers = (paramset.tn_is_paired && (0x1 & paramset.tn_flag));
    const bool is_assay_UMI = (dflag & 0x1);
    
    const auto indel_len = UNSIGN2SIGN(indel_len_arg);
    const uvc1_refgpos_t bias_thres_veryhighBQ = paramset.bias_thres_highBQ; // veryhigh overriden by high. TODO: check if it makes sense?
    
    const auto rend = bam_endpos(aln);
    
    const uvc1_qual_t seg_l_baq = baq_offsetarr.getByPos(rpos) - baq_offsetarr.getByPos(aln->core.pos) + 1;
    const uvc1_qual_t _seg_r_baq = baq_offsetarr.getByPos(rend-1) - baq_offsetarr.getByPos(rpos) + 1;
    const uvc1_qual_t seg_r_baq = (isGap ? MIN(_seg_r_baq, baq_offsetarr2.getByPos(rend-1) - baq_offsetarr2.getByPos(rpos) + 7) : _seg_r_baq);
    
    const uvc1_refgpos_t seg_l_nbases = (rpos - aln->core.pos + 1);
    const uvc1_refgpos_t seg_r_nbases = (bam_endpos(aln) - rpos);
    const uvc1_refgpos_t frag_pos_L = MIN(aln->core.pos, aln->core.mpos);
    const uvc1_refgpos_t frag_pos_R = frag_pos_L + abs(aln->core.isize);
    const uvc1_refgpos_t frag_l_nbases2 = ((aln->core.isize != 0) ? MIN(rpos - frag_pos_L + 1, MAX_INSERT_SIZE) : (MAX_INSERT_SIZE));
    const uvc1_refgpos_t frag_r_nbases2 = ((aln->core.isize != 0) ? MIN(frag_pos_R - rpos + 0, MAX_INSERT_SIZE) : (MAX_INSERT_SIZE));
    
    const bool is_normal = ((aln->core.isize != 0) || (0 == (aln->core.flag & 0x1)));
    const bool isrc = ((aln->core.flag & 0x10) == 0x10);
    const bool isr2 = ((aln->core.flag & 0x80) == 0x80 && (aln->core.flag & 0x1) == 0x1);
    const bool strand = (isrc ^ isr2);
    
    const auto a1BQ = (isrc ? VQ_a1BQr : VQ_a1BQf);
    const auto a2BQ = (isrc ? VQ_a2BQr : VQ_a2BQf);
    symbol_to_VQ_format_tag_set[a1BQ] += bq;
    symbol_to_VQ_format_tag_set[a2BQ] += bq * bq / SQR_QUAL_DIV;
    
    symbol_to_seg_format_depth_set.seginfo_aMQs += aln->core.qual;
    auto & symbol_to_seg_aDP_depth_set = (strand
        ? (isrc ? symbol_to_seg_format_depth_set.seginfo_aDPrr : symbol_to_seg_format_depth_set.seginfo_aDPrf) 
        : (isrc ? symbol_to_seg_format_depth_set.seginfo_aDPfr : symbol_to_seg_format_depth_set.seginfo_aDPff));
    
    symbol_to_seg_aDP_depth_set += 1;
    if (MIN3(dist_to_interfering_indel, seg_l_nbases, seg_r_nbases) >= paramset.bias_thres_interfering_indel) {
        symbol_to_seg_format_depth_set.seginfo_aP3 += 1;
    }
    
    symbol_to_seg_format_depth_set.seginfo_aLPT += seg_l_nbases;
    symbol_to_seg_format_depth_set.seginfo_aRPT += seg_r_nbases;
    if (isrc) {
        symbol_to_seg_format_depth_set.seginfo_aLIT += ((aln->core.isize != 0) ? frag_l_nbases2 : 0);
    } else {
        symbol_to_seg_format_depth_set.seginfo_aRIT += ((aln->core.isize != 0) ? frag_r_nbases2 : 0);
    }
    
#if COMPILATION_ENABLE_XMGOT
    const auto const_XM1T = seg_format_thres_set.segthres_aXM1T;
    const auto const_XM2T = seg_format_thres_set.segthres_aXM2T;
    const auto const_GO1T = seg_format_thres_set.segthres_aGO1T;
    const auto const_GO2T = seg_format_thres_set.segthres_aGO2T;
#endif
    
    auto _const_LPxT = seg_format_thres_set.segthres_aLPxT; 
    const auto const_RPxT = seg_format_thres_set.segthres_aRPxT; 
    const auto const_LPxT = (isGap ? _const_LPxT : MIN(_const_LPxT, const_RPxT));
    
    const bool is_far_from_edge = (seg_l_nbases + ((BAM_CINS == cigar_op) ? non_neg_minus(indel_len, paramset.microadjust_nobias_pos_indel_maxlen) : 0) >= const_LPxT) && (seg_r_nbases >= const_RPxT);
    const auto bias_thres_highBAQ = paramset.bias_thres_highBAQ + (isGap ? 0 : 3);
    const bool is_unaffected_by_edge = (seg_l_baq >= bias_thres_highBAQ && seg_r_baq >= bias_thres_highBAQ);
    
    const auto min_dist2iend = ((aln->core.flag & 0x1) ? MIN(frag_l_nbases2, frag_r_nbases2) : (isrc ? seg_r_nbases : seg_l_nbases));
    if (is_far_from_edge && is_unaffected_by_edge && (min_dist2iend > paramset.primerlen2 || !is_assay_amplicon)) {
        symbol_to_seg_format_depth_set.seginfo_aP1 += 1;
    }
    if (is_assay_UMI || !is_assay_amplicon) {
        symbol_to_seg_format_depth_set.seginfo_aP2 += 1;
    }
    
    uvc1_readnum100x_t ampfact1, ampfact2;
    
    if (isGap) {
        ampfact1 = 100;
#if COMPILATION_ENABLE_XMGOT
        if (xm1500 > const_XM1T || go1500 > const_GO1T) {
            ampfact1 = MIN(
                    100 * mathsquare(const_XM1T) / MAX(1, mathsquare(xm1500)), 
                    100 * mathsquare(const_GO1T) / MAX(1, mathsquare(go1500)));
        }
#endif
        ampfact2 = 100;
        if (bq < paramset.bias_thres_PFBQ1) {
            ampfact2 = 100 * mathsquare(bq) / mathsquare(paramset.bias_thres_PFBQ1);
        }
        
        symbol_to_seg_format_depth_set.seginfo_aPF1 += MIN(ampfact1, ampfact2);
        
        ampfact1 = 100;
#if COMPILATION_ENABLE_XMGOT
        if ((xm1500 > const_XM2T || go1500 > const_GO2T)) {
            ampfact1 = MIN(
                    100 * mathsquare(const_XM2T) / MAX(1, mathsquare(xm1500)),
                    100 * mathsquare(const_GO2T) / MAX(1, mathsquare(go1500)));
        }
#endif
        ampfact2 = 100;
        if (bq < paramset.bias_thres_PFBQ2) {
            ampfact2 = 100 * mathsquare(bq) / mathsquare(paramset.bias_thres_PFBQ2);
        }
        
        symbol_to_seg_format_depth_set.seginfo_aPF2 += MIN(ampfact1, ampfact2);
        
    } else {
        ampfact1 = 100;
#if COMPILATION_ENABLE_XMGOT
        if (xm1500 > const_XM1T) {
            ampfact1 = 100 * mathsquare(const_XM1T) / mathsquare(xm1500);
        }
#endif
        ampfact2 = 100;
        if (bq < paramset.bias_thres_PFBQ1) {
            ampfact2 = 100 * mathsquare(bq) / mathsquare(paramset.bias_thres_PFBQ1);
        }

        symbol_to_seg_format_depth_set.seginfo_aPF1 += (ampfact1 * ampfact2 / (100));
        
        ampfact1 = 100;
#if COMPILATION_ENABLE_XMGOT
        if (xm1500 > const_XM2T) {
            ampfact1 = 100 * mathsquare(const_XM2T) / mathsquare(xm1500);
        }
#endif
        ampfact2 = 100; 
        if (bq < paramset.bias_thres_PFBQ2) {
            ampfact2 = 100 * mathsquare(bq) / mathsquare(paramset.bias_thres_PFBQ2);
        }
        
        symbol_to_seg_format_depth_set.seginfo_aPF2 += (ampfact1 * ampfact2 / (100));
        
        symbol_to_seg_format_depth_set.seginfo_a2XM2 += (xm1500 > 20 ? (100 * mathsquare(20) / mathsquare(xm1500)) : 100);
        symbol_to_seg_format_depth_set.seginfo_a2BM2 += (bm1500 > 20 ? (100 * mathsquare(20) / mathsquare(bm1500)) : 100);
    }
    
    if (((!isGap) && bq >= paramset.bias_thres_highBQ) || (isGap && dist_to_interfering_indel >= paramset.bias_thres_interfering_indel)) {
        const bool is_l1_unbiased = (seg_l_nbases + indel_len >= seg_format_thres_set.segthres_aLP1t);
        const bool is_l2_unbiased = (seg_l_nbases + indel_len >= seg_format_thres_set.segthres_aLP2t);
        const bool is_r1_unbiased = (seg_r_nbases >= seg_format_thres_set.segthres_aRP1t);
        const bool is_r2_unbiased = (seg_r_nbases >= seg_format_thres_set.segthres_aRP2t);
        if (is_far_from_edge) {
            if (is_l1_unbiased) {
                symbol_to_seg_format_depth_set.seginfo_aLP1 += 1;
            }
            if (is_l2_unbiased && (isGap || bq >= bias_thres_veryhighBQ)) {
                symbol_to_seg_format_depth_set.seginfo_aLP2 += 1;
            }
            if (is_r1_unbiased) {
                symbol_to_seg_format_depth_set.seginfo_aRP1 += 1;
            }
            if (is_r2_unbiased && (isGap || bq >= bias_thres_veryhighBQ)) {
                symbol_to_seg_format_depth_set.seginfo_aRP2 += 1;
            }
            symbol_to_seg_format_depth_set.seginfo_aLPL += seg_l_nbases;
            symbol_to_seg_format_depth_set.seginfo_aRPL += seg_r_nbases;
        }
        if (is_unaffected_by_edge) {
            if (seg_l_baq >= paramset.bias_thres_BAQ1) {
                symbol_to_seg_format_depth_set.seginfo_aLB1 += 1;
            }
            if (seg_l_baq >= paramset.bias_thres_BAQ2 && (isGap || bq >= bias_thres_veryhighBQ)) {
                symbol_to_seg_format_depth_set.seginfo_aLB2 += 1;
            }
            if (seg_r_baq >= paramset.bias_thres_BAQ1) {
                symbol_to_seg_format_depth_set.seginfo_aRB1 += 1;
            }
            if (seg_r_baq >= paramset.bias_thres_BAQ2 && (isGap || bq >= bias_thres_veryhighBQ)) {
                symbol_to_seg_format_depth_set.seginfo_aRB2 += 1;
            }
            symbol_to_seg_format_depth_set.seginfo_aLBL += seg_l_baq;
            symbol_to_seg_format_depth_set.seginfo_aRBL += seg_r_baq;
        }
        symbol_to_seg_format_depth_set.seginfo_aBQ2 += 1;
    }
    
    const bool is_l_nonbiased = (((0 == (aln->core.flag & 0x8)) || (0 == (aln->core.flag & 0x1))) && seg_l_nbases > seg_r_nbases);
    const bool is_r_nonbiased = (((0 == (aln->core.flag & 0x8)) || (0 == (aln->core.flag & 0x1))) && seg_l_nbases < seg_r_nbases);
    // rc : test if the left-side of the insert is biased, and vice versa.
    const bool is_pos_good_for_bias_calc = ((!is_assay_amplicon) || (!is_normal_used_to_filter_vars_on_primers) || ((is_far_from_edge && is_unaffected_by_edge))); 
    if (isrc) {
        auto dist2iend = frag_l_nbases2;
        if ((dist2iend >= seg_format_thres_set.segthres_aLI1t) && (dist2iend <= seg_format_thres_set.segthres_aLI1T || isGap) 
                && (is_normal || (isGap && is_l_nonbiased))) {
            symbol_to_seg_format_depth_set.seginfo_aLI1 += 1; // aLI1
        } 
        if ((dist2iend >= seg_format_thres_set.segthres_aLI2t) && (dist2iend <= seg_format_thres_set.segthres_aLI2T || isGap) 
                && (is_normal || (isGap && is_l_nonbiased))) {
            if (is_pos_good_for_bias_calc) { symbol_to_seg_format_depth_set.seginfo_aLI2 += 1; } // aLI2
        }
        if (is_pos_good_for_bias_calc) { symbol_to_seg_format_depth_set.seginfo_aLIr += 1; }
    } else {
        auto dist2iend = frag_r_nbases2;
        if ((dist2iend >= seg_format_thres_set.segthres_aRI1t) && (dist2iend <= seg_format_thres_set.segthres_aRI1T || isGap) 
                && (is_normal || (isGap && is_r_nonbiased))) {
            symbol_to_seg_format_depth_set.seginfo_aRI1 += 1; // aRI1
        } 
        if ((dist2iend >= seg_format_thres_set.segthres_aRI2t) && (dist2iend <= seg_format_thres_set.segthres_aRI2T || isGap) 
                && (is_normal || (isGap && is_r_nonbiased))) {
            if (is_pos_good_for_bias_calc) { symbol_to_seg_format_depth_set.seginfo_aRI2 += 1; } // aRI2
        }
        if (is_pos_good_for_bias_calc) { symbol_to_seg_format_depth_set.seginfo_aRIf += 1; }
    }
    
    return 0;
}

template <class TSymbol2Count>
class GenericSymbol2CountCoverage : public CoveredRegion<TSymbol2Count> {
public:
    GenericSymbol2CountCoverage() : CoveredRegion<TSymbol2Count>(0, 0, 1) { };
    GenericSymbol2CountCoverage(auto tid, auto beg, auto end) : CoveredRegion<TSymbol2Count>(tid, beg, end) {}
    
    void
    assertUpdateIsLegal(const GenericSymbol2CountCoverage<TSymbol2Count> & other) const {
        assert(this->tid == other.tid);
        assert(this->getIncluBegPosition() <= other.getIncluBegPosition() || !fprintf(stderr, "%d <= %d failed!", this->getIncluBegPosition(), other.getIncluBegPosition()));
        assert(this->getExcluEndPosition() >= other.getExcluEndPosition() || !fprintf(stderr, "%d >= %d failed!", this->getExcluEndPosition(), other.getExcluEndPosition())); 
    }
    // mainly for grouping reads sharing the same UMI into one family
    template <bool TIsIncVariable = true>
    void
    updateByRepresentative(const GenericSymbol2CountCoverage<TSymbol2Count> & other, uvc1_readnum_t incvalue = 1,
            const bool update_pos2indel2count = true, const bool update_idx2symbol2data = true) {
        this->assertUpdateIsLegal(other);
        
        if (update_idx2symbol2data) {
            for (uvc1_refgpos_t epos = UNSIGN2SIGN(other.getIncluBegPosition()); epos < UNSIGN2SIGN(other.getExcluEndPosition()); epos++) {
                AlignmentSymbol consymbol = this->getRefByPos(epos).template updateByRepresentative<TIsIncVariable>(other.getByPos(epos), incvalue);
                if (update_pos2indel2count) {
                    if (isSymbolIns(consymbol)) {
                        posToIndelToCount_updateByRepresentative<TIsIncVariable>(this->getRefPosToIseqToData(consymbol), other.getPosToIseqToData(consymbol), epos, incvalue);
                    } else if (isSymbolDel(consymbol)) {
                        posToIndelToCount_updateByRepresentative<TIsIncVariable>(this->getRefPosToDlenToData(consymbol), other.getPosToDlenToData(consymbol), epos, incvalue);
                    }
                }
            }
        }
    }

    // mainly for merging R1 and R2 into one read
    template<bool TIndelIsMajor = false> 
    void
    updateByConsensus(const GenericSymbol2CountCoverage<TSymbol2Count> &other, uvc1_readnum_t incvalue = 1,
            const bool update_pos2indel2count = true, const bool update_idx2symbol2data = true) {
        this->assertUpdateIsLegal(other);
        
        if (update_idx2symbol2data) {
            for (uvc1_refgpos_t epos = (other.getIncluBegPosition()); epos < UNSIGN2SIGN(other.getExcluEndPosition()); epos++) {
                const std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols = this->getRefByPos(epos).template updateByConsensus<TIndelIsMajor>(other.getByPos(epos), incvalue);
                if (update_pos2indel2count) {
                    if (isSymbolIns(consymbols[1])) {
                        posToIndelToCount_updateByConsensus(this->getRefPosToIseqToData(consymbols[1]), other.getPosToIseqToData(consymbols[1]), epos, incvalue);
                    } else if (isSymbolDel(consymbols[1])) {
                        posToIndelToCount_updateByConsensus(this->getRefPosToDlenToData(consymbols[1]), other.getPosToDlenToData(consymbols[1]), epos, incvalue);
                    }
                }
            }
        }
    }
    
    // Add read supports to a bigger family, while excluding read supports that did not pass the threshold.
    int
    updateByFiltering(
            const GenericSymbol2CountCoverage<TSymbol2Count> &other, 
            const std::array<uvc1_qual_t, NUM_SYMBOL_TYPES> thres,
            uvc1_readnum_t incvalue = 1,
            const bool update_pos2indel2count = true) {
        this->assertUpdateIsLegal(other);
        int num_updated_pos = 0;
        // other.assertUpdateIsLegal(thres); // may not hold because threshold is delimited by bed whereas this is not delimited by bed.
        auto incluBegPos = other.getIncluBegPosition();
        auto excluEndPos = other.getExcluEndPosition();
        std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols; 
        for (auto epos = incluBegPos; epos < excluEndPos; epos++) {
            int updateresult = this->getRefByPos(epos).updateByFiltering(
                    consymbols, 
                    other.getByPos(epos), 
                    thres,
                    incvalue);
            if (update_pos2indel2count) {
                if (isSymbolIns(consymbols[1])) {
                    posToIndelToCount_updateByConsensus(this->getRefPosToIseqToData(consymbols[1]), other.getPosToIseqToData(consymbols[1]), epos, incvalue);
                } else if (isSymbolDel(consymbols[1])) {
                    posToIndelToCount_updateByConsensus(this->getRefPosToDlenToData(consymbols[1]), other.getPosToDlenToData(consymbols[1]), epos, incvalue);
                }
            }
            if (updateresult) { num_updated_pos++; }
        }
        return num_updated_pos;
    };
    
    int
    updateByMajorMinusMinor(
            const GenericSymbol2CountCoverage<TSymbol2Count> &other,
            const bool update_pos2indel2count = true) {
        this->assertUpdateIsLegal(other);
        int num_updated_pos = 0;
        // other.assertUpdateIsLegal(thres); // may not hold because threshold is delimited by bed whereas this is not delimited by bed.
        std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols; 
        std::array<uvc1_readnum_t, NUM_SYMBOL_TYPES> conmmmcnts;
        for (auto epos = other.getIncluBegPosition(); epos < other.getExcluEndPosition(); epos++) {
            int updateresult = this->getRefByPos(epos).updateByMajorMinusMinor(
                    consymbols,
                    conmmmcnts,
                    other.getByPos(epos));
            if (update_pos2indel2count) {
                if (isSymbolIns(consymbols[LINK_SYMBOL])) {
                    posToIndelToCount_updateByConsensus(this->getRefPosToIseqToData(consymbols[1]), other.getPosToIseqToData(consymbols[1]), 
                            epos, conmmmcnts[LINK_SYMBOL]);
                } else if (isSymbolDel(consymbols[LINK_SYMBOL])) {
                    posToIndelToCount_updateByConsensus(this->getRefPosToDlenToData(consymbols[1]), other.getPosToDlenToData(consymbols[1]), 
                            epos, conmmmcnts[LINK_SYMBOL]);
                }
            }
            if (updateresult) { num_updated_pos++; }
        }
        return num_updated_pos;
    };

    template <ValueType TUpdateType2>
    void // GenericSymbol2CountCoverage<TSymbol2Count>::
    inc(const uvc1_refgpos_t epos, const AlignmentSymbol symbol, const uvc1_readnum_t incvalue = 1, const bam1_t *bam = NULL) {
        auto &r = this->getRefByPos(epos, bam);
        r.template incSymbolCount<(TUpdateType2)>(symbol, incvalue);
    };

    void // GenericSymbol2CountCoverage<TSymbol2Count>::
    incIns(const uvc1_refgpos_t epos, const std::string & iseq, const AlignmentSymbol symbol, const uvc1_readnum_t incvalue = 1) {
        assert (incvalue > 0); 
        assert (iseq.size() > 0);
        size_t ipos = epos;
        posToIndelToCount_inc(this->getRefPosToIseqToData(symbol), ipos, iseq, incvalue);
    };

    void // GenericSymbol2CountCoverage<TSymbol2Count>::
    incDel(const uvc1_refgpos_t epos, const uvc1_refgpos_t dlen, const AlignmentSymbol symbol, const uvc1_readnum_t incvalue = 1) {
        assert (incvalue > 0);
        assert (dlen > 0);
        size_t ipos = epos;
        posToIndelToCount_inc(this->getRefPosToDlenToData(symbol), ipos, dlen, incvalue);
    };
    
    template<bool TIsProton, ValueType TUpdateType, bool TIsBiasUpdated>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByAln(
            const bam1_t *const aln, 
            
            const uvc1_refgpos_t region_offset,
            const auto & region_symbolvec, 
            const std::vector<RegionalTandemRepeat> & region_repeatvec,
            const CoveredRegion<uvc1_qual_big_t> & baq_offsetarr,
            const CoveredRegion<uvc1_qual_big_t> & baq_offsetarr2,

            auto & seg_format_depth_sets,
            auto & symbol_to_VQ_format_tag_sets,
            const auto & seg_format_prep_sets,
            const auto & seg_format_thres_sets,
            
            const uvc1_flag_t dflag,
            const CommandLineArgs & paramset,
            const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
        
        static_assert(BASE_QUALITY_MAX == TUpdateType || SYMBOL_COUNT_SUM == TUpdateType);
        assert(this->tid == SIGN2UNSIGN(aln->core.tid));
        assert(this->getIncluBegPosition() <= SIGN2UNSIGN(aln->core.pos)   || !fprintf(stderr, "%d <= %d failed", this->getIncluBegPosition(), aln->core.pos));
        assert(this->getExcluEndPosition() >= SIGN2UNSIGN(bam_endpos(aln)) || !fprintf(stderr, "%d >= %d failed", this->getExcluEndPosition(), bam_endpos(aln)));
        
        const bool is_assay_amplicon = (dflag & 0x4);

        const auto n_cigar = aln->core.n_cigar;
        const auto *cigar = bam_get_cigar(aln);
        const auto *bseq = bam_get_seq(aln);
        const auto rend = bam_endpos(aln);
        const std::array<uvc1_qual_t, NUM_SYMBOL_TYPES> symboltype2addPhred = {{paramset.bq_phred_added_misma, paramset.bq_phred_added_indel}};
        
        uvc1_readnum_t nge_cnt = 0; // number of gap extensions.
        uvc1_readnum_t ngo_cnt = 0; // number of gap opens.
        for (uvc1_unsigned_int_t i = 0; i < n_cigar; i++) {
            const auto c = cigar[i];
            const auto cigar_op = bam_cigar_op(c);
            if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) { 
                nge_cnt += bam_cigar_oplen(c); 
                ngo_cnt++; 
            }
        }

        const auto *bam_aux_data = bam_aux_get(aln, "NM");
        const uvc1_readnum_t nm_cnt = ((bam_aux_data != NULL) ? bam_aux2i(bam_aux_data) : nge_cnt);
        assert (nm_cnt >= nge_cnt);
        const uvc1_base1500x_t xm_cnt = nm_cnt - nge_cnt;
        const uvc1_base1500x_t xm1500 = xm_cnt * 1500 / (rend - aln->core.pos);
        const uvc1_base1500x_t go1500 = ngo_cnt * 1500 / (rend - aln->core.pos);
        
        std::vector<uvc1_refgpos_t> indel_rposs = {{ 0 }};
        size_t indel_rposs_idx = 0;
        std::array<uvc1_readpos_t, NUM_ALIGNMENT_SYMBOLS> bm_cnts = {{ 0 }}; // mismatch of the same base type
        {
            uvc1_refgpos_t qpos = 0;
            uvc1_refgpos_t rpos = aln->core.pos;
            for (uint32_t i = 0; i < n_cigar; i++) {
                const auto c = cigar[i];
                const auto cigar_op = bam_cigar_op(c);
                const auto cigar_oplen = bam_cigar_oplen(c);
                if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                    for (uint32_t i2 = 0; i2 < cigar_oplen; i2++) {
                        const auto base4bit = bam_seqi(bseq, qpos);
                        const auto base3bit = seq_nt16_int[base4bit];
                        const auto refbase = region_symbolvec[rpos - region_offset];
                        if (refbase != base3bit) { 
                            AlignmentSymbol symbol = AlignmentSymbol(base3bit);
                            bm_cnts[symbol] += 1;
                        }
                        qpos++;
                        rpos++;
                    }
                } else if (cigar_op == BAM_CINS) {
                    bool is_lowBQ = false;
                    for (uvc1_refgpos_t qpos2 = qpos - MIN(qpos, 1); qpos2 < MIN(qpos + UNSIGN2SIGN(cigar_oplen) + 1, rend); qpos2++) {
                        if (BAM_PHREDI(aln, qpos2) < paramset.bias_thres_interfering_indel_BQ) { is_lowBQ = true; }
                    }
                    if (is_lowBQ) { 
                        indel_rposs.push_back(rpos); 
                    }
                    qpos += cigar_oplen;
                } else if (cigar_op == BAM_CDEL) {
                    bool is_lowBQ = (MIN(BAM_PHREDI(aln, MAX(1, qpos) - 1), BAM_PHREDI(aln, qpos)) <= paramset.bias_thres_interfering_indel_BQ);
                    if (is_lowBQ) { 
                        indel_rposs.push_back(rpos); 
                    }
                    rpos += cigar_oplen;
                } else {
                    process_cigar(qpos, rpos, cigar_op, cigar_oplen);
                }
            }
            indel_rposs.push_back(INT32_MAX); 
        }
        std::array<uvc1_base1500x_t, NUM_ALIGNMENT_SYMBOLS> bm1500s = {{ 0 }};
        for (size_t i = 0; i < bm_cnts.size(); i++) {
            bm1500s[i] = bm_cnts[i] * 1500 / (rend - aln->core.pos);
        }
        
        const bool isrc = ((aln->core.flag & 0x10) == 0x10);
        
        const bool is_normal_used_to_filter_vars_on_primers = (paramset.tn_is_paired && (0x1 & paramset.tn_flag));
        const uvc1_readpos_t primerlen_lside = paramset.primerlen;
        const uvc1_readpos_t primerlen_rside = paramset.primerlen;
        // NOTE: if insert size is zero, then assume that the whole insert including the two primers are sequenced
        // A correct pair of primers should not generate zero insert size
        const uvc1_readpos_t ibeg = ((aln->core.isize != 0) ? (MIN(aln->core.pos, aln->core.mpos) + primerlen_lside) 
                : ((isrc && (0x0 == (0x1 & aln->core.flag))) ? 0 : (aln->core.pos + primerlen_lside)));
        const uvc1_readpos_t iend = ((aln->core.isize != 0) ? non_neg_minus(MIN(aln->core.pos, aln->core.mpos) + abs(aln->core.isize), primerlen_rside) 
                : ((isrc && (0x0 == (0x1 & aln->core.flag))) ? non_neg_minus(rend, primerlen_rside) : INT32_MAX));
        
        uvc1_readpos_t qpos = 0;
        uvc1_readpos_t rpos = aln->core.pos;
        uvc1_qual_t incvalue = 1; // can also be uvc1_readnum_t or uvc1_readpos_t
        uvc1_readpos_t lclip_len = ((n_cigar > 0 && bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) ? bam_cigar_oplen(cigar[0]) : 0);
        uvc1_readpos_t rclip_len = ((n_cigar > 0 && bam_cigar_op(cigar[n_cigar - 1]) == BAM_CSOFT_CLIP) ? bam_cigar_oplen(cigar[n_cigar - 1]) : 0);
        const uvc1_qual_t micro_indel_penal_by_clip = MAX(lclip_len, rclip_len) / 6; // / paramset.indel_str_repeatsize_max;
        const uvc1_qual_t micro_indel_penal_by_nm = (xm1500 + go1500) / (30); // paramset.bias_thres_aXM1T_add;
        const uvc1_qual_t micro_indel_penal = MIN(1, micro_indel_penal_by_nm + micro_indel_penal_by_clip);
        const uvc1_qual_t micro_nogap_penal = MIN(4, micro_indel_penal_by_nm + micro_indel_penal_by_clip) + 1;
        for (uint32_t i = 0; i < n_cigar; i++) {
            const auto c = cigar[i];
            const auto cigar_op = bam_cigar_op(c);
            const auto cigar_oplen = bam_cigar_oplen(c);
            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                for (uint32_t i2 = 0; i2 < cigar_oplen; i2++) {
                    assert((rpos >= SIGN2UNSIGN(aln->core.pos) && rpos < SIGN2UNSIGN(bam_endpos(aln)))
                            || !fprintf(stderr, "Bam line with QNAME %s has rpos %d which is not the in range (%d - %d)", 
                            bam_get_qname(aln), rpos, aln->core.pos, bam_endpos(aln)));
if ((is_normal_used_to_filter_vars_on_primers || !is_assay_amplicon) || (ibeg <= rpos && rpos < iend)) {

                    uvc1_refgpos_t dist_to_interfering_indel = 10000;
                    if (TIsBiasUpdated && nge_cnt > 0) {
                        if (indel_rposs[indel_rposs_idx] <= rpos) {
                            indel_rposs_idx++;
                        }
                        uvc1_refgpos_t prev_indel_rpos = indel_rposs[indel_rposs_idx - 1];
                        uvc1_refgpos_t next_indel_rpos = indel_rposs[indel_rposs_idx];
                        
                        const auto & rtr1 = region_repeatvec[MAX(rpos-region_offset, paramset.indel_adj_tracklen_dist) - paramset.indel_adj_tracklen_dist];
                        const auto & rtr2 = region_repeatvec[MIN(rpos-region_offset + paramset.indel_adj_tracklen_dist, UNSIGN2SIGN(region_repeatvec.size()) - 1)];
                        const auto prevlen = non_neg_minus(rpos - prev_indel_rpos, MAX(
                                rpos - (region_offset + rtr1.begpos),
                                seg_format_thres_sets.getByPos(rpos).segthres_aLP1t));
                        const auto nextlen = non_neg_minus(next_indel_rpos - rpos, MAX(
                                (region_offset + rtr2.begpos + rtr2.tracklen) - rpos,
                                seg_format_thres_sets.getByPos(rpos).segthres_aRP1t));
                        dist_to_interfering_indel = MIN(prevlen, nextlen);
                        assert (rpos >= (region_offset + rtr1.begpos));
                        assert ((region_offset + rtr2.begpos + rtr2.tracklen) >= rpos);
                    }
                    
                    if (i2 > 0) {
                        const auto noindel_phredvalue = MIN(
                                region_repeatvec[rpos-region_offset - 1].indelphred,
                                region_repeatvec[rpos-region_offset    ].indelphred);
                       const uvc1_qual_t qfromBQ2 = (TIsProton ? MIN(BAM_PHREDI(aln, qpos-1), BAM_PHREDI(aln, qpos)) : 80);
                       incvalue = non_neg_minus(MIN(qfromBQ2, noindel_phredvalue), micro_nogap_penal) + 1;
                        this->template inc<TUpdateType>(rpos, LINK_M, incvalue, aln);
                        if (TIsBiasUpdated) {
                            dealwith_segbias<true>(
                                    incvalue,
                                    rpos,
                                    seg_format_depth_sets.getRefByPos(rpos)[LINK_M],
                                    symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[LINK_M],
                                    seg_format_thres_sets.getByPos(rpos),
                                    aln,
                                    xm1500,
                                    go1500,
                                    bm1500s[LINK_M],
                                    baq_offsetarr,
                                    baq_offsetarr2,
                                    
                                    cigar_op,
                                    0,
                                    dist_to_interfering_indel,
                                    dflag,
                                    
                                    paramset,
                                    0);
                        }
                    }
                    const auto base4bit = bam_seqi(bseq, qpos);
                    const auto base3bit = seq_nt16_int[base4bit];
                    AlignmentSymbol symbol = AlignmentSymbol(base3bit);
                    if (TIsProton && ((0 == i2) || (cigar_oplen - 1 == i2))) {
                        const auto prev_cigar = (0 < i ? cigar[i - 1] : -1);
                        const auto next_cigar = (i + 1 < n_cigar ? cigar[i + 1] : -1);
                        if (((cigar_oplen - 1 == i2) && (BAM_CMATCH != next_cigar) && (BAM_CEQUAL != next_cigar) && (BAM_CDIFF != next_cigar))
                                       || ((0 == i2) && (BAM_CMATCH != prev_cigar) && (BAM_CEQUAL != prev_cigar) && (BAM_CDIFF != prev_cigar))) {
                            
                            const bool isrc2 = (0 != i2);
                            uvc1_qual_t prev_base_phred = 1;
                            if ((isrc2) && (qpos + 1 < aln->core.l_qseq)) {
                                prev_base_phred = BAM_PHREDI(aln, qpos + 1);
                            }
                            if ((!isrc2) && (qpos > 0)) {
                                prev_base_phred = BAM_PHREDI(aln, qpos - 1);
                            }
                            
                            incvalue = MIN(BAM_PHREDI(aln, qpos), prev_base_phred) + MIN(symboltype2addPhred[BASE_SYMBOL], symboltype2addPhred[LINK_SYMBOL]);
                        } else {
                            incvalue = BAM_PHREDI(aln, qpos) + symboltype2addPhred[BASE_SYMBOL]; 
                        }
                    } else {
                        incvalue = BAM_PHREDI(aln, qpos) + symboltype2addPhred[BASE_SYMBOL];
                    }
                    this->template inc<TUpdateType>(rpos, AlignmentSymbol(base3bit), incvalue, aln);
                    if (TIsBiasUpdated) {
                        dealwith_segbias<false>(
                                incvalue,
                                rpos,
                                seg_format_depth_sets.getRefByPos(rpos)[symbol],
                                symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                                seg_format_thres_sets.getByPos(rpos),
                                aln,
                                xm1500,
                                go1500,
                                bm1500s[symbol],
                                baq_offsetarr,
                                baq_offsetarr2,
                                
                                cigar_op,
                                0,
                                dist_to_interfering_indel,
                                dflag,

                                paramset,
                                0);
                    }
                    
}
                    rpos += 1;
                    qpos += 1;
                }
            } else if (cigar_op == BAM_CINS) {
if ((is_normal_used_to_filter_vars_on_primers || !is_assay_amplicon) || (ibeg <= rpos && rpos < iend)) {
                const auto nbases2end = MIN(qpos, UNSIGN2SIGN(aln->core.l_qseq) - UNSIGN2SIGN(qpos + UNSIGN2SIGN(cigar_oplen)));
                const bool is_ins_at_read_end = (nbases2end <= 0);
                uvc1_readpos_t inslen = cigar_oplen;
                if (is_ins_at_read_end) {
                    if (UNSIGN2SIGN(cigar_oplen) > paramset.debug_warn_min_read_end_ins_cigar_oplen) {
                        LOG(logWARNING) << "Query " << bam_get_qname(aln) << " has insertion of legnth " << cigar_oplen << " at " << qpos
                                << " which is not exclusively between 0 and " << aln->core.l_qseq << " aligned to tid " << aln->core.tid << " and position " << rpos;
                    }
                    incvalue = (0 != qpos ? BAM_PHREDI(aln, qpos-1) : 
                            ((qpos + UNSIGN2SIGN(cigar_oplen) < SIGN2UNSIGN(aln->core.l_qseq)) ? 
                            BAM_PHREDI(aln, qpos + SIGN2UNSIGN(cigar_oplen)) : 1)) 
                            + (symboltype2addPhred[LINK_SYMBOL]);
                } else {
                    uvc1_qual_t phredvalue = ref_to_phredvalue(
                            inslen, 
                            region_symbolvec, 
                            rpos - region_offset,
                            paramset.indel_BQ_max,
                            paramset.indel_polymerase_slip_rate,
                            cigar_oplen, 
                            cigar_op, 
                            paramset.indel_str_repeatsize_max,
                            paramset.indel_del_to_ins_err_ratio,
                            0);
                    uvc1_qual_t phredinc = (uvc1_qual_t)round(2 * numstates2phred((double)seg_format_prep_sets.getByPos(rpos).segprep_a_dp
                            / (double)(1.0 + non_neg_minus(seg_format_prep_sets.getByPos(rpos).segprep_a_dp, 
                                    seg_format_prep_sets.getByPos(rpos).segprep_a_at_ins_dp + seg_format_prep_sets.getByPos(rpos).segprep_a_at_del_dp))));
                    const auto & p = seg_format_prep_sets.getByPos(rpos);
                    const uvc1_readnum_t qfromBQ2_ratiothres = (NOT_PROVIDED == paramset.vcf_tumor_fname ? 2 : 4);
                    const bool is_multiallelic_ins = (p.segprep_a_near_ins_pow2len * qfromBQ2_ratiothres > MAX(1, p.segprep_a_near_ins_dp) * UNSIGN2SIGN(cigar_oplen * 3));
                    if (1 == inslen && !is_multiallelic_ins) { phredvalue += BETWEEN(phredinc - 3, 0, 4); }
                    const uvc1_readnum_t thisdp = (seg_format_prep_sets.getByPos(rpos).segprep_a_at_ins_dp); 
                    const uvc1_readnum_t neardp = (MAX(seg_format_prep_sets.getByPos(rpos).segprep_a_near_ins_dp, seg_format_prep_sets.getByPos(rpos).segprep_a_near_RTR_ins_dp));
                    uvc1_qual_t insbase_minphred = 80;
                    for (uvc1_refgpos_t qpos2 = qpos; qpos2 < qpos + UNSIGN2SIGN(cigar_oplen); qpos2++) {
                        UPDATE_MIN(insbase_minphred, BAM_PHREDI(aln, qpos2));
                    }
                    uvc1_qual_t ancbase_minphred = 80;
                    if (qpos > 0) {
                        ancbase_minphred = MIN(ancbase_minphred, BAM_PHREDI(aln, qpos - 1)); 
                    }
                    if (qpos + UNSIGN2SIGN(cigar_oplen) + 1 < aln->core.l_qseq) { 
                        ancbase_minphred = MIN(ancbase_minphred, BAM_PHREDI(aln, qpos + cigar_oplen + 1)); 
                    }
                    // IonTorrent may generate erroneous indels along with true indels, so be more lenient for IonTorrent sequencers.
                    uvc1_refgpos_t qfromBQ1 = (TIsProton ? ancbase_minphred : MIN(ancbase_minphred, insbase_minphred)); 
                    uvc1_refgpos_t qfromBQ2 = ((thisdp * qfromBQ2_ratiothres <= neardp || (1 == cigar_oplen && 
                                (xm1500 >= paramset.microadjust_xm || 
                                    ((lclip_len + paramset.microadjust_cliplen >= rpos - aln->core.pos) && isrc) 
                                 || ((rclip_len + paramset.microadjust_cliplen >= rend - aln->core.pos) && !isrc)))) 
                            ? qfromBQ1 : (TIsProton ? MIN(qfromBQ1 + proton_cigarlen2phred(cigar_oplen), MAX(3, qfromBQ1) * UNSIGN2SIGN(cigar_oplen)) : 80));
                    incvalue = non_neg_minus(MIN(qfromBQ2, phredvalue + symboltype2addPhred[LINK_SYMBOL]), micro_indel_penal) + 1;
                }
                if (nbases2end >= paramset.indel_filter_edge_dist) {
                    const auto symbol = insLenToSymbol(inslen, aln);
                    this->template inc<TUpdateType>(rpos, symbol, MAX(SIGN2UNSIGN(1), incvalue), aln);
                    if (TIsBiasUpdated) {
                        dealwith_segbias<true>(
                                MAX(SIGN2UNSIGN(1), incvalue),
                                rpos,
                                seg_format_depth_sets.getRefByPos(rpos)[symbol],
                                symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                                seg_format_thres_sets.getByPos(rpos),
                                aln,
                                xm1500,
                                go1500,
                                bm1500s[symbol],
                                baq_offsetarr,
                                baq_offsetarr2,
                                
                                cigar_op,
                                cigar_oplen,
                                10000,
                                dflag,
                                
                                paramset,
                                0);
                    }
                    std::string iseq;
                    iseq.reserve(cigar_oplen);
                    uvc1_qual_t incvalue2 = incvalue;
                    for (uint32_t i2 = 0; i2 < cigar_oplen; i2++) {
                        const auto base4bit = bam_seqi(bseq, qpos+i2);
                        const auto base8bit = seq_nt16_str[base4bit];
                        iseq.push_back(base8bit);
                        incvalue2 = MIN(incvalue2, SIGN2UNSIGN(bam_seqi(bseq, qpos+i2)))
                                + (symboltype2addPhred[LINK_SYMBOL]);
                    }
                    this->incIns(rpos, iseq, insLenToSymbol(inslen, aln), MAX(SIGN2UNSIGN(1), incvalue2));
                }
}
                qpos += cigar_oplen;
            } else if (cigar_op == BAM_CDEL) {
if ((is_normal_used_to_filter_vars_on_primers || !is_assay_amplicon) || (ibeg <= rpos && rpos < iend)) {
                const auto nbases2end = MIN(qpos, UNSIGN2SIGN(aln->core.l_qseq) - UNSIGN2SIGN(qpos));
                const bool is_del_at_read_end = (nbases2end <= 0);
                uvc1_readpos_t dellen = cigar_oplen;
                if (is_del_at_read_end) {
                    LOG(logWARNING) << "Query " << bam_get_qname(aln) << " has deletion of legnth " << cigar_oplen << " at " << qpos
                            << " which is not exclusively between 0 and " << aln->core.l_qseq << " aligned to tid " << aln->core.tid << " and position " << rpos; 
                    incvalue = (0 != qpos ? BAM_PHREDI(aln, qpos-1) : 
                            ((qpos < SIGN2UNSIGN(aln->core.l_qseq)) ? 
                            BAM_PHREDI(aln, qpos) : 1))
                            + (symboltype2addPhred[LINK_SYMBOL]);
                } else {
                    uvc1_qual_t phredvalue = ref_to_phredvalue(
                            dellen, 
                            region_symbolvec, 
                            rpos - region_offset, 
                            paramset.indel_BQ_max,
                            paramset.indel_polymerase_slip_rate, 
                            cigar_oplen, 
                            cigar_op, 
                            paramset.indel_str_repeatsize_max,
                            paramset.indel_del_to_ins_err_ratio,
                            0);
                    uvc1_qual_t phredinc = round(2 * numstates2phred((double)seg_format_prep_sets.getByPos(rpos).segprep_a_dp
                            / (double)(1.0 + non_neg_minus(seg_format_prep_sets.getByPos(rpos).segprep_a_dp, 
                                    seg_format_prep_sets.getByPos(rpos).segprep_a_at_ins_dp + seg_format_prep_sets.getByPos(rpos).segprep_a_at_del_dp))));
                    if (1 == dellen) { phredvalue += BETWEEN(phredinc - 3, 0, 4); }
                    uvc1_readnum_t thisdp = (seg_format_prep_sets.getByPos(rpos).segprep_a_at_del_dp);
                    uvc1_readnum_t neardp = (MAX(seg_format_prep_sets.getByPos(rpos).segprep_a_near_del_dp, seg_format_prep_sets.getByPos(rpos).segprep_a_near_RTR_del_dp));
                    uvc1_qual_t qfromBQ1 = MIN(BAM_PHREDI(aln, qpos), BAM_PHREDI(aln, qpos-1));
                    const uvc1_readnum_t qfromBQ2_ratiothres = (NOT_PROVIDED == paramset.vcf_tumor_fname ? 2 : 4);
                    uvc1_qual_t qfromBQ2 = ((thisdp * qfromBQ2_ratiothres <= neardp) ? non_neg_minus(qfromBQ1, 1) : (TIsProton ? MIN(qfromBQ1 + proton_cigarlen2phred(cigar_oplen), MAX(3, qfromBQ1) * UNSIGN2SIGN(cigar_oplen)) : 80));
                    double delFA = ((double)(thisdp + 0.5) / (double)(seg_format_prep_sets.getByPos(rpos).segprep_a_dp + 1));
                    uvc1_qual_t delFAQ = MAX(0, paramset.microadjust_delFAQmax + (uvc1_qual_t)round(paramset.powlaw_exponent * numstates2phred(delFA)));
                    
                    uint32_t prev_cidx = i;
                    uvc1_refgpos_t prev_rpos = rpos;
                    while ((0 != prev_cidx) && (BAM_CINS != bam_cigar_op(cigar[prev_cidx]) || cigar_oplen != bam_cigar_oplen(cigar[prev_cidx]))) { 
                        prev_cidx--;
                        if (BAM_CMATCH == bam_cigar_op(cigar[prev_cidx]) || BAM_CEQUAL == bam_cigar_op(cigar[prev_cidx]) || BAM_CDIFF == bam_cigar_op(cigar[prev_cidx]) 
                                || BAM_CDEL == bam_cigar_op(cigar[prev_cidx]) || BAM_CREF_SKIP == bam_cigar_op(cigar[prev_cidx])) {
                            prev_rpos -= bam_cigar_oplen(cigar[prev_cidx]);
                        }
                    }
                    uint32_t next_cidx = i;
                    uvc1_refgpos_t next_rpos = rpos + cigar_oplen;
                    while ((n_cigar - 1 != next_cidx) && (BAM_CINS != bam_cigar_op(cigar[next_cidx]) || cigar_oplen != bam_cigar_oplen(cigar[next_cidx]))) { 
                        next_cidx++;
                        if (BAM_CMATCH == bam_cigar_op(cigar[next_cidx]) || BAM_CEQUAL == bam_cigar_op(cigar[next_cidx]) || BAM_CDIFF == bam_cigar_op(cigar[next_cidx])
                                || BAM_CDEL == bam_cigar_op(cigar[next_cidx]) || BAM_CREF_SKIP == bam_cigar_op(cigar[next_cidx])) {
                            next_rpos += bam_cigar_oplen(cigar[next_cidx]);
                        }
                    }
                    uvc1_qual_t qfromBAQl = baq_offsetarr.getByPos(rpos) - baq_offsetarr.getByPos(prev_rpos);
                    uvc1_qual_t qfromBAQr = baq_offsetarr.getByPos(next_rpos) - baq_offsetarr.getByPos(rpos + cigar_oplen);
                    uvc1_qual_t qfromBAQ = MAX3(delFAQ, qfromBQ1, MIN(qfromBAQl, qfromBAQr));
                    incvalue = non_neg_minus(MIN3(qfromBQ2, qfromBAQ, phredvalue + symboltype2addPhred[LINK_SYMBOL]), micro_indel_penal) + 1;
                }
                if (nbases2end >= paramset.indel_filter_edge_dist) {
                    AlignmentSymbol symbol = delLenToSymbol(dellen, aln);
                    this->template inc<TUpdateType>(rpos, symbol, MAX(SIGN2UNSIGN(1), incvalue), aln);
                    if (TIsBiasUpdated) {
                        dealwith_segbias<true>(
                                MAX(SIGN2UNSIGN(1), incvalue),
                                rpos,
                                seg_format_depth_sets.getRefByPos(rpos)[symbol],
                                symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                                seg_format_thres_sets.getByPos(rpos),
                                aln,
                                xm1500,
                                go1500,
                                bm1500s[symbol],
                                baq_offsetarr,
                                baq_offsetarr2,
                                
                                cigar_op,
                                cigar_oplen,
                                10000,
                                dflag,
                                
                                paramset,
                                0);
                    }
                    this->incDel(rpos, cigar_oplen, delLenToSymbol(dellen, aln), MAX(SIGN2UNSIGN(1), incvalue));
#if 1 
// The definition of non-ref for deletion is not clearly defined, this piece of code makes sure that deletion identity conforms to the VCF specs.
                    for (uvc1_refgpos_t rpos2 = rpos; rpos2 < MIN(rpos + UNSIGN2SIGN(cigar_oplen), rend); rpos2++) {
                        for (AlignmentSymbol s : std::array<AlignmentSymbol, 2> {{ BASE_NN, LINK_NN }} ) {
                            const auto p = ((BASE_NN == s) ? (rpos2) : (rpos2 + 1));
                            if (p >= rend) { continue ; }
                            this->template inc<TUpdateType>(p, s, MAX(SIGN2UNSIGN(1), incvalue), aln);
                            if (TIsBiasUpdated) {
                                if (indel_rposs[indel_rposs_idx] <= rpos) {
                                    indel_rposs_idx++;
                                }
                                uint32_t prev_indel_rpos = indel_rposs[indel_rposs_idx-1];
                                uint32_t next_indel_rpos = indel_rposs[indel_rposs_idx];
                                dealwith_segbias<true>(
                                        MAX(SIGN2UNSIGN(1), incvalue),
                                        p,
                                        seg_format_depth_sets.getRefByPos(p)[s],
                                        symbol_to_VQ_format_tag_sets.getRefByPos(p)[s],
                                        seg_format_thres_sets.getByPos(p),
                                        aln,
                                        xm1500,
                                        go1500,
                                        bm1500s[s],
                                        baq_offsetarr,
                                        baq_offsetarr2,
                                        
                                        cigar_op,
                                        cigar_oplen,
                                        MIN(rpos - prev_indel_rpos, next_indel_rpos - rpos),
                                        dflag,

                                        paramset,
                                        0);
                            }
                        }
                    }
                }
#endif
}
                rpos += cigar_oplen;
            } else {
                process_cigar(qpos, rpos, cigar_op, cigar_oplen);
            } 
        }
        return 0;
    }
    
    template <ValueType TUpdateType = BASE_QUALITY_MAX, bool TIsBiasUpdated = false>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByRead1Aln(
            const std::vector<bam1_t *> & aln_vec,
            
            uvc1_refgpos_t region_offset,
            const auto & region_symbolvec,
            const auto & region_repeatvec,
            const auto & baq_offsetarr,
            const auto & baq_offsetarr2,
            
            auto & symbol_to_seg_format_depth_sets,
            auto & symbol_to_VQ_format_tag_sets,
            const auto & seg_format_prep_sets,
            const auto & seg_format_thres_sets,
            
            const auto dflag,
            const auto & paramset,
            const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
        for (const bam1_t *aln : aln_vec) {
            if (SEQUENCING_PLATFORM_IONTORRENT == paramset.inferred_sequencing_platform) {
                this->template updateByAln<true, TUpdateType, TIsBiasUpdated>(
                        aln, 
                        region_offset, 
                        region_symbolvec, 
                        region_repeatvec, 
                        baq_offsetarr,
                        baq_offsetarr2,
                        
                        symbol_to_seg_format_depth_sets,
                        symbol_to_VQ_format_tag_sets,
                        seg_format_prep_sets,
                        seg_format_thres_sets,

                        dflag,
                        paramset,
                        0);
            } else {
                this->template updateByAln<false, TUpdateType, TIsBiasUpdated>(
                        aln, 
                        region_offset, 
                        region_symbolvec, 
                        region_repeatvec, 
                        baq_offsetarr,
                        baq_offsetarr2,
                        
                        symbol_to_seg_format_depth_sets,
                        symbol_to_VQ_format_tag_sets,
                        seg_format_prep_sets,
                        seg_format_thres_sets,
                        
                        dflag,
                        paramset,
                        0);
            }
        }
        return 0;
    }
};

typedef GenericSymbol2CountCoverage<Symbol2Count> Symbol2CountCoverage; 
typedef GenericSymbol2CountCoverage<std::array<std::string, NUM_ALIGNMENT_SYMBOLS>> Symbol2CountCoverageString; 

struct Symbol2CountCoverageSet {
    uvc1_refgpos_t tid;
    uvc1_refgpos_t incluBegPosition;
    uvc1_refgpos_t excluEndPosition;
    std::string refstring;
    
    SegFormatPrepSets seg_format_prep_sets;
    SegFormatThresSets seg_format_thres_sets;
    Symbol2SegFormatDepthSets symbol_to_seg_format_depth_sets;
    std::array<Symbol2FragFormatDepthSets, 2> symbol_to_frag_format_depth_sets; 
    std::array<Symbol2FamFormatDepthSets, 2> symbol_to_fam_format_depth_sets_2strand;
    Symbol2DuplexFormatDepthSets symbol_to_duplex_format_depth_sets;
    Symbol2VQFormatTagSets symbol_to_VQ_format_tag_sets;
    
    std::array<Symbol2Bucket2CountCoverage, 2> dedup_ampDistr;
    Symbol2CountCoverageString additional_note;
    
    Symbol2CountCoverageSet(uvc1_refgpos_t t, uvc1_refgpos_t beg, uvc1_refgpos_t end):
        tid(t), 
        incluBegPosition(beg), 
        excluEndPosition(end),
        seg_format_prep_sets(SegFormatPrepSets(t, beg, end)),
        seg_format_thres_sets(SegFormatThresSets(t, beg, end)),
        symbol_to_seg_format_depth_sets(Symbol2SegFormatDepthSets(t, beg, end)),
        symbol_to_frag_format_depth_sets({{Symbol2FragFormatDepthSets(t, beg, end), Symbol2FragFormatDepthSets(t, beg, end)}}),
        symbol_to_fam_format_depth_sets_2strand({{Symbol2FamFormatDepthSets(t, beg, end), Symbol2FamFormatDepthSets(t, beg, end)}}),
        symbol_to_duplex_format_depth_sets(Symbol2DuplexFormatDepthSets(t, beg, end)),
        symbol_to_VQ_format_tag_sets(Symbol2VQFormatTagSets(t, beg, end)),
        
        dedup_ampDistr({Symbol2Bucket2CountCoverage(t, beg, end), Symbol2Bucket2CountCoverage(t, beg, end)}),
        additional_note(Symbol2CountCoverageString(t, beg, end)) 
    {
        assert(beg < end);
    };
    
    uvc1_refgpos_t
    getUnifiedIncluBegPosition() const {
        return incluBegPosition;
    };
    uvc1_refgpos_t
    getUnifiedExcluEndPosition() const {
        return excluEndPosition;
    };
    
    int 
    updateByAlns3UsingBQ(
            std::map<std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>>, std::array<uvc1_readnum_t, 2>> & mutform2count4map,
            const std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, uvc1_flag_t>> & alns3, 
            
            const std::basic_string<AlignmentSymbol> & region_symbolvec,
            auto & region_repeatvec,
            const auto & baq_offsetarr,
            const auto & baq_offsetarr2,

            const CommandLineArgs & paramset,
            const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
        
        Symbol2CountCoverage bg_seg_bqsum_conslogo(tid, this->getUnifiedIncluBegPosition(), this->getUnifiedExcluEndPosition());
        for (const auto & alns2pair2dflag : alns3) {
            for (int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair2dflag.first[strand];
                for (const auto & alns1 : alns2) {
                    for (const bam1_t *aln : alns1) {
                        update_seg_format_prep_sets_by_aln(
                                this->seg_format_prep_sets,
                                aln,
                                region_repeatvec,
                                baq_offsetarr,
                                this->getUnifiedIncluBegPosition(),
                                alns2pair2dflag.second,
                                region_symbolvec,
                                
                                paramset,
                                0);
                    }
                }
            }
        }
        update_seg_format_thres_from_prep_sets(
                region_repeatvec,
                this->seg_format_thres_sets,
                this->seg_format_prep_sets,
                
                paramset,
                0);
        int n_updates = 0;
        for (const auto & alns2pair2dflag : alns3) {
            for (int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair2dflag.first[strand];
                for (const auto & alns1 : alns2) {
                    n_updates += alns1.size();
                    bg_seg_bqsum_conslogo.updateByRead1Aln<SYMBOL_COUNT_SUM, true>(
                            alns1, 
                            
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec, 
                            region_repeatvec,
                            baq_offsetarr,
                            baq_offsetarr2,

                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,
                            alns2pair2dflag.second,
                            
                            paramset,
                            0);
                }
            }
        }
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            for (int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                for (const auto & alns1 : alns2) {
                    uvc1_refgpos_t tid2, beg2, end2;
                    fillTidBegEndFromAlns1(tid2, beg2, end2, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid, beg2, end2);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln(
                            alns1, 
                            
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec,
                            region_repeatvec,
                            baq_offsetarr,
                            baq_offsetarr2,

                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,
                            
                            alns2pair2dflag.second,
                            paramset,
                            0);
                    uvc1_qual_t normMQ = 0;
                    for (const bam1_t *aln : alns1) {
                        normMQ = MAX(normMQ, aln->core.qual);
                    }
                    std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>> pos_symbol_string;

                    size_t tlen = read_ampBQerr_fragWithR1R2.getExcluEndPosition() - read_ampBQerr_fragWithR1R2.getIncluBegPosition();
                    // 1 means is covered, 2 means has mut, 4 means is near mut.
                    std::vector<int8_t> cov_mut_vec(tlen, 0);
                    std::vector<AlignmentSymbol> cov_mut_base_symbol_vec(tlen, END_ALIGNMENT_SYMBOLS);
                    std::vector<AlignmentSymbol> cov_mut_link_symbol_vec(tlen, END_ALIGNMENT_SYMBOLS);
                    
                    for (auto epos = read_ampBQerr_fragWithR1R2.getIncluBegPosition(); epos < read_ampBQerr_fragWithR1R2.getExcluEndPosition(); epos++) {
                        for (SymbolType symboltype : SYMBOL_TYPES_IN_VCF_ORDER) {
                            AlignmentSymbol con_symbol;
                            uvc1_qual_t con_count, tot_count;
                            if (LINK_SYMBOL == symboltype) {
                                read_ampBQerr_fragWithR1R2.getByPos(epos).template fillConsensusCounts<true >(con_symbol, con_count, tot_count, symboltype);
                            } else {
                                read_ampBQerr_fragWithR1R2.getByPos(epos).template fillConsensusCounts<false>(con_symbol, con_count, tot_count, symboltype); 
                            }
                            assert (con_count * 2 >= tot_count);
                            if (0 == tot_count) { continue; }
                            
                            uvc1_qual_t max_qual = 8 + get_avgBQ(bg_seg_bqsum_conslogo, symbol_to_seg_format_depth_sets, epos, con_symbol);
                            uvc1_qual_t phredlike = MIN(con_count * 2 - tot_count, max_qual);
                            
                            int pbucket = max_qual - phredlike;
                            if (pbucket < -8) {
                                std::string qnames;
                                for (const auto *aln : alns1) { 
                                    qnames += std::string("/") + bam_get_qname(aln) + "/" + std::to_string(aln->core.tid) + "/" + std::to_string(aln->core.pos); 
                                }
                                LOG(logWARNING) << "The qname " << qnames << " has base quality " << phredlike << " at position " << epos << " which is higher than " << max_qual;
                            }
                            pbucket = MAX(0, pbucket);
                            if (pbucket < NUM_BUCKETS) {
                                this->dedup_ampDistr[0].getRefByPos(epos).incSymbolBucketCount(con_symbol, pbucket, 1);
                            }
                            
                            this->symbol_to_frag_format_depth_sets[strand].getRefByPos(epos)[con_symbol][FRAG_bDP] += 1;
                            this->symbol_to_VQ_format_tag_sets.getRefByPos(epos)[con_symbol][VQ_bMQ] += (normMQ * normMQ) / SQR_QUAL_DIV;
                            
                            if (isSymbolIns(con_symbol)) {
                                posToIndelToCount_updateByConsensus(this->symbol_to_frag_format_depth_sets[strand].getRefPosToIseqToData(con_symbol), 
                                        read_ampBQerr_fragWithR1R2.getPosToIseqToData(con_symbol), epos, 1);
                            }
                            if (isSymbolDel(con_symbol)) {
                                posToIndelToCount_updateByConsensus(this->symbol_to_frag_format_depth_sets[strand].getRefPosToDlenToData(con_symbol),
                                        read_ampBQerr_fragWithR1R2.getPosToDlenToData(con_symbol), epos, 1);
                            }
                            AlignmentSymbol refsymbol = region_symbolvec[epos-this->getUnifiedIncluBegPosition()]; 
                            cov_mut_vec[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()] |= 0x1;
                            if (areSymbolsMutated(refsymbol, con_symbol) && (LINK_SYMBOL == symboltype || phredlike >= paramset.bias_thres_highBQ)) {
                                pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
                                cov_mut_vec[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()] |= 0x2;
                            }
                            if (LINK_SYMBOL == symboltype) {
                                cov_mut_link_symbol_vec[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()] = con_symbol;
                            } else {
                                cov_mut_base_symbol_vec[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()] = con_symbol;
                            }
                        }
                    }
                    if (pos_symbol_string.size() > 1) {
                        mutform2count4map.insert(std::make_pair(pos_symbol_string, std::array<uvc1_readnum_t, 2>({0, 0})));
                        mutform2count4map[pos_symbol_string][strand]++;
                    }
                    for (size_t i = 0; i < cov_mut_vec.size(); i++) {
                        if (cov_mut_vec[i] & 0x2) {
                            for (size_t j = i - paramset.syserr_mut_region_n_bases; j < i + paramset.syserr_mut_region_n_bases + 1; j++) {
                                if (j < cov_mut_vec.size()) {
                                    cov_mut_vec[j] |= 0x4;
                                }
                            }
                        }
                    }
                    uvc1_base_t n_cov_positions = 0;
                    uvc1_base_t n_near_mut_positions = 0;
                    for (size_t i = 0; i < cov_mut_vec.size(); i++) {
                        if ((cov_mut_vec[i]) & 0x1) {
                            n_cov_positions++;
                            if ((cov_mut_vec[i]) & 0x4) { n_near_mut_positions++; }
                        }
                    }
                    const auto b10xSeqTlen = n_cov_positions;
                    const auto b10xSeqTNevents = n_near_mut_positions;
                    for (auto epos = read_ampBQerr_fragWithR1R2.getIncluBegPosition(); epos < read_ampBQerr_fragWithR1R2.getExcluEndPosition(); epos++) {
                        const auto ivec = epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition();
                        for (const auto con_symbol : std::array<AlignmentSymbol, 2> {{ cov_mut_base_symbol_vec[ivec], cov_mut_link_symbol_vec[ivec] }} ) {
                            if (con_symbol != END_ALIGNMENT_SYMBOLS) {
                                this->symbol_to_frag_format_depth_sets[strand].getRefByPos(epos)[con_symbol][FRAG_bTA] += b10xSeqTlen;
                                this->symbol_to_frag_format_depth_sets[strand].getRefByPos(epos)[con_symbol][FRAG_bTB] += b10xSeqTNevents;
                            }
                        }
                    }
                }
            }
        }
        assert(this->symbol_to_seg_format_depth_sets.getIncluBegPosition() == this->dedup_ampDistr[0].getIncluBegPosition());
        assert(this->symbol_to_seg_format_depth_sets.getExcluEndPosition() == this->dedup_ampDistr[0].getExcluEndPosition());
        for (auto epos = this->getUnifiedIncluBegPosition(); epos < this->getUnifiedExcluEndPosition(); epos++) {
            for (SymbolType symboltype : SYMBOL_TYPE_ARR) {
                for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symboltype]) {
                    const auto totDP = 
                            formatSumBySymbolType(this->symbol_to_frag_format_depth_sets[0].getByPos(epos), 
                            symboltype, FRAG_bDP) + 
                            formatSumBySymbolType(this->symbol_to_frag_format_depth_sets[1].getByPos(epos), 
                            symboltype, FRAG_bDP); 
                    uvc1_qual_t max_qual = 8 + get_avgBQ(bg_seg_bqsum_conslogo, symbol_to_seg_format_depth_sets, epos, symbol);
                    uvc1_qual_t maxvqual = 0;
                    uvc1_readnum_t argmaxAD = 0;
                    uvc1_qual_t argmaxBQ = 0;
                    infer_max_qual_assuming_independence(
                            maxvqual,
                            argmaxAD,
                            argmaxBQ,
                            max_qual,
                            1, 
                            this->dedup_ampDistr[0].getByPos(epos).getSymbolCounts(symbol),
                            totDP,
                            0);
                    this->symbol_to_VQ_format_tag_sets.getRefByPos(epos)[symbol][VQ_bIAQb] += maxvqual;
                    this->symbol_to_VQ_format_tag_sets.getRefByPos(epos)[symbol][VQ_bIADb] += argmaxAD;
                    this->symbol_to_VQ_format_tag_sets.getRefByPos(epos)[symbol][VQ_bIDQb] += argmaxBQ;
                }
            }
            this->dedup_ampDistr[0].getRefByPos(epos).clearSymbolBucketCount();
        }
        return 0;
    };
    
    int
    updateByAlns3UsingFQ(
            std::map<std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>>, std::array<uvc1_readnum_t, 2>> & mutform2count4map,
            const auto & alns3, 
            
            const std::basic_string<AlignmentSymbol> & region_symbolvec,
            const auto & region_repeatvec,
            const auto & baq_offsetarr,
            const auto & baq_offsetarr2,
            
            const CommandLineArgs & paramset,
            const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
        
        const PhredMutationTable sscs_mut_table(
                paramset.fam_phred_sscs_transition_CG_TA,
                paramset.fam_phred_sscs_transition_AT_GC,
                paramset.fam_phred_sscs_transversion_CG_AT,
                paramset.fam_phred_sscs_transversion_other,
                paramset.fam_phred_sscs_indel_open,
                paramset.fam_phred_sscs_indel_ext);
        
        size_t niters = 0;
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            niters++;
            assert (alns2pair[0].size() != 0 || alns2pair[1].size() != 0);
            for (int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                if (alns2.size() == 0) { continue; }
                uvc1_refgpos_t tid2, beg2, end2;
                fillTidBegEndFromAlns2(tid2, beg2, end2, alns2);
                Symbol2CountCoverage read_family_con_ampl(tid2, beg2, end2); 
                for (const auto & alns1 : alns2) {
                    uvc1_refgpos_t tid1, beg1, end1;
                    fillTidBegEndFromAlns1(tid1, beg1, end1, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid1, beg1, end1);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln(
                            alns1, 
                            
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec,
                            region_repeatvec,
                            baq_offsetarr,
                            baq_offsetarr2,
                            
                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,
                            alns2pair2dflag.second,

                            paramset,
                            0);
                    
                    read_family_con_ampl.updateByFiltering(read_ampBQerr_fragWithR1R2, std::array<uvc1_qual_t, NUM_SYMBOL_TYPES> {{ 
                            paramset.fam_thres_highBQ_snv, 0
                    }});
                }
                for (auto epos = read_family_con_ampl.getIncluBegPosition(); epos < read_family_con_ampl.getExcluEndPosition(); epos++) {
                    const auto & con_ampl_symbol2count = read_family_con_ampl.getByPos(epos);
                    for (SymbolType symboltype : SYMBOL_TYPE_ARR) {
                        AlignmentSymbol con_symbol;
                        uvc1_qual_t con_count, tot_count;
                        con_ampl_symbol2count.fillConsensusCounts(con_symbol, con_count, tot_count, symboltype);
                        if (0 == tot_count) { continue ; }
                        this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP12] += 1;
                        if (1 == tot_count) {
                            this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_c1DP] += 1;
                        }
                        if (paramset.fam_thres_dup1add <= tot_count && (con_count * 100 >= tot_count * paramset.fam_thres_dup1perc)) {
                            this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP2] += 1;
                        }
                        if (paramset.fam_thres_dup2add <= tot_count && (con_count * 100 >= tot_count * paramset.fam_thres_dup2perc)) {
                            this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP3] += 1;
                        }
                        
                        if (isSymbolIns(con_symbol)) {
                            posToIndelToCount_updateByConsensus(
                                    this->symbol_to_fam_format_depth_sets_2strand[strand].getRefPosToIseqToData(con_symbol),
                                    read_family_con_ampl.getPosToIseqToData(con_symbol), epos, 1);
                        }
                        if (isSymbolDel(con_symbol)) {
                            posToIndelToCount_updateByConsensus(
                                    this->symbol_to_fam_format_depth_sets_2strand[strand].getRefPosToDlenToData(con_symbol),
                                    read_family_con_ampl.getPosToDlenToData(con_symbol), epos, 1);
                        }
                        // This is one round of bootstrapping for EM. We can use a full EM algorithm to improve this estimator, but it is probably overkill
                        // 0 means no coverage, 1 means no error correction, 2 means low quality family if the symbols disagree with each other, 
                        // 3 means up to 1 erroneous basecall is tolerated, 4 means up to 2 erroneous basecalls are tolerated
                        const auto fam_thres_emperr_all_flat = (isSymbolSubstitution(con_symbol) ? 
                                paramset.fam_thres_emperr_all_flat_snv : paramset.fam_thres_emperr_all_flat_indel);
                        const auto fam_thres_emperr_con_perc = (isSymbolSubstitution(con_symbol) ? 
                                paramset.fam_thres_emperr_con_perc_snv : paramset.fam_thres_emperr_con_perc_indel);
                        if (tot_count < fam_thres_emperr_all_flat) { continue; } 
                        if (con_count * 100 < tot_count * fam_thres_emperr_con_perc) { continue; }
                        for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symboltype]) {
                            if (con_symbol != symbol) {
                                this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDPm] += 
                                        con_ampl_symbol2count.getSymbolCount(symbol);
                                this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDPM] += tot_count;
                            }
                        }
                    }
                }
            }
        }
        niters = 0;
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            niters++;
            uvc1_refgpos_t tid2, beg2, end2;
            tid2 = 0; beg2 = INT32_MAX; end2 = 0; bool initialized = false;
            assert (alns2pair[0].size() != 0 || alns2pair[1].size() != 0);
            if (alns2pair[0].size() > 0) { fillTidBegEndFromAlns2(tid2, beg2, end2, alns2pair[0], initialized); initialized = true; }
            if (alns2pair[1].size() > 0) { fillTidBegEndFromAlns2(tid2, beg2, end2, alns2pair[1], initialized); initialized = true; }
            Symbol2CountCoverage read_duplex_amplicon(tid2, beg2, end2);
            for (int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                if (alns2.size() == 0) { continue; }
                uvc1_refgpos_t tid2, beg2, end2;
                fillTidBegEndFromAlns2(tid2, beg2, end2, alns2);
                Symbol2CountCoverage read_family_mmm_ampl(tid2, beg2, end2);
                Symbol2CountCoverage read_family_con_ampl(tid2, beg2, end2);
                for (const auto & aln_vec : alns2) {
                    uvc1_refgpos_t tid1, beg1, end1;
                    fillTidBegEndFromAlns1(tid1, beg1, end1, aln_vec);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid1, beg1, end1);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln(
                            aln_vec,
                            
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec, 
                            region_repeatvec,
                            baq_offsetarr,      
                            baq_offsetarr2,

                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,
                            alns2pair2dflag.second,

                            paramset,
                            0);
                    // The line below is similar to : read_family_mmm_ampl.updateByConsensus<SYMBOL_COUNT_SUM>(read_ampBQerr_fragWithR1R2);
                    read_family_con_ampl.updateByFiltering(read_ampBQerr_fragWithR1R2, std::array<uvc1_qual_t, NUM_SYMBOL_TYPES> {{ 
                            paramset.fam_thres_highBQ_snv, 0
                    }});
                    read_family_mmm_ampl.updateByMajorMinusMinor(read_ampBQerr_fragWithR1R2, false); 
                }
                if ((0x2 == (alns2pair2dflag.second & 0x2)) && alns2pair[0].size() > 0 && alns2pair[1].size() > 0) { // is duplex
                    read_duplex_amplicon.template updateByConsensus<false>(read_family_con_ampl);
                }
                std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>> pos_symbol_string;
                for (auto epos = read_family_mmm_ampl.getIncluBegPosition(); 
                        epos < read_family_mmm_ampl.getExcluEndPosition(); 
                        epos++) {
                    for (SymbolType symboltype : SYMBOL_TYPES_IN_VCF_ORDER) {
                        
                        AlignmentSymbol con_symbol; // assign with the value AlignmentSymbol(NUM_ALIGNMENT_SYMBOLS) to flag for error
                        uvc1_qual_t con_sumBQs, tot_sumBQs;
                        read_family_mmm_ampl.getRefByPos(epos).fillConsensusCounts(con_symbol, con_sumBQs, tot_sumBQs, symboltype);
                        if (0 == tot_sumBQs) { continue; }
                        uvc1_readnum_t con_nfrags = read_family_con_ampl.getByPos(epos).getSymbolCount(con_symbol);
                        uvc1_readnum_t tot_nfrags = read_family_con_ampl.getByPos(epos).sumBySymbolType(symboltype);
                        
                        this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP1] += 1; // in some rare cases, cDP1 > cDP0 
                        
                        const uvc1_qual_t avgBQ = ((0 == tot_nfrags) ? 1 : (con_sumBQs / tot_nfrags));
                        const uvc1_readnum_t majorcount = this->symbol_to_fam_format_depth_sets_2strand[strand].getByPos(epos).at(con_symbol).at(FAM_cDPM);
                        const uvc1_readnum_t minorcount = this->symbol_to_fam_format_depth_sets_2strand[strand].getByPos(epos).at(con_symbol).at(FAM_cDPm);
                        const double prior_weight = 1.0 / (minorcount + 1.0);
                        const double realphred = prob2realphred((minorcount + prior_weight) / 
                                (majorcount + minorcount + prior_weight / phred2prob(avgBQ)));
                        // Compute empirical error assuming statistical independence and majority vote by frag count.
                        // If no majority is obtained, then confam_qual is set to one.
                        const uvc1_qual_t indep_frag_phred = round(((con_nfrags * 2) - tot_nfrags) * realphred);
                        uvc1_qual_t confam_qual = 0;
                        if (LINK_SYMBOL == symboltype) {
                            confam_qual = MAX(1, (MIN(indep_frag_phred,
                                    // PCR error of the first cycle + error in consensus selection
                                    (uvc1_qual_t)paramset.fam_phred_indel_inc_before_barcode_labeling + (uvc1_qual_t)round(realphred))));
                        } else {
                            confam_qual = MAX(1, MIN(indep_frag_phred,
                                    ((con_sumBQs * 2) - tot_sumBQs)));
                        }
                        
                        AlignmentSymbol ref_symbol = region_symbolvec[epos - this->getUnifiedIncluBegPosition()]; 
                        if (areSymbolsMutated(ref_symbol, con_symbol) && (LINK_SYMBOL == symboltype || confam_qual >= paramset.bias_thres_highBQ)) {
                            pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
                        }
                        uvc1_qual_t max_qual = sscs_mut_table.toPhredErrRate(ref_symbol, con_symbol);
                        uvc1_qual_t confam_qual2 = MIN(confam_qual, max_qual);
                        if (tot_nfrags >= paramset.fam_thres_dup1add) {
                            int pbucket = (max_qual - confam_qual2 + 2) / 4;
                            this->dedup_ampDistr[strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pbucket, 1);
                        }
                    }
                }
                
                if (pos_symbol_string.size() > 1) {
                    mutform2count4map.insert(std::make_pair(pos_symbol_string, std::array<uvc1_readnum_t, 2>({0, 0})));
                    mutform2count4map[pos_symbol_string][strand]++;
                }
            }
            if ((0x2 == (alns2pair2dflag.second & 0x2)) && alns2pair[0].size() > 0 && alns2pair[1].size() > 0) { // is duplex
                for (auto epos = read_duplex_amplicon.getIncluBegPosition(); epos < read_duplex_amplicon.getExcluEndPosition(); epos++) {
                    for (SymbolType symboltype : SYMBOL_TYPE_ARR) {
                        AlignmentSymbol con_symbol;
                        uvc1_readnum_t con_count, tot_count;
                        read_duplex_amplicon.getRefByPos(epos).fillConsensusCounts(con_symbol, con_count, tot_count, symboltype);
                        assert (tot_count <= 2 || !fprintf(stderr, "%d <= 2 failed for duplex family, a duplex family is supported by two single-strand families!\n", tot_count));
                        if (0 < tot_count) {
                            this->symbol_to_duplex_format_depth_sets.getRefByPos(epos)[con_symbol][DUPLEX_dDP1] += 1;
                        }
                        if (1 < tot_count) {
                            this->symbol_to_duplex_format_depth_sets.getRefByPos(epos)[con_symbol][DUPLEX_dDP2] += 1;
                        }
                    }
                }
            }
        }
        for (int strand = 0; strand < 2; strand++) {
            assert(this->symbol_to_fam_format_depth_sets_2strand[strand].getExcluEndPosition() == this->dedup_ampDistr[strand].getExcluEndPosition());
            assert(this->symbol_to_fam_format_depth_sets_2strand[strand].getIncluBegPosition() == this->dedup_ampDistr[strand].getIncluBegPosition()); 
            auto VQ_cIAQ = (strand ? VQ_cIAQr : VQ_cIAQf);
            auto VQ_cIAD = (strand ? VQ_cIADr : VQ_cIADf);
            auto VQ_cIDQ = (strand ? VQ_cIDQr : VQ_cIDQf);
            for (uvc1_refgpos_t epos = this->getUnifiedIncluBegPosition(); epos < this->getUnifiedExcluEndPosition(); epos++) {
                for (SymbolType symboltype : SYMBOL_TYPE_ARR) {
                    AlignmentSymbol ref_symbol = region_symbolvec[epos - this->getUnifiedIncluBegPosition()];
                    const uvc1_readnum_t totDP = formatSumBySymbolType(this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos), symboltype, FAM_cDP1);
                    for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symboltype]) {
                        const uvc1_qual_t max_qual = sscs_mut_table.toPhredErrRate(ref_symbol, symbol);
                        uvc1_qual_t maxvqual = 0;
                        uvc1_readnum_t argmaxAD = 0;
                        uvc1_qual_t argmaxBQ = 0;
                        infer_max_qual_assuming_independence(
                                maxvqual,
                                argmaxAD,
                                argmaxBQ,
                                max_qual,
                                4,
                                dedup_ampDistr[strand].getByPos(epos).getSymbolCounts(symbol),
                                totDP,
                                0);
                        this->symbol_to_VQ_format_tag_sets.getRefByPos(epos)[symbol][VQ_cIAQ] += maxvqual;
                        this->symbol_to_VQ_format_tag_sets.getRefByPos(epos)[symbol][VQ_cIAD] += argmaxAD;
                        this->symbol_to_VQ_format_tag_sets.getRefByPos(epos)[symbol][VQ_cIDQ] += argmaxBQ;
                        if (paramset.should_add_note) {
                            this->additional_note.getRefByPos(epos)[symbol] += std::string("//fq-distr/");
                            for (size_t idx = 0; idx < dedup_ampDistr[strand].getByPos(epos).getSymbolCounts(symbol).size(); idx++) {
                                this->additional_note.getRefByPos(epos)[symbol] += std::to_string(idx) + "/";
                                this->additional_note.getRefByPos(epos)[symbol] += std::to_string(dedup_ampDistr[strand].getByPos(epos).getSymbolCounts(symbol)[idx]) + "/";
                            }
                            this->additional_note.getRefByPos(epos)[symbol] += std::string("//");
                        }
                    }
                }
                this->dedup_ampDistr[strand].getRefByPos(epos).clearSymbolBucketCount();
            }
        }
        return 0;
    };
    
    int 
    updateHapMap(std::map<std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>>, std::array<uvc1_readnum_t, 2>> & mutform2count4map, 
            const auto & tsum_depth, const auto read_count_field_id, int max_ploidy = 2+1) {
        for (auto it = mutform2count4map.begin(); it != mutform2count4map.end();) {
            std::basic_string<std::pair<uvc1_readnum_t, AlignmentSymbol>> mutform = it->first;
            auto counts = it->second;
            std::vector<uvc1_readnum_t> dsADs;
            dsADs.reserve(mutform.size() + 1);
            dsADs.push_back(0);
            for (std::pair<uvc1_refgpos_t, AlignmentSymbol> simplemut : mutform) {
                uvc1_readnum_t ad0 = tsum_depth.at(0).getByPos(simplemut.first)[simplemut.second][read_count_field_id];
                uvc1_readnum_t ad1 = tsum_depth.at(1).getByPos(simplemut.first)[simplemut.second][read_count_field_id];
                dsADs.push_back(ad0 + ad1);
            }
            const auto med = MEDIAN(dsADs);
            if ((counts[0] +counts[1]) * (max_ploidy) <= med) {
                mutform2count4map.erase(it++);
            } else {
                it++;
            }
        }
        return 0;
    };
    
    int 
    updateByRegion3Aln(
            std::map<std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>>, std::array<uvc1_readnum_t, 2>> & mutform2count4map_bq,
            std::map<std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>>, std::array<uvc1_readnum_t, 2>> & mutform2count4map_fq,
            const std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, uvc1_flag_t>> & alns3, 
            
            const std::string & refstring,
            std::vector<RegionalTandemRepeat> & region_repeatvec,
            const CoveredRegion<uvc1_qual_big_t> & baq_offsetarr,
            const CoveredRegion<uvc1_qual_big_t> & baq_offsetarr2,

            const auto & paramset,
            const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
        
        std::basic_string<AlignmentSymbol> ref_symbol_string = string2symbolseq(refstring);
        
        updateByAlns3UsingBQ(
                mutform2count4map_bq, 
                alns3, 
                ref_symbol_string,

                region_repeatvec,
                baq_offsetarr,
                baq_offsetarr2,
                
                paramset,
                0);
        updateHapMap(mutform2count4map_bq, this->symbol_to_frag_format_depth_sets, FRAG_bDP);
        updateByAlns3UsingFQ(
                mutform2count4map_fq, 
                alns3,
                
                ref_symbol_string, 
                region_repeatvec,
                baq_offsetarr,
                baq_offsetarr2,
                
                paramset,
                0);
        updateHapMap(mutform2count4map_fq, this->symbol_to_fam_format_depth_sets_2strand, FAM_cDP1);
        return 0;
    };
};

void 
fill_symboltype_fmt(
        auto & fmtDP,
        const auto & symbol_to_abcd_format_depth_sets,
        const auto format_field,
        const auto refpos,
        const SymbolType symboltype, 
        const AlignmentSymbol refsymbol IGNORE_UNUSED_PARAM) {

    const auto symbolNN = SYMBOL_TYPE_TO_AMBIG[symboltype];
    fmtDP[0] = formatSumBySymbolType(symbol_to_abcd_format_depth_sets.getByPos(refpos), symboltype, format_field);
    fmtDP[1] = symbol_to_abcd_format_depth_sets.getByPos(refpos)[symbolNN][format_field];
};

#define filla_symboltype_fmt(fmtDP, symbol_to_abcd_format_depth_sets, format_field, refpos, symboltype, refsymbol) { \
    const auto symbolNN = SYMBOL_TYPE_TO_AMBIG[symboltype]; \
    int64_t ret = 0; \
    for (const auto symbol : SYMBOL_TYPE_TO_SYMBOLS[symboltype]) { ret += (int64_t)symbol_to_abcd_format_depth_sets.getByPos(refpos)[symbol] . format_field; } \
    fmtDP[0] = ret; \
    fmtDP[1] = symbol_to_abcd_format_depth_sets.getByPos(refpos)[symbolNN] . format_field; \
};

void
fill_symbol_fmt(
        auto & fmtAD,
        const auto & symbol_to_abcd_format_depth_sets,
        const auto format_field,
        const auto refpos,
        const auto symbol,
        const int allele_idx) {
    if (0 == allele_idx) {
        fmtAD.clear();
    } else {
        assert(fmtAD.size() > 0);
    }
    const auto cnt = symbol_to_abcd_format_depth_sets.getByPos(refpos)[symbol][format_field];
    fmtAD.push_back(cnt);
};

#define filla_symbol_fmt(fmtAD, symbol_to_abcd_format_depth_sets, format_field, refpos, symbol, allele_idx) { \
    if (0 == allele_idx) { fmtAD.clear(); } \
    else { assert(fmtAD.size() > 0); } \
    const auto cnt = symbol_to_abcd_format_depth_sets.getByPos(refpos)[symbol] . format_field; \
    fmtAD.push_back(cnt); \
};

int
fill_symbol_VQ_fmts(
        auto & fmt,
        const auto & symbol_to_VQ_format_tag_sets,
        const uvc1_refgpos_t refpos,
        const AlignmentSymbol symbol,
        uvc1_qual_t minABQ,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    const int a = 0;
    
    const auto a2BQf = symbol_to_VQ_format_tag_sets.getByPos(refpos)[symbol][VQ_a2BQf];
    const auto a2BQr = symbol_to_VQ_format_tag_sets.getByPos(refpos)[symbol][VQ_a2BQr];
    uvc1_readnum_t aDPf = fmt.aDPff[a] + fmt.aDPrf[a];
    uvc1_readnum_t aDPr = fmt.aDPfr[a] + fmt.aDPrr[a];
    uvc1_readnum_t ADP = fmt.ADPff[0] + fmt.ADPrf[0] + fmt.ADPfr[0] + fmt.ADPrr[0];
    const uvc1_qual_t rssDPfBQ = (aDPf * sqrt(int64mul(a2BQf, SQR_QUAL_DIV) / MAX(1, aDPf)));
    const uvc1_qual_t rssDPrBQ = (aDPr * sqrt(int64mul(a2BQr, SQR_QUAL_DIV) / MAX(1, aDPr)));
    
    const uvc1_qual_t rssDPbBQ = ((aDPf + aDPr) * sqrt((a2BQf + a2BQr) * SQR_QUAL_DIV / MAX(1, aDPf + aDPr)));
    
    assert ((aDPf + aDPr) * 100 >= LAST(fmt.a2XM2));
    assert ((aDPf + aDPr) * 100 >= LAST(fmt.a2BM2));
    // TODO: what is the exact theory behind the following line of code to rescue hetero and homalt?
    uvc1_qual_t minABQa = minABQ - (uvc1_qual_t)(5 * 10.0 * mathsquare(MAX(0, ((aDPf + aDPr + 0.5) * 2.0 / (ADP + 1.0) - 1.0))));
    const uvc1_readnum_t dp10pc = 10;
    double sbratio = (double)(MAX(aDPf, aDPr) * 10 + dp10pc) / (double)(MIN(aDPf, aDPr) * 10 + dp10pc);
    minABQa += BETWEEN((uvc1_qual_t)mathsquare(sbratio) - paramset.syserr_BQ_sbratio_q_add, 0, paramset.syserr_BQ_sbratio_q_max);
    const uvc1_readnum_t xmratio = (paramset.syserr_BQ_xmratio_q_max * 10 * (aDPf + aDPr) / MAX(1, LAST(fmt.a2XM2)));
    const uvc1_readnum_t bmratio = (paramset.syserr_BQ_bmratio_q_max * 10 * (aDPf + aDPr) / MAX(1, LAST(fmt.a2BM2)));
    minABQa += BETWEEN(xmratio - paramset.syserr_BQ_xmratio_q_add, 0, paramset.syserr_BQ_xmratio_q_max) 
             + BETWEEN(bmratio - paramset.syserr_BQ_bmratio_q_add, 0, paramset.syserr_BQ_bmratio_q_max);
    if (paramset.should_add_note) {
        fmt.note += std::string("//minABQa/") 
                + std::to_string(minABQa) + "/" 
                + std::to_string(sbratio) + "/" 
                + std::to_string(xmratio) + "/" 
                + std::to_string(bmratio) + "//";
    }
    const uvc1_readnum_t m = paramset.syserr_BQ_strand_favor_mul;
    uvc1_qual_t a_BQ_syserr_qual_fw = (rssDPfBQ * m - minABQa * aDPf * m / 10 + rssDPrBQ - minABQa * aDPr / 10) / m;
    uvc1_qual_t a_BQ_syserr_qual_rv = (rssDPrBQ * m - minABQa * aDPr * m / 10 + rssDPfBQ - minABQa * aDPf / 10) / m;
    uvc1_qual_t a_BQ_syserr_qual_2d = (rssDPbBQ) - minABQa * (aDPf + aDPr) / 10;
    const uvc1_qual_t a_rmsBQ = (rssDPbBQ) / MAX(1, aDPf + aDPr);
    fill_symbol_fmt(fmt.bMQ,  symbol_to_VQ_format_tag_sets,  VQ_bMQ,  refpos, symbol, a);
    fmt.bMQ[a] = round(sqrt(fmt.bMQ[a] * SQR_QUAL_DIV / MAX(fmt.bDPf[a] + fmt.bDPr[a], 1)) + (double)(1.0 - FLT_EPSILON));
    
    const auto aBQQ = MAX(a_rmsBQ, paramset.syserr_BQ_prior + MAX3(a_BQ_syserr_qual_2d, a_BQ_syserr_qual_fw, a_BQ_syserr_qual_rv));
    clear_push(fmt.a2BQf, rssDPfBQ, a);
    clear_push(fmt.a2BQr, rssDPrBQ, a);
    clear_push(fmt.aBQ, a_rmsBQ, a);
    clear_push(fmt.aBQQ, aBQQ, a);
    
    fill_symbol_fmt(fmt.bIAQb, symbol_to_VQ_format_tag_sets, VQ_bIAQb, refpos, symbol, a);
    fill_symbol_fmt(fmt.bIADb, symbol_to_VQ_format_tag_sets, VQ_bIADb, refpos, symbol, a);
    fill_symbol_fmt(fmt.bIDQb, symbol_to_VQ_format_tag_sets, VQ_bIDQb, refpos, symbol, a);
    fill_symbol_fmt(fmt.cIAQf, symbol_to_VQ_format_tag_sets, VQ_cIAQf, refpos, symbol, a);
    fill_symbol_fmt(fmt.cIADf, symbol_to_VQ_format_tag_sets, VQ_cIADf, refpos, symbol, a);
    fill_symbol_fmt(fmt.cIDQf, symbol_to_VQ_format_tag_sets, VQ_cIDQf, refpos, symbol, a);
    fill_symbol_fmt(fmt.cIAQr, symbol_to_VQ_format_tag_sets, VQ_cIAQr, refpos, symbol, a);
    fill_symbol_fmt(fmt.cIADr, symbol_to_VQ_format_tag_sets, VQ_cIADr, refpos, symbol, a);
    fill_symbol_fmt(fmt.cIDQr, symbol_to_VQ_format_tag_sets, VQ_cIDQr, refpos, symbol, a);
    
    clear_push(fmt.VTD,  SYMBOL_TO_DESC_ARR[symbol], a);
    clear_push(fmt.VTI,  (int32_t)symbol, a);
    return 0;
}

std::array<uvc1_readnum_t, 2>
BcfFormat_symboltype_init(bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet12, 
        uvc1_refgpos_t refpos, 
        const SymbolType symboltype,
        const AlignmentSymbol refsymbol,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    const auto & p = symbol2CountCoverageSet12.seg_format_prep_sets.getByPos(refpos);
    fmt.APDP = {{ 
        p.segprep_a_dp, 
        
        p.segprep_a_near_ins_dp, 
        p.segprep_a_near_del_dp, 
        p.segprep_a_near_RTR_ins_dp, 
        p.segprep_a_near_RTR_del_dp, 
        
        p.segprep_a_pcr_dp,
        p.segprep_a_snv_dp,
        p.segprep_a_dnv_dp,
        p.segprep_a_highBQ_dp,
        
        p.segprep_a_near_pcr_clip_dp,
        p.segprep_a_near_long_clip_dp,
    }};
    
    fmt.APXM = {{ 
        p.segprep_a_XM1500,
        p.segprep_a_GO1500,
        p.segprep_a_qlen,

        p.segprep_a_GAPLEN,
        p.segprep_a_near_ins_pow2len, 
        p.segprep_a_near_del_pow2len,
        // no longer used afterwards
        
        p.segprep_a_near_ins_inv100len,
        p.segprep_a_near_del_inv100len,
    }};
    
    fmt.APLRID = {{
        p.segprep_a_near_ins_l_pow2len,
        p.segprep_a_near_ins_r_pow2len,
        p.segprep_a_near_del_l_pow2len,
        p.segprep_a_near_del_r_pow2len,
    }};

    fmt.APLRI = {{ p.segprep_a_LI, p.segprep_a_LIDP, p.segprep_a_RI, p.segprep_a_RIDP }};
    fmt.APLRP = {{ p.segprep_a_l_dist_sum, p.segprep_a_r_dist_sum, p.segprep_a_inslen_sum, p.segprep_a_dellen_sum }};
    
    const auto & t = symbol2CountCoverageSet12.seg_format_thres_sets.getByPos(refpos);

#if COMPILATION_TRY_HIGH_DEPTH_POS_BIAS // can be useful at very high depth
    fmt.APPB  = {{
        p.segprep_aa_l_ins_dist_x_wei,
        p.segprep_aa_l_ins_weight,
        p.segprep_aa_r_ins_dist_x_wei,
        p.segprep_aa_r_ins_weight,
        
        p.segprep_aa_l_del_dist_x_wei,
        p.segprep_aa_l_del_weight,
        p.segprep_aa_r_del_dist_x_wei,
        p.segprep_aa_r_del_weight,
    }};
#endif
#if COMPILATION_ENABLE_XMGOT 
    fmt.AXMT = {{ t.segthres_aXM2T, t.segthres_aXM2T }};
#endif

    fmt.ALRPxT = {{ t.segthres_aLPxT, t.segthres_aRPxT }};
    
    fmt.ALRIT = {{ t.segthres_aLI1T, t.segthres_aLI2T, t.segthres_aRI1T, t.segthres_aRI2T }};
    fmt.ALRIt = {{ t.segthres_aLI1t, t.segthres_aLI2t, t.segthres_aRI1t, t.segthres_aRI2t }};
    fmt.ALRPt = {{ t.segthres_aLP1t, t.segthres_aLP2t, t.segthres_aRP1t, t.segthres_aRP2t }};
    fmt.ALRBt = {{ t.segthres_aLB1t, t.segthres_aLB2t, t.segthres_aRB1t, t.segthres_aRB2t }};
    
    const auto & symbol_to_VQ_format_tag_sets = symbol2CountCoverageSet12.symbol_to_VQ_format_tag_sets;
    const auto & symbol_to_seg_format_depth_sets = symbol2CountCoverageSet12.symbol_to_seg_format_depth_sets;
    
    fill_symboltype_fmt(fmt.A1BQf, symbol_to_VQ_format_tag_sets,    VQ_a1BQf,  refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.A1BQr, symbol_to_VQ_format_tag_sets,    VQ_a1BQr,  refpos, symboltype, refsymbol);
    // this MQ can be used for addition and substraction.
    filla_symboltype_fmt(fmt.AMQs,  symbol_to_seg_format_depth_sets, seginfo_aMQs,  refpos, symboltype, refsymbol); 
    
    filla_symboltype_fmt(fmt.AP1,   symbol_to_seg_format_depth_sets, seginfo_aP1,   refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.AP2,   symbol_to_seg_format_depth_sets, seginfo_aP2,   refpos, symboltype, refsymbol);

    filla_symboltype_fmt(fmt.ADPff, symbol_to_seg_format_depth_sets, seginfo_aDPff, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ADPfr, symbol_to_seg_format_depth_sets, seginfo_aDPfr, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ADPrf, symbol_to_seg_format_depth_sets, seginfo_aDPrf, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ADPrr, symbol_to_seg_format_depth_sets, seginfo_aDPrr, refpos, symboltype, refsymbol);
    
    filla_symboltype_fmt(fmt.ALP1,  symbol_to_seg_format_depth_sets, seginfo_aLP1, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ALP2,  symbol_to_seg_format_depth_sets, seginfo_aLP2, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ALPL,  symbol_to_seg_format_depth_sets, seginfo_aLPL, refpos, symboltype, refsymbol);
    
    filla_symboltype_fmt(fmt.ARP1,  symbol_to_seg_format_depth_sets, seginfo_aRP1, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ARP2,  symbol_to_seg_format_depth_sets, seginfo_aRP2, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ARPL,  symbol_to_seg_format_depth_sets, seginfo_aRPL, refpos, symboltype, refsymbol);
    
    filla_symboltype_fmt(fmt.ALB1,  symbol_to_seg_format_depth_sets, seginfo_aLB1, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ALB2,  symbol_to_seg_format_depth_sets, seginfo_aLB2, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ALBL,  symbol_to_seg_format_depth_sets, seginfo_aLBL, refpos, symboltype, refsymbol);
    
    filla_symboltype_fmt(fmt.ARB1,  symbol_to_seg_format_depth_sets, seginfo_aRB1, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ARB2,  symbol_to_seg_format_depth_sets, seginfo_aRB2, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ARBL,  symbol_to_seg_format_depth_sets, seginfo_aRBL, refpos, symboltype, refsymbol);
    
    filla_symboltype_fmt(fmt.ABQ2,  symbol_to_seg_format_depth_sets, seginfo_aBQ2, refpos, symboltype, refsymbol);
    
    filla_symboltype_fmt(fmt.APF1,  symbol_to_seg_format_depth_sets, seginfo_aPF1, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.APF2,  symbol_to_seg_format_depth_sets, seginfo_aPF2, refpos, symboltype, refsymbol);
    
    filla_symboltype_fmt(fmt.ALI1,  symbol_to_seg_format_depth_sets, seginfo_aLI1, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ALI2,  symbol_to_seg_format_depth_sets, seginfo_aLI2, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ARIf,  symbol_to_seg_format_depth_sets, seginfo_aRIf, refpos, symboltype, refsymbol);

    filla_symboltype_fmt(fmt.ARI1,  symbol_to_seg_format_depth_sets, seginfo_aRI1, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ARI2,  symbol_to_seg_format_depth_sets, seginfo_aRI2, refpos, symboltype, refsymbol);
    filla_symboltype_fmt(fmt.ALIr,  symbol_to_seg_format_depth_sets, seginfo_aLIr, refpos, symboltype, refsymbol);
    
    const auto & symbol_to_frag_format_depth_sets = symbol2CountCoverageSet12.symbol_to_frag_format_depth_sets;
    
    fill_symboltype_fmt(fmt.BDPf,  symbol_to_frag_format_depth_sets[0], FRAG_bDP, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.BTAf,  symbol_to_frag_format_depth_sets[0], FRAG_bTA, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.BTBf,  symbol_to_frag_format_depth_sets[0], FRAG_bTB, refpos, symboltype, refsymbol);
    
    fill_symboltype_fmt(fmt.BDPr,  symbol_to_frag_format_depth_sets[1], FRAG_bDP, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.BTAr,  symbol_to_frag_format_depth_sets[1], FRAG_bTA, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.BTBr,  symbol_to_frag_format_depth_sets[1], FRAG_bTB, refpos, symboltype, refsymbol);
    
    const auto & symbol_to_fam_format_depth_sets = symbol2CountCoverageSet12.symbol_to_fam_format_depth_sets_2strand;
    
    fill_symboltype_fmt(fmt.CDP1f, symbol_to_fam_format_depth_sets[0], FAM_cDP1, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDP12f,symbol_to_fam_format_depth_sets[0], FAM_cDP12,refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDP2f, symbol_to_fam_format_depth_sets[0], FAM_cDP2, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDP3f, symbol_to_fam_format_depth_sets[0], FAM_cDP3, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.C1DPf, symbol_to_fam_format_depth_sets[0], FAM_c1DP, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDPMf, symbol_to_fam_format_depth_sets[0], FAM_cDPM, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDPmf, symbol_to_fam_format_depth_sets[0], FAM_cDPm, refpos, symboltype, refsymbol);
    
    fill_symboltype_fmt(fmt.CDP1r, symbol_to_fam_format_depth_sets[1], FAM_cDP1, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDP12r,symbol_to_fam_format_depth_sets[1], FAM_cDP12,refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDP2r, symbol_to_fam_format_depth_sets[1], FAM_cDP2, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDP3r, symbol_to_fam_format_depth_sets[1], FAM_cDP3, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.C1DPr, symbol_to_fam_format_depth_sets[1], FAM_c1DP, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDPMr, symbol_to_fam_format_depth_sets[1], FAM_cDPM, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.CDPmr, symbol_to_fam_format_depth_sets[1], FAM_cDPm, refpos, symboltype, refsymbol);

    const auto & symbol_to_duplex_format_depth_sets = symbol2CountCoverageSet12.symbol_to_duplex_format_depth_sets;
    fill_symboltype_fmt(fmt.DDP1,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP1, refpos, symboltype, refsymbol);
    fill_symboltype_fmt(fmt.DDP2,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP2, refpos, symboltype, refsymbol);
    
    resetBcfFormatD(fmt);
    return {{ fmt.BDPf[0] + fmt.BDPr[0], MAX(fmt.CDP1f[0], fmt.CDP12f[0]) + MAX(fmt.CDP1r[0], fmt.CDP12r[0]) }};
};

const auto 
vectorsum(const auto v1, const auto v2) {
    assert(v1.size() == v2.size());
    auto ret = v1;
    for (size_t i = 0; i < v1.size(); i++) {
        ret[i] += v2[i];
    }
    return ret;
}

std::array<uvc1_readnum_t, 2>
BcfFormat_symbol_init(
        bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet12, 
        uvc1_refgpos_t refpos, 
        const AlignmentSymbol symbol,
        const auto & mutform2count4vec_bq,
        const auto & indices_bq,
        const auto & mutform2count4vec_fq,
        const auto & indices_fq,
        const uvc1_readnum_t bDPa,
        const uvc1_readnum_t cDP0a,
        const std::string & gapSa,
        uvc1_qual_t minABQ,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    const int a = 0;
    
    const auto & symbol_to_VQ_format_tag_sets = symbol2CountCoverageSet12.symbol_to_VQ_format_tag_sets;
    const auto & symbol_to_seg_format_depth_sets = symbol2CountCoverageSet12.symbol_to_seg_format_depth_sets;
    
    fill_symbol_fmt(fmt.a1BQf, symbol_to_VQ_format_tag_sets,    VQ_a1BQf, refpos, symbol, a);
    fill_symbol_fmt(fmt.a1BQr, symbol_to_VQ_format_tag_sets,    VQ_a1BQr, refpos, symbol, a);
    filla_symbol_fmt(fmt.aMQs,  symbol_to_seg_format_depth_sets, seginfo_aMQs,  refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aP1,  symbol_to_seg_format_depth_sets, seginfo_aP1,  refpos, symbol, a);
    filla_symbol_fmt(fmt.aP2,  symbol_to_seg_format_depth_sets, seginfo_aP2,  refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aDPff, symbol_to_seg_format_depth_sets, seginfo_aDPff, refpos, symbol, a);
    filla_symbol_fmt(fmt.aDPfr, symbol_to_seg_format_depth_sets, seginfo_aDPfr, refpos, symbol, a);
    filla_symbol_fmt(fmt.aDPrf, symbol_to_seg_format_depth_sets, seginfo_aDPrf, refpos, symbol, a);
    filla_symbol_fmt(fmt.aDPrr, symbol_to_seg_format_depth_sets, seginfo_aDPrr, refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aLP1, symbol_to_seg_format_depth_sets, seginfo_aLP1, refpos, symbol, a);
    filla_symbol_fmt(fmt.aLP2, symbol_to_seg_format_depth_sets, seginfo_aLP2, refpos, symbol, a);
    filla_symbol_fmt(fmt.aLPL, symbol_to_seg_format_depth_sets, seginfo_aLPL, refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aRP1, symbol_to_seg_format_depth_sets, seginfo_aRP1, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRP2, symbol_to_seg_format_depth_sets, seginfo_aRP2, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRPL, symbol_to_seg_format_depth_sets, seginfo_aRPL, refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aLB1, symbol_to_seg_format_depth_sets, seginfo_aLB1, refpos, symbol, a);
    filla_symbol_fmt(fmt.aLB2, symbol_to_seg_format_depth_sets, seginfo_aLB2, refpos, symbol, a);
    filla_symbol_fmt(fmt.aLBL, symbol_to_seg_format_depth_sets, seginfo_aLBL, refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aRB1, symbol_to_seg_format_depth_sets, seginfo_aRB1, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRB2, symbol_to_seg_format_depth_sets, seginfo_aRB2, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRBL, symbol_to_seg_format_depth_sets, seginfo_aRBL, refpos, symbol, a);
    
    // extra
    filla_symbol_fmt(fmt.a2XM2,symbol_to_seg_format_depth_sets, seginfo_a2XM2,refpos, symbol, a);
    filla_symbol_fmt(fmt.a2BM2,symbol_to_seg_format_depth_sets, seginfo_a2BM2,refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aBQ2, symbol_to_seg_format_depth_sets, seginfo_aBQ2, refpos, symbol, a);

    filla_symbol_fmt(fmt.aPF1, symbol_to_seg_format_depth_sets, seginfo_aPF1, refpos, symbol, a);
    filla_symbol_fmt(fmt.aPF2, symbol_to_seg_format_depth_sets, seginfo_aPF2, refpos, symbol, a);

    filla_symbol_fmt(fmt.aLI1, symbol_to_seg_format_depth_sets, seginfo_aLI1, refpos, symbol, a);
    filla_symbol_fmt(fmt.aLI2, symbol_to_seg_format_depth_sets, seginfo_aLI2, refpos, symbol, a);
    filla_symbol_fmt(fmt.aLIr, symbol_to_seg_format_depth_sets, seginfo_aLIr, refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aRI1, symbol_to_seg_format_depth_sets, seginfo_aRI1, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRI2, symbol_to_seg_format_depth_sets, seginfo_aRI2, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRIf, symbol_to_seg_format_depth_sets, seginfo_aRIf, refpos, symbol, a);
    
    const auto & symbol_to_frag_format_depth_sets = symbol2CountCoverageSet12.symbol_to_frag_format_depth_sets;
    
    fill_symbol_fmt(fmt.bDPf, symbol_to_frag_format_depth_sets[0], FRAG_bDP, refpos, symbol, a);
    fill_symbol_fmt(fmt.bTAf, symbol_to_frag_format_depth_sets[0], FRAG_bTA, refpos, symbol, a);
    fill_symbol_fmt(fmt.bTBf, symbol_to_frag_format_depth_sets[0], FRAG_bTB, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.bDPr, symbol_to_frag_format_depth_sets[1], FRAG_bDP, refpos, symbol, a);
    fill_symbol_fmt(fmt.bTAr, symbol_to_frag_format_depth_sets[1], FRAG_bTA, refpos, symbol, a);
    fill_symbol_fmt(fmt.bTBr, symbol_to_frag_format_depth_sets[1], FRAG_bTB, refpos, symbol, a);
    
    const auto & symbol_to_fam_format_depth_sets = symbol2CountCoverageSet12.symbol_to_fam_format_depth_sets_2strand;
    fill_symbol_fmt(fmt.cDP1f, symbol_to_fam_format_depth_sets[0], FAM_cDP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP12f,symbol_to_fam_format_depth_sets[0], FAM_cDP12,refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP2f, symbol_to_fam_format_depth_sets[0], FAM_cDP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP3f, symbol_to_fam_format_depth_sets[0], FAM_cDP3, refpos, symbol, a);
    fill_symbol_fmt(fmt.c1DPf, symbol_to_fam_format_depth_sets[0], FAM_c1DP, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPMf, symbol_to_fam_format_depth_sets[0], FAM_cDPM, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPmf, symbol_to_fam_format_depth_sets[0], FAM_cDPm, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.cDP1r, symbol_to_fam_format_depth_sets[1], FAM_cDP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP12r,symbol_to_fam_format_depth_sets[1], FAM_cDP12,refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP2r, symbol_to_fam_format_depth_sets[1], FAM_cDP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP3r, symbol_to_fam_format_depth_sets[1], FAM_cDP3, refpos, symbol, a);
    fill_symbol_fmt(fmt.c1DPr, symbol_to_fam_format_depth_sets[1], FAM_c1DP, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPMr, symbol_to_fam_format_depth_sets[1], FAM_cDPM, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPmr, symbol_to_fam_format_depth_sets[1], FAM_cDPm, refpos, symbol, a);
    
    const auto & symbol_to_duplex_format_depth_sets = symbol2CountCoverageSet12.symbol_to_duplex_format_depth_sets;
    fill_symbol_fmt(fmt.dDP1,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.dDP2,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP2, refpos, symbol, a);
    
    filla_symbol_fmt(fmt.aLPT, symbol_to_seg_format_depth_sets, seginfo_aLPT, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRPT, symbol_to_seg_format_depth_sets, seginfo_aRPT, refpos, symbol, a);
    filla_symbol_fmt(fmt.aLIT, symbol_to_seg_format_depth_sets, seginfo_aLIT, refpos, symbol, a);
    filla_symbol_fmt(fmt.aRIT, symbol_to_seg_format_depth_sets, seginfo_aRIT, refpos, symbol, a);
    filla_symbol_fmt(fmt.aP3,  symbol_to_seg_format_depth_sets, seginfo_aP3,  refpos, symbol, a);
    
    fmt.DP = fmt.CDP1f[0] + fmt.CDP1r[0];
    fmt.AD = vectorsum(fmt.cDP1f, fmt.cDP1r);
    fmt.bDP = fmt.BDPf[0] + fmt.BDPr[0];
    fmt.bAD = vectorsum(fmt.bDPf, fmt.bDPr);
    fmt.c2DP = fmt.CDP2f[0] + fmt.CDP2r[0];
    fmt.c2AD = vectorsum(fmt.cDP2f, fmt.cDP2r);
    
    fill_symbol_VQ_fmts(
            fmt,
            symbol_to_VQ_format_tag_sets,
            refpos,
            symbol,
            minABQ,
            paramset,
            a);
    
    fmt.bHap = mutform2count4map_to_phase(mutform2count4vec_bq, indices_bq);
    fmt.cHap = mutform2count4map_to_phase(mutform2count4vec_fq, indices_fq);
    clear_push(fmt.bDPa, bDPa, a);
    clear_push(fmt.cDP0a, cDP0a, a);
    clear_push(fmt.gapSa, gapSa, a);
    
    fmt.note = symbol2CountCoverageSet12.additional_note.getByPos(refpos)[symbol];
    return std::array<uvc1_readnum_t, 2>({{fmt.bDPf[a] + fmt.bDPr[a], fmt.cDP1f[a] + fmt.cDP1r[a]}});
};

double 
calc_normFA_from_rawFA_refbias(double FA, double refbias) {
    return (FA + FA * refbias) / (FA + (1.0 - FA) / (1.0 + refbias) + FA * refbias);
}

int
fmt_bias_push(auto & vecFA, const auto refFA, const auto biasFA, const double thresFAratio, std::string & fts, const std::string & ft) {
    vecFA.push_back(-numstates2deciphred(biasFA));
    if (biasFA * thresFAratio < refFA) {
        if (fts.size() > 0) { fts += std::string("&"); }
        fts += ft;
        return 1;
    }
    return 0;
}

int
BcfFormat_symbol_calc_DPv(
        bcfrec::BcfFormat & fmt,
        const RegionalTandemRepeat & rtr1,
        const RegionalTandemRepeat & rtr2,
        const double tpfa,
        const AlignmentSymbol refsymbol,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    const double unbias_ratio = ((NOT_PROVIDED == paramset.vcf_tumor_fname) ? 1.0 : sqrt(2.0));
    const double unbias_qualadd = ((NOT_PROVIDED == paramset.vcf_tumor_fname) ? 0 : 3);
    const uvc1_qual_t allbias_allprior = ((NOT_PROVIDED == paramset.vcf_tumor_fname) ? 0 : 31);
    
    const bool is_strong_amplicon = (fmt.APDP[5] * 100 > fmt.APDP[0] * 60);
    const bool is_weak_amplicon = (fmt.APDP[5] * 100 > fmt.APDP[0] * 40);
    
    const int a = 0;
    const bool is_rescued = (tpfa >= 0);
    const double pfa = (is_rescued ? tpfa : 0.5);
    
    // Non-UMI universality-based prequal allele fraction
    const uvc1_readnum_t ADP = MAX(fmt.ADPff[0] + fmt.ADPfr[0] + fmt.ADPrf[0] + fmt.ADPrr[0], fmt.APDP[9]);
    const uvc1_readnum_t aDP = (fmt.aDPff[a] + fmt.aDPfr[a] + fmt.aDPrf[a] + fmt.aDPrr[a]);
    
    const auto & f = fmt;
    const auto symbol = AlignmentSymbol(LAST(f.VTI));
    
    double _counterbias_P_FA = 1e-9;
    double _counterbias_BQ_FA = 1e-9;
    double _dir_bias_div = 1.0;
    const bool is_nmore_amplicon = ((NOT_PROVIDED == paramset.vcf_tumor_fname) ? is_strong_amplicon : is_weak_amplicon);
    if ((is_nmore_amplicon && (0x2 == (0x2 & paramset.nobias_flag))) || ((!is_nmore_amplicon) && (0x1 == (0x1 & paramset.nobias_flag)))) {
        // counter bias : position 23-10 bp and base quality 13
        // not-pos-bias < pos-bias
        const double using_bias_aFA = (aDP - LAST(fmt.aP1) + 0.5) / (ADP - fmt.AP1[0] + 1.0);
        const double using_all_aFA = (aDP + 0.5) / (ADP + 1.0);
        
        // for ERP015684 T610 17_37871708_A_G, is_pos_counterbias becomes true to rescue the TP call by eliminating the primer-induced position bias and orientation bias
        const bool is_pos_counterbias = (
                   (using_bias_aFA * 4  < using_all_aFA * (unbias_ratio - DBL_EPSILON))
                && (LAST(fmt.aP1) * (unbias_ratio - DBL_EPSILON) > aDP - LAST(fmt.aP1))
                && ((ADP - fmt.AP1[0]) * 4 * (unbias_ratio - DBL_EPSILON)  > fmt.AP1[0])
                && ((0 == paramset.primerlen && 0 != paramset.primerlen2) || !isSymbolSubstitution(symbol)));
        if (is_pos_counterbias) {
            UPDATE_MAX(_counterbias_P_FA, (LAST(fmt.aP1) + 0.5) / (MAX(fmt.AP1[0], fmt.APDP[9]) + 1.0));
        } else {
            UPDATE_MAX(_counterbias_P_FA, 2e-9);
        }
        if (isSymbolSubstitution(symbol)) {
            const bool is_f_good_cov = ((fmt.ADPfr[0] + fmt.ADPrr[0]) + 150 <= (fmt.ADPff[0] + fmt.ADPrf[0]) * 5 * unbias_ratio);
            const bool is_r_good_cov = ((fmt.ADPff[0] + fmt.ADPrf[0]) + 150 <= (fmt.ADPfr[0] + fmt.ADPrr[0]) * 5 * unbias_ratio);
            const uvc1_qual_t avg_f_aBQ = (LAST(fmt.a1BQf) / MAX(1, fmt.aDPff[a] + fmt.aDPrf[a]));
            const uvc1_qual_t avg_r_aBQ = (LAST(fmt.a1BQr) / MAX(1, fmt.aDPfr[a] + fmt.aDPrr[a]));
            const uvc1_qual_t avg_f_ABQ = (fmt.A1BQf[0] / MAX(1, fmt.ADPff[0] + fmt.ADPrf[0]));
            const uvc1_qual_t avg_r_ABQ = (fmt.A1BQr[0] / MAX(1, fmt.ADPfr[0] + fmt.ADPrr[0]));
            
            const bool is_f_BQ_counterbias = (
                    (LAST(fmt.a1BQf) >= LAST(fmt.a1BQr)) 
                    && (is_f_good_cov && is_r_good_cov)
                    && (avg_f_aBQ + unbias_qualadd >= avg_r_ABQ + 14) 
                    && (avg_r_ABQ <= 14 + unbias_qualadd));
            if (is_f_BQ_counterbias) {
                UPDATE_MAX(_counterbias_BQ_FA, (fmt.aDPff[a] + fmt.aDPrf[a] + 0.5) / (fmt.ADPff[0] + fmt.ADPrf[0] + 1.0));
            }
            const bool is_r_BQ_counterbias = (
                    (LAST(fmt.a1BQr) >= LAST(fmt.a1BQf))
                    && (is_f_good_cov && is_r_good_cov)
                    && (avg_r_aBQ + unbias_qualadd >= avg_f_ABQ + 14)
                    && (avg_f_ABQ <= 14 + unbias_qualadd));
            if (is_r_BQ_counterbias) {
                UPDATE_MAX(_counterbias_BQ_FA, (fmt.aDPfr[a] + fmt.aDPrr[a] + 0.5) / (fmt.ADPfr[0] + fmt.ADPrr[0] + 1.0));
            }
        } else {
            _dir_bias_div = (1.0 + LAST(fmt.gapSa).size() / paramset.indel_str_repeatsize_max); 
        }
    }
    const double counterbias_P_FA = _counterbias_P_FA;
    const double counterbias_BQ_FA = _counterbias_BQ_FA;
    const double dir_bias_div = _dir_bias_div;
    
    // const double aDPdel = (fmt.ADPff[1] + fmt.ADPfr[1] + fmt.ADPrf[1] + fmt.ADPrr[1]); // APDP
    const auto aDPgap = non_neg_minus(MAX(fmt.APDP[1], fmt.APDP[2]), fmt.aP3[a]);
    const double aDPFAgap = ((rtr1.tracklen + rtr2.tracklen < paramset.indel_str_repeatsize_max) ? 1.0 : ((fmt.aP3[a] + pfa) / (aDPgap + 1.0)));
    const double aDPFA1 = ((aDP + pfa) / (ADP + 1.0));
    const double labelFA = (fmt.aP2[a] + 1.5 + fmt.aP2[a]) / (fmt.AP2[0] + 2.0 + fmt.aP2[a]); // assay-type bias
    const double aDPFA = MIN(
            (isSymbolSubstitution(symbol) ? MIN(aDPFA1, MAX(aDPFA1 / 3, aDPFAgap)) : (aDPFA1)), // the number 3 is magic
            labelFA * (ADP + 1.0) / (fmt.AP2[0] + 0.5) * unbias_ratio);
    // substitution in indel region, substitution in indel-prone region or indel, other cases 
    uvc1_readnum_t aDPplus = (isSymbolSubstitution(symbol) ? 0 : ((aDP + 1) * paramset.bias_prior_DPadd_perc / 100));
    double dp_coef = ((symbol == LINK_M) ? MAX(paramset.contam_any_mul_frac, 1.0 - MAX(rtr1.tracklen, rtr2.tracklen) / (MAX3(1, f.ALPL[0], f.ARPL[0]) / MAX(1.0/150.0, f.ABQ2[0]))) : 1.0);
    double _aPpriorfreq = paramset.bias_priorfreq_pos; // probability that there is one indel at a pos
    double _aBpriorfreq = _aPpriorfreq;
    const bool is_in_indel_read = (   (f.APXM[1]) / 15.0 * paramset.microadjust_bias_pos_indel_fold * (paramset.bias_prior_var_DP_mul) > (aDP + aDPplus) * dp_coef);
    const bool is_in_indel_len  = (MAX(f.APDP[1],  f.APDP[2]) * (paramset.bias_prior_var_DP_mul) > (aDP + aDPplus) * dp_coef);
    const bool is_in_indel_rtr  = (MAX(f.APDP[3],  f.APDP[4]) * (paramset.bias_prior_var_DP_mul) > (aDP + aDPplus) * dp_coef);
    const bool is_in_rtr = (MAX(rtr1.tracklen, rtr2.tracklen) > round(paramset.indel_polymerase_size));
    
    const bool is_in_dnv_read = ((SEQUENCING_PLATFORM_IONTORRENT == paramset.inferred_sequencing_platform) && (f.APDP[7] * 2 > f.APDP[6]));
    
    // is_outlier_del is set to false to prevent potential false negative calls on Illumina, it may be useful for some rare corner-cases in BGI data.
    const bool is_outlier_del = false ; 
    // (is_in_indel_read && (isSymbolDel(symbol)) && (fmt.APXM[3] * 2 > UNSIGN2SIGN(3 * LAST(fmt.gapSa).size())) && (MIN(fmt.aLPT[a], fmt.aRPT[a]) < aDP * 20)); 
    // read level
    if (is_in_indel_read || is_in_dnv_read || 
            (((isSymbolIns(symbol) || isSymbolDel(symbol)) 
                && (fmt.APXM[0] > fmt.APXM[1] * paramset.microadjust_bias_pos_indel_misma_to_indel_ratio)))) {
        _aPpriorfreq -= paramset.bias_priorfreq_indel_in_read_div;
        _aBpriorfreq -= paramset.bias_priorfreq_indel_in_read_div;
    }
    // regional level
    if (LINK_M != symbol && LINK_NN != symbol) { 
        double maxpf = 0;
        if (is_in_indel_len) { UPDATE_MAX(maxpf, paramset.bias_priorfreq_indel_in_var_div2); }
        if (is_in_indel_rtr) { UPDATE_MAX(maxpf, paramset.bias_priorfreq_indel_in_str_div2); }
        if (is_in_rtr)       { UPDATE_MAX(maxpf, paramset.bias_priorfreq_var_in_str_div2); }
        if (is_outlier_del) {
            UPDATE_MAX(maxpf, 10);
        }
        _aBpriorfreq -= maxpf;
        _aPpriorfreq -= maxpf;
    }
    
    const double aPpriorfreq = _aPpriorfreq + allbias_allprior;
    const double aBpriorfreq = _aBpriorfreq + allbias_allprior;
    fmt.nPF = {{ (uvc1_qual_t)round(aPpriorfreq), (uvc1_qual_t)round(aBpriorfreq) }};
    
    const double aIpriorfreq = (isSymbolSubstitution(symbol) ? paramset.bias_priorfreq_ipos_snv : paramset.bias_priorfreq_ipos_indel) 
            + allbias_allprior;
    const auto homopol_len = ((1 == rtr1.unitlen) ? rtr1.tracklen : 0) + ((1 == rtr2.unitlen) ? rtr2.tracklen : 0);
    const double aSBpriorfreq = (isSymbolSubstitution(symbol) 
            ? (MIN(
                // The IonTorrent platform is very susceptible to strand bias in homopolymer regions
                non_neg_minus(fmt.aBQ[a], (
                        ((SEQUENCING_PLATFORM_IONTORRENT == paramset.inferred_sequencing_platform) 
                        && (homopol_len > 0)
                        && (isSymbolSubstitution(symbol) || (LINK_D1 == symbol) || (LINK_I1 == symbol)))
                    ? MIN(5 * homopol_len, 20) : 0)),
                fmt.bMQ[a])
              + paramset.bias_priorfreq_strand_snv_base)
            : (paramset.bias_priorfreq_strand_indel))
            + allbias_allprior;
    
    const auto aLPFAx2 = dp4_to_pcFA<false>(f.aLP1[a], aDP, f.ALP2[0] + f.aLP1[a] - f.aLP2[a], ADP, paramset.powlaw_exponent, phred2nat(aPpriorfreq),
            MAX(1, f.aLPL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ALPL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? paramset.bias_FA_pseudocount_indel_in_read : 0.5)); 
            // Here we should not use f.aBQ2[a] instead of aDP because biased read count is not in aBQ2.
    const auto aRPFAx2 = dp4_to_pcFA<false>(f.aRP1[a], aDP, f.ARP2[0] + f.aRP1[a] - f.aRP2[a], ADP, paramset.powlaw_exponent, phred2nat(aPpriorfreq),
            MAX(1, f.aRPL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ARPL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? paramset.bias_FA_pseudocount_indel_in_read : 0.5));
    double aLPFA = aLPFAx2[0];
    double aRPFA = aRPFAx2[0];
    
    const auto aLBFAx2 = dp4_to_pcFA<false>(f.aLB1[a], aDP, f.ALB2[0] + f.aLB1[a] - f.aLB2[a], ADP, paramset.powlaw_exponent, phred2nat(aBpriorfreq),
            MAX(1, f.aLBL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ALBL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? paramset.bias_FA_pseudocount_indel_in_read : 0.5));
    const auto aRBFAx2 = dp4_to_pcFA<false>(f.aRB1[a], aDP, f.ARB2[0] + f.aRB1[a] - f.aRB2[a], ADP, paramset.powlaw_exponent, phred2nat(aBpriorfreq),
            MAX(1, f.aRBL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ARBL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? paramset.bias_FA_pseudocount_indel_in_read : 0.5));
    double aLBFA = aLBFAx2[0];
    double aRBFA = aRBFAx2[0];
    const bool is_tmore_amplicon = ((NOT_PROVIDED == paramset.vcf_tumor_fname) ? is_weak_amplicon : is_strong_amplicon);
    
    std::array<double, 2> _aLIFAx2 = {{ 0, 0 }};
    {
        double ALpd = (f.ALI2[0] +                 0.5) / (f.ADPfr[0] + f.ADPrr[0] - f.ALI2[0] + 0.5);
        double aLpd = (f.aLI1[a] + ALpd / (1.0 + ALpd)) / (f.aDPfr[a] + f.aDPrr[a] - f.aLI1[a] + 1.0 / (1.0 + ALpd)); 
        _aLIFAx2 = dp4_to_pcFA<false>(
                (f.aLI1[a]), (f.aDPfr[a] + f.aDPrr[a]), 
                (f.ALI2[0] + f.aLI1[a] - f.aLI2[a]), (f.ADPfr[0] + f.ADPrr[0]), 
                paramset.powlaw_exponent, phred2nat(aIpriorfreq), 
                aLpd, ALpd, 0.25, 0.5);
        
    }
    double aLIFA = _aLIFAx2[0] * ((is_tmore_amplicon) ? (dir_bias_div) : MAX(dir_bias_div, aDPFA / _aLIFAx2[1]));
    
    std::array<double, 2> _aRIFAx2 = {{ 0, 0 }};
    {
        double ARpd = (f.ARI2[0] +                 0.5) / (f.ADPff[0] + f.ADPrf[0] - f.ARI2[0] + 0.5);
        double aRpd = (f.aRI1[a] + ARpd / (1.0 + ARpd)) / (f.aDPff[a] + f.aDPrf[a] - f.aRI1[a] + 1.0 / (1.0 + ARpd));
        _aRIFAx2 = dp4_to_pcFA<false>(
                (f.aRI1[a]), (f.aDPff[a] + f.aDPrf[a]), 
                (f.ARI2[0] + f.aRI1[a] - f.aRI2[a]), (f.ADPff[0] + f.ADPrf[0]), 
                paramset.powlaw_exponent, phred2nat(aIpriorfreq),
                aRpd, ARpd, 0.25, 0.5);
    }
    double aRIFA = _aRIFAx2[0] * ((is_tmore_amplicon) ? (dir_bias_div) : MAX(dir_bias_div, aDPFA / _aRIFAx2[1]));
    
    const double aSIFA = MAX( 
        (f.aLI1[a] + 0.5) / (f.ALI2[0] + f.aLI1[a] - f.aLI2[a] + 1.0), 
        (f.aRI1[a] + 0.5) / (f.ARI2[0] + f.aRI1[a] - f.aRI2[a] + 1.0));
    if (isSymbolIns(symbol) || isSymbolDel(symbol)) {
        const auto & indelstring = LAST(fmt.gapSa);
        const double indel_multialleles_coef = MAX(1, fmt.bDPa[a]) / (double)MAX(1, fmt.bDPf[a] + fmt.bDPr[a]);
        const bool is_in_indel_major_reg = ((MAX(f.APDP[1], f.APDP[3]) + MAX(f.APDP[2], f.APDP[4])) * 0.5 * (1.0+(double)FLT_EPSILON) < aDP * indel_multialleles_coef);
        // The 16 is from germ-HG002_het_fp_0-1_13_nochr1_173328669_A_ATTCAAGGACTTTCTTTTTACCAGCTGT
        if (    (MIN(UNSIGN2SIGN(indelstring.size()), paramset.microadjust_nobias_pos_indel_maxlen) * aDPFA * indel_multialleles_coef 
                    >= paramset.nobias_pos_indel_lenfrac_thres) || 
                (MAX(rtr1.tracklen, rtr2.tracklen) >= paramset.nobias_pos_indel_str_track_len
                    && is_in_indel_major_reg 
                    && (!is_outlier_del)
                    && !(fmt.APXM[0] > fmt.APXM[1] * paramset.microadjust_nobias_pos_indel_misma_to_indel_ratio))) {
            aLPFA += 2.0;
            aRPFA += 2.0;
            aLBFA += 2.0;
            aRBFA += 2.0;
        }
        if (LAST(fmt.bMQ) >= paramset.microadjust_nobias_pos_indel_bMQ && LAST(fmt.a2XM2) * 100 >= aDP * 100 * paramset.microadjust_nobias_pos_indel_perc) { 
            aLIFA += 2.0; 
            aRIFA += 2.0; 
        }
    } else if (LINK_M == symbol || LINK_NN == symbol) {
        double pc = paramset.bias_FA_pseudocount_indel_in_read; 
        aLBFA = (double)MIN(aLBFA, (pc + LAST(fmt.aLB1)) / (double)(pc * 2 + ADP));
        aRBFA = (double)MIN(aRBFA, (pc + LAST(fmt.aRB1)) / (double)(pc * 2 + ADP));
    } else if (refsymbol == symbol) {
        aLIFA = aRIFA = MAX(aLIFA, aRIFA); // reference error or long indel on either the left or right frag side does not affect the ref SNP allele.
    }
    
    const auto avg_sqr_indel_len = MAX(fmt.APXM[4] / MAX(1, fmt.APDP[1]), fmt.APXM[5] / MAX(1, fmt.APDP[2]));
    if ((!isSymbolSubstitution(symbol)) 
            && (mathsquare(paramset.microadjust_nobias_pos_indel_maxlen) < avg_sqr_indel_len)
            && (LINK_M == symbol || LINK_NN == symbol || (UNSIGN2SIGN(mathsquare(LAST(fmt.gapSa).size() * 2)) < avg_sqr_indel_len))) {
        double pc = paramset.bias_FA_pseudocount_indel_in_read; 
        aLPFA = (double)MIN(aLPFA, (pc + LAST(fmt.aLP1)) / (double)(pc * 2 + fmt.ALP1[0]));
        aRPFA = (double)MIN(aRPFA, (pc + LAST(fmt.aRP1)) / (double)(pc * 2 + fmt.ALP1[0]));
    }
    
    if (NOT_PROVIDED != paramset.vcf_tumor_fname) {
        aLIFA = aRIFA = MAX(aLIFA, aRIFA);
    }
    
    const double aPFFA = (fmt.aPF1[a] + pfa * 100.0) / (fmt.APF2[0] + (fmt.aPF1[a] - fmt.aPF2[a]) + 100.0);

    const auto aSSFAx2 = dp4_to_pcFA<true>(
            fmt.aRIf[a], fmt.aLIr[a], fmt.ARIf[0], fmt.ALIr[0],
            paramset.powlaw_exponent, phred2nat(aSBpriorfreq));
    const auto bias_priorfreq_orientation_base = (isSymbolSubstitution(symbol) 
            ? paramset.bias_priorfreq_orientation_snv_base : paramset.bias_priorfreq_orientation_indel_base) + allbias_allprior;
    
    auto _cROFA1x2 = dp4_to_pcFA<true>(fmt.cDP1f[a], fmt.cDP1r[a], fmt.CDP1f[0], fmt.CDP1r[0], paramset.powlaw_exponent, 
                log(mathsquare(aDPFA)) + phred2nat(bias_priorfreq_orientation_base));
    if (paramset.bias_is_orientation_artifact_mixed_with_sequencing_error) {
        // this is useful for dealing with heavy FFPE artifact and heavy sequencing error
        const auto cROFA10x2 = dp4_to_pcFA<true>(fmt.cDP1f[a], fmt.cDP1r[a], fmt.CDP1f[0], fmt.CDP1r[0], paramset.powlaw_exponent, 
                log(mathsquare(aDPFA)) + phred2nat(bias_priorfreq_orientation_base));
        const auto cROFA12x2 = dp4_to_pcFA<true>(fmt.cDP12f[a], fmt.cDP12r[a], fmt.CDP12f[0], fmt.CDP12r[0], paramset.powlaw_exponent, 
                log(mathsquare(aDPFA)) + phred2nat(bias_priorfreq_orientation_base));
        _cROFA1x2 = (
                    ((fmt.ADPff[0] * 8 >= ADP) && (fmt.ADPfr[0] * 8 >= ADP) 
                  && (fmt.ADPrf[0] * 8 >= ADP) && (fmt.ADPrr[0] * 8 >= ADP)) 
                ? cROFA12x2 : cROFA10x2);
    }
    const auto cROFA1x2 = _cROFA1x2; 
    const auto cROFA2x2 = dp4_to_pcFA<true>(fmt.cDP2f[a], fmt.cDP2r[a], fmt.CDP2f[0], fmt.CDP2r[0], paramset.powlaw_exponent, 
            log(mathsquare(aDPFA)) + phred2nat(bias_priorfreq_orientation_base));
    
    double aSSFA = aSSFAx2[0] * dir_bias_div;
    double cROFA1 = cROFA1x2[0] * dir_bias_div;
    double cROFA2 = cROFA2x2[0] * dir_bias_div;
    
    double bFA = (fmt.bDPa[a] + pfa) / (fmt.BDPf[0] + fmt.BDPr[0] + 1.0);
    double cFA0 = (fmt.cDP0a[a] + pfa * (does_fmt_imply_short_frag(fmt) ? 0.1 : 1.0)) / (fmt.CDP1f[0] + fmt.CDP1r[0] + 1.0);

    const bool is_strand_r_weak = ((f.ADPfr[0] + f.ADPrr[0]) * paramset.microadjust_nobias_strand_all_fold < (f.ADPff[0] + f.ADPrf[0]) * unbias_ratio);
    const bool is_strand_f_weak = ((f.ADPff[0] + f.ADPrf[0]) * paramset.microadjust_nobias_strand_all_fold < (f.ADPfr[0] + f.ADPrr[0]) * unbias_ratio);
    if (is_strand_r_weak) {
        aLIFA += 4.0;
        aSSFA += 4.0;
    }
    if (is_strand_f_weak) {
        aRIFA += 4.0;
        aSSFA += 4.0;
    }
    
    const double aLPFA2 = MAX(aDPFA * 0.01, aLPFA);
    const double aRPFA2 = MAX(aDPFA * 0.01, aRPFA);
    const double aLBFA2 = MAX(aDPFA * 0.01, aLBFA);
    const double aRBFA2 = MAX(aDPFA * 0.01, aRBFA);
    const double aLIFA2 = MAX(aDPFA * 0.01, aLIFA);
    const double aRIFA2 = MAX(aDPFA * 0.01, aRIFA);
    const double aSSFA2 = MAX(aDPFA * 0.05, aSSFA);
    
    cROFA1 = MAX(aDPFA * 1e-4, cROFA1);
    cROFA2 = MAX(aDPFA * 1e-4, cROFA2);
    
    // compute systematic error here
    
    const auto fBTA = (double)(fmt.BTAf[0] + fmt.BTAr[0] + 200);
    const auto fBTB = (double)(fmt.BTBf[0] + fmt.BTBr[0] + 6);
    const auto fbTA = (double)(fmt.bTAf[a] + fmt.bTAr[a] + 100);
    const auto fbTB = (double)(fmt.bTBf[a] + fmt.bTBr[a] + 3);
    const double frag_sidelen_frac = 1.0 - MIN(
        BETWEEN(fmt.aLIT[a] / MAX(1, fmt.aDPfr[a] + fmt.aDPrr[a]) - paramset.microadjust_longfrag_sidelength_min, 0, paramset.microadjust_longfrag_sidelength_max), 
        BETWEEN(fmt.aRIT[a] / MAX(1, fmt.aDPff[a] + fmt.aDPrf[a]) - paramset.microadjust_longfrag_sidelength_min, 0, paramset.microadjust_longfrag_sidelength_max))
        / paramset.microadjust_longfrag_sidelength_zeroMQpenalty;
    
    const double _alt_frac_mut_affected_tpos = fbTB / fbTA; // is low by default
    const double alt_frac_mut_affected_tpos = (is_nmore_amplicon ? (MAX(0, _alt_frac_mut_affected_tpos - 0.2) * 1.25) : _alt_frac_mut_affected_tpos);
    const double nonalt_frac_mut_affected_tpos = (fBTB + paramset.contam_any_mul_frac * fbTB - fbTB) / (fBTA + paramset.contam_any_mul_frac * fbTA - fbTA); // is same as alt by default
    const double frac_mut_affected_pos = MAX(paramset.syserr_MQ_NMR_expfrac, 
              paramset.syserr_MQ_NMR_altfrac_coef    * alt_frac_mut_affected_tpos * frag_sidelen_frac
            - paramset.syserr_MQ_NMR_nonaltfrac_coef * nonalt_frac_mut_affected_tpos);
    const uvc1_qual_t bNMQ = round(numstates2phred(pow(frac_mut_affected_pos / paramset.syserr_MQ_NMR_expfrac, (paramset.syserr_MQ_NMR_pl_exponent))) * (frac_mut_affected_pos)); 
    
    clear_push(fmt.bNMa, round(100 * alt_frac_mut_affected_tpos));
    clear_push(fmt.bNMb, round(100 * nonalt_frac_mut_affected_tpos));
    clear_push(fmt.bNMQ, bNMQ);
    
    // end of computation of systematic error
    
    auto min_aFA_vec = std::vector<double>{{
            aDPFA * BETWEEN(1.0 + aDPFA - alt_frac_mut_affected_tpos, 0.1, 1.0),
            aLPFA2,
            aRPFA2,
            aLBFA2,
            aRBFA2,
            aLIFA2,
            aRIFA2,
            aSSFA2,
            aPFFA * aSSFA2 / MAX(aSSFA2, aSSFAx2[1]),
            MAX(aDPFA * 0.01, aSIFA),
            }};
    
    double cFA2 = (fmt.cDP2f[a] + fmt.cDP2r[a] + pfa) / (fmt.CDP2f[0] + fmt.CDP2r[0] + 1.0);
    // The following code without + 1 pseudocount in both nominator and denominator can result in false negative calls at low allele fraction (it is rare but can happen).
    // double cFA3 = (fmt.cDP3f[a] + fmt.cDP3r[a] + pfa + 1) / (fmt.CDP3f[0] + fmt.CDP3r[0] + 1.0 + 1);
    double cFA3 = (fmt.cDP3f[a] + fmt.cDP3r[a] + pfa) / (fmt.CDP3f[0] + fmt.CDP3r[0] + 1.0);
    
    assert( cFA2 > 0 );
    assert( cFA2 < 1 );
    assert( cFA3 > 0 );
    assert( cFA3 < 1 );
    
    fmt.nNFA.push_back(-numstates2deciphred(counterbias_P_FA));
    fmt.nNFA.push_back(-numstates2deciphred(counterbias_BQ_FA));
    
    fmt.nNFA.push_back(-numstates2deciphred(aDPFA));
    fmt_bias_push(fmt.nAFA,  aDPFA,  aLBFA2, paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::BB1L]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aRBFA2, paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::BB1R]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aLPFA2, paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::PB1L]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aRPFA2, paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::PB1R]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aLIFA2, paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::PB2L]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aRIFA2, paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::PB2R]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aSSFA2, paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::SB1]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aPFFA,  paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::PFB]);
    fmt_bias_push(fmt.nAFA,  aDPFA,  aSIFA,  paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::ASI]);
    
    fmt_bias_push(fmt.nBCFA, bFA,    cFA0,   paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::DB1]);
    fmt_bias_push(fmt.nBCFA, cFA0,   bFA,    paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::DB2]);
    
    fmt.nBCFA.push_back(-numstates2deciphred(cFA2));
    fmt.nBCFA.push_back(-numstates2deciphred(cFA3));
    
    fmt_bias_push(fmt.nBCFA, cFA0,  cROFA1,  paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::ROB1]);
    fmt_bias_push(fmt.nBCFA, cFA2,  cROFA2,  paramset.bias_thres_FTS_FA, fmt.FTS, bcfrec::FILTER_IDS[bcfrec::ROB2]);
    
    if (fmt.FTS.size() == 0) {
        fmt.FTS = "PASS";
    }
    
    // self-rescue of the normal by high allele frac
    const double counterbias_normalgerm_FA = ((NOT_PROVIDED == paramset.vcf_tumor_fname || !does_fmt_imply_short_frag(fmt))
        ? 1e-9 : BETWEEN(aDPFA * aDPFA * 6, aDPFA / 6, aDPFA));
    double min_aFA = MAX4(MINVEC(min_aFA_vec), counterbias_P_FA, counterbias_BQ_FA, counterbias_normalgerm_FA);
    
    double dedup_FA = ((NOT_PROVIDED == paramset.vcf_tumor_fname) ? MIN(bFA, cFA0) : (MAX(bFA, cFA0)));
    double min_bcFA = MAX3(MIN3(bFA, cFA0, cROFA1), counterbias_P_FA, counterbias_BQ_FA);
    
    auto min_cFA23_vec = std::vector<double> {{ cFA2, cFA3 }};
    
    double umi_cFA = MINVEC(min_cFA23_vec);
    
    // TODO: find the theory for the following three lines of code.
    double frac_umi2seg1 = MIN(1.0, umi_cFA / aDPFA);
    double frac_umi2seg2 = MIN(1.0, cFA0 / umi_cFA);
    double frac_umi2seg = (frac_umi2seg1 * frac_umi2seg2);
    
    double refbias = 0;
    if ((isSymbolIns(symbol) || isSymbolDel(symbol)) && is_rescued) {
        const auto & indelstring = LAST(fmt.gapSa);
        const uvc1_readpos_t indel_noinfo_nbases = (UNSIGN2SIGN(indelstring.size()) * (isSymbolIns(symbol) ? 2 : 1) 
                + MAX3(UNSIGN2SIGN(indelstring.size()), rtr1.tracklen, rtr2.anyTR_tracklen));
        refbias = (double)(indel_noinfo_nbases) / ((double)(MIN(f.ALPL[0], f.ARPL[0]) * 2 + indel_noinfo_nbases) / (double)(f.ABQ2[0] + 0.5));
        refbias = MIN(refbias, paramset.microadjust_refbias_indel_max);
    }
    // non-UMI push
    double min_abcFA_v = MAX3(MIN(min_aFA, min_bcFA), counterbias_P_FA, counterbias_BQ_FA);
    clear_push(fmt.cDP1v, (uvc1_readnum100x_t)(calc_normFA_from_rawFA_refbias(min_abcFA_v, refbias) * (fmt.CDP1f[0] +fmt.CDP1r[0]) * 100), a);
    double min_abcFA_w = MAX3(MINVEC(std::vector<double>{{
            aLPFA2, aRPFA2,
            aLBFA2, aRBFA2,
            bFA }}), counterbias_P_FA, counterbias_BQ_FA);
    clear_push(fmt.cDP1w, (uvc1_readnum100x_t)(calc_normFA_from_rawFA_refbias(min_abcFA_w, refbias) * (fmt.CDP1f[0] +fmt.CDP1r[0]) * 100), a);
    double min_abcFA_x = MINVEC(std::vector<double>{{ aPFFA, dedup_FA }});
    if (NOT_PROVIDED != paramset.vcf_tumor_fname) {
        min_abcFA_x = MAX3(min_abcFA_x, counterbias_P_FA, counterbias_BQ_FA);
    }
    clear_push(fmt.cDP1x, 1+(uvc1_readnum100x_t)                                       (min_abcFA_x * (fmt.CDP1f[0] +fmt.CDP1r[0]) * 100), a);
    
    // UMI push
    clear_push(fmt.cDP2v, (uvc1_readnum100x_t)(calc_normFA_from_rawFA_refbias(MIN(min_aFA * frac_umi2seg, cROFA2), refbias) 
            * (fmt.CDP2f[0] + fmt.CDP2r[0]) * 100), a);
    double umi_cFA_w = MINVEC(std::vector<double>{{
            aLPFA2 * frac_umi2seg, aRPFA2 * frac_umi2seg,
            aLBFA2 * frac_umi2seg, aRBFA2 * frac_umi2seg,
            aDPFA * frac_umi2seg,
            cROFA2
            }});
    clear_push(fmt.cDP2w, (uvc1_readnum100x_t)(calc_normFA_from_rawFA_refbias(umi_cFA_w, refbias) 
            * (fmt.CDP2f[0] + fmt.CDP2r[0]) * 100), a);
    double umi_cFA_x = MINVEC(std::vector<double>{{ 
            aPFFA * frac_umi2seg,
            aDPFA * frac_umi2seg,
            cROFA2
            }});
    clear_push(fmt.cDP2x, 1+(uvc1_readnum100x_t)(MIN(umi_cFA_x 
            * (fmt.CDP2f[0] + fmt.CDP2r[0]), fmt.cDP2f[a] +fmt.cDP2r[a]) * 100), a);
    
    return 0;
};

enum ReductionType {
    REDUCTION_DP1v,
    REDUCTION_DP1w,
    REDUCTION_DP1x,
    REDUCTION_DP2v,
    REDUCTION_DP2w,
    REDUCTION_DP2x,
    NUM_REDUCTIONS
};

const std::array<ReductionType, NUM_REDUCTIONS> REDUCTION_DPS = {{
    REDUCTION_DP1v,
    REDUCTION_DP1w,
    REDUCTION_DP1x,
    REDUCTION_DP2v,
    REDUCTION_DP2w,
    REDUCTION_DP2x
}};

auto 
getptr_cCDPxvV2(bcfrec::BcfFormat & f, const ReductionType fidx) {
    if (REDUCTION_DP1v == fidx) {
        return std::make_tuple(&(f.CDP1v), &(f.cDP1v));
    }
    if (REDUCTION_DP1w == fidx) {
        return std::make_tuple(&(f.CDP1w), &(f.cDP1w));
    }
    if (REDUCTION_DP1x == fidx) {
        return std::make_tuple(&(f.CDP1x), &(f.cDP1x));
    }
    if (REDUCTION_DP2v == fidx) {
        return std::make_tuple(&(f.CDP2v), &(f.cDP2v));
    }
    if (REDUCTION_DP2w == fidx) {
        return std::make_tuple(&(f.CDP2w), &(f.cDP2w));
    }
    if (REDUCTION_DP2x == fidx) {
        return std::make_tuple(&(f.CDP2x), &(f.cDP2x));
    }
    abort();
}

int
BcfFormat_symbol_sum_DPv(auto & fmts) {
    std::array<uvc1_readnum_t, NUM_REDUCTIONS> cDPxx1s = {{0}};
    std::array<uvc1_readnum_t, NUM_REDUCTIONS> cDPxx2s = {{0}};
    for (auto & fmt : fmts) {
        for (const auto i : REDUCTION_DPS) { cDPxx1s[i] += LAST(*std::get<1>(getptr_cCDPxvV2(std::get<0>(fmt), i))); }
        if (BASE_NN == LAST(std::get<0>(fmt).VTI) || LINK_NN == LAST(std::get<0>(fmt).VTI)) {
            for (const auto i : REDUCTION_DPS) { cDPxx2s[i] = LAST(*std::get<1>(getptr_cCDPxvV2(std::get<0>(fmt), i))); }
        }
    }
    for (auto & fmt : fmts) {
        for (const auto i : REDUCTION_DPS) {
            (*std::get<0>(getptr_cCDPxvV2(std::get<0>(fmt), i)))[0] = cDPxx1s[i];
            (*std::get<0>(getptr_cCDPxvV2(std::get<0>(fmt), i)))[1] = cDPxx2s[i]; 
        }
    }
    return 0;
};

int
BcfFormat_symbol_calc_qual(
        bcfrec::BcfFormat & fmt,
        
        const uvc1_readnum_t ins_cdepth,
        const uvc1_readnum_t del_cdepth,
        
        const uvc1_readnum_t ins1_cdepth,
        const uvc1_readnum_t del1_cdepth,
        
        const std::string & repeatunit,
        const uvc1_readpos_t repeatnum,
        const bool is_rescued,
        const auto & rtr1,
        const auto & rtr2,
        uvc1_refgpos_t tid,
        uvc1_refgpos_t refpos,
        AlignmentSymbol refsymbol,
        double tpfa,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    const int a = 0;
    const auto symbol = AlignmentSymbol(LAST(fmt.VTI));
    
    const PhredMutationTable sscs_mut_table(
            paramset.fam_phred_sscs_transition_CG_TA,
            paramset.fam_phred_sscs_transition_AT_GC,
            paramset.fam_phred_sscs_transversion_CG_AT,
            paramset.fam_phred_sscs_transversion_other,
            paramset.fam_phred_sscs_indel_open,
            paramset.fam_phred_sscs_indel_ext);
    
    const uvc1_qual_t powlaw_sscs_phrederr = sscs_mut_table.toPhredErrRate(refsymbol, symbol);
    const uvc1_qual_t powlaw_sscs_inc1 = (powlaw_sscs_phrederr  - (isSymbolSubstitution(symbol) 
            ? (((BASE_A == refsymbol && BASE_T == symbol) || (BASE_T == refsymbol && BASE_A == symbol))
                ? paramset.fam_phred_pow_sscs_transversion_AT_TA_origin : paramset.fam_phred_pow_sscs_snv_origin)
            : paramset.fam_phred_pow_sscs_indel_origin));
    
    const double t2n_contam_frac = (tpfa > 0 ? tpfa : 0) * paramset.contam_t2n_mul_frac;
    
    const double contamfrac = paramset.contam_any_mul_frac + (1.0 - paramset.contam_any_mul_frac) * t2n_contam_frac;
    
    const uvc1_readnum_t aDP = (fmt.aDPff[a] + fmt.aDPfr[a] + fmt.aDPrf[a] + fmt.aDPrr[a]);
    const uvc1_readnum_t ADP = (fmt.ADPff[0] + fmt.ADPrf[0] + fmt.ADPfr[0] + fmt.ADPrr[0]);
    const uvc1_readnum_t cDP0 = (fmt.cDP1f[a] + fmt.cDP1r[a]);
    const uvc1_readnum_t CDP0 = (fmt.CDP1f[0] + fmt.CDP1r[0]);
    const uvc1_readnum_t cDP2 = (fmt.cDP2f[a] + fmt.cDP2r[a]);
    const uvc1_readnum_t CDP2 = (fmt.CDP2f[0] + fmt.CDP2r[0]);
    
    const uvc1_qual_t aavgMQ = (uvc1_qual_t)(fmt.aMQs[a] / MAX(1, aDP));
    const auto diffAaMQs = (uvc1_qual_t)((fmt.AMQs[0] - fmt.aMQs[a]) / MAX(1, ADP - aDP)) - aavgMQ;
    const auto tn_q_inc_max = paramset.tn_q_inc_max;
    
    const auto noUMI_bias_inc = MIN(paramset.bias_FA_powerlaw_noUMI_phred_inc_snv, aDP/2); 
    const auto pl_noUMI_phred_inc = paramset.powlaw_anyvar_base 
            + (isSymbolSubstitution(symbol) ? noUMI_bias_inc : paramset.bias_FA_powerlaw_noUMI_phred_inc_indel);
    const auto withUMI_bias_inc = MIN(paramset.bias_FA_powerlaw_withUMI_phred_inc_snv - paramset.bias_FA_powerlaw_noUMI_phred_inc_snv, cDP2/2) + noUMI_bias_inc; 
    const auto pl_withUMI_phred_inc = paramset.powlaw_anyvar_base 
            + (isSymbolSubstitution(symbol) ? withUMI_bias_inc : paramset.bias_FA_powerlaw_withUMI_phred_inc_indel);
    
    const bool is_cytosine_deanim_CT = ((BASE_C == refsymbol && BASE_T == symbol) || (BASE_G == refsymbol && BASE_A == symbol));
    const bool is_cytosine_deanim_FA = (int64mul(cDP2, paramset.microadjust_fam_lowfreq_invFA) < (int64_t)CDP2);
    const bool is_cysotine_deanimation = (is_cytosine_deanim_CT && is_cytosine_deanim_FA);
    uvc1_qual_big_t non_duplex_binom_dec_x10 = ((is_cysotine_deanimation)
            ? (10L * non_neg_minus(paramset.fam_phred_sscs_transition_CG_TA, paramset.fam_phred_dscs_all / 2L) * ((uvc1_qual_big_t)CDP2 - int64mul(cDP2, paramset.microadjust_fam_lowfreq_invFA)) / CDP2) : 0);
    double prior_weight = 1.0 / (fmt.cDPmf[a] + fmt.cDPmr[a] + 1.0);
    const uvc1_qual_t fam_thres_highBQ = (isSymbolSubstitution(symbol) ? paramset.fam_thres_highBQ_snv : paramset.fam_thres_highBQ_indel); 
    const uvc1_qual_t cMmQ = (uvc1_qual_t)round(numstates2phred((fmt.cDPMf[a] + fmt.cDPmf[a] + fmt.cDPMr[a] + fmt.cDPmr[a] + pow(10, fam_thres_highBQ / 10.0) * prior_weight) 
            / (fmt.cDPmf[a] + fmt.cDPmr[a] + prior_weight)));
    
    uvc1_readnum100x_t nbases_x100_1 = fmt.bIADb[a] * 100 + 1;
    uvc1_readnum100x_t nbases_x100_2 = MIN(nbases_x100_1, fmt.cDP1v[a] + 1);
    uvc1_qual_big_t perbase_likeratio_q_x10_1 = 10 * fmt.bIAQb[a] / MAX(1, fmt.bIADb[a]);
    uvc1_qual_big_t perbase_likeratio_q_x10_2 = perbase_likeratio_q_x10_1 + (uvc1_qual_big_t)round(10 * numstates2phred((double)nbases_x100_2 / (double)nbases_x100_1));
    uvc1_qual_big_t duped_frag_binom_qual = ((isSymbolIns(symbol) || isSymbolDel(symbol)) ? perbase_likeratio_q_x10_1 : perbase_likeratio_q_x10_2)  * nbases_x100_2 / (10 * 100);
    uvc1_qual_big_t contam_frag_withmin_qual = (uvc1_qual_big_t)round(calc_binom_10log10_likeratio(t2n_contam_frac, cDP0, CDP0 - cDP0)) + 9 - 3;
    
    uvc1_qual_t phred_het3al_chance_inc_snp = MAX(0, 2 * paramset.germ_phred_hetero_snp - paramset.germ_phred_het3al_snp - TIN_CONTAM_MICRO_VQ_DELTA);
    uvc1_qual_t phred_het3al_chance_inc_indel = MAX(0, 2 * paramset.germ_phred_hetero_indel - paramset.germ_phred_het3al_indel - TIN_CONTAM_MICRO_VQ_DELTA);
    // systematic error affects both T and N upstream to the same extent, whereas T-in-N contamination affects only the N from T.
    uvc1_qual_t contam_syserr_phred_bypassed = (isSymbolSubstitution(symbol) ? phred_het3al_chance_inc_snp : phred_het3al_chance_inc_indel);
    const uvc1_readnum_t normCDP1 = (fmt.CDP12f[0] + fmt.CDP12r[0] + 1);
    const uvc1_readnum_t normBDP = (fmt.BDPf[0] + fmt.BDPr[0] + 1);
    const uvc1_qual_big_t sscs_dec1a = 2 * int64mul(non_neg_minus(paramset.fam_min_n_copies, normCDP1),                             (powlaw_sscs_inc1 + 3))
            / (paramset.fam_min_n_copies);
    const uvc1_qual_big_t sscs_dec1b = 2 * int64mul(non_neg_minus(normCDP1 * (paramset.fam_min_overseq_perc + 100) / 100, normBDP), (powlaw_sscs_inc1 + 3))
            / MAX(1, normCDP1 * (paramset.fam_min_overseq_perc / 2 + 100) / 100);
    const uvc1_qual_big_t sscs_dec1 = MAX(sscs_dec1a, sscs_dec1b);
    
    // how to reduce qual by bias: div first, then minus
    const uvc1_qual_big_t sscs_dec2 = non_neg_minus(fam_thres_highBQ, cMmQ);
    const uvc1_qual_big_t cIADnormcnt = int64mul(fmt.cIADf[a] + fmt.cIADr[a], 100) + 1;
    const uvc1_qual_big_t cIADmincnt = MIN(cIADnormcnt, 
            int64mul((fmt.CDP2f[0] + fmt.CDP2r[0]) * 100UL + 1, fmt.cDP1v[a] + 1) 
            / (fmt.CDP1f[0] + fmt.CDP1r[0] + 1) + 1);
    
    const uvc1_qual_big_t sscs_binom_qual_fw = fmt.cIAQf[a] + int64mul(fmt.cIAQr[a], MIN(paramset.fam_phred_dscs_all - fmt.cIDQf[a], fmt.cIDQr[a])) / MAX(fmt.cIDQr[a], 1);
    const uvc1_qual_big_t sscs_binom_qual_rv = fmt.cIAQr[a] + int64mul(fmt.cIAQf[a], MIN(paramset.fam_phred_dscs_all - fmt.cIDQr[a], fmt.cIDQf[a])) / MAX(fmt.cIDQf[a], 1);
    
    const uvc1_qual_big_t contam_sscs_withmin_qual = (uvc1_qual_big_t)round(calc_binom_10log10_likeratio(t2n_contam_frac, cDP2, CDP2 - cDP2)) + 9 - 3;
    
    uvc1_qual_big_t sscs_binom_qual = int64mul(MAX(sscs_binom_qual_fw, sscs_binom_qual_rv), cIADmincnt) / (cIADnormcnt) - (non_duplex_binom_dec_x10) * cIADmincnt;
    if (MAX(sscs_binom_qual_fw, sscs_binom_qual_rv) > paramset.microadjust_fam_binom_qual_halving_thres && isSymbolSubstitution(symbol)) {
        sscs_binom_qual = MIN(sscs_binom_qual, 
                paramset.microadjust_fam_binom_qual_halving_thres 
                + (MAX(sscs_binom_qual_fw, sscs_binom_qual_rv) - paramset.microadjust_fam_binom_qual_halving_thres) / 2); // intrinsic bias
    }
    sscs_binom_qual -= sscs_dec1 + sscs_dec2;
    double min_bcFA_v = (((double)(fmt.cDP1v[a]) + 0.5) / (double)(fmt.CDP1f[0] * 100 + fmt.CDP1r[0] * 100 + 1.0));
    uvc1_qual_t dedup_frag_powlaw_qual_v = round(paramset.powlaw_exponent * numstates2phred(min_bcFA_v) + (pl_noUMI_phred_inc));
    
    double min_bcFA_w = (((double)(fmt.cDP1w[a]) + 0.5) / (double)(fmt.CDP1f[0] * 100 + fmt.CDP1r[0] * 100 + 1.0));
    uvc1_qual_t dedup_frag_powlaw_qual_w = round(paramset.powlaw_exponent * numstates2phred(min_bcFA_w) 
            + (pl_noUMI_phred_inc) + tn_q_inc_max);
    
    uvc1_qual_t ds_vq_inc_powlaw = round(10/log(10)*MIN(log((fmt.cDP12f[a] + 0.5) / (fmt.CDP12f[0] + 1.0)), log((fmt.cDP12r[a] + 0.5) / (fmt.CDP12r[0] + 1.0)))) 
            + (powlaw_sscs_phrederr);
    uvc1_qual_t ds_vq_inc_binom = 3 * MIN(fmt.cDP2f[a], fmt.cDP2r[a]);
    
    uvc1_qual_t powlaw_sscs_inc2 = MAX(0, MIN5(sscs_binom_qual_fw, sscs_binom_qual_rv, ds_vq_inc_powlaw, ds_vq_inc_binom, 3));
    double umi_cFA =  (((double)(fmt.cDP2v[a]) + 0.5) / ((double)(fmt.CDP2f[0] * 100 + fmt.CDP2r[0] * 100) + 1.0));
    
    uvc1_qual_t sscs_base_2 = pl_withUMI_phred_inc + powlaw_sscs_inc1 + powlaw_sscs_inc2 - sscs_dec1 - sscs_dec2;
    uvc1_qual_t sscs_powlaw_qual_v = round((paramset.powlaw_exponent * numstates2phred(umi_cFA)    + sscs_base_2));
    uvc1_qual_t sscs_powlaw_qual_w = round((paramset.powlaw_exponent * numstates2phred(min_bcFA_w) + sscs_base_2)
            + tn_q_inc_max);
        
    double dFA = (double)(fmt.dDP2[a] + 0.5) / (double)(fmt.DDP1[0] + 1.0);
    double dSNR = (double)(fmt.dDP2[a] + 0.5) / (double)(fmt.dDP1[0] + 1.0);
    double dnormFA = dFA * pow(dSNR, 1.0 / paramset.powlaw_exponent);
    uvc1_qual_big_t dFA_vq_binom = ((uvc1_qual_big_t)paramset.fam_phred_dscs_all - (uvc1_qual_big_t)round(numstates2phred(1.0 / (dnormFA))))
            * (uvc1_qual_big_t)fmt.dDP2[a] * (uvc1_qual_big_t)cIADmincnt / (uvc1_qual_big_t)cIADnormcnt;
    uvc1_qual_t dFA_vq_powlaw = paramset.powlaw_anyvar_base + (paramset.fam_phred_dscs_all - paramset.fam_phred_pow_dscs_all_origin) 
            + (uvc1_qual_t)round(numstates2phred((dnormFA) * MIN(1.0, (double)((fmt.cDP1v[a]) + 0.5) / (double)(fmt.CDP1f[0] * 100 + fmt.CDP1r[0] * 100 + 1.0))));
    
// This code is supposed to work given highly accurate estimation of InDel error.
// However, it seems that in practice, the estimation of InDel errors is not sufficiently accurate.
// Perhaps we can enable this code once an accurate method for estimating InDel errors is found.
#if 0 
    if (!isSymbolSubstitution(symbol) && (aDP / (double)(ADP + 0.5)) < 0.5) {
        int sscs_base_3 = powlaw_sscs_inc2 - sscs_dec1;
        const auto phred_label = paramset.fam_phred_indel_inc_before_barcode_labeling;
        sscs_powlaw_qual_v = (int)(paramset.powlaw_exponent * (numstates2phred(umi_cFA)    + cMmQ + phred_label) + sscs_base_3);
        sscs_powlaw_qual_w = (int)(paramset.powlaw_exponent * (numstates2phred(min_bcFA_w) + cMmQ + phred_label) + sscs_base_3) + (int)paramset.tn_q_inc_max;
    }
#endif
    clear_push(fmt.cMmQ, cMmQ);
    
    const double eps = (double)FLT_EPSILON;
    const bool is_indel_penal_applied = ((SEQUENCING_PLATFORM_IONTORRENT == paramset.inferred_sequencing_platform) && (NOT_PROVIDED == paramset.vcf_tumor_fname));
    const uvc1_qual_t indel_penal_base = (is_indel_penal_applied 
                ? ((uvc1_qual_t)round(paramset.indel_multiallele_samepos_penal / log(2) 
                * log((double)MAX3(aDP + eps, fmt.APDP[1], fmt.APDP[2]) / (double)(aDP + eps)))) 
                : 0);
    if (paramset.should_add_note) {
        fmt.note += std::string("/pb/") + std::to_string(indel_penal_base) + "/" 
                + std::to_string(paramset.inferred_sequencing_platform) + "/(" + paramset.vcf_tumor_fname + ")/";
    }
    uvc1_qual_t indel_penal4multialleles = 0;
    uvc1_qual_t indel_penal4multialleles_g = 0;
    uvc1_qual_t indel_penal4multialleles_soma = 0;
    const std::string & indelstring = fmt.gapSa[a];
    
    uvc1_qual_t indel_UMI_penal = 0;
    
    if (indelstring.size() > 0 && fmt.cDP0a[a] > 0) {
        const double indel_pq = (double)MIN(indel_phred(paramset.indel_polymerase_slip_rate, 
                UNSIGN2SIGN(repeatunit.size()), repeatnum), 24) + 2 - (double)10;
        const auto eff_tracklen1 = (UNSIGN2SIGN(repeatunit.size()) * UNSIGN2SIGN(MAX(1, repeatnum)) - UNSIGN2SIGN(repeatunit.size()));
        const auto eff_tracklen2 = (MAX(rtr1.tracklen - rtr1.unitlen, rtr2.tracklen - rtr2.unitlen) / 3);
        const double indel_ic = numstates2phred((double)MAX(indelstring.size() + (isSymbolIns(symbol) ? INS_N_ANCHOR_BASES : 0), (size_t)1) / (double)(MAX(UNSIGN2SIGN(eff_tracklen1), eff_tracklen2) + 1)) 
                + (isSymbolIns(symbol) ? (numstates2phred(paramset.indel_del_to_ins_err_ratio) * MIN(200, fmt.cDP0a[a]) / 200) : 0);
        auto indelcdepth = (isSymbolIns(symbol) ? ins_cdepth : del_cdepth);
        if (LINK_D1 == symbol) {
            indelcdepth += ins1_cdepth;
        }
        if (LINK_I1 == symbol) {
            indelcdepth += del1_cdepth / paramset.indel_del_to_ins_err_ratio;
        }
        const uvc1_readnum_t nearInDelDP = (isSymbolIns(symbol) ? fmt.APDP[1] : fmt.APDP[2]);
        
        assert (nearInDelDP >= aDP || !fprintf(stderr, "nearInDelDP >= aDP failed (%d >= %d failed) at tid %d pos %d!\n", nearInDelDP, aDP, tid, refpos));
                
        const auto indel_penal4multialleles1 = (uvc1_qual_t)round(paramset.indel_multiallele_samepos_penal / log(2.0) 
                * log((double)(indelcdepth + eps) / (double)(fmt.cDP0a[a] + eps)));
        const auto indel_penal4multialleles2 = (uvc1_qual_t)round(paramset.indel_multiallele_diffpos_penal / log(2.0) 
                * log((double)(nearInDelDP + eps) / (double)(aDP + eps)));
        indel_penal4multialleles_g = (uvc1_qual_t)round(paramset.indel_tetraallele_germline_penal_value /log(2.0) * log((double)(ins_cdepth + del_cdepth + eps) / (double)(fmt.cDP0a[a] + eps))) - paramset.indel_tetraallele_germline_penal_thres;
        
        if (isSymbolIns(symbol)) {
            indel_penal4multialleles = (indel_penal4multialleles1 * paramset.indel_ins_penal_pseudocount / 
                    (uvc1_qual_t)(paramset.indel_ins_penal_pseudocount + UNSIGN2SIGN(indelstring.size()))); 
            indel_penal4multialleles_soma = (indel_penal4multialleles1 * paramset.indel_ins_penal_pseudocount / 
                    (uvc1_qual_t)(paramset.indel_ins_penal_pseudocount + UNSIGN2SIGN(indelstring.size())));
        } else {
            indel_penal4multialleles = MAX(indel_penal4multialleles1, indel_penal4multialleles2);
            indel_penal4multialleles_soma = indel_penal4multialleles1;
        }
        dedup_frag_powlaw_qual_v += round(indel_ic);
        dedup_frag_powlaw_qual_w += round(indel_ic);
        duped_frag_binom_qual += round(indel_pq);
        
        // assume that first PCR-cycle InDel error rate approx equals to in-vivo error rate over multiple cell generations.
        const double sscs_indel_ic = numstates2phred((double)mathsquare(MAX(indelstring.size(), 1U)) / (double)(MAX(UNSIGN2SIGN(eff_tracklen1), eff_tracklen2) + 1));
        const uvc1_qual_t sscs_ins_vs_del_inc = round(paramset.powlaw_exponent * numstates2phred(paramset.indel_del_to_ins_err_ratio));
        const uvc1_qual_t extra_reward = non_neg_minus(
                sscs_ins_vs_del_inc,
                sscs_indel_ic * (isSymbolIns(symbol) ? 0 : MAX(UNSIGN2SIGN(eff_tracklen1), eff_tracklen2)) / round(paramset.indel_polymerase_size)) 
            - sscs_ins_vs_del_inc / 2;
        sscs_powlaw_qual_v += round(sscs_indel_ic) + extra_reward;
        sscs_powlaw_qual_w += round(sscs_indel_ic) + extra_reward;
        sscs_binom_qual += round(indel_pq) + extra_reward;
        indel_UMI_penal = non_neg_minus((fmt.BDPf[0] + fmt.BDPr[0] + 1.0) / (double)(fmt.CDP1f[0] + fmt.CDP1r[0] + 1.0) * paramset.fam_indel_nonUMI_phred_dec_per_fold_overseq,
                (paramset.fam_thres_emperr_all_flat_indel + 1) * paramset.fam_indel_nonUMI_phred_dec_per_fold_overseq);
    }
    
    clear_push(fmt.aAaMQ, diffAaMQs);
    
    const uvc1_qual_t readlenMQcap = (fmt.APXM[2]) / MAX(1, fmt.APDP[0]) - 17;
    const uvc1_qual_t diffMQ = (MAX(0, diffAaMQs));
    const bool is_aln_extra_accurate = (paramset.inferred_maxMQ > 60);
    const uvc1_qual_t _systematicMQVQadd = (uvc1_qual_t)((symbol == refsymbol) 
            ? 0 : (MIN(paramset.germ_phred_homalt_snp, ADP * 3)));
    const bool is_MQ_unadjusted = (is_aln_extra_accurate || (!isSymbolSubstitution(symbol)) || (aDP > ADP * 3/4));
    const uvc1_qual_t _systematicMQVQminus =  
             (is_MQ_unadjusted ? 0 : (non_neg_minus((60- 30), aavgMQ) * 2 / 5))
          + ((is_MQ_unadjusted || (refsymbol != symbol)) ? 0 : non_neg_minus(MIN(15, diffMQ), aavgMQ));
    
    uvc1_qual_t diffMQ2 = diffMQ;
    if (fmt.bMQ[a] < 20 && (NOT_PROVIDED == paramset.vcf_tumor_fname)) {
        const auto aDPxf = (fmt.aDPff[a] + fmt.aDPrf[a] + 0.5);
        const auto aDPxr = (fmt.aDPfr[a] + fmt.aDPrr[a] + 0.5);
        const auto ADPxf = (fmt.ADPff[0] + fmt.ADPrf[0] + 1.0);
        const auto ADPxr = (fmt.ADPfr[0] + fmt.ADPrr[0] + 1.0);
        if ((aDPxr / ADPxr) * 2 < (aDPxf / ADPxf) || (aDPxf / ADPxf) * 2 < (aDPxr / ADPxr)
                || (fmt.aLI1[a] + 0.5) / (fmt.ALI2[0] + 1.0) * (2 * (1.0+DBL_EPSILON)) < (aDPxr) / (ADPxr)
                || (fmt.aRI1[a] + 0.5) / (fmt.ARI2[0] + 1.0) * (2 * (1.0+DBL_EPSILON)) < (aDPxf) / (ADPxf)) {
            UPDATE_MAX(diffMQ2, 20 - MIN(fmt.bMQ[a], 20));
        }
    }
    
    const uvc1_qual_t _systematicMQ = (((refsymbol == symbol) && (ADP > aDP * 2))
                ? fmt.bMQ[a] // small
                : (fmt.bMQ[a] * (paramset.syserr_MQ_max - paramset.syserr_MQ_nonref_base) / paramset.syserr_MQ_max + paramset.syserr_MQ_nonref_base))
            - (uvc1_qual_t)(diffMQ2)
            - (uvc1_qual_t)(fmt.bNMQ[a])
            - (uvc1_qual_t)(numstates2phred((ADP + 1.0) / (aDP + 0.5)));
    auto systematicMQVQ1 = MIN((MAX(_systematicMQ, paramset.syserr_MQ_min) + _systematicMQVQadd), readlenMQcap);
    
    const uvc1_qual_t systematicBQVQ = (
            ((SEQUENCING_PLATFORM_IONTORRENT != paramset.inferred_sequencing_platform) && isSymbolSubstitution(AlignmentSymbol(LAST(fmt.VTI)))) 
            ? fmt.aBQQ[a] : (200));
    
    // begin ad-hoc
    const bool is_strong_amplicon = (fmt.APDP[5] * 100 > fmt.APDP[0] * 60);
    const bool is_weak_amplicon = (fmt.APDP[5] * 100 > fmt.APDP[0] * 40);
    const bool is_tmore_amplicon = ((NOT_PROVIDED == paramset.vcf_tumor_fname) ? is_weak_amplicon : is_strong_amplicon);
    
    if (is_tmore_amplicon && (isSymbolIns(symbol) || isSymbolDel(symbol)) && (systematicBQVQ > 70)
            && (fmt.APXM[1] / MAX(fmt.APDP[0], 1) > 20)) {
        systematicMQVQ1 = 70 + ((systematicMQVQ1 - 70) * 5 / (fmt.APXM[1] / MAX(fmt.APDP[0], 1) - 15));
    }
    uvc1_qual_t indel_penal_base_add = 0;
    if (NOT_PROVIDED == paramset.vcf_tumor_fname) {
        const auto delAPDP = MAX(fmt.APDP[2], fmt.APDP[4]);
        if ((fmt.APDP[0] < 3 * delAPDP) && (fmt.APDP[0] < 3 * fmt.APDP[6]) && (aDP * 3 < delAPDP) && (aDP * 3 < fmt.APDP[6]) 
                && isSymbolSubstitution(symbol) && (rtr2.tracklen >= 8 * rtr2.unitlen)) {
            indel_penal_base_add = 9;
        }
        if (is_tmore_amplicon && (isSymbolDel(symbol))) {
            if (aDP * 4 < fmt.APDP[2]) {
                indel_penal_base_add = 5;
            } else if (fmt.cDP0a[a] * 3 < 2 * (del_cdepth)) {
                indel_penal_base_add = 2;
            }
        }
    }
    // end ad-hoc

    const auto systematicMQVQ = MAX(0, systematicMQVQ1);
    const auto indel_penal_base2 = indel_penal_base + indel_penal_base_add;
    
    const auto tn_syserr_q = systematicMQVQ + paramset.tn_q_inc_max;
    clear_push(fmt.bMQQ, systematicMQVQ);
    
    clear_push(fmt.bIAQ, duped_frag_binom_qual - indel_penal_base2, a);
    clear_push(fmt.cIAQ, sscs_binom_qual - indel_penal_base, a);
    
    clear_push(fmt.cPCQ1, MIN(dedup_frag_powlaw_qual_w - indel_penal_base2, tn_syserr_q), a);
    clear_push(fmt.cPLQ1, dedup_frag_powlaw_qual_v - indel_penal_base2, a);
    
    clear_push(fmt.cPCQ2, MIN(sscs_powlaw_qual_w - indel_penal_base, tn_syserr_q), a);
    clear_push(fmt.cPLQ2, sscs_powlaw_qual_v - indel_penal_base, a);
    
    clear_push(fmt.bTINQ, contam_frag_withmin_qual + contam_syserr_phred_bypassed, a);
    clear_push(fmt.cTINQ, contam_sscs_withmin_qual + contam_syserr_phred_bypassed, a);
    
    const uvc1_readnum_t aDPpc = ((refsymbol == symbol) ? 1 : 0);
    const uvc1_qual_t penal4BQerr = (isSymbolSubstitution(symbol) ? (5 + UNSIGN2SIGN(((int64_t)paramset.penal4lowdep) / (int64_t)mathsquare((int64_t)MAX(1, aDP + aDPpc)))) : 0);
    
    const uvc1_qual_t indel_q_inc = ((((!isSymbolIns(symbol)) && (!isSymbolDel(symbol))) || is_rescued) ? 0 : indel_len_rusize_phred(indelstring.size(), repeatnum));
    clear_push(fmt.gVQ1, MAX(
            0, 
            indel_q_inc + MIN3(
                MIN(systematicBQVQ, non_neg_minus(systematicMQVQ, _systematicMQVQminus)),
                LAST(fmt.bIAQ) - penal4BQerr, 
                LAST(fmt.cPLQ1)) - 2 * MAX3(
                    0, 
                    indel_penal4multialleles - paramset.indel_multiallele_soma_penal_thres, 
                    indel_penal4multialleles_g)),
    a);
    const uvc1_qual_t systematicVQsomatic_minus = (is_rescued ? 0 : (15 - MIN3(ADP * 15 / 100, aDP, 15)));
    const uvc1_qual_t systematicVQsomatic = non_neg_minus(MIN(systematicBQVQ, systematicMQVQ), systematicVQsomatic_minus);
    const auto bcVQ1 = MIN3(
                systematicVQsomatic,
                LAST(fmt.bIAQ) - (is_rescued ? 0 : penal4BQerr), 
                LAST(fmt.cPLQ1)) - indel_penal4multialleles_soma;
    clear_push(fmt.cVQ1, MAX(0, MIN(bcVQ1, LAST(fmt.bTINQ)) - indel_UMI_penal), a);
    
    // This adjustment makes sure that variant with UMI support is ranked before variant without UMI suppport.
    const std::array<uvc1_qual_t, 10> cysotine_deanim_to_score = {{0,1,1, 1,1,2, 2,2,2, 3}};
    const std::array<uvc1_qual_t, 5> other_to_score = {{0,2,3, 3,4}};
    uvc1_qual_t mincVQ2 = (is_cytosine_deanim_CT 
        ? cysotine_deanim_to_score[MIN(cDP2, UNSIGN2SIGN(cysotine_deanim_to_score.size()) - 1)] 
        : other_to_score[MIN(cDP2, UNSIGN2SIGN(other_to_score.size()) - 1)]); 
    
    if (isSymbolIns(symbol) || isSymbolDel(symbol)) {
        const uvc1_qual_t sscs_floor_qual_v = MIN(paramset.germ_phred_homalt_indel + numstates2phred(umi_cFA), fmt.cDP2v[a] * 3 / 100) + ((isSymbolIns(symbol) ? INS_N_ANCHOR_BASES : 0) - INS_N_ANCHOR_BASES) * 3;
        UPDATE_MAX(mincVQ2, sscs_floor_qual_v);
    }
    
    const uvc1_qual_big_t dVQinc = MIN(dFA_vq_binom, dFA_vq_powlaw) - MAX(0, MIN(LAST(fmt.cIAQ), LAST(fmt.cPLQ2)));
    clear_push(fmt.dVQinc, dVQinc, a);
    
    const uvc1_qual_t cVQ2 = MIN3(systematicVQsomatic,
            LAST(fmt.cIAQ) + MAX(0, dVQinc),
            LAST(fmt.cPLQ2) + MAX(0, dVQinc)) - indel_penal4multialleles;
    clear_push(fmt.cVQ2, MAX(mincVQ2, MIN(cVQ2, LAST(fmt.cTINQ))), a);
    
    // TODO; check if reducing all allele read count to increase alt allele frac in case of ref bias makes more sense
    const auto & cDP1y = (is_rescued ? fmt.cDP1x : fmt.cDP1v); // it was always not is_rescued before.
    const auto & CDP1y = (is_rescued ? fmt.CDP1x : fmt.CDP1v);
    auto binom_contam_LODQ = calc_binom_10log10_likeratio(contamfrac, cDP1y[a], CDP1y[0]); // it was cDP1v before
    auto power_contam_LODQ = round(10.0/log(10.0) * paramset.powlaw_exponent * 
            MAX(logit2((cDP1y[a] + 1) / (double)(CDP1y[0] + 1), contamfrac), 0.0)); // it was min_bcFA_v before
    clear_push(fmt.CONTQ, MIN(binom_contam_LODQ, power_contam_LODQ), a);
    return 0;
};

#define INDEL_ID 1
#include "instcode.hpp"
#undef INDEL_ID
#define INDEL_ID 2
#include "instcode.hpp"
std::array<uvc1_readnum_t, 2>
fill_by_indel_info(
        bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet, 
        const int strand, 
        const uvc1_refgpos_t refpos, 
        const AlignmentSymbol symbol, 
        const std::string & refstring, 
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {

    assert(isSymbolIns(symbol) || isSymbolDel(symbol));
    if (isSymbolIns(symbol)) {
        return fill_by_indel_info2_1(fmt, symbol2CountCoverageSet, strand, refpos, symbol,
                symbol2CountCoverageSet.symbol_to_frag_format_depth_sets.at(strand).getPosToIseqToData(symbol),
                symbol2CountCoverageSet.symbol_to_fam_format_depth_sets_2strand.at(strand).getPosToIseqToData(symbol),
                refstring, 
                0);
    } else {
        return fill_by_indel_info2_2(fmt, symbol2CountCoverageSet, strand, refpos, symbol,
                symbol2CountCoverageSet.symbol_to_frag_format_depth_sets.at(strand).getPosToDlenToData(symbol),
                symbol2CountCoverageSet.symbol_to_fam_format_depth_sets_2strand.at(strand).getPosToDlenToData(symbol),
                refstring, 
                0);
    }
};

std::string 
mutform2count4map_to_phase(const auto & mutform2count4vec, const auto & indices, uvc1_readnum_t pseudocount = 1) {
    std::string phase_string;
    for (auto idx : indices) {
        auto mutform2count4pair = mutform2count4vec.at(idx);
        auto counts = mutform2count4pair.second;
        if ((counts[0] + counts[1]) > pseudocount) {
            phase_string +="(";
            for (auto pos2symbol4pair : mutform2count4pair.first) {
                AlignmentSymbol symbol = pos2symbol4pair.second;
                uvc1_refgpos_t mutpos = pos2symbol4pair.first + (isSymbolSubstitution(symbol) ? 1 : 0);
                phase_string += std::string("(") + std::to_string(mutpos) + "&" + SYMBOL_TO_DESC_ARR[symbol] + ")";
            }
            phase_string += std::string("&") + std::to_string(counts[0]) + "&" + std::to_string(counts[1]) + ")";
        }
    }
    return phase_string;
}

const std::vector<std::pair<std::array<uvc1_readnum_t, 2>, std::string>> 
indel_get_majority(
        const bcfrec::BcfFormat & fmt, 
        const char *tname, 
        uvc1_refgpos_t refpos, 
        const AlignmentSymbol symbol,
        bool is_warning_disabled,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    std::vector<std::pair<std::array<uvc1_readnum_t, 2>, std::string>> indelstrings;
    if ((SUMVEC(fmt.gapNf) + SUMVEC(fmt.gapNr)) == 0) {
        if (!is_warning_disabled) {
            std::cerr << "Invalid indel detected (invalid mutation) : " << tname << ", " << refpos << ", " << SYMBOL_TO_DESC_ARR[symbol] << std::endl;
            std::string msg;
            bcfrec::streamAppendBcfFormat(msg, fmt);
            std::cerr << msg << "\n";
        }
        indelstrings.push_back(std::make_pair(std::array<uvc1_readnum_t, 2>({{0, 0}}), SYMBOL_TO_DESC_ARR[symbol]));
    } else {
        std::map<std::string, std::array<uvc1_readnum_t, 2>> indelmap;
        assert(UNSIGN2SIGN(fmt.gapSeq.size()) == (SUMVEC(fmt.gapNf) + SUMVEC(fmt.gapNr)));
        assert(fmt.gapbAD1.size() == fmt.gapSeq.size());
        for (size_t i = 0; i < fmt.gapSeq.size(); i++) {
            const auto & it = indelmap.find(fmt.gapSeq[i]);
            if (it == indelmap.end()) {
                indelmap.insert(std::make_pair(fmt.gapSeq[i], std::array<uvc1_readnum_t, 2>({{fmt.gapbAD1[i], fmt.gapcAD1[i]}})));
            } else {
                indelmap[fmt.gapSeq[i]][0] += fmt.gapbAD1[i];
                indelmap[fmt.gapSeq[i]][1] += fmt.gapcAD1[i];
            }
        }
        uvc1_readnum_t max_bAD1 = 0;
        for (auto & indelit : indelmap) {
            UPDATE_MAX(max_bAD1, indelit.second[0]);
        }
        for (auto & indelit : indelmap) {
            if (indelit.second[0] >= (max_bAD1 + 3) / 4) {
                indelstrings.push_back(std::make_pair(indelit.second, indelit.first));
            }
        }
        struct InDelLess {
            bool operator() (
                    const std::pair<std::array<uvc1_readnum_t, 2>, std::string>& x, 
                    const std::pair<std::array<uvc1_readnum_t, 2>, std::string>& y) const {
                return ((mathsquare((int64_t)x.first[0]) * (int64_t)x.second.size()) < (mathsquare((int64_t)y.first[0]) * (int64_t)y.second.size()));
            }
        } InDelLess;
        std::sort(indelstrings.rbegin(), indelstrings.rend(), InDelLess);
    }
    return indelstrings;
}

uvc1_qual_t
hetLODQ(double allele1count, double allele2count, double expfrac, double powlaw_exponent = 3.0) {
    uvc1_qual_t binomLODQ = calc_binom_10log10_likeratio(expfrac, allele1count, allele2count);
    uvc1_qual_t powerLODQ = round(10.0/log(10.0) * powlaw_exponent * MAX(logit2((allele1count + 0.5) * 0.5 / expfrac, (allele2count + 0.5) * 0.5 / (1.0 - expfrac)), 0.0));
    return MIN(binomLODQ, powerLODQ);
}

struct {
    template<class T1, class T2>
    bool operator()(std::pair<T1, T2> a, std::pair<T1, T2> b) const {   
        return (a.second < b.second) || (a.second == b.second && a.first < b.first);
    }
} PairSecondLess;

double 
compute_norm_ad(const bcfrec::BcfFormat *fmtp,
    const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    return fmtp->cDP1v[fmtp->cDP1v.size()-1] / 100.0;
}

const auto
ALODQ(const auto x) {
    return x->gVQ1[x->gVQ1.size() - 1];
};

auto
output_germline(
        std::string & out_string,
        AlignmentSymbol refsymbol, 
        std::vector<std::pair<AlignmentSymbol, bcfrec::BcfFormat*>> symbol_format_vec,
        const char *tname,
        const std::string & refstring,
        uvc1_refgpos_t refpos,
        uvc1_refgpos_t extended_inclu_beg_pos,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    const bool is_rescued = (NOT_PROVIDED != paramset.vcf_tumor_fname);
    assert(symbol_format_vec.size() >= 4 || 
            !fprintf(stderr, " The variant-type %s:%d %d has symbol_format_vec of length %lu", 
            tname, refpos, refsymbol, symbol_format_vec.size()));
    uvc1_refgpos_t regionpos = refpos- extended_inclu_beg_pos;
    struct {
        bool operator()(std::pair<AlignmentSymbol, bcfrec::BcfFormat*> & p1, std::pair<AlignmentSymbol, bcfrec::BcfFormat*> & p2) const {
            return ALODQ(p1.second) < ALODQ(p2.second);
        }
    } SymbolBcfFormatPairLess;
    std::sort(symbol_format_vec.rbegin(), symbol_format_vec.rend(), SymbolBcfFormatPairLess);
    std::array<std::pair<AlignmentSymbol, bcfrec::BcfFormat*>, 4> ref_alt1_alt2_alt3 = {{std::make_pair(END_ALIGNMENT_SYMBOLS, (bcfrec::BcfFormat*)NULL)}};
    int allele_idx = 1;
    uvc1_qual_t ref_alodq = INT_MIN;
    for (auto symb_fmt : symbol_format_vec) {
        const bool isref = (refsymbol == symb_fmt.first || BASE_NN == symb_fmt.first || LINK_NN == symb_fmt.first); 
        if (isref && ALODQ(symb_fmt.second) > ref_alodq) {
            ref_alt1_alt2_alt3[0] = symb_fmt;
            ref_alodq = ALODQ(symb_fmt.second);
        }
        if ((!isref) && allele_idx <= 3) {
            ref_alt1_alt2_alt3[allele_idx] = symb_fmt;
            allele_idx++;
        }
    }
    assert(ref_alt1_alt2_alt3[0].second != NULL);
    assert(ref_alt1_alt2_alt3[3].second != NULL || !fprintf(stderr, "%d %d %d %d is invalid!\n", 
        ref_alt1_alt2_alt3[0].first, 
        ref_alt1_alt2_alt3[1].first,
        ref_alt1_alt2_alt3[2].first, 
        ref_alt1_alt2_alt3[3].first
        ));
    
    uvc1_qual_t a0LODQ = ALODQ(ref_alt1_alt2_alt3[0].second);
    uvc1_qual_t a1LODQ = ALODQ(ref_alt1_alt2_alt3[1].second);
    uvc1_qual_t a2LODQ = ALODQ(ref_alt1_alt2_alt3[2].second);
    uvc1_qual_t a3LODQ = ALODQ(ref_alt1_alt2_alt3[3].second);
    
    auto fmtptr0 = ref_alt1_alt2_alt3[0].second;
    auto fmtptr1 = ref_alt1_alt2_alt3[1].second;
    auto fmtptr2 = ref_alt1_alt2_alt3[2].second;
    const bool isSubst = isSymbolSubstitution(refsymbol);
    const AlignmentSymbol symbolNN = ((isSubst || !is_rescued) ? BASE_NN : LINK_NN);
    double ad0norm = compute_norm_ad(fmtptr0, 0);
    double ad1norm = compute_norm_ad(fmtptr1, 0);
    double ad2norm = compute_norm_ad(fmtptr2, 0);
    if (symbolNN == ref_alt1_alt2_alt3[1].first) {
        ad0norm += ad1norm;
        ad1norm = 0;
    }
    if (symbolNN == ref_alt1_alt2_alt3[2].first) {
        ad0norm += ad2norm;
        ad2norm = 0;
    }
    
    uvc1_qual_t a0a1LODQ = hetLODQ(ad0norm, ad1norm, 1.0 - paramset.germ_hetero_FA, paramset.powlaw_exponent);
    uvc1_qual_t a1a0LODQ = hetLODQ(ad1norm, ad0norm, paramset.germ_hetero_FA, paramset.powlaw_exponent);
    uvc1_qual_t a1a2LODQ = hetLODQ(ad1norm, ad2norm, 0.5, paramset.powlaw_exponent);
    uvc1_qual_t a2a1LODQ = hetLODQ(ad2norm, ad1norm, 0.5, paramset.powlaw_exponent);
    
    std::array<std::string, 4> GTidx2GT {{
        "0/0",
        "0/1",
        "1/1",
        "1/2"
    }};
    
    uvc1_qual_t phred_homref = 0;
    uvc1_qual_t phred_hetero = (isSubst ? paramset.germ_phred_hetero_snp : paramset.germ_phred_hetero_indel);
    uvc1_qual_t phred_homalt = (isSubst ? paramset.germ_phred_homalt_snp : paramset.germ_phred_homalt_indel);
    uvc1_qual_t phred_tri_al = (isSubst ? paramset.germ_phred_het3al_snp : paramset.germ_phred_het3al_indel); 
    // https://www.genetics.org/content/184/1/233 and https://www.fsigenetics.com/article/S1872-4973(20)30003-X/fulltext : triallelic-SNP-phred = 29*2-3+5
    
    const uvc1_qual_t contqadd = 0; // (phred_hetero - paramset.germ_phred_hetero_snp); // this should be enabled only if the alignment is highly accurate.
    if (is_rescued) {
        UPDATE_MIN(a0LODQ, ref_alt1_alt2_alt3[0].second->CONTQ[1] + contqadd);
        UPDATE_MIN(a1LODQ, ref_alt1_alt2_alt3[1].second->CONTQ[1]);
        UPDATE_MIN(a2LODQ, ref_alt1_alt2_alt3[2].second->CONTQ[1]);
        UPDATE_MIN(a3LODQ, ref_alt1_alt2_alt3[3].second->CONTQ[1]);
    } else {
        UPDATE_MIN(a0LODQ, ref_alt1_alt2_alt3[0].second->CONTQ[1] + contqadd);
    }
    // Two principles: 
    //  1. Every unexpected signal given a genotype reduces the likelihood of this genotype
    //  2. Signals are correlated with each other with a phred-scaled correlation of phred_hetero
    const auto a2penal = MAX(a2LODQ - (phred_tri_al - phred_hetero), 0);
    const auto a3penal = MAX(a3LODQ - phred_hetero, 0);
    const auto a01hetp = MAX(MAX(a0a1LODQ, a1a0LODQ) - (0-0), 0);
    const auto a12hetp = MAX(MAX(a1a2LODQ, a2a1LODQ) - (3-0), 0);
    const auto a03trip = MAX(a0LODQ, a3LODQ);
    
    uvc1_qual_t tri_al_penal = 0;
    const auto symb1 = AlignmentSymbol(LAST(ref_alt1_alt2_alt3[1].second->VTI));
    const auto symb2 = AlignmentSymbol(LAST(ref_alt1_alt2_alt3[2].second->VTI));
    if (isSymbolIns(symb1) && isSymbolIns(symb2)) {
        tri_al_penal += 3;
        if (LAST(ref_alt1_alt2_alt3[1].second->VTI) == LAST(ref_alt1_alt2_alt3[2].second->VTI)) {
            tri_al_penal += 3;
            if (LINK_I3P == AlignmentSymbol(LAST(ref_alt1_alt2_alt3[1].second->VTI))) {
                tri_al_penal += 3;
            }
        }
    }
    
    {
        const auto n1 = SYMBOL_TO_INDEL_N_UNITS[symb1];
        const auto n2 = SYMBOL_TO_INDEL_N_UNITS[symb2];
        if (n1 != 0 && n2 != 0) {
            tri_al_penal -= BETWEEN(abs(n1 - n2) * 3 - 5, 0, 9);
        }
    }
    std::array<std::pair<size_t, uvc1_qual_t>, 4> GL4raw = {{
        std::make_pair(0,     (-phred_homref - a1LODQ                - a2penal - a3penal)),
        std::make_pair(1,     (-phred_hetero - MAX(a01hetp, a2LODQ)  - MAX(MIN(a01hetp, a2LODQ)  - phred_hetero, 0) - a3penal)),
        std::make_pair(2,     (-phred_homalt - MAX( a0LODQ, a2LODQ)  - MAX(MIN( a0LODQ, a2LODQ)  - phred_hetero, 0) - a3penal)),
        std::make_pair(3,     (-phred_tri_al - MAX(a12hetp, a03trip) - MAX(MIN(a12hetp, a03trip) - phred_hetero, 0) 
            - MAX(MIN(a12hetp, MIN(a0LODQ, a3LODQ)) - phred_hetero, 0) - tri_al_penal))
    }};
    auto ret = GL4raw[0].second - MAX3(GL4raw[1].second, GL4raw[2].second, GL4raw[3].second);
    if (0 == ((OUTVAR_GERMLINE) & paramset.outvar_flag)) {
        return std::make_tuple(ret, fmtptr1, fmtptr2);
    }
    
    auto GL4 = GL4raw;
    std::sort(GL4.rbegin(), GL4.rend(), PairSecondLess);
    
    size_t GLidx = GL4[0].first;
    if (0 == GLidx && (!paramset.should_output_all_germline) && MAX(LAST(fmtptr1->cDP0a), LAST(fmtptr2->cDP0a)) <= 2) {
        return std::make_tuple(ret, fmtptr1, fmtptr2);
    }
    std::array<std::string, 3> ref_alt1_alt2_vcfstr_arr = {{
        SYMBOL_TO_DESC_ARR[ref_alt1_alt2_alt3[0].first],
        SYMBOL_TO_DESC_ARR[ref_alt1_alt2_alt3[1].first],
        SYMBOL_TO_DESC_ARR[ref_alt1_alt2_alt3[2].first]
    }};
    
    std::string vcfref = "";
    std::string vcfalt = "";
    uvc1_qual_t alt1_uniallelic_phred = 200;
    bool is_alt1_uniallelic = true;
    if (isSymbolSubstitution(refsymbol)) {
        assert(refstring.substr(regionpos, 1) == ref_alt1_alt2_vcfstr_arr[0]);
        vcfref = std::string(ref_alt1_alt2_vcfstr_arr[0]);
        vcfalt = std::string(ref_alt1_alt2_vcfstr_arr[1]);
        if (3 == GLidx) {
            vcfalt += std::string(",") + std::string(ref_alt1_alt2_vcfstr_arr[2]);
        }
    } else {
        ref_alt1_alt2_vcfstr_arr[0] = (regionpos > 0 ? refstring.substr(regionpos-1, 1) : "n");
        AlignmentSymbol s1 = ref_alt1_alt2_alt3[1].first;
        std::vector<std::pair<std::array<uvc1_readnum_t, 2>, std::string>> indelstrings1 = indel_get_majority(
                *(ref_alt1_alt2_alt3[1].second), 
                tname, 
                refpos, 
                s1, 
                is_rescued,
                0);
        const std::string & indelstring1 = indelstrings1[0].second;
        if (indelstrings1.size() > 1) {
            alt1_uniallelic_phred = (phred_tri_al - phred_hetero) *  log(1.0 + (double)indelstrings1[0].first[1] / (double)indelstrings1[1].first[1]) / log(2.0);
        }
        if (3 != GLidx && is_alt1_uniallelic) {
            vcfref = (regionpos > 0 ? refstring.substr(regionpos-1, 1) : "n");
            if (indelstring1.size() == 0 || indelstring1[0] == '<') { 
                vcfalt = SYMBOL_TO_DESC_ARR[s1]; 
            } else {
                vcfalt = vcfref;
                if (isSymbolIns(s1)) {
                    vcfalt += indelstring1;
                } else if (isSymbolDel(s1)) {
                    vcfref += indelstring1;
                } else {
                    vcfalt = SYMBOL_TO_DESC_ARR[s1];
                }
            }
        } else {
            AlignmentSymbol s2 = END_ALIGNMENT_SYMBOLS;
            std::string indelstring2;
            if (!is_alt1_uniallelic) {
                s2 = ref_alt1_alt2_alt3[1].first;
                indelstring2 = indelstrings1[1].second;
            } else {
                s2 = ref_alt1_alt2_alt3[2].first;
                std::vector<std::pair<std::array<uvc1_readnum_t, 2>, std::string>> indelstrings2 = indel_get_majority(
                        *(ref_alt1_alt2_alt3[2].second), 
                        tname, 
                        refpos, 
                        ref_alt1_alt2_alt3[2].first, 
                        is_rescued, 
                        0);
                indelstring2 = indelstrings2[0].second;
                if (s2 == s1) {
                    if (indelstrings2.size() < 2) {
                        fprintf(stderr, "Runtime error: indel at %u is invalid!\n", refpos);
                        abort();
                    }
                    indelstring2 = indelstrings2[1].second;
                }
            }
            auto vcfref1 = (regionpos > 0 ? refstring.substr(regionpos-1, 1) : "n");
            auto vcfalt1 = vcfref1;
            vcfref = vcfref1;
            vcfalt = vcfalt1;
            if (indelstring1.size() == 0 || indelstring1[0] == '<' || indelstring2.size() == 0 || indelstring2[0] == '<') {
                vcfalt = std::string(SYMBOL_TO_DESC_ARR[s1]) + "," + SYMBOL_TO_DESC_ARR[s2]; 
            } else {
                if (isSymbolIns(s1) && isSymbolIns(s2)) {
                    vcfalt = vcfref1 + indelstring1 + "," + vcfref1 + indelstring2;
                } else if (isSymbolDel(s1) && isSymbolDel(s2)) {
                    assert(indelstring1.size() != indelstring2.size());
                    assert(indelstring1.substr(0, MIN(indelstring1.size(), indelstring2.size())) 
                        == indelstring2.substr(0, MIN(indelstring1.size(), indelstring2.size())));
                    if (indelstring1.size() > indelstring2.size()) {
                        vcfref = vcfref1 + indelstring1;
                        vcfalt = vcfalt1 + "," + vcfalt1 + indelstring1.substr(indelstring2.size());
                    } else {
                        vcfref = vcfref1 + indelstring2;
                        vcfalt = vcfalt1 + indelstring2.substr(indelstring1.size()) + "," + vcfalt1;
                    }
                } else if (isSymbolIns(s1) && isSymbolDel(s2)) {
                    vcfalt = vcfref1 + indelstring1 + indelstring2 + "," + vcfref1;
                    vcfref = vcfref1 + indelstring2;
                } else if (isSymbolDel(s1) && isSymbolIns(s2)) {
                    vcfalt = vcfref1 + "," + vcfref1 + indelstring2 + indelstring1;
                    vcfref = vcfref1 + indelstring1;
                } else {
                    vcfalt = std::string(SYMBOL_TO_DESC_ARR[s1]) + "," + SYMBOL_TO_DESC_ARR[s2];  
                }
            }
        }
    }
    auto vcfpos = refpos + (isSymbolSubstitution(refsymbol) ? 1 : 0);
    std::string germ_GT = GTidx2GT[(is_alt1_uniallelic ? GLidx : 3)];
    uvc1_qual_t germ_GQ =  (is_alt1_uniallelic ? (GL4[0].second - GL4[1].second) : MIN(alt1_uniallelic_phred, GL4[0].second - GL4[1].second));
    std::vector<uvc1_readnum_t> germ_ADR;
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[0].second->cDP0a, 1, 0));
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[1].second->cDP0a, 1, 0));
    if (3 == GLidx) {
        germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[2].second->cDP0a, 1, 0));
    }
    
    std::string bcfline = string_join(std::array<std::string, 10>
    {{
        std::string(tname), 
        std::to_string(vcfpos), 
        std::string("."), 
        vcfref, 
        vcfalt, 
        std::to_string(germ_GQ), 
        std::string("PASS"), 
        std::string("GERMLINE"),
        std::string("GT:GQ:HQ:FT:CDP1:cDP1:GL4:GST:note"),
        string_join(std::array<std::string, 9>
        {{
            germ_GT, 
            std::to_string(germ_GQ), 
            std::string("0,0"), 
            "PASS",
            other_join(std::array<uvc1_readnum_t, 2> {{
                (symbol_format_vec[0].second->CDP12f[0] + symbol_format_vec[0].second->CDP12r[0]), 
                (symbol_format_vec[0].second->CDP12f[1] + symbol_format_vec[0].second->CDP12r[1]) }},
                ","
            ),
            other_join(germ_ADR, std::string(",")), 
            other_join(std::array<uvc1_readnum_t, 4>
            {{ 
                GL4raw[0].second, GL4raw[1].second, GL4raw[2].second, GL4raw[3].second, 
            }}, ","),
            other_join(std::array<uvc1_readnum_t, 8>
            {{
                a0LODQ, a1LODQ, a2LODQ, a3LODQ,
                a0a1LODQ, a1a0LODQ, a1a2LODQ, a2a1LODQ
            }}, ","),
            ref_alt1_alt2_alt3[0].second->note
        }}, ":")
    }}, "\t") + "\n";
    out_string += bcfline;
    return std::make_tuple(ret, fmtptr1, fmtptr2);
}

#include "version.h"
std::string 
generate_vcf_header(
        int argc,
        const char *const *argv,
        int32_t n_targets,
        const char *const *target_name,
        const uint32_t *target_len,
        const char *const tumor_sampleName,
        const CommandLineArgs & paramset) {
    time_t rawtime;
    time(&rawtime);
    char timestring[80];
    strftime(timestring, 80, "%F %T", localtime(&rawtime));
    std::string ret = "";
    ret += std::string("") + "##fileformat=VCFv4.2" + "\n" ;
    ret += std::string("") + "##fileDate=" + timestring + "\n" ;
    ret += std::string("") + "##variantCallerVersion=" + VERSION_DETAIL + "\n";
    ret += std::string("") + "##variantCallerCommand=";
    for (int i = 0; i < argc; i++) {
        ret += std::string("") + std::string(argv[i]) + "  ";
    }
    ret += "\n";
    ret += std::string("") + "##variantCallerInferredParameters=<" 
            + "inferred_sequencing_platform=" +  SEQUENCING_PLATFORM_TO_NAME.at(paramset.inferred_sequencing_platform)
            + ",central_readlen=" + std::to_string(paramset.central_readlen) 
            + ">\n";
    ret += std::string("") + "##reference=" + paramset.fasta_ref_fname + "\n";
    for (int i = 0; i < n_targets; i++) {
        ret += std::string("") + "##contig=<ID=" + target_name[i] + ",length=" + std::to_string(target_len[i]) + ">\n";
    }
    ret += std::string("") + "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location, where POS (start position) is one-based inclusive. "
        + "CAVEAT: this VCF line record is similar to a GVCF block but does not conform to the GVCF specifications. \">\n";
    
    for (size_t i = 0; i < bcfrec::FILTER_NUM; i++) {
        ret += std::string("") + bcfrec::FILTER_LINES[i] + "\n";
    }
    
    ret += "##INFO=<ID=ANY_VAR,Number=0,Type=Flag,Description=\"Any type of variant which may be caused by germline polymorphism and/or somatic mutation\">\n";
    ret += "##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"germline variant\">\n";
    ret += "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic variant\">\n";
    ret += "##INFO=<ID=MGVCF_BLOCK,Number=0,Type=Flag,Description=\"Multi-sample GVCF-like genomic regions consisting of " + std::to_string(MGVCF_REGION_MAX_SIZE) + " consecutive positions. " 
        + "MGVCF is modified from GVCF to allow for easy comparison of sequencing depths of multiple samples at any arbitrary position. "
        + "More detail is described in FORMAT/POS_VT_BDP_CDP_HomRefQ. \">\n";
    ret += "##INFO=<ID=ADDITIONAL_INDEL_CANDIDATE,Number=0,Type=Flag,Description=\"Position with an abnormally high number of (soft/hard)-clipped sequences adjacent to this position (which can be caused by long InDel, copy-number variation (CNV), structural variation (SV), etc.) or with a high STR track length after it\">\n";

    ret += ("##INFO=<ID=SomaticQ,Number=A,Type=Float,Description=\"Somatic quality of the variant, the PHRED-scale probability that this variant is not somatic. "
          "CAVEAT: if only tumor bam file is provided, then this quality usually cannot reach 60 even with the help of a very big germline database because "
          "germline and somatic variants share similar characteristics in the tumor. "
          "Therefore, a matched normal is absolutely required to confidently determine the germline-vs-somatic origin of a biological variant. \">\n");
    ret += "##INFO=<ID=TLODQ,Number=A,Type=Float,Description=\"Tumor log-of-data-likelihood quality, the PHRED-scale probability that this variant is not of biological origin (i.e., artifactual). \">\n";
    ret += "##INFO=<ID=NLODQ,Number=A,Type=Float,Description=\"Normal log-of-data-likelihood quality, the PHRED-scale probability that this variant is of germline origin. \">\n";
    ret += "##INFO=<ID=NLODV,Number=A,Type=String,Description=\"The variant symbol that minimizes NLODQ. \">\n";

    ret += "##INFO=<ID=TNBQF,Number=4,Type=Float,Description=\"Binomial reward, power-law reward, systematic-error penalty, and normal-adjusted tumor variant quality computed using deduplicated read fragments. \">\n";
    ret += "##INFO=<ID=TNCQF,Number=4,Type=Float,Description=\"Binomial reward, power-law reward, systematic-error penalty, and normal-adjusted tumor variant quality computed using consensus families of read fragments. \">\n";
    
    ret += "##INFO=<ID=tDP,Number=1,Type=Integer,Description=\"Tumor total deduped depth (deprecated, please see CDP1f and CDP1r). \">\n";
    ret += "##INFO=<ID=tADR,Number=R,Type=Integer,Description=\"Tumor deduped depth of each allele (deprecated, please see cDP1f and cDP1r). \">\n";
    ret += "##INFO=<ID=nDP,Number=1,Type=Integer,Description=\"Normal total deduped depth (deprecated, please see CDP1f and CDP1r). \">\n";
    ret += "##INFO=<ID=nADR,Number=R,Type=Integer,Description=\"Normal deduped depth of each allele (deprecated, please see cDP1f and cDP1r). \">\n";
    ret += "##INFO=<ID=tDPC,Number=1,Type=Integer,Description=\"Tumor total UMI-barcoded-family depth (deprecated, please see CDP2f and CDP2r). \">\n";
    ret += "##INFO=<ID=tADCR,Number=R,Type=Integer,Description=\"Tumor UMI-barcoded-family depth of each allele (deprecated, please see cDP2f and cDP2r). \">\n";
    ret += "##INFO=<ID=tbDP,Number=1,Type=Integer,Description=\"Tumor total non-deduped depth (deprecated, please see BDPf and BDPr). \">\n";
    
    ret += "##INFO=<ID=RU,Number=1,Type=String,Description=\"The shortest repeating unit in the reference\">\n";
    ret += "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"The number of non-interrupted RUs in the reference\">\n";
    ret += "##INFO=<ID=R3X2,Number=6,Type=Integer,Description=\"Repeat start position, repeat track length, and repeat unit size at the two positions before and after this VCF position. \">\n"; 
    
    for (size_t i = 0; i < bcfrec::FORMAT_NUM; i++) {
        ret += std::string("") + bcfrec::FORMAT_LINES[i] + "\n";
    }
    
    ret += std::string("") + "##FORMAT=<ID=GL4,Number=4,Type=Integer,Description=\"The four genotype likelihoods for 0/0, 0/1, 1/1, and 1/2\">\n";
    ret += std::string("") + "##FORMAT=<ID=GST,Number=.,Type=Integer,Description=\"The genotype statistics\">\n";
    ret += std::string("") + "##FORMAT=<ID=CDP1,Number=2,Type=Integer,Description=\"CDP1f + CDP1r\">\n";
    ret += std::string("") + "##FORMAT=<ID=cDP1,Number=2,Type=Integer,Description=\"cDP1f + cDP1r\">\n";
    
    ret += std::string("") + "##FORMAT=<ID=POS_VT_BDP_CDP_HomRefQ,Number=.,Type=Integer,Description=\"Summary of multiple GVCF regions in a line with INFO/MGVCF. "
            "This field conforms to the following regular expression: ((<pos>,<postype>,<.>,<dup>,<dedup>,<dedupBQ>,<homrefQ>,<.>)+<endpos>) "
            "where (x)+ means one or more occurrence of the expression x. "
            "The integer <pos> denotes position (coordinate on the reference sequence) that separates adjacent regions on the reference sequence. "
            "The integer <postype> denotes position type, where 1 and 2 mean SNV and InDel sub-positions, respectively. "
            "The missing integer represented by the dot symbol <.> is a sentinel value that delimits region separators (aka positions) and region information. "
            "The integer <dup> is the minimum non-deduplicated fragment depth of the region. "
            "The integer <dedup> is the minimum deduplicated fragment depth (with duplicated fragments counted only once). "
            "The integer <dedupBQ> is similar to <dedup> but is computed using only support with R1R2-adjusted BQ passing the threshold set by the command-line parameter --fam-thres-highBQ. "
            "The integer <homrefQ> is the minimum likelihood of the homozygous-reference (homref) genotype (GT) in this region. "
            "The integer <endpos> denotes the SNV ending sub-position of the set of regions on this VCF line, and <endpos> is the last number in this field. "
            "The (inclusive) begin position of the current region is the (exclusive) end position of the previous region. "
            "Each genomic position (e.g., chr1:99) is divided into (a) one SNV sub-position and (b) one InDel sub-position that is right after the SNV sub-position. "
            "The SNV prior of homref GT is used here. "
            "Thus, the actual InDel likelihood of homref GT is the one shown here plus "
            + std::to_string(paramset.germ_phred_hetero_indel - paramset.germ_phred_hetero_snp) + ". " +    
            "CAVEAT: HomRefQ is computed by a very fast but imprecise algorithm, so it is not as accurate as GQ. \">\n";
    
    ret += std::string("") + "##FORMAT=<ID=clipDP,Number=2,Type=Integer,Description=\"Total segment depth and segment depth with adjacent long clips "
            "(for the " + SYMBOL_TO_DESC_ARR[ADDITIONAL_INDEL_CANDIDATE_SYMBOL] + " symbolic ALT allele indicating that this position has a lot of long (soft/hard) clips nearby) or that this position is at the beginning of a long STR track\">\n";

    ret += std::string("") + "##phasing=partial\n";
    ret += std::string("") + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" 
            + paramset.sample_name + ((tumor_sampleName != NULL && paramset.is_tumor_format_retrieved) ? (std::string("\t") + tumor_sampleName) : std::string("")) + "\n";
    return ret;
}

std::string
bcf1_to_string_2(const bcf_hdr_t *tki_bcf1_hdr, const bcf1_t *bcf1_record) {
    kstring_t ks = { 0, 0, NULL };
    vcf_format(tki_bcf1_hdr, bcf1_record, &ks);
    assert (ks.l > 2);
    std::string ret = ks.s;
    if (ks.s != NULL) {
        free(ks.s);
    }
    return ret;
}

std::string
bcf1_to_string(const bcf_hdr_t *tki_bcf1_hdr, const bcf1_t *bcf1_record) {
    kstring_t ks = { 0, 0, NULL };
    vcf_format(tki_bcf1_hdr, bcf1_record, &ks);
    assert (ks.l > 2); 
    size_t idx = ks.l - 1;
    for (;idx != 0 && ks.s[idx] != '\t'; idx--) {
    }
    std::string ret = std::string(ks.s, idx, ks.l - idx - 1);
    if (ks.s != NULL) {
        free(ks.s);
    }
    return ret;
}

int 
fill_tki(auto & tki, const auto & fmt, size_t a = 1) {
    tki.BDP = fmt.BDPf.at(0) +fmt.BDPr.at(0);
    tki.bDP = fmt.bDPf.at(a) +fmt.bDPr.at(a);
    
    tki.CDP1x = fmt.CDP1x.at(0);
    tki.cDP1x = fmt.cDP1x.at(a);
    tki.cVQ1  = fmt.cVQ1.at(a);
    tki.cPCQ1 = fmt.cPCQ1.at(a);

    tki.CDP2x = fmt.CDP2x.at(0);
    tki.cDP2x = fmt.cDP2x.at(a);
    tki.cVQ2  = fmt.cVQ2.at(a);
    tki.cPCQ2 = fmt.cPCQ2.at(a);
    
    tki.bNMQ = fmt.bNMQ.at(a);
    tki.vHGQ = fmt.vHGQ;
    return 0;
};

template <bool TIsFmtTumor>
int 
fill_conditional_tki(auto & tki, const auto & fmt) {
    // extra code for backward compatibility
    if (TIsFmtTumor) {
        tki.tDP = (fmt.CDP1f[0] + fmt.CDP1r[0]);
        tki.tADR = {{ fmt.cDP1f[0] + fmt.cDP1r[0], LAST(fmt.cDP1f) + LAST(fmt.cDP1r) }};
        tki.nDP = 0;
        tki.nADR = {{ 0 }};
        tki.tDPC = (fmt.CDP2f[0] + fmt.CDP2r[0]);
        tki.tADCR = {{ fmt.cDP2f[0] + fmt.cDP2r[0], LAST(fmt.cDP2f) + LAST(fmt.cDP2r) }};
    } else {
        tki.nDP = (fmt.CDP1f[0] + fmt.CDP1r[0]);
        tki.nADR = {{ fmt.cDP1f[0] + fmt.cDP1r[0], LAST(fmt.cDP1f) + LAST(fmt.cDP1r) }};
    }
    return TIsFmtTumor;
}

const auto
calc_binom_powlaw_syserr_normv_quals(
        auto tAD, 
        auto tDP, 
        uvc1_qual_t tVQ, 
        uvc1_qual_t tnVQcap,
        auto nAD, 
        auto nDP, 
        uvc1_qual_t nVQ, 
        const double penal_dimret_coef,
        const uvc1_qual_t prior_phred,
        const uvc1_qual_t tn_dec_by_xm,
        const double powlaw_exponent,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    uvc1_qual_t binom_b10log10like = calc_binom_10log10_likeratio((tDP - tAD) / (tDP), nDP - nAD, nAD);
    const auto nADplus = nAD * BETWEEN(nDP / tDP - 1.0, 0, 1);
    double bjpfrac = ((tAD + 0.5) / (tDP + 1.0)) / ((nAD + 0.5 + nADplus) / (nDP + 1.0 + nADplus));
    uvc1_qual_t powlaw_b10log10like = round(powlaw_exponent * numstates2phred(bjpfrac));
    
    uvc1_qual_t tnVQinc = MAX3(-prior_phred, (-(uvc1_qual_t)nAD)*3, MIN(binom_b10log10like - prior_phred, powlaw_b10log10like - prior_phred));
    uvc1_qual_t tnVQdec = MAX(0, nVQ - MAX(0, MIN(
            binom_b10log10like - prior_phred, 
            (uvc1_qual_t)(mathsquare(log(MAX(bjpfrac, 1.001)) / log(2)) * penal_dimret_coef))));
    UPDATE_MAX(tnVQdec, MIN(nVQ + 9, tn_dec_by_xm));
    
    uvc1_qual_t tnVQ = MIN(tnVQcap, tVQ + tnVQinc) - tnVQdec;
    return std::array<uvc1_qual_t, 4> {{binom_b10log10like, powlaw_b10log10like, tnVQdec, tnVQ }};
};

const auto
calc_binom_powlaw_syserr_normv_quals2(
        auto tAD, 
        auto tDP, 
        uvc1_qual_t tVQ, 
        uvc1_qual_t tnVQcap,
        auto nAD, 
        auto nDP, 
        uvc1_qual_t nVQ,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    uvc1_qual_t binom_b10log10like = calc_binom_10log10_likeratio((tDP - tAD) / (tDP), nDP - nAD, nAD);
    uvc1_qual_t powlaw_b10log10like = (nAD <= 3 ? binom_b10log10like : (uvc1_qual_t)round(binom_b10log10like * 3 / nAD));
    uvc1_qual_t tnVQ = BETWEEN(tVQ + MAX3(MIN(binom_b10log10like, powlaw_b10log10like) - TVN_MICRO_VQ_DELTA, - TVN_MICRO_VQ_DELTA * nAD, -TVN_MICRO_VQ_DELTA) - nVQ, 0, tnVQcap);
    return std::array<uvc1_qual_t, 4> {{binom_b10log10like, powlaw_b10log10like, nVQ, tnVQ }};
};

int
append_vcf_record(
        std::string & out_string,
        const char *tname,
        const uvc1_refgpos_t refpos,
        const uvc1_refgpos_t region_offset,
        
        const std::string & refstring,
        const std::vector<RegionalTandemRepeat> & region_repeatvec,
        const std::string & repeatunit,
        const uvc1_readpos_t repeatnum,
        
        const AlignmentSymbol refsymbol,
        const AlignmentSymbol symbol,
        const bcfrec::BcfFormat & fmt,
        auto & tki,

        const uvc1_qual_t nlodq1,
        const AlignmentSymbol argmin_nlodq_symbol,
        const bool should_output_ref_allele,
        const bcf_hdr_t *g_bcf_hdr,
        const auto & baq_offsetarr,
        
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    
    const std::string & indelstring = LAST(fmt.gapSa);
    
    const bcfrec::BcfFormat & nfm = (tki.ref_alt.size() > 0 ? fmt : FORMAT_UNCOV);
    if (tki.ref_alt.size() == 0) {
        fill_tki(tki, fmt);
        fill_conditional_tki<true>(tki, fmt);
    } else {
        fill_conditional_tki<false>(tki, fmt);
    }
    const auto regionpos = refpos - region_offset;
    uvc1_refgpos_t vcfpos;
    std::string vcfref, vcfalt;
    if (indelstring.size() > 0) {
        vcfpos = refpos; 
        vcfref = (regionpos > 0 ? refstring.substr(regionpos-1, 1) : "n");
        vcfalt = vcfref;
        if (indelstring.size() == 0) {
            vcfalt = SYMBOL_TO_DESC_ARR[symbol];
        } else if ('<' == indelstring[0]) {
            vcfalt = indelstring;
        } else if (isSymbolIns(symbol)) {
            vcfalt += indelstring;
        } else {
            vcfref += indelstring;
        }
    } else {
        if (isSymbolSubstitution(symbol)) {
            vcfpos = refpos+1;
            vcfref = refstring.substr(regionpos, 1);
        } else {
            vcfpos = refpos;
            vcfref = (regionpos > 0 ? refstring.substr(regionpos-1, 1) : "n");
        }
        vcfalt = SYMBOL_TO_DESC_ARR[symbol];
    }
    
    assert (2 == fmt.cDP1v.size());
    assert (2 == fmt.cDP2v.size());
    assert (2 == fmt.CDP1v.size());
    assert (2 == fmt.CDP2v.size());
    assert (2 == fmt.cVQ1.size());
    assert (2 == fmt.cVQ2.size());
    
    assert (tki.cDP1x <= tki.CDP1x || !fprintf(stderr, "%d <= %d failed for cDP1x tname %s pos %d aln-symbol %d\n", tki.cDP1x, tki.CDP1x, tname, refpos, symbol));
    assert (tki.cDP2x <= tki.CDP2x || !fprintf(stderr, "%d <= %d failed for cDP2x tname %s pos %d aln-symbol %d\n", tki.cDP2x, tki.CDP2x, tname, refpos, symbol));
    const auto & rtr1 =  region_repeatvec.at(MAX(refpos - region_offset, paramset.indel_adj_tracklen_dist) - paramset.indel_adj_tracklen_dist);
    const auto & rtr2 =  region_repeatvec.at(MIN(refpos - region_offset + paramset.indel_adj_tracklen_dist, UNSIGN2SIGN(region_repeatvec.size()) - paramset.indel_adj_tracklen_dist));
    const auto rtr1_tpos = ((0 == rtr1.tracklen) ? 0 : (region_offset + rtr1.begpos));
    const auto rtr2_tpos = ((0 == rtr2.tracklen) ? 0 : (region_offset + rtr2.begpos));

    const auto nfm_cDP1x = collectget(nfm.cDP1x, 1);
    const auto nfm_CDP1x = collectget(nfm.CDP1x, 0);
    const auto nfm_cDP2x = collectget(nfm.cDP2x, 1);
    const auto nfm_CDP2x = collectget(nfm.CDP2x, 0);
    
    uvc1_refgpos_t phred_het3al_chance_inc_snp = MAX(0, 2 * paramset.germ_phred_hetero_snp - paramset.germ_phred_het3al_snp - TIN_CONTAM_MICRO_VQ_DELTA);
    uvc1_refgpos_t phred_het3al_chance_inc_indel = MAX(0, 2 * paramset.germ_phred_hetero_indel - paramset.germ_phred_het3al_indel - TIN_CONTAM_MICRO_VQ_DELTA);
    
    const auto tn_dec_by_xm = BETWEEN(MIN(LAST(fmt.bNMQ), tki.bNMQ),
        paramset.microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual_min, 
        paramset.microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual_max) 
        - paramset.microadjust_syserr_MQ_NMR_tn_syserr_no_penal_qual_min; // guaranteed least penalty
    
    double nfm_cDP1x_1add = 0;
    double nfm_cDP2x_1add = 0;
    if (tki.ref_alt.size() > 0 && does_fmt_imply_short_frag(fmt)) {
        // heuristics: if the assay is not WGS, then hybrid-capture or PCR-amplicon is used, so there is greater variation in ALT read counts
        nfm_cDP1x_1add = 0.5 * nfm_cDP1x / 100.0;
        nfm_cDP2x_1add = 0.5 * nfm_cDP2x / 100.0;
    }
    
    auto tn_dec_both_tlodq_nlodq = 0;
    if (tki.ref_alt.size() > 0 && tki.tDP > 500 && tki.nDP > 500 && (isSymbolDel(symbol)) && nfm.APDP[2] * 3 > nfm.APDP[0]) {
        tn_dec_both_tlodq_nlodq = MIN(non_neg_minus(collectget(nfm.cVQ1, 1), 31), 9);
    }
    const uvc1_qual_t prior_phred = ((SEQUENCING_PLATFORM_IONTORRENT == paramset.inferred_sequencing_platform) 
            ? (3+8) : (3));
    const std::array<uvc1_qual_t, 4> b_binom_powlaw_syserr_normv_q4filter = (paramset.tn_syserr_norm_devqual >= 0 
        ? calc_binom_powlaw_syserr_normv_quals(
            (tki.cDP1x + 0.5) / 100.0 + 0.0, 
            (tki.CDP1x + 1.0) / 100.0 + 0.0, 
            (tki.cVQ1),
            (tki.cPCQ1),
            (nfm_cDP1x + 0.5) / 100.0 + 0.0 + nfm_cDP1x_1add,
            (nfm_CDP1x + 1.0) / 100.0 + 0.0 + nfm_cDP1x_1add,
            non_neg_minus(collectget(nfm.cVQ1, 1), (isSymbolSubstitution(symbol) ? phred_het3al_chance_inc_snp : phred_het3al_chance_inc_indel)),
            paramset.tn_syserr_norm_devqual,
            prior_phred,
            tn_dec_by_xm,
            paramset.powlaw_exponent,
            0)
        : calc_binom_powlaw_syserr_normv_quals2(
            (tki.cDP1x + 0.5) / 100.0 + 0.0, 
            (tki.CDP1x + 1.0) / 100.0 + 0.0, 
            (tki.cVQ1),
            (tki.cPCQ1),
            (nfm_cDP1x + 0.5) / 100.0 + 0.0 + nfm_cDP1x_1add,
            (nfm_CDP1x + 1.0) / 100.0 + 0.0 + nfm_cDP1x_1add, 
            non_neg_minus(collectget(nfm.cVQ1, 1), (isSymbolSubstitution(symbol) ? phred_het3al_chance_inc_snp : phred_het3al_chance_inc_indel)),
            0));
    
    const auto c_binom_powlaw_syserr_normv_q4 = (paramset.tn_syserr_norm_devqual >= 0 
        ?  calc_binom_powlaw_syserr_normv_quals(
            (tki.cDP2x + 0.5) / 100.0 + 0.0, 
            (tki.CDP2x + 1.0) / 100.0 + 0.0, 
            (tki.cVQ2),
            (tki.cPCQ2),
            (nfm_cDP2x + 0.0) / 100.0 + 0.5 + nfm_cDP2x_1add,
            (nfm_CDP2x + 0.0) / 100.0 + 1.0 + nfm_cDP2x_1add,
            non_neg_minus(collectget(nfm.cVQ2, 1), (isSymbolSubstitution(symbol) ? phred_het3al_chance_inc_snp : phred_het3al_chance_inc_indel)),
            paramset.tn_syserr_norm_devqual,
            prior_phred,
            tn_dec_by_xm,
            paramset.powlaw_exponent,
            0)
        : calc_binom_powlaw_syserr_normv_quals2(
            (tki.cDP2x + 0.5) / 100.0 + 0.0, 
            (tki.CDP2x + 1.0) / 100.0 + 0.0, 
            (tki.cVQ2),
            (tki.cPCQ2),
            (nfm_cDP2x + 0.0) / 100.0 + 0.5 + nfm_cDP2x_1add,
            (nfm_CDP2x + 0.0) / 100.0 + 1.0 + nfm_cDP2x_1add,
            non_neg_minus(collectget(nfm.cVQ2, 1), (isSymbolSubstitution(symbol) ? phred_het3al_chance_inc_snp : phred_het3al_chance_inc_indel)),
            0));
    
    uvc1_qual_t tlodq1 = MAX(b_binom_powlaw_syserr_normv_q4filter[3], c_binom_powlaw_syserr_normv_q4[3]);
    float lowestVAQ = MIN(tlodq1, 4) + 0.5 + (prob2realphred(1 / (double)(tki.bDP + 1)) * (tki.bDP) / (double)(tki.BDP + + tki.bDP + 1));
    
    const auto tlodq = tlodq1 - tn_dec_both_tlodq_nlodq;
    const auto nlodq = nlodq1 - tn_dec_both_tlodq_nlodq;
    uvc1_qual_t somaticq = MIN(tlodq, nlodq);
    float vcfqual = calc_non_negative((tki.ref_alt.size() > 0) ? ((float)somaticq) : MAX((float)tlodq, lowestVAQ));
    std::string infostring = std::string(tki.ref_alt.size() > 0 ? "SOMATIC" : "ANY_VAR");
    infostring += std::string(";SomaticQ=") + std::to_string(somaticq);
    infostring += std::string(";TLODQ=") + std::to_string(tlodq);
    infostring += std::string(";NLODQ=") + std::to_string(nlodq);
    infostring += std::string(";NLODV=") + SYMBOL_TO_DESC_ARR[argmin_nlodq_symbol];
    infostring += std::string(";TNBQF=") + other_join(b_binom_powlaw_syserr_normv_q4filter);
    infostring += std::string(";TNCQF=") + other_join(c_binom_powlaw_syserr_normv_q4);
    
    infostring += std::string(";tDP=") + std::to_string(tki.tDP);
    infostring += std::string(";tADR=") + other_join(tki.tADR, ",");
    infostring += std::string(";nDP=") + std::to_string(tki.nDP);
    infostring += std::string(";nADR=") + other_join(tki.nADR, ",");
    infostring += std::string(";tDPC=") + std::to_string(tki.tDPC);
    infostring += std::string(";tADCR=") + other_join(tki.tADCR, ",");
    infostring += std::string(";tbDP=") + std::to_string(tki.BDP);
    
    infostring += std::string(";RU=") + repeatunit + ";RC=" + std::to_string(repeatnum);
    // This can be useful for debugging only.
    if ((paramset.debug_note_flag & DEBUG_NOTE_FLAG_BITMASK_BAQ_OFFSETARR)) { 
        infostring += std::string(";RBAQ=") + std::to_string(baq_offsetarr.getByPos(refpos)); 
    }
    infostring += std::string(";R3X2=") + other_join(std::array<int32_t, 6>{{rtr1_tpos, rtr1.tracklen, rtr1.unitlen, rtr2_tpos, rtr2.tracklen, rtr2.unitlen}});
    
    std::string vcffilter = "";
    if (vcfqual < 10) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::Q10]) + ";");
    } else if (vcfqual < 20) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::Q20]) + ";");
    } else if (vcfqual < 30) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::Q30]) + ";");
    } else if (vcfqual < 40) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::Q40]) + ";");
    } else if (vcfqual < 50) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::Q50]) + ";");
    } else if (vcfqual < 60) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::Q60]) + ";");
    } else {
        vcffilter += (std::string("PASS") + ";");
    }
    vcffilter.pop_back();
    
    const auto vad1curr = LAST(fmt.aBQ2);
    const auto vdp1curr = fmt.ABQ2[0];
    const auto vad2curr = tki.bDP;
    const auto vdp2curr = tki.BDP;
    
    const bool keep_var = (((vcfqual >= paramset.vqual)
            || ((NOT_PROVIDED == paramset.vcf_tumor_fname) && 
                    ((vad1curr >= paramset.vad1 && vdp1curr >= paramset.vdp1 && (vdp1curr * paramset.vfa1) <= vad1curr)
                  || (vad2curr >= paramset.vad2 && vdp2curr >= paramset.vdp2 && (vdp2curr * paramset.vfa2) <= vad2curr))))
            && (symbol != refsymbol || (should_output_ref_allele)));
    const auto min_ad = ((symbol == refsymbol) ? paramset.min_r_ad : paramset.min_a_ad);
    if (keep_var && tki.bDP >= min_ad) {
        out_string += string_join(std::array<std::string, 9>{{
                std::string(tname), std::to_string(vcfpos), ".", vcfref, vcfalt, std::to_string(vcfqual), vcffilter, 
                infostring, bcfrec::FORMAT_STRING_PER_REC}}, "\t") + "\t";
        bcfrec::streamAppendBcfFormat(out_string, fmt);
        out_string += ((tki.ref_alt.size() > 0 && paramset.is_tumor_format_retrieved) ? bcf1_to_string(g_bcf_hdr, tki.bcf1_record) : std::string("")) + "\n";
    }
    return 0;
};

#endif
