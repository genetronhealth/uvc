#ifndef main_hpp_INCLUDED
#define main_hpp_INCLUDED

#include "bcf_formats.step1.c" // auto-generated

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

#define THE_BAQ_MAX 150 // 200
//#define INDEL_MUL_PER_BAQ 2.0
//#define SNV_MUL_PER_BAQ 3.0
// T=uint32_t for deletion and T=string for insertion
template <class T>
uint32_t
posToIndelToData_get(const std::map<uint32_t, std::map<T, uint32_t>> & pos2indel2data, const uint32_t refpos, const T & indel, uint32_t default1 = 0, uint32_t default2 = 0) {
    if (pos2indel2data.find(refpos) == pos2indel2data.end()) { return default1; }
    if (pos2indel2data.at(refpos).find(indel) == pos2indel2data.at(refpos).end()) { return default2; }
    return pos2indel2data.at(refpos).at(indel);
}

template <class T>
void
posToIndelToCount_inc(std::map<uint32_t, std::map<T, uint32_t>> & pos2indel2count, uint32_t pos, const T indel, uint32_t incvalue = 1) {
    assert(incvalue > 0);
    auto pos2indel2count4it = pos2indel2count.insert(std::make_pair(pos, std::map<T, uint32_t>()));
    auto indel2count4it = pos2indel2count4it.first->second.insert(std::make_pair(indel, 0));
    indel2count4it.first->second += incvalue;
    assert (posToIndelToData_get(pos2indel2count, pos, indel) > 0);
}

template <class T>
T
posToIndelToCount_updateByConsensus(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src, uint32_t epos, uint32_t incvalue = 1) {
    const auto pos2indel2count4it = src.find(epos);
    assert(pos2indel2count4it != src.end());

    auto indel2count = pos2indel2count4it->second;
    // if (indel2count.size() > 1) { return T(); }

    const T & src_indel = (indel2count.size() > 1 ? T() : indel2count.begin()->first);
    // const uint32_t src_count = indel2count.begin()->second;
    assert (indel2count.begin()->second > 0);
    posToIndelToCount_inc<T>(dst, epos, src_indel, incvalue);
    return src_indel;
    //return (int)src_count;
}

template <bool TIsIncVariable, class T>
int
posToIndelToCount_updateByRepresentative(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src, uint32_t epos, uint32_t incvalue = 1) {
    const auto pos2indel2count4it = src.find(epos);
    assert(pos2indel2count4it != src.end());

    T max_indel;
    uint32_t max_count = 0;
    uint32_t sum_count = 0;
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
    return (int)max_count;
}

template <class T>
void
posToIndelToCount_updateBySummation(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src) {
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
posToIndelToCount_DlenToDseq(std::map<uint32_t, std::map<std::string, uint32_t>> & dst, const std::map<uint32_t, std::map<uint32_t, uint32_t>> & src,
        const std::string & refchars, uint32_t incluBegPos) {
    for (auto src_pos2indel2count4it : src) {
        uint32_t src_pos = src_pos2indel2count4it.first;
        auto src_del2count = src_pos2indel2count4it.second;
        for (auto src_del2count4it : src_del2count) {
            std::string dseq = refchars.substr(src_pos - incluBegPos, src_del2count4it.first);
            posToIndelToCount_inc(dst, src_pos, dseq, src_del2count4it.second);
        }
    }
}

#define bam_phredi(b, i) (bam_get_qual((b))[(i)])

const bcfrec::BcfFormat FORMAT_UNCOV = bcfrec::BcfFormat();

enum ValueType {
    SYMBOL_COUNT_SUM,
    BASE_QUALITY_MAX,
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
insLenToSymbol(unsigned int len) {
    //assert(len > 0);
    return 1 == len ? LINK_I1 : (2 == len ? LINK_I2 : LINK_I3P);
}

constexpr AlignmentSymbol 
delLenToSymbol(unsigned int len) {
    //assert(len > 0);
    return 1 == len ? LINK_D1 : (2 == len ? LINK_D2 : LINK_D3P);
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


struct PhredMutationTable {
    const unsigned int transition_CG_TA;
    const unsigned int transition_TA_CG;
    const unsigned int transversion_any;
    const unsigned int indel_open;
    const unsigned int indel_ext;
    PhredMutationTable(unsigned int CG_TA, unsigned int TA_CG, unsigned int transversionany, unsigned int idopen, unsigned int idext)
            : transition_CG_TA(CG_TA), transition_TA_CG(TA_CG), transversion_any(transversionany), indel_open(idopen), indel_ext(idext) 
            {
    }
    const unsigned int toPhredErrRate(const AlignmentSymbol con_symbol, const AlignmentSymbol alt_symbol) const {
        if (isSymbolIns(con_symbol) || isSymbolDel(con_symbol)) {
            return indel_open;
        } else if (con_symbol == LINK_M) {
            if (LINK_D1 == alt_symbol || LINK_I1 == alt_symbol) {
                return indel_open; // + indel_ext * 1;
            } else if (LINK_D2== alt_symbol || LINK_I2 == alt_symbol) {
                return indel_open + indel_ext * 1;
            } else {
                return indel_open + indel_ext * 2;
            }
        } else if ((con_symbol == BASE_C && alt_symbol == BASE_T) || (con_symbol == BASE_G && alt_symbol == BASE_A)) {
            return transition_CG_TA;
        } else if ((con_symbol == BASE_T && alt_symbol == BASE_C) || (con_symbol == BASE_A && alt_symbol == BASE_G)) {
            return transition_TA_CG;
        }  else {
            return transversion_any;
        }
    }
};

const unsigned int SYMBOL_TO_INDEL_N_UNITS[] = {
    [BASE_A] = 0, [BASE_C] = 0, [BASE_G] = 0, [BASE_T] = 0, [BASE_N] = 0,
    [BASE_NN] = 0, 
    [LINK_M] = 0, 
    [LINK_D3P] = 3, [LINK_D2] = 2, [LINK_D1] = 1,
    [LINK_I3P] = 3, [LINK_I2] = 2, [LINK_I1] = 1,
    [LINK_NN] = 0,
};


const char* SYMBOL_TO_DESC_ARR[] = {
    [BASE_A] = "A", [BASE_C] = "C", [BASE_G] = "G", [BASE_T] = "T", [BASE_N] = "N",
    [BASE_NN] = "*", 
    [LINK_M] = "<LR>", 
    [LINK_D3P] = "<LD3P>", [LINK_D2] = "<LD2>", [LINK_D1] = "<LD1>",
    [LINK_I3P] = "<LI3P>", [LINK_I2] = "<LI2>", [LINK_I1] = "<LI1>",
    [LINK_NN] = "*",
    [END_ALIGNMENT_SYMBOLS] = "<ALN_SYMBOL_END>",
};

const std::map<std::string, AlignmentSymbol>
_generateDescToSymbolMap() {
    std::map<std::string, AlignmentSymbol> ret;
    for (AlignmentSymbol s = AlignmentSymbol(0); s < END_ALIGNMENT_SYMBOLS; s = AlignmentSymbol(1+(unsigned int)s)) {
        ret[SYMBOL_TO_DESC_ARR[s]] = s;
    }
    return ret;
}

const std::map<std::string, AlignmentSymbol> DESC_TO_SYMBOL_MAP = _generateDescToSymbolMap();

template <class T>
class TDistribution {
protected:
    std::array<T, NUM_ALIGNMENT_SYMBOLS> symbol2data; // [NUM_ALIGNMENT_SYMBOLS];
};

typedef uint32_t molcount_t;

template <class TB2C>
class GenericSymbol2Bucket2Count : TDistribution<TB2C> {
public:    
    molcount_t
    getSymbolBucketCount(AlignmentSymbol symbol, unsigned int bucket) const {
        // assert(bucket < NUM_BUCKETS);
        return this->symbol2data[symbol][bucket];
    };
    
    const TB2C &
    getSymbolCounts(AlignmentSymbol symbol) const {
        return this->symbol2data[symbol];
    };

    void
    incSymbolBucketCount(AlignmentSymbol symbol, unsigned int bucket, unsigned int increment) {
        // assert(bucket < NUM_BUCKETS || !fprintf(stderr, "%d < %d failed!", bucket, NUM_BUCKETS));
        this->symbol2data[symbol].at(bucket) += increment;
    };
    
    const TB2C
    vectorsumBySymbolType(const SymbolType symbolType) const {
        TB2C ret = {0};
        for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; 
                symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; 
                symbol = AlignmentSymbol(1+((unsigned int)symbol))) {
            for (size_t i = 0; i < this->symbol2data[symbol].size(); i++) {
                ret[i] += this->getSymbolBucketCount(symbol, i);
            }
        }
        return ret;
    };

    void
    clearSymbolBucketCount() { // AlignmentSymbol symbol, unsigned int bucket, unsigned int increment) {
        //char errmsg[100];
        //sprintf(errmsg, "%d * %d != %d\n", sizeof(TB2C), NUM_BUCKETS, sizeof(symbol2data));
        assert(sizeof(TB2C) * NUM_ALIGNMENT_SYMBOLS == sizeof(this->symbol2data) || !fprintf(stderr, "%lu * %u != %lu\n", sizeof(TB2C), NUM_ALIGNMENT_SYMBOLS, sizeof(this->symbol2data)));
        // memset(this->symbol2data, 0, sizeof(this->symbol2data));
        
        for (unsigned int i = 0; i < NUM_ALIGNMENT_SYMBOLS; i++) {
            for (unsigned int j = 0; j < this->symbol2data[0].size(); j++) {
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
    incSymbolCount(const AlignmentSymbol symbol, const TInteger increment, const unsigned int update_max_inc = 0) {
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
        for (AlignmentSymbol symb = beg; symb <= end; symb = AlignmentSymbol(((unsigned int)symb) + 1)) {
            alpha_sum += this->symbol2data[symb];
        }
        return alpha_sum;
    };
    
    const TInteger
    sumBySymbolType(const SymbolType symbolType) const {
        if (symbolType == BASE_SYMBOL) {
            return this->_sumBySymbolType(BASE_A, BASE_NN);
        } else if (symbolType == LINK_SYMBOL) {
            return this->_sumBySymbolType(LINK_M, LINK_NN);
        } else {
            abort();
            return -1;
        }
    };
    /*
    const TInteger
    sumBySymbolType2(const SymbolType symbolType) const {
        return sumBySymbolType(symbolType);
        if (symbolType == BASE_SYMBOL) {
            return this->_sumBySymbolType(BASE_A, AlignmentSymbol(((unsigned int)BASE_NN) - 1));
        } else if (symbolType == LINK_SYMBOL) {
            return this->_sumBySymbolType(LINK_M, AlignmentSymbol(((unsigned int)LINK_NN) - 1));
        } else {
            abort();
            return -1;
        }
    };
    
    const TInteger
    _maxBySymbolType(AlignmentSymbol beg, AlignmentSymbol end) const {
        assert (beg <= end);
        
        TInteger alpha_sum = 0;
        for (AlignmentSymbol symb = beg; symb <= end; symb = AlignmentSymbol(((unsigned int)symb) + 1)) {
            alpha_sum = MAX(alpha_sum, this->symbol2data[symb]);
        }
        return alpha_sum;
    };
    
    const TInteger
    maxBySymbolType(const SymbolType symbolType) const {
        if (symbolType == BASE_SYMBOL) {
            return this->_maxBySymbolType(BASE_A, BASE_NN);
        } else if (symbolType == LINK_SYMBOL) {
            return this->_maxBySymbolType(LINK_M, LINK_NN);
        } else {
            abort();
            return -1;
        }
    };
    */

    template <bool TIndelIsMajor>
    int
    _fillConsensusCounts(
            AlignmentSymbol & count_argmax, unsigned int & count_max, unsigned int & count_sum,
            AlignmentSymbol incluBeg, AlignmentSymbol incluEnd) const {
        assert (incluBeg <= incluEnd);
        
        count_argmax = incluEnd; // AlignmentSymbol(NUM_ALIGNMENT_SYMBOLS); // TODO:investigate further
        count_max = 0;
        count_sum = 0;
        for (AlignmentSymbol symb = incluBeg; symb <= incluEnd; symb = AlignmentSymbol(((unsigned int)symb) + 1)) {
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
    const int
    fillConsensusCounts(
            AlignmentSymbol & count_argmax, unsigned int & count_max, unsigned int & count_sum,
            const SymbolType symbolType) const {
        if (symbolType == BASE_SYMBOL) {
            return this->template _fillConsensusCounts<false        >(count_argmax, count_max, count_sum, BASE_A, BASE_NN);
        } else if (symbolType == LINK_SYMBOL) {
            return this->template _fillConsensusCounts<TIndelIsMajor>(count_argmax, count_max, count_sum, LINK_M, LINK_NN);
        } else {
            abort();
            return -1;
        }
    };
    
    template</*ValueType T_SymbolCountType,*/ bool TIndelIsMajor>
    const AlignmentSymbol
    _updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count,
            const SymbolType symbolType, const AlignmentSymbol ambig_pos, unsigned int incvalue2) {
        AlignmentSymbol argmax_count = END_ALIGNMENT_SYMBOLS; // AlignmentSymbol(0) is not fully correct
        unsigned int max_count = 0;
        unsigned int sum_count = 0;
        thatSymbol2Count.template fillConsensusCounts<TIndelIsMajor>(argmax_count, max_count, sum_count, symbolType);
        unsigned int incvalue;
        if (1 /*T_SymbolCountType == SYMBOL_COUNT_SUM*/) {
            incvalue = incvalue2;
        //} else if (T_SymbolCountType == BASE_QUALITY_MAX) {
            //assert(max_count < 96);
            //incvalue = THE_PHRED_TO_ERROR_PROBABILITY.over2pow16[max_count];
        //
        } else {
            assert(false);
        }
        
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

    template</*ValueType T_SymbolCountType, */
        bool TIndelIsMajor = false>
    std::array<AlignmentSymbol, 2>
    updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count, unsigned int incvalue = 1) {
        AlignmentSymbol baseSymb = this->template _updateByConsensus</*T_SymbolCountType,*/ false        >(thatSymbol2Count, BASE_SYMBOL, BASE_NN, incvalue);
        AlignmentSymbol linkSymb = this->template _updateByConsensus</*T_SymbolCountType,*/ TIndelIsMajor>(thatSymbol2Count, LINK_SYMBOL, LINK_NN, incvalue);
        return {baseSymb, linkSymb};
    };
    
    template <bool TIsIncVariable = true>
    AlignmentSymbol
    updateByRepresentative(const GenericSymbol2Count<TInteger> & other, unsigned int incvalue = 1) {
        AlignmentSymbol consalpha; 
        unsigned int countalpha, totalalpha;
        for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1 + (unsigned int)symbolType)) { 
            other.fillConsensusCounts(consalpha, countalpha, totalalpha, symbolType);
            if (countalpha > 0) {
                this->symbol2data[consalpha] += (TIsIncVariable ? totalalpha : incvalue);
            }
        }
        return consalpha;
    };
    /*
    int
    updateByFiltering(
            std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> & con_symbols, 
            const GenericSymbol2Count<TInteger> & other, 
            const unsigned int thres = 0, // const GenericSymbol2Count<TInteger> & thres, 
            const unsigned int incvalue = 1) {
        int ret = 0;
        AlignmentSymbol consalpha;
        unsigned int countalpha, totalalpha;
        for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1 + (unsigned int)symbolType)) {
            if (LINK_SYMBOL == symbolType) {
                other.template fillConsensusCounts<true >(consalpha, countalpha, totalalpha, symbolType);
            } else {
                other.template fillConsensusCounts<false>(consalpha, countalpha, totalalpha, symbolType);
            }
            auto adjcount = MAX(countalpha * 2, totalalpha) - totalalpha;
            if (adjcount >= thres && adjcount > 0) {
                this->symbol2data[consalpha] += incvalue;
                ret++;
            }
            con_symbols[symbolType] = consalpha;
        }
        return ret;
    };
    */
    int
    updateByMajorMinusMinor(
            std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> & con_symbols,
            std::array<unsigned int, NUM_SYMBOL_TYPES> & con_counts,
            const GenericSymbol2Count<TInteger> & other) {
        int ret = 0;
        AlignmentSymbol consalpha;
        unsigned int countalpha, totalalpha;
        for (SymbolType symbolType : SYMBOL_TYPE_ARR) {
            if (LINK_SYMBOL == symbolType) {
                other.template fillConsensusCounts<true >(consalpha, countalpha, totalalpha, symbolType);
            } else {
                other.template fillConsensusCounts<false>(consalpha, countalpha, totalalpha, symbolType);
            }
            auto adjcount = MAX(countalpha * 2, totalalpha) - totalalpha;
            if (adjcount > 0) {
                this->symbol2data[consalpha] += adjcount;
                ret++;
            }
            con_symbols[symbolType] = consalpha;
        }
        return ret;
    };
};

template<class T>
class CoveredRegion {
         
protected:
    
    std::vector<T> idx2symbol2data;
    std::array<std::map<uint32_t, std::map<uint32_t   , uint32_t>>, NUM_INS_SYMBOLS> pos2dlen2data;
    std::array<std::map<uint32_t, std::map<std::string, uint32_t>>, NUM_DEL_SYMBOLS> pos2iseq2data;
public:
    
    const uint32_t tid;
    const uint32_t incluBegPosition; // end_pos = incluBegPosition + idx2symbol2data

    CoveredRegion() {};    
    CoveredRegion(uint32_t tid, unsigned int beg, unsigned int end): tid(tid), incluBegPosition(beg)  {
        assert (beg < end || !fprintf(stderr, "assertion %d < %d failed!\n", beg, end));
        this->idx2symbol2data = std::vector<T>(end-beg); // end should be real end position plus one
    };
    
    T &
    getRefByPos(const unsigned int pos, const bam1_t *bam = NULL) {
        assert(pos >= this->incluBegPosition || !fprintf(stderr, "%u >= %u failed for qname %s !!\n", pos, this->incluBegPosition, (NULL != bam ? bam_get_qname(bam) : "?")));
        unsigned int pos2 = pos - this->incluBegPosition;
        assert(pos2 < idx2symbol2data.size() || !fprintf(stderr, "%u  < %lu failed for qname %s !!\n", pos, this->incluBegPosition + idx2symbol2data.size(), (NULL != bam ? bam_get_qname(bam) : "?")));
        return this->idx2symbol2data[pos2];
    };
    
    const T &
    getByPos(const unsigned int pos, const bam1_t *bam = NULL) const {
        assert(pos >= this->incluBegPosition || !fprintf(stderr, "%u >= %u failed for qname %s !\n", pos, this->incluBegPosition, (NULL != bam ? bam_get_qname(bam) : "?")));
        unsigned int pos2 = pos - this->incluBegPosition;
        assert(pos2 < idx2symbol2data.size() || !fprintf(stderr, "%u  < %lu failed for qname %s !\n", pos, this->incluBegPosition + idx2symbol2data.size(), (NULL != bam ? bam_get_qname(bam) : "?")));
        return this->idx2symbol2data[pos2];
    };
    
    const size_t 
    getIncluBegPosition() const {
        return this->incluBegPosition;
    };
    const size_t 
    getExcluEndPosition() const {
        return this->incluBegPosition + idx2symbol2data.size();
    };
    
    const std::map<uint32_t, std::map<uint32_t   , uint32_t>> & 
    getPosToDlenToData(const AlignmentSymbol s) const { 
        unsigned int idx = (LINK_D1 == s ? 0 : ((LINK_D2 == s) ? 1: 2));
        return pos2dlen2data[idx];
    };
    const std::map<uint32_t, std::map<std::string, uint32_t>> & 
    getPosToIseqToData(const AlignmentSymbol s) const {
        unsigned int idx = (LINK_I1 == s ? 0 : ((LINK_I2 == s) ? 1: 2));
        return pos2iseq2data[idx];
    };
    
    std::map<uint32_t, std::map<uint32_t   , uint32_t>> & 
    getRefPosToDlenToData(const AlignmentSymbol s) {
        unsigned int idx = (LINK_D1 == s ? 0 : ((LINK_D2 == s) ? 1: 2));
        return pos2dlen2data[idx];
    };
    std::map<uint32_t, std::map<std::string, uint32_t>> & 
    getRefPosToIseqToData(const AlignmentSymbol s) {
        unsigned int idx = (LINK_I1 == s ? 0 : ((LINK_I2 == s) ? 1: 2));
        return pos2iseq2data[idx];
    };
};

template <class TSymbol2Bucket2Count>
class GenericSymbol2Bucket2CountCoverage : public CoveredRegion<TSymbol2Bucket2Count> {
    public:
    GenericSymbol2Bucket2CountCoverage() : CoveredRegion<TSymbol2Bucket2Count>(0, 0, 1) { /*assert(false);*/ } // just to fix compiling error
    GenericSymbol2Bucket2CountCoverage(auto tid, auto beg, auto end) : CoveredRegion<TSymbol2Bucket2Count>(tid, beg, end) {}
};

typedef std::array<molcount_t, NUM_BUCKETS> Bucket2Count;

typedef GenericSymbol2Bucket2Count<Bucket2Count> Symbol2Bucket2Count;
//typedef GenericSymbol2Bucket2Count<Bucket2CountEdgeDist> Symbol2Bucket2CountEdgeDist;
//typedef GenericSymbol2Bucket2Count<Bucket2CountNumMisma> Symbol2Bucket2CountNumMisma;

typedef GenericSymbol2Count<uint32_t> Symbol2Count;
typedef GenericSymbol2Count<uint64_t> Symbol2CountUint64;

typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2Count> Symbol2Bucket2CountCoverage;
//typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2CountEdgeDist> Symbol2Bucket2CountCoverageEdgeDist;
//typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2CountNumMisma> Symbol2Bucket2CountCoverageNumMisma;

typedef CoveredRegion<std::array<molcount_t, NUM_SEG_FORMAT_PREP_SETS>> SegFormatPrepSets;
typedef CoveredRegion<std::array<molcount_t, NUM_SEG_FORMAT_THRES_SETS>> SegFormatThresSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_SEG_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2SegFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_FRAG_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2FragFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_FAM_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2FamFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_DUPLEX_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2DuplexFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_VQ_FORMAT_TAG_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2VQFormatTagSets;

// template <class TFormatYSet>
unsigned int
formatSumBySymbolType(
        const auto & xFormatYSets,
        const auto symbolType,
        const auto formatSet) {
    unsigned int ret = 0;
    for (auto symbol : SYMBOL_TYPE_TO_SYMBOLS[symbolType]) {
        ret += xFormatYSets[symbol][formatSet];
    }
    return ret;
};

void
initTidBegEnd(uint32_t & tid, uint32_t & inc_beg, uint32_t & exc_end) {
   tid = UINT32_MAX;
   inc_beg = UINT32_MAX;
   exc_end = 0;
};

int 
fillTidBegEndFromAlns1(uint32_t & tid, uint32_t & inc_beg, uint32_t & exc_end, const std::vector<bam1_t *> & alns1, bool initialized=false) {
    assert(alns1.size() > 0);
    if (!initialized) {
        initTidBegEnd(tid, inc_beg, exc_end);
    }
    for (const bam1_t *aln : alns1) {
        assert (tid == UINT32_MAX || SIGN2UNSIGN(aln->core.tid) == tid);
        tid = aln->core.tid;
        inc_beg = MIN(inc_beg, SIGN2UNSIGN(aln->core.pos));
        exc_end = MAX(exc_end, SIGN2UNSIGN(bam_endpos(aln))) + 1; // accounts for insertion at the end 
    }
    assert (tid != UINT32_MAX);
    assert (inc_beg < exc_end);
    return 0;
};

int 
fillTidBegEndFromAlns2(uint32_t  & tid, uint32_t & inc_beg, uint32_t & exc_end, const std::vector<std::vector<bam1_t *>> &alns2, bool initialized=false) {
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
fillTidBegEndFromAlns3(uint32_t & tid, uint32_t & inc_beg, uint32_t & exc_end, const std::vector<std::vector<std::vector<bam1_t *>>> & alns3, bool initialized=false) {
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
is_indel_context_more_STR(unsigned int rulen1, unsigned int rc1, unsigned int rulen2, unsigned int rc2) {
    if (rulen1 > 6 || rulen2 > 6) {
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

int 
indelpos_to_context(
        std::string & repeatunit, unsigned int & max_repeatnum,
        const std::string & refstring, unsigned int refpos) {
    max_repeatnum = 0;
    if (refpos >= refstring.size()) {
        repeatunit = "";
        return -1;
    }
    unsigned int repeatsize_at_max_repeatnum = 0;
    for (unsigned int repeatsize = 1; repeatsize <= 6; repeatsize++) {
        unsigned int qidx = refpos;
        while (qidx + repeatsize < refstring.size() && refstring[qidx] == refstring[qidx+repeatsize]) {
            qidx++;
        }
        unsigned int repeatnum = (qidx - refpos) / repeatsize + 1;
        if (is_indel_context_more_STR(repeatsize, repeatnum, repeatsize_at_max_repeatnum, max_repeatnum)) {
            max_repeatnum = repeatnum;
            repeatsize_at_max_repeatnum = repeatsize;
        }
    }
    repeatunit = refstring.substr(refpos, repeatsize_at_max_repeatnum);
    return 0;
}

unsigned int
indel_len_rusize_phred(unsigned int indel_len, unsigned int repeatunit_size) {
    assert (indel_len > 0 && repeatunit_size > 0);
    // python code: for i in range(1, 20): print("{}  {}".format(  int( round(10.0/log(10.0)*log(i)) ) , i  )) 
    const std::array<unsigned int, 18+1> n_units_to_phred = {{
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
        return n_units_to_phred[MIN(n_units   , n_units_to_phred.size()-1)];
    } else {
        return n_units_to_phred[MIN(indel_len, n_units_to_phred.size()-1)];
        // return 7;
    }
}

unsigned int
indel_phred(double ampfact, unsigned int cigar_oplen, unsigned int repeatsize_at_max_repeatnum, unsigned int max_repeatnum) {
    unsigned int region_size = repeatsize_at_max_repeatnum * max_repeatnum;
    double num_slips = (region_size > 64 ? (double)(region_size - 8) : log1p(exp((double)region_size - (double)8))) 
            * ampfact / ((double)(repeatsize_at_max_repeatnum * repeatsize_at_max_repeatnum)); //  / indel_n_units;
    return prob2phred((1.0 - DBL_EPSILON) / (num_slips + 1.0));
    // AC AC AC : repeatsize_at_max_repeatnum = 2, indel_n_units = 3
}

std::vector<RegionalTandemRepeat>
refstring2repeatvec(
        const auto & refstring, 
        unsigned int specialflag = 0) {
    std::vector<RegionalTandemRepeat> region_repeatvec(refstring.size());
    for (auto & rtr : region_repeatvec) {
        rtr.indelphred = 40;
    }
    for (int refpos = 0; refpos < refstring.size();) {
        unsigned int repeatsize_at_max_repeatnum = 0;
        unsigned int max_repeatnum = 0;
        unsigned int repeat_endpos = refpos;
        for (unsigned int repeatsize = 1; repeatsize <= 6; repeatsize++) {
            unsigned int qidx = refpos;
            while (qidx + repeatsize < refstring.size() && refstring[qidx] == refstring[qidx+repeatsize]) {
                qidx++;
            }
            unsigned int repeatnum = (qidx - refpos) / repeatsize + 1;
            if (is_indel_context_more_STR(repeatsize, repeatnum, repeatsize_at_max_repeatnum, max_repeatnum)) {
                repeatsize_at_max_repeatnum = repeatsize;
                max_repeatnum = repeatnum;
                repeat_endpos = qidx + repeatsize;
            }
        }
        assert(repeat_endpos > refpos);
        unsigned int tl = MIN(repeat_endpos, refstring.size()) - refpos;
        const unsigned int decphred = indel_phred(8.0*4.0, repeatsize_at_max_repeatnum, repeatsize_at_max_repeatnum, tl / repeatsize_at_max_repeatnum);
        for (unsigned int i = refpos; i != MIN(repeat_endpos, refstring.size()); i++) {
            if (tl > region_repeatvec[i].tracklen) {
                region_repeatvec[i].begpos = refpos;
                region_repeatvec[i].tracklen = tl;
                region_repeatvec[i].unitlen = repeatsize_at_max_repeatnum;
                region_repeatvec[i].indelphred = 40+3 - MIN(40+3-1, decphred);
            }
        }
        refpos += repeatsize_at_max_repeatnum * max_repeatnum;
    }
    return region_repeatvec;
}

void
repeatvec_LOG(const auto & region_repeatvec, unsigned int region_offset) {
    for (size_t i = 0; i < region_repeatvec.size(); i++) {
        LOG(logINFO) << "RegionalTandemRepeat\t" 
                << (region_offset + i) << "\t"
                << region_offset << "\t"
                << region_repeatvec[i].begpos << "\t" 
                << region_repeatvec[i].tracklen << "\t" 
                << region_repeatvec[i].unitlen;
    }
}

unsigned int
ref_to_phredvalue(unsigned int & n_units, const auto & refstring, const size_t refpos, 
        const unsigned int max_phred, double ampfact, const unsigned int cigar_oplen, const auto cigar_op, const bam1_t *b) {
    unsigned int max_repeatnum = 0;
    unsigned int repeatsize_at_max_repeatnum = 0;
    for (unsigned int repeatsize = 1; repeatsize <= 6; repeatsize++) {
        unsigned int qidx = refpos;
        while (qidx + repeatsize < refstring.size() && refstring[qidx] == refstring[qidx+repeatsize]) {
            qidx++;
        }
        unsigned int repeatnum = (qidx - refpos) / repeatsize + 1;
        if (is_indel_context_more_STR(repeatsize, repeatnum, repeatsize_at_max_repeatnum, max_repeatnum)) {
            max_repeatnum = repeatnum;
            repeatsize_at_max_repeatnum = repeatsize;
        }
    }
    if (cigar_oplen == repeatsize_at_max_repeatnum && cigar_op == BAM_CDEL) {
        ampfact *= 4.0;
    }
    // Because of a lower number of PCR cycles, it is higher than the one set in Fig. 3 at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/
    unsigned int decphred = indel_phred(ampfact, cigar_oplen, repeatsize_at_max_repeatnum, max_repeatnum);
    n_units = ((0 == cigar_oplen % repeatsize_at_max_repeatnum) ? (cigar_oplen / repeatsize_at_max_repeatnum) : 0);
    return max_phred - MIN(max_phred, decphred) + indel_len_rusize_phred(cigar_oplen, repeatsize_at_max_repeatnum); 
}

/*
unsigned int 
bam_to_phredvalue(unsigned int & n_units, const bam1_t *b, unsigned int qpos, unsigned int max_phred, double ampfact, const unsigned int cigar_oplen, const auto cigar_op) {
    unsigned int max_repeatnum = 0;
    unsigned int repeatsize_at_max_repeatnum = 0;
    for (unsigned int repeatsize = 1; repeatsize <= 6; repeatsize++) {
        unsigned int qidx = qpos;
        while (qidx + repeatsize < SIGN2UNSIGN(b->core.l_qseq) && bam_seqi(bam_get_seq(b), qidx) == bam_seqi(bam_get_seq(b), qidx+repeatsize)) {
            qidx++;
        }
        unsigned int repeatnum = (qidx - qpos) / repeatsize + 1;
        if (repeatnum > max_repeatnum) {
            max_repeatnum = repeatnum;
            repeatsize_at_max_repeatnum = repeatsize;
        }
    }
    if (cigar_oplen == repeatsize_at_max_repeatnum && cigar_op == BAM_CDEL) {
        ampfact *= 2.5;
    }
    // Because of a lower number of PCR cycles, it is higher than the one set in Fig. 3 at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/
    unsigned int decphred = indel_phred(ampfact, cigar_oplen, repeatsize_at_max_repeatnum, max_repeatnum);
    n_units = ((0 == cigar_oplen % repeatsize_at_max_repeatnum) ? (cigar_oplen / repeatsize_at_max_repeatnum) : 0);
    return max_phred - MIN(max_phred, decphred); 
}
*/

int
update_seg_format_prep_sets_by_aln(
        SegFormatPrepSets & seg_format_prep_sets,
        const bam1_t *aln,
        unsigned int specialflag) {
    
    const unsigned int rend = bam_endpos(aln);
    const auto cigar = bam_get_cigar(aln);
    const auto *bseq = bam_get_seq(aln);
    unsigned int nge_cnt = 0; // number of gap extensions.
    unsigned int ngo_cnt = 0; // number of gap opens.
    for (unsigned int i = 0; i < aln->core.n_cigar; i++) {
        int32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) { 
            nge_cnt += bam_cigar_oplen(c); 
            ngo_cnt++; 
        }
    }
    const uint8_t *bam_aux_data = bam_aux_get(aln, "NM");
    const unsigned int nm_cnt = ((bam_aux_data != NULL) ? bam_aux2i(bam_aux_data) : nge_cnt);
    assert (nm_cnt >= nge_cnt);
    const unsigned int xm_cnt = nm_cnt - nge_cnt;
    const unsigned int xm150 = xm_cnt * 150 / (rend - aln->core.pos);
    const unsigned int frag_pos_L = ((aln->core.isize != 0) ? MIN(aln->core.pos, aln->core.mpos) : aln->core.pos);
    const unsigned int frag_pos_R = ((aln->core.isize != 0) ? (frag_pos_L + aln->core.isize) : rend);
    const bool isrc = ((aln->core.flag & 0x10) == 0x10); 
    unsigned int qpos = 0;
    unsigned int rpos = aln->core.pos;
    
    for (unsigned int i = 0; i < aln->core.n_cigar; i++) {
        uint32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        unsigned int cigar_oplen = bam_cigar_oplen(c);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            for (unsigned int j = 0; j < cigar_oplen; j++) {
                // a_total_link_dp[rpos - roffset] += 1;
                seg_format_prep_sets.getRefByPos(rpos)[SEG_a_DP] += 1;
                seg_format_prep_sets.getRefByPos(rpos)[SEG_a_XM] += xm150;
                if (isrc) { 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_LI] += MIN(rpos - frag_pos_L + 1, 1200); 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_LIDP] += 1;
                } else { 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_RI] += MIN(frag_pos_R - rpos    , 1200); 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_RIDP] += 1;
                }
                qpos++;
                rpos++;
            }
        } else if (cigar_op == BAM_CINS) {
            unsigned int rtotlen = 0;
            for (int rpos2 = MAX((int)rpos - (int)cigar_oplen, aln->core.pos); rpos2 < rpos; rpos2++) {
                rtotlen++;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_INS_R] += rtotlen;
            }           
            unsigned int ltotlen = cigar_oplen;
            for (int rpos2 = (int)rpos; MIN(rpos2 < rpos + cigar_oplen, rend); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_INS_L] += ltotlen;
                ltotlen--;
            }
            qpos += cigar_oplen;
        } else if (cigar_op == BAM_CDEL) {
            unsigned int ltotlen = 0;
            unsigned int rtotlen = 0;
            for (int rpos2 = MAX((int)rpos - (int)cigar_oplen, aln->core.pos); rpos2 < rpos; rpos2++) {
                rtotlen++;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DEL_R] += rtotlen;
            }
            for (int rpos2 = (int)rpos; rpos2 < rpos + cigar_oplen; rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DEL_L] += ltotlen;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DEL_R] += rtotlen;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DP] += 1;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_XM] += xm150;
                if (isrc) { 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_LI] += MIN(rpos - frag_pos_L + 1, 1200); 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_LIDP] += 1;
                } else { 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_RI] += MIN(frag_pos_R - rpos    , 1200); 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_RIDP] += 1;
                }
                rtotlen--;
                ltotlen++;
            }
            for (int rpos2 = (int)rpos + (int)cigar_oplen; rpos2 < MIN(rpos + 2*cigar_oplen, rend); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DEL_L] += ltotlen;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DP] += 1;
                ltotlen--;
            }
            rpos += cigar_oplen;
        } else {
            process_cigar(qpos, rpos, cigar_op, cigar_oplen);
        }
    }
    assert(bam_endpos(aln) == rpos || !fprintf(stderr, "%u == %u failed for bam %s at tid %d position %d", 
            bam_endpos(aln), rpos, bam_get_qname(aln), aln->core.tid, aln->core.pos)); 
}

int // GenericSymbol2CountCoverage<TSymbol2Count>::
update_seg_format_depth_sets(
        auto & seg_format_depth_sets,
        auto & symbol_to_VQ_format_tag_sets,
        const auto & seg_format_thres_sets,
        const bam1_t *const aln,
        uint32_t specialflag) {
    
    const uint32_t n_cigar = aln->core.n_cigar;
    const uint32_t *cigar =  bam_get_cigar(aln);
    const uint8_t *bseq = bam_get_seq(aln);
    const auto rend = bam_endpos(aln);
    
    unsigned int nge_cnt = 0; // number of gap extensions.
    unsigned int ngo_cnt = 0; // number of gap opens.
    for (unsigned int i = 0; i < n_cigar; i++) {
        int32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) { 
            nge_cnt += bam_cigar_oplen(c); 
            ngo_cnt++; 
        }
    }
    const uint8_t *bam_aux_data = bam_aux_get(aln, "NM");
    const unsigned int nm_cnt = ((bam_aux_data != NULL) ? bam_aux2i(bam_aux_data) : nge_cnt);
    assert (nm_cnt >= nge_cnt);
    const unsigned int xm_cnt = nm_cnt - nge_cnt;
    unsigned int xm150 = xm_cnt * 150 / aln->core.l_qseq;
    // std::vector<uint16_t> qr_baq_vec = compute_qr_baq_vec(aln, region_repeatvec, region_offset, baq_per_aligned_base, 1);
    
    unsigned int qpos = 0;
    unsigned int rpos = aln->core.pos;
    unsigned int incvalue = 1;

    for (unsigned int i = 0; i < n_cigar; i++) {
        uint32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        unsigned int cigar_oplen = bam_cigar_oplen(c);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                if (i2 > 0) {
                    //assert (rpos >= region_repeatvec[rpos-roffset-1].begpos);
                    dealwith_seg_regbias(
                            (bam_phredi(aln, qpos - 1) + bam_phredi(aln, qpos)) / 2,
                            rpos,
                            seg_format_depth_sets.getRefByPos(rpos)[LINK_M],
                            symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[LINK_M],
                            seg_format_thres_sets.getByPos(rpos),
                            aln,
                            xm150,
                            0);
                }
                unsigned int base4bit = bam_seqi(bseq, qpos);
                unsigned int base3bit = seq_nt16_int[base4bit];
                AlignmentSymbol symbol = AlignmentSymbol(base3bit);
                dealwith_seg_regbias(
                        bam_phredi(aln, qpos),
                        rpos,
                        seg_format_depth_sets.getRefByPos(rpos)[symbol],
                        symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                        seg_format_thres_sets.getByPos(rpos),
                        aln,
                        xm150,
                        0);
                rpos += 1;
                qpos += 1;
            }
        } else if (cigar_op == BAM_CINS) {
            const bool is_ins_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(aln->core.l_qseq));
            if (!is_ins_at_read_end) {
                AlignmentSymbol symbol = insLenToSymbol(SIGN2UNSIGN(cigar_oplen));
                dealwith_seg_regbias(
                        (bam_phredi(aln, qpos) + bam_phredi(aln, qpos + cigar_oplen - 1)) / 2,
                        rpos,
                        seg_format_depth_sets.getRefByPos(rpos)[symbol],
                        symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                        seg_format_thres_sets.getByPos(rpos),
                        aln,
                        xm150,
                        0);
            }
            qpos += cigar_oplen;
        } else if (cigar_op == BAM_CDEL) {
            const bool is_del_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(aln->core.l_qseq));
            if (!is_del_at_read_end) {
                AlignmentSymbol symbol = delLenToSymbol(SIGN2UNSIGN(cigar_oplen));
                dealwith_seg_regbias(
                        (bam_phredi(aln, qpos) + bam_phredi(aln, qpos - 1)) / 2,
                        rpos,
                        seg_format_depth_sets.getRefByPos(rpos)[symbol],
                        symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                        seg_format_thres_sets.getByPos(rpos),
                        aln,
                        xm150,
                        0);
            }
            rpos += cigar_oplen;
        } else {
            process_cigar(qpos, rpos, cigar_op, cigar_oplen);
        }
    }
    return 0;
};

int
update_seg_format_thres_from_prep_sets(
        SegFormatThresSets & seg_format_thres_sets, 
        const SegFormatPrepSets & seg_format_prep_sets, 
        size_t region_offset,
        const std::vector<RegionalTandemRepeat> & region_repeatvec,
        const unsigned int specialflag) {
    assert(seg_format_thres_sets.getIncluBegPosition() == seg_format_prep_sets.getIncluBegPosition());
    assert(seg_format_thres_sets.getExcluEndPosition() == seg_format_prep_sets.getExcluEndPosition());
    for (size_t epos = seg_format_prep_sets.getIncluBegPosition(); epos != seg_format_prep_sets.getExcluEndPosition(); epos++) {
        const unsigned int seg_a_dp = MAX(seg_format_prep_sets.getByPos(epos)[SEG_a_DP], 1U);
        const unsigned int indel_len_per_DP = MAX4(
                seg_format_prep_sets.getByPos(epos)[SEG_a_INS_L], 
                seg_format_prep_sets.getByPos(epos)[SEG_a_INS_R], 
                seg_format_prep_sets.getByPos(epos)[SEG_a_DEL_L], 
                seg_format_prep_sets.getByPos(epos)[SEG_a_DEL_R] 
                ) / seg_a_dp;
        const auto potential_indel_len = MIN3(
                region_repeatvec[MAX(epos-region_offset, 1) - 1].tracklen, 
                region_repeatvec[MIN(epos-region_offset, region_repeatvec.size()-2) + 1].tracklen,
                region_repeatvec[epos-region_offset].tracklen);
        
        assert(seg_format_prep_sets.getByPos(epos)[SEG_a_LIDP] + seg_format_prep_sets.getByPos(epos)[SEG_a_RIDP] > 0);
        auto segLIDP = MAX(seg_format_prep_sets.getByPos(epos)[SEG_a_LIDP], 1);
        auto segRIDP = MAX(seg_format_prep_sets.getByPos(epos)[SEG_a_RIDP], 1);
        seg_format_thres_sets.getRefByPos(epos)[SEG_aEP1t] = MAX(indel_len_per_DP * 2, potential_indel_len) + 20;
        seg_format_thres_sets.getRefByPos(epos)[SEG_aEP2t] = MAX(indel_len_per_DP * 2, potential_indel_len) + 10; // the lower the strong the bias
        seg_format_thres_sets.getRefByPos(epos)[SEG_aXM1T] = seg_format_prep_sets.getByPos(epos)[SEG_a_XM] / seg_a_dp + 4;
        seg_format_thres_sets.getRefByPos(epos)[SEG_aXM2T] = seg_format_prep_sets.getByPos(epos)[SEG_a_XM] / seg_a_dp + 8; // the higher the stronger the bias
        seg_format_thres_sets.getRefByPos(epos)[SEG_aLI1T] = (seg_format_prep_sets.getByPos(epos)[SEG_a_LI] * 2 / segLIDP);
        seg_format_thres_sets.getRefByPos(epos)[SEG_aLI2T] = (seg_format_prep_sets.getByPos(epos)[SEG_a_LI] * 3 / segLIDP); // higher > stronger bias
        seg_format_thres_sets.getRefByPos(epos)[SEG_aRI1T] = (seg_format_prep_sets.getByPos(epos)[SEG_a_RI] * 2 / segRIDP);
        seg_format_thres_sets.getRefByPos(epos)[SEG_aRI2T] = (seg_format_prep_sets.getByPos(epos)[SEG_a_RI] * 3 / segRIDP); // higher > stronger bias
        seg_format_thres_sets.getRefByPos(epos)[SEG_aLI1t] = (seg_format_prep_sets.getByPos(epos)[SEG_a_LI] / 2 / segLIDP);
        seg_format_thres_sets.getRefByPos(epos)[SEG_aLI2t] = (seg_format_prep_sets.getByPos(epos)[SEG_a_LI] / 3 / segLIDP); // lower > stronger bias
        seg_format_thres_sets.getRefByPos(epos)[SEG_aRI1t] = (seg_format_prep_sets.getByPos(epos)[SEG_a_RI] / 2 / segRIDP);
        seg_format_thres_sets.getRefByPos(epos)[SEG_aRI2t] = (seg_format_prep_sets.getByPos(epos)[SEG_a_RI] / 3 / segRIDP); // lower > stronger bias
    }
}

inline const
int
dealwith_seg_regbias(
        const auto bq,
        const unsigned int rpos,
        auto & symbol_to_seg_format_depth_set,
        auto & symbol_to_VQ_format_tag_set,
        const auto & seg_format_thres_set,
        const bam1_t *aln, 
        const auto xm150,
        const unsigned int specialflag) {
    
    const auto rend = bam_endpos(aln); 
    const size_t seg_l_nbases = (rpos - aln->core.pos + 1);
    const size_t seg_r_nbases = (bam_endpos(aln) - aln->core.pos);
    const size_t frag_pos_L = ((aln->core.isize != 0) ? MIN(aln->core.pos, aln->core.mpos) : aln->core.pos);
    const size_t frag_pos_R = ((aln->core.isize != 0) ? (frag_pos_L + aln->core.isize) : rend);
    const size_t frag_l_nbases = (rpos - frag_pos_L + 1); 
    const size_t frag_r_nbases = (frag_pos_R - rpos);
    const bool isrc = ((aln->core.flag & 0x10) == 0x10);
    const bool isr2 = ((aln->core.flag & 0x80) == 0x80 && (aln->core.flag & 0x1) == 0x1);
    const bool strand = (isrc ^ isr2);
    
    const auto aSBQ = (isrc ? VQ_aSBQr : VQ_aSBQf);
    symbol_to_VQ_format_tag_set[aSBQ] += bq;
    
    auto & symbol_to_seg_aDP_depth_set = (strand
        ? (isrc ? symbol_to_seg_format_depth_set[SEG_aDPrr] : symbol_to_seg_format_depth_set[SEG_aDPrf]) 
        : (isrc ? symbol_to_seg_format_depth_set[SEG_aDPfr] : symbol_to_seg_format_depth_set[SEG_aDPff]));
    symbol_to_seg_aDP_depth_set += 1;
    
    // const auto aBQidx = (isrc ? SEG_aBQa : SEG_aBQb);
    // symbol_to_seg_format_depth_set[SEG_aBQa] += 1;
    
    if ((seg_l_nbases < seg_format_thres_set[SEG_aEP1t]) || (seg_r_nbases < seg_format_thres_set[SEG_aEP1t])) {
        symbol_to_seg_format_depth_set[SEG_aEP1] += 1; // unmerged position quality
    }
    if ((seg_l_nbases < seg_format_thres_set[SEG_aEP2t]) || (seg_r_nbases < seg_format_thres_set[SEG_aEP2t])) {
        symbol_to_seg_format_depth_set[SEG_aEP2] += 1;
    }
    
    if (xm150 >= seg_format_thres_set[SEG_aXM1T]) {
        symbol_to_seg_format_depth_set[SEG_aXM1] += 1; // unmerged base quality
    } 
    if (xm150 >= seg_format_thres_set[SEG_aXM2T]) {
        symbol_to_seg_format_depth_set[SEG_aXM2] += 1;
    }
    
    if (isrc) {
        auto dist2iend = aln->core.pos - frag_pos_L;
        if (dist2iend <= seg_format_thres_set[SEG_aLI1T] && dist2iend >= seg_format_thres_set[SEG_aLI1t]) {
            symbol_to_seg_format_depth_set[SEG_aLI1] += 1; // unmerged base quality
        } 
        if (dist2iend <= seg_format_thres_set[SEG_aLI2T] && dist2iend >= seg_format_thres_set[SEG_aLI2t]) {
            symbol_to_seg_format_depth_set[SEG_aLI2] += 1;
        }
    } else {
        auto dist2iend = frag_pos_R - aln->core.pos;
        if (dist2iend <= seg_format_thres_set[SEG_aRI1T] && dist2iend >= seg_format_thres_set[SEG_aRI1t]) {
            symbol_to_seg_format_depth_set[SEG_aRI1] += 1; // unmerged base quality
        } 
        if (dist2iend <= seg_format_thres_set[SEG_aRI2T] && dist2iend >= seg_format_thres_set[SEG_aRI2t]) {
            symbol_to_seg_format_depth_set[SEG_aRI2] += 1;
        }
    }
    //auto baq = MIN(THE_BAQ_MAX, qr_baq_vec.at(rpos - b->core.pos));
    //symbol_to_seg_format_depth_sets[SEG_SQR_BAQ_A].template inc<SYMBOL_COUNT_SUM>(rpos, symbol, 1, baq);
    
    return 0;
}

template <class TSymbol2Count>
class GenericSymbol2CountCoverage : public CoveredRegion<TSymbol2Count> {
public:
    GenericSymbol2CountCoverage() : CoveredRegion<TSymbol2Count>(0, 0, 1) { /* assert(false); */ };
    GenericSymbol2CountCoverage(auto tid, auto beg, auto end) : CoveredRegion<TSymbol2Count>(tid, beg, end) {}
    
    void
    assertUpdateIsLegal(const GenericSymbol2CountCoverage<TSymbol2Count> & other) const {
        assert(this->tid == other.tid);
        assert(this->getIncluBegPosition() <= other.getIncluBegPosition() || !fprintf(stderr, "%lu <= %lu failed!", this->getIncluBegPosition(), other.getIncluBegPosition()));
        assert(this->getExcluEndPosition() >= other.getExcluEndPosition() || !fprintf(stderr, "%lu >= %lu failed!", this->getExcluEndPosition(), other.getExcluEndPosition())); 
    }
    // mainly for merging reads in one family
    template <bool TIsIncVariable = true>
    void
    updateByRepresentative(const GenericSymbol2CountCoverage<TSymbol2Count> & other, unsigned int incvalue = 1,
            const bool update_pos2indel2count = true, const bool update_idx2symbol2data = true) {
        this->assertUpdateIsLegal(other);
        
        if (update_idx2symbol2data) {
            for (size_t epos = other.getIncluBegPosition(); epos < other.getExcluEndPosition(); epos++) {
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
    template</*ValueType T_ConsensusType,*/ bool TIndelIsMajor = false> 
    void
    updateByConsensus(const GenericSymbol2CountCoverage<TSymbol2Count> &other, unsigned int incvalue = 1,
            const bool update_pos2indel2count = true, const bool update_idx2symbol2data = true) {
        this->assertUpdateIsLegal(other);
        
        if (update_idx2symbol2data) {
            for (size_t epos = other.getIncluBegPosition(); epos < other.getExcluEndPosition(); epos++) {
                const std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols = this->getRefByPos(epos).template updateByConsensus</*T_ConsensusType,*/ TIndelIsMajor>(other.getByPos(epos), incvalue);
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
    /*
    int
    updateByFiltering(
            const GenericSymbol2CountCoverage<TSymbol2Count> &other, 
            const unsigned int thres, // const GenericSymbol2CountCoverage<TSymbol2Count> &thres, 
            unsigned int incvalue = 1,
            const bool update_pos2indel2count = true) {
        this->assertUpdateIsLegal(other);
        int num_updated_pos = 0;
        // other.assertUpdateIsLegal(thres); // may not hold because threshold is delimited by bed whereas this is not delimited by bed.
        auto incluBegPos = other.getIncluBegPosition(); // MAX(other.getIncluBegPosition(), thres.getIncluBegPosition()); 
        auto excluEndPos = other.getExcluEndPosition(); // MIN(other.getExcluEndPosition(), thres.getExcluEndPosition());
        std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols; 
        for (auto epos = incluBegPos; epos < excluEndPos; epos++) {
            int updateresult = this->getRefByPos(epos).updateByFiltering(
                    consymbols, 
                    other.getByPos(epos), 
                    thres, // thres.getByPos(epos), 
                    incvalue, 
                    epos);
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
    */
    int
    updateByMajorMinusMinor(
            const GenericSymbol2CountCoverage<TSymbol2Count> &other,
            const bool update_pos2indel2count = true) {
        this->assertUpdateIsLegal(other);
        int num_updated_pos = 0;
        // other.assertUpdateIsLegal(thres); // may not hold because threshold is delimited by bed whereas this is not delimited by bed.
        auto incluBegPos = other.getIncluBegPosition();
        auto excluEndPos = other.getExcluEndPosition();
        std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols; 
        std::array<unsigned int   , NUM_SYMBOL_TYPES> conmmmcnts;
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
    inc(const unsigned int epos, const AlignmentSymbol symbol, const unsigned int incvalue = 1, const bam1_t *bam = NULL) {
        auto &r = this->getRefByPos(epos, bam);
        r.template incSymbolCount<(TUpdateType2)>(symbol, incvalue);
    };

    void // GenericSymbol2CountCoverage<TSymbol2Count>::
    incIns(const unsigned int epos, const std::string & iseq, const AlignmentSymbol symbol, const uint32_t incvalue = 1) {
        assert (incvalue > 0); 
        assert (iseq.size() > 0);
        size_t ipos = epos; //  _extern2intern4pos(epos);
        posToIndelToCount_inc(this->getRefPosToIseqToData(symbol), ipos, iseq, incvalue);
    };

    void // GenericSymbol2CountCoverage<TSymbol2Count>::
    incDel(const unsigned int epos, const uint32_t dlen, const AlignmentSymbol symbol, const uint32_t incvalue = 1) {
        assert (incvalue > 0);
        assert (dlen > 0);
        size_t ipos = epos; // _extern2intern4pos(epos);
        posToIndelToCount_inc(this->getRefPosToDlenToData(symbol), ipos, dlen, incvalue);
    };
    
    std::vector<uint16_t> 
    compute_qr_baq_vec(const bam1_t *b, const auto & region_repeatvec, const unsigned int region_offset, unsigned int baq_per_aligned_base, unsigned int baq_per_additional_base) {
        const unsigned int rbeg = b->core.pos;
        const unsigned int rend = bam_endpos(b);
        const uint8_t *bseq = bam_get_seq(b);
        const uint32_t n_cigar = b->core.n_cigar;
        const uint32_t *cigar = bam_get_cigar(b);
        
        std::vector<uint16_t> refBAQvec(rend-rbeg);
        {
            unsigned int qpos = 0;
            unsigned int rpos = b->core.pos;
            unsigned int baq = 0;
            unsigned int old_baq = 0;
            unsigned int prev_repeat_begpos = region_repeatvec.at(rpos - region_offset).begpos; // UINT_MAX;
            // unsigned int prev_rpos = rbeg;
            // unsigned int baq_per_aligned_repeatregion = baq_per_aligned_base; // baq_per_additional_base;
            for (int i = 0; i != n_cigar; i++) {
                const uint32_t c = cigar[i];
                const unsigned int cigar_op = bam_cigar_op(c);
                const unsigned int cigar_oplen = bam_cigar_oplen(c);
                if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                    for (unsigned int i2 = 0; i2 != cigar_oplen; i2++) {
                        unsigned int repeat_begpos = region_repeatvec.at(rpos - region_offset).begpos;
                        assert(repeat_begpos >= prev_repeat_begpos);
                        if (prev_repeat_begpos != repeat_begpos) {
                            baq += baq_per_aligned_base; 
                            old_baq = baq;
                        } else {
                            baq = MIN(baq + baq_per_additional_base, old_baq + baq_per_aligned_base);
                        }
                        refBAQvec[rpos - rbeg] = baq;
                        prev_repeat_begpos = repeat_begpos;
                        qpos++;
                        rpos++;
                    }
                } else if (cigar_op == BAM_CINS) {
                    qpos += cigar_oplen;
                } else if (cigar_op == BAM_CDEL) {
                    for (unsigned int i2 = 0; i2 != cigar_oplen; i2++) {
                        refBAQvec[rpos - rbeg] = baq;
                        rpos++;
                    }
                }  else if (cigar_op == BAM_CSOFT_CLIP) {
                    qpos += cigar_oplen;
                }
            }
            assert(b->core.l_qseq == qpos || !fprintf(stderr, "%u == %u failed for bam %s", b->core.l_qseq, qpos, bam_get_qname(b)));
            assert(rend == rpos || !fprintf(stderr, "%u == %u failed for bam %s", rend, rpos, bam_get_qname(b))); 
        }
#if 0
        if (b->core.pos % 1000 == 0) {
            LOG(logINFO) << "For readname " << bam_get_qname(b) << " at tid " << b->core.tid << " position " <<  b->core.pos;
            std::string tmpstring = "";
            for (auto a : refBAQvec) {
                tmpstring += std::to_string(a) + ",";
            }
            LOG(logINFO) << tmpstring;
            tmpstring = "";
            for (unsigned int i = b->core.pos; i < bam_endpos(b); i++) {
                tmpstring += std::to_string(region_repeatvec.at(i - region_offset).begpos) + ",";
            }
            LOG(logINFO) << tmpstring;
        }
#endif
        {
            unsigned int qpos = b->core.l_qseq - 1;
            unsigned int rpos = rend - 1;
            unsigned int baq = 0;
            unsigned int old_baq = 0;
            unsigned int prev_repeat_begpos = region_repeatvec.at(rpos - region_offset).begpos; // UINT_MAX;
            for (int i = n_cigar - 1; i != -1; i--) {
                const uint32_t c = cigar[i];
                const unsigned int cigar_op = bam_cigar_op(c);
                const unsigned int cigar_oplen = bam_cigar_oplen(c);
                if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                    for (unsigned int i2 = 0; i2 != cigar_oplen; i2++) {
                        unsigned int repeat_begpos = region_repeatvec.at(rpos - region_offset).begpos;
                        assert(repeat_begpos <= prev_repeat_begpos);
                        if (repeat_begpos == 0) {
                            // LOG(logWARNING) << "repeat_begpos = 0 for readname " << bam_get_qname(b) << " at tid " << b->core.tid << " position " <<  b->core.pos;
                        }
                        if (prev_repeat_begpos != repeat_begpos) { 
                            baq += baq_per_aligned_base;
                            old_baq = baq;
                        } else {
                            baq = MIN(baq + baq_per_additional_base, old_baq + baq_per_aligned_base);
                        }
                        UPDATE_MIN(refBAQvec[rpos - rbeg], baq);
                        prev_repeat_begpos = repeat_begpos;
                        qpos--;
                        rpos--;
                    }
                } else if (cigar_op == BAM_CINS) {
                    qpos -= cigar_oplen;
                } else if (cigar_op == BAM_CDEL) {
                    for (unsigned int i2 = 0; i2 != cigar_oplen; i2++) {
                        UPDATE_MIN(refBAQvec[rpos - rbeg], baq);
                        rpos--;
                    }
                } else if (cigar_op == BAM_CSOFT_CLIP) {
                    qpos -= cigar_oplen;
                }
            }
            assert(-1 == qpos || !fprintf(stderr, "-1 == %u failed for bam %s", qpos, bam_get_qname(b)));
            assert(rbeg - 1 == rpos || !fprintf(stderr, "%u - 1 == %u failed for bam %s", rbeg, rpos, bam_get_qname(b))); 
        }
        return refBAQvec;
    }
    
    template<bool TIsProton>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByAln(
            const bam1_t *const aln, 
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhredArg, 
            unsigned int frag_indel_basemax,
            const unsigned int region_offset,
            const auto & region_symbolvec, 
            const std::vector<RegionalTandemRepeat> & region_repeatvec,
            uint32_t sepcialflag) {
        
        // static_assert(BASE_QUALITY_MAX == TUpdateType || SYMBOL_COUNT_SUM == TUpdateType);
        assert(this->tid == SIGN2UNSIGN(aln->core.tid));
        assert(this->getIncluBegPosition() <= SIGN2UNSIGN(aln->core.pos)   || !fprintf(stderr, "%lu <= %d failed", this->getIncluBegPosition(), aln->core.pos));
        assert(this->getExcluEndPosition() >= SIGN2UNSIGN(bam_endpos(aln)) || !fprintf(stderr, "%lu >= %d failed", this->getExcluEndPosition(), bam_endpos(aln)));
        
        const uint32_t n_cigar = aln->core.n_cigar;
        const uint32_t *cigar = bam_get_cigar(aln);
        const uint8_t *bseq = bam_get_seq(aln);
        const auto symbolType2addPhred = symbolType2addPhredArg; // std::array({0, 0});
        // const auto roffset = bq_baq_sum.at(0).getIncluBegPosition();
        
        /*
        unsigned int nge_cnt = 0; // number of gap extensions.
        unsigned int ngo_cnt = 0; // number of gap opens.
        for (unsigned int i = 0; i < n_cigar; i++) {
            int32_t c = cigar[i];
            unsigned int cigar_op = bam_cigar_op(c);
            if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) { 
                nge_cnt += bam_cigar_oplen(c); 
                ngo_cnt++; 
            }
        }
        const uint8_t *bam_aux_data = bam_aux_get(aln, "NM");
        const unsigned int nm_cnt = ((bam_aux_data != NULL) ? bam_aux2i(bam_aux_data) : nge_cnt);
        assert (nm_cnt >= nge_cnt);
        const unsigned int xm_cnt = nm_cnt - nge_cnt;
        const unsigned int max_noindel_phred = MIN(40, prob2phred((xm_cnt + 1) / (double)(aln->core.l_qseq * 4 * 2 + 2)));
        */
        unsigned int qpos = 0;
        unsigned int rpos = aln->core.pos;
        unsigned int incvalue = 1;
        for (unsigned int i = 0; i < n_cigar; i++) {
            uint32_t c = cigar[i];
            unsigned int cigar_op = bam_cigar_op(c);
            unsigned int cigar_oplen = bam_cigar_oplen(c);
            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                    assert((rpos >= SIGN2UNSIGN(aln->core.pos) && rpos < SIGN2UNSIGN(bam_endpos(aln)))
                            || !fprintf(stderr, "Bam line with QNAME %s has rpos that is not within the range (%d - %d)", bam_get_qname(aln), aln->core.pos, bam_endpos(aln)));
                    if (i2 > 0) {
                        const auto noindel_phred = MIN3(
                                region_repeatvec[rpos-region_offset - 1].indelphred,
                                region_repeatvec[rpos-region_offset].indelphred,
                                40);
                        incvalue = ((TIsProton) ? (MIN3(noindel_phred, bam_phredi(aln, qpos-1), bam_phredi(aln, qpos))) : noindel_phred); 
                        this->template inc<BASE_QUALITY_MAX>(rpos, LINK_M, incvalue, aln);
                    }
                    unsigned int base4bit = bam_seqi(bseq, qpos);
                    unsigned int base3bit = seq_nt16_int[base4bit];
                    incvalue = bam_phredi(aln, qpos) + symbolType2addPhred[BASE_SYMBOL];
                    this->template inc<BASE_QUALITY_MAX>(rpos, AlignmentSymbol(base3bit), incvalue, aln);
                    rpos += 1;
                    qpos += 1;
                }
            } else if (cigar_op == BAM_CINS) {
                const bool is_ins_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(aln->core.l_qseq));
                unsigned int inslen = SIGN2UNSIGN(cigar_oplen);
                if (is_ins_at_read_end) {
                    LOG(logWARNING) << "Query " << bam_get_qname(aln) << " has insertion of legnth " << cigar_oplen << " at " << qpos
                            << " which is not exclusively between 0 and " << aln->core.l_qseq << " aligned to tid " << aln->core.tid << " and position " << rpos;
                    incvalue = (0 != qpos ? bam_phredi(aln, qpos-1) : 
                            ((qpos + cigar_oplen < SIGN2UNSIGN(aln->core.l_qseq)) ? 
                            bam_phredi(aln, qpos + SIGN2UNSIGN(cigar_oplen)) : 1)) 
                            + (symbolType2addPhred[LINK_SYMBOL]); // + addidq; // 
                } else {
                    unsigned int phredvalue = ref_to_phredvalue(inslen, region_symbolvec, rpos - region_offset,
                            frag_indel_basemax, 8.0, cigar_oplen, cigar_op, aln);
                    incvalue = MIN(MIN(bam_phredi(aln, qpos-1), bam_phredi(aln, qpos + cigar_oplen)) + (TIsProton ? proton_cigarlen2phred(cigar_oplen) : 0), phredvalue)
                            + (symbolType2addPhred[LINK_SYMBOL]); // + addidq; 
                }
                if (!is_ins_at_read_end) {
                    this->template inc<BASE_QUALITY_MAX>(rpos, insLenToSymbol(inslen), MAX(SIGN2UNSIGN(1), incvalue), aln);
                    std::string iseq;
                    iseq.reserve(cigar_oplen);
                    unsigned int incvalue2 = incvalue;
                    for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                        unsigned int base4bit = bam_seqi(bseq, qpos+i2);
                        const char base8bit = seq_nt16_str[base4bit];
                        iseq.push_back(base8bit);
                        incvalue2 = MIN(incvalue2, SIGN2UNSIGN(bam_seqi(bseq, qpos+i2)))
                                + (symbolType2addPhred[LINK_SYMBOL]); // + symbolType2addPhred[LINK_SYMBOL];
                    }
                    this->incIns(rpos, iseq, insLenToSymbol(inslen), MAX(SIGN2UNSIGN(1), incvalue2));
                }
                qpos += cigar_oplen;
            } else if (cigar_op == BAM_CDEL) {
                const bool is_del_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(aln->core.l_qseq));
                unsigned int dellen = SIGN2UNSIGN(cigar_oplen);
                if (is_del_at_read_end) {
                    LOG(logWARNING) << "Query " << bam_get_qname(aln) << " has deletion of legnth " << cigar_oplen << " at " << qpos
                            << " which is not exclusively between 0 and " << aln->core.l_qseq << " aligned to tid " << aln->core.tid << " and position " << rpos; 
                    incvalue = (0 != qpos ? bam_phredi(aln, qpos-1) : 
                            ((qpos + cigar_oplen < SIGN2UNSIGN(aln->core.l_qseq)) ? 
                            bam_phredi(aln, qpos + SIGN2UNSIGN(cigar_oplen)) : 1))
                            + (symbolType2addPhred[LINK_SYMBOL]); // + addidq;
                } else {
                    unsigned int phredvalue = ref_to_phredvalue(dellen, region_symbolvec, rpos - region_offset, 
                            frag_indel_basemax, 8.0, cigar_oplen, cigar_op, aln);
                    incvalue = MIN(MIN(bam_phredi(aln, qpos), bam_phredi(aln, qpos-1)) + (TIsProton ? proton_cigarlen2phred(cigar_oplen) : 0), phredvalue)
                            + (symbolType2addPhred[LINK_SYMBOL]);
                }
                if (!is_del_at_read_end) {
                    this->template inc<BASE_QUALITY_MAX>(rpos, delLenToSymbol(dellen), MAX(SIGN2UNSIGN(1), incvalue), aln);
                    this->incDel(rpos, cigar_oplen, delLenToSymbol(dellen), MAX(SIGN2UNSIGN(1), incvalue));
                }
#if 1 
// The definition of non-ref for indel is not clearly defined, this piece of code can result in germline risk that is not in the ground truth.
// However, if we consider any position covered by the boundary of a germline indel to be non-ref, then this code should be enabled.
                unsigned int endpos = SIGN2UNSIGN(bam_endpos(aln));
                for (unsigned int p = rpos+1; p < MIN(rpos + cigar_oplen + 1, endpos); p++) {
                    this->template inc<BASE_QUALITY_MAX>(p, LINK_NN, MAX(SIGN2UNSIGN(1), incvalue), aln);
                }
                for (unsigned int p = rpos; p < MIN(rpos + cigar_oplen, endpos); p++) {
                    this->template inc<BASE_QUALITY_MAX>(p, BASE_NN, MAX(SIGN2UNSIGN(1), incvalue), aln);
                }
#endif
                rpos += cigar_oplen;
            } else {
                process_cigar(qpos, rpos, cigar_op, cigar_oplen);
            } 
        }
        return 0;
    }
    
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByRead1Aln(
            const std::vector<bam1_t *> & aln_vec,
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred, 
            const unsigned int frag_indel_basemax, 
            const bool is_proton, 
            const unsigned int region_offset,
            const auto & region_symbolvec,
            const auto & region_repeatvec,
            unsigned int specialflag) {
        for (const bam1_t *aln : aln_vec) {
            if (is_proton) {
                this->template updateByAln<true>(
                        aln, 
                        symbolType2addPhred, 
                        frag_indel_basemax, 
                        region_offset, 
                        region_symbolvec, 
                        region_repeatvec, 
                        0);
            } else {
                this->template updateByAln<false>(
                        aln, 
                        symbolType2addPhred, 
                        frag_indel_basemax, 
                        region_offset, 
                        region_symbolvec, 
                        region_repeatvec, 
                        0);
            } 
        }
        return 0;
    }
};

typedef GenericSymbol2CountCoverage<Symbol2Count> Symbol2CountCoverage; 
typedef GenericSymbol2CountCoverage<Symbol2CountUint64> Symbol2CountCoverageUint64; 
typedef GenericSymbol2CountCoverage<std::array<std::string, NUM_ALIGNMENT_SYMBOLS>> Symbol2CountCoverageString; 

struct Symbol2CountCoverageSet {
    const uint32_t tid;
    const uint32_t incluBegPosition;
    const uint32_t excluEndPosition;
    std::string refstring;
    
    SegFormatPrepSets seg_format_prep_sets;
    SegFormatThresSets seg_format_thres_sets;
    Symbol2SegFormatDepthSets symbol_to_seg_format_depth_sets;
    std::array<Symbol2FragFormatDepthSets, 2> symbol_to_frag_format_depth_sets; 
    std::array<Symbol2FamFormatDepthSets, 2> symbol_to_fam_format_depth_sets_2strand;
    Symbol2DuplexFormatDepthSets symbol_to_duplex_format_depth_sets;
    Symbol2VQFormatTagSets symbol_to_VQ_format_tag_sets;
    
    std::array<Symbol2Bucket2CountCoverage, 2> dedup_ampDistr; // dedup, important for error correction by barcoding, need major and minor, second pass
    Symbol2CountCoverageString additional_note;
    
    Symbol2CountCoverageSet(unsigned int t, unsigned int beg, unsigned int end):
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
    
    const size_t
    getUnifiedIncluBegPosition() const {
        return incluBegPosition;
    };
    const size_t
    getUnifiedExcluEndPosition() const {
        return excluEndPosition;
    };
    
    int 
    updateByAlns3UsingBQ(
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map,
            const std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> & alns3, 
            const std::basic_string<AlignmentSymbol> & region_symbolvec,
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred,
            bool should_add_note, 
            unsigned int frag_indel_basemax,
            const PhredMutationTable & phred_max_table, 
            const bool is_proton, 
            //const unsigned int baq_per_aligned_base,
            //const unsigned int sidereg_nbases,
            const auto & region_repeatvec, 
            unsigned int specialflag) {

        for (const auto & alns2pair2dflag : alns3) {
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair2dflag.first[strand];
                for (const auto & alns1 : alns2) {
                    for (const bam1_t *aln : alns1) {
                        update_seg_format_prep_sets_by_aln(
                                this->seg_format_prep_sets,
                                aln,
                                0);
                    }
                }
            }
        }
        update_seg_format_thres_from_prep_sets(
                this->seg_format_thres_sets,
                this->seg_format_prep_sets,
                this->getUnifiedIncluBegPosition(),
                region_repeatvec,
                0);
        for (const auto & alns2pair2dflag : alns3) {
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair2dflag.first[strand];
                for (const auto & alns1 : alns2) {
                    for (const bam1_t *aln : alns1) {
                        update_seg_format_depth_sets(
                                this->symbol_to_seg_format_depth_sets,
                                this->symbol_to_VQ_format_tag_sets,
                                this->seg_format_thres_sets,
                                aln, 
                                0);
                    }
                }
            }
        }
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                for (const auto & alns1 : alns2) {
                    uint32_t tid2, beg2, end2;
                    fillTidBegEndFromAlns1(tid2, beg2, end2, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid, beg2, end2);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln(
                            alns1, 
                            symbolType2addPhred, 
                            frag_indel_basemax, 
                            is_proton, 
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec, 
                            //this->a_edge_dist_l,
                            //this->a_edge_dist_r,
                            region_repeatvec,
                            0);
                    unsigned int normMQ = 0;
                    for (const bam1_t *aln : alns1) {
                        normMQ = MAX(normMQ, (unsigned int)aln->core.qual);
                    }
                    std::basic_string<std::pair<unsigned int, AlignmentSymbol>> pos_symbol_string;
                    for (auto epos = read_ampBQerr_fragWithR1R2.getIncluBegPosition(); epos < read_ampBQerr_fragWithR1R2.getExcluEndPosition(); epos++) {
                        for (SymbolType symbolType : SYMBOL_TYPES_IN_VCF_ORDER) {
                            AlignmentSymbol con_symbol;
                            unsigned int con_count, tot_count;
                            if (LINK_SYMBOL == symbolType) {
                                read_ampBQerr_fragWithR1R2.getByPos(epos).template fillConsensusCounts<true >(con_symbol, con_count, tot_count, symbolType);
                            } else {
                                read_ampBQerr_fragWithR1R2.getByPos(epos).template fillConsensusCounts<false>(con_symbol, con_count, tot_count, symbolType); 
                            }
                            assert (con_count * 2 >= tot_count);
                            if (0 == tot_count) { continue; }
                            
                            int max_qual = 8 + seg_format_get_avgBQ(
                                    this->symbol_to_seg_format_depth_sets.getByPos(epos)[con_symbol], 
                                    this->symbol_to_VQ_format_tag_sets.getByPos(epos)[con_symbol]);
                            unsigned int phredlike = MIN(con_count * 2 - tot_count, max_qual);
                            
                            // shared between BQ and FQ
                            // con_symbols_vec[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()][symbolType] = con_symbol;
                            
                            int pbucket = max_qual - phredlike;
                            if (pbucket < NUM_BUCKETS) {
                                this->dedup_ampDistr[0].getRefByPos(epos).incSymbolBucketCount(con_symbol, pbucket, 1);
                            }
                            // unsigned int frag_bAD_idx = ((0 == strand) ? FRAG_bDP1 : FRAG_bDP2);
                            
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
                            if (areSymbolsMutated(refsymbol, con_symbol) && (LINK_SYMBOL == symbolType || phredlike >= 20)) {
                                pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
                            }
                        }
                    }
                    if (pos_symbol_string.size() > 1) {
                        mutform2count4map.insert(std::make_pair(pos_symbol_string, std::array<unsigned int, 2>({0, 0})));
                        mutform2count4map[pos_symbol_string][strand]++;
                    }
                }
            }
        }
        assert(this->symbol_to_seg_format_depth_sets.getExcluEndPosition() - this->symbol_to_seg_format_depth_sets.getIncluBegPosition() == dedup_ampDistr.size());
        for (unsigned int epos = this->getUnifiedIncluBegPosition(); epos < this->getUnifiedExcluEndPosition(); epos++) {
            for (SymbolType symbolType : SYMBOL_TYPE_ARR) {
                for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symbolType]) {
                    const auto totDP = 
                            formatSumBySymbolType(this->symbol_to_frag_format_depth_sets[0].getByPos(epos), 
                            symbolType, FRAG_bDP) + 
                            formatSumBySymbolType(this->symbol_to_frag_format_depth_sets[1].getByPos(epos), 
                            symbolType, FRAG_bDP); 
                    const int max_qual = 8 + seg_format_get_avgBQ(
                            this->symbol_to_seg_format_depth_sets.getByPos(epos)[symbol],
                            this->symbol_to_VQ_format_tag_sets.getByPos(epos)[symbol]);
                    int maxvqual = 0;
                    unsigned int argmaxAD = 0;
                    unsigned int argmaxBQ = 0;
                    infer_max_qual_assuming_independence(
                            maxvqual,
                            argmaxAD,
                            argmaxBQ,
                            max_qual,
                            1, 
                            dedup_ampDistr[0].getByPos(epos).getSymbolCounts(symbol),
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
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map,
            const auto & alns3, 
            const std::basic_string<AlignmentSymbol> & region_symbolvec,
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred, 
            bool should_add_note,
            const unsigned int frag_indel_basemax,
            const PhredMutationTable & phred_max_table,
            const bool is_proton, 
            const auto & region_repeatvec,
            const unsigned int phred_indel_error_before_barcode_labeling,
            unsigned int specialflag = 0) {
        bool is_loginfo_enabled = true; 
        
        unsigned int niters = 0;
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            niters++;
            // bool log_alns2 = (is_loginfo_enabled && ispowerof2(niters));
            assert (alns2pair[0].size() != 0 || alns2pair[1].size() != 0);
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                if (alns2.size() == 0) { continue; }
                uint32_t tid2, beg2, end2;
                fillTidBegEndFromAlns2(tid2, beg2, end2, alns2);
                Symbol2CountCoverage read_family_con_ampl(tid2, beg2, end2); 
                for (const auto & alns1 : alns2) {
                    uint32_t tid1, beg1, end1;
                    fillTidBegEndFromAlns1(tid1, beg1, end1, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid1, beg1, end1);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln(
                            alns1, 
                            symbolType2addPhred, 
                            frag_indel_basemax, 
                            is_proton, 
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec, 
                            region_repeatvec,
                            0);
                    read_family_con_ampl.template updateByConsensus</*SYMBOL_COUNT_SUM,*/ true>(read_ampBQerr_fragWithR1R2);
                }
                for (size_t epos = read_family_con_ampl.getIncluBegPosition(); epos < read_family_con_ampl.getExcluEndPosition(); epos++) {
                    const auto & con_ampl_symbol2count = read_family_con_ampl.getByPos(epos);
                    for (SymbolType symbolType : SYMBOL_TYPE_ARR) {
                        AlignmentSymbol con_symbol; // = END_ALIGNMENT_SYMBOLS;
                        unsigned int con_count, tot_count;
                        con_ampl_symbol2count.fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        if (0 == tot_count) { continue ; }
                        this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP1] += 1;
                        if (1 == tot_count) {
                            this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_c1DP] += 1;
                        }
                        if (1 < tot_count && (con_count * 5 >= tot_count * 4)) {
                            this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP2] += 1;
                        }
                        if (7 < tot_count && (con_count * 5 >= tot_count * 4)) {
                            this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP3] += 1;
                        }
                        // this->fam_total_dep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1); // for genomic region, so add DP (every type of  DP is unfiltered)
                        /*if (1 == tot_count) {
                            this->fam_size1_dep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        } else if ((con_count * 5 < tot_count * 4)) { 
                            this->fam_nocon_dep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        }
                        */
                        //AlignmentSymbol con_symbol; // = END_ALIGNMENT_SYMBOLS;
                        //unsigned int con_count, tot_count;
                        //con_ampl_symbol2count.fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        // This is one round of bootstrapping for EM. We can use a full EM algorithm to improve this estimator, but it is probably overkill
                        // 0 means no coverage, 1 means no error correction, 2 means low quality family if the symbols disagree with each other, 
                        // 3 means up to 1 erroneous basecall is tolerated, 4 means up to 2 erroneous basecalls are tolerated
                        if (tot_count < 4) { continue; } 
                        if (con_count * 3 < tot_count * 2) { continue; }
                        for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symbolType]) {
                            if (con_symbol != symbol) {
                            // this->minor_amplicon[strand].getRefByPos(epos).incSymbolCount(symbol, con_ampl_symbol2count.getSymbolCount(symbol));
                            // this->major_amplicon[strand].getRefByPos(epos).incSymbolCount(symbol, tot_count);
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
            bool log_alns2 = (is_loginfo_enabled && ispowerof2(niters));
            uint32_t tid2, beg2, end2;
            tid2 = 0; beg2 = UINT32_MAX; end2 = 0; bool initialized = false;
            assert (alns2pair[0].size() != 0 || alns2pair[1].size() != 0);
            if (alns2pair[0].size() > 0) { fillTidBegEndFromAlns2(tid2, beg2, end2, alns2pair[0], initialized); initialized = true; }
            if (alns2pair[1].size() > 0) { fillTidBegEndFromAlns2(tid2, beg2, end2, alns2pair[1], initialized); initialized = true; }
            Symbol2CountCoverage read_duplex_amplicon(tid2, beg2, end2);
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                if (alns2.size() == 0) { continue; }
                uint32_t tid2, beg2, end2;
                fillTidBegEndFromAlns2(tid2, beg2, end2, alns2);
                Symbol2CountCoverage read_family_amplicon(tid2, beg2, end2);
                Symbol2CountCoverage read_family_con_ampl(tid2, beg2, end2);
                for (const auto & aln_vec : alns2) {
                    uint32_t tid1, beg1, end1;
                    fillTidBegEndFromAlns1(tid1, beg1, end1, aln_vec);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid1, beg1, end1);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln(
                            aln_vec,
                            symbolType2addPhred, 
                            frag_indel_basemax, 
                            is_proton, 
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec, 
                            region_repeatvec,
                            0);
                    // The line below is similar to : read_family_amplicon.updateByConsensus<SYMBOL_COUNT_SUM>(read_ampBQerr_fragWithR1R2);
                    read_family_con_ampl.template updateByConsensus</*SYMBOL_COUNT_SUM,*/ true>(read_ampBQerr_fragWithR1R2);
                    //if (is_proton) { 
                    read_family_amplicon.updateByMajorMinusMinor(read_ampBQerr_fragWithR1R2); 
                    //}
                }
                if ((2 == alns2pair2dflag.second) && alns2pair[0].size() > 0 && alns2pair[1].size() > 0) { // is duplex
                    read_duplex_amplicon.template updateByConsensus</*SYMBOL_COUNT_SUM,*/ false>(read_family_amplicon);
                }
                std::basic_string<std::pair<unsigned int, AlignmentSymbol>> pos_symbol_string;
                //std::vector<std::array<AlignmentSymbol, NUM_SYMBOL_TYPES>> con_symbols_vec(
                //        read_family_amplicon.getExcluEndPosition() - read_family_amplicon.getIncluBegPosition(),
                //        std::array<AlignmentSymbol, NUM_SYMBOL_TYPES>({END_ALIGNMENT_SYMBOLS, END_ALIGNMENT_SYMBOLS}));
                for (size_t epos = read_family_amplicon.getIncluBegPosition(); 
                        epos < read_family_amplicon.getExcluEndPosition(); 
                        epos++) {
                    for (SymbolType symbolType : SYMBOL_TYPES_IN_VCF_ORDER) {
                        AlignmentSymbol con_symbol; // = END_ALIGNMENT_SYMBOLS;
                        unsigned int con_count, tot_count;
                        read_family_amplicon.getRefByPos(epos).fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        if (0 == tot_count) { continue; }
                        
                        unsigned int majorcount = this->symbol_to_fam_format_depth_sets_2strand[strand].getByPos(epos).at(con_symbol).at(FAM_cDPM);
                        unsigned int minorcount = this->symbol_to_fam_format_depth_sets_2strand[strand].getByPos(epos).at(con_symbol).at(FAM_cDPm);
                        // unsigned int majorcount = this->major_amplicon[strand].getByPos(epos).getSymbolCount(con_symbol);
                        // unsigned int minorcount = this->minor_amplicon[strand].getByPos(epos).getSymbolCount(con_symbol);
                        unsigned int confam_qual = MAX(con_count * 2, tot_count + 1) - tot_count;
                        if (LINK_SYMBOL == symbolType) {
                            unsigned int con_count = read_family_con_ampl.getByPos(epos).getSymbolCount(con_symbol);
                            unsigned int tot_count = read_family_con_ampl.getByPos(epos).sumBySymbolType(symbolType);
                            double effective_nreads = log((double)(con_count + 1) / (double)(tot_count - con_count + 1)) / log(2); // PCR errors
                            double prior_weight = 1.0 / (minorcount + 1.0);
                            double realphred = prob2realphred((minorcount + prior_weight) / (majorcount + minorcount + prior_weight / phred2prob(confam_qual)));
                            confam_qual = (unsigned int)MIN(effective_nreads * realphred, phred_indel_error_before_barcode_labeling + realphred);
                        } else if (is_proton) {
                            double effective_nreads = (double)read_family_con_ampl.getByPos(epos).getSymbolCount(con_symbol);
                            double prior_weight = 1.0 / (minorcount + 1.0);
                            double realphred = prob2realphred((minorcount + prior_weight) / (majorcount + minorcount + prior_weight / phred2prob(confam_qual)));
                            confam_qual = (unsigned int)(effective_nreads * realphred);
                        }
                        
                        AlignmentSymbol ref_symbol = region_symbolvec[epos - this->getUnifiedIncluBegPosition()]; 
                        if (areSymbolsMutated(ref_symbol, con_symbol) && (LINK_SYMBOL == symbolType || confam_qual >= 20)) {
                            pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
                        }
                        unsigned int max_qual = phred_max_table.toPhredErrRate(ref_symbol, con_symbol);
                        unsigned int confam_qual2 = MIN(confam_qual, max_qual);
                        
                        int pbucket = (max_qual - confam_qual2 + 2) / 4;
                        this->dedup_ampDistr[strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pbucket, 1);
                        
                        // this->symbol_to_fam_format_depth_sets[strand].getRefByPos(epos).incSymbol

                        //if (BASE_N == con_symbol) { phredlike = MIN(phredlike, phred_thres); }
                        if (isSymbolIns(con_symbol)) {
                            posToIndelToCount_updateByConsensus(
                                    this->symbol_to_fam_format_depth_sets_2strand[strand].getRefPosToIseqToData(con_symbol),
                                    read_family_amplicon.getPosToIseqToData(con_symbol), epos, 1);
                        }
                        if (isSymbolDel(con_symbol)) {
                            posToIndelToCount_updateByConsensus(
                                    this->symbol_to_fam_format_depth_sets_2strand[strand].getRefPosToDlenToData(con_symbol),
                                    read_family_amplicon.getPosToDlenToData(con_symbol), epos, 1);
                        }
                    }
                }
                
                if (pos_symbol_string.size() > 1) {
                    mutform2count4map.insert(std::make_pair(pos_symbol_string, std::array<unsigned int, 2>({0, 0})));
                    mutform2count4map[pos_symbol_string][strand]++;
                }
            }
            if ((2 == alns2pair2dflag.second) && alns2pair[0].size() > 0 && alns2pair[1].size() > 0) { // is duplex
                for (size_t epos = read_duplex_amplicon.getIncluBegPosition(); epos < read_duplex_amplicon.getExcluEndPosition(); epos++) {
                    for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+(unsigned int)symbolType)) {
                        AlignmentSymbol con_symbol;
                        unsigned int con_count, tot_count;
                        read_duplex_amplicon.getRefByPos(epos).fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        assert (tot_count <= 2 || !fprintf(stderr, "%d <= 2 failed for duplex family, a duplex family is supported by two single-strand families!\n", tot_count));
                        if (0 < tot_count) {
                            //this->duplex_tsum_depth.getRefByPos(epos).incSymbolCount(con_symbol, 1);
                            this->symbol_to_duplex_format_depth_sets.getRefByPos(epos)[con_symbol][DUPLEX_dDP1] += 1;
                        }
                        if (1 < tot_count) {
                            //this->duplex_pass_depth.getRefByPos(epos).incSymbolCount(con_symbol, 1);
                            this->symbol_to_duplex_format_depth_sets.getRefByPos(epos)[con_symbol][DUPLEX_dDP2] += 1;
                        }
                    }
                }
            }
        }
        for (unsigned int strand = 0; strand < 2; strand++) {
            assert(this->symbol_to_fam_format_depth_sets_2strand[strand].getExcluEndPosition() == dedup_ampDistr[strand].getExcluEndPosition());
            assert(this->symbol_to_fam_format_depth_sets_2strand[strand].getIncluBegPosition() == dedup_ampDistr[strand].getIncluBegPosition()); 
            auto VQ_cIAQ = (strand ? VQ_cIAQr : VQ_cIAQf);
            auto VQ_cIAD = (strand ? VQ_cIADr : VQ_cIADf);
            auto VQ_cIDQ = (strand ? VQ_cIDQr : VQ_cIDQf);
            for (unsigned int epos = this->getUnifiedIncluBegPosition(); epos < this->getUnifiedExcluEndPosition(); epos++) {
                for (SymbolType symbolType : SYMBOL_TYPE_ARR) {
                    AlignmentSymbol ref_symbol = region_symbolvec[epos - this->getUnifiedIncluBegPosition()];
                    const auto totDP = formatSumBySymbolType(this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos), symbolType, FAM_cDP1);
                    for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symbolType]) {
                        const unsigned int max_qual = phred_max_table.toPhredErrRate(ref_symbol, symbol);
                        int maxvqual = 0;
                        unsigned int argmaxAD = 0;
                        unsigned int argmaxBQ = 0;
                        AlignmentSymbol refsymbol = region_symbolvec[epos - this->getUnifiedIncluBegPosition()];                   
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
                    }
                }
                this->dedup_ampDistr[strand].getRefByPos(epos).clearSymbolBucketCount();
            }
        }
    };
    
    int 
    updateHapMap(std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map, 
            const auto & tsum_depth, const auto read_count_field_id, unsigned int max_ploidy = 2+1) {
        for (auto it = mutform2count4map.begin(); it != mutform2count4map.end();) {
            std::basic_string<std::pair<unsigned int, AlignmentSymbol>> mutform = it->first;
            auto counts = it->second;
            std::vector<unsigned int> dsADs; 
            dsADs.reserve(mutform.size() + 1);
            dsADs.push_back(0);
            for (std::pair<unsigned int, AlignmentSymbol> simplemut : mutform) {
                unsigned int ad0 = tsum_depth.at(0).getByPos(simplemut.first)[simplemut.second][read_count_field_id];
                unsigned int ad1 = tsum_depth.at(1).getByPos(simplemut.first)[simplemut.second][read_count_field_id];
                dsADs.push_back(ad0 + ad1);
            }
            auto med = MEDIAN(dsADs);
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
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map_bq,
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map_fq,
            const std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> & alns3, 
            const std::string & refstring,
            const std::vector<RegionalTandemRepeat> & region_repeatvec,
            unsigned int bq_phred_added_misma, 
            unsigned int bq_phred_added_indel, 
            bool should_add_note, 
            // const unsigned int frag_indel_ext, 
            const unsigned int frag_indel_basemax,
            const PhredMutationTable & phred_max_sscs_table, 
            //unsigned int phred_thres,
            //const double ess_georatio_dedup, 
            //const double ess_georatio_duped_pcr,
            bool use_deduplicated_reads, 
            // bool is_loginfo_enabled, 
            // unsigned int thread_id, 
            // unsigned int fixedthresBQ, 
            // unsigned int nogap_phred,
            // unsigned int highqual_thres_snv, 
            // unsigned int highqual_thres_indel, 
            // unsigned int uni_bias_r_max,
            const bool is_proton, 
            // const AssayType assay_type, // unused?
            const unsigned int phred_indel_error_before_barcode_labeling,
            const unsigned int baq_per_aligned_base,
            // const unsigned int regside_nbases,
            // const uint32_t bias_flag_snv,
            // const uint32_t bias_flag_indel,
            const unsigned int specialflag = 0) {
        
        const std::array<unsigned int, NUM_SYMBOL_TYPES> symbolType2addPhred = {{bq_phred_added_misma, bq_phred_added_indel}};
        std::basic_string<AlignmentSymbol> ref_symbol_string = string2symbolseq(refstring);
        
        updateByAlns3UsingBQ(
                mutform2count4map_bq, 
                alns3, 
                ref_symbol_string, 
                symbolType2addPhred, 
                should_add_note, 
                frag_indel_basemax, 
                phred_max_sscs_table, 
                is_proton, 
                region_repeatvec, 
                0); // base qualities
        updateHapMap(mutform2count4map_bq, this->symbol_to_frag_format_depth_sets, FRAG_bDP);
        if (use_deduplicated_reads) {
            updateByAlns3UsingFQ(
                    mutform2count4map_fq, 
                    alns3,
                    ref_symbol_string, 
                    symbolType2addPhred, 
                    should_add_note, 
                    frag_indel_basemax, 
                    phred_max_sscs_table, 
                    is_proton, 
                    region_repeatvec,
                    phred_indel_error_before_barcode_labeling,
                    0); // family qualities
            updateHapMap(mutform2count4map_fq, this->symbol_to_fam_format_depth_sets_2strand, FAM_cDP1);
        }
        return 0;
    };
};

// fmt.DDP1,  fmt.dDP,   symbol_to_duplex_format_depth_sets, DUPLEX_dDP1, refpos, symbolType, refsymbol
void 
fill_symboltype_fmt(
        auto & fmtDP,
        auto & fmtAD,
        const auto & symbol_to_abcd_format_depth_sets,
        const auto format_field,
        const auto refpos,
        const SymbolType symbolType, 
        const AlignmentSymbol refsymbol) {
    const auto symbolNN = SYMBOL_TYPE_TO_AMBIG[symbolType];
    fmtDP[0] = formatSumBySymbolType(symbol_to_abcd_format_depth_sets.getByPos(refpos), symbolType, format_field);
    fmtDP[1] = symbol_to_abcd_format_depth_sets.getByPos(refpos)[symbolNN][format_field];
    // clear_push(fmtAD, symbol_to_abcd_format_depth_sets.getByPos(refpos)[refsymbol][format_field]);
};
// fmt.aDPff, symbol_to_seg_format_depth_sets, SEG_aDPff, refpos, symbo
void
fill_symbol_fmt(
        auto & fmtAD,
        const auto & symbol_to_abcd_format_depth_sets,
        const auto format_field,
        const auto refpos,
        const auto symbol, //,
        // const auto cnt, 
        //const unsigned int newsize = 1,
        const unsigned int allele_idx
        ) {
    //fmtAD.resize(newsize);
    // assert(fmtAD.size() == 1);
    if (0 == allele_idx) {
        fmtAD.clear();
    } else {
        assert(fmtAD.size() > 0);
    }
    auto cnt = symbol_to_abcd_format_depth_sets.getByPos(refpos)[symbol][format_field];
    fmtAD.push_back(cnt);
};

int
fill_symbol_VQ_fmts(
        auto & fmt,
        const auto & symbol_to_VQ_format_tag_sets,
        const size_t refpos,
        const AlignmentSymbol symbol,
        const unsigned int specialflag) {
    const unsigned int a = 0;

    fill_symbol_fmt(fmt.aSBQf, symbol_to_VQ_format_tag_sets, VQ_aSBQf, refpos, symbol, a);
    fill_symbol_fmt(fmt.aSBQr, symbol_to_VQ_format_tag_sets, VQ_aSBQr, refpos, symbol, a);
    
    int a_BQ_syserr_qual_fw = (int)(fmt.aSBQf[a] * 2 + fmt.aSBQr[a]) - (int)(22 * ((fmt.aDPff[a] + fmt.aDPrf[a]) * 2 + (fmt.aDPfr[a] + fmt.aDPrr[a])));
    int a_BQ_syserr_qual_rv = (int)(fmt.aSBQf[a] + fmt.aSBQr[a] * 2) - (int)(22 * ((fmt.aDPff[a] + fmt.aDPrf[a]) + (fmt.aDPfr[a] + fmt.aDPrr[a]) * 2));
    int a_BQ_avg_qual = (int)(fmt.aSBQf[a] + fmt.aSBQr[a]) - (int)(22 * (fmt.aDPff[a] + fmt.aDPrf[a] + fmt.aDPfr[a] + fmt.aDPrr[a]));
    int a_BQ_syserr_qual = MAX3(a_BQ_syserr_qual_fw, a_BQ_syserr_qual_rv, a_BQ_avg_qual);
    clear_push(fmt.aBQQ, (MAX(a_BQ_syserr_qual_fw, a_BQ_syserr_qual_rv)), a);
    fill_symbol_fmt(fmt.bMQ,  symbol_to_VQ_format_tag_sets,  VQ_bMQ,  refpos, symbol, a);
    fmt.bMQ[a] = fmt.bMQ[a] / SQR_QUAL_DIV / MAX(fmt.bDPf[a] + fmt.bDPr[a], 1);
    
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
    clear_push(fmt.VTI,  (unsigned int)symbol, a);
}

std::array<unsigned int, 2>
BcfFormat_symboltype_init(bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet12, 
        unsigned int refpos, 
        const SymbolType symbolType,
        bool use_deduplicated_reads, 
        const AlignmentSymbol refsymbol,
        const unsigned int specialflag) {
    
    const auto symbolNN = SYMBOL_TYPE_TO_AMBIG[symbolType]; 
    const auto & symbol_to_seg_format_depth_sets = symbol2CountCoverageSet12.symbol_to_seg_format_depth_sets;
    
    fill_symboltype_fmt(fmt.ADPff, fmt.aDPff, symbol_to_seg_format_depth_sets, SEG_aDPff, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ADPfr, fmt.aDPfr, symbol_to_seg_format_depth_sets, SEG_aDPfr, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ADPrf, fmt.aDPrf, symbol_to_seg_format_depth_sets, SEG_aDPrf, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ADPrr, fmt.aDPrr, symbol_to_seg_format_depth_sets, SEG_aDPrr, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.AEP1, fmt.aEP1, symbol_to_seg_format_depth_sets, SEG_aEP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.AEP2, fmt.aEP2, symbol_to_seg_format_depth_sets, SEG_aEP2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.AXM1, fmt.aXM1, symbol_to_seg_format_depth_sets, SEG_aXM1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.AXM2, fmt.aXM2, symbol_to_seg_format_depth_sets, SEG_aXM2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALI1, fmt.aLI1, symbol_to_seg_format_depth_sets, SEG_aLI1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALI2, fmt.aLI2, symbol_to_seg_format_depth_sets, SEG_aLI2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARI1, fmt.aRI1, symbol_to_seg_format_depth_sets, SEG_aRI1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARI2, fmt.aRI2, symbol_to_seg_format_depth_sets, SEG_aRI2, refpos, symbolType, refsymbol);
    
    const auto & symbol_to_frag_format_depth_sets = symbol2CountCoverageSet12.symbol_to_frag_format_depth_sets;
    
    fill_symboltype_fmt(fmt.BDPf,  fmt.bDPf,  symbol_to_frag_format_depth_sets[0], FRAG_bDP, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.BDPr,  fmt.bDPr,  symbol_to_frag_format_depth_sets[1], FRAG_bDP, refpos, symbolType, refsymbol);
    //fill_symboltype_fmt(fmt.BDPBQ, fmt.bDPBQ, symbol_to_frag_format_depth_sets, FRAG_bBQ,  refpos, symbolType, refsymbol);
    //fill_symboltype_fmt(fmt.BDPMQ, fmt.bDPMQ, symbol_to_frag_format_depth_sets, FRAG_bMQ,  refpos, symbolType, refsymbol);
    
    const auto & symbol_to_fam_format_depth_sets = symbol2CountCoverageSet12.symbol_to_fam_format_depth_sets_2strand;
    
    fill_symboltype_fmt(fmt.CDP1f, fmt.cDP1f, symbol_to_fam_format_depth_sets[0], FAM_cDP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP2f, fmt.cDP1f, symbol_to_fam_format_depth_sets[0], FAM_cDP2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP3f, fmt.cDP2f, symbol_to_fam_format_depth_sets[0], FAM_cDP3, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.C1DPf, fmt.cDP3f, symbol_to_fam_format_depth_sets[0], FAM_c1DP, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPMf, fmt.cDPMf, symbol_to_fam_format_depth_sets[0], FAM_cDPM, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPmf, fmt.cDPmf, symbol_to_fam_format_depth_sets[0], FAM_cDPm, refpos, symbolType, refsymbol);

    fill_symboltype_fmt(fmt.CDP1r, fmt.cDP1r, symbol_to_fam_format_depth_sets[1], FAM_cDP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP2r, fmt.cDP1r, symbol_to_fam_format_depth_sets[1], FAM_cDP2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP3r, fmt.cDP2r, symbol_to_fam_format_depth_sets[1], FAM_cDP3, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.C1DPr, fmt.cDP3r, symbol_to_fam_format_depth_sets[1], FAM_c1DP, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPMr, fmt.cDPMr, symbol_to_fam_format_depth_sets[1], FAM_cDPM, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPmr, fmt.cDPmr, symbol_to_fam_format_depth_sets[1], FAM_cDPm, refpos, symbolType, refsymbol);

    const auto & symbol_to_duplex_format_depth_sets = symbol2CountCoverageSet12.symbol_to_duplex_format_depth_sets;
    fill_symboltype_fmt(fmt.DDP1,  fmt.dDP1,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.DDP2,  fmt.dDP2,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP2, refpos, symbolType, refsymbol);
    
    //fill_symboltype_fmt(fmt.AABQf, fmt.aSBQf, symbol_to_VQ_format_tag_sets, VQ_ASBQf, refpos, symbolType, refsymbol);
    //fill_symboltype_fmt(fmt.AABQr, fmt.aSBQr, symbol_to_VQ_format_tag_sets, VQ_ASBQr, refpos, symbolType, refsymbol);
    //fmt.AABQf /= MAX(1, fmt.ADPff);
    
    fmt.gapNf.clear();
    fmt.gapNr.clear();
    fmt.gapSeq.clear();
    fmt.gapbAD1.clear();
    fmt.gapcAD1.clear();
    fmt.bHap.clear();
    fmt.cHap.clear();
    return {fmt.BDPf[0] + fmt.BDPr[0], fmt.CDP1f[0] + fmt.CDP1r[0]};
};


std::array<unsigned int, 2>
BcfFormat_symbol_init(
        bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet12, 
        unsigned int refpos, 
        const AlignmentSymbol symbol,
        const auto & mutform2count4vec_bq,
        const auto & indices_bq,
        const auto & mutform2count4vec_fq,
        const auto & indices_fq,
        const unsigned int bDPa,
        const unsigned int cDP1a,
        const std::string & gapSa,
        const unsigned int specialflag) {
    const unsigned int a = 0;
    
    const auto & symbol_to_seg_format_depth_sets = symbol2CountCoverageSet12.symbol_to_seg_format_depth_sets;
    
    fill_symbol_fmt(fmt.aDPff, symbol_to_seg_format_depth_sets, SEG_aDPff, refpos, symbol, a);
    fill_symbol_fmt(fmt.aDPfr, symbol_to_seg_format_depth_sets, SEG_aDPfr, refpos, symbol, a);
    fill_symbol_fmt(fmt.aDPrf, symbol_to_seg_format_depth_sets, SEG_aDPrf, refpos, symbol, a);
    fill_symbol_fmt(fmt.aDPrr, symbol_to_seg_format_depth_sets, SEG_aDPrr, refpos, symbol, a);
    
    //fill_symbol_fmt(fmt.aBQ1, symbol_to_seg_format_depth_sets, SEG_aBQ1, refpos, symbol);
    //fill_symbol_fmt(fmt.aBQ2, symbol_to_seg_format_depth_sets, SEG_aBQ2, refpos, symbol);
    fill_symbol_fmt(fmt.aEP1, symbol_to_seg_format_depth_sets, SEG_aEP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aEP2, symbol_to_seg_format_depth_sets, SEG_aEP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aXM1, symbol_to_seg_format_depth_sets, SEG_aXM1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aXM2, symbol_to_seg_format_depth_sets, SEG_aXM2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLI1, symbol_to_seg_format_depth_sets, SEG_aLI1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLI2, symbol_to_seg_format_depth_sets, SEG_aLI2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRI1, symbol_to_seg_format_depth_sets, SEG_aRI1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRI2, symbol_to_seg_format_depth_sets, SEG_aRI2, refpos, symbol, a);
    
    const auto & symbol_to_frag_format_depth_sets = symbol2CountCoverageSet12.symbol_to_frag_format_depth_sets;
    
    fill_symbol_fmt(fmt.bDPf, symbol_to_frag_format_depth_sets[0], FRAG_bDP, refpos, symbol, a);
    fill_symbol_fmt(fmt.bDPr, symbol_to_frag_format_depth_sets[1], FRAG_bDP, refpos, symbol, a);
    // fill_symbol_fmt(fmt.b2BQ, symbol_to_frag_format_depth_sets, FRAG_b2BQ, refpos, symbol);
    
    const auto & symbol_to_fam_format_depth_sets = symbol2CountCoverageSet12.symbol_to_fam_format_depth_sets_2strand;
    fill_symbol_fmt(fmt.cDP1f, symbol_to_fam_format_depth_sets[0], FAM_cDP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP2f, symbol_to_fam_format_depth_sets[0], FAM_cDP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP3f, symbol_to_fam_format_depth_sets[0], FAM_cDP3, refpos, symbol, a);
    fill_symbol_fmt(fmt.c1DPf, symbol_to_fam_format_depth_sets[0], FAM_c1DP, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPMr, symbol_to_fam_format_depth_sets[0], FAM_cDPM, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPmr, symbol_to_fam_format_depth_sets[0], FAM_cDPm, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP1r, symbol_to_fam_format_depth_sets[1], FAM_cDP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP2r, symbol_to_fam_format_depth_sets[1], FAM_cDP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP3r, symbol_to_fam_format_depth_sets[1], FAM_cDP3, refpos, symbol, a);
    fill_symbol_fmt(fmt.c1DPr, symbol_to_fam_format_depth_sets[1], FAM_c1DP, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPMr, symbol_to_fam_format_depth_sets[1], FAM_cDPM, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPmr, symbol_to_fam_format_depth_sets[1], FAM_cDPm, refpos, symbol, a);
    
    const auto & symbol_to_duplex_format_depth_sets = symbol2CountCoverageSet12.symbol_to_duplex_format_depth_sets;
    fill_symbol_fmt(fmt.dDP1,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.dDP2,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP2, refpos, symbol, a);
    
    const auto & symbol_to_VQ_format_tag_sets = symbol2CountCoverageSet12.symbol_to_VQ_format_tag_sets;
    
    fill_symbol_VQ_fmts(
            fmt,
            symbol_to_VQ_format_tag_sets,
            refpos,
            symbol,
            a);
    if (a > 0) {
        fmt.bHap = mutform2count4map_to_phase(mutform2count4vec_bq, indices_bq);
        fmt.cHap = mutform2count4map_to_phase(mutform2count4vec_fq, indices_fq);
    }
    fmt.bDPa.push_back(bDPa);
    fmt.cDP1a.push_back(cDP1a);
    fmt.gapSa.push_back(gapSa);
    return std::array<unsigned int, 2>({{fmt.bDPf[a] + fmt.bDPr[a], fmt.cDP1f[a] + fmt.cDP1r[a]}});
    /*
    std::vector<std::map<std::string, uint32_t>> iseq2cnt_vec;
    std::vector<std::map<uint32_t,    uint32_t>> dlen2cnt_vec;
    for (size_t strand = 0; strand < 2; strand++) {
        for (const auto ins_symb : INS_SYMBOLS) {
            auto & m = symbol2CountCoverageSet12.fq_tsum_depth.at(strand).getPosToIseqToData(ins_symb);
            if (m.find(refpos) != m.end()) {
                iseq2cnt_vec.push_back(m.at(refpos));
            }
        }
        for (const auto del_symb : DEL_SYMBOLS) {
            auto & m = symbol2CountCoverageSet12.fq_tsum_depth.at(strand).getPosToDlenToData(del_symb);
            if (m.find(refpos) != m.end()) {
                dlen2cnt_vec.push_back(m.at(refpos));
            }
        }
    }*/
};

int
BcfFormat_symbol_calc_DPv(
        bcfrec::BcfFormat & fmt,
        //unsigned int bdepth_allele,
        //unsigned int cdepth_allele,
        unsigned int specialflag = 0) {
    
    const unsigned int a = 0;
    //clear_push(fmt.bDPa, bdepth_allele, a);
    //clear_push(fmt.cDPa, cdepth_allele, a);
    
    // Non-UMI universality-based prequal allele fraction
    unsigned int ADP = fmt.ADPff[0] + fmt.ADPfr[0] + fmt.ADPrf[0] + fmt.ADPrr[0];
    unsigned int aDP = fmt.aDPff[a] + fmt.aDPfr[a] + fmt.aDPrf[a] + fmt.aDPrr[a];
    double aEPFA = (fmt.aEP1[a] + 0.5) / (ADP - fmt.aEP1[a] - fmt.aEP2[a] - (aDP - fmt.aEP1[a]) + 1.0);
    double aXMFA = (fmt.aXM1[a] + 0.5) / (ADP - fmt.aXM1[a] - fmt.aXM2[a] - (aDP - fmt.aXM1[a]) + 1.0);
    double aLIFA = (fmt.aLI1[a] + 0.5) / (ADP - fmt.aLI1[a] - fmt.aLI2[a] - (aDP - fmt.aLI1[a]) + 1.0);
    double aRIFA = (fmt.aRI1[a] + 0.5) / (ADP - fmt.aRI1[a] - fmt.aRI2[a] - (aDP - fmt.aRI1[a]) + 1.0);
    
    double aSSFA = dp4_to_pcFA(fmt.ADPff[0] + fmt.ADPrf[0], fmt.ADPfr[0] + fmt.ADPrr[0], fmt.aDPff[a] + fmt.aDPrf[a], fmt.aDPrf[a] + fmt.aDPff[a]);
    double aROFA = dp4_to_pcFA(fmt.ADPff[0] + fmt.ADPrr[0], fmt.ADPfr[0] + fmt.ADPrf[0], fmt.aDPff[a] + fmt.aDPrr[a], fmt.aDPfr[a] + fmt.aDPrf[a]);
    
    double bFA = (fmt.bDPa[a] + 0.5) / (fmt.BDPf[0] + fmt.BDPr[0] + 1.0);
    double cFA = (fmt.cDP1a[a] + 0.5) / (fmt.CDP1f[0] + fmt.CDP1r[0] + 1.0);
    
    double min_aFA = MINVEC(std::array<double, 6>{{aEPFA, aXMFA, aLIFA, aRIFA, aSSFA, aROFA}});
    double min_bcFA = MIN3(min_aFA, bFA, cFA);
    
    // UMI universality-based prequal allele fraction
    double cFA2 = ((fmt.cDP2f[a] + fmt.cDP2r[a]) + 0.5) / (fmt.CDP2f[0] + fmt.CDP2r[0] + 1.0);
    double cFA3 = ((fmt.cDP3f[a] + fmt.cDP3r[a]) + 0.5) / (fmt.CDP3f[0] + fmt.CDP3r[0] + 1.0);
    double umi_cFA = MIN3(min_bcFA, cFA2, cFA3);
    
    clear_push(fmt.cDP1v, (int)(min_bcFA * (fmt.CDP1f[0] +fmt.cDP1r[0]) * 100), a);
    clear_push(fmt.cDP2v, (int)(umi_cFA * (fmt.CDP2f[0] +fmt.cDP2r[0]) * 100), a);
    return 0;
};

int
BcfFormat_symbol_sum_DPv(auto & fmts) {
    unsigned int CDP1v = 0;
    unsigned int CDP2v = 0;
    for (auto & fmt : fmts) {
        CDP1v += std::get<0>(fmt).cDP1v[0];
        CDP2v += std::get<0>(fmt).cDP2v[0];
    }
    for (auto & fmt : fmts) {
        std::get<0>(fmt).CDP1v[0] = CDP1v;
        std::get<0>(fmt).CDP2v[0] = CDP2v;
    }
    return 0;
};

int
BcfFormat_symbol_calc_qual(
        bcfrec::BcfFormat & fmt,
        //const AlignmentSymbol symbol,
        //const unsigned int indelbdepth,
        //const unsigned int indelcdepth,
        //const std::string & indelstring,
        const std::string & repeatunit,
        const unsigned int repeatnum,
        // const unsigned int a = 0,
        double powlaw_anyvar_base = 90,
        double powlaw_exponent = 3,
        double powlaw_sscs_inc = 30,
        const unsigned int phred_varcall_err_per_map_err_per_base = 20,
        const unsigned int specialflag = 0) {
    
    const unsigned int a = 0;
    int duped_frag_binom_qual = fmt.bIAQb[a];
    
    int sscs_binom_qual_fw = fmt.cIAQf[a] + fmt.cIAQr[a] * MIN(60 - fmt.cIDQf[a], fmt.cIDQr[a]) / fmt.cIDQr[a];
    int sscs_binom_qual_rv = fmt.cIAQr[a] + fmt.cIAQf[a] * MIN(60 - fmt.cIDQr[a], fmt.cIDQf[a]) / fmt.cIDQf[a];
    int sscs_binom_qual = MAX(sscs_binom_qual_fw, sscs_binom_qual_rv); 
    
    double min_bcFA = (((double)(fmt.cDP1v[a])) / ((double)(fmt.CDP1v[0]))); 
    int dedup_frag_powlaw_qual = (unsigned int)(powlaw_exponent * 10.0 / log(10.0) * log(min_bcFA) + powlaw_anyvar_base + 10);
    double umi_cFA = (((double)(fmt.cDP2v[a])) / ((double)(fmt.CDP2v[0]))); 
    int sscs_powlaw_qual = (int)(powlaw_exponent * 10.0 / log(10.0) * log(umi_cFA) + powlaw_anyvar_base + powlaw_sscs_inc);
    
    const std::string & indelstring = fmt.gapSa[a];
    if (indelstring.size() > 0) {
        // const double symbol_to_allele_frac = 1.0 - pow((isSymbolIns(symbol) ? 0.9 : (isSymbolDel(symbol) ? 0.95 : 1.0)), indelstring.size());
        // fmt.cADR[1] * ((double)(fmt.gapDP4[2] + 1) / (double)(fmt.gapDP4[0] + 1));
        // const double tDP0a = (fmt.cDP1f[a] + fmt.cDP1r[a]) * (double)(indelbdepth + 0.5) / (double)(fmt.bDPf[a] + fmt.bDPr[a] + 1);
        const double indel_pq = MIN(indel_phred(8.0, indelstring.size(), repeatunit.size(), repeatnum), 24) + 2 - 10; // phred_snv_to_indel_ratio;
        const double indel_ic = 10.0 / log(10.0) * log((double)MAX(indelstring.size(), 1U) / (double)(repeatunit.size() * (MAX(1, repeatnum) - 1) + 1));
        //const double indel_penal4multialleles = log((double)(1 + fmt.CDP1f[0] + fmt.CDP1r[0] - fmt.cDP1f[0] - fmt.cDP1r[0]) / (indelcdepth + 1.0)) / log(2.0) * 8.0;
        dedup_frag_powlaw_qual += (int)(indel_ic); // - indel_penal4multialleles);
        duped_frag_binom_qual  -= (int)(indel_pq); // - indel_penal4multialleles);
    }

    clear_push(fmt.bIAQ, duped_frag_binom_qual, a); 
    clear_push(fmt.cIAQ, sscs_binom_qual, a);
    
    clear_push(fmt.cPLQ1, dedup_frag_powlaw_qual, a); // phred_varcall_err_per_map_err_per_base
    clear_push(fmt.cPLQ2, sscs_powlaw_qual, a);
        
    const int syserr_q = MIN(fmt.aBQQ[a], fmt.bMQ[a] + phred_varcall_err_per_map_err_per_base);
    clear_push(fmt.cVQ1, MIN3(syserr_q, fmt.bIAQ[a], fmt.cPLQ1[a]), a);
    clear_push(fmt.cVQ2, MIN(syserr_q, MIN(fmt.cIAQ[a], fmt.cPLQ2[a])), a);
    // TODO; check if reducing all allele read count to increase alt allele frac in case of ref bias makes more sense
};
    
/*      
    fmt.bHap = mutform2count4map_to_phase(mutform2count4vec_bq, indices_bq);
    fmt.cHap = mutform2count4map_to_phase(mutform2count4vec_fq, indices_fq);
    bool is_novar; // = (symbol == LINK_M || (isSymbolSubstitution(symbol) && vcfref == vcfalt));
    
    // auto fmtAD = SIGN2UNSIGN(0);
    // fmt.MQ = (unsigned int)sqrt((double)bq_qsum_sqrMQ_tot / (DBL_MIN + (double)(fmt.bAD1[0] + fmt.bAD1[1])));
    
        /*fill_TN_germline(
                fmt, 
                symbol,
                powlaw_exponent,
                refmul,
                altmul,
                any_mul_contam_frac,
                ref_mul_contam_frac, 
                t2n_add_contam_transfrac,
                prev_is_tumor, 
                tki, 
                is_novar,
                somaticGT,
                0); 
    } else {
        fmt.GT = "./.";
        fmt.GQ = 0;
        fmt.RefBias = 0;
        fmt.GTa = ".";
        fmt.GTb = ".";
        fmt.GQa = 0;
        fmt.GQb = 0;
        fmt.GLa = {0};
        fmt.GLb = {0};
        fmt.GSTa = {0};
        fmt.GSTb = {0};
    }
    *
    fmt.HQ[0] = 0; 
    fmt.HQ[1] = 0;
    

* */

// END of new

#define INDEL_ID 1
#include "instcode.hpp"
#undef INDEL_ID
#define INDEL_ID 2
#include "instcode.hpp"
std::array<unsigned int, 2>
fill_by_indel_info(
        bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet, 
        const unsigned int strand, 
        const unsigned int refpos, 
        const AlignmentSymbol symbol, 
        const std::string & refstring, 
        const std::string & repeatunit, 
        unsigned int repeatnum) {
    assert(isSymbolIns(symbol) || isSymbolDel(symbol));
    if (isSymbolIns(symbol)) {
        return fill_by_indel_info2_1(fmt, symbol2CountCoverageSet, strand, refpos, symbol,
                symbol2CountCoverageSet.symbol_to_frag_format_depth_sets.at(strand).getPosToIseqToData(symbol),
                symbol2CountCoverageSet.symbol_to_fam_format_depth_sets_2strand.at(strand).getPosToIseqToData(symbol),
                refstring, repeatunit, repeatnum);
    } else {
        return fill_by_indel_info2_2(fmt, symbol2CountCoverageSet, strand, refpos, symbol,
                symbol2CountCoverageSet.symbol_to_frag_format_depth_sets.at(strand).getPosToDlenToData(symbol),
                symbol2CountCoverageSet.symbol_to_fam_format_depth_sets_2strand.at(strand).getPosToDlenToData(symbol),
                refstring, repeatunit, repeatnum);
    }
};

std::string 
mutform2count4map_to_phase(const auto & mutform2count4vec, const auto & indices, unsigned int pseudocount = 1) {
    std::string phase_string;
    for (auto idx : indices) {
        auto mutform2count4pair = mutform2count4vec.at(idx);
        auto counts = mutform2count4pair.second;
        if ((counts[0] + counts[1]) > pseudocount) {
            phase_string +="(";
            for (auto pos2symbol4pair : mutform2count4pair.first) {
                AlignmentSymbol symbol = pos2symbol4pair.second;
                unsigned int mutpos = pos2symbol4pair.first + (isSymbolSubstitution(symbol) ? 1 : 0);
                phase_string += std::string("(") + std::to_string(mutpos) + "&" + SYMBOL_TO_DESC_ARR[symbol] + ")";
            }
            phase_string += std::string("&") + std::to_string(counts[0]) + "&" + std::to_string(counts[1]) + ")";
        }
    }
    return phase_string;
}

const std::vector<std::pair<std::array<unsigned int, 2>, std::string>> 
indel_get_majority(
        const bcfrec::BcfFormat & fmt, 
        const char *tname, 
        unsigned int refpos, 
        const AlignmentSymbol symbol, 
        bool is_warning_generated = true) {
    std::vector<std::pair<std::array<unsigned int, 2>, std::string>> indelstrings;
    if ((SUMVEC(fmt.gapNf) + SUMVEC(fmt.gapNr)) == 0) {
        if (is_warning_generated) {
            std::cerr << "Invalid indel detected (invalid mutation) : " << tname << ", " << refpos << ", " << SYMBOL_TO_DESC_ARR[symbol] << std::endl;
            std::string msg;
            bcfrec::streamAppendBcfFormat(msg, fmt);
            std::cerr << msg << "\n";
        }
        indelstrings.push_back(std::make_pair(std::array<unsigned int, 2>({{0, 0}}), SYMBOL_TO_DESC_ARR[symbol]));
    } else {
        std::map<std::string, std::array<unsigned int, 2>> indelmap;
        assert(fmt.gapSeq.size() ==  SUMVEC(fmt.gapNf) + SUMVEC(fmt.gapNr));
        assert(fmt.gapbAD1.size() ==  fmt.gapSeq.size());
        for (unsigned int i = 0; i < fmt.gapSeq.size(); i++) {
            const auto & it = indelmap.find(fmt.gapSeq[i]);
            if (it == indelmap.end()) {
                indelmap.insert(std::make_pair(fmt.gapSeq[i], std::array<unsigned int, 2>({{fmt.gapbAD1[i], fmt.gapcAD1[i]}})));
            } else {
                indelmap[fmt.gapSeq[i]][0] += fmt.gapbAD1[i];
                indelmap[fmt.gapSeq[i]][1] += fmt.gapcAD1[i];
            }
        }
        unsigned int max_bAD1 = 0;
        for (auto & indelit : indelmap) {
            UPDATE_MAX(max_bAD1, indelit.second[0]);
        }
        for (auto & indelit : indelmap) {
            if (indelit.second[0] >= (max_bAD1 + 3) / 4) {
                indelstrings.push_back(std::make_pair(indelit.second, indelit.first));
            }
        }
        std::sort(indelstrings.rbegin(), indelstrings.rend());
    }
    return indelstrings;
}

// higher allele1count (lower allele2count) results in higher LODQ if FA is above frac, meaning allele1 is more likely to be homo, and vice versa
/*
unsigned int hetLODQ(double allele1count, double allele2count, double frac, double powlaw_exponent = 3.0) {
    auto binomLODQ = (int)calc_binom_10log10_likeratio(frac, allele1count, allele2count);
    auto powerLODQ = (int)(10.0/log(10.00) * powlaw_exponent * MAX(logit2(allele1count / (frac * 2.0) + 0.5, allele2count / ((1.0 - frac) * 2.0) + 0.5), 0.0));
    return MIN(binomLODQ, powerLODQ);
}*/

unsigned int hetLODQ(double allele1count, double allele2count, double powlaw_exponent = 3.0) {
    auto binomLODQ = (int)calc_binom_10log10_likeratio(0.5, allele1count, allele2count);
    auto powerLODQ = (int)(10.0/log(10.00) * powlaw_exponent * MAX(logit2(allele1count + 0.5, allele2count + 0.5), 0.0));
    return MIN(binomLODQ, powerLODQ);
}

struct {
    template<class T1, class T2>
    bool operator()(std::pair<T1, T2> a, std::pair<T1, T2> b) const {   
        return (a.second < b.second) || (a.second == b.second && a.first < b.first);
    }
} PairSecondLess;

double 
compute_norm_ad(const bcfrec::BcfFormat *fmtp, const bool isSubst, const bool isNN = false) {
    if (isNN) {
        return (double)(fmtp->cDP1f[fmtp->cDP1f.size()-1] + fmtp->cDP1r[fmtp->cDP1r.size()-1]) / 100.0;
    } else {
        return fmtp->cDP1v[fmtp->cDP1v.size()-1] / 100.0;
    }
    // return fmtp->c1FA[1] * (double)(fmtp->CAD1f[1] + fmtp->CAD1r[1]);
}

const auto
ALODQ(const auto x) {
    return x->cVQ1[x->cVQ1.size() - 1];
};

int
output_germline(
        std::string & out_string,
        AlignmentSymbol refsymbol, 
        std::vector<std::pair<AlignmentSymbol, bcfrec::BcfFormat*>> symbol_format_vec,
        const char *tname,
        const std::string & refstring,
        unsigned int refpos,
        unsigned int extended_inclu_beg_pos, // regionpos,
        unsigned int central_readlen,
        // const auto & tki, 
        unsigned int outputflag,
        unsigned int specialflag) {
    
    assert(symbol_format_vec.size() >= 4 || 
            !fprintf(stderr, " The variant-type %s:%u %u has symbol_format_vec of length %u", 
            tname, refpos, refsymbol, symbol_format_vec.size()));
    unsigned int regionpos = refpos- extended_inclu_beg_pos;
    struct {
        bool operator()(std::pair<AlignmentSymbol, bcfrec::BcfFormat*> & p1, std::pair<AlignmentSymbol, bcfrec::BcfFormat*> & p2) const {
            return ALODQ(p1.second) < ALODQ(p2.second);
        }
    } SymbolBcfFormatPairLess;
    std::sort(symbol_format_vec.rbegin(), symbol_format_vec.rend(), SymbolBcfFormatPairLess);
    std::array<std::pair<AlignmentSymbol, bcfrec::BcfFormat*>, 4> ref_alt1_alt2_alt3 = {{std::make_pair(END_ALIGNMENT_SYMBOLS, (bcfrec::BcfFormat*)NULL)}};
    unsigned int allele_idx = 1;
    for (auto symb_fmt : symbol_format_vec) {
        if (refsymbol == symb_fmt.first) {
            ref_alt1_alt2_alt3[0] = symb_fmt;
        } else if (allele_idx <= 3) {
            ref_alt1_alt2_alt3[allele_idx] = symb_fmt;
            allele_idx++;
        }
    }
    assert(ref_alt1_alt2_alt3[0].second != NULL);
    assert(ref_alt1_alt2_alt3[3].second != NULL);
    
    //int a0LODQA = ref_alt1_alt2_alt3[0].second->ALODQ;
    // int a0LODQB = ref_alt1_alt2_alt3[0].second->BLODQ;
    // int a0LODQ = MIN(a0LODQA, a0LODQB) + 1;
    int a0LODQ = ALODQ(ref_alt1_alt2_alt3[0].second);
    int a1LODQ = ALODQ(ref_alt1_alt2_alt3[1].second);
    int a2LODQ = ALODQ(ref_alt1_alt2_alt3[2].second);
    int a3LODQ = ALODQ(ref_alt1_alt2_alt3[3].second);
    
    /*
    unsigned int readlen = MAX(30U, central_readlen);
    const double alt1frac = (double)(readlen 
            // - MIN(readlen - 30, ref_alt1_alt2_alt3[1].second->RefBias / 2.0 )
            ) / (double)readlen / 2.0;
    const double alt2frac = (double)(readlen 
            // - MIN(readlen - 30, ref_alt1_alt2_alt3[2].second->RefBias / 2.0)
            ) / (double)readlen / 2.0;
    */
    auto fmtptr0 = ref_alt1_alt2_alt3[0].second;
    auto fmtptr1 = ref_alt1_alt2_alt3[1].second;
    auto fmtptr2 = ref_alt1_alt2_alt3[2].second;
    const bool isSubst = isSymbolSubstitution(refsymbol);
    const AlignmentSymbol symbolNN = BASE_NN; // (isSubst ? BASE_NN : LINK_NN);
    double ad0norm = compute_norm_ad(fmtptr0, isSubst);
    double ad1norm = compute_norm_ad(fmtptr1, isSubst, symbolNN == ref_alt1_alt2_alt3[1].first);
    double ad2norm = compute_norm_ad(fmtptr2, isSubst, symbolNN == ref_alt1_alt2_alt3[2].first);
    if (symbolNN == ref_alt1_alt2_alt3[1].first) {
        ad0norm += ad1norm;
        ad1norm = 0;
    }
    if (symbolNN == ref_alt1_alt2_alt3[2].first) {
        ad0norm += ad2norm;
        ad2norm = 0;
    }
    int a0a1LODQ = hetLODQ(ad0norm, ad1norm); //, 1.0 - alt1frac);
    int a1a0LODQ = hetLODQ(ad1norm, ad0norm); //, alt1frac);
    int a1a2LODQ = hetLODQ(ad1norm, ad2norm); //, alt1frac / (alt1frac + alt2frac));
    int a2a1LODQ = hetLODQ(ad2norm, ad1norm); //, alt2frac / (alt1frac + alt2frac));
    
    std::array<std::string, 4> GTidx2GT {{
        "0/0",
        "0/1",
        "1/1",
        "1/2"
    }};
    
    int phred_homref = 0; // (isSubst ? 31 : 41);
    int phred_hetero = (isSubst ? 31 : 41-1);
    int phred_homalt = (isSubst ? 33 : 43-1);
    int phred_tri_al = (isSubst ? 60 : 41-1); // 49 // https://www.genetics.org/content/184/1/233 : triallelic-SNP-phred = 29*2-3
    // tri_al for InDels is lower than expected because indels were already penalized for tri-allelelity (aka triallelic InDels) in their TLODQs
    // const double qfrac = (isSubst ? 1.0 : 0.25);
    std::array<std::pair<int, int>, 4> GL4raw = {{
        std::make_pair(0,     (-phred_homref - a1LODQ              - MAX(a2LODQ - phred_hetero - 3, 0))),
        std::make_pair(1,  MIN(-phred_hetero - MAX(a0a1LODQ, a1a0LODQ), -a2LODQ)),
        std::make_pair(2,     (-phred_homalt - MAX(a0LODQ, a2LODQ) - MAX(MIN(a0LODQ, a2LODQ) - phred_hetero - 3, 0))),
        std::make_pair(3,  MIN(-phred_tri_al - MAX(a1a2LODQ, a2a1LODQ), -MAX(a0LODQ, a3LODQ)))
    }};
    auto ret = GL4raw[0].second - MAX3(GL4raw[1].second, GL4raw[2].second, GL4raw[3].second);
    if ((0x3 & outputflag) == 0) {
        return ret;
    }
    
    auto GL4 = GL4raw;
    std::sort(GL4.rbegin(), GL4.rend(), PairSecondLess);
    
    size_t GLidx = GL4[0].first;
    
    std::array<std::string, 3> ref_alt1_alt2_vcfstr_arr = {{
        SYMBOL_TO_DESC_ARR[ref_alt1_alt2_alt3[0].first],
        SYMBOL_TO_DESC_ARR[ref_alt1_alt2_alt3[1].first],
        SYMBOL_TO_DESC_ARR[ref_alt1_alt2_alt3[2].first]
    }};
    
    std::string vcfref = "";
    std::string vcfalt = "";
    unsigned int alt1_uniallelic_phred = 200;
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
        std::vector<std::pair<std::array<unsigned int, 2>, std::string>> indelstrings1 = indel_get_majority(*(ref_alt1_alt2_alt3[1].second), 
            // false, tki, false, 
            tname, refpos, s1, false);
        const std::string & indelstring1 = indelstrings1[0].second;
        if (indelstrings1.size() > 1) {
            alt1_uniallelic_phred = (phred_tri_al - phred_hetero) *  log(1.0 + (double)indelstrings1[0].first[1] / (double)indelstrings1[1].first[1]) / log(2.0);
            if (a1LODQ > MAX(a0LODQ, a2LODQ) + alt1_uniallelic_phred) {
                // is_alt1_uniallelic = false;
            }
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
                std::vector<std::pair<std::array<unsigned int, 2>, std::string>> indelstrings2 = indel_get_majority(*(ref_alt1_alt2_alt3[2].second), 
                    // false, tki,false, 
                    tname, refpos, ref_alt1_alt2_alt3[2].first, false);
                indelstring2 = indelstrings2[0].second;
                if (s2 == s1) {
                    if (indelstrings2.size() < 2) {
                        fprintf(stderr, "Runtime error: indel at %u is invalid!\b", refpos);
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
                    auto minsize = MIN(indelstring1.size(), indelstring2.size());
                    assert(indelstring1.substr(0, minsize) == indelstring2.substr(0, minsize));
                    if (indelstring1.size() > indelstring2.size()) {
                        vcfref = vcfref1 + indelstring1;
                        vcfalt = vcfalt1 + "," + vcfalt1 + indelstring1.substr(indelstring2.size());
                    } else {
                        vcfref = vcfref1 + indelstring2;
                        vcfalt = vcfalt1 + indelstring2.substr(indelstring1.size()) + "," + vcfalt1;
                    }
                    // vcfref = vcfref1 + indelstring1 + "," + vcfref1 + indelstring2;
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
    int germ_GQ =  (is_alt1_uniallelic ? (GL4[0].second - GL4[1].second) : MIN(alt1_uniallelic_phred, GL4[0].second - GL4[1].second));
    std::vector<int> germ_ADR;
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[0].second->cDP1a, 1, 0));
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[1].second->cDP1a, 1, 0));
    if (3 == GLidx) {
        germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[2].second->cDP1a, 1, 0));
    }
    std::vector<std::string> germ_FT;
    //germ_FT.push_back(ref_alt1_alt2_alt3[0].second->FT);
    //germ_FT.push_back(ref_alt1_alt2_alt3[1].second->FT);
    if (3 == GLidx) {
        // germ_FT.push_back(ref_alt1_alt2_alt3[2].second->FT);
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
        std::string("GT:GQ:HQ:FT:DP:cADR:GLa:GSTa:note"),
        string_join(std::array<std::string, 9>
        {{
            germ_GT, 
            std::to_string(germ_GQ), 
            std::string("0,0"), 
            string_join(germ_FT, "/"),
            std::to_string(symbol_format_vec[0].second->DP),
            other_join(germ_ADR, std::string(",")), 
            other_join(std::array<int, 3>
            {{ 
                GL4raw[0].second, GL4raw[1].second, GL4raw[2].second 
            }}, ","),
            other_join(std::array<int, 10>
            {{
                GL4raw[3].second,
                // a0LODQA, a0LODQB, 
                a0LODQ, a1LODQ, a2LODQ, a3LODQ,
                a0a1LODQ, a1a0LODQ, a1a2LODQ, a2a1LODQ, 0
            }}, ","),
            ref_alt1_alt2_alt3[0].second->note
        }}, ":")
    }}, "\t") + "\n";
    out_string += bcfline;
    return ret;
}

/*
int
fill_TN_germline(
        bcfrec::BcfFormat & fmt, 
        AlignmentSymbol symbol,
        double powlaw_exponent,
        double refmul,
        double altmul,
        double any_mul_contam_frac,
        double ref_mul_contam_frac, 
        double t2n_add_contam_transfrac,
        bool prev_is_tumor, 
        const auto & tki, 
        bool is_novar,
        bool somaticGT,
        const unsigned int specialflag) {
    
    const bool isInDel = (isSymbolIns(symbol) || isSymbolDel(symbol));
    // likelihood that the reads are generated by tumor contam, other contam, hetero genotype, homo-alt genotype, (and homo-ref genotype)
    std::array<unsigned int, 2> pranks = {{0, 0}}; // number of alleles for REF vs ALT+other and ALT vs REF+other.
    for (unsigned int t = 0; t < 2; t++) {
        auto & fmtGT = (0 == t ? fmt.GTa : fmt.GTb);
        auto & fmtGQ = (0 == t ? fmt.GQa : fmt.GQb);
        auto & fmtGL = (0 == t ? fmt.GLa : fmt.GLb);
        auto & fmtGST = (0 == t ? fmt.GSTa : fmt.GSTb);
        auto & prank = pranks[t];
        
        auto fmt_AD = MAX(0.0, fmt.FA) * fmt.DP;
        auto fmt_OD = MAX(0.0, 1.0 - (fmt.FA + (is_novar ? 0 : fmt.FR))) * fmt.DP;
        auto fmt_RD = MAX(0.0, fmt.FR) * fmt.DP;
        auto pre1AD = ((0 == t) ? fmt_AD : fmt_OD);
        auto pre2AD = ((0 == t) ? MAX(fmt_RD, fmt_OD) : MAX(fmt_RD, fmt_AD));
        
        auto fmt_bAD = SUM2(fmt.bAD1);
        assert(SUM2(fmt.bDP1) >= (SUM2(fmt.bAD1) + (is_novar ? 0 : SUM2(fmt.bRD1))));
        auto fmt_bOD = SUM2(fmt.bDP1) - (SUM2(fmt.bAD1) + (is_novar ? 0 : SUM2(fmt.bRD1)));
        auto fmt_bRD = SUM2(fmt.bRD1);
        auto pre1bAD = ((0 == t) ? fmt_bAD : fmt_bOD);
        auto pre2bAD = ((0 == t) ? MAX(fmt_bRD, fmt_bOD) : MAX(fmt_bRD, fmt_bAD));
        
        auto fmt_bABQ = SUM2(fmt.bAltBQ);
        assert(SUM2(fmt.bAllBQ) >= (SUM2(fmt.bAltBQ) + (is_novar ? 0 : SUM2(fmt.bRefBQ))));
        auto fmt_bOBQ = SUM2(fmt.bAllBQ) - (SUM2(fmt.bAltBQ) + (is_novar ? 0 : SUM2(fmt.bRefBQ)));
        auto fmt_bRBQ = SUM2(fmt.bRefBQ);
        int pre1bBQ = ((0 == t) ? fmt_bABQ : fmt_bOBQ);
        int pre2bBQ = ((0 == t) ? MAX(fmt_bRBQ, fmt_bOBQ) : MAX(fmt_bRBQ, fmt_bABQ));
        
        int allAvgBQ = SUM2(fmt.bAllBQ) / (SUM2(fmt.bDP1) + 1);
        assert(1 + pre1bAD != 0);
        assert(1 + pre2bAD != 0);
        int pre1bBQAvg = (int)((pre1bBQ + allAvgBQ) / (1 + pre1bAD));
        int pre2bBQAvg = (int)((pre2bBQ + allAvgBQ) / (1 + pre2bAD));
        if (isInDel) {
            pre1bBQAvg = 0;
            pre2bBQAvg = 0;
        }

        // const double preFrac = MAX(0, (0 == t) ? (fmt.FA + MAX(fmt.FR, fmt_FO)) : (1.0 - MIN(fmt.FA, fmt.FR)));
        // const double fa1 = MAX(0.0, ((0 == t) ? (fmtFA / (fmtFA + fmtFR)) : ((1.0 - fmtFA - fmtFR) / ()) ));
        const double fa1 = (double)pre1AD / (double)(pre1AD + pre2AD); // / preFrac;
        
        const double fa_l = (pre1AD + 0.5) / (pre1AD + pre2AD + 1.0);
        //const double da_l = fa_l * fmt.DP;
        const double fr_l = 1.0 - fa_l;
        //const double dr_l = fr_l * fmt.DP;
        assert (fa_l > 0 && fa_l < 1 || !fprintf(stderr, "The fraction %lf is not a valid fa_l\n", fa_l));
        assert (fr_l > 0 && fr_l < 1 || !fprintf(stderr, "The fraction %lf is not a valid fr_l\n", fr_l));
        
        const double fa_v = (pre1AD +  DBL_EPSILON) / (pre1AD + pre2AD + 2 * DBL_EPSILON);
        const double da_v = fa_v * (pre1AD + pre2AD);
        const double fr_v = 1.0 - fa_v;
        const double dr_v = fr_v * (pre1AD + pre2AD);
        
        double t2n_conref_frac = 0.0;
        double t2n_conalt_frac = 0.0;
        int homref_likecon1 = -200;
        int homalt_likecon1 = -200;
        
        if (prev_is_tumor) {
            double tkiFA = ((!isSymbolSubstitution(symbol) || true) ? tki.FA : (
                         (double)(tki.cAltBQ) / (         (double)(tki.cAllBQ) + DBL_MIN)));
            double tkiFR = ((!isSymbolSubstitution(symbol) || true) ? tki.FR : (
                         (double)(tki.cRefBQ) / (         (double)(tki.cAllBQ) + DBL_MIN)));
            const double tki_fa1 = MAX(0.0, ((0 == t) ? (tkiFA) : (1.0 - tkiFA - tkiFR)));
            const double tki_fa_l = (tki_fa1 * (double)tki.DP + 0.5) / (double)(tki.DP + 1.0);
            const double tki_fr_l = 1.0 - tki_fa_l;
            homref_likecon1 = -(int)calc_binom_10log10_likeratio(t2n_add_contam_transfrac, fa_l * fmt.DP, tki_fa_l * tki.DP);
            homalt_likecon1 = -(int)calc_binom_10log10_likeratio(t2n_add_contam_transfrac, fr_l * fmt.DP, tki_fr_l * tki.DP);
            t2n_conref_frac += (tki_fr_l * any_mul_contam_frac + tki_fr_l * any_mul_contam_frac * MIN(5.0, tki.DP / (double)(fmt.DP + DBL_MIN)));
            t2n_conalt_frac += (tki_fa_l * any_mul_contam_frac + tki_fa_l * any_mul_contam_frac * MIN(5.0, tki.DP / (double)(fmt.DP + DBL_MIN)));
        }
        
        // two models (additive and multiplicative)
        // two alleles (REF and ALT)
        // two sources of stochasticity (heterozygosity and contamination)
        // two genotyes (REF-with-ALT genotype and REF-with-(all-minus-REF-minus-ALT) genotype)
        // = 16 combinations in total
        // ContEst https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3167057/ Fig. 1 TCGA Ovarian 
                    
        // assuming classical preferential attachment, deviation from its theoretical distribution is translated into a phred-scaled error probability.
        // zero minus the phred-scale probability that ALT is generated by a stochastic process of power-law error given that ref is good, higher means ref is of better quality (zero is best)
        // logit2(theoretical-deviation, observed-deviation), higher-than-expected deviation <=> more-negative-score
        int hetREF_likelim1 =  (int)(10.0/log(10.00) * powlaw_exponent * MIN(logit2(fr_l / refmul * (1.0 + 0.0),                                    fa_l / altmul), 0.0));
        int hetALT_likelim1 =  (int)(10.0/log(10.00) * powlaw_exponent * MIN(logit2(fa_l / altmul * (1.0 + 0.0),                                    fr_l / refmul), 0.0));         
        int homref_likelim1 =  (int)(10.0/log(10.00) * powlaw_exponent * MIN(logit2(fr_l / refmul * any_mul_contam_frac + t2n_conalt_frac / altmul, fa_l / altmul), 0.0));
        int homalt_likelim1 =  (int)(10.0/log(10.00) * powlaw_exponent * MIN(logit2(fa_l / altmul * any_mul_contam_frac + t2n_conref_frac / refmul, fr_l / refmul), 0.0));
        
        // assuming statistical independence of reads, kl-divergence is translated into a phred-scaled error probability.
        // binom_10log10_likeratio(theoretical-deviation-rate, number-of-deviation-signals, number-of-all-signals), higher-than-expected deviation <=> more-negative-score
        int hetREF_likeval1 = -(int)calc_binom_10log10_likeratio(0.500 * altmul, da_v, dr_v);                                 // het-ref to ALT add error phred
        int hetALT_likeval1 = -(int)calc_binom_10log10_likeratio(0.500 * refmul, dr_v, da_v);                                 // het-alt to REF add error phred
        int homref_likeval1 = -(int)calc_binom_10log10_likeratio(any_mul_contam_frac * altmul + t2n_conalt_frac, da_v, dr_v); // hom-alt to REF add error phred by contamination
        double any_mul_contam_frac2 = any_mul_contam_frac * refmul + t2n_conref_frac;
        int homalt_likeval1 = -(int)calc_binom_10log10_likeratio(MAX(any_mul_contam_frac2, ref_mul_contam_frac), dr_v, da_v); // hom-ref to ALT add error phred by contamination
        
        int homref_likelim2 = MAX(homref_likelim1, homref_likecon1);
        int homref_likeval2 = MAX(homref_likeval1, homref_likecon1);
        
        int homalt_likelim2 = MAX(homalt_likelim1, homalt_likecon1);
        int homalt_likeval2 = MAX(homalt_likeval1, homalt_likecon1);
        
        fmtGST = {
                -homref_likelim1, -homref_likeval1, -homref_likecon1,
                -hetREF_likelim1, -hetALT_likeval1,
                -hetALT_likelim1, -hetREF_likeval1,
                -homalt_likelim1, -homalt_likeval1, -homalt_likecon1
        };
        
        // An important note: two models must generate different alleles to perform model selection
        fmtGL[0] =     MAX(homref_likelim2, homref_likeval2)                                         + pre2bBQAvg;
        fmtGL[1] = MIN(MAX(hetREF_likelim1, hetREF_likeval1), MAX(hetALT_likelim1, hetALT_likeval1)) + MIN(pre1bBQAvg, pre2bBQAvg);
        fmtGL[2] =     MAX(homalt_likelim2, homalt_likeval2)                                         + pre1bBQAvg;
        
        std::array<int, 3> likes = {{ fmtGL[0], fmtGL[1], fmtGL[2] }};
        std::sort(likes.rbegin(), likes.rend());
        const auto & gt_homref = (somaticGT ? TT_HOMREF : GT_HOMREF);
        const auto & gt_homalt = (somaticGT ? TT_HOMALT : GT_HOMALT);
        const auto & gt_hetero = (somaticGT ? TT_HETERO : GT_HETERO);
        if        (likes[0] == fmtGL[0]) {
            fmtGT = (is_novar ? gt_homalt[t] : gt_homref[t]);
            prank = (is_novar ? 2 : 0);
        } else if (likes[0] == fmtGL[1]) {
            fmtGT = (is_novar ? gt_hetero[t] : gt_hetero[t]);
            prank = (is_novar ? 1 : 1);
        } else if (likes[0] == fmtGL[2]) {
            fmtGT = (is_novar ? gt_homref[t] : gt_homalt[t]);
            prank = (is_novar ? 0 : 2);
        } else {
            abort(); // should not happen
        }
        fmtGQ = likes[0] - likes[1];
    }
    // homalt (highest priority) > hetero > homref (lowest priority), tiebreak with GTa > GTb
    if (pranks[0] > pranks[1] || (pranks[0] == pranks[1] && fmt.GQa >= fmt.GQb)) {
        fmt.GT = fmt.GTa;
        fmt.GQ = fmt.GQa;
    } else {
        fmt.GT = fmt.GTb;
        fmt.GQ = fmt.GQb;
    }
}
*/

#include "version.h"
std::string 
generate_vcf_header(const char *ref_fasta_fname, 
        const char *platform, 
        unsigned int central_readlen,
        //const unsigned int minABQ_pcr_snv, 
        //const unsigned int minABQ_pcr_indel, 
        //const unsigned int minABQ_cap_snv, 
        //const unsigned int minABQ_cap_indel, 
        unsigned int argc,
        const char *const *argv,
        unsigned int n_targets,
        const char *const *target_name,
        const uint32_t *target_len,
        const char *const sampleName,
        const char *const tumor_sampleName,
        bool is_tumor_format_retrieved) {
    time_t rawtime;
    time(&rawtime);
    char timestring[80];
    strftime(timestring, 80, "%F %T", localtime(&rawtime));
    std::string ret = "";
    ret += std::string("") + "##fileformat=VCFv4.2" + "\n" ;
    ret += std::string("") + "##fileDate=" + timestring + "\n" ;
    ret += std::string("") + "##variantCallerVersion=" + VERSION_DETAIL + "\n";
    ret += std::string("") + "##variantCallerCommand=";
    for (unsigned int i = 0; i < argc; i++) {
        ret += std::string("") + std::string(argv[i]) + "  ";
    }
    ret += "\n";
    ret += std::string("") + "##variantCallerInferredParameters=<" 
            + "inferred_sequencing_platform=" + platform 
            + ",central_readlen=" + std::to_string(central_readlen) 
            // + ",minABQs=(" + std::to_string(minABQ_pcr_snv) + "x" + std::to_string(minABQ_pcr_indel) + "x" +  std::to_string(minABQ_cap_snv) + "x" + std::to_string(minABQ_cap_indel) 
            // + ")" 
            + ">\n";
    ret += std::string("") + "##reference=" + ref_fasta_fname + "\n";
    for (size_t i = 0; i < n_targets; i++) {
        ret += std::string("") + "##contig=<ID=" + target_name[i] + ",length=" + std::to_string(target_len[i]) + ">\n";
    }
    ret += std::string("") + "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location, where POS (start position) is one-based inclusive.\">\n";
    
    for (unsigned int i = 0; i < bcfrec::FILTER_NUM; i++) {
        ret += std::string("") + bcfrec::FILTER_LINES[i] + "\n";
    }
    
    ret += "##INFO=<ID=ANY_VAR,Number=0,Type=Flag,Description=\"Any type of variant which may be caused by germline polymorphism and/or somatic mutation\">\n";
    ret += "##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"germline variant\">\n";
    ret += "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic variant\">\n";
    ret += "##INFO=<ID=SomaticQ,Number=A,Type=Float,Description=\"Somatic quality of the variant, the PHRED-scale probability that this variant is not somatic.\">\n";
    ret += "##INFO=<ID=TLODQ,Number=A,Type=Float,Description=\"Tumor log-of-data-likelihood quality, the PHRED-scale probability that this variant is not of biological origin (i.e., artifactual).\">\n";
    ret += "##INFO=<ID=NLODQ,Number=A,Type=Float,Description=\"Normal log-of-data-likelihood quality, the PHRED-scale probability that this variant is of germline origin.\">\n";
    ret += "##INFO=<ID=RU,Number=1,Type=String,Description=\"The shortest repeating unit in the reference\">\n";
    ret += "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"The number of non-interrupted RUs in the reference\">\n";
    ret += "##INFO=<ID=R3X2,Number=6,Type=Integer,Description=\"Repeat start position, repeat track length, and repeat unit size at the two positions before and after this VCF position.\">\n"; 
    
    for (unsigned int i = 0; i < bcfrec::FORMAT_NUM; i++) {
        ret += std::string("") + bcfrec::FORMAT_LINES[i] + "\n";
    }
    ret += std::string("") + "##FORMAT=<ID=gbDP,Number=1,Type=Integer,Description=\"Minimum duped   fragment depths in the genomic block for SNV and InDel\">\n";
    ret += std::string("") + "##FORMAT=<ID=gcDP,Number=1,Type=Integer,Description=\"Minimum deduped fragment depths in the genomic block for SNV and InDel\">\n";
    ret += std::string("") + "##FORMAT=<ID=gSTS,Number=2,Type=Integer,Description=\"Variant types for start and end positions, where 0 means SNV and 1 means InDel.\">\n";
    ret += std::string("") + "##FORMAT=<ID=gBEG,Number=1,Type=Integer,Description=\"Begin position of the genomic block (one-based inclusive)\">\n";
    ret += std::string("") + "##FORMAT=<ID=gEND,Number=1,Type=Integer,Description=\"End position of the genomic block (one-based inclusive)\">\n";
    ret += std::string("") + "##phasing=partial\n";
    ret += std::string("") + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" 
            + sampleName + ((tumor_sampleName != NULL && is_tumor_format_retrieved) ? (std::string("\t") + tumor_sampleName) : std::string("")) + "\n";
    return ret;
}

/*
auto
fmtFTSupdate(auto & maxval, std::string & ft, std::vector<unsigned int> & ftv, const char *fkey , const auto fthres, const auto fval) {
    maxval = MAX(maxval, fval);
    if ((unsigned int)fthres <= (unsigned int)fval) {
        ft  += (std::string(fkey) + "&"); // amperstand-separated list of filter strings
        ftv.push_back((unsigned int)fval);
    }
    return fval;
}
*/

std::string
bcf1_to_string(const bcf_hdr_t *tki_bcf1_hdr, const bcf1_t *bcf1_record) {
    kstring_t ks = { 0, 0, NULL };
    vcf_format(tki_bcf1_hdr, bcf1_record, &ks);
    assert (ks.l > 2); 
    size_t idx = ks.l - 1;
    // fprintf(stderr, "kstring_t.size = %d\n", idx);
    for (;idx != 0 && ks.s[idx] != '\t'; idx--) {
    }
    std::string ret = std::string(ks.s, idx, ks.l - idx - 1);
    if (ks.s != NULL) {
        free(ks.s);
    }
    // fprintf(stderr, "tumor-vcf-format=%s\n", ret.c_str());
    return ret;
}

int 
fill_tki(auto & tki, const auto & fmt, unsigned int a = 1) {
    tki.BDP = fmt.BDPf.at(0) +fmt.BDPr.at(0); // 
    tki.bDP = fmt.bDPf.at(a) +fmt.bDPr.at(a); // 
    tki.CDP1v = fmt.CDP1v.at(0);
    tki.cDP1v = fmt.cDP1v.at(a);
    tki.cVQ1  = fmt.cVQ1.at(a);
    tki.CDP2v = fmt.CDP2v.at(0);
    tki.cDP2v = fmt.cDP2v.at(a);
    tki.cVQ2  = fmt.cVQ2.at(a);
};

const auto
calc_binom_powlaw_syserr_normv_quals(
        auto tAD, auto tDP, auto tVQ,
        auto nAD, auto nDP, auto nVQ, bool is_penal_applied) {
    int binom_b10log10like = (int)calc_binom_10log10_likeratio((tAD + 0.5) / (tDP + 1.0), nAD, nDP);
    double bjpfrac = ((tAD + 0.5) / (tDP + 1.0)) / ((nAD + 0.5) / (nDP + 1.0));
    int powlaw_b10log10like = (int)(3 * 10 / log(10) * log(bjpfrac));
    int syserr_b10log10like = (int)MAX(0, nVQ - 12.5 * mathsquare(MAX(0, bjpfrac - 1.0)));
    int tnVQ = tVQ + MIN(22, CENTER(binom_b10log10like, powlaw_b10log10like));
    if (is_penal_applied) { tnVQ -= syserr_b10log10like; }
    return std::array<int, 4> {{binom_b10log10like, powlaw_b10log10like, syserr_b10log10like, tnVQ }};
};

int
append_vcf_record(
        std::string & out_string,
        const char *tname,
        const unsigned int refpos,
        const unsigned int region_offset,
        const std::string & refstring,
        const std::vector<RegionalTandemRepeat> & region_repeatvec,
        const std::string & repeatunit,
        const unsigned int repeatnum,
        const AlignmentSymbol refsymbol,
        const AlignmentSymbol symbol,
        // const std::string & indelstring,
        const bcfrec::BcfFormat & fmt,
        auto & tki,
        const bool is_tumor_format_retrieved,
        const int nlodq,
        const double vcfqual_thres,
        const unsigned int all_allele_depth_thres,
        const unsigned int var_allele_depth_thres,
        const bool should_output_ref_allele,
        const bcf_hdr_t *g_bcf_hdr,
        unsigned int specialflag) {
    
    const std::string & indelstring = LAST(fmt.gapSa);
    const bcfrec::BcfFormat & nfm = (is_tumor_format_retrieved ? fmt : FORMAT_UNCOV);
    if (!is_tumor_format_retrieved) {
        fill_tki(tki, fmt);
    }
    const auto regionpos = refpos - region_offset;
    unsigned int vcfpos;
    std::string vcfref, vcfalt;
    if (indelstring.size() > 0) {
        vcfpos = refpos; // refpos > 0?
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
    
    const auto b_binom_powlaw_syserr_normv_q4 = calc_binom_powlaw_syserr_normv_quals(tki.cDP1v, tki.CDP1v, tki.cVQ1,
            nfm.cDP1v[1], nfm.CDP1v[1], nfm.cVQ1[1], (indelstring.size() > 0));
    const auto c_binom_powlaw_syserr_normv_q4 = calc_binom_powlaw_syserr_normv_quals(tki.cDP2v, tki.CDP2v, tki.cVQ2,
            nfm.cDP2v[1], nfm.CDP2v[1], nfm.cVQ2[1], true);
    
    int tlodq = MAX(b_binom_powlaw_syserr_normv_q4[3], c_binom_powlaw_syserr_normv_q4[3]);
    int somaticq = MIN(tlodq, nlodq);
    int vcfqual = (is_tumor_format_retrieved ? somaticq : tlodq);
    std::string infostring = std::string(is_tumor_format_retrieved ? "SOMATIC" : "ANY_VAR");
    infostring += std::string(";SomaticQ=") + std::to_string(somaticq);
    infostring += std::string(";TLODQ=") + std::to_string(tlodq);
    infostring += std::string(";NLODQ=") + std::to_string(nlodq);
    infostring += std::string(";TNBQ4=") + other_join(b_binom_powlaw_syserr_normv_q4);
    infostring += std::string(";TNCQ4=") + other_join(c_binom_powlaw_syserr_normv_q4);
    infostring += std::string(";RU=") + repeatunit + ";RC=" + std::to_string(repeatnum);
    const auto & rtr1 =  region_repeatvec.at(MAX(refpos - region_offset, 1) - 1);
    const auto & rtr2 =  region_repeatvec.at(MIN(refpos - region_offset + 1, region_repeatvec.size() - 1));
    const auto rtr1_tpos = ((0 == rtr1.tracklen) ? 0 : (region_offset + rtr1.begpos));
    const auto rtr2_tpos = ((0 == rtr2.tracklen) ? 0 : (region_offset + rtr2.begpos));
    infostring += std::string(";R3X2=") + other_join(std::array<unsigned int, 6>{{rtr1_tpos, rtr1.tracklen, rtr1.unitlen, rtr2_tpos, rtr2.tracklen, rtr2.unitlen}});
    
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
    
    const bool keep_var = (((vcfqual >= vcfqual_thres) || (tki.BDP >= all_allele_depth_thres) || (tki.bDP >= var_allele_depth_thres)) && (symbol != refsymbol || should_output_ref_allele));
    if (keep_var) {
        out_string += string_join(std::array<std::string, 9>{{
                std::string(tname), std::to_string(vcfpos), ".", vcfref, vcfalt, std::to_string(vcfqual), vcffilter, 
                infostring, bcfrec::FORMAT_STR_PER_REC}});
        bcfrec::streamAppendBcfFormat(out_string, fmt);
        out_string += ((g_bcf_hdr != NULL && is_tumor_format_retrieved) ? bcf1_to_string(g_bcf_hdr, tki.bcf1_record) : std::string("")) + "\n";
    }
    return 0;
};

#endif
