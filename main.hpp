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
    const T & src_indel = (indel2count.size() > 1 ? T() : indel2count.begin()->first);
    assert (indel2count.begin()->second > 0);
    posToIndelToCount_inc<T>(dst, epos, src_indel, incvalue);
    return src_indel;
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

#define INS_N_ANCHOR_BASES 1 //  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC149199/

#define bam_phredi(b, i) (bam_get_qual((b))[(i)])

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

// Warning: misuse can cause arr-out-of-bound error
const AlignmentSymbol SYMBOL_TYPE_TO_NO_VAR_SYMBOL[NUM_SYMBOL_TYPES] = {
    [BASE_SYMBOL] = BASE_NO_SNV,
    [LINK_SYMBOL] = LINK_NO_INDEL,
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
                return indel_open;
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

unsigned int is_mut_transition(const AlignmentSymbol con_symbol, const AlignmentSymbol alt_symbol) {
    return ((con_symbol == BASE_C && alt_symbol == BASE_T) || (con_symbol == BASE_G && alt_symbol == BASE_A)
        || (con_symbol == BASE_T && alt_symbol == BASE_C) || (con_symbol == BASE_A && alt_symbol == BASE_G)
    );
}

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
    [BASE_NO_SNV] = "<NO_SNV>",
    [LINK_NO_INDEL] = "<NO_INDEL>"
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
    std::array<T, NUM_ALIGNMENT_SYMBOLS> symbol2data;
};

typedef uint32_t molcount_t;

template <class TB2C>
class GenericSymbol2Bucket2Count : TDistribution<TB2C> {
public:    
    molcount_t
    getSymbolBucketCount(AlignmentSymbol symbol, unsigned int bucket) const {
        assert(bucket < NUM_BUCKETS);
        return this->symbol2data[symbol][bucket];
    };
    
    const TB2C &
    getSymbolCounts(AlignmentSymbol symbol) const {
        return this->symbol2data[symbol];
    };

    void
    incSymbolBucketCount(AlignmentSymbol symbol, unsigned int bucket, unsigned int increment) {
        assert(bucket < NUM_BUCKETS || !fprintf(stderr, "%d < %d failed!", bucket, NUM_BUCKETS));
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
    clearSymbolBucketCount() {
        assert(sizeof(TB2C) * NUM_ALIGNMENT_SYMBOLS == sizeof(this->symbol2data) || !fprintf(stderr, "%lu * %u != %lu\n", sizeof(TB2C), NUM_ALIGNMENT_SYMBOLS, sizeof(this->symbol2data)));
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
    
    template <bool TIndelIsMajor>
    int
    _fillConsensusCounts(
            AlignmentSymbol & count_argmax, unsigned int & count_max, unsigned int & count_sum,
            AlignmentSymbol incluBeg, AlignmentSymbol incluEnd) const {
        assert (incluBeg <= incluEnd);
        
        count_argmax = incluEnd; // assign with the value AlignmentSymbol(NUM_ALIGNMENT_SYMBOLS) to flag for error
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
    
    template<bool TIndelIsMajor>
    const AlignmentSymbol
    _updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count,
            const SymbolType symbolType, const AlignmentSymbol ambig_pos, unsigned int incvalue) {
        AlignmentSymbol argmax_count = END_ALIGNMENT_SYMBOLS; // AlignmentSymbol(0) is not fully correct
        unsigned int max_count = 0;
        unsigned int sum_count = 0;
        thatSymbol2Count.template fillConsensusCounts<TIndelIsMajor>(argmax_count, max_count, sum_count, symbolType);
        
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
    updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count, unsigned int incvalue = 1) {
        AlignmentSymbol baseSymb = this->template _updateByConsensus<false        >(thatSymbol2Count, BASE_SYMBOL, BASE_NN, incvalue);
        AlignmentSymbol linkSymb = this->template _updateByConsensus<TIndelIsMajor>(thatSymbol2Count, LINK_SYMBOL, LINK_NN, incvalue);
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
    
    int
    updateByFiltering(
            std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> & con_symbols, 
            const GenericSymbol2Count<TInteger> & other, 
            const std::array<unsigned int, NUM_SYMBOL_TYPES> thres = {{0}},
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
            if (adjcount >= thres[symbolType] && adjcount > 0) {
                this->symbol2data[consalpha] += incvalue;
                ret++;
            }
            con_symbols[symbolType] = consalpha;
        }
        return ret;
    };
    
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
    const uint32_t incluBegPosition; // where end_pos = incluBegPosition + idx2symbol2data.size()

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
    GenericSymbol2Bucket2CountCoverage() : CoveredRegion<TSymbol2Bucket2Count>(0, 0, 1) { }
    GenericSymbol2Bucket2CountCoverage(auto tid, auto beg, auto end) : CoveredRegion<TSymbol2Bucket2Count>(tid, beg, end) {}
};

typedef std::array<molcount_t, NUM_BUCKETS> Bucket2Count;

typedef GenericSymbol2Bucket2Count<Bucket2Count> Symbol2Bucket2Count;

typedef GenericSymbol2Count<uint32_t> Symbol2Count;
typedef GenericSymbol2Count<uint64_t> Symbol2CountUint64;

typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2Count> Symbol2Bucket2CountCoverage;

typedef CoveredRegion<std::array<molcount_t, NUM_SEG_FORMAT_PREP_SETS>> SegFormatPrepSets;
typedef CoveredRegion<std::array<molcount_t, NUM_SEG_FORMAT_THRES_SETS>> SegFormatThresSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_SEG_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2SegFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_FRAG_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2FragFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_FAM_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2FamFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_DUPLEX_FORMAT_DEPTH_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2DuplexFormatDepthSets;
typedef CoveredRegion<std::array<std::array<molcount_t, NUM_VQ_FORMAT_TAG_SETS>, NUM_ALIGNMENT_SYMBOLS>> Symbol2VQFormatTagSets;

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
    if (rulen2 * rc2 == 0) {
        return true;
    }
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
    // derived from the python code: for i in range(1, 20): print("{}  {}".format(  int( round(10.0/log(10.0)*log(i)) ) , i  )) 
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
        return n_units_to_phred[MIN(n_units  , n_units_to_phred.size()-1)];
    } else {
        return n_units_to_phred[MIN(indel_len, n_units_to_phred.size()-1)];
    }
}

unsigned int
indel_phred(double ampfact, unsigned int cigar_oplen, unsigned int repeatsize_at_max_repeatnum, unsigned int max_repeatnum) {
    unsigned int region_size = repeatsize_at_max_repeatnum * max_repeatnum;
    double num_slips = (region_size > 64 ? (double)(region_size - 8) : log1p(exp((double)region_size - (double)8))) 
            * ampfact / ((double)(repeatsize_at_max_repeatnum * repeatsize_at_max_repeatnum));
    return prob2phred((1.0 - DBL_EPSILON) / (num_slips + 1.0));
    // e.g. for the short tadem repeat ACACAC : repeatsize_at_max_repeatnum = 2, indel_n_units = 3, repeatunit = AC
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
        const unsigned int decphred = indel_phred(8.0*(4.0+1.0), repeatsize_at_max_repeatnum, repeatsize_at_max_repeatnum, tl / repeatsize_at_max_repeatnum);
        for (unsigned int i = refpos; i != MIN(repeat_endpos, refstring.size()); i++) {
            if (tl > region_repeatvec[i].tracklen) {
                region_repeatvec[i].begpos = refpos;
                region_repeatvec[i].tracklen = tl;
                region_repeatvec[i].unitlen = repeatsize_at_max_repeatnum;
                region_repeatvec[i].indelphred = 30+18 - MIN(30+18-1, decphred);
            }
        }
        refpos += MAX(repeatsize_at_max_repeatnum * max_repeatnum, 6+1) - 6;
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
ref_to_phredvalue(unsigned int & n_units, 
        const auto & refstring, 
        const size_t refpos, 
        const unsigned int max_phred, 
        double ampfact, 
        const unsigned int cigar_oplen, 
        const auto cigar_op, 
        const bam1_t *b) {

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

int
update_seg_format_prep_sets_by_aln(
        SegFormatPrepSets & seg_format_prep_sets,
        const bam1_t *aln,
        const std::vector<RegionalTandemRepeat> & rtr_vec,
        const CoveredRegion<int64_t> & baq_offsetarr,
        const unsigned int region_offset,
        unsigned int specialflag) {
    
    const unsigned int rend = bam_endpos(aln);
    const auto cigar = bam_get_cigar(aln);
    const auto *bseq = bam_get_seq(aln);
    unsigned int nge_cnt = 0; // number of gap extensions.
    unsigned int ngo_cnt = 0; // number of gap opens.
    
    unsigned int insbaq_sum = 0;
    unsigned int delbaq_sum = 0;
    unsigned int inslen_sum = 0;
    unsigned int dellen_sum = 0;
    
    unsigned int qpos = 0;
    unsigned int rpos = aln->core.pos;
    for (unsigned int i = 0; i < aln->core.n_cigar; i++) {
        int32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        unsigned int cigar_oplen = bam_cigar_oplen(c);
        if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) { 
            nge_cnt += bam_cigar_oplen(c);
            ngo_cnt++;
            if (BAM_CINS == cigar_op) {
                insbaq_sum += baq_offsetarr.getByPos(MIN(rpos + cigar_oplen, baq_offsetarr.getExcluEndPosition() - 1)) - baq_offsetarr.getByPos(rpos);
                inslen_sum += bam_cigar_oplen(c);
                qpos += cigar_oplen;
            } else if (BAM_CDEL == cigar_op) {
                delbaq_sum += baq_offsetarr.getByPos(MIN(rpos + cigar_oplen, baq_offsetarr.getExcluEndPosition() - 1)) - baq_offsetarr.getByPos(rpos);
                dellen_sum += bam_cigar_oplen(c);
                rpos += cigar_oplen;
            }
        } else if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            qpos += cigar_oplen;
            rpos += cigar_oplen;
        } else {
            process_cigar(qpos, rpos, cigar_op, cigar_oplen);
        }
    }
    const uint8_t *bam_aux_data = bam_aux_get(aln, "NM");
    const unsigned int nm_cnt = ((bam_aux_data != NULL) ? bam_aux2i(bam_aux_data) : nge_cnt);
    assert (nm_cnt >= nge_cnt);
    const unsigned int xm_cnt = nm_cnt - nge_cnt;
    const unsigned int xm150 = xm_cnt * 1500 / (rend - aln->core.pos);
    const unsigned int go150 = ngo_cnt * 1500 / (rend - aln->core.pos);
    const unsigned int frag_pos_L = ((aln->core.isize != 0) ? MIN(aln->core.pos, aln->core.mpos) : aln->core.pos);
    const unsigned int frag_pos_R = ((aln->core.isize != 0) ? (frag_pos_L + aln->core.isize) : rend);
    const bool isrc = ((aln->core.flag & 0x10) == 0x10); 
    qpos = 0;
    rpos = aln->core.pos;
    
    for (unsigned int i = 0; i < aln->core.n_cigar; i++) {
        uint32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        unsigned int cigar_oplen = bam_cigar_oplen(c);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            for (unsigned int j = 0; j < cigar_oplen; j++) {
                seg_format_prep_sets.getRefByPos(rpos)[SEG_a_DP] += 1;
                seg_format_prep_sets.getRefByPos(rpos)[SEG_a_XM] += xm150;
                seg_format_prep_sets.getRefByPos(rpos)[SEG_a_GO] += go150;
                if (isrc) { 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_LI] += MIN(rpos - frag_pos_L + 1, MAX_INSERT_SIZE); 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_LIDP] += 1;
                } else { 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_RI] += MIN(frag_pos_R - rpos    , MAX_INSERT_SIZE); 
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_RIDP] += 1;
                }
                unsigned int base4bit = bam_seqi(bseq, qpos);
                unsigned int base3bit = seq_nt16_int[base4bit];
                if (bam_phredi(aln, qpos) >= 20) {
                    unsigned int ldist = rpos - aln->core.pos + 1;
                    unsigned int rdist = rend - rpos;
                    
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_L_DIST_SUM] += ldist;
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_R_DIST_SUM] += rdist;
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_INSLEN_SUM] += inslen_sum;
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_DELLEN_SUM] += dellen_sum;
                    
                    unsigned int lbaq = baq_offsetarr.getByPos(rpos) - baq_offsetarr.getByPos(aln->core.pos) + 1;
                    unsigned int rbaq = baq_offsetarr.getByPos(rend-1) - baq_offsetarr.getByPos(rpos) + 1;
                    
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_L_BAQ_SUM] += lbaq;
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_R_BAQ_SUM] += rbaq;
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_INSBAQ_SUM] += insbaq_sum;
                    seg_format_prep_sets.getRefByPos(rpos)[SEG_a_DELBAQ_SUM] += delbaq_sum;

                    seg_format_prep_sets.getRefByPos(rpos)[SEG_aBQ2_DP] += 1;
                }
                qpos++;
                rpos++;
            }
        } else if (cigar_op == BAM_CINS) {
            for (int rpos2 = MAX((int)rpos - (int)cigar_oplen - 1 - 1, aln->core.pos); rpos2 < MIN(rpos + cigar_oplen + 1 + 1, rend); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_INS_DP] += 1;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_INS_LEN] += mathsquare(cigar_oplen);
            }
            
            const auto & rtr1 = rtr_vec[MAX(1, rpos - region_offset) - 1];
            const auto & rtr2 = rtr_vec[MIN(rtr_vec.size() - 1, rpos - region_offset + 1)];
            for (int rpos2 = MAX((int)rpos - (int)rtr1.tracklen - 1 - 1, aln->core.pos); rpos2 < rpos; rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_RTR_INS_DP] += 1;
            }           
            for (int rpos2 = (int)rpos; rpos2 < MIN(rpos + rtr2.tracklen + 1 + 1, rend); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_RTR_INS_DP] += 1;
            }
            seg_format_prep_sets.getRefByPos(rpos)[SEG_a_AT_INS_DP] += 1;
            qpos += cigar_oplen;
        } else if (cigar_op == BAM_CDEL) {
            for (int rpos2 = (int)rpos; rpos2 < rpos + cigar_oplen; rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DP] += 1;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_aBQ2_DP] += 1;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_XM] += xm150;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_GO] += go150;
                if (isrc) { 
                    seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_LI] += MIN(rpos - frag_pos_L + 1, MAX_INSERT_SIZE); 
                    seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_LIDP] += 1;
                } else { 
                    seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_RI] += MIN(frag_pos_R - rpos    , MAX_INSERT_SIZE); 
                    seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_RIDP] += 1;
                }
                unsigned int ldist = rpos - aln->core.pos + 1;
                unsigned int rdist = rend - rpos;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_L_DIST_SUM] += ldist;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_R_DIST_SUM] += rdist;
               
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_INSLEN_SUM] += inslen_sum;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DELLEN_SUM] += dellen_sum;

                unsigned int lbaq = baq_offsetarr.getByPos(rpos) - baq_offsetarr.getByPos(aln->core.pos) + 1;
                unsigned int rbaq = baq_offsetarr.getByPos(rend-1) - baq_offsetarr.getByPos(rpos) + 1;
                seg_format_prep_sets.getRefByPos(rpos)[SEG_a_L_BAQ_SUM] += lbaq;
                seg_format_prep_sets.getRefByPos(rpos)[SEG_a_R_BAQ_SUM] += rbaq;
                
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_INSBAQ_SUM] += insbaq_sum;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_DELBAQ_SUM] += delbaq_sum;
            }
            
            for (int rpos2 = MAX((int)rpos - (int)cigar_oplen - 1, aln->core.pos); rpos2 < MIN(rpos + cigar_oplen + 1, rend); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_DEL_DP] += 1;
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_DEL_LEN] += mathsquare(cigar_oplen);
            }
            const auto & rtr1 = rtr_vec[MAX(1, rpos - region_offset) - 1];
            const auto & rtr2 = rtr_vec[MIN(rtr_vec.size() - 1, rpos - region_offset + 1)];
            for (int rpos2 = MAX((int)rpos - (int)rtr1.tracklen - 1, aln->core.pos); rpos2 < rpos; rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_RTR_DEL_DP] += 1;
            }
            for (int rpos2 = (int)rpos; rpos2 < MIN(rpos + rtr2.tracklen + 1, rend); rpos2++) {
                seg_format_prep_sets.getRefByPos(rpos2)[SEG_a_NEAR_RTR_DEL_DP] += 1;
            }
            seg_format_prep_sets.getRefByPos(rpos)[SEG_a_AT_DEL_DP] += 1;
            rpos += cigar_oplen;
        } else {
            process_cigar(qpos, rpos, cigar_op, cigar_oplen);
        }
    }
    assert(bam_endpos(aln) == rpos || !fprintf(stderr, "%u == %u failed for bam %s at tid %d position %d", 
            bam_endpos(aln), rpos, bam_get_qname(aln), aln->core.tid, aln->core.pos)); 
}

int
update_seg_format_thres_from_prep_sets(
        SegFormatThresSets & seg_format_thres_sets, 
        const SegFormatPrepSets & seg_format_prep_sets, 
        const unsigned int specialflag) {
    assert(seg_format_thres_sets.getIncluBegPosition() == seg_format_prep_sets.getIncluBegPosition());
    assert(seg_format_thres_sets.getExcluEndPosition() == seg_format_prep_sets.getExcluEndPosition());
    for (size_t epos = seg_format_prep_sets.getIncluBegPosition(); epos != seg_format_prep_sets.getExcluEndPosition(); epos++) {
        const unsigned int seg_a_dp = MAX(seg_format_prep_sets.getByPos(epos)[SEG_a_DP], 1U);
        const auto & p = seg_format_prep_sets.getByPos(epos);
        auto & t = seg_format_thres_sets.getRefByPos(epos);
        
        auto segLIDP = MAX(seg_format_prep_sets.getByPos(epos)[SEG_a_LIDP], 1);
        auto segRIDP = MAX(seg_format_prep_sets.getByPos(epos)[SEG_a_RIDP], 1);
        const auto ins_border_len = (unsigned int)ceil(sqrt((p[SEG_a_NEAR_INS_LEN]) / MAX(p[SEG_a_NEAR_INS_DP], 1)));
        const auto del_border_len = (unsigned int)ceil(sqrt((p[SEG_a_NEAR_DEL_LEN]) / MAX(p[SEG_a_NEAR_DEL_DP], 1)));
        t[SEG_aLPxT] =    (ins_border_len)                 * 5/4 + 5;
        t[SEG_aRPxT] = MAX(ins_border_len, del_border_len) * 5/4 + 5;

        t[SEG_aXM1T] = (p[SEG_a_XM] * 5 / (seg_a_dp*10) + (30+5)); // easier to pass
        t[SEG_aXM2T] = (p[SEG_a_XM] * 7 / (seg_a_dp*10) + (15+5)); // the higher the stronger the bias
        t[SEG_aGO1T] = (p[SEG_a_GO] * 5 / (seg_a_dp*10) + (20+5)); // easier to pass
        t[SEG_aGO2T] = (p[SEG_a_GO] * 7 / (seg_a_dp*10) + (10+5)); // the higher the stronger the bias
        
        t[SEG_aLI1T] = (unsigned int)((uint64_t)p[SEG_a_LI] * 250/100 / segLIDP);
        t[SEG_aLI2T] = (unsigned int)((uint64_t)p[SEG_a_LI] * 200/100 / segLIDP); // higher > stronger bias
        t[SEG_aRI1T] = (unsigned int)((uint64_t)p[SEG_a_RI] * 250/100 / segRIDP);
        t[SEG_aRI2T] = (unsigned int)((uint64_t)p[SEG_a_RI] * 200/100 / segRIDP); // higher > stronger bias
        t[SEG_aLI1t] = (unsigned int)((uint64_t)p[SEG_a_LI] *  50/100 / segLIDP);
        t[SEG_aLI2t] = (unsigned int)((uint64_t)p[SEG_a_LI] *  67/100 / segLIDP); // lower > stronger bias
        t[SEG_aRI1t] = (unsigned int)((uint64_t)p[SEG_a_RI] *  50/100 / segRIDP);
        t[SEG_aRI2t] = (unsigned int)((uint64_t)p[SEG_a_RI] *  67/100 / segRIDP); // lower > stronger bias
        
        t[SEG_aLP1t] = (unsigned int)non_neg_minus(p[SEG_a_L_DIST_SUM] / MAX(1, p[SEG_aBQ2_DP]), p[SEG_a_DELLEN_SUM] / MAX(1, p[SEG_aBQ2_DP]) + 10);
        t[SEG_aLP2t] = (unsigned int)non_neg_minus(p[SEG_a_L_DIST_SUM] / MAX(1, p[SEG_aBQ2_DP]), 5);
        t[SEG_aRP1t] = (unsigned int)non_neg_minus(p[SEG_a_R_DIST_SUM] / MAX(1, p[SEG_aBQ2_DP]), p[SEG_a_DELLEN_SUM] / MAX(1, p[SEG_aBQ2_DP]) + 10);
        t[SEG_aRP2t] = (unsigned int)non_neg_minus(p[SEG_a_R_DIST_SUM] / MAX(1, p[SEG_aBQ2_DP]), 5);
        
        t[SEG_aLB1t] = (unsigned int)non_neg_minus(p[SEG_a_L_BAQ_SUM] / MAX(1, p[SEG_aBQ2_DP]), p[SEG_a_DELBAQ_SUM] / MAX(1, p[SEG_aBQ2_DP]) + 10*5);
        t[SEG_aLB2t] = (unsigned int)non_neg_minus(p[SEG_a_L_BAQ_SUM] / MAX(1, p[SEG_aBQ2_DP]), 5*5);
        t[SEG_aRB1t] = (unsigned int)non_neg_minus(p[SEG_a_R_BAQ_SUM] / MAX(1, p[SEG_aBQ2_DP]), p[SEG_a_DELBAQ_SUM] / MAX(1, p[SEG_aBQ2_DP]) + 10*5);
        t[SEG_aRB2t] = (unsigned int)non_neg_minus(p[SEG_a_R_BAQ_SUM] / MAX(1, p[SEG_aBQ2_DP]), 5*5);
    }
}

template <bool isGap>
inline const
int
dealwith_segbias(
        const auto bq,
        const unsigned int rpos,
        auto & symbol_to_seg_format_depth_set,
        auto & symbol_to_VQ_format_tag_set,
        const auto & seg_format_thres_set,
        const bam1_t *aln,
        const auto xm150,
        const auto go150,
        const auto bm150,
        const CoveredRegion<int64_t> & baq_offsetarr,
        const unsigned int dist_to_interfering_indel,
        const unsigned int specialflag) {
    
    const auto rend = bam_endpos(aln); 
    
    const size_t seg_l_baq = baq_offsetarr.getByPos(rpos) - baq_offsetarr.getByPos(aln->core.pos) + 1;
    const size_t seg_r_baq = baq_offsetarr.getByPos(rend-1) - baq_offsetarr.getByPos(rpos) + 1;
    
    const size_t seg_l_nbases = (rpos - aln->core.pos + 1);
    const size_t seg_r_nbases = (bam_endpos(aln) - rpos);
    const size_t seg_region_size = ((bam_endpos(aln) - aln->core.pos) + 4) / 10;
    const size_t frag_pos_L = ((aln->core.isize != 0) ? MIN(aln->core.pos, aln->core.mpos) : aln->core.pos);
    const size_t frag_pos_R = ((aln->core.isize != 0) ? (frag_pos_L + aln->core.isize) : rend);
    const size_t frag_l_nbases1 = rpos - frag_pos_L + 1;
    const size_t frag_l_nbases2 = MIN(frag_l_nbases1, MAX_INSERT_SIZE);
    const size_t frag_r_nbases1 = frag_pos_R - rpos;
    const size_t frag_r_nbases2 = MIN(frag_r_nbases1, MAX_INSERT_SIZE);
    
    const bool is_normal = ((aln->core.isize != 0) || (0 == (aln->core.flag & 0x1)));
    const bool isrc = ((aln->core.flag & 0x10) == 0x10);
    const bool isr2 = ((aln->core.flag & 0x80) == 0x80 && (aln->core.flag & 0x1) == 0x1);
    const bool strand = (isrc ^ isr2);
    
    const auto a1BQ = (isrc ? VQ_a1BQr : VQ_a1BQf);
    const auto a2BQ = (isrc ? VQ_a2BQr : VQ_a2BQf);
    symbol_to_VQ_format_tag_set[a1BQ] += bq; 
    symbol_to_VQ_format_tag_set[a2BQ] += bq * bq / SQR_QUAL_DIV;
    symbol_to_VQ_format_tag_set[VQ_a1XM] += xm150; 
    
    auto & symbol_to_seg_aDP_depth_set = (strand
        ? (isrc ? symbol_to_seg_format_depth_set[SEG_aDPrr] : symbol_to_seg_format_depth_set[SEG_aDPrf]) 
        : (isrc ? symbol_to_seg_format_depth_set[SEG_aDPfr] : symbol_to_seg_format_depth_set[SEG_aDPff]));
    symbol_to_seg_aDP_depth_set += 1;
    
    const unsigned int const_XM1T = seg_format_thres_set[SEG_aXM1T];
    const unsigned int const_XM2T = seg_format_thres_set[SEG_aXM2T];
    const unsigned int const_GO1T = seg_format_thres_set[SEG_aGO1T];
    const unsigned int const_GO2T = seg_format_thres_set[SEG_aGO2T];
    
    const unsigned int const_LPxT = seg_format_thres_set[SEG_aLPxT]; 
    const unsigned int const_RPxT = seg_format_thres_set[SEG_aRPxT]; 
    if (isGap) {
        
        if (xm150 <= const_XM1T && go150 <= const_GO1T && bq >= 0*25) {
            symbol_to_seg_format_depth_set[SEG_aPF1] += 100;
        } else {
            symbol_to_seg_format_depth_set[SEG_aPF1] += MIN3(
                    100 * mathsquare(const_XM1T) / MAX(1, mathsquare(xm150)), 
                    100 * mathsquare(const_GO1T) / MAX(1, mathsquare(go150)),
                    100 * mathsquare(bq) / (25*25));
        } 
        if (xm150 <= const_XM2T && go150 <= const_GO2T && bq >= 0*30) {
            symbol_to_seg_format_depth_set[SEG_aPF2] += 100;
        } else {
            symbol_to_seg_format_depth_set[SEG_aPF2] += MIN3(
                    100 * mathsquare(const_XM2T) / MAX(1, mathsquare(xm150)), 
                    100 * mathsquare(const_GO2T) / MAX(1, mathsquare(go150)),
                    100 * mathsquare(bq) / (30*30));
        }
        symbol_to_seg_format_depth_set[SEG_aXM1] += (xm150 > 30 ? (100 * mathsquare(30) / mathsquare(xm150)) : 100);
    } else {
        unsigned int ampfact1 = 100;
        unsigned int ampfact2 = 100;

        if (xm150 > const_XM2T) {
            ampfact1 = 100 * mathsquare(const_XM2T) / mathsquare(xm150);
        } else { 
            ampfact1 = 100; 
        }
        if (bq < 25) {
            ampfact2 = 100 * mathsquare(bq) / (25*25);
        } else { 
            ampfact2 = 100; 
        }
        symbol_to_seg_format_depth_set[SEG_aPF1] += (ampfact1 * ampfact2 / (100));
        
        if (xm150 > const_XM1T) {
            ampfact1 = 100 * mathsquare(const_XM1T) / mathsquare(xm150);
        } else { 
            ampfact1 = 100; 
        }
        if (bq < 30) {
            ampfact2 = 100 * mathsquare(bq) / (30*30);
        } else { 
            ampfact2 = 100; 
        }
        symbol_to_seg_format_depth_set[SEG_aPF2] += (ampfact1 * ampfact2 / (100));
        
        symbol_to_seg_format_depth_set[SEG_aXM1] += (xm150 > 30 ? (100 * mathsquare(30) / mathsquare(xm150)) : 100);
        symbol_to_seg_format_depth_set[SEG_aXM2] += (xm150 > 20 ? (100 * mathsquare(20) / mathsquare(xm150)) : 100);
        symbol_to_seg_format_depth_set[SEG_aBM2] += (bm150 > 20 ? (100 * mathsquare(20) / mathsquare(bm150)) : 100);
    }
    if (((!isGap) && bq >= 20) || (isGap && dist_to_interfering_indel >= 20)) {
        const bool is_l1_unbiased = (seg_l_nbases >= seg_format_thres_set[SEG_aLP1t]);
        const bool is_l2_unbiased = (seg_l_nbases >= seg_format_thres_set[SEG_aLP2t]);
        const bool is_r1_unbiased = (seg_r_nbases >= seg_format_thres_set[SEG_aRP1t]);
        const bool is_r2_unbiased = (seg_r_nbases >= seg_format_thres_set[SEG_aRP2t]);
        if ((seg_l_nbases >= const_LPxT) && (seg_r_nbases >= const_RPxT)) {
            if (is_l1_unbiased) {
                symbol_to_seg_format_depth_set[SEG_aLP1] += 1;
            }
            if (is_l2_unbiased) {
                symbol_to_seg_format_depth_set[SEG_aLP2] += 1;
            }
            if (is_r1_unbiased) {
                symbol_to_seg_format_depth_set[SEG_aRP1] += 1;
            }
            if (is_r2_unbiased) {
                symbol_to_seg_format_depth_set[SEG_aRP2] += 1;
            }
            symbol_to_seg_format_depth_set[SEG_aLPL] += seg_l_nbases;
            symbol_to_seg_format_depth_set[SEG_aRPL] += seg_r_nbases;
        }
        if (seg_l_baq >= 23 && seg_r_baq >= 23) {
            if (seg_l_baq >= 45-10) {
                symbol_to_seg_format_depth_set[SEG_aLB1] += 1;
            }
            if (seg_l_baq >= 60-10) {
                symbol_to_seg_format_depth_set[SEG_aLB2] += 1;
            }
            if (seg_r_baq >= 45-10) {
                symbol_to_seg_format_depth_set[SEG_aRB1] += 1;
            }
            if (seg_r_baq >= 60-10) {
                symbol_to_seg_format_depth_set[SEG_aRB2] += 1;
            }
            symbol_to_seg_format_depth_set[SEG_aLBL] += seg_l_baq;
            symbol_to_seg_format_depth_set[SEG_aRBL] += seg_r_baq;
        }
        symbol_to_seg_format_depth_set[SEG_aBQ2] += 1;
    }
    
    if (0 == strand) {
        symbol_to_seg_format_depth_set[SEG_aLILf] += frag_l_nbases2;
        symbol_to_seg_format_depth_set[SEG_aRILf] += frag_r_nbases2;
    } else {
        symbol_to_seg_format_depth_set[SEG_aLILr] += frag_l_nbases2;
        symbol_to_seg_format_depth_set[SEG_aRILr] += frag_r_nbases2;
    }
    
    const bool is_l_nonbiased = ((0 != (aln->core.flag & 0x8)) && seg_l_nbases > seg_r_nbases);
    const bool is_r_nonbiased = ((0 != (aln->core.flag & 0x8)) && seg_l_nbases < seg_r_nbases);
    if (isrc) {
        auto dist2iend = frag_l_nbases2;
        if ((dist2iend <= seg_format_thres_set[SEG_aLI1T] || isGap) && dist2iend >= seg_format_thres_set[SEG_aLI1t] && (is_normal || isGap) || is_l_nonbiased) {
            symbol_to_seg_format_depth_set[SEG_aLI1] += 1;
        } 
        if ((dist2iend <= seg_format_thres_set[SEG_aLI2T] || isGap) && dist2iend >= seg_format_thres_set[SEG_aLI2t] && (is_normal || isGap)) {
            symbol_to_seg_format_depth_set[SEG_aLI2] += 1;
        }
    } else {
        auto dist2iend = frag_r_nbases2;
        if ((dist2iend <= seg_format_thres_set[SEG_aRI1T] || isGap) && dist2iend >= seg_format_thres_set[SEG_aRI1t] && (is_normal || isGap) || is_r_nonbiased) {
            symbol_to_seg_format_depth_set[SEG_aRI1] += 1;
        } 
        if ((dist2iend <= seg_format_thres_set[SEG_aRI2T] || isGap) && dist2iend >= seg_format_thres_set[SEG_aRI2t] && (is_normal || isGap)) {
            symbol_to_seg_format_depth_set[SEG_aRI2] += 1;
        }
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
        assert(this->getIncluBegPosition() <= other.getIncluBegPosition() || !fprintf(stderr, "%lu <= %lu failed!", this->getIncluBegPosition(), other.getIncluBegPosition()));
        assert(this->getExcluEndPosition() >= other.getExcluEndPosition() || !fprintf(stderr, "%lu >= %lu failed!", this->getExcluEndPosition(), other.getExcluEndPosition())); 
    }
    // mainly for grouping reads sharing the same UMI into one family
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
    template<bool TIndelIsMajor = false> 
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
    int
    updateByFiltering(
            const GenericSymbol2CountCoverage<TSymbol2Count> &other, 
            const std::array<unsigned int, NUM_SYMBOL_TYPES> thres,
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
        size_t ipos = epos;
        posToIndelToCount_inc(this->getRefPosToIseqToData(symbol), ipos, iseq, incvalue);
    };

    void // GenericSymbol2CountCoverage<TSymbol2Count>::
    incDel(const unsigned int epos, const uint32_t dlen, const AlignmentSymbol symbol, const uint32_t incvalue = 1) {
        assert (incvalue > 0);
        assert (dlen > 0);
        size_t ipos = epos;
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
    
    template<bool TIsProton, ValueType TUpdateType, bool TIsBiasUpdated>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByAln(
            const bam1_t *const aln, 
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhredArg, 
            unsigned int frag_indel_basemax,
            const unsigned int region_offset,
            const auto & region_symbolvec, 
            const std::vector<RegionalTandemRepeat> & region_repeatvec,
            const CoveredRegion<int64_t> & baq_offsetarr,
            
            auto & seg_format_depth_sets,
            auto & symbol_to_VQ_format_tag_sets,
            const auto & seg_format_prep_sets,
            const auto & seg_format_thres_sets,
            
            unsigned int primerlen,
            uint32_t sepcialflag) {
        
        static_assert(BASE_QUALITY_MAX == TUpdateType || SYMBOL_COUNT_SUM == TUpdateType);
        assert(this->tid == SIGN2UNSIGN(aln->core.tid));
        assert(this->getIncluBegPosition() <= SIGN2UNSIGN(aln->core.pos)   || !fprintf(stderr, "%lu <= %d failed", this->getIncluBegPosition(), aln->core.pos));
        assert(this->getExcluEndPosition() >= SIGN2UNSIGN(bam_endpos(aln)) || !fprintf(stderr, "%lu >= %d failed", this->getExcluEndPosition(), bam_endpos(aln)));
        
        const uint32_t n_cigar = aln->core.n_cigar;
        const uint32_t *cigar = bam_get_cigar(aln);
        const uint8_t *bseq = bam_get_seq(aln);
        const uint32_t rend = bam_endpos(aln);
        const auto symbolType2addPhred = symbolType2addPhredArg;
        // const auto roffset = bq_baq_sum.at(0).getIncluBegPosition();
        
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
        const unsigned int xm150 = xm_cnt * 1500 / (rend - aln->core.pos);
        const unsigned int go150 = ngo_cnt * 1500 / (rend - aln->core.pos);
        
        std::vector<unsigned int> indel_rposs = {{ 0 }};
        unsigned int indel_rposs_idx = 0;
        std::array<unsigned int, NUM_ALIGNMENT_SYMBOLS> bm_cnts = {{ 0 }}; // mismatch of the same base type
        {
            unsigned int qpos = 0;
            unsigned int rpos = aln->core.pos;
            for (unsigned int i = 0; i < n_cigar; i++) {
                uint32_t c = cigar[i];
                unsigned int cigar_op = bam_cigar_op(c);
                unsigned int cigar_oplen = bam_cigar_oplen(c);
                if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                    for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                        unsigned int base4bit = bam_seqi(bseq, qpos);
                        unsigned int base3bit = seq_nt16_int[base4bit];
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
                    for (unsigned int qpos2 = qpos; qpos2 < qpos + cigar_oplen; qpos2++) {
                        if (bam_phredi(aln, qpos2) < 20) { is_lowBQ = true; }
                    }
                    if (is_lowBQ) { indel_rposs.push_back(rpos); }
                    qpos += cigar_oplen;
                } else if (cigar_op == BAM_CDEL) {
                    bool is_lowBQ = (MIN(bam_phredi(aln, MAX(1, qpos) - 1), bam_phredi(aln, qpos)) < 20);
                    if (is_lowBQ) { indel_rposs.push_back(rpos); }
                    rpos += cigar_oplen;
                } else {
                    process_cigar(qpos, rpos, cigar_op, cigar_oplen);
                }
            }
            indel_rposs.push_back(INT32_MAX);
        }
        std::array<unsigned int, NUM_ALIGNMENT_SYMBOLS> bm150s = {{ 0 }};
        for (unsigned int i = 0; i < NUM_ALIGNMENT_SYMBOLS; i++) {
            bm150s[i] = bm_cnts[i] * 1500 / (rend - aln->core.pos);
        }
        
        const bool isrc = ((aln->core.flag & 0x10) == 0x10);
        const unsigned int ibeg = ((aln->core.isize != 0) ? (MIN(aln->core.pos, aln->core.mpos) + primerlen) : (isrc ? 0 : (aln->core.pos + primerlen)));
        const unsigned int iend = ((aln->core.isize != 0) ? non_neg_minus(MIN(aln->core.pos, aln->core.mpos) + abs(aln->core.isize),  primerlen) 
                : (isrc ? non_neg_minus(rend, primerlen) : INT32_MAX));
        
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
                            || !fprintf(stderr, "Bam line with QNAME %s has rpos %d which is not the in range (%d - %d)", 
                            bam_get_qname(aln), rpos, aln->core.pos, bam_endpos(aln)));
if ((0 == primerlen) || (ibeg <= rpos && rpos < iend)) {
                    if (i2 > 0) {
                        const auto noindel_phredvalue = MIN3(
                                region_repeatvec[rpos-region_offset - 1].indelphred,
                                region_repeatvec[rpos-region_offset].indelphred,
                                40);
                        incvalue = MIN(MIN(bam_phredi(aln, qpos-1), bam_phredi(aln, qpos)) + (symbolType2addPhred[LINK_SYMBOL]), noindel_phredvalue) + 1;
                        this->template inc<TUpdateType>(rpos, LINK_M, incvalue, aln);
                        if (TIsBiasUpdated) {
                            if (indel_rposs[indel_rposs_idx] <= rpos) {
                                indel_rposs_idx++;
                            }
                            unsigned int prev_indel_rpos = indel_rposs[indel_rposs_idx-1];
                            unsigned int next_indel_rpos = indel_rposs[indel_rposs_idx];
                            dealwith_segbias<true>(
                                    incvalue,
                                    rpos,
                                    seg_format_depth_sets.getRefByPos(rpos)[LINK_M],
                                    symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[LINK_M],
                                    seg_format_thres_sets.getByPos(rpos),
                                    aln,
                                    xm150,
                                    go150,
                                    bm150s[LINK_M],
                                    baq_offsetarr,
                                    MIN(rpos - prev_indel_rpos, next_indel_rpos - rpos),
                                    0);
                        }
                    }
                    unsigned int base4bit = bam_seqi(bseq, qpos);
                    unsigned int base3bit = seq_nt16_int[base4bit];
                    AlignmentSymbol symbol = AlignmentSymbol(base3bit);
                    incvalue = bam_phredi(aln, qpos) + symbolType2addPhred[BASE_SYMBOL];
                    this->template inc<TUpdateType>(rpos, AlignmentSymbol(base3bit), incvalue, aln);
                    if (TIsBiasUpdated) {
                        dealwith_segbias<false>(
                                incvalue,
                                rpos,
                                seg_format_depth_sets.getRefByPos(rpos)[symbol],
                                symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                                seg_format_thres_sets.getByPos(rpos),
                                aln,
                                xm150,
                                go150,
                                bm150s[symbol],
                                baq_offsetarr,
                                10000,
                                0);
                    }
}
                    rpos += 1;
                    qpos += 1;
                }
            } else if (cigar_op == BAM_CINS) {
if ((0 == primerlen) || (ibeg <= rpos && rpos < iend)) {
                const bool is_ins_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(aln->core.l_qseq));
                unsigned int inslen = SIGN2UNSIGN(cigar_oplen);
                if (is_ins_at_read_end) {
                    LOG(logWARNING) << "Query " << bam_get_qname(aln) << " has insertion of legnth " << cigar_oplen << " at " << qpos
                            << " which is not exclusively between 0 and " << aln->core.l_qseq << " aligned to tid " << aln->core.tid << " and position " << rpos;
                    incvalue = (0 != qpos ? bam_phredi(aln, qpos-1) : 
                            ((qpos + cigar_oplen < SIGN2UNSIGN(aln->core.l_qseq)) ? 
                            bam_phredi(aln, qpos + SIGN2UNSIGN(cigar_oplen)) : 1)) 
                            + (symbolType2addPhred[LINK_SYMBOL]);
                } else {
                    unsigned int phredvalue = ref_to_phredvalue(
                            inslen, 
                            region_symbolvec, 
                            rpos - region_offset,
                            frag_indel_basemax, 
                            8.0, 
                            cigar_oplen, 
                            cigar_op, 
                            aln);
                    unsigned int thisdp = (seg_format_prep_sets.getByPos(rpos)[SEG_a_AT_INS_DP]); 
                    unsigned int neardp = (MAX(seg_format_prep_sets.getByPos(rpos)[SEG_a_NEAR_INS_DP], seg_format_prep_sets.getByPos(rpos)[SEG_a_NEAR_RTR_INS_DP]));
                    unsigned int insbase_minphred = 80;
                    for (unsigned int qpos2 = qpos; qpos2 < qpos + cigar_oplen; qpos2++) {
                        UPDATE_MIN(insbase_minphred, bam_phredi(aln, qpos2));
                    }
                    unsigned int qfromBQ1 = insbase_minphred; //  + (TIsProton ? proton_cigarlen2phred(cigar_oplen) : 0);
                    unsigned int qfromBQ2 = ((thisdp * 2 < neardp) ? qfromBQ1 : (TIsProton ? (qfromBQ1 + proton_cigarlen2phred(cigar_oplen)) : 80));
                    incvalue = MIN(qfromBQ2, phredvalue + symbolType2addPhred[LINK_SYMBOL]) + 1;
                }
                if (!is_ins_at_read_end) {
                    const auto symbol = insLenToSymbol(inslen);
                    this->template inc<TUpdateType>(rpos, symbol, MAX(SIGN2UNSIGN(1), incvalue), aln);
                    if (TIsBiasUpdated) {
                        dealwith_segbias<true>(
                                MAX(SIGN2UNSIGN(1), incvalue),
                                rpos,
                                seg_format_depth_sets.getRefByPos(rpos)[symbol],
                                symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                                seg_format_thres_sets.getByPos(rpos),
                                aln,
                                xm150,
                                go150,
                                bm150s[symbol],
                                baq_offsetarr,
                                10000,
                                0);
                    }
                    std::string iseq;
                    iseq.reserve(cigar_oplen);
                    unsigned int incvalue2 = incvalue;
                    for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                        unsigned int base4bit = bam_seqi(bseq, qpos+i2);
                        const char base8bit = seq_nt16_str[base4bit];
                        iseq.push_back(base8bit);
                        incvalue2 = MIN(incvalue2, SIGN2UNSIGN(bam_seqi(bseq, qpos+i2)))
                                + (symbolType2addPhred[LINK_SYMBOL]);
                    }
                    if (867113 <= rpos && rpos <= 867115) {
                        LOG(logINFO) << " rpos " << rpos << " has ins-seq " << iseq;
                    }
                    this->incIns(rpos, iseq, insLenToSymbol(inslen), MAX(SIGN2UNSIGN(1), incvalue2));
                }
}
                qpos += cigar_oplen;
            } else if (cigar_op == BAM_CDEL) {
if ((0 == primerlen) || (ibeg <= rpos && rpos < iend)) {
                const bool is_del_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(aln->core.l_qseq));
                unsigned int dellen = SIGN2UNSIGN(cigar_oplen);
                if (is_del_at_read_end) {
                    LOG(logWARNING) << "Query " << bam_get_qname(aln) << " has deletion of legnth " << cigar_oplen << " at " << qpos
                            << " which is not exclusively between 0 and " << aln->core.l_qseq << " aligned to tid " << aln->core.tid << " and position " << rpos; 
                    incvalue = (0 != qpos ? bam_phredi(aln, qpos-1) : 
                            ((qpos + cigar_oplen < SIGN2UNSIGN(aln->core.l_qseq)) ? 
                            bam_phredi(aln, qpos + SIGN2UNSIGN(cigar_oplen)) : 1))
                            + (symbolType2addPhred[LINK_SYMBOL]);
                } else {
                    unsigned int phredvalue = ref_to_phredvalue(
                            dellen, 
                            region_symbolvec, 
                            rpos - region_offset, 
                            frag_indel_basemax, 
                            8.0, 
                            cigar_oplen, 
                            cigar_op, 
                            aln);
                    unsigned int thisdp = (seg_format_prep_sets.getByPos(rpos)[SEG_a_AT_DEL_DP]);
                    unsigned int neardp = (MAX(seg_format_prep_sets.getByPos(rpos)[SEG_a_NEAR_DEL_DP], seg_format_prep_sets.getByPos(rpos)[SEG_a_NEAR_RTR_DEL_DP]));
                    unsigned int qfromBQ1 = MIN(bam_phredi(aln, qpos), bam_phredi(aln, qpos-1));
                    unsigned int qfromBQ2 = ((thisdp * 2 < neardp) ? qfromBQ1 : (TIsProton ? (qfromBQ1 + proton_cigarlen2phred(cigar_oplen)) : 80));
                    incvalue = MIN(qfromBQ2, phredvalue + symbolType2addPhred[LINK_SYMBOL]) + 1;
                }
                if (!is_del_at_read_end) {
                    AlignmentSymbol symbol = delLenToSymbol(dellen);
                    this->template inc<TUpdateType>(rpos, symbol, MAX(SIGN2UNSIGN(1), incvalue), aln);
                    if (TIsBiasUpdated) {
                        dealwith_segbias<true>(
                                MAX(SIGN2UNSIGN(1), incvalue),
                                rpos,
                                seg_format_depth_sets.getRefByPos(rpos)[symbol],
                                symbol_to_VQ_format_tag_sets.getRefByPos(rpos)[symbol],
                                seg_format_thres_sets.getByPos(rpos),
                                aln,
                                xm150,
                                go150,
                                bm150s[symbol],
                                baq_offsetarr,
                                10000,
                                0);
                    }
                    this->incDel(rpos, cigar_oplen, delLenToSymbol(dellen), MAX(SIGN2UNSIGN(1), incvalue));
#if 1 
// The definition of non-ref for deletion is not clearly defined, this piece of code makes sure that deletion identity conforms to the VCF specs.
                    for (unsigned int rpos2 = rpos; rpos2 < MIN(rpos + cigar_oplen, rend); rpos2++) {
                        for (AlignmentSymbol s : std::array<AlignmentSymbol, 2> {{ BASE_NN, LINK_NN }} ) {
                            const auto p = ((BASE_NN == s) ? (rpos2) : (rpos2 + 1));
                            if (p >= rend) { continue ; }
                            this->template inc<TUpdateType>(p, s, MAX(SIGN2UNSIGN(1), incvalue), aln);
                            if (TIsBiasUpdated) {
                                if (indel_rposs[indel_rposs_idx] <= rpos) {
                                    indel_rposs_idx++;
                                }
                                unsigned int prev_indel_rpos = indel_rposs[indel_rposs_idx-1];
                                unsigned int next_indel_rpos = indel_rposs[indel_rposs_idx];
                                dealwith_segbias<true>(
                                        MAX(SIGN2UNSIGN(1), incvalue),
                                        p,
                                        seg_format_depth_sets.getRefByPos(p)[s],
                                        symbol_to_VQ_format_tag_sets.getRefByPos(p)[s],
                                        seg_format_thres_sets.getByPos(p),
                                        aln,
                                        xm150,
                                        go150,
                                        bm150s[s],
                                        baq_offsetarr,
                                        MIN(rpos - prev_indel_rpos, next_indel_rpos - rpos),
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
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred, 
            const unsigned int frag_indel_basemax, 
            const bool is_proton, 
            const unsigned int region_offset,
            const auto & region_symbolvec,
            const auto & region_repeatvec,
            const auto & baq_offsetarr,
            
            auto & symbol_to_seg_format_depth_sets,
            auto & symbol_to_VQ_format_tag_sets,
            const auto & seg_format_prep_sets,
            const auto & seg_format_thres_sets,
            
            unsigned int primerlen,
            unsigned int specialflag) {
        for (const bam1_t *aln : aln_vec) {
            if (is_proton) {
                this->template updateByAln<true, TUpdateType, TIsBiasUpdated>(
                        aln, 
                        symbolType2addPhred, 
                        frag_indel_basemax, 
                        region_offset, 
                        region_symbolvec, 
                        region_repeatvec, 
                        baq_offsetarr,

                        symbol_to_seg_format_depth_sets,
                        symbol_to_VQ_format_tag_sets,
                        seg_format_prep_sets,
                        seg_format_thres_sets,
                        
                        primerlen,
                        0);
            } else {
                this->template updateByAln<false, TUpdateType, TIsBiasUpdated>(
                        aln, 
                        symbolType2addPhred, 
                        frag_indel_basemax, 
                        region_offset, 
                        region_symbolvec, 
                        region_repeatvec, 
                        baq_offsetarr,

                        symbol_to_seg_format_depth_sets,
                        symbol_to_VQ_format_tag_sets,
                        seg_format_prep_sets,
                        seg_format_thres_sets,
                        
                        primerlen,
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
    
    std::array<Symbol2Bucket2CountCoverage, 2> dedup_ampDistr;
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
            const auto & region_repeatvec,
            const auto & baq_offsetarr,

            unsigned int primerlen,
            unsigned int specialflag) {
        
        Symbol2CountCoverage bg_seg_bqsum_conslogo(tid, this->getUnifiedIncluBegPosition(), this->getUnifiedExcluEndPosition());
        for (const auto & alns2pair2dflag : alns3) {
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair2dflag.first[strand];
                for (const auto & alns1 : alns2) {
                    for (const bam1_t *aln : alns1) {
                        update_seg_format_prep_sets_by_aln(
                                this->seg_format_prep_sets,
                                aln,
                                region_repeatvec,
                                baq_offsetarr,
                                this->getUnifiedIncluBegPosition(),
                                0);
                    }
                }
            }
        }
        update_seg_format_thres_from_prep_sets(
                this->seg_format_thres_sets,
                this->seg_format_prep_sets,
                0);
        for (const auto & alns2pair2dflag : alns3) {
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair2dflag.first[strand];
                for (const auto & alns1 : alns2) {
                    bg_seg_bqsum_conslogo.updateByRead1Aln<SYMBOL_COUNT_SUM, true>(
                            alns1, 
                            symbolType2addPhred, 
                            frag_indel_basemax, 
                            is_proton, 
                            this->getUnifiedIncluBegPosition(), 
                            region_symbolvec, 
                            region_repeatvec,
                            baq_offsetarr,

                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,

                            primerlen,
                            0);
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
                            region_repeatvec,
                            baq_offsetarr,                            

                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,
                            
                            primerlen,
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
                            
                            int max_qual = 8 + get_avgBQ(bg_seg_bqsum_conslogo, symbol_to_seg_format_depth_sets, epos, con_symbol);
                            unsigned int phredlike = MIN(con_count * 2 - tot_count, max_qual);
                            
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
        assert(this->symbol_to_seg_format_depth_sets.getIncluBegPosition() == this->dedup_ampDistr[0].getIncluBegPosition());
        assert(this->symbol_to_seg_format_depth_sets.getExcluEndPosition() == this->dedup_ampDistr[0].getExcluEndPosition());
        for (unsigned int epos = this->getUnifiedIncluBegPosition(); epos < this->getUnifiedExcluEndPosition(); epos++) {
            for (SymbolType symbolType : SYMBOL_TYPE_ARR) {
                for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symbolType]) {
                    const auto totDP = 
                            formatSumBySymbolType(this->symbol_to_frag_format_depth_sets[0].getByPos(epos), 
                            symbolType, FRAG_bDP) + 
                            formatSumBySymbolType(this->symbol_to_frag_format_depth_sets[1].getByPos(epos), 
                            symbolType, FRAG_bDP); 
                    int max_qual = 8 + get_avgBQ(bg_seg_bqsum_conslogo, symbol_to_seg_format_depth_sets, epos, symbol);
                    int maxvqual = 0;
                    unsigned int argmaxAD = 0;
                    unsigned int argmaxBQ = 0;
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
                    if (2358481 == epos) {
                        LOG(logINFO) << "Symbol " << symbol << " bq distributions:";
                        for (auto x : this->dedup_ampDistr[0].getByPos(epos).getSymbolCounts(symbol)) {
                            LOG(logINFO) << "\t" << x;
                        }
                    }
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
            const auto & baq_offsetarr,
            unsigned int primerlen,
            const unsigned int phred_indel_error_before_barcode_labeling,
            unsigned int specialflag = 0) {
        bool is_loginfo_enabled = true; 
        
        unsigned int niters = 0;
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            niters++;
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
                            baq_offsetarr,
                            
                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,
                            
                            primerlen,
                            0);
                    read_family_con_ampl.updateByFiltering(read_ampBQerr_fragWithR1R2, std::array<unsigned int, NUM_SYMBOL_TYPES> {{ 25, 0 }});
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
                        if (1 < tot_count && (con_count * 100 >= tot_count * 80)) {
                            this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP2] += 1;
                        }
                        if (2 < tot_count && (con_count * 100 >= tot_count * 85)) {
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
                        if (tot_count < 4) { continue; } 
                        if (con_count * 3 < tot_count * 2) { continue; }
                        for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symbolType]) {
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
                Symbol2CountCoverage read_family_mmm_ampl(tid2, beg2, end2);
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
                            baq_offsetarr,      
                            
                            this->symbol_to_seg_format_depth_sets,
                            this->symbol_to_VQ_format_tag_sets,
                            this->seg_format_prep_sets,
                            this->seg_format_thres_sets,
                            
                            primerlen,
                            0);
                    // The line below is similar to : read_family_mmm_ampl.updateByConsensus<SYMBOL_COUNT_SUM>(read_ampBQerr_fragWithR1R2);
                    read_family_con_ampl.updateByFiltering(read_ampBQerr_fragWithR1R2, std::array<unsigned int, NUM_SYMBOL_TYPES> {{ 25, 0 }});
                    read_family_mmm_ampl.updateByMajorMinusMinor(read_ampBQerr_fragWithR1R2, false); 
                }
                if ((2 == alns2pair2dflag.second) && alns2pair[0].size() > 0 && alns2pair[1].size() > 0) { // is duplex
                    read_duplex_amplicon.template updateByConsensus<false>(read_family_mmm_ampl);
                }
                std::basic_string<std::pair<unsigned int, AlignmentSymbol>> pos_symbol_string;
                for (size_t epos = read_family_mmm_ampl.getIncluBegPosition(); 
                        epos < read_family_mmm_ampl.getExcluEndPosition(); 
                        epos++) {
                    for (SymbolType symbolType : SYMBOL_TYPES_IN_VCF_ORDER) {
                        AlignmentSymbol con_symbol; // assign with the value AlignmentSymbol(NUM_ALIGNMENT_SYMBOLS) to flag for error
                        unsigned int con_count, tot_count;
                        read_family_mmm_ampl.getRefByPos(epos).fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        if (0 == tot_count) { continue; }
                        this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos)[con_symbol][FAM_cDP0] += 1; // in some rare cases, cDP1 > cDP0 
                        
                        unsigned int majorcount = this->symbol_to_fam_format_depth_sets_2strand[strand].getByPos(epos).at(con_symbol).at(FAM_cDPM);
                        unsigned int minorcount = this->symbol_to_fam_format_depth_sets_2strand[strand].getByPos(epos).at(con_symbol).at(FAM_cDPm);
                        unsigned int confam_qual = MAX(con_count * 2, tot_count + 1) - tot_count;
                        if (LINK_SYMBOL == symbolType) {
                            unsigned int con_count = read_family_con_ampl.getByPos(epos).getSymbolCount(con_symbol);
                            unsigned int tot_count = read_family_con_ampl.getByPos(epos).sumBySymbolType(symbolType);
                            double effective_nreads = log((double)(con_count + 1) / (double)(tot_count - con_count + 1)) / log(2); // PCR errors
                            double prior_weight = 1.0 / (minorcount + 1.0);
                            double realphred = prob2realphred((minorcount + prior_weight) / (majorcount + minorcount + prior_weight / phred2prob(confam_qual)));
                            confam_qual = (unsigned int)MIN(effective_nreads * realphred, phred_indel_error_before_barcode_labeling + realphred);
                        } else {
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
                            this->symbol_to_duplex_format_depth_sets.getRefByPos(epos)[con_symbol][DUPLEX_dDP1] += 1;
                        }
                        if (1 < tot_count) {
                            this->symbol_to_duplex_format_depth_sets.getRefByPos(epos)[con_symbol][DUPLEX_dDP2] += 1;
                        }
                    }
                }
            }
        }
        for (unsigned int strand = 0; strand < 2; strand++) {
            assert(this->symbol_to_fam_format_depth_sets_2strand[strand].getExcluEndPosition() == this->dedup_ampDistr[strand].getExcluEndPosition());
            assert(this->symbol_to_fam_format_depth_sets_2strand[strand].getIncluBegPosition() == this->dedup_ampDistr[strand].getIncluBegPosition()); 
            auto VQ_cIAQ = (strand ? VQ_cIAQr : VQ_cIAQf);
            auto VQ_cIAD = (strand ? VQ_cIADr : VQ_cIADf);
            auto VQ_cIDQ = (strand ? VQ_cIDQr : VQ_cIDQf);
            for (unsigned int epos = this->getUnifiedIncluBegPosition(); epos < this->getUnifiedExcluEndPosition(); epos++) {
                for (SymbolType symbolType : SYMBOL_TYPE_ARR) {
                    AlignmentSymbol ref_symbol = region_symbolvec[epos - this->getUnifiedIncluBegPosition()];
                    const auto totDP = formatSumBySymbolType(this->symbol_to_fam_format_depth_sets_2strand[strand].getRefByPos(epos), symbolType, FAM_cDP0);
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
            const CoveredRegion<int64_t> & baq_offsetarr,
            unsigned int bq_phred_added_misma, 
            unsigned int bq_phred_added_indel, 
            bool should_add_note, 
            const unsigned int frag_indel_basemax,
            const PhredMutationTable & phred_max_sscs_table, 
            bool use_deduplicated_reads, 
            const bool is_proton, 
            const unsigned int primerlen,
            const unsigned int phred_indel_error_before_barcode_labeling,
            const unsigned int baq_per_aligned_base,
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
                baq_offsetarr,
                primerlen,
                0);
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
                    baq_offsetarr,
                    primerlen,
                    phred_indel_error_before_barcode_labeling,
                    0);
            updateHapMap(mutform2count4map_fq, this->symbol_to_fam_format_depth_sets_2strand, FAM_cDP0);
        }
        return 0;
    };
};

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
};

void
fill_symbol_fmt(
        auto & fmtAD,
        const auto & symbol_to_abcd_format_depth_sets,
        const auto format_field,
        const auto refpos,
        const auto symbol,
        const unsigned int allele_idx) {
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
        int minABQ,
        const unsigned int specialflag) {
    const unsigned int a = 0;
    
    fill_symbol_fmt(fmt.a1XM,  symbol_to_VQ_format_tag_sets, VQ_a1XM,  refpos, symbol, a);
    
    fill_symbol_fmt(fmt.a1BQf, symbol_to_VQ_format_tag_sets, VQ_a1BQf, refpos, symbol, a);
    fill_symbol_fmt(fmt.a1BQr, symbol_to_VQ_format_tag_sets, VQ_a1BQr, refpos, symbol, a);
    
    const auto a2BQf = symbol_to_VQ_format_tag_sets.getByPos(refpos)[symbol][VQ_a2BQf];
    const auto a2BQr = symbol_to_VQ_format_tag_sets.getByPos(refpos)[symbol][VQ_a2BQr];
    int aDPf = fmt.aDPff[a] + fmt.aDPrf[a];
    int aDPr = fmt.aDPfr[a] + fmt.aDPrr[a];
    int ADP = fmt.ADPff[0] + fmt.ADPrf[0] + fmt.ADPfr[0] + fmt.ADPrr[0];
    const int rssDPfBQ = (int)(aDPf * sqrt(a2BQf * SQR_QUAL_DIV / MAX(1, aDPf)));
    const int rssDPrBQ = (int)(aDPr * sqrt(a2BQr * SQR_QUAL_DIV / MAX(1, aDPr)));
    
    const int rssDPbBQ = (int)((aDPf + aDPr) * sqrt((a2BQf + a2BQr) * SQR_QUAL_DIV / MAX(1, aDPf + aDPr)));
    
    assert ((aDPf + aDPr) * 100 >= LAST(fmt.aXM1));
    assert ((aDPf + aDPr) * 100 >= LAST(fmt.aXM2));
    assert ((aDPf + aDPr) * 100 >= LAST(fmt.aBM2));
    int64_t minABQa = (int)minABQ - (int)(50.0 * mathsquare(MAX(0, ((aDPf + aDPr + 0.5) * 2.0 / (ADP + 1.0) - 1.0))));
    const int dp10pc = 10;
    double sbratio = (double)(MAX(aDPf, aDPr) * 10 + dp10pc) / (double)(MIN(aDPf, aDPr) * 10 + dp10pc);
    minABQa += (int)BETWEEN((int)mathsquare(sbratio) - 5, 0, 40);
    const int xmratio = (int)(4*100 * (int)(aDPf + aDPr) / (int)MAX(1, LAST(fmt.aXM2)) - 0);
    const int bmratio = (int)(4*100 * (int)(aDPf + aDPr) / (int)MAX(1, LAST(fmt.aBM2)) - 0);
    minABQa += (int)BETWEEN((int)xmratio - 5, 0, 30) + (int)BETWEEN((int)bmratio - 5, 0, 40);
    fmt.note += std::string("//minABQa/") 
            + std::to_string(minABQa) + "/" 
            + std::to_string(sbratio) + "/" 
            + std::to_string(xmratio) + "/" 
            + std::to_string(bmratio) + "//";
    int a_BQ_syserr_qual_fw = (rssDPfBQ * 3 - minABQa * aDPf * 3 / 10 + rssDPrBQ - minABQa * aDPr / 10) / 3;
    int a_BQ_syserr_qual_rv = (rssDPrBQ * 3 - minABQa * aDPr * 3 / 10 + rssDPfBQ - minABQa * aDPf / 10) / 3;
    int a_BQ_syserr_qual_2d = (rssDPbBQ) - minABQa * (aDPf + aDPr) / 10;
    const int a_rmsBQ = (rssDPbBQ) / MAX(1, aDPf + aDPr);
    fill_symbol_fmt(fmt.bMQ,  symbol_to_VQ_format_tag_sets,  VQ_bMQ,  refpos, symbol, a);
    fmt.bMQ[a] = (unsigned int)floor(sqrt(fmt.bMQ[a] * SQR_QUAL_DIV / MAX(fmt.bDPf[a] + fmt.bDPr[a], 1)) + (double)(1.0 - FLT_EPSILON));
    
    clear_push(fmt.a2BQf, rssDPfBQ, a);
    clear_push(fmt.a2BQr, rssDPrBQ, a);
    clear_push(fmt.aBQ, a_rmsBQ);
    const int sys_seq_err_phred_per_ref_bp = 30; // https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451
    clear_push(fmt.aBQQ, MAX(a_rmsBQ, sys_seq_err_phred_per_ref_bp + MAX3(a_BQ_syserr_qual_2d, a_BQ_syserr_qual_fw, a_BQ_syserr_qual_rv)), a);
    
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
    
    const auto & p = symbol2CountCoverageSet12.seg_format_prep_sets.getByPos(refpos);
    fmt.APDP = {{ 
        p[SEG_a_DP], 
        p[SEG_a_NEAR_INS_DP], p[SEG_a_NEAR_DEL_DP], 
        p[SEG_a_NEAR_RTR_INS_DP], p[SEG_a_NEAR_RTR_DEL_DP], 
        p[SEG_aBQ2_DP],
        p[SEG_a_NEAR_INS_LEN], p[SEG_a_NEAR_DEL_LEN] }};
    
    fmt.APXM  = {{ p[SEG_a_XM], p[SEG_a_GO] }};
    fmt.APLRI = {{ p[SEG_a_LI], p[SEG_a_LIDP], p[SEG_a_RI], p[SEG_a_RIDP] }};
    fmt.APLRP = {{ p[SEG_a_L_DIST_SUM], p[SEG_a_R_DIST_SUM], p[SEG_a_INSLEN_SUM], p[SEG_a_DELLEN_SUM] }};
    
    const auto & t = symbol2CountCoverageSet12.seg_format_thres_sets.getByPos(refpos);
    fmt.AXMT =  {{ t[SEG_aXM2T], t[SEG_aXM2T] }};
    fmt.ALRIT = {{ t[SEG_aLI1T], t[SEG_aLI2T], t[SEG_aRI1T], t[SEG_aRI2T] }};
    fmt.ALRIt = {{ t[SEG_aLI1t], t[SEG_aLI2t], t[SEG_aRI1t], t[SEG_aRI2t] }};
    fmt.ALRPt = {{ t[SEG_aLP1t], t[SEG_aLP2t], t[SEG_aRP1t], t[SEG_aRP2t] }};
    fmt.ALRBt = {{ t[SEG_aLB1t], t[SEG_aLB2t], t[SEG_aRB1t], t[SEG_aRB2t] }};
    
    const auto & symbol_to_seg_format_depth_sets = symbol2CountCoverageSet12.symbol_to_seg_format_depth_sets;
    
    fill_symboltype_fmt(fmt.ADPff, fmt.aDPff, symbol_to_seg_format_depth_sets, SEG_aDPff, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ADPfr, fmt.aDPfr, symbol_to_seg_format_depth_sets, SEG_aDPfr, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ADPrf, fmt.aDPrf, symbol_to_seg_format_depth_sets, SEG_aDPrf, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ADPrr, fmt.aDPrr, symbol_to_seg_format_depth_sets, SEG_aDPrr, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.ALP1, fmt.aLP1, symbol_to_seg_format_depth_sets, SEG_aLP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALP2, fmt.aLP2, symbol_to_seg_format_depth_sets, SEG_aLP2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALPL, fmt.aLPL, symbol_to_seg_format_depth_sets, SEG_aLPL, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.ARP1, fmt.aRP1, symbol_to_seg_format_depth_sets, SEG_aRP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARP2, fmt.aRP2, symbol_to_seg_format_depth_sets, SEG_aRP2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARPL, fmt.aRPL, symbol_to_seg_format_depth_sets, SEG_aRPL, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.ALB1, fmt.aLB1, symbol_to_seg_format_depth_sets, SEG_aLB1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALB2, fmt.aLB2, symbol_to_seg_format_depth_sets, SEG_aLB2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALBL, fmt.aLBL, symbol_to_seg_format_depth_sets, SEG_aLBL, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.ARB1, fmt.aRB1, symbol_to_seg_format_depth_sets, SEG_aRB1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARB2, fmt.aRB2, symbol_to_seg_format_depth_sets, SEG_aRB2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARBL, fmt.aRBL, symbol_to_seg_format_depth_sets, SEG_aRBL, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.ABQ2, fmt.aBQ2, symbol_to_seg_format_depth_sets, SEG_aBQ2, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.APF1, fmt.aPF1, symbol_to_seg_format_depth_sets, SEG_aPF1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.APF2, fmt.aPF2, symbol_to_seg_format_depth_sets, SEG_aPF2, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.ALI1, fmt.aLI1, symbol_to_seg_format_depth_sets, SEG_aLI1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALI2, fmt.aLI2, symbol_to_seg_format_depth_sets, SEG_aLI2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALILf,fmt.aLILf,symbol_to_seg_format_depth_sets, SEG_aLILf,refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ALILr,fmt.aLILr,symbol_to_seg_format_depth_sets, SEG_aLILr,refpos, symbolType, refsymbol);

    fill_symboltype_fmt(fmt.ARI1, fmt.aRI1, symbol_to_seg_format_depth_sets, SEG_aRI1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARI2, fmt.aRI2, symbol_to_seg_format_depth_sets, SEG_aRI2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARILf,fmt.aRILf,symbol_to_seg_format_depth_sets, SEG_aRILf,refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.ARILr,fmt.aRILr,symbol_to_seg_format_depth_sets, SEG_aRILr,refpos, symbolType, refsymbol);
    
    const auto & symbol_to_frag_format_depth_sets = symbol2CountCoverageSet12.symbol_to_frag_format_depth_sets;
    
    fill_symboltype_fmt(fmt.BDPf,  fmt.bDPf,  symbol_to_frag_format_depth_sets[0], FRAG_bDP, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.BDPr,  fmt.bDPr,  symbol_to_frag_format_depth_sets[1], FRAG_bDP, refpos, symbolType, refsymbol);
    
    const auto & symbol_to_fam_format_depth_sets = symbol2CountCoverageSet12.symbol_to_fam_format_depth_sets_2strand;
    
    fill_symboltype_fmt(fmt.CDP0f, fmt.cDP0f, symbol_to_fam_format_depth_sets[0], FAM_cDP0, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP1f, fmt.cDP1f, symbol_to_fam_format_depth_sets[0], FAM_cDP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP2f, fmt.cDP2f, symbol_to_fam_format_depth_sets[0], FAM_cDP2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP3f, fmt.cDP3f, symbol_to_fam_format_depth_sets[0], FAM_cDP3, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.C1DPf, fmt.c1DPf, symbol_to_fam_format_depth_sets[0], FAM_c1DP, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPMf, fmt.cDPMf, symbol_to_fam_format_depth_sets[0], FAM_cDPM, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPmf, fmt.cDPmf, symbol_to_fam_format_depth_sets[0], FAM_cDPm, refpos, symbolType, refsymbol);
    
    fill_symboltype_fmt(fmt.CDP0r, fmt.cDP0r, symbol_to_fam_format_depth_sets[1], FAM_cDP0, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP1r, fmt.cDP1r, symbol_to_fam_format_depth_sets[1], FAM_cDP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP2r, fmt.cDP2r, symbol_to_fam_format_depth_sets[1], FAM_cDP2, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDP3r, fmt.cDP3r, symbol_to_fam_format_depth_sets[1], FAM_cDP3, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.C1DPr, fmt.c1DPr, symbol_to_fam_format_depth_sets[1], FAM_c1DP, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPMr, fmt.cDPMr, symbol_to_fam_format_depth_sets[1], FAM_cDPM, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.CDPmr, fmt.cDPmr, symbol_to_fam_format_depth_sets[1], FAM_cDPm, refpos, symbolType, refsymbol);

    const auto & symbol_to_duplex_format_depth_sets = symbol2CountCoverageSet12.symbol_to_duplex_format_depth_sets;
    fill_symboltype_fmt(fmt.DDP1,  fmt.dDP1,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP1, refpos, symbolType, refsymbol);
    fill_symboltype_fmt(fmt.DDP2,  fmt.dDP2,  symbol_to_duplex_format_depth_sets, DUPLEX_dDP2, refpos, symbolType, refsymbol);
    
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
        const unsigned int cDP0a,
        const std::string & gapSa,
        unsigned int minABQ,
        const unsigned int specialflag) {
    const unsigned int a = 0;
    
    const auto & symbol_to_seg_format_depth_sets = symbol2CountCoverageSet12.symbol_to_seg_format_depth_sets;
    
    fill_symbol_fmt(fmt.aDPff, symbol_to_seg_format_depth_sets, SEG_aDPff, refpos, symbol, a);
    fill_symbol_fmt(fmt.aDPfr, symbol_to_seg_format_depth_sets, SEG_aDPfr, refpos, symbol, a);
    fill_symbol_fmt(fmt.aDPrf, symbol_to_seg_format_depth_sets, SEG_aDPrf, refpos, symbol, a);
    fill_symbol_fmt(fmt.aDPrr, symbol_to_seg_format_depth_sets, SEG_aDPrr, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.aLP1, symbol_to_seg_format_depth_sets, SEG_aLP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLP2, symbol_to_seg_format_depth_sets, SEG_aLP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLPL, symbol_to_seg_format_depth_sets, SEG_aLPL, refpos, symbol, a);

    fill_symbol_fmt(fmt.aRP1, symbol_to_seg_format_depth_sets, SEG_aRP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRP2, symbol_to_seg_format_depth_sets, SEG_aRP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRPL, symbol_to_seg_format_depth_sets, SEG_aRPL, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.aLB1, symbol_to_seg_format_depth_sets, SEG_aLB1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLB2, symbol_to_seg_format_depth_sets, SEG_aLB2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLBL, symbol_to_seg_format_depth_sets, SEG_aLBL, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.aRB1, symbol_to_seg_format_depth_sets, SEG_aRB1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRB2, symbol_to_seg_format_depth_sets, SEG_aRB2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRBL, symbol_to_seg_format_depth_sets, SEG_aRBL, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.aXM1, symbol_to_seg_format_depth_sets, SEG_aXM1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aXM2, symbol_to_seg_format_depth_sets, SEG_aXM2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aBM2, symbol_to_seg_format_depth_sets, SEG_aBM2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aBQ2, symbol_to_seg_format_depth_sets, SEG_aBQ2, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.aPF1, symbol_to_seg_format_depth_sets, SEG_aPF1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aPF2, symbol_to_seg_format_depth_sets, SEG_aPF2, refpos, symbol, a);

    fill_symbol_fmt(fmt.aLI1, symbol_to_seg_format_depth_sets, SEG_aLI1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLI2, symbol_to_seg_format_depth_sets, SEG_aLI2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aLILf,symbol_to_seg_format_depth_sets, SEG_aLILf,refpos, symbol, a);
    fill_symbol_fmt(fmt.aLILr,symbol_to_seg_format_depth_sets, SEG_aLILr,refpos, symbol, a);
    
    fill_symbol_fmt(fmt.aRI1, symbol_to_seg_format_depth_sets, SEG_aRI1, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRI2, symbol_to_seg_format_depth_sets, SEG_aRI2, refpos, symbol, a);
    fill_symbol_fmt(fmt.aRILf,symbol_to_seg_format_depth_sets, SEG_aRILf,refpos, symbol, a);
    fill_symbol_fmt(fmt.aRILr,symbol_to_seg_format_depth_sets, SEG_aRILr,refpos, symbol, a);

    const auto & symbol_to_frag_format_depth_sets = symbol2CountCoverageSet12.symbol_to_frag_format_depth_sets;
    
    fill_symbol_fmt(fmt.bDPf, symbol_to_frag_format_depth_sets[0], FRAG_bDP, refpos, symbol, a);
    fill_symbol_fmt(fmt.bDPr, symbol_to_frag_format_depth_sets[1], FRAG_bDP, refpos, symbol, a);
    
    const auto & symbol_to_fam_format_depth_sets = symbol2CountCoverageSet12.symbol_to_fam_format_depth_sets_2strand;
    fill_symbol_fmt(fmt.cDP0f, symbol_to_fam_format_depth_sets[0], FAM_cDP0, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP1f, symbol_to_fam_format_depth_sets[0], FAM_cDP1, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP2f, symbol_to_fam_format_depth_sets[0], FAM_cDP2, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDP3f, symbol_to_fam_format_depth_sets[0], FAM_cDP3, refpos, symbol, a);
    fill_symbol_fmt(fmt.c1DPf, symbol_to_fam_format_depth_sets[0], FAM_c1DP, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPMf, symbol_to_fam_format_depth_sets[0], FAM_cDPM, refpos, symbol, a);
    fill_symbol_fmt(fmt.cDPmf, symbol_to_fam_format_depth_sets[0], FAM_cDPm, refpos, symbol, a);
    
    fill_symbol_fmt(fmt.cDP0r, symbol_to_fam_format_depth_sets[1], FAM_cDP0, refpos, symbol, a);
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
            minABQ,
            a);
    
    fmt.bHap = mutform2count4map_to_phase(mutform2count4vec_bq, indices_bq);
    fmt.cHap = mutform2count4map_to_phase(mutform2count4vec_fq, indices_fq);
    clear_push(fmt.bDPa, bDPa, a);
    clear_push(fmt.cDP0a, cDP0a, a);
    clear_push(fmt.gapSa, gapSa, a);
    return std::array<unsigned int, 2>({{fmt.bDPf[a] + fmt.bDPr[a], fmt.cDP0f[a] + fmt.cDP0r[a]}});
};

double calc_normFA_from_rawFA_refbias(double FA, double refbias) {
    return (FA + FA * refbias) / (FA + (1.0 - FA) / (1.0 + refbias) + FA * refbias);
}

int
BcfFormat_symbol_calc_DPv(
        bcfrec::BcfFormat & fmt,
        unsigned int refpos,
        const RegionalTandemRepeat & rtr1,
        const RegionalTandemRepeat & rtr2,
        const double tpfa,
        const bool is_amplicon,
        unsigned int specialflag) {
    
    const unsigned int a = 0;
    const bool is_rescued = (tpfa >= 0);
    const double pfa = (is_rescued ? tpfa : 0.5);
    
    // Non-UMI universality-based prequal allele fraction
    unsigned int ADP = (fmt.ADPff[0] + fmt.ADPfr[0] + fmt.ADPrf[0] + fmt.ADPrr[0]);
    unsigned int aDP = (fmt.aDPff[a] + fmt.aDPfr[a] + fmt.aDPrf[a] + fmt.aDPrr[a]);
    double aDPFA = (aDP + pfa) / (double)(ADP + 1.0);
    
    const auto & f = fmt;
    const auto symbol = AlignmentSymbol(LAST(f.VTI));
    
    // substitution in indel region, substitution in indel-prone region or indel, other cases 
    unsigned int aDPplus = (isSymbolSubstitution(symbol) ? 0 : ((aDP + 1) / 2));
    double dp_coef = ((symbol == LINK_M) ?  MAX(0.02, 1.0 - MAX(rtr1.tracklen, rtr2.tracklen) / (MAX3(1, f.ALPL[0], f.ARPL[0]) / MAX(1.0/150.0, f.ABQ2[0]))) : 1.0);
    double _aPprior = 1024*9/2 / (is_rescued ? phred2prob(8) : 1); // probability that there is one indel at a pos
    const bool is_in_indel_read = (   (f.APXM[1]) / 15.0 * 2  * (1.25+1e-7) > (aDP + aDPplus) * dp_coef);
    const bool is_in_indel_len  = (MAX(f.APDP[1],  f.APDP[2]) * (1.25+1e-7) > (aDP + aDPplus) * dp_coef);
    const bool is_in_indel_rtr  = (MAX(f.APDP[3],  f.APDP[4]) * (1.25+1e-7) > (aDP + aDPplus) * dp_coef);
    const bool is_in_rtr = (MAX(rtr1.tracklen, rtr2.tracklen) > 8);
    if (is_in_indel_read)      { _aPprior /= (64+32); }
    if (is_in_indel_len)       { _aPprior /= 24; }
    else if (is_in_indel_rtr)  { _aPprior /=  8; }
    else if (is_in_rtr)        { _aPprior /=  3; }
    const double aPprior = _aPprior;
    fmt.note += std::string("//aPprior/") +std::to_string(aPprior) + "//";
    const double aIprior = log(isSymbolSubstitution(symbol) ? 1e5 : 3e2) + 1;
    const double aSBprior = log(isSymbolSubstitution(symbol) ? (pow(10.0, fmt.aBQ[a] / 10.0) + 1) : 3e2);
    
    const auto aLPFAx2 = dp4_to_pcFA<false>(f.aLP1[a], aDP, f.ALP2[0] + f.aLP1[a] - f.aLP2[a], ADP, 3.0, log(aPprior),
            MAX(1, f.aLPL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ALPL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? 0.10 : 0.5)); 
            // use f.aBQ2[a] instead of aDP ???
    const auto aRPFAx2 = dp4_to_pcFA<false>(f.aRP1[a], aDP, f.ARP2[0] + f.aRP1[a] - f.aRP2[a], ADP, 3.0, log(aPprior),
            MAX(1, f.aRPL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ARPL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? 0.10 : 0.5));
    double aLPFA = aLPFAx2[0];
    double aRPFA = aRPFAx2[0];
    
    const auto aLBFAx2 = dp4_to_pcFA<false>(f.aLB1[a], aDP, f.ALB2[0] + f.aLB1[a] - f.aLB2[a], ADP, 3.0, log(aPprior),
            MAX(1, f.aLBL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ALBL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? 0.10 : 0.5));
    const auto aRBFAx2 = dp4_to_pcFA<false>(f.aRB1[a], aDP, f.ARB2[0] + f.aRB1[a] - f.aRB2[a], ADP, 3.0, log(aPprior),
            MAX(1, f.aRBL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ARBL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? 0.10 : 0.5));
    double aLBFA = aLBFAx2[0];
    double aRBFA = aRBFAx2[0];
    fmt.note += std::string("symbol/") + SYMBOL_TO_DESC_ARR[symbol] + "/";
    fmt.note += std::string("aLBFAdetail//") + other_join2(std::vector<double> {{ 
            (double)f.aRB1[a], (double)aDP, (double)(f.ARB2[0] + f.aRB1[a] - f.aRB2[a]), (double)ADP, 3.0, log(aPprior),
            MAX(1, f.aRBL[a]) / (double)MAX(1, f.aBQ2[a]), MAX(1, f.ARBL[0]) / (double)MAX(1, f.ABQ2[0]), ((is_in_indel_read) ? 0.10 : 0.5) }} , "/") + "//";

    std::array<double, 2> _aLIFAx2 = {{ 0, 0 }};
    {
        double ALpd = (f.ALI2[0] +                 0.5) / (f.ADPfr[0] + f.ADPrr[0] - f.ALI2[0] + 0.5);
        double aLpd = (f.aLI1[a] + ALpd / (1.0 + ALpd)) / (f.aDPfr[a] + f.aDPrr[a] - f.aLI1[a] + 1.0 / (1.0 + ALpd)); 
        _aLIFAx2 = dp4_to_pcFA<false>(
                (f.aLI1[a]), (f.aDPfr[a] + f.aDPrr[a]), 
                (f.ALI2[0] + f.aLI1[a] - f.aLI2[a]), (f.ADPfr[0] + f.ADPrr[0]), 
                3.0, log(aIprior), 
                aLpd, ALpd, 0.25, 0.5);
        fmt.note += std::string("ALpd/aLpd/") + std::to_string(ALpd) + "/" + std::to_string(aLpd) + "//";
    }
    double aLIFA = _aLIFAx2[0] * (is_amplicon ? 1.0 : MAX(1.0, aDPFA / _aLIFAx2[1]));
    
    std::array<double, 2> _aRIFAx2 = {{ 0, 0 }};
    {
        double ARpd = (f.ARI2[0] +                 0.5) / (f.ADPff[0] + f.ADPrf[0] - f.ARI2[0] + 0.5);
        double aRpd = (f.aRI1[a] + ARpd / (1.0 + ARpd)) / (f.aDPff[a] + f.aDPrf[a] - f.aRI1[a] + 1.0 / (1.0 + ARpd));
        _aRIFAx2 = dp4_to_pcFA<false>(
                (f.aRI1[a]), (f.aDPff[a] + f.aDPrf[a]), 
                (f.ARI2[0] + f.aRI1[a] - f.aRI2[a]), (f.ADPff[0] + f.ADPrf[0]), 
                3.0, log(aIprior),
                aRpd, ARpd, 0.25, 0.5);
        fmt.note += std::string("ARpd/aRpd/") + std::to_string(ARpd) + "/" + std::to_string(aRpd) + "//";
    }
    double aRIFA = _aRIFAx2[0] * (is_amplicon ? 1.0 : MAX(1.0, aDPFA / _aRIFAx2[1]));
    
    if (isSymbolIns(symbol) || isSymbolDel(symbol)) {
        const auto & indelstring = LAST(fmt.gapSa);
        const bool is_in_indel_major_reg = ((MAX(f.APDP[1], f.APDP[3]) + MAX(f.APDP[2], f.APDP[4])) * 0.5 * (1.0+(double)FLT_EPSILON) < aDP);
        if ((indelstring.size() * aDPFA >= 2) || 
                (MAX(rtr1.tracklen, rtr2.tracklen) >= 16
                && is_in_indel_major_reg)) {
            aLPFA = 2.0;
            aRPFA = 2.0;
            aLBFA = 2.0;
            aRBFA = 2.0;
        } else {
            fmt.note += std::string("rtr-tracklen/") + std::to_string(rtr1.tracklen) + "/" + std::to_string(rtr2.tracklen) + "/";
        }
        aLIFA = aRIFA = 2.0;
    } else if (LINK_M == symbol || LINK_NN == symbol) {
        aLBFA *= MAX(1, fmt.ALB1[0]) / (double)MAX(ADP, 1);
        aRBFA *= MAX(1, fmt.ARB1[0]) / (double)MAX(ADP, 1);
    } 
    
    const double aPFFA = (fmt.aPF1[a] + pfa * (is_rescued ? 100.0 : (double)(fmt.aXM2[a] / MAX(1, aDP)))) / (fmt.APF2[0] + (fmt.aPF1[a] - fmt.aPF2[a]) + 1.0*100);
    const auto aSSFAx2 = dp4_to_pcFA<true>(fmt.aDPff[a] + fmt.aDPrf[a], fmt.aDPfr[a] + fmt.aDPrr[a], fmt.ADPff[0] + fmt.ADPrf[0], fmt.ADPfr[0] + fmt.ADPrr[0], 3.0, aSBprior);
    const auto cROFA0x2 = dp4_to_pcFA<true>(fmt.cDP0f[a], fmt.cDP0r[a], fmt.CDP0f[0], fmt.CDP0r[0], 3.0, log(mathsquare(aDPFA) * 1e6 + 2));
    const auto cROFA2x2 = dp4_to_pcFA<true>(fmt.cDP2f[a], fmt.cDP2r[a], fmt.CDP2f[0], fmt.CDP2r[0], 3.0, log(mathsquare(aDPFA) * 1e6 + 2));
    const double aSSFA = aSSFAx2[0];
    double cROFA0 = cROFA0x2[0];
    double cROFA2 = cROFA2x2[0];
    if ((fmt.cDP1r[a] * fmt.CDP1f[0] < fmt.cDP1f[a] * fmt.CDP1r[0]) 
            && (MIN(fmt.aLILf[a], fmt.aRILf[a]) / MAX(fmt.aDPff[a] + fmt.aDPfr[a], 1) > 25)
            && (MIN(fmt.ALILr[0], fmt.ARILr[0]) / MAX(fmt.ADPrf[0] + fmt.ADPrr[0], 1) < 25)) { // fw allele frac > rv allele frac
        cROFA0 = cROFA2 = 2.0;
    }
    if ((fmt.cDP1f[a] * fmt.CDP1r[0] < fmt.cDP1r[a] * fmt.CDP1f[0]) 
            && (MIN(fmt.aLILr[a], fmt.aRILr[a]) / MAX(fmt.aDPrr[a] + fmt.aDPrr[a], 1) > 25)
            && (MIN(fmt.ALILf[0], fmt.ARILf[0]) / MAX(fmt.ADPff[0] + fmt.ADPfr[0], 1) < 25)) { // fw allele frac > rv allele frac
        cROFA0 = cROFA2 = 2.0;
    }
    double bFA = (fmt.bDPa[a] + pfa) / (fmt.BDPf[0] + fmt.BDPr[0] + 1.0);
    double cFA0 = (fmt.cDP0a[a] + pfa) / (fmt.CDP0f[0] + fmt.CDP0r[0] + 1.0);

    auto min_aFA_vec = std::vector<double>{{
            aDPFA,
            aLPFA, aRPFA,
            aLBFA, aRBFA,
            aLIFA, aRIFA, 
            aSSFA, aPFFA * aSSFAx2[0] / MAX(aSSFAx2[0], aSSFAx2[1]) }};
    
    double min_aFA = MINVEC(min_aFA_vec);    
    fmt.note += std::string("raw/pos4/ss/xm//") + other_join(min_aFA_vec, "/") + "//";
    
    double min_bcFA = MIN3(bFA, cFA0, cROFA0);
    fmt.note += std::string("b/c/o1/o2//") + other_join(std::vector<double> {{ bFA, cFA0, cROFA0, cROFA2 }}, "/") + "//";
    
    double cFA2 = (fmt.cDP2f[a] + fmt.cDP2r[a] + pfa) / (fmt.CDP2f[0] + fmt.CDP2r[0] + 1.0);
    double cFA3 = (fmt.cDP3f[a] + fmt.cDP3r[a] + pfa) / (fmt.CDP3f[0] + fmt.CDP3r[0] + 1.0);
    assert( cFA2 > 0);
    assert( cFA2 < 1);
    assert( cFA3 > 0);
    assert( cFA3 < 1);
    
    auto min_cFA23_vec = std::vector<double> {{ cFA2, cFA3 }};
    fmt.note += std::string("2/3//") + other_join(min_cFA23_vec, "/") + "//";
    double umi_cFA = MINVEC(min_cFA23_vec);
    
    double frac_umi2seg = MIN(1.0, umi_cFA / aDPFA * BETWEEN(2.0 * mathsquare(min_aFA / aDPFA), 0.5, 1.0));
    fmt.note += std::string("frac_umi2seg//") + std::to_string(frac_umi2seg) + "//";
    double frac_c1_b = BETWEEN(2.5 * mathsquare(cFA0 / bFA), 0.4, 1.0);
    
    double refbias = 0;
    if ((isSymbolIns(symbol) || isSymbolDel(symbol)) && is_rescued) {
        const auto & indelstring = LAST(fmt.gapSa);
        const unsigned int indel_noinfo_nbases = (indelstring.size() * (isSymbolIns(symbol) ? 2 : 1) + MAX3(indelstring.size(), rtr1.tracklen, rtr2.tracklen));
        refbias = (double)(indel_noinfo_nbases) / ((double)(MIN(f.ALPL[0], f.ARPL[0]) * 2 + indel_noinfo_nbases) / (double)(f.ABQ2[0] + 0.5));
        refbias = MIN(refbias, 2.0);
    }
    // non-UMI push
    clear_push(fmt.cDP1v, (int)(calc_normFA_from_rawFA_refbias(MIN(min_aFA, min_bcFA), refbias) * (fmt.CDP0f[0] +fmt.CDP0r[0]) * 100), a);
    double min_abcFA_w = MINVEC(std::vector<double>{{
            aLPFA, aRPFA,
            aLBFA, aRBFA,
            bFA }});
    clear_push(fmt.cDP1w, (int)(calc_normFA_from_rawFA_refbias(min_abcFA_w, refbias) * (fmt.CDP0f[0] +fmt.CDP0r[0]) * 100), a);
    double min_abcFA_x = MINVEC(std::vector<double>{{ aPFFA, cFA0 }});
    clear_push(fmt.cDP1x, 1+(int)(min_abcFA_x * (fmt.CDP0f[0] +fmt.CDP0r[0]) * 100), a);
    
    // UMI push
    clear_push(fmt.cDP2v, (int)(calc_normFA_from_rawFA_refbias(MIN(min_aFA * frac_umi2seg, cROFA2) * frac_c1_b, refbias) * (fmt.CDP2f[0] +fmt.CDP2r[0]) * 100), a);
    double umi_cFA_w = MINVEC(std::vector<double>{{
            aLPFA * frac_umi2seg, aRPFA * frac_umi2seg,
            aLBFA * frac_umi2seg, aRBFA * frac_umi2seg,
            aDPFA * frac_umi2seg,
            cROFA2
            }});
    clear_push(fmt.cDP2w, (int)(calc_normFA_from_rawFA_refbias(umi_cFA_w, refbias) * (fmt.CDP2f[0] +fmt.CDP2r[0]) * 100), a);
    double umi_cFA_x = MINVEC(std::vector<double>{{ 
            aPFFA * frac_umi2seg,
            aDPFA * frac_umi2seg,
            cROFA2
            }});
    clear_push(fmt.cDP2x, 1+(int)(MIN(umi_cFA_x * (fmt.CDP2f[0] +fmt.CDP2r[0]), fmt.cDP2f[a] +fmt.cDP2r[a]) * 100), a);
    
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
BcfFormat_symbol_sum_DPv(auto & fmts, const unsigned int a = 0) {
    std::array<unsigned int, NUM_REDUCTIONS> cDPxx1s = {{0}};
    std::array<unsigned int, NUM_REDUCTIONS> cDPxx2s = {{0}};
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
        
        const unsigned int var_bdepth,
        const unsigned int var_cdepth,
        
        const unsigned int ins_bdepth,
        const unsigned int ins_cdepth,
        const unsigned int del_bdepth,
        const unsigned int del_cdepth,
        
        const unsigned int ins1_bdepth,
        const unsigned int ins1_cdepth,
        const unsigned int del1_bdepth,
        const unsigned int del1_cdepth,
        
        const std::string & repeatunit,
        const unsigned int repeatnum,
        double powlaw_anyvar_base, // is approx (double)90,
        double powlaw_exponent, // is approx (double)3,
        double powlaw_sscs_phrederr,
        double powlaw_sscs_inc, // is approx 30,
        const unsigned int phred_varcall_err_per_map_err_per_base, // is approx 20,
        const bool is_rescued,
        const auto & rtr1,
        const auto & rtr2,
        int tid,
        int refpos,
        AlignmentSymbol refsymbol,
        double tpfa,
        const unsigned int specialflag) {
    
    const unsigned int a = 0;
    const auto symbol = AlignmentSymbol(LAST(fmt.VTI));
    
    const unsigned int aDP = (fmt.aDPff[a] + fmt.aDPfr[a] + fmt.aDPrf[a] + fmt.aDPrr[a]);
    const unsigned int ADP = (fmt.ADPff[0] + fmt.ADPrf[0] + fmt.ADPfr[0] + fmt.ADPrr[0]);
    const unsigned int bDP = (fmt.bDPf[a] + fmt.bDPr[a]);
    const unsigned int BDP = (fmt.BDPf[0] + fmt.BDPr[0]);
    
    double prior_weight;
    prior_weight = 1.0 / (fmt.cDPmf[a] + 1.0);
    double cMmQf = 10/log(10) * log((fmt.cDPMf[a] + fmt.cDPmf[a] + pow(10, 25/10.0) * prior_weight) 
            / (fmt.cDPmf[a] + prior_weight));
    prior_weight = 1.0 / (fmt.cDPmr[a] + 1.0);
    double cMmQr = 10/log(10) * log((fmt.cDPMr[a] + fmt.cDPmr[a] + pow(10, 25/10.0) * prior_weight) 
            / (fmt.cDPmr[a] + prior_weight));
    
    int cMmQ = sqrt((mathsquare(cMmQf) * fmt.cDP3f[a] + mathsquare(cMmQr) * fmt.cDP3r[a]) / (fmt.cDP3f[a] + fmt.cDP3r[a]));
    
    const int bIADnormcnt = (fmt.bIADb[a] * 100 + 1);
    const unsigned int bIAQinc = ((isSymbolIns(symbol) || isSymbolDel(symbol)) ? 
            (unsigned int)round(10/log(10)*log(MIN(5, (double)(BDP + var_bdepth + 0.5) / (double)(BDP - var_bdepth + 0.5)))): 0);
    int duped_frag_binom_qual = (fmt.bIAQb[a] + bIAQinc) * MIN(bIADnormcnt, fmt.cDP1v[a] + 1) / bIADnormcnt;
    if (fmt.cDP1v[a] > bIADnormcnt) {
        duped_frag_binom_qual += (int)(10/log(10) * log((double)fmt.cDP1v[a] / (double)bIADnormcnt) * fmt.bIADb[a]);
    }
    
    const int sscs_dec = (int)(mathsquare((double)(fmt.CDP1f[0] + fmt.CDP1r[0] + 300) / (double)(fmt.BDPf[0] + fmt.BDPr[0] + 300)) * 60);
    
    const int cIADnormcnt = (fmt.cIADf[a] * 100 + 1);
    int sscs_binom_qual_fw = fmt.cIAQf[a] + fmt.cIAQr[a] * MIN(60 - fmt.cIDQf[a], fmt.cIDQr[a]) / MAX(fmt.cIDQr[a], 1);
    int sscs_binom_qual_rv = fmt.cIAQr[a] + fmt.cIAQf[a] * MIN(60 - fmt.cIDQr[a], fmt.cIDQf[a]) / MAX(fmt.cIDQf[a], 1);
    int dedup_sscs_binom_qual = MAX(sscs_binom_qual_fw, sscs_binom_qual_rv) * MIN(cIADnormcnt, fmt.cDP2v[a] + 1) / cIADnormcnt - sscs_dec;
    
    const int tn_reward = 10;
    double min_bcFA_v = (((double)(fmt.cDP1v[a]) + 0.5) / (double)(fmt.CDP1f[0] * 100 + fmt.CDP1r[0] * 100 + 1.0));
    int dedup_frag_powlaw_qual_v = (int)(powlaw_exponent * 10.0 / log(10.0) * log(min_bcFA_v) + (powlaw_anyvar_base + 7));
    
    double min_bcFA_w = ((double)(fmt.cDP1w[a] + 0.5) / (double)(fmt.CDP1f[0] * 100 + fmt.CDP1r[0] * 100 + 1.0));
    int dedup_frag_powlaw_qual_w = (int)(powlaw_exponent * 10.0 / log(10.0) * log(min_bcFA_w) + (powlaw_anyvar_base + 7 + tn_reward));
    
    int ds_vq_inc_powlaw = (int)(10/log(10)*MIN(log((fmt.cDP1f[a] + 0.5) / (fmt.CDP1f[0] + 1.0)), log((fmt.cDP1r[a] + 0.5) / (fmt.CDP1r[0] +1.0)))) 
            + (powlaw_sscs_phrederr);
    int ds_vq_inc_binom = 3 * MIN(fmt.cDP2f[a], fmt.cDP2r[a]);
    
    powlaw_sscs_inc += MAX(0, MIN5(sscs_binom_qual_fw, sscs_binom_qual_rv, ds_vq_inc_powlaw, ds_vq_inc_binom, 4)) + MIN(0, (int)cMmQ - 25);
    double umi_cFA =  (((double)(fmt.cDP2v[a]) + 0.5) / ((double)(fmt.CDP2f[0] * 100 + fmt.CDP2r[0] * 100) + 1.0 + 300));
    int sscs_base_2 = powlaw_anyvar_base + 7 + 7 + powlaw_sscs_inc - sscs_dec;
    int sscs_powlaw_qual_v = (int)((powlaw_exponent * 10.0 / log(10.0) * log(umi_cFA)    + sscs_base_2));
    int sscs_powlaw_qual_w = (int)((powlaw_exponent * 10.0 / log(10.0) * log(min_bcFA_w) + sscs_base_2));
    
    fmt.note += std::string("min_bcFA_v(") + std::to_string(min_bcFA_v) + ")cDP1v(" + std::to_string(fmt.cDP1v[a]) + ")" + "CDP1v(" + std::to_string(fmt.CDP1v[0]) + ")";
    fmt.note += std::string("umi_cFA_v(") + std::to_string(umi_cFA) + ")cDP2v(" + std::to_string(fmt.cDP2v[a]) + ")" + "CDP2v(" + std::to_string(fmt.CDP2v[0]) + ")" + "cMmQ(" + std::to_string(cMmQ) + ")";
    
    int indel_penal4multialleles = 0;
    int indel_penal4multialleles_g = 0;
    int indel_penal4multialleles_soma = 0;
    const std::string & indelstring = fmt.gapSa[a];
    
    if (indelstring.size() > 0 && fmt.cDP0a[a] > 0) {
        const double indel_pq = (double)MIN((int)indel_phred(8.0, indelstring.size(), repeatunit.size(), repeatnum), (int)24) + 2 - (double)10;
        const auto eff_tracklen1 = (repeatunit.size() * MAX(1, repeatnum) - repeatunit.size());
        const auto eff_tracklen2 = (MAX(rtr1.tracklen - rtr1.unitlen, rtr2.tracklen - rtr2.unitlen) / 3);
        const double indel_ic = 10.0 / log(10.0) * log((double)MAX(indelstring.size() + (isSymbolIns(symbol) ? INS_N_ANCHOR_BASES : 0), 1U) / (double)(MAX(eff_tracklen1, eff_tracklen2) + 1));
        auto indelcdepth = (isSymbolIns(symbol) ? ins_cdepth : del_cdepth);
        if (LINK_D1 == symbol) {
            indelcdepth += ins1_cdepth;
        }
        if (LINK_I1 == symbol) {
            indelcdepth += del1_cdepth / 4.0;
        }
        const unsigned int nearInDelDP = (isSymbolIns(symbol) ? fmt.APDP[1] : fmt.APDP[2]);
        
        assert (nearInDelDP >= aDP || !fprintf(stderr, "nearInDelDP >= aDP failed ( %d >= %d failed) at tid %d pos %d!\n", nearInDelDP, aDP, tid, refpos));
        
        const auto indel_penal4multialleles1 = (int)round(11.0/log(2.0) * log((double)(indelcdepth + 1e-3) / (double)(fmt.cDP0a[a] + 1e-3)));
        const auto indel_penal4multialleles2 = (int)round( 8.0/log(2.0) * log((double)(nearInDelDP + 1e-3) / (double)(aDP + 1e-3)));
        indel_penal4multialleles_g = (int)round(2 * 8.0/log(2.0) * log((double)(ins_cdepth + del_cdepth + 1e-3) / (double)(fmt.cDP0a[a] + 1e-3))) - 22;
        
        if (isSymbolIns(symbol)) {
            indel_penal4multialleles = (indel_penal4multialleles1 * 16 / (int)(16 + indelstring.size())); 
            indel_penal4multialleles_soma = (indel_penal4multialleles1 * 16 / (int)(16 + indelstring.size()));
        } else {
            indel_penal4multialleles = MAX(indel_penal4multialleles1, indel_penal4multialleles2);
            indel_penal4multialleles_soma = indel_penal4multialleles1;
        }
        fmt.note += std::string("/indelpq/") + std::to_string(indel_pq) + "/" + indelstring + "/" + repeatunit + "/" + std::to_string(repeatnum) 
                + std::string("/indelic/") + std::to_string(indel_ic)
                + std::string("/indelcdepth/") + std::to_string(indelcdepth) 
                + std::string("/indelpm/") + std::to_string(indel_penal4multialleles) + "/" + std::to_string(indel_penal4multialleles_g)
                ;
        dedup_frag_powlaw_qual_v += (int)(indel_ic);
        dedup_frag_powlaw_qual_w += (int)(indel_ic);
        duped_frag_binom_qual  += (int)(indel_pq);
    }
    
    clear_push(fmt.bIAQ, duped_frag_binom_qual, a);
    clear_push(fmt.cIAQ, dedup_sscs_binom_qual, a);
    
    clear_push(fmt.cPCQ1, dedup_frag_powlaw_qual_w, a);
    clear_push(fmt.cPLQ1, dedup_frag_powlaw_qual_v, a);
    
    clear_push(fmt.cPCQ2, sscs_powlaw_qual_w, a);
    clear_push(fmt.cPLQ2, sscs_powlaw_qual_v, a);
    
    const auto bMQinc = (int)MAX(phred_varcall_err_per_map_err_per_base, 31 + 2 + (int)(2 * 10.0/log(10.0) * log((bDP + 0.5) / (double)(BDP + 1.0))));
    const int syserr_q = MIN(
            (int)((fmt.bMQ[a]) + MIN(bMQinc, bDP * 3)),
            (isSymbolSubstitution(AlignmentSymbol(LAST(fmt.VTI))) ? (fmt.aBQQ[a]) : (1000*1000)));
    
    const int PENAL4LOWDEP = 37;
    const int aDPpc = ((refsymbol == symbol) ? 1 : 0);
    int penal4BQerr = (isSymbolSubstitution(symbol) ? (5 + (int)(PENAL4LOWDEP / (int)mathsquare((int64_t)MAX(1, aDP + aDPpc)))) : 0);
    
    int indel_q_inc = ((isSymbolSubstitution(symbol) || is_rescued) ? 0 : indel_len_rusize_phred(indelstring.size(), repeatnum));
    clear_push(fmt.gVQ1, MAX(
            (int)0, 
            indel_q_inc + MIN3(
                syserr_q,
                (int)LAST(fmt.bIAQ) - penal4BQerr, 
                (int)LAST(fmt.cPLQ1)) - 2 * MAX3(
                    (int)0, 
                    (int)indel_penal4multialleles - (int)11, 
                    (int)indel_penal4multialleles_g)),
    a);
    
    clear_push(fmt.cVQ1, MAX(
            (int)0, 
            MIN3(
                syserr_q,
                (int)LAST(fmt.bIAQ) - penal4BQerr, // Does this correspond to in-silico mixing artifact ???
                (int)LAST(fmt.cPLQ1)) - indel_penal4multialleles_soma),
    a);
    
    clear_push(fmt.cVQ2, MIN3(syserr_q,
            (int)LAST(fmt.cIAQ),
            (int)LAST(fmt.cPLQ2)) - indel_penal4multialleles, 
    a);
    
    // TODO; check if reducing all allele read count to increase alt allele frac in case of ref bias makes more sense
    double base_contamfrac = 0.02; // ContEst https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3167057/ Fig. 1 TCGA Ovarian 
    double contamfrac = base_contamfrac + (1 - base_contamfrac) * (tpfa > 0 ? tpfa : 0) * 0.05;
    double errfrac = contamfrac + (1.0 - contamfrac) * (refsymbol == symbol ? phred2prob(MIN(rtr1.indelphred, rtr2.indelphred)) : 0);
    auto binom_contam_LODQ = (int)calc_binom_10log10_likeratio(errfrac, fmt.cDP1v[a], fmt.CDP1v[0]);
    auto power_contam_LODQ = (int)(10.0/log(10.00) * powlaw_exponent * MAX(logit2(min_bcFA_v, errfrac), 0.0));
    clear_push(fmt.CONTQ, MIN(binom_contam_LODQ, power_contam_LODQ), a);
};

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

unsigned int hetLODQ(double allele1count, double allele2count, double expfrac, double powlaw_exponent = 3.0) {
    auto binomLODQ = (int)calc_binom_10log10_likeratio(expfrac, allele1count, allele2count);
    auto powerLODQ = (int)(10.0/log(10.00) * powlaw_exponent * MAX(logit2((allele1count + 0.5) * 0.5 / expfrac, (allele2count + 0.5) * 0.5 / (1.0 - expfrac)), 0.0));
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
    // return (double)(fmtp->cDP1f[fmtp->cDP1f.size()-1] + fmtp->cDP1r[fmtp->cDP1r.size()-1]) / 100.0;
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
        unsigned int refpos,
        unsigned int extended_inclu_beg_pos,
        unsigned int central_readlen,
        unsigned int outputflag,
        const bool should_output_all_germline,
        const bool is_rescued,
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
    int ref_alodq = INT_MIN;
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
    
    int a0LODQ = ALODQ(ref_alt1_alt2_alt3[0].second);
    int a1LODQ = ALODQ(ref_alt1_alt2_alt3[1].second);
    int a2LODQ = ALODQ(ref_alt1_alt2_alt3[2].second);
    int a3LODQ = ALODQ(ref_alt1_alt2_alt3[3].second);
    
    auto fmtptr0 = ref_alt1_alt2_alt3[0].second;
    auto fmtptr1 = ref_alt1_alt2_alt3[1].second;
    auto fmtptr2 = ref_alt1_alt2_alt3[2].second;
    const bool isSubst = isSymbolSubstitution(refsymbol);
    const AlignmentSymbol symbolNN = ((isSubst || !is_rescued) ? BASE_NN : LINK_NN);
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
    
    int a0a1LODQ = hetLODQ(ad0norm, ad1norm, 1.0 - 0.47);
    int a1a0LODQ = hetLODQ(ad1norm, ad0norm, 0.47);
    int a1a2LODQ = hetLODQ(ad1norm, ad2norm, 0.5);
    int a2a1LODQ = hetLODQ(ad2norm, ad1norm, 0.5);
    
    std::array<std::string, 4> GTidx2GT {{
        "0/0",
        "0/1",
        "1/1",
        "1/2"
    }};
    
    int phred_homref = 0;
    int phred_hetero = (isSubst ? 31     : 41-1);
    int phred_homalt = (isSubst ? 33     : 43-1);
    int phred_tri_al = (isSubst ? (54+5) : 50-1); 
    // https://www.genetics.org/content/184/1/233 and https://www.fsigenetics.com/article/S1872-4973(20)30003-X/fulltext : triallelic-SNP-phred = 29*2-3+5
    
    if (is_rescued) {
        UPDATE_MIN(a0LODQ, ref_alt1_alt2_alt3[0].second->CONTQ[1]);
        UPDATE_MIN(a1LODQ, ref_alt1_alt2_alt3[1].second->CONTQ[1]);
        UPDATE_MIN(a2LODQ, ref_alt1_alt2_alt3[2].second->CONTQ[1]);
        UPDATE_MIN(a3LODQ, ref_alt1_alt2_alt3[3].second->CONTQ[1]);
    } else {
        UPDATE_MIN(a0LODQ, ref_alt1_alt2_alt3[0].second->CONTQ[1]);
    }
    // Two principles: 
    //  1. Every unexpected signal given a genotype reduces the likelihood of this genotype
    //  2. Signals are correlated with each other with a phred-scaled correlation of phred_hetero
    const auto a2penal = MAX(a2LODQ - phred_hetero, 0);
    const auto a3penal = MAX(a3LODQ - phred_hetero, 0);
    const auto a01hetp = MAX(MAX(a0a1LODQ, a1a0LODQ) - (0-0), 0);
    const auto a12hetp = MAX(MAX(a1a2LODQ, a2a1LODQ) - (0-0), 0);
    const auto a03trip = MAX(a0LODQ, a3LODQ);
    
    std::array<std::pair<int, int>, 4> GL4raw = {{
        std::make_pair(0,     (-phred_homref - a1LODQ                - a2penal - a3penal)),
        std::make_pair(1,     (-phred_hetero - MAX(a01hetp, a2LODQ)  - MAX(MIN(a01hetp, a2LODQ)  - phred_hetero, 0) - a3penal)),
        std::make_pair(2,     (-phred_homalt - MAX( a0LODQ, a2LODQ)  - MAX(MIN( a0LODQ, a2LODQ)  - phred_hetero, 0) - a3penal)),
        std::make_pair(3,     (-phred_tri_al - MAX(a12hetp, a03trip) - MAX(MIN(a12hetp, a03trip) - phred_hetero, 0) 
            - MAX(MIN(a12hetp, MIN(a0LODQ, a3LODQ)) - phred_hetero, 0)))
    }};
    auto ret = GL4raw[0].second - MAX3(GL4raw[1].second, GL4raw[2].second, GL4raw[3].second);
    if (0 == ((OUTVAR_GERMLINE) & outputflag)) {
        return std::make_tuple(ret, fmtptr1, fmtptr2);
    }
    
    auto GL4 = GL4raw;
    std::sort(GL4.rbegin(), GL4.rend(), PairSecondLess);
    
    size_t GLidx = GL4[0].first;
    if (0 == GLidx && (!should_output_all_germline) && MAX(LAST(fmtptr1->cDP0a), LAST(fmtptr2->cDP0a)) <= 2) {
        return std::make_tuple(ret, fmtptr1, fmtptr2);
    }
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
            tname, refpos, s1, false);
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
                std::vector<std::pair<std::array<unsigned int, 2>, std::string>> indelstrings2 = indel_get_majority(*(ref_alt1_alt2_alt3[2].second), 
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
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[0].second->cDP0a, 1, 0));
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[1].second->cDP0a, 1, 0));
    if (3 == GLidx) {
        germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[2].second->cDP0a, 1, 0));
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
        std::string("GT:GQ:HQ:FT:CDP1:cDP1:GL4:GST:note"),
        string_join(std::array<std::string, 9>
        {{
            germ_GT, 
            std::to_string(germ_GQ), 
            std::string("0,0"), 
            "PASS", // string_join(germ_FT, "/"),
            other_join(std::array<int, 2> {{
                (int)(symbol_format_vec[0].second->CDP1f[0] + symbol_format_vec[0].second->CDP1r[0]), 
                (int)(symbol_format_vec[0].second->CDP1f[1] + symbol_format_vec[0].second->CDP1r[1]) }},
                ","
            ),
            other_join(germ_ADR, std::string(",")), 
            other_join(std::array<int, 4>
            {{ 
                GL4raw[0].second, GL4raw[1].second, GL4raw[2].second, GL4raw[3].second, 
            }}, ","),
            other_join(std::array<int, 8>
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
generate_vcf_header(const char *ref_fasta_fname, 
        const char *platform, 
        unsigned int central_readlen,
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
    ret += "##INFO=<ID=TNBQF,Number=4,Type=Float,Description=\"Binomial reward, power-law reward, systematic-error penalty, and normal-adjusted tumor variant quality computed using deduplicated read fragments.\">\n";
    ret += "##INFO=<ID=TNBQN,Number=4,Type=Float,Description=\"TNBQF using non-bias-adjusted deduplicated read fragments.\">\n";
    ret += "##INFO=<ID=TNBQP,Number=1,Type=Float,Description=\"Bias-induced penalty to deduplicated read fragments after comparing tumor bias with normal bias.\">\n";
    ret += "##INFO=<ID=TNCQF,Number=4,Type=Float,Description=\"Binomial reward, power-law reward, systematic-error penalty, and normal-adjusted tumor variant quality computed using consensus families of read fragments.\">\n";
    ret += "##INFO=<ID=RU,Number=1,Type=String,Description=\"The shortest repeating unit in the reference\">\n";
    ret += "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"The number of non-interrupted RUs in the reference\">\n";
    ret += "##INFO=<ID=R3X2,Number=6,Type=Integer,Description=\"Repeat start position, repeat track length, and repeat unit size at the two positions before and after this VCF position.\">\n"; 
    
    for (unsigned int i = 0; i < bcfrec::FORMAT_NUM; i++) {
        ret += std::string("") + bcfrec::FORMAT_LINES[i] + "\n";
    }
    
    ret += std::string("") + "##FORMAT=<ID=CDP1,Number=2,Type=Integer,Description=\"Total deduped depth and deletion deduped depth\">\n";
    ret += std::string("") + "##FORMAT=<ID=cDP1,Number=R,Type=Integer,Description=\"The deduped depth of each allele\">\n";
    ret += std::string("") + "##FORMAT=<ID=GL4,Number=4,Type=Integer,Description=\"The four genotype likelihoods for 0/0, 0/1, 1/1, and 1/2\">\n";
    ret += std::string("") + "##FORMAT=<ID=GST,Number=.,Type=Integer,Description=\"The genotype statistics\">\n";

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
fill_tki(auto & tki, const auto & fmt, unsigned int a = 1) {
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
};

const auto
calc_binom_powlaw_syserr_normv_quals(
        auto tAD, auto tDP, int tVQ, int tnVQcap,
        auto nAD, auto nDP, int nVQ, const double penal_dimret_coef = 15) {
    int binom_b10log10like = (int)calc_binom_10log10_likeratio((tDP - tAD) / (tDP), nDP - nAD, nAD);
    double bjpfrac = ((tAD + 0.5) / (tDP + 1.0)) / ((nAD + 0.5) / (nDP + 1.0));
    int powlaw_b10log10like = (int)(3 * 10 / log(10) * log(bjpfrac));
    int prior_phred = ((penal_dimret_coef < 25) ? 8 : 1);
    int tnVQinc = MAX3(-prior_phred, -nAD*3, MIN(binom_b10log10like - prior_phred, powlaw_b10log10like - prior_phred));
    int tnVQdec = MAX(0, nVQ - MAX(0, MIN(binom_b10log10like - prior_phred, (int)(mathsquare(log(MAX(bjpfrac, 1.01)) / log(2)) * penal_dimret_coef))));
    
    int tnVQ = tVQ + tnVQinc; 
    if (penal_dimret_coef < 1e6) {
        tnVQ -= tnVQdec;
    }
    tnVQ = MIN(tnVQcap, tnVQ);
    return std::array<int, 4> {{binom_b10log10like, powlaw_b10log10like, tnVQdec, tnVQ }};
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
        const bcfrec::BcfFormat & fmt,
        auto & tki,
        const bool is_tumor_format_retrieved,
        const int nlodq,
        const double vcfqual_thres,
        const unsigned int all_allele_depth_thres,
        const unsigned int var_allele_depth_thres,
        const bool should_output_ref_allele,
        const bcf_hdr_t *g_bcf_hdr,
        const auto & baq_offsetarr,
        unsigned int specialflag) {
    
    const std::string & indelstring = LAST(fmt.gapSa);
    
    const bcfrec::BcfFormat & nfm = (tki.ref_alt.size() > 0 ? fmt : FORMAT_UNCOV);
    if (tki.ref_alt.size() == 0) {
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
    
    assert (2 == fmt.cDP1v.size());
    assert (2 == fmt.cDP2v.size());
    assert (2 == fmt.CDP1v.size());
    assert (2 == fmt.CDP2v.size());
    assert (2 == fmt.cVQ1.size());
    assert (2 == fmt.cVQ2.size());
    
    assert (tki.cDP1x <= tki.CDP1x || !fprintf(stderr, "%d <= %d failed for cDP1x tname %s pos %d aln-symbol %d\n", tki.cDP1x, tki.CDP1x, tname, refpos, symbol));
    assert (tki.cDP2x <= tki.CDP2x || !fprintf(stderr, "%d <= %d failed for cDP2x tname %s pos %d aln-symbol %d\n", tki.cDP2x, tki.CDP2x, tname, refpos, symbol));
    const auto & rtr1 =  region_repeatvec.at(MAX(refpos - region_offset, 3) - 3);
    const auto & rtr2 =  region_repeatvec.at(MIN(refpos - region_offset + 3, region_repeatvec.size() - 3));
    const auto rtr1_tpos = ((0 == rtr1.tracklen) ? 0 : (region_offset + rtr1.begpos));
    const auto rtr2_tpos = ((0 == rtr2.tracklen) ? 0 : (region_offset + rtr2.begpos));

    double nfm_cDP1 = collectget(nfm.cDP1f, 1) + collectget(nfm.cDP1r, 1);
    double nfm_CDP1 = collectget(nfm.CDP1f, 0) + collectget(nfm.CDP1r, 0);
    const auto nfm_cDP1x = collectget(nfm.cDP1x, 1);
    const auto nfm_CDP1x = collectget(nfm.CDP1x, 0);

    const auto b_binom_powlaw_syserr_normv_q4filter = calc_binom_powlaw_syserr_normv_quals(
            (tki.cDP1x + 0.5) / 100.0 + 0.0, 
            (tki.CDP1x + 1.0) / 100.0 + 0.0, 
            (tki.cVQ1),
            (tki.cPCQ1),
            (nfm_cDP1x + 0.5) / 100.0 + 0.0, 
            (nfm_CDP1x + 1.0) / 100.0 + 0.0, 
            non_neg_minus(collectget(nfm.cVQ1, 1), (indelstring.size() == 0 ? 0 : 40)),
            15.0);
    
    const auto c_binom_powlaw_syserr_normv_q4 = calc_binom_powlaw_syserr_normv_quals(
            (tki.cDP2x + 0.5) / 100.0 + 0.0, 
            (tki.CDP2x + 1.0) / 100.0 + 0.0, 
            (tki.cVQ2),
            (tki.cPCQ1) + 30,
            (collectget(nfm.cDP2x, 1) + 0.0) / 100.0 + 0.5, 
            (collectget(nfm.CDP2x, 0) + 0.0) / 100.0 + 1.0, 
            (collectget(nfm.cVQ2, 1)), 
            (15.0));
    
    int tlodq = MAX(b_binom_powlaw_syserr_normv_q4filter[3], c_binom_powlaw_syserr_normv_q4[3]);
    
    int somaticq = MIN(tlodq, nlodq);
    int vcfqual = ((tki.ref_alt.size() > 0) ? somaticq : tlodq);
    std::string infostring = std::string(tki.ref_alt.size() > 0 ? "SOMATIC" : "ANY_VAR");
    infostring += std::string(";SomaticQ=") + std::to_string(somaticq);
    infostring += std::string(";TLODQ=") + std::to_string(tlodq);
    infostring += std::string(";NLODQ=") + std::to_string(nlodq);
    infostring += std::string(";TNBQF=") + other_join(b_binom_powlaw_syserr_normv_q4filter);
    infostring += std::string(";TNCQF=") + other_join(c_binom_powlaw_syserr_normv_q4);
    infostring += std::string(";RU=") + repeatunit + ";RC=" + std::to_string(repeatnum);
    // This can be useful for debugging only.
    // if (isSymbolSubstitution(symbol)) { infostring += std::string(";RBAQ=") + std::to_string(baq_offsetarr.getByPos(refpos)); }
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
                infostring, bcfrec::FORMAT_STR_PER_REC}}, "\t") + "\t";
        bcfrec::streamAppendBcfFormat(out_string, fmt);
        out_string += ((tki.ref_alt.size() > 0 && is_tumor_format_retrieved) ? bcf1_to_string(g_bcf_hdr, tki.bcf1_record) : std::string("")) + "\n";
    }
    return 0;
};

#endif
