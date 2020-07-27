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

#define THE_BAQ_MAX 200
#define INDEL_MUL_PER_BAQ 2.0
#define SNV_MUL_PER_BAQ 3.0
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
posToIndelToCount_updateByConsensus_old(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src, uint32_t incvalue = 1) {
    for (auto src_pos2indel2count4it : src) {
        auto src_pos = src_pos2indel2count4it.first;
        auto src_indel2count = src_pos2indel2count4it.second;
        if (src_indel2count.size() == 1) {
            assert (src_indel2count.begin()->second >= incvalue
                    || !fprintf(stderr, "Condition %d >= %d failed!\n", src_indel2count.begin()->second, incvalue)
                    || !(std::cerr << "value is " << src_indel2count.begin()->first << std::endl));
            posToIndelToCount_inc<T>(dst, src_pos, src_indel2count.begin()->first, incvalue);
        } else if (src_indel2count.size() > 1) {
            // auto insresult2 = insresult->insert(std::make_pair(T(), 0)); insresult2-second++; // do nothing
        } else { assert(false); }
    }
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

#define NM_MULT_NORM_COEF (100.0 * (1.0 + DBL_EPSILON))

#define N_MODELS 1 // previously used 23 & 11 models
#define QUAL_PRE_ADD 1
// #define bam_phredi(b, i) (phred2bucket(bam_get_qual((b))[(i)]))
#define bam_phredi(b, i) (bam_get_qual((b))[(i)])

const bcfrec::BcfFormat FORMAT_UNCOV = bcfrec::BcfFormat();

enum ValueType {
    SYMBOL_COUNT_SUM,
    BASE_QUALITY_MAX,
    _NUM_VALUE_TYPES,
};

// base and link between bases.
// partition many insertion sequences into very few categories 
// to prevent overfitting and 
// to prevent false negative due to many erroneous inserted sequences.
enum AlignmentSymbol {
    BASE_A,  //   = 0,
    BASE_C,  //   = 1,
    BASE_G,  //   = 2,
    BASE_T,  //   = 3,
    BASE_N,  //   = 4, // ambigous in the original sequencing data
    BASE_NN, //   = 5, // ambigous base after collapsing different reads
    //BASE_P = 6; // padded in deleted sequence
    LINK_M,  //   = 6, // absence of any gap
    LINK_D3P,// = 7, // deletion of length 3 or plus
    LINK_D2, //  = 8,  // deletion of length 2
    LINK_D1, //  = 9,
    LINK_I3P,//  = 10, // insertion of length 1 // where the inserted sequence is not a repeat
    LINK_I2, //  = 11, 
    LINK_I1, // = 12, 
    LINK_NN, //  = 13, // ambiguous link between bases
    // LINK_P = 13; // padded in deleted sequence
    END_ALIGNMENT_SYMBOLS,
};

#define NUM_INS_SYMBOLS 3
const std::array<AlignmentSymbol, NUM_INS_SYMBOLS> INS_SYMBOLS = {{LINK_I1, LINK_I2, LINK_I3P}};

#define NUM_DEL_SYMBOLS 3
const std::array<AlignmentSymbol, NUM_DEL_SYMBOLS> DEL_SYMBOLS = {{LINK_D1, LINK_D2, LINK_D3P}};

const std::array<AlignmentSymbol, (NUM_INS_SYMBOLS+NUM_DEL_SYMBOLS)> INDEL_SYMBOLS = {{LINK_I1, LINK_I2, LINK_I3P, LINK_D1, LINK_D2, LINK_D3P}};

bool 
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

const AlignmentSymbol 
insLenToSymbol(unsigned int len) {
    //assert(len > 0);
    return 1 == len ? LINK_I1 : (2 == len ? LINK_I2 : LINK_I3P);
}

const AlignmentSymbol 
delLenToSymbol(unsigned int len) {
    //assert(len > 0);
    return 1 == len ? LINK_D1 : (2 == len ? LINK_D2 : LINK_D3P);
}

#define NUM_ALIGNMENT_SYMBOLS 14
static_assert(NUM_ALIGNMENT_SYMBOLS == END_ALIGNMENT_SYMBOLS);

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
    [BASE_NN] = "<BN>", 
    [LINK_M] = "<LR>", 
    [LINK_D3P] = "<LD3P>", [LINK_D2] = "<LD2>", [LINK_D1] = "<LD1>",
    [LINK_I3P] = "<LI3P>", [LINK_I2] = "<LI2>", [LINK_I1] = "<LI1>",
    [LINK_NN] = "<LN>",
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

typedef \
uint32_t \
molcount_t;

// edge distance bucket
#define NUM_EDBUCKS (SIGN2UNSIGN(11+1)) // (11*12/2+1)
#define NUM_NMBUCKS (SIGN2UNSIGN(12))
// #define EDBUCK_SIZE 4

const std::array<unsigned int, (NUM_EDBUCKS*(NUM_EDBUCKS-1)/2+1+1)> 
DIST_TO_EDBUCK = {{0, 
    0, 
    1,  1,                                  // 3
    2,  2,  2, 
    3,  3,  3,  3,                          // 10 
    4,  4,  4,  4,  4, 
    5,  5,  5,  5,  5,  5,                  // 21
    6,  6,  6,  6,  6,  6,  6,    
    7,  7,  7,  7,  7,  7,  7,  7,          // 36
    8,  8,  8,  8,  8,  8,  8,  8,  8,      // 45
    9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  // 55
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, // 66
    11
}};
const std::array<unsigned int, NUM_EDBUCKS>
EDBUCK_TO_DIST = {{
    1*2/2,
    2*3/2,
    3*4/2,
    4*5/2,
    5*6/2,
    6*7/2,
    7*8/2,
    8*9/2,
    9*10/2,
    10*11/2,
    11*12/2,
    12*13/2
}};

unsigned int 
pos2edbuck(unsigned int pos) {
    //return MIN(pos / EDBUCK_SIZE, NUM_EDBUCKS - 1);
    return DIST_TO_EDBUCK[MIN(pos, DIST_TO_EDBUCK.size()-1)];
}

unsigned int 
edbuck2pos(unsigned int edbuck) {
    assert(EDBUCK_TO_DIST.size() > edbuck);
    // return edbuck * EDBUCK_SIZE;
    return EDBUCK_TO_DIST[edbuck];
}

typedef std::array<molcount_t, NUM_BUCKETS> Bucket2Count;
typedef std::array<molcount_t, NUM_EDBUCKS> Bucket2CountEdgeDist;
typedef std::array<molcount_t, NUM_NMBUCKS> Bucket2CountNumMisma;

const int 
_print_Bucket2CountEdgeDist(const Bucket2CountEdgeDist & arg) {
    for (size_t i = 0; i < NUM_EDBUCKS; i++) {
        LOG(logINFO) << arg.at(i) << "\t";
    }
    return 0;
}

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

typedef GenericSymbol2Bucket2Count<Bucket2Count> Symbol2Bucket2Count;
typedef GenericSymbol2Bucket2Count<Bucket2CountEdgeDist> Symbol2Bucket2CountEdgeDist;
typedef GenericSymbol2Bucket2Count<Bucket2CountNumMisma> Symbol2Bucket2CountNumMisma;

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
    
    template<ValueType T_SymbolCountType, bool TIndelIsMajor>
    const AlignmentSymbol
    _updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count,
            const SymbolType symbolType, const AlignmentSymbol ambig_pos, unsigned int incvalue2) {
        AlignmentSymbol argmax_count = END_ALIGNMENT_SYMBOLS; // AlignmentSymbol(0) is not fully correct
        unsigned int max_count = 0;
        unsigned int sum_count = 0;
        thatSymbol2Count.template fillConsensusCounts<TIndelIsMajor>(argmax_count, max_count, sum_count, symbolType);
        unsigned int incvalue;
        if (T_SymbolCountType == SYMBOL_COUNT_SUM) {
            incvalue = incvalue2;
        } else if (T_SymbolCountType == BASE_QUALITY_MAX) {
            assert(max_count < 96);
            incvalue = THE_PHRED_TO_ERROR_PROBABILITY.over2pow16[max_count];
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

    template<ValueType T_SymbolCountType, bool TIndelIsMajor = false>
    std::array<AlignmentSymbol, 2>
    updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count, unsigned int incvalue = 1) {
        AlignmentSymbol baseSymb = this->template _updateByConsensus<T_SymbolCountType, false        >(thatSymbol2Count, BASE_SYMBOL, BASE_NN, incvalue);
        AlignmentSymbol linkSymb = this->template _updateByConsensus<T_SymbolCountType, TIndelIsMajor>(thatSymbol2Count, LINK_SYMBOL, LINK_NN, incvalue);
        return {baseSymb, linkSymb};
    };
    
    template <bool TIsIncVariable = false>
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
    updateByFiltering(std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> & con_symbols, const GenericSymbol2Count<TInteger> & other, const GenericSymbol2Count<TInteger> & thres, 
            unsigned int incvalue, unsigned int TPos, unsigned int TStrand) {
        int ret = 0;
        AlignmentSymbol consalpha;
        unsigned int countalpha, totalalpha;
        for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1 + (unsigned int)symbolType)) {
            // consalpha = END_ALIGNMENT_SYMBOLS;
            if (LINK_SYMBOL == symbolType) {
                other.template fillConsensusCounts<true >(consalpha, countalpha, totalalpha, symbolType);
            } else {
                other.template fillConsensusCounts<false>(consalpha, countalpha, totalalpha, symbolType);
            }
            auto adjcount = MAX(countalpha * 2, totalalpha) - totalalpha;
            if (adjcount >= (thres.getSymbolCount(consalpha)) && adjcount > 0) {
                this->symbol2data[consalpha] += incvalue;
                ret++;
                //LOG(logINFO) << " The value " << countalpha << " can indeed pass the base quality threshold " << thres.getSymbolCount(consalpha) << " for symbol " << SYMBOL_TO_DESC_ARR[consalpha];
            } else {
                //LOG(logINFO) << " The value " << countalpha << " cannot pass the base quality threshold " << thres.getSymbolCount(consalpha) << " for symbol " << SYMBOL_TO_DESC_ARR[consalpha];
            }
            con_symbols[symbolType] = consalpha;
        }
        return ret;
    };
};

typedef GenericSymbol2Count<uint32_t> Symbol2Count;
typedef GenericSymbol2Count<uint64_t> Symbol2CountUint64;

template<class T>
class CoveredRegion {
         
protected:
    std::vector<T> idx2symbol2data;
    std::array<std::map<uint32_t, std::map<uint32_t   , uint32_t>>, NUM_INS_SYMBOLS> pos2dlen2data;
    std::array<std::map<uint32_t, std::map<std::string, uint32_t>>, NUM_DEL_SYMBOLS> pos2iseq2data;
    
    const size_t 
    _extern2intern4pos(size_t extern_ref_pos) {
        assert(extern_ref_pos >= incluBegPosition);
        assert(extern_ref_pos <  incluBegPosition + idx2symbol2data.size() || !fprintf(stderr, "%d is not within (%d - %d)", extern_ref_pos, incluBegPosition, incluBegPosition + idx2symbol2data.size()));
        return extern_ref_pos - incluBegPosition;
    };
public:
    
    const uint32_t tid;
    const uint32_t incluBegPosition; // end_pos = incluBegPosition + idx2symbol2data

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
    
    template <LinkType TLinkType>
    const auto &
    getPosToIndelToData(const AlignmentSymbol s) const {
        static_assert(INS_LINK == TLinkType || DEL_LINK == TLinkType);
        if (TLinkType == INS_LINK) {
            return getRefPosToIseqToData(s);
        } else {
            return getRefPosToDlenToData(s);
        }
    };
};

template <class TSymbol2Bucket2Count>
class GenericSymbol2Bucket2CountCoverage : public CoveredRegion<TSymbol2Bucket2Count> {
    public:
    GenericSymbol2Bucket2CountCoverage() : CoveredRegion<TSymbol2Bucket2Count>(0, 0, 1) { /*assert(false);*/ } // just to fix compiling error
    GenericSymbol2Bucket2CountCoverage(auto tid, auto beg, auto end) : CoveredRegion<TSymbol2Bucket2Count>(tid, beg, end) {}
};

typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2Count> Symbol2Bucket2CountCoverage;
typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2CountEdgeDist> Symbol2Bucket2CountCoverageEdgeDist;
typedef GenericSymbol2Bucket2CountCoverage<Symbol2Bucket2CountNumMisma> Symbol2Bucket2CountCoverageNumMisma;

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
    // const unsigned int rank_STR[6+1] = {0, 65, 32, 33, 34, 35, 36}; // monomer is ranked first, followed by 2xdimer, hexamer, pentamer, ..., dimer, and zero
    // const int rank_STR[6+1] = {0, 201, 196, 197, 198, 199, 200};
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
    /*
    if (0 == (cigar_oplen % repeatsize_at_max_repeatnum)) {
        indel_n_units = cigar_oplen / repeatsize_at_max_repeatnum;
    } else {
        indel_n_units = max_repeatnum;
    }
    */
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
        for (unsigned int i = refpos; i != MIN(repeat_endpos, refstring.size()); i++) {
            if (tl > region_repeatvec[i].tracklen) {
                region_repeatvec[i].begpos = refpos;
                region_repeatvec[i].tracklen = tl;
                region_repeatvec[i].unitlen = repeatsize_at_max_repeatnum;
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

int
dealwith_seg_regbias(unsigned int qpos1, unsigned int qpos2, const bam1_t *b, unsigned int seq_pos_bias_nbases, 
        auto & bq_dirs_count, unsigned int strand, unsigned int isrc, 
        auto & bq_regs_count, unsigned int rpos, AlignmentSymbol symbol, auto & bq_nms_count, unsigned int nm,
        auto & bq_baq_sum, const auto & qr_baq_vec, 
        auto & bq_dirs_bqsum,
        unsigned int specialflag = 0) {
    bq_dirs_count[strand*2+isrc].template inc<SYMBOL_COUNT_SUM>(rpos, symbol, 1, b);
    
    const bool is_lbias = (qpos1 + 1 <= seq_pos_bias_nbases);
    const bool is_rbias = (b->core.l_qseq - qpos2 <= seq_pos_bias_nbases);
    int idx = -1;
    if (is_lbias && is_rbias) {
        idx = 0;
    } else if (is_lbias) {
        idx = 1;
    } else if (is_rbias) {
        idx = 2;
    }
    if (idx >= 0) {
        bq_regs_count[idx].template inc<SYMBOL_COUNT_SUM>(rpos, symbol, 1, b);
    }
    bq_nms_count[strand*2+isrc].template inc<SYMBOL_COUNT_SUM>(rpos, symbol, nm, b);
    auto baq = MIN(THE_BAQ_MAX, qr_baq_vec.at(rpos - b->core.pos));
    bq_baq_sum[0].template inc<SYMBOL_COUNT_SUM>(rpos, symbol, mathsquare(baq) / THE_BAQ_MAX, b);
    auto qualphred = (bam_phredi(b, qpos1) + bam_phredi(b, qpos2)) / 2;
    bq_dirs_bqsum[strand*2+isrc].template inc<SYMBOL_COUNT_SUM>(rpos, symbol, qualphred, b);
    return idx;
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
    
    std::vector<unsigned int> 
    computeZeroBasedPosToInsLenVec(unsigned int & totInsLen) {
        std::vector<unsigned int> ret(this->getExcluEndPosition() - this->getIncluBegPosition(), 0);
        for (AlignmentSymbol s : INS_SYMBOLS) {
            for (auto & refPosToIseqToData : this->getRefPosToIseqToData(s)) {
                uint32_t refPos = refPosToIseqToData.first;
                std::map<std::string, uint32_t> iseqToData = refPosToIseqToData.second;
                uint32_t insCount = 0;
                uint32_t insSumSize = 0; 
                for (auto iseqToDataPair : iseqToData) {
                    std::string iseq = iseqToDataPair.first;
                    uint32_t count = iseqToDataPair.second;
                    insCount += count;
                    insSumSize = iseq.size() * count;
                }
                if (insCount > this->getByPos(refPos).getSymbolCount(LINK_M)) {
                    ret[refPos-this->getIncluBegPosition()] = insSumSize / insCount;
                    totInsLen += insSumSize / insCount;
                }
            }
        }
        return ret;
    }
    
    // mainly for merging reads in one family
    template <bool TIsIncVariable = false>
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
    template<ValueType T_ConsensusType, bool TIndelIsMajor = false> 
    void
    updateByConsensus(const GenericSymbol2CountCoverage<TSymbol2Count> &other, unsigned int incvalue = 1,
            const bool update_pos2indel2count = true, const bool update_idx2symbol2data = true) {
        this->assertUpdateIsLegal(other);
        
        if (update_idx2symbol2data) {
            for (size_t epos = other.getIncluBegPosition(); epos < other.getExcluEndPosition(); epos++) {
                const std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols = this->getRefByPos(epos).template updateByConsensus<T_ConsensusType, TIndelIsMajor>(other.getByPos(epos), incvalue);
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
    updateByFiltering(const GenericSymbol2CountCoverage<TSymbol2Count> &other, const GenericSymbol2CountCoverage<TSymbol2Count> &thres, unsigned int incvalue = 1,
            const bool update_pos2indel2count = true, unsigned int TStrand=0) {
        this->assertUpdateIsLegal(other);
        int num_updated_pos = 0;
        // other.assertUpdateIsLegal(thres); // may not hold because threshold is delimited by bed whereas this is not delimited by bed.
        auto incluBegPos = MAX(other.getIncluBegPosition(), thres.getIncluBegPosition()); 
        auto excluEndPos = MIN(other.getExcluEndPosition(), thres.getExcluEndPosition());
        std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols; 
        for (size_t epos = incluBegPos; epos < excluEndPos; epos++) {
            int updateresult = this->getRefByPos(epos).updateByFiltering(consymbols, other.getByPos(epos), thres.getByPos(epos), incvalue, epos, TStrand);
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
        // return (int)excluEndPos - (int)incluBegPos;
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
    compute_qr_baq_vec(const bam1_t *b, const auto & region_repeatvec, const unsigned int region_offset, unsigned int baq_per_aligned_base) {
        const unsigned int rbeg = b->core.pos;
        const unsigned int rend = bam_endpos(b);
        const uint8_t *bseq = bam_get_seq(b);
        const uint32_t n_cigar = b->core.n_cigar;
        const uint32_t *cigar = bam_get_cigar(b);
        /*
        unsigned int baq_init = baq_per_aligned_base;
        for (unsigned int i = 0; i != n_cigar; i++) {
            const uint32_t c = cigar[i];
            const unsigned int cigar_op = bam_cigar_op(c);
            if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) {
                baq_init = 0;
            }
        }
        */
        std::vector<uint16_t> refBAQvec(rend-rbeg);
        {
            unsigned int qpos = 0;
            unsigned int rpos = b->core.pos;
            unsigned int baq = 0;
            unsigned int prev_repeat_begpos = 0; // UINT_MAX;
            unsigned int prev_rpos = rbeg;
            // unsigned int baq_per_aligned_repeatregion = baq_per_aligned_base;
            for (int i = 0; i != n_cigar; i++) {
                const uint32_t c = cigar[i];
                const unsigned int cigar_op = bam_cigar_op(c);
                const unsigned int cigar_oplen = bam_cigar_oplen(c);
                if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                    for (unsigned int i2 = 0; i2 != cigar_oplen; i2++) {
                        unsigned int repeat_begpos = region_repeatvec.at(rpos - region_offset).begpos;
                        assert(repeat_begpos >= prev_repeat_begpos);
                        if (prev_repeat_begpos != repeat_begpos) {
                            /*
                            unsigned int rrsize = rpos - prev_rpos;
                            if (1 < rrsize) {
                                baq += MIN(baq_per_aligned_base * rrsize, baq_per_aligned_repeatregion);
                                baq_per_aligned_repeatregion += 1;
                            } else {
                                baq += baq_per_aligned_base;
                            }
                            prev_rpos = rpos;
                            */
                            baq += baq_per_aligned_base;
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
            unsigned int prev_repeat_begpos = UINT_MAX;
            // unsigned int prev_rpos = rend;
            // unsigned int baq_per_aligned_repeatregion = baq_per_aligned_base;
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
                            /*
                            unsigned int rrsize = prev_rpos - rpos;
                            if (1 < rrsize) {
                                baq += MIN(baq_per_aligned_base * rrsize, baq_per_aligned_repeatregion);
                                baq_per_aligned_repeatregion += 1;
                            } else {
                                baq += baq_per_aligned_base;
                            }
                            prev_rpos = rpos;
                            */
                            baq += baq_per_aligned_base;
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
#if 0
        if (b->core.pos % 10000 == 0) {
            LOG(logINFO) << "For readname " << bam_get_qname(b) << " at tid " << b->core.tid << " position " <<  b->core.pos;
            std::string refBAQvecStr = "";
            for (auto a : refBAQvec) {
                refBAQvecStr += std::to_string(a) + ",";
            }
            LOG(logINFO) << refBAQvecStr;
        }
#endif
        return refBAQvec;
    }
    
    unsigned int
    proton_cigarlen2phred(unsigned int cigarlen) {
        unsigned int oplen2cigar[8] = {0, 7, 13, 18, 22, 25, 27};
        return oplen2cigar[MIN(cigarlen, 6)];
    }
    
    template<ValueType TUpdateType, bool TIsProton, bool THasDups, bool TFillSeqDir, unsigned int TIndelAddPhred = 0*29>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByAln(const bam1_t *const b, unsigned int frag_indel_ext, 
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhredArg, 
            unsigned int frag_indel_basemax, 
            unsigned int nogap_phred, // this is obsolete and replace by frag_indel_basemax
            const auto & region_symbolvec, const unsigned int region_offset,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 4> & bq_dirs_count,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 3> & bq_regs_count,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 4> & bq_nms_count,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 1> & bq_baq_sum,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 4> & bq_dirs_bqsum,
            unsigned int sidereg_nbases,
            const std::vector<RegionalTandemRepeat> & region_repeatvec,
            unsigned int baq_per_aligned_base,
            uint32_t sepcialflag = 0) {
        static_assert(BASE_QUALITY_MAX == TUpdateType || SYMBOL_COUNT_SUM == TUpdateType);
        assert(this->tid == SIGN2UNSIGN(b->core.tid));
        assert(this->getIncluBegPosition() <= SIGN2UNSIGN(b->core.pos)   || !fprintf(stderr, "%lu <= %d failed", this->getIncluBegPosition(), b->core.pos));
        assert(this->getExcluEndPosition() >= SIGN2UNSIGN(bam_endpos(b)) || !fprintf(stderr, "%lu >= %d failed", this->getExcluEndPosition(), bam_endpos(b)));
        const auto symbolType2addPhred = symbolType2addPhredArg; // std::array({0, 0});

        const bool isrc = ((b->core.flag & 0x10) == 0x10);
        const bool isr2 = ((b->core.flag & 0x80) == 0x80 && (b->core.flag & 0x1) == 0x1);
        const unsigned int strand = (isrc ^ isr2);
        unsigned int qpos = 0;
        unsigned int rpos = b->core.pos;
        const uint32_t n_cigar = b->core.n_cigar;
        const uint32_t *cigar =  bam_get_cigar(b);
        const uint8_t *bseq = bam_get_seq(b);
        unsigned int incvalue = 1;
        
        unsigned int nm_var = 0;
        std::vector<uint16_t> qr_baq_vec;
        if (TFillSeqDir) {
            uint8_t * bam_aux_data = bam_aux_get(b, "NM");
            if (bam_aux_data != NULL) {
                nm_var = bam_aux2i(bam_aux_data);
            } else {
                for (unsigned int i = 0; i < n_cigar; i++) {
                    int32_t c = cigar[i];
                    unsigned int cigar_op = bam_cigar_op(c);
                    unsigned int cigar_oplen = bam_cigar_oplen(c);
                    if (cigar_op == BAM_CDIFF) { nm_var++ ; }
                    else if (BAM_CINS == cigar_op || BAM_CDEL == cigar_op) { nm_var += cigar_oplen; }
                }
            }
            nm_var = (unsigned int)floor(NM_MULT_NORM_COEF * sqrt((double)nm_var / (double)MAX(30, b->core.l_qseq) + DBL_EPSILON));
            qr_baq_vec = compute_qr_baq_vec(b, region_repeatvec, region_offset, baq_per_aligned_base);
        }
        const unsigned int nm = nm_var;
        
        for (unsigned int i = 0; i < n_cigar; i++) {
            uint32_t c = cigar[i];
            unsigned int cigar_op = bam_cigar_op(c);
            unsigned int cigar_oplen = bam_cigar_oplen(c);
            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                    assert((rpos >= SIGN2UNSIGN(b->core.pos) && rpos < SIGN2UNSIGN(bam_endpos(b)))
                            || !fprintf(stderr, "Bam line with QNAME %s has rpos that is not within the range (%d - %d)", bam_get_qname(b), b->core.pos, bam_endpos(b)));
                    if (i2 > 0) {
                        const bool T_update_ref_nogap = TFillSeqDir;
                        if (TUpdateType == BASE_QUALITY_MAX) {
                            const auto noindel_phred = (T_update_ref_nogap ? (MIN(qr_baq_vec[rpos - b->core.pos] / 2, 17)) : 17); // frag_indel_basemax
                            incvalue = ((TIsProton || T_update_ref_nogap) ? (MIN3(noindel_phred, bam_phredi(b, qpos-1), bam_phredi(b, qpos))) : noindel_phred); 
                            // + symbolType2addPhred[LINK_SYMBOL];
                        }
                        this->template inc<TUpdateType>(rpos, LINK_M, incvalue, b);
                        if (TFillSeqDir) {
                            dealwith_seg_regbias(qpos, qpos, b, sidereg_nbases, bq_dirs_count, strand, isrc, bq_regs_count, rpos, LINK_M, bq_nms_count, nm, 
                                    bq_baq_sum, qr_baq_vec, bq_dirs_bqsum, 0);
                        }
                    }
                    unsigned int base4bit = bam_seqi(bseq, qpos);
                    unsigned int base3bit = seq_nt16_int[base4bit];
                    if (TUpdateType == BASE_QUALITY_MAX) {
                        incvalue = bam_phredi(b, qpos) 
                                + (QUAL_PRE_ADD ? symbolType2addPhred[BASE_SYMBOL] : 0);
                    }
                    this->template inc<TUpdateType>(rpos, AlignmentSymbol(base3bit), incvalue, b);
                    if (TFillSeqDir) {
                        dealwith_seg_regbias(qpos, qpos, b, sidereg_nbases, bq_dirs_count, strand, isrc, bq_regs_count, rpos, AlignmentSymbol(base3bit), bq_nms_count, nm,
                                bq_baq_sum, qr_baq_vec, bq_dirs_bqsum, 0);
                    }
                    rpos += 1;
                    qpos += 1;
                }
            } else if (cigar_op == BAM_CINS) {
                const bool is_ins_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(b->core.l_qseq));
                unsigned int inslen = SIGN2UNSIGN(cigar_oplen);
                if (TUpdateType == BASE_QUALITY_MAX) {
                    if (TIndelAddPhred) {
                        auto addidq = (THasDups ? 0 : (MIN(cigar_oplen - 1, SIGN2UNSIGN(3 - 1)) * frag_indel_ext));
                        incvalue = TIndelAddPhred + addidq;
                    } else if (is_ins_at_read_end) {
                        LOG(logWARNING) << "Query " << bam_get_qname(b) << " has insertion of legnth " << cigar_oplen << " at " << qpos
                                << " which is not exclusively between 0 and " << b->core.l_qseq << " aligned to tid " << b->core.tid << " and position " << rpos;
                        incvalue = (0 != qpos ? bam_phredi(b, qpos-1) : 
                                ((qpos + cigar_oplen < SIGN2UNSIGN(b->core.l_qseq)) ? 
                                bam_phredi(b, qpos + SIGN2UNSIGN(cigar_oplen)) : 1)) 
                                + (QUAL_PRE_ADD ? symbolType2addPhred[LINK_SYMBOL] : 0); // + addidq; // 
                    } else {
                        unsigned int phredvalue = ref_to_phredvalue(inslen, region_symbolvec, rpos - region_offset,
                                frag_indel_basemax, 8.0, cigar_oplen, cigar_op, b);
                        incvalue = MIN(MIN(bam_phredi(b, qpos-1), bam_phredi(b, qpos + cigar_oplen)) + (TIsProton ? proton_cigarlen2phred(cigar_oplen) : 0), phredvalue)
                                + (QUAL_PRE_ADD ? symbolType2addPhred[LINK_SYMBOL] : 0); // + addidq; 
                    }
                }
                if (!is_ins_at_read_end) {
                    this->template inc<TUpdateType>(rpos, insLenToSymbol(inslen), MAX(SIGN2UNSIGN(1), incvalue), b);
                    if (TFillSeqDir) {
                        dealwith_seg_regbias(qpos + cigar_oplen, qpos, b, sidereg_nbases, bq_dirs_count, strand, isrc, bq_regs_count, rpos, insLenToSymbol(inslen), bq_nms_count, nm,
                                bq_baq_sum, qr_baq_vec, bq_dirs_bqsum, 0);
                    }
                    std::string iseq;
                    iseq.reserve(cigar_oplen);
                    unsigned int incvalue2 = incvalue;
                    for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                        unsigned int base4bit = bam_seqi(bseq, qpos+i2);
                        const char base8bit = seq_nt16_str[base4bit];
                        iseq.push_back(base8bit);
                        if (TUpdateType == BASE_QUALITY_MAX) {
                            incvalue2 = MIN(incvalue2, SIGN2UNSIGN(bam_seqi(bseq, qpos+i2)))
                                    + (QUAL_PRE_ADD ? symbolType2addPhred[LINK_SYMBOL] : 0); // + symbolType2addPhred[LINK_SYMBOL];
                        }
                    }
                    this->incIns(rpos, iseq, insLenToSymbol(inslen), MAX(SIGN2UNSIGN(1), incvalue2));
                }
                qpos += cigar_oplen;
            } else if (cigar_op == BAM_CDEL) {
                const bool is_del_at_read_end = (0 == qpos || qpos + SIGN2UNSIGN(cigar_oplen) >= SIGN2UNSIGN(b->core.l_qseq));
                unsigned int dellen = SIGN2UNSIGN(cigar_oplen);
                if (TUpdateType == BASE_QUALITY_MAX) {
                    if (TIndelAddPhred) {
                        unsigned int addidq = (THasDups ? 0 : SIGN2UNSIGN(MIN(cigar_oplen - 1, SIGN2UNSIGN(3 - 1)) * frag_indel_ext));
                        incvalue = TIndelAddPhred + addidq;
                    } else if (is_del_at_read_end) {
                        LOG(logWARNING) << "Query " << bam_get_qname(b) << " has deletion of legnth " << cigar_oplen << " at " << qpos
                                << " which is not exclusively between 0 and " << b->core.l_qseq << " aligned to tid " << b->core.tid << " and position " << rpos; 
                        incvalue = (0 != qpos ? bam_phredi(b, qpos-1) : 
                                ((qpos + cigar_oplen < SIGN2UNSIGN(b->core.l_qseq)) ? 
                                bam_phredi(b, qpos + SIGN2UNSIGN(cigar_oplen)) : 1))
                                + (QUAL_PRE_ADD ? symbolType2addPhred[LINK_SYMBOL] : 0); // + addidq;
                    } else {
                        // double afa = ((cigar_oplen <= 2) ? 18.0 : 6.0);
                        unsigned int phredvalue = ref_to_phredvalue(dellen, region_symbolvec, rpos - region_offset, 
                                frag_indel_basemax, 8.0, cigar_oplen, cigar_op, b);
                        // unsigned int phredvalue = bam_to_phredvalue(dellen, b, qpos, frag_indel_basemax, 6.0, cigar_oplen, cigar_op); // THasDups is not used here
                        incvalue = MIN(MIN(bam_phredi(b, qpos), bam_phredi(b, qpos-1)) + (TIsProton ? proton_cigarlen2phred(cigar_oplen) : 0), phredvalue) // + addidq; 
                                + (QUAL_PRE_ADD ? symbolType2addPhred[LINK_SYMBOL] : 0);
                        // + symbolType2addPhred[LINK_SYMBOL];
                    }
                }
                if (!is_del_at_read_end) {
                    this->template inc<TUpdateType>(rpos, delLenToSymbol(dellen), MAX(SIGN2UNSIGN(1), incvalue), b);
                    if (TFillSeqDir) {
                        dealwith_seg_regbias(qpos, qpos, b, sidereg_nbases, bq_dirs_count, strand, isrc, bq_regs_count, rpos, delLenToSymbol(dellen), bq_nms_count, nm,
                                bq_baq_sum, qr_baq_vec, bq_dirs_bqsum, 0);
                    }
                    this->incDel(rpos, cigar_oplen, delLenToSymbol(dellen), MAX(SIGN2UNSIGN(1), incvalue));
                }
#if 1 
// The definition of non-ref for indel is not clearly defined, this piece of code can result in germline risk that is not in the ground truth, so it it commented out.
// However, if we consider any position covered by the boundary of a germline indel to be non-ref, then this code should be enabled.
                unsigned int endpos = SIGN2UNSIGN(bam_endpos(b));
                for (unsigned int p = rpos+1; p < MIN(rpos + cigar_oplen + 1, endpos); p++) {
                    this->template inc<TUpdateType>(p, LINK_NN, MAX(SIGN2UNSIGN(1), incvalue), b);
                }
                for (unsigned int p = rpos; p < MIN(rpos + cigar_oplen, endpos); p++) {
                    this->template inc<TUpdateType>(p, BASE_NN, MAX(SIGN2UNSIGN(1), incvalue), b);
                }
#endif
                rpos += cigar_oplen;
            } else if (cigar_op == BAM_CREF_SKIP) {
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
        return 0;
    }

    template<ValueType TUpdateType, bool TFillSeqDir>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByRead1Aln(
            std::vector<bam1_t *> aln_vec, 
            unsigned int frag_indel_ext, 
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred, 
            const unsigned int alns2size, 
            const unsigned int frag_indel_basemax, 
            unsigned int dflag, 
            unsigned int nogap_phred, 
            const bool is_proton, 
            const auto & region_symbolvec,
            const unsigned int region_offset,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 4> & bq_dirs_count, 
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 3> & bq_regs_count,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 4> & bq_nms_count,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 1> & bq_baq_sum,
            std::array<GenericSymbol2CountCoverage<Symbol2Count>, 4> & bq_dirs_bqsum,
            unsigned int sidereg_nbases,
            const auto & region_repeatvec,
            unsigned int baq_per_aligned_base,
            unsigned int specialflag) {
        for (bam1_t *aln : aln_vec) {
            if (alns2size > 1 && dflag > 0) { // is barcoded and not singleton
                if (is_proton) {
                    this->template updateByAln<TUpdateType, true , true , TFillSeqDir>(aln, frag_indel_ext, symbolType2addPhred, frag_indel_basemax, 
                            nogap_phred, region_symbolvec, region_offset, bq_dirs_count, bq_regs_count, bq_nms_count, bq_baq_sum, bq_dirs_bqsum, sidereg_nbases, region_repeatvec, baq_per_aligned_base);
                } else {
                    this->template updateByAln<TUpdateType, false, true , TFillSeqDir>(aln, frag_indel_ext, symbolType2addPhred, frag_indel_basemax, 
                            nogap_phred, region_symbolvec, region_offset, bq_dirs_count, bq_regs_count, bq_nms_count, bq_baq_sum, bq_dirs_bqsum, sidereg_nbases, region_repeatvec, baq_per_aligned_base);
                }
            } else {
                if (is_proton) {
                    this->template updateByAln<TUpdateType, true , false, TFillSeqDir>(aln, frag_indel_ext, symbolType2addPhred, frag_indel_basemax, 
                            nogap_phred, region_symbolvec, region_offset, bq_dirs_count, bq_regs_count, bq_nms_count, bq_baq_sum, bq_dirs_bqsum, sidereg_nbases, region_repeatvec, baq_per_aligned_base);
                } else {
                    this->template updateByAln<TUpdateType, false, false, TFillSeqDir>(aln, frag_indel_ext, symbolType2addPhred, frag_indel_basemax, 
                            nogap_phred, region_symbolvec, region_offset, bq_dirs_count, bq_regs_count, bq_nms_count, bq_baq_sum, bq_dirs_bqsum, sidereg_nbases, region_repeatvec, baq_per_aligned_base);
                }
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
    
    std::array<Symbol2CountCoverage, 4> bq_dirs_count;
    std::array<Symbol2CountCoverage, 4> bq_nms_count;
    std::array<Symbol2CountCoverage, 1> bq_baq_sum;
    std::array<Symbol2CountCoverage, 3> bq_regs_count;
    std::array<Symbol2CountCoverage, 4> bq_bias_sedir;
    std::array<Symbol2CountCoverage, 4> bq_dirs_bqsum;
    
    // std::array<Symbol2CountCoverage, 2> bq_qsum_rawMQ;
    std::array<Symbol2CountCoverageUint64, 2> bq_qsum_sqrMQ;
    std::array<Symbol2CountCoverage, 2> bq_qual_p1sum;
    std::array<Symbol2CountCoverageUint64, 2> bq_qual_p2sum;
    std::array<Symbol2CountCoverage, 2> bq_tsum_LQdep;
    // std::array<Symbol2CountCoverage, 2> bq_n_edit_ops;
    
    std::array<Symbol2CountCoverage, 2> du_bias_dedup;  
    
    std::array<Symbol2CountCoverage, 2> bq_amax_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> bq_bias_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> bq_amax_rdist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> bq_bias_rdist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> bq_amax_nvars; // number of mismatches
    std::array<Symbol2CountCoverage, 2> bq_bias_nvars; // number of mismatches

    std::array<Symbol2CountCoverage, 2> bq_bsum_ldist;
    std::array<Symbol2CountCoverage, 2> bq_bsum_rdist;
    std::array<Symbol2CountCoverage, 2> bq_bias_1stra;
    std::array<Symbol2CountCoverage, 2> bq_bias_2stra;
    std::array<Symbol2CountCoverage, 2> bq_n_seq_bases;

    std::array<Symbol2CountCoverage, 2> bq_tsum_depth;
    std::array<Symbol2CountCoverage, 2> bq_pass_thres;
    std::array<Symbol2CountCoverage, 2> bq_pass_depth;
    std::array<Symbol2CountCoverage, 2> bq_vars_thres;
    std::array<Symbol2CountCoverage, 2> bq_vars_depth;
    std::array<Symbol2CountCoverage, 2> bq_vars_badep;
    std::array<Symbol2CountCoverage, 2> bq_vars_vqual;
        
    std::array<Symbol2CountCoverage, 2> major_amplicon; // duped, is needed by distr
    std::array<Symbol2CountCoverage, 2> minor_amplicon; // duped, is needed by distr
    std::array<Symbol2CountCoverage, 2> fam_total_dep; // deduped
    std::array<Symbol2CountCoverage, 2> fam_size1_dep; // deduped 
    std::array<Symbol2CountCoverage, 2> fam_nocon_dep; // deduped 
    
    std::array<Symbol2CountCoverage, 2> fq_qual_p1sum;
    // std::array<Symbol2CountCoverageUint64, 2> fq_qual_p2sum;
    std::array<Symbol2CountCoverage, 2> fq_hiqual_dep;

    std::array<Symbol2CountCoverage, 2> fq_amax_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> fq_bias_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> fq_amax_rdist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> fq_bias_rdist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> fq_amax_nvars; // number of mismatches
    std::array<Symbol2CountCoverage, 2> fq_bias_nvars; // number of mismatches

    std::array<Symbol2CountCoverage, 2> fq_bsum_ldist;
    std::array<Symbol2CountCoverage, 2> fq_bsum_rdist;
    std::array<Symbol2CountCoverage, 2> fq_bias_1stra;
    std::array<Symbol2CountCoverage, 2> fq_bias_2stra;
    std::array<Symbol2CountCoverage, 2> fq_n_seq_bases;

    std::array<Symbol2CountCoverage, 2> fq_tsum_depth;
    std::array<Symbol2CountCoverage, 2> fq_pass_thres;
    std::array<Symbol2CountCoverage, 2> fq_pass_depth;
    std::array<Symbol2CountCoverage, 2> fq_vars_thres;
    std::array<Symbol2CountCoverage, 2> fq_vars_depth;
    std::array<Symbol2CountCoverage, 2> fq_vars_badep;
    std::array<Symbol2CountCoverage, 2> fq_vars_vqual;
   
    Symbol2CountCoverage duplex_pass_depth;
    Symbol2CountCoverage duplex_tsum_depth;
    
    std::array<Symbol2Bucket2CountCoverage, 2> dedup_ampDistr; // dedup, important for error correction by barcoding, need major and minor, second pass
    std::array<Symbol2Bucket2CountCoverageEdgeDist, 2> pb_dist_lpart;
    std::array<Symbol2Bucket2CountCoverageEdgeDist, 2> pb_dist_rpart;
    std::array<Symbol2Bucket2CountCoverageNumMisma, 2> pb_dist_nvars;
    
    Symbol2CountCoverageString additional_note;
    
    Symbol2CountCoverageSet(unsigned int t, unsigned int beg, unsigned int end):
        tid(t), incluBegPosition(beg), excluEndPosition(end)
        , bq_dirs_count({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end),
                         Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_nms_count ({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end),
                         Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_baq_sum   ({{Symbol2CountCoverage(t, beg, end)}})
        , bq_regs_count({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_sedir({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end),
                         Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_dirs_bqsum({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end),
                         Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        // , bq_qsum_rawMQ({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_qsum_sqrMQ({Symbol2CountCoverageUint64(t, beg, end), Symbol2CountCoverageUint64(t, beg, end)})
        , bq_qual_p1sum({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_qual_p2sum({Symbol2CountCoverageUint64(t, beg, end), Symbol2CountCoverageUint64(t, beg, end)})
        , bq_tsum_LQdep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        // , bq_n_edit_ops({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , du_bias_dedup({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , bq_amax_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_amax_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_amax_nvars({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_nvars({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})

        , bq_bsum_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bsum_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_1stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_2stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_n_seq_bases({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})

        , bq_tsum_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_pass_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_pass_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_badep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_vqual({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , major_amplicon({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , minor_amplicon({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fam_total_dep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fam_size1_dep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fam_nocon_dep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , fq_qual_p1sum({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        // , fq_qual_p2sum({Symbol2CountCoverageUint64(t, beg, end), Symbol2CountCoverageUint64(t, beg, end)})
        , fq_hiqual_dep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})

        , fq_amax_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_amax_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_amax_nvars({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_nvars({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
 
        , fq_bsum_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bsum_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_1stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_2stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_n_seq_bases({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , fq_tsum_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_pass_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_pass_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_badep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_vqual({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , duplex_pass_depth(Symbol2CountCoverage(t, beg, end))
        , duplex_tsum_depth(Symbol2CountCoverage(t, beg, end))
        , dedup_ampDistr({Symbol2Bucket2CountCoverage(t, beg, end), Symbol2Bucket2CountCoverage(t, beg, end)})
        , pb_dist_lpart({Symbol2Bucket2CountCoverageEdgeDist(t, beg, end), Symbol2Bucket2CountCoverageEdgeDist(t, beg, end)})
        , pb_dist_rpart({Symbol2Bucket2CountCoverageEdgeDist(t, beg, end), Symbol2Bucket2CountCoverageEdgeDist(t, beg, end)})
        , pb_dist_nvars({Symbol2Bucket2CountCoverageNumMisma(t, beg, end), Symbol2Bucket2CountCoverageNumMisma(t, beg, end)})
        , additional_note(Symbol2CountCoverageString(t, beg, end))
    {
        assert(beg < end);
        assert(this->bq_tsum_depth[0].getIncluBegPosition() == beg);
        assert(this->bq_tsum_depth[1].getIncluBegPosition() == beg);
        Symbol2CountCoverage cov(t, beg, end);
    };
    
    template <bool TIsFilterStrong=false>
    int 
    getbest(auto & max_pqual, auto & best_phred, auto & best_count,
            const auto & ampDistrByPos, const double symbolTypeSum, const AlignmentSymbol symbol, const unsigned int bias_adjusted_mincount, double imba_fact,
            const unsigned int phred_max, const unsigned int add_phred, double ess_georatio_dedup, const double homogeneity = 0) const {
        max_pqual = 0;
        best_phred = 0;
        best_count = 0;
        unsigned int tot_count = 0;
        double prev_adj_tot_count = DBL_MAX;
        for (unsigned int rev_buc_idx = 0; rev_buc_idx < NUM_BUCKETS; rev_buc_idx++) {
            unsigned int bucket = NUM_BUCKETS - 1 - rev_buc_idx;
            unsigned int count = ampDistrByPos.getSymbolBucketCount(symbol, bucket);
            tot_count += count;
            unsigned int phred = MIN(bucket2phred(bucket), phred_max);
            auto tot_pqual = 0;
            assert(tot_count <= symbolTypeSum || !fprintf(stderr, "%d <= %f failed for symbol %d and bucket %d !!!\n", tot_count, symbolTypeSum, symbol, bucket));
            if (0 < count) {
                if (TIsFilterStrong) {
                    const double adj_tot_count = MIN((double)tot_count, bias_adjusted_mincount / 10.0);
                    if (adj_tot_count != prev_adj_tot_count) {
                        // NOTE: add_phred (addPhred) is disabled here, it should be enabled when reading bam
                        if (phred > 0) {
                            double pr = phred2prob(phred + (QUAL_PRE_ADD ? 0 : add_phred));
                            tot_pqual = MAX(tot_pqual, h01_to_phredlike<false>(pr, 1 + DBL_EPSILON, 
                                    adj_tot_count, symbolTypeSum, 1.0, ess_georatio_dedup));
                        }
                        prev_adj_tot_count = adj_tot_count;
                    }
#if 0 // QUESTION: is this piece of code useless? It was only run half of the time without its intended effect due to a bug and yet results generated by the bug look very good
                    bool some_is_biased = (tot_count - count > (bias_adjusted_mincount / 10.0));
                    if (some_is_biased) {
                        double pen = (tot_count - count - bias_adjusted_mincount / 10.0) * 10.0/log(10.0) * log(MAX(1.0, imba_fact)); // penetratation into the bias
                        double phred = phred + (QUAL_PRE_ADD ? 0 : add_phred) - pen;
                        if (phred > 0) {
                            double pr = phred2prob(phred);
                            tot_pqual = MAX(tot_pqual, h01_to_phredlike<false>(pr, 1 + DBL_EPSILON, 
                                    tot_count, symbolTypeSum, 1.0, ess_georatio_dedup));
                        }
                    }
#endif
                } else {
                    tot_pqual = tot_count * phred;
                }
                if (max_pqual < tot_pqual) { // must be less than
                    max_pqual = tot_pqual;
                    best_phred = phred;
                    best_count = tot_count;
                }
            }
        }
        return 0;
    };
    
    template <bool TSmallerIsOutlier>
    std::pair<unsigned int, unsigned int>
    adabias(const auto & t0v, const auto & t1v, const double pseudocount, unsigned int gapdist) {
        static_assert(t0v.size() == t1v.size());
        static_assert(t0v.size() >= 2);
        unsigned int argmax = 0; 
        const size_t usize = t0v.size();
        double sum0 = 0;
        double sum1 = 0;
        for (size_t i = 0; i < usize; i++) {
            sum0 += t0v.at(i);
            sum1 += t1v.at(i);
        }
        double cur0 = 0;
        double cur1 = 0;
        unsigned int max_biasfact100 = 0;
        std::vector<unsigned int> prev_biasfact100s(gapdist, 0);
        for (size_t i = 0; i < usize - 1; i++) {
            cur0 += (double)t0v[i];
            cur1 += (double)t1v[i];
            unsigned int curr_biasfact100 = (TSmallerIsOutlier
                    ? any4_to_biasfact100(sum0 - cur0, cur0, sum1 - cur1, cur1, false, pseudocount) // position bias
                    : any4_to_biasfact100(cur0, sum0 - cur0, cur1, sum1 - cur1, false, pseudocount)); // mismatch bias
            // LOG(logINFO) << "At " << i << " : " << cur0 << " , " << sum0 - cur0 << " , " << cur1 << " , " << sum1 - cur1;
            unsigned int norm_biasfact100 = curr_biasfact100;
            for (unsigned int prev_bf100 : prev_biasfact100s) {
                norm_biasfact100 = MIN(norm_biasfact100, prev_bf100);
            }
            if (norm_biasfact100 > max_biasfact100) {
                max_biasfact100 = norm_biasfact100;
                argmax = i+1;
            }
            for (int j = ((int)gapdist) - 1; j > 0; j--) {
                prev_biasfact100s[j] = prev_biasfact100s[j-1];
            }
            if (0 < gapdist) {  prev_biasfact100s[0] = curr_biasfact100; }
        }
        return std::make_pair(argmax, max_biasfact100);
    }
    
    template<bool TUsePrev = false>
    int 
    adafilter(
            auto & du_bias_dedup,
            const auto & bq_qual_p1sum, const auto & bq_tsum_depth,
            auto & amax_ldist, auto & amax_rdist, auto & bias_ldist, auto & bias_rdist, 
            auto & amax_nvars, auto & bias_nvars, 
            const auto & bsum_ldist, const auto & bsum_rdist, auto & bias_1stra, auto & bias_2stra,
            const auto & curr_tsum_depth, 
            auto & pass_thres, auto & pass_depth,
            auto & vars_thres, auto & vars_depth, auto & vars_badep, auto & vars_vqual,
            auto & dedup_ampDistr, const auto & prev_tsum_depth,
            auto & pb_dist_lpart, auto & pb_dist_rpart, auto & pb_dist_nvars, auto & additional_note, bool should_add_note, 
            const auto & bq_dirs_count, auto & bq_bias_sedir, const auto & q_n_seq_bases, const AssayType assay_type,
            const PhredMutationTable & phred_max_table, const auto & symbolType2addPhred, const double ess_georatio_dedup, 
            const unsigned int uni_bias_r_max, uint32_t bias_flag_snv, uint32_t bias_flag_indel) {
        
        assert(dedup_ampDistr.at(0).getIncluBegPosition() == dedup_ampDistr.at(1).getIncluBegPosition());
        assert(dedup_ampDistr.at(0).getExcluEndPosition() == dedup_ampDistr.at(1).getExcluEndPosition());
        for (unsigned int strand = 0; strand < 2; strand++) {
            for (auto pos = dedup_ampDistr.at(strand).getIncluBegPosition(); pos < dedup_ampDistr.at(strand).getExcluEndPosition(); pos++) {
                for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+(unsigned int)symbolType)) {
                    const auto bias_flag_symb = (BASE_SYMBOL == symbolType ? bias_flag_snv : bias_flag_indel);
                    // prepare duplication bias
                    const auto prev_depth_typesum = (TUsePrev ? prev_tsum_depth[strand].getByPos(pos).sumBySymbolType(symbolType) : 0);
                    const auto curr_depth_typesum = curr_tsum_depth[0+strand].getByPos(pos).sumBySymbolType(symbolType); 
                    const auto curr_deprv_typesum = curr_tsum_depth[1-strand].getByPos(pos).sumBySymbolType(symbolType);
 
                    // prepare positional bias
                    Bucket2CountEdgeDist vsum_pb_dist_lpart = pb_dist_lpart[strand].getByPos(pos).vectorsumBySymbolType(symbolType);
                    Bucket2CountEdgeDist vsum_pb_dist_rpart = pb_dist_rpart[strand].getByPos(pos).vectorsumBySymbolType(symbolType);
                    Bucket2CountNumMisma vsum_pb_dist_nvars = pb_dist_nvars[strand].getByPos(pos).vectorsumBySymbolType(symbolType);

                    // prepare strand bias
                    auto typesum_uqual_v0 = bq_qual_p1sum[1-strand].getByPos(pos).sumBySymbolType(symbolType);
                    auto typesum_depth_v0 = bq_tsum_depth[1-strand].getByPos(pos).sumBySymbolType(symbolType);
                    auto typesum_uqual_v1 = bq_qual_p1sum[0+strand].getByPos(pos).sumBySymbolType(symbolType);
                    auto typesum_depth_v1 = bq_tsum_depth[0+strand].getByPos(pos).sumBySymbolType(symbolType);
                    double typesum_uqual_v0_avg = typesum_uqual_v0/ (double)(typesum_depth_v0 + DBL_MIN);
                    
                    auto typebsum_ldist_v0 = bsum_ldist[1-strand].getByPos(pos).sumBySymbolType(symbolType);
                    auto typebsum_rdist_v0 = bsum_rdist[1-strand].getByPos(pos).sumBySymbolType(symbolType); 
                    
                    const auto dp0 = curr_tsum_depth.at(1-strand).getByPos(pos).sumBySymbolType(symbolType);
                    const auto dp1 = curr_tsum_depth.at(0+strand).getByPos(pos).sumBySymbolType(symbolType);
                    
                    // prepare sequencing-segment biases
                    auto bq_dir_s0 = bq_dirs_count.at(0*2+0).getByPos(pos).sumBySymbolType(symbolType)
                                   + bq_dirs_count.at(1*2+0).getByPos(pos).sumBySymbolType(symbolType);
                    auto bq_dir_s1 = bq_dirs_count.at(0*2+1).getByPos(pos).sumBySymbolType(symbolType)
                                   + bq_dirs_count.at(1*2+1).getByPos(pos).sumBySymbolType(symbolType);
                    auto side_stotal= bq_dir_s0 + bq_dir_s1;
                    auto only_sboth = bq_regs_count.at(0).getByPos(pos).sumBySymbolType(symbolType);
                    auto only_sleft = bq_regs_count.at(1).getByPos(pos).sumBySymbolType(symbolType);
                    auto only_sright= bq_regs_count.at(2).getByPos(pos).sumBySymbolType(symbolType);
                    auto side_sleft = only_sleft + only_sboth;
                    auto side_sright= only_sright+ only_sboth;
                    auto side_sboth = only_sleft + only_sright+ only_sboth;

                    for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType];
                            symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType];
                            symbol = AlignmentSymbol(1+((unsigned int)symbol))) {
                        auto curr_depth_symbsum = curr_tsum_depth[0+strand].getByPos(pos).getSymbolCount(symbol);
                        auto curr_deprv_symbsum = curr_tsum_depth[1-strand].getByPos(pos).getSymbolCount(symbol);
                        unsigned int max_imba_depth = (MAX_IMBA_DEP); // magic number meaning no limit on imba depth
                        double imba_fact = 1.0;
if (SYMBOL_TYPE_TO_AMBIG[symbolType] != symbol 
        && ((curr_depth_symbsum * 5 < curr_depth_typesum * 4 && curr_depth_symbsum > 0)
         || (curr_deprv_symbsum * 5 < curr_deprv_typesum * 4 && curr_deprv_symbsum > 0))) {
                        const unsigned int add1count = 1;
                        const double pseudocount = (double)add1count;
                        // compute duplication bias
                        double dup_imba = 1;
                        if (TUsePrev) {
                            auto prev_depth_symbsum = prev_tsum_depth[strand].getByPos(pos).getSymbolCount(symbol);
                            auto db100 = any4_to_biasfact100(
                                    MAX(prev_depth_typesum, curr_depth_typesum) - curr_depth_typesum + add1count,
                                    curr_depth_typesum, 
                                    MAX(prev_depth_symbsum, curr_depth_symbsum) - curr_depth_symbsum + add1count,
                                    curr_depth_symbsum, 
                                    false, pseudocount / 2.0);
                            
                            du_bias_dedup[strand].getRefByPos(pos).incSymbolCount(symbol, db100);
                            dup_imba = biasfact100_to_imba(db100);
                        }
                        
                        auto ad0 = curr_tsum_depth.at(1-strand).getByPos(pos).getSymbolCount(symbol);
                        auto ad1 = curr_tsum_depth.at(0+strand).getByPos(pos).getSymbolCount(symbol);
                        
                        unsigned int ad01 = ad0 + ad1;
                        double fa = (double)(ad0 + ad1) / (double)(dp0 + dp1 + 1);
                        // unsigned int nvars_gapdist = ((fa < 0.025) ? 4 : ((fa < 0.025 * 5) ? 5 : 6));
                        unsigned int nvars_gapdist = q_n_seq_bases[strand].getByPos(pos).getSymbolCount(symbol) / MAX(ad1, 1) / 90
                                + (ad01 < 4 ? 1 : (ad01 < 16 ? 2 : 3));
                        // compute positional bias
                        auto pb_ldist_pair = adabias<true >(vsum_pb_dist_lpart, pb_dist_lpart[strand].getByPos(pos).getSymbolCounts(symbol), pseudocount / 2.0, 2);
                        amax_ldist[strand].getRefByPos(pos).incSymbolCount(symbol, edbuck2pos(pb_ldist_pair.first));
                        bias_ldist[strand].getRefByPos(pos).incSymbolCount(symbol, pb_ldist_pair.second);
                        auto pb_ldist_imba = biasfact100_to_imba(bias_ldist[strand].getRefByPos(pos).getSymbolCount(symbol));
                        
                        auto pb_rdist_pair = adabias<true >(vsum_pb_dist_rpart, pb_dist_rpart[strand].getByPos(pos).getSymbolCounts(symbol), pseudocount / 2.0, 2);
                        amax_rdist[strand].getRefByPos(pos).incSymbolCount(symbol, edbuck2pos(pb_rdist_pair.first));
                        bias_rdist[strand].getRefByPos(pos).incSymbolCount(symbol, pb_rdist_pair.second);
                        auto pb_rdist_imba = biasfact100_to_imba(bias_rdist[strand].getRefByPos(pos).getSymbolCount(symbol));
                        
                        auto pb_nvars_pair = adabias<false>(vsum_pb_dist_nvars, pb_dist_nvars[strand].getByPos(pos).getSymbolCounts(symbol), pseudocount / 2.0, nvars_gapdist);
                        amax_nvars[strand].getRefByPos(pos).incSymbolCount(symbol, pb_nvars_pair.first);
                        bias_nvars[strand].getRefByPos(pos).incSymbolCount(symbol, pb_nvars_pair.second);
                        auto pb_nvars_imba = biasfact100_to_imba(bias_nvars[strand].getRefByPos(pos).getSymbolCount(symbol));
                        
                        if (should_add_note) {
                            unsigned int allcurr, altcurr, allrest, altrest;
                            this->additional_note.getRefByPos(pos).at(symbol) += "//(";
                            
                            allcurr = altcurr = allrest = altrest = 0;
                            for (unsigned int i = 0; i < vsum_pb_dist_lpart.size(); i++) {
                                allrest += vsum_pb_dist_lpart[i];
                                altrest += pb_dist_lpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                            }
                            for (unsigned int i = 0; i < vsum_pb_dist_lpart.size(); i++) {
                                allcurr += vsum_pb_dist_lpart[i];
                                altcurr += pb_dist_lpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                                this->additional_note.getRefByPos(pos).at(symbol) += std::to_string(i) + "(" + std::to_string(edbuck2pos(i)) + "/" 
                                        + std::to_string(allrest-allcurr) + "/" + std::to_string(allcurr) + "/" + std::to_string(altrest-altcurr) + "/" + std::to_string(altcurr)  + ")";
                            }

                            allcurr = altcurr = allrest = altrest = 0;
                            for (unsigned int i = 0; i < vsum_pb_dist_rpart.size(); i++) {
                                allrest += vsum_pb_dist_rpart[i];
                                altrest += pb_dist_rpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                            }
                            for (unsigned int i = 0; i < vsum_pb_dist_rpart.size(); i++) {
                                allcurr += vsum_pb_dist_rpart[i];
                                altcurr += pb_dist_rpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                                this->additional_note.getRefByPos(pos).at(symbol) += std::to_string(i) + "(" + std::to_string(edbuck2pos(i)) + "/" 
                                        + std::to_string(allrest-allcurr) + "/" + std::to_string(allcurr) + "/" + std::to_string(altrest-altcurr) + "/" + std::to_string(altcurr)  + ")";
                            }
                            
                            allcurr = altcurr = allrest = altrest = 0;
                            for (unsigned int i = 0; i < vsum_pb_dist_nvars.size(); i++) {
                                allrest += vsum_pb_dist_nvars[i];
                                altrest += pb_dist_nvars[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                            }
                            for (unsigned int i = 0; i < vsum_pb_dist_nvars.size(); i++) {
                                allcurr += vsum_pb_dist_nvars[i];
                                altcurr += pb_dist_nvars[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                                this->additional_note.getRefByPos(pos).at(symbol) += std::to_string(i) + "(" + std::to_string((i)) + "/" 
                                        + std::to_string(allrest-allcurr) + "/" + std::to_string(allcurr) + "/" + std::to_string(altrest-altcurr) + "/" + std::to_string(altcurr)  + ")";
                            }
                            
                            this->additional_note.getRefByPos(pos).at(symbol) += ")//";
                            this->additional_note.getRefByPos(pos).at(symbol) += "ampDistr:"+std::to_string(this->dedup_ampDistr[strand].getByPos(pos).getSymbolCounts(symbol).size())+":";
                            for (size_t k = 0; k < this->dedup_ampDistr[strand].getByPos(pos).getSymbolCounts(symbol).size(); k++) {
                                this->additional_note.getRefByPos(pos).at(symbol) += std::to_string(this->dedup_ampDistr[strand].getByPos(pos).getSymbolCounts(symbol).at(k)) + ",";
                            }
                        }

                        // compute strand bias
                        auto symbval_uqual_v0 = bq_qual_p1sum[1-strand].getByPos(pos).getSymbolCount(symbol);
                        auto symbval_depth_v0 = bq_tsum_depth[1-strand].getByPos(pos).getSymbolCount(symbol);
                        auto symbval_uqual_v1 = bq_qual_p1sum[0+strand].getByPos(pos).getSymbolCount(symbol);
                        auto symbval_depth_v1 = bq_tsum_depth[0+strand].getByPos(pos).getSymbolCount(symbol);
                        double symbval_uqual_v1_avg = symbval_uqual_v1 / (double)(symbval_depth_v1 + DBL_MIN);
                        auto uqual_avg_imba = pow((double)10, (symbval_uqual_v1_avg - typesum_uqual_v0_avg) / (double)10);
                        
                        double symbval_uqual_va_avg = (symbval_uqual_v0 + symbval_uqual_v1) / (double)(symbval_depth_v0 + symbval_depth_v1 + DBL_MIN);
                        double typesum_uqual_va_avg = (typesum_uqual_v0 + typesum_uqual_v1) / (double)(typesum_depth_v0 + typesum_depth_v1 + DBL_MIN);
                        
                        auto bsum_ldist_v1 = bsum_ldist[0+strand].getByPos(pos).getSymbolCount(symbol);
                        auto bsum_rdist_v1 = bsum_rdist[0+strand].getByPos(pos).getSymbolCount(symbol);
                        auto bsum_dist_imba0 = ((bsum_ldist_v1 + 1) / (double)(ad1 + 1)) / ((typebsum_rdist_v0 + 1) / (double)(dp0 + 1));
                        auto bsum_dist_imba1 = ((bsum_rdist_v1 + 1) / (double)(ad1 + 1)) / ((typebsum_ldist_v0 + 1) / (double)(dp0 + 1));
                        
                        auto sb100raw = any4_to_biasfact100((double)dp0, (double)dp1, (double)ad0, (double)ad1, false, pseudocount / 2.0 + MIN(2.0, fa * 100.0));
                        bias_1stra[strand].getRefByPos(pos).incSymbolCount(symbol, sb100raw);
                        assert(bsum_dist_imba0 + bsum_dist_imba1 > 0 || !fprintf(stderr, "%g + %g > 0 failed! (will encounter division by zero)\n", bsum_dist_imba0, bsum_dist_imba1));
                        
                        auto sb100fin = (unsigned int)(sb100raw / MAX(uqual_avg_imba, MAX(bsum_dist_imba0, bsum_dist_imba1)));
                        bias_2stra[strand].getRefByPos(pos).incSymbolCount(symbol, sb100fin);
                        auto str_imba = biasfact100_to_imba(sb100fin);
                        
                        auto bq_dir_v0 = bq_dirs_count.at(0*2+0).getByPos(pos).getSymbolCount(symbol) 
                                       + bq_dirs_count.at(1*2+0).getByPos(pos).getSymbolCount(symbol);
                        auto bq_dir_v1 = bq_dirs_count.at(0*2+1).getByPos(pos).getSymbolCount(symbol)
                                       + bq_dirs_count.at(1*2+1).getByPos(pos).getSymbolCount(symbol);
                        unsigned int segbias100 = 100;
                        if (!TUsePrev && 0 == strand) {
                            const bool test1bias = (bq_dir_v0 * bq_dir_s1 < bq_dir_v1 * bq_dir_s0);
                            const auto bq_dir_vmin = (test1bias ? bq_dir_v0 : bq_dir_v1);
                            const auto bq_dir_vmax = (test1bias ? bq_dir_v1 : bq_dir_v0);
                            const auto bq_dir_smin = (test1bias ? bq_dir_s0 : bq_dir_s1);
                            const auto bq_dir_smax = (test1bias ? bq_dir_s1 : bq_dir_s0);
                            const double pseudocount_1 = ((double)pseudocount) / 4.0 + (double)(symbolType == BASE_SYMBOL ? 
                                    (MIN(2.0, 2.0 * pow(10.0, (3 + symbval_uqual_va_avg - typesum_uqual_va_avg) / 10.0))) : 1);
                            auto segbias100strand = any4_to_biasfact100((double)bq_dir_smin, (double)bq_dir_smax, (double)bq_dir_vmin, (double)bq_dir_vmax, 
                                    false, pseudocount_1);
                            bq_bias_sedir[0].getRefByPos(pos).incSymbolCount(symbol, segbias100strand); 
                            
                            auto side_vtotal= bq_dir_v0 + bq_dir_v1;
                            auto only_vboth = bq_regs_count.at(0).getByPos(pos).getSymbolCount(symbol);
                            auto only_vleft = bq_regs_count.at(1).getByPos(pos).getSymbolCount(symbol);
                            auto only_vright= bq_regs_count.at(2).getByPos(pos).getSymbolCount(symbol);
                            auto side_vleft = only_vleft + only_vboth;
                            auto side_vright= only_vright+ only_vboth;
                            auto side_vboth = only_vleft + only_vright+ only_vboth;
                            
                            auto segbias100bside = any4_to_biasfact100(
                                    (double)(side_stotal - side_sboth), (double)side_sboth,
                                    (double)(side_vtotal - side_vboth), (double)side_vboth, false, pseudocount);
                            bq_bias_sedir[1].getRefByPos(pos).incSymbolCount(symbol, segbias100bside);
                            auto segbias100lside = any4_to_biasfact100(
                                    (double)(side_stotal - side_sleft), (double)side_sleft, 
                                    (double)(side_vtotal - side_vleft), (double)side_vleft, false, pseudocount);
                            bq_bias_sedir[2].getRefByPos(pos).incSymbolCount(symbol, segbias100lside);
                            auto segbias100rside = any4_to_biasfact100(
                                    (double)(side_stotal - side_sright),(double)side_sright, 
                                    (double)(side_vtotal - side_vright),(double)side_vright,false, pseudocount * 2.0);
                            bq_bias_sedir[3].getRefByPos(pos).incSymbolCount(symbol, segbias100rside);
                        }
                        if (bias_flag_symb & BIAS_SSEG_STR) {UPDATE_MAX(segbias100, bq_bias_sedir.at(0).getByPos(pos).getSymbolCount(symbol)); }
                        if (bias_flag_symb & BIAS_SSEG_END) {UPDATE_MAX(segbias100, bq_bias_sedir.at(1).getByPos(pos).getSymbolCount(symbol)); }
                        if (bias_flag_symb & BIAS_SSEG_POS) {UPDATE_MAX(segbias100, MAX(
                                bq_bias_sedir.at(2).getByPos(pos).getSymbolCount(symbol), 
                                bq_bias_sedir.at(3).getByPos(pos).getSymbolCount(symbol))); }
                        
                        imba_fact = biasfact100_to_imba(segbias100);
                        if (bias_flag_symb & BIAS_FRAG_DUP) { UPDATE_MAX(imba_fact, dup_imba); }
                        if (bias_flag_symb & BIAS_FRAG_POS) { UPDATE_MAX(imba_fact, MAX(pb_ldist_imba, pb_rdist_imba)); }
                        if (bias_flag_symb & BIAS_FRAG_STR) { UPDATE_MAX(imba_fact, str_imba); }
                        if (bias_flag_symb & BIAS_FRAG_MIS) { UPDATE_MAX(imba_fact, pb_nvars_imba); }
                        
                        max_imba_depth = (unsigned int)ceil(0.5 * 10.0 + 10.0 * curr_depth_symbsum / 
                                MIN(((double)uni_bias_r_max) / 100.0, imba_fact) / (1 + DBL_EPSILON));
                        
                        if (should_add_note) {
                            this->additional_note.getRefByPos(pos).at(symbol) += "//(" +
                                    std::to_string(uqual_avg_imba) + "/" + 
                                    std::to_string(bsum_dist_imba0) + "/" + 
                                    std::to_string(bsum_dist_imba1) + "/" + 
                                    ")//";
                        }
}
                        auto phred_max = 0;
                        // previously it was: if (LINK_SYMBOL == symbolType) { phred_max = phred_max_table.indel_any; } else ...
                        {
                            AlignmentSymbol con_symbol;
                            unsigned int con_count, tot_count;
                            curr_tsum_depth.at(0+strand).getByPos(pos).fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                            phred_max = phred_max_table.toPhredErrRate(con_symbol, symbol);
                        }
                        // find best cutoff from families
                        double max_pqual = 0;
                        unsigned int best_phred = 0;
                        unsigned int best_count = 0;
                        if (curr_depth_symbsum > 0) {
                            getbest<false>(
                                    max_pqual, best_phred, best_count,
                                    dedup_ampDistr[strand].getRefByPos(pos), curr_depth_typesum, symbol, max_imba_depth, imba_fact, phred_max, 0, ess_georatio_dedup);
                        } else {
                            best_phred = 0;
                            best_count = 0;
                        }
                        pass_thres[strand].getRefByPos(pos).incSymbolCount(symbol, best_phred);
                        pass_depth[strand].getRefByPos(pos).incSymbolCount(symbol, best_count);
                        
                        if (curr_depth_symbsum > 0) {
                            getbest<true> (
                                    max_pqual, best_phred, best_count,
                                    dedup_ampDistr[strand].getRefByPos(pos), curr_depth_typesum, symbol, max_imba_depth, imba_fact, phred_max, symbolType2addPhred[symbolType], ess_georatio_dedup);
                        } else {
                            max_pqual = 0;
                            best_phred = 0;
                            best_count = 0;
                        }
                        vars_thres[strand].getRefByPos(pos).incSymbolCount(symbol, best_phred);
                        vars_depth[strand].getRefByPos(pos).incSymbolCount(symbol, best_count);
                        vars_badep[strand].getRefByPos(pos).incSymbolCount(symbol, max_imba_depth); // bias-adjusted-depth
                        vars_vqual[strand].getRefByPos(pos).incSymbolCount(symbol, (unsigned int)max_pqual);
                    }
                }
                dedup_ampDistr[strand].getRefByPos(pos).clearSymbolBucketCount();
                pb_dist_lpart[strand].getRefByPos(pos).clearSymbolBucketCount();
                pb_dist_rpart[strand].getRefByPos(pos).clearSymbolBucketCount();
                pb_dist_nvars[strand].getRefByPos(pos).clearSymbolBucketCount(); 
            }
        }
        return 0;
    };
    
    int 
    updateByAlns3UsingBQ(
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map,
            const std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> & alns3, 
            const std::basic_string<AlignmentSymbol> & region_symbolvec,
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred,
            bool should_add_note, 
            unsigned int frag_indel_ext, 
            unsigned int frag_indel_basemax,
            const PhredMutationTable & phred_max_table, 
            unsigned int phred_thres,
            double ess_georatio_dedup,
            double ess_georatio_duped_pcr,
            unsigned int fixedthresBQ, 
            unsigned int nogap_phred, 
            unsigned int uni_bias_r_max,
            const bool is_proton, 
            const AssayType assay_type,
            const unsigned int baq_per_aligned_base,
            const unsigned int sidereg_nbases,
            uint32_t bias_flag_snv,
            uint32_t bias_flag_indel,
            const auto & region_repeatvec, 
            unsigned int specialflag = 0) {

        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                for (const auto & alns1 : alns2) {
                    uint32_t tid2, beg2, end2;
                    fillTidBegEndFromAlns1(tid2, beg2, end2, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid, beg2, end2);
                    read_ampBQerr_fragWithR1R2.template updateByRead1Aln<BASE_QUALITY_MAX, true>(
                            alns1, 
                            frag_indel_ext, 
                            symbolType2addPhred, 
                            alns2.size(), 
                            frag_indel_basemax, 
                            alns2pair2dflag.second, 
                            nogap_phred, 
                            is_proton, 
                            region_symbolvec, 
                            this->dedup_ampDistr.at(strand).getIncluBegPosition(), 
                            this->bq_dirs_count, 
                            this->bq_regs_count, 
                            this->bq_nms_count,
                            this->bq_baq_sum,
                            this->bq_dirs_bqsum,
                            sidereg_nbases,
                            region_repeatvec,
                            baq_per_aligned_base,
                            0);
                    unsigned int normMQ = 0;
                    for (const bam1_t * b : alns1) {
                        normMQ = MAX(normMQ, b->core.qual);
                    }
                    std::basic_string<std::pair<unsigned int, AlignmentSymbol>> pos_symbol_string;
                    unsigned int ldist_inc = 0;
                    unsigned int rdist_inc = 0;
                    std::vector<unsigned int> posToInsertLen = read_ampBQerr_fragWithR1R2.computeZeroBasedPosToInsLenVec(rdist_inc); // extend
                    unsigned int n_vars = 0;
                    unsigned int n_ops = 0;
                    unsigned int n_seq_bases = 0;
                    std::vector<std::array<AlignmentSymbol, NUM_SYMBOL_TYPES>> con_symbols_vec(
                            read_ampBQerr_fragWithR1R2.getExcluEndPosition() - read_ampBQerr_fragWithR1R2.getIncluBegPosition(),
                            std::array<AlignmentSymbol, NUM_SYMBOL_TYPES>({END_ALIGNMENT_SYMBOLS, END_ALIGNMENT_SYMBOLS}));
                    for (auto epos = read_ampBQerr_fragWithR1R2.getIncluBegPosition(); epos < read_ampBQerr_fragWithR1R2.getExcluEndPosition(); epos++) {
                        unsigned int ldist = 1 + epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition();
                        unsigned int rdist = read_ampBQerr_fragWithR1R2.getExcluEndPosition() - epos;
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
                            unsigned int phredlike = (con_count * 2 - tot_count);
                            
                            // this->bq_qsum_rawMQ [strand].getRefByPos(epos).incSymbolCount(con_symbol, normMQ);
                            this->bq_qsum_sqrMQ [strand].getRefByPos(epos).incSymbolCount(con_symbol, normMQ * normMQ); 
                           
                            // shared between BQ and FQ
                            con_symbols_vec[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()][symbolType] = con_symbol;
                            this->bq_tsum_depth [strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                            // unsigned int edge_baq = MIN(ldist, rdist) * baq_per_aligned_base;
                            unsigned int overallq = phredlike; // MIN(edge_baq, phredlike);
                            this->bq_qual_p1sum [strand].getRefByPos(epos).incSymbolCount(con_symbol, overallq);
                            this->bq_qual_p2sum [strand].getRefByPos(epos).incSymbolCount(con_symbol, mathsquare(overallq)); // specific to BQ
                            unsigned int pbucket = phred2bucket(overallq);
                            assert (pbucket < NUM_BUCKETS || !fprintf(stderr, "%u < %u failed at position %lu and con_symbol %u symboltype %u plusbucket %u\n", 
                                    pbucket,  NUM_BUCKETS, epos, con_symbol, symbolType, SIGN2UNSIGN(symbolType2addPhred[symbolType])));
                            if (isSymbolIns(con_symbol)) {
                                posToIndelToCount_updateByConsensus(this->bq_tsum_depth[strand].getRefPosToIseqToData(con_symbol), read_ampBQerr_fragWithR1R2.getPosToIseqToData(con_symbol), epos, 1);
                            }
                            if (isSymbolDel(con_symbol)) {
                                posToIndelToCount_updateByConsensus(this->bq_tsum_depth[strand].getRefPosToDlenToData(con_symbol), read_ampBQerr_fragWithR1R2.getPosToDlenToData(con_symbol), epos, 1);
                            }
                            AlignmentSymbol refsymbol = region_symbolvec[epos-this->dedup_ampDistr.at(strand).getIncluBegPosition()]; 
                            if (areSymbolsMutated(refsymbol, con_symbol)) {
                                pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
                                if (symbolType == BASE_SYMBOL && phredlike >= phred_thres) {
                                    n_vars++;
                                }
                                if (symbolType == BASE_SYMBOL) {
                                    n_ops += 1; 
                                } else {
                                    n_ops += 2; // gap-opening is one edit operation and gap-extension is another edit operation
                                }
                            }
                            if (symbolType == BASE_SYMBOL) {
                                n_seq_bases++;
                            }
                            this->dedup_ampDistr[strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pbucket, 1);
                            ldist_inc += posToInsertLen[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()];
                            this->pb_dist_lpart [strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pos2edbuck(ldist + ldist_inc), 1);
                            this->pb_dist_rpart [strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pos2edbuck(rdist + rdist_inc), 1);
                            this->bq_bsum_ldist [strand].getRefByPos(epos).incSymbolCount(con_symbol, ldist + ldist_inc);
                            this->bq_bsum_rdist [strand].getRefByPos(epos).incSymbolCount(con_symbol, rdist + rdist_inc);
                            rdist_inc -= posToInsertLen[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()];
                            
                            if (overallq < fixedthresBQ) {
                                this->bq_tsum_LQdep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1); 
                            }
                        }
                    }
                    n_vars = MIN(n_vars, NUM_NMBUCKS - 1);
                    for (auto epos = read_ampBQerr_fragWithR1R2.getIncluBegPosition(); epos < read_ampBQerr_fragWithR1R2.getExcluEndPosition(); epos++) {
                        for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+((unsigned int)symbolType))) {
                            auto con_symbol = con_symbols_vec.at(epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()).at(symbolType);
                            if (END_ALIGNMENT_SYMBOLS != con_symbol) {
                                this->pb_dist_nvars[strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, n_vars, 1);
                                //this->bq_n_edit_ops[strand].getRefByPos(epos).incSymbolCount(con_symbol, 
                                //        (unsigned int)floor(0.5 + 10.0 / log(10.0) * log(n_ops + 1.0)));
                                this->bq_n_seq_bases[strand].getRefByPos(epos).incSymbolCount(con_symbol, n_seq_bases);
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
        adafilter<false>(
                this->du_bias_dedup, 
                this->bq_qual_p1sum, this->bq_tsum_depth, 
                this->bq_amax_ldist, this->bq_amax_rdist, this->bq_bias_ldist, this->bq_bias_rdist,
                this->bq_amax_nvars, this->bq_bias_nvars,
                this->bq_bsum_ldist, this->bq_bsum_rdist, this->bq_bias_1stra, this->bq_bias_2stra,
                this->bq_tsum_depth, 
                this->bq_pass_thres, this->bq_pass_depth,
                this->bq_vars_thres, this->bq_vars_depth, this->bq_vars_badep, this->bq_vars_vqual,
                this->dedup_ampDistr,this->bq_tsum_depth,
                this->pb_dist_lpart, this->pb_dist_rpart, this->pb_dist_nvars, this->additional_note, should_add_note, 
                this->bq_dirs_count, this->bq_bias_sedir, this->bq_n_seq_bases, assay_type,
                phred_max_table, symbolType2addPhred, ess_georatio_dedup, uni_bias_r_max, bias_flag_snv, bias_flag_indel);
        return 0;
    }
    
    int
    updateByAlns3UsingFQ(
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map,
            const auto & alns3, 
            const std::basic_string<AlignmentSymbol> & region_symbolvec,
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred, 
            bool should_add_note,
            const unsigned int frag_indel_ext, 
            const unsigned int frag_indel_basemax, 
            const PhredMutationTable & phred_max_table, 
            unsigned int phred_thres,
            const double ess_georatio_dedup, 
            const double ess_georatio_duped_pcr,
            bool is_loginfo_enabled, 
            unsigned int thread_id, 
            unsigned int nogap_phred,
            unsigned int highqual_thres_snv, 
            unsigned int highqual_thres_indel, 
            unsigned int uni_bias_r_max,
            const bool is_proton, 
            const AssayType assay_type,
            const unsigned int phred_indel_error_before_barcode_labeling,
            const unsigned int baq_per_aligned_base,
            const unsigned int sidereg_nbases,
            uint32_t bias_flag_snv,
            uint32_t bias_flag_indel,
            const auto & region_repeatvec,
            unsigned int specialflag = 0) {
        
        unsigned int niters = 0;
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            niters++;
            bool log_alns2 = (is_loginfo_enabled && ispowerof2(niters));
            if (log_alns2) {
                //LOG(logINFO) << "Thread " << thread_id << " is processing the " << niters << "/" << alns3.size() << " th duplex-aware families which has " << alns2pair[0].size() << " and " << alns2pair[1].size() << " members.\n" ;
            }
            assert (alns2pair[0].size() != 0 || alns2pair[1].size() != 0);
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                if (alns2.size() == 0) { continue; }
                uint32_t tid2, beg2, end2;
                fillTidBegEndFromAlns2(tid2, beg2, end2, alns2);
                Symbol2CountCoverage read_family_con_ampl(tid2, beg2, end2); 
                Symbol2CountCoverage read_family_amplicon(tid2, beg2, end2); 
                if (log_alns2) {
                    // LOG(logINFO) << "    has " << alns2.size() << " fragments on strand " << strand;
                }
                for (const auto & alns1 : alns2) {
                    uint32_t tid1, beg1, end1;
                    fillTidBegEndFromAlns1(tid1, beg1, end1, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid1, beg1, end1);
                    read_ampBQerr_fragWithR1R2.template updateByRead1Aln<BASE_QUALITY_MAX, false>(
                            alns1, 
                            frag_indel_ext, 
                            symbolType2addPhred, 
                            alns2.size(), 
                            frag_indel_basemax, 
                            alns2pair2dflag.second, 
                            nogap_phred, 
                            is_proton, 
                            region_symbolvec, 
                            this->dedup_ampDistr.at(strand).getIncluBegPosition(), 
                            this->bq_dirs_count, 
                            this->bq_regs_count,
                            this->bq_nms_count,
                            this->bq_baq_sum,
                            this->bq_dirs_bqsum,
                            sidereg_nbases,
                            region_repeatvec,
                            baq_per_aligned_base,
                            0);
                    read_family_con_ampl.template updateByConsensus<SYMBOL_COUNT_SUM, true>(read_ampBQerr_fragWithR1R2);
                    read_family_amplicon.updateByFiltering(read_ampBQerr_fragWithR1R2, this->bq_pass_thres[strand], 1, true, strand);
                    if (log_alns2) {
                        // LOG(logINFO) << "        has " << alns1.size() << " sequenced read templates.";
                        // LOG(logDEBUG) << "num-updates = " << updateresult;
                    }
                }
                for (size_t epos = read_family_amplicon.getIncluBegPosition(); epos < read_family_amplicon.getExcluEndPosition(); epos++) {
                    const auto & con_ampl_symbol2count = read_family_con_ampl.getByPos(epos);
                    for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+(unsigned int)symbolType)) {
                        AlignmentSymbol con_symbol; // = END_ALIGNMENT_SYMBOLS;
                        unsigned int con_count, tot_count;
                        con_ampl_symbol2count.fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        if (0 == tot_count) { continue ; }
                        this->fam_total_dep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1); // for genomic region, so add DP (every type of  DP is unfiltered)
                        if (1 == tot_count) {
                            this->fam_size1_dep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        } else if ((con_count * 5 < tot_count * 4)) { 
                            this->fam_nocon_dep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        }
                    }
                    const auto & amplicon_symbol2count = read_family_amplicon.getByPos(epos);
                    for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+(unsigned int)symbolType)) {
                        AlignmentSymbol con_symbol; // = END_ALIGNMENT_SYMBOLS;
                        unsigned int con_count, tot_count;
                        amplicon_symbol2count.fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        // This is one round of bootstrapping for EM. We can use a full EM algorithm to improve this estimator, but it is probably overkill
                        // 0 means no coverage, 1 means no error correction, 2 means low quality family if the symbols disagree with each other, 
                        // 3 means up to 1 erroneous basecall is tolerated, 4 means up to 2 erroneous basecalls are tolerated
                        if (tot_count < 4) { continue; } 
                        if ((con_count * 3 < tot_count * 2)) { continue; }
                        for (AlignmentSymbol 
                                symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; 
                                symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; 
                                symbol = AlignmentSymbol(1+(unsigned int)symbol)) {
                            if (con_symbol != symbol) {
                                this->minor_amplicon[strand].getRefByPos(epos).incSymbolCount(symbol, amplicon_symbol2count.getSymbolCount(symbol));
                                this->major_amplicon[strand].getRefByPos(epos).incSymbolCount(symbol, tot_count);
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
                for (const auto & aln_vec : alns2) {
                    uint32_t tid1, beg1, end1;
                    fillTidBegEndFromAlns1(tid1, beg1, end1, aln_vec);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid1, beg1, end1);
                    read_ampBQerr_fragWithR1R2.template updateByRead1Aln<BASE_QUALITY_MAX, false>(
                            aln_vec, 
                            frag_indel_ext, 
                            symbolType2addPhred, 
                            alns2.size(), 
                            frag_indel_basemax, 
                            alns2pair2dflag.second, 
                            nogap_phred, 
                            is_proton, 
                            region_symbolvec, 
                            this->dedup_ampDistr.at(strand).getIncluBegPosition(), 
                            this->bq_dirs_count, 
                            this->bq_regs_count,
                            this->bq_nms_count,
                            this->bq_baq_sum,
                            this->bq_dirs_bqsum,
                            sidereg_nbases,
                            region_repeatvec,
                            baq_per_aligned_base,
                            0);
                    // The line below is similar to : read_family_amplicon.updateByConsensus<SYMBOL_COUNT_SUM>(read_ampBQerr_fragWithR1R2);
                    read_family_amplicon.updateByFiltering(read_ampBQerr_fragWithR1R2, this->bq_pass_thres[strand], 1, true, strand);
                    if (log_alns2) {
                        // LOG(logDEBUG) << "num-updates = " << updateresult << " (second time)";
                    }
                }
                if ((2 == alns2pair2dflag.second) && alns2pair[0].size() > 0 && alns2pair[1].size() > 0) { // is duplex
                    read_duplex_amplicon.template updateByConsensus<SYMBOL_COUNT_SUM>(read_family_amplicon);
                }
                std::basic_string<std::pair<unsigned int, AlignmentSymbol>> pos_symbol_string;
                unsigned int ldist_inc = 0;
                unsigned int rdist_inc = 0;
                std::vector<unsigned int> posToInsertLen = read_family_amplicon.computeZeroBasedPosToInsLenVec(rdist_inc); // extend
                unsigned int n_vars = 0;
                unsigned int n_seq_bases = 0;
                std::vector<std::array<AlignmentSymbol, NUM_SYMBOL_TYPES>> con_symbols_vec(
                        read_family_amplicon.getExcluEndPosition() - read_family_amplicon.getIncluBegPosition(),
                        std::array<AlignmentSymbol, NUM_SYMBOL_TYPES>({END_ALIGNMENT_SYMBOLS, END_ALIGNMENT_SYMBOLS}));
                for (size_t epos = read_family_amplicon.getIncluBegPosition(); epos < read_family_amplicon.getExcluEndPosition(); epos++) {
                    unsigned int ldist = 1 + epos - read_family_amplicon.getIncluBegPosition();
                    unsigned int rdist = read_family_amplicon.getExcluEndPosition() - epos;
                    for (SymbolType symbolType : SYMBOL_TYPES_IN_VCF_ORDER) {
                        AlignmentSymbol con_symbol; // = END_ALIGNMENT_SYMBOLS;
                        unsigned int con_count, tot_count;
                        read_family_amplicon.getRefByPos(epos).fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        if (0 == tot_count) { continue; }
                        
                        unsigned int majorcount = this->major_amplicon[strand].getByPos(epos).getSymbolCount(con_symbol);
                        unsigned int minorcount = this->minor_amplicon[strand].getByPos(epos).getSymbolCount(con_symbol);
                        auto con_bq_pass_thres = this->bq_pass_thres[strand].getByPos(epos).getSymbolCount(con_symbol);
                        double con_bq_pass_prob = phred2prob(con_bq_pass_thres) * (1 - DBL_EPSILON); // < 1
                        assert (con_bq_pass_prob >= pow(10, (-(double)NUM_BUCKETS)/10) 
                                || !fprintf(stderr, "%lf >= phred51 failed at position %lu and symbol %u!\n", con_bq_pass_prob, epos, con_symbol));
                        double prior_weight = 1.0 / (minorcount + 1.0);
                        double phredlike_db = h01_to_phredlike<true>(minorcount + prior_weight, majorcount + minorcount + prior_weight / con_bq_pass_prob, 
                                con_count * 2 - MIN(con_count * 2, tot_count), tot_count, 1.0, (ess_georatio_duped_pcr));
                        unsigned int phredlike = (unsigned int)MAX(0.0, phredlike_db);
                        if (BASE_N == con_symbol) { phredlike = MIN(phredlike, phred_thres); }
                        phredlike = MIN(phredlike, NUM_BUCKETS - SIGN2UNSIGN(1));
                        if (LINK_SYMBOL == symbolType) { 
                            unsigned int lim_phred = prob2phred((double)(minorcount + prior_weight) 
                                    / (double)(majorcount + minorcount + prior_weight / con_bq_pass_prob));
                            phredlike = MIN(phredlike, lim_phred + phred_indel_error_before_barcode_labeling); 
                        } // assuming errors equivalent to 200 pre-UMI-attachment PCR cycles after UMI attachment and one such PCR cycle before UMI attachment. QUESTION-TODO: justify
                        // no base quality stuff
                        
                        con_symbols_vec[epos - read_family_amplicon.getIncluBegPosition()][symbolType] = con_symbol;
                        this->fq_tsum_depth [strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        // unsigned int edge_baq = MIN(ldist, rdist) * baq_per_aligned_base;
                        unsigned int overallq = phredlike; // MIN(edge_baq, phredlike);
                        this->fq_qual_p1sum [strand].getRefByPos(epos).incSymbolCount(con_symbol, overallq);
                        // this->fq_qual_p2sum [strand].getRefByPos(epos).incSymbolCount(con_symbol, mathsquare(overallq));
                        if (overallq >= ((BASE_SYMBOL == symbolType) ? highqual_thres_snv : ((LINK_SYMBOL == symbolType) ? highqual_thres_indel : 0))) {
                            this->fq_hiqual_dep[strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        }
                        unsigned int pbucket = phred2bucket(overallq);
                        assert (pbucket < NUM_BUCKETS || !fprintf(stderr, "%u < %u failed at position %lu and con_symbol %u symboltype %u plusbucket %u\n", 
                                 pbucket,  NUM_BUCKETS, epos, con_symbol, symbolType, symbolType2addPhred[symbolType]));
                        if (isSymbolIns(con_symbol)) {
                            posToIndelToCount_updateByConsensus(this->fq_tsum_depth[strand].getRefPosToIseqToData(con_symbol), read_family_amplicon.getPosToIseqToData(con_symbol), epos, 1);
                        }
                        if (isSymbolDel(con_symbol)) {
                            posToIndelToCount_updateByConsensus(this->fq_tsum_depth[strand].getRefPosToDlenToData(con_symbol), read_family_amplicon.getPosToDlenToData(con_symbol), epos, 1);
                        }
                        AlignmentSymbol refsymbol = region_symbolvec[epos - this->dedup_ampDistr.at(strand).getIncluBegPosition()]; 
                        if (areSymbolsMutated(refsymbol, con_symbol)) {
                            pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
                            if (symbolType == BASE_SYMBOL && phredlike >= phred_thres) {
                                n_vars++;
                            }
                        }
                        if (symbolType == BASE_SYMBOL) {
                            n_seq_bases++;
                        }
                        this->dedup_ampDistr[strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pbucket, 1);
                        ldist_inc += posToInsertLen[epos - read_family_amplicon.getIncluBegPosition()];
                        this->pb_dist_lpart [strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pos2edbuck(ldist + ldist_inc), 1);
                        this->pb_dist_rpart [strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pos2edbuck(rdist + rdist_inc), 1);
                        this->fq_bsum_ldist [strand].getRefByPos(epos).incSymbolCount(con_symbol, ldist + ldist_inc);
                        this->fq_bsum_rdist [strand].getRefByPos(epos).incSymbolCount(con_symbol, rdist + rdist_inc);
                        rdist_inc -= posToInsertLen[epos - read_family_amplicon.getIncluBegPosition()];
                    }
                }
                n_vars = MIN(n_vars, NUM_NMBUCKS - 1);
                for (auto epos = read_family_amplicon.getIncluBegPosition(); epos < read_family_amplicon.getExcluEndPosition(); epos++) {
                    for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+((unsigned int)symbolType))) {
                        auto con_symbol = con_symbols_vec.at(epos - read_family_amplicon.getIncluBegPosition()).at(symbolType);
                        if (END_ALIGNMENT_SYMBOLS != con_symbol) {
                            this->pb_dist_nvars[strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, n_vars, 1);
                            this->fq_n_seq_bases[strand].getRefByPos(epos).incSymbolCount(con_symbol, n_seq_bases);
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
                            this->duplex_tsum_depth.getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        }
                        if (1 < tot_count) {
                            this->duplex_pass_depth.getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        }
                    }
                }
            }
        }
        adafilter<true>(
                this->du_bias_dedup,
                this->bq_qual_p1sum, this->bq_tsum_depth, //  this->fq_qual_p2sum, // this->fq_imba_depth,
                this->fq_amax_ldist, this->fq_amax_rdist, this->fq_bias_ldist, this->fq_bias_rdist, 
                this->fq_amax_nvars, this->fq_bias_nvars,
                this->fq_bsum_ldist, this->fq_bsum_rdist, this->fq_bias_1stra, this->fq_bias_2stra,
                this->fq_tsum_depth, 
                this->fq_pass_thres, this->fq_pass_depth,
                this->fq_vars_thres, this->fq_vars_depth, this->fq_vars_badep, this->fq_vars_vqual,
                this->dedup_ampDistr,this->bq_tsum_depth,
                this->pb_dist_lpart, this->pb_dist_rpart, this->pb_dist_nvars, this->additional_note, should_add_note, 
                this->bq_dirs_count, this->bq_bias_sedir, this->fq_n_seq_bases, assay_type,
                phred_max_table, symbolType2addPhred, ess_georatio_dedup, uni_bias_r_max, bias_flag_snv, bias_flag_indel);
        return 0;
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
    
    int 
    updateHapMap(std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map, 
            const auto & tsum_depth, unsigned int max_ploidy = 2+1) {
        for (auto it = mutform2count4map.begin(); it != mutform2count4map.end();) {
            std::basic_string<std::pair<unsigned int, AlignmentSymbol>> mutform = it->first;
            auto counts = it->second;
            std::vector<unsigned int> dsADs; 
            dsADs.reserve(mutform.size() + 1);
            dsADs.push_back(0);
            for (std::pair<unsigned int, AlignmentSymbol> simplemut : mutform) {
                unsigned int ad0 = tsum_depth.at(0).getByPos(simplemut.first).getSymbolCount(simplemut.second);
                unsigned int ad1 = tsum_depth.at(1).getByPos(simplemut.first).getSymbolCount(simplemut.second);
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
            const unsigned int frag_indel_ext, 
            const unsigned int frag_indel_basemax,
            const PhredMutationTable & phred_max_sscs_table, 
            unsigned int phred_thres,
            const double ess_georatio_dedup, 
            const double ess_georatio_duped_pcr,
            bool use_deduplicated_reads, 
            bool is_loginfo_enabled, 
            unsigned int thread_id, 
            unsigned int fixedthresBQ, 
            unsigned int nogap_phred,
            unsigned int highqual_thres_snv, 
            unsigned int highqual_thres_indel, 
            unsigned int uni_bias_r_max,
            const bool is_proton, 
            const AssayType assay_type, // unused?
            const unsigned int phred_indel_error_before_barcode_labeling,
            const unsigned int baq_per_aligned_base,
            const unsigned int regside_nbases,
            const uint32_t bias_flag_snv,
            const uint32_t bias_flag_indel,
            const unsigned int specialflag = 0) {
        
        const std::array<unsigned int, NUM_SYMBOL_TYPES> symbolType2addPhred = {{bq_phred_added_misma, bq_phred_added_indel}};
        std::basic_string<AlignmentSymbol> ref_symbol_string = string2symbolseq(refstring);
        // std::vector<RegionalTandemRepeat> region_repeatvec = refstring2repeatvec(refstring);

        updateByAlns3UsingBQ(mutform2count4map_bq, 
                alns3, 
                ref_symbol_string, 
                symbolType2addPhred, 
                should_add_note, 
                frag_indel_ext, 
                frag_indel_basemax, 
                phred_max_sscs_table, 
                phred_thres, 
                ess_georatio_dedup, 
                ess_georatio_duped_pcr, 
                fixedthresBQ, 
                nogap_phred, 
                uni_bias_r_max, is_proton, 
                assay_type,
                baq_per_aligned_base,
                regside_nbases,
                bias_flag_snv,
                bias_flag_indel,
                region_repeatvec, 
                0); // base qualities
        updateHapMap(mutform2count4map_bq, this->bq_tsum_depth);
        if (use_deduplicated_reads) {
            updateByAlns3UsingFQ(mutform2count4map_fq, 
                    alns3,
                    ref_symbol_string, 
                    symbolType2addPhred, 
                    should_add_note, 
                    frag_indel_ext, 
                    frag_indel_basemax, 
                    phred_max_sscs_table, 
                    phred_thres,
                    ess_georatio_dedup, 
                    ess_georatio_duped_pcr,
                    is_loginfo_enabled, 
                    thread_id, 
                    nogap_phred,
                    highqual_thres_snv, 
                    highqual_thres_indel, 
                    uni_bias_r_max,
                    is_proton, 
                    assay_type,
                    phred_indel_error_before_barcode_labeling,
                    baq_per_aligned_base,
                    regside_nbases,
                    bias_flag_snv,
                    bias_flag_indel,
                    region_repeatvec,
                    0); // family qualities
            updateHapMap(mutform2count4map_fq, this->fq_tsum_depth);
        }
        return 0;
    };
};

std::array<unsigned int, 2>
BcfFormat_init(bcfrec::BcfFormat & fmt, 
        const Symbol2CountCoverageSet & symbolDistrSets12, 
        unsigned int refpos, 
        SymbolType symbolType, 
        bool use_deduplicated_reads, 
        const AlignmentSymbol refsymbol) {
    
    for (unsigned int strand = 0; strand < 2; strand++) {
        fmt.bAllBQ[strand] = (symbolDistrSets12.bq_qual_p1sum.at(strand).getByPos(refpos).sumBySymbolType2(symbolType)); 
        // fmt.bAllBQ[strand] = div_by_20(symbolDistrSets12.bq_qual_p2sum.at(strand).getByPos(refpos).sumBySymbolType2(symbolType)); 

        if (use_deduplicated_reads) {
            fmt.cAllBQ[strand] = (symbolDistrSets12.fq_qual_p1sum.at(strand).getByPos(refpos).sumBySymbolType2(symbolType)); 
            //fmt.cAllBQ2[strand] = div_by_20(symbolDistrSets12.fq_qual_p2sum.at(strand).getByPos(refpos).sumBySymbolType2(symbolType)); 
        } else {
            fmt.cAllBQ[strand] = fmt.bAllBQ[strand]; 
            //fmt.cAllBQ2[strand] = div_by_20(symbolDistrSets12.bq_qual_p2sum.at(strand).getByPos(refpos).sumBySymbolType2(symbolType)); 
        }
        fmt.cAllHD[strand] = symbolDistrSets12.fq_hiqual_dep.at(strand).getByPos(refpos).sumBySymbolType2(symbolType);
        
        fmt.bRefBQ[strand] = (symbolDistrSets12.bq_qual_p1sum.at(strand).getByPos(refpos).getSymbolCount(refsymbol)); 
        // fmt.bRefBQ[strand] = div_by_20(symbolDistrSets12.bq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(refsymbol)); 
        if (use_deduplicated_reads) {
            fmt.cRefBQ[strand] = (symbolDistrSets12.fq_qual_p1sum.at(strand).getByPos(refpos).getSymbolCount(refsymbol)); 
            //fmt.cRefBQ2[strand] = div_by_20(symbolDistrSets12.fq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(refsymbol)); 
        } else {
            fmt.cRefBQ[strand] = fmt.bRefBQ[strand]; 
            //fmt.cRefBQ2[strand] = div_by_20(symbolDistrSets12.bq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(refsymbol)); 
        }
        fmt.cRefHD[strand] = symbolDistrSets12.fq_hiqual_dep.at(strand).getByPos(refpos).getSymbolCount(refsymbol); 
        
        fmt.bDP1[strand] = symbolDistrSets12.bq_tsum_depth.at(strand).getByPos(refpos).sumBySymbolType2(symbolType);
        fmt.cDP1[strand] = symbolDistrSets12.fq_tsum_depth.at(strand).getByPos(refpos).sumBySymbolType2(symbolType);
        fmt.bDPLQ[strand] = symbolDistrSets12.bq_tsum_LQdep.at(strand).getByPos(refpos).sumBySymbolType2(symbolType);
       
        fmt.cDPTT[strand] = symbolDistrSets12.fam_total_dep.at(strand).getByPos(refpos).sumBySymbolType2(symbolType);
        fmt.cDPT1[strand] = symbolDistrSets12.fam_size1_dep.at(strand).getByPos(refpos).sumBySymbolType2(symbolType);
        auto fmt_cDPTN = symbolDistrSets12.fam_nocon_dep.at(strand).getByPos(refpos).sumBySymbolType2(symbolType);
        fmt.cDPTC[strand] = fmt.cDPTT[strand] - fmt.cDPT1[strand] - fmt_cDPTN;
        
        fmt.bRD1[strand] = symbolDistrSets12.bq_tsum_depth.at(strand).getByPos(refpos).getSymbolCount(refsymbol);
        fmt.cRD1[strand] = symbolDistrSets12.fq_tsum_depth.at(strand).getByPos(refpos).getSymbolCount(refsymbol);
        fmt.bRDLQ[strand] = symbolDistrSets12.bq_tsum_LQdep.at(strand).getByPos(refpos).getSymbolCount(refsymbol);

        fmt.cRDTT[strand] = symbolDistrSets12.fam_total_dep.at(strand).getByPos(refpos).getSymbolCount(refsymbol);
        fmt.cRDT1[strand] = symbolDistrSets12.fam_size1_dep.at(strand).getByPos(refpos).getSymbolCount(refsymbol);
        auto fmt_cRDTN = symbolDistrSets12.fam_nocon_dep.at(strand).getByPos(refpos).getSymbolCount(refsymbol); 
        fmt.cRDTC[strand] = fmt.cRDTT[strand] - fmt.cRDT1[strand] - fmt_cRDTN;
        
        fmt.aDP[strand*2+0] = symbolDistrSets12.bq_dirs_count.at(strand*2+0).getByPos(refpos).sumBySymbolType2(symbolType);
        fmt.aDP[strand*2+1] = symbolDistrSets12.bq_dirs_count.at(strand*2+1).getByPos(refpos).sumBySymbolType2(symbolType);
        
        fmt.aRD[strand*2+0] = symbolDistrSets12.bq_dirs_count.at(strand*2+0).getByPos(refpos).getSymbolCount(refsymbol);
        fmt.aRD[strand*2+1] = symbolDistrSets12.bq_dirs_count.at(strand*2+1).getByPos(refpos).getSymbolCount(refsymbol);
        
        fmt.aNMRD[strand*2+0] = symbolDistrSets12.bq_nms_count.at(strand*2+0).getByPos(refpos).getSymbolCount(refsymbol);
        fmt.aNMRD[strand*2+1] = symbolDistrSets12.bq_nms_count.at(strand*2+1).getByPos(refpos).getSymbolCount(refsymbol);
        
        // fmt.bEDRD[strand] = symbolDistrSets12.bq_n_edit_ops.at(strand).getByPos(refpos).getSymbolCount(refsymbol);
    }
    fmt.aBAQDP = symbolDistrSets12.bq_baq_sum.at(0).getByPos(refpos).sumBySymbolType2(symbolType);
    fmt.aBAQADR = {{symbolDistrSets12.bq_baq_sum.at(0).getByPos(refpos).getSymbolCount(refsymbol), 0}};

    assert(fmt.aPBDP.size() == symbolDistrSets12.bq_regs_count.size());
    for (unsigned int region = 0; region < symbolDistrSets12.bq_regs_count.size(); region++) {
        fmt.aPBDP[region] = symbolDistrSets12.bq_regs_count.at(region).getByPos(refpos).sumBySymbolType2(symbolType);
    }
    
    fmt.gapSeq.clear();
    fmt.gapbAD1.clear();
    fmt.gapcAD1.clear();
    fmt.gapNum[0] = 0;
    fmt.gapNum[1] = 0;

    fmt.dDP1 = symbolDistrSets12.duplex_tsum_depth.getByPos(refpos).sumBySymbolType2(symbolType);
    return {fmt.bDP1[0] + fmt.bDP1[1], fmt.cDPTT[0] + fmt.cDPTT[1]};
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
                symbol2CountCoverageSet.bq_tsum_depth.at(strand).getPosToIseqToData(symbol),
                symbol2CountCoverageSet.fq_tsum_depth.at(strand).getPosToIseqToData(symbol),
                refstring, repeatunit, repeatnum);
    } else {
        return fill_by_indel_info2_2(fmt, symbol2CountCoverageSet, strand, refpos, symbol,
                symbol2CountCoverageSet.bq_tsum_depth.at(strand).getPosToDlenToData(symbol),
                symbol2CountCoverageSet.fq_tsum_depth.at(strand).getPosToDlenToData(symbol),
                refstring, repeatunit, repeatnum);
    }
};

double 
qthres_ad_dp_to_qtotal(unsigned int QP, unsigned int ADP, unsigned int DPT,
        double positive_pseudocount = 1, double negative_pseudocount = 1) {
    auto observed_unit_phred = log((double)ADP / (double)DPT) * (10 / log(10));
    return (observed_unit_phred - QP) * log(ADP+positive_pseudocount) / log(positive_pseudocount + negative_pseudocount);
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

const std::vector<std::pair<unsigned int, std::string>> 
indel_get_majority(
        const bcfrec::BcfFormat & fmt, 
        // const bool prev_is_tumor, const auto & tki, bool is_rescued, 
        const char *tname, 
        unsigned int refpos, 
        const AlignmentSymbol symbol, 
        bool is_warning_generated = true) {
    // std::string indelstring = "";
    std::vector<std::pair<unsigned int, std::string>> indelstrings;
    /*
    if (prev_is_tumor && tki.ref_alt.size() > 3) {
        size_t midtry = tki.ref_alt.size() - 2;
        if ('\t' == tki.ref_alt[1] && tki.ref_alt[0] == tki.ref_alt[2]) {
            indelstring = tki.ref_alt.substr(3, tki.ref_alt.size() - 3);
        } else if ('\t' == tki.ref_alt[midtry] && tki.ref_alt[0] == tki.ref_alt[midtry+1]) {
            indelstring = tki.ref_alt.substr(1, tki.ref_alt.size() - 3);
        } else {
            if (is_warning_generated) {
                fprintf(stderr, "ThisIsAVerySevereWarning: the ref_alt %s is not valid!\n", tki.ref_alt.c_str());
                std::string msg;
                bcfrec::streamAppendBcfFormat(msg, fmt);
                std::cerr << msg << "\n";
                // abort();
            }
            indelstring = SYMBOL_TO_DESC_ARR[symbol];
        }
    } else */
    if (fmt.gapNum[0] <= 0 && fmt.gapNum[1] <= 0) {
        //if (!is_rescued) {
            if (is_warning_generated) {
                std::cerr << "Invalid indel detected (invalid mutation) : " << tname << ", " << refpos << ", " << SYMBOL_TO_DESC_ARR[symbol] << std::endl;
                std::string msg;
                bcfrec::streamAppendBcfFormat(msg, fmt);
                std::cerr << msg << "\n";
                // assert(false);
            }
            indelstrings.push_back(std::make_pair(0, SYMBOL_TO_DESC_ARR[symbol]));
        // }
    } else {
        std::map<std::string, unsigned int> indelmap;
        assert(fmt.gapSeq.size() ==  fmt.gapNum[0] + fmt.gapNum[1]);
        assert(fmt.gapbAD1.size() ==  fmt.gapNum[0] + fmt.gapNum[1]);
        for (unsigned int i = 0; i < fmt.gapNum[0] + fmt.gapNum[1]; i++) {
            const auto & it = indelmap.find(fmt.gapSeq[i]);
            if (it == indelmap.end()) {
                indelmap.insert(std::make_pair(fmt.gapSeq[i], fmt.gapbAD1[i]));
            } else {
                indelmap[fmt.gapSeq[i]] = indelmap[fmt.gapSeq[i]] + fmt.gapbAD1[i];
            }
        }
        unsigned int max_bAD1 = 0;
        for (auto & indelit : indelmap) {
            UPDATE_MAX(max_bAD1, indelit.second);
        }
        for (auto & indelit : indelmap) {
            if (indelit.second >= (max_bAD1 + 3) / 4) {
                indelstrings.push_back(std::make_pair(indelit.second, indelit.first));
            }
        }
        std::sort(indelstrings.rbegin(), indelstrings.rend());
    }
    /*
    if (fmt.gapNum[0] <= 0 || fmt.gapNum[1] <= 0) {
        indelstring = fmt.gapSeq[0];
    } else {
        indelstring = (fmt.gapbAD1[0] > fmt.gapbAD1[fmt.gapNum[0]] ? fmt.gapSeq[0] : fmt.gapSeq.at(fmt.gapNum[0]));
    }
    */
    return indelstrings;
}

template <class T = int>
struct RepNumCluster {
    T mode;
    unsigned int cnt0;  // centroid count
    unsigned int cnt1m; // 1st left count
    unsigned int cnt1p; // 1st right count
    unsigned int cnt2m; // 2nd left count
    unsigned int cnt2p; // 2nd right count
};

// 
int
indel_fill_rep_num_clusters(std::array<RepNumCluster<int>, RCC_NUM> & rep_num_clusters, 
        const unsigned int refcnt,
        const std::string & repeatunit,
        const std::vector<std::map<std::string, uint32_t>> & iseq2cnt_vec,
        const std::vector<std::map<uint32_t   , uint32_t>> & dlen2cnt_vec) {
    
    for (size_t i = 0; i < rep_num_clusters.size(); i++) {
        rep_num_clusters[i].mode = 0;
        rep_num_clusters[i].cnt0 = 0;
        rep_num_clusters[i].cnt1m = 0;
        rep_num_clusters[i].cnt1p = 0;
        rep_num_clusters[i].cnt2m = 0;
        rep_num_clusters[i].cnt2p = 0;

    }
    
    unsigned int max_ilen = 0;
    for (auto & iseq2cnt : iseq2cnt_vec) { 
        for (auto & iseq_cnt : iseq2cnt) {
            max_ilen = MAX(max_ilen, iseq_cnt.first.size() / repeatunit.size());
        }
    }
    unsigned int max_dlen = 0;
    for (auto & dlen2cnt : dlen2cnt_vec) { 
        for (auto & dlen_cnt : dlen2cnt) {
            max_dlen = MAX(max_dlen, dlen_cnt.first        / repeatunit.size());
        }
    }
    
    // at index max_dlen is the ref non-indel // ... (STR-2,STR-1] (STR-1,REF) REF (REF,STR+1) [STR+1,STR+2) ...
    std::vector<unsigned int> idx2cnt(max_dlen + max_ilen + 1, 0);
    for (auto & iseq2cnt : iseq2cnt_vec) {
        for (auto & iseq_cnt : iseq2cnt) {
            auto idx = max_dlen + iseq_cnt.first.size() / repeatunit.size();
            idx2cnt[idx] += iseq_cnt.second;
        }
    }
    for (auto & dlen2cnt : dlen2cnt_vec) {
        for (auto & dlen_cnt : dlen2cnt) {
            auto idx = max_dlen - dlen_cnt.first       / repeatunit.size();
            idx2cnt[idx] += dlen_cnt.second;
        }
    }
    idx2cnt[max_dlen] += refcnt;
    
    std::vector<std::pair<unsigned int, unsigned int>> cnt_len_vec; 
    for (size_t idx = 0; idx < idx2cnt.size(); idx++) {
        if (idx2cnt[idx] > 0 || max_dlen == idx) {
            cnt_len_vec.push_back(std::make_pair(idx2cnt[idx], idx));
        }
    }
    std::sort(cnt_len_vec.rbegin(), cnt_len_vec.rend());
    
    for (size_t i = 0; i < rep_num_clusters.size(); i++) {
        size_t idx = cnt_len_vec[MIN(i, cnt_len_vec.size()-1)].second;
        rep_num_clusters[i].mode  = (int)idx - (int)max_dlen;
        rep_num_clusters[i].cnt0  = idx2cnt[idx]; // (int)cnt_len_vec[0].first;
        rep_num_clusters[i].cnt1m = ((                  idx<1) ? 0 : idx2cnt[idx-1]);
        rep_num_clusters[i].cnt1p = ((idx2cnt.size() <= idx+1) ? 0 : idx2cnt[idx+1]);
        rep_num_clusters[i].cnt2m = ((                  idx<2) ? 0 : idx2cnt[idx-2]);
        rep_num_clusters[i].cnt2p = ((idx2cnt.size() <= idx+2) ? 0 : idx2cnt[idx+2]);
    }

    return 0;
}

double
penal_indel_2(double AD0a, int dst_str_units, const auto & RCC, const unsigned int phred_triallelic_indel, const double snr_data2model = 2.0, const double snr_model2data = 4.0) {
    double max_noise = 0.0;
    for (int c = 0; c < RCC_NUM; c++) {
        int peakidx = c*RCC_NFS;
        double src_str_units = RCC[peakidx];
        if (dst_str_units == src_str_units) { continue; }
        //const double peak_height1 = (double)RCC[peakidx+3];
        const double deldev2 = (double)RCC[peakidx+1] * (snr_data2model * snr_data2model); // peak_height1;
        const double deldev1 = (double)RCC[peakidx+2] * (snr_data2model); // peak_height1;
        const double insdev1 = (double)RCC[peakidx+4] * (snr_data2model); // peak_height1;
        const double insdev2 = (double)RCC[peakidx+5] * (snr_data2model * snr_data2model); // peak_height1;
        double _noise = 1.0;
        if        (src_str_units + 2 == dst_str_units) { // ins of two units wrt to peak
            _noise = MAX3(deldev2 / 2.0, deldev1 / 4.0, insdev1 / 2.0                ) / (snr_model2data * snr_model2data);
        } else if (src_str_units + 1 == dst_str_units) {
            _noise = MAX3(deldev2 / 2.0, deldev1 / 4.0,                insdev2 / 2.0 ) / (snr_model2data);
        } else if (src_str_units - 1 == dst_str_units) {
            _noise = MAX3(deldev2      ,                insdev1 * 2.0, insdev2       ) / (snr_model2data);
        } else if (src_str_units - 2 == dst_str_units) {
            _noise = MAX3(               deldev1,       insdev1 * 2.0, insdev2       ) / (snr_model2data * snr_model2data);
        } else if (src_str_units      < dst_str_units) { // is other types of ins
            _noise = MAX4(deldev2 / 2.0, deldev1 / 4.0, insdev1 / 2.0, insdev2 / 2.0 ) / pow(snr_model2data, dst_str_units - src_str_units);
        } else {
            _noise = MAX4(deldev2,       deldev1,       insdev1 * 2.0, insdev2       ) / pow(snr_model2data, src_str_units - dst_str_units);
        }
        max_noise = MAX(max_noise, _noise);
    }
    return ((double)(phred_triallelic_indel)) / log(2.0) * MIN(2.0 * log(2.0), log((AD0a + DBL_MIN + max_noise) / (AD0a + DBL_MIN)));
}

// higher allele1count (lower allele2count) results in higher LODQ if FA is above frac, meaning allele1 is more likely to be homo, and vice versa
unsigned int hetLODQ(double allele1count, double allele2count, double frac, double powlaw_exponent = 3.0) {
    auto binomLODQ = (int)calc_binom_10log10_likeratio(frac, allele1count, allele2count);
    auto powerLODQ = (int)(10.0/log(10.00) * powlaw_exponent * MAX(logit2(allele1count / (frac * 2.0) + 0.5, allele2count / ((1.0 - frac) * 2.0) + 0.5), 0.0));
    return MIN(binomLODQ, powerLODQ);
}

struct {
    template<class T1, class T2>
    bool operator()(std::pair<T1, T2> a, std::pair<T1, T2> b) const {   
        return (a.second < b.second) || (a.second == b.second && a.first < b.first);
    }
} PairSecondLess;

double compute_norm_ad(const bcfrec::BcfFormat *fmtp, const bool isSubst) {
    if (isSubst) {
        double fa_baq = fmtp->aBAQADR[1] * INDEL_MUL_PER_BAQ / ((double)(fmtp->aBAQDP + fmtp->aBAQADR[1] * (INDEL_MUL_PER_BAQ - 1)) + DBL_EPSILON);
        double fa_bq  = SUM2(fmtp->bAltBQ) / ((double)SUM2(fmtp->bAllBQ) + DBL_EPSILON);
        return MIN(MIN(fa_baq, fa_bq) * SUM2(fmtp->cDPTT), (double)SUM2(fmtp->cADTT));
    } else {
        double fa_baq = fmtp->aBAQADR[1] * SNV_MUL_PER_BAQ / ((double)(fmtp->aBAQDP + fmtp->aBAQADR[1] * (SNV_MUL_PER_BAQ - 1)) + DBL_EPSILON);
        return MIN(fa_baq * SUM2(fmtp->cDPTT), (double)SUM2(fmtp->cADTT));
    }
}

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
        unsigned int specialflag) {
    
    assert(symbol_format_vec.size() >= 4 || 
            !fprintf(stderr, " The variant-type %s:%u %u has symbol_format_vec of length %u", 
            tname, refpos, refsymbol, symbol_format_vec.size()));
    unsigned int regionpos = refpos- extended_inclu_beg_pos;
    struct {
        bool operator()(std::pair<AlignmentSymbol, bcfrec::BcfFormat*> & p1, std::pair<AlignmentSymbol, bcfrec::BcfFormat*> & p2) const {
            return p1.second->ALODQ < p2.second->ALODQ;
        }
    } SymbolBcfFormatPairLess;
    std::sort(symbol_format_vec.rbegin(), symbol_format_vec.rend(), SymbolBcfFormatPairLess);
    std::array<std::pair<AlignmentSymbol, bcfrec::BcfFormat*>, 4> ref_alt1_alt2_alt3 = {{std::make_pair(AlignmentSymbol(0), (bcfrec::BcfFormat*)NULL)}};
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
    
    int a0LODQA = ref_alt1_alt2_alt3[0].second->ALODQ;
    int a0LODQB = ref_alt1_alt2_alt3[0].second->BLODQ;
    int a0LODQ = MIN(a0LODQA, a0LODQB) + 1;
    int a1LODQ = (ref_alt1_alt2_alt3[1].second->ALODQ);
    int a2LODQ = (ref_alt1_alt2_alt3[2].second->ALODQ);
    int a3LODQ = (ref_alt1_alt2_alt3[3].second->ALODQ);
    
    unsigned int readlen = MAX(30U, central_readlen);
    const double alt1frac = (double)(readlen - MIN(readlen - 30, ref_alt1_alt2_alt3[1].second->RefBias / 2.0)) / (double)readlen / 2.0;
    const double alt2frac = (double)(readlen - MIN(readlen - 30, ref_alt1_alt2_alt3[2].second->RefBias / 2.0)) / (double)readlen / 2.0;
    
    auto fmtptr0 = ref_alt1_alt2_alt3[0].second;
    auto fmtptr1 = ref_alt1_alt2_alt3[1].second;
    auto fmtptr2 = ref_alt1_alt2_alt3[2].second;
    bool isSubst = isSymbolSubstitution(refsymbol);
    double ad0 = compute_norm_ad(fmtptr0, isSubst);
    double ad1 = compute_norm_ad(fmtptr1, isSubst);
    double ad2 = compute_norm_ad(fmtptr2, isSubst);
    int a0a1LODQ = hetLODQ(ad0, ad1, 1.0 - alt1frac);
    int a1a0LODQ = hetLODQ(ad1, ad0, alt1frac);
    int a1a2LODQ = hetLODQ(ad1, ad2, alt1frac / (alt1frac + alt2frac));
    int a2a1LODQ = hetLODQ(ad2, ad1, alt2frac / (alt1frac + alt2frac));

    std::array<std::string, 4> GTidx2GT {{
        "0/0",
        "0/1",
        "1/1",
        "1/2"
    }};
    
    int phred_homref = 0; // (isSubst ? 31 : 41);
    int phred_hetero = (isSubst ? 31 : 41);
    int phred_homalt = (isSubst ? 33 : 43);
    int phred_tri_al = (isSubst ? 55 : 49-2); // https://www.genetics.org/content/184/1/233 : triallelic-SNP-phred = 29*2-3
    // tri_al for InDels is lower than expected because indels were already penalized for tri-allelelity in their TLODQs
    const double qfrac = (isSubst ? 1.0 : 0.25);
    std::array<std::pair<int, int>, 4> GL4raw = {{
        std::make_pair(0,     (-phred_homref - a1LODQ              - MAX(a2LODQ - phred_tri_al, 0))),
        std::make_pair(1,  MIN(-phred_hetero - MAX(a0a1LODQ, a1a0LODQ), -a2LODQ)),
        std::make_pair(2,     (-phred_homalt - MAX(a0LODQ, a2LODQ) - MAX(MIN(a0LODQ, a2LODQ) - phred_tri_al, 0))),
        std::make_pair(3,  MIN(-phred_tri_al - MAX(a1a2LODQ, a2a1LODQ), -MAX(a0LODQ, a3LODQ)))
    }};
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
        std::vector<std::pair<unsigned int, std::string>> indelstrings1 = indel_get_majority(*(ref_alt1_alt2_alt3[1].second), 
            // false, tki, false, 
            tname, refpos, s1, false);
        const std::string & indelstring1 = indelstrings1[0].second;
        if (indelstrings1.size() > 1) {
            alt1_uniallelic_phred = (phred_tri_al - phred_hetero) *  log(1.0 + (double)indelstrings1[0].first / (double)indelstrings1[1].first) / log(2.0);
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
                std::vector<std::pair<unsigned int, std::string>> indelstrings2 = indel_get_majority(*(ref_alt1_alt2_alt3[2].second), 
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
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[0].second->cADR, 1, 0));
    germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[1].second->cADR, 1, 0));
    if (3 == GLidx) {
        germ_ADR.push_back(collectget(ref_alt1_alt2_alt3[2].second->cADR, 1, 0));
    }
    std::vector<std::string> germ_FT;
    germ_FT.push_back(ref_alt1_alt2_alt3[0].second->FT);
    germ_FT.push_back(ref_alt1_alt2_alt3[1].second->FT);
    if (3 == GLidx) {
        germ_FT.push_back(ref_alt1_alt2_alt3[2].second->FT);
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
                a0LODQA, a0LODQB, a1LODQ, a2LODQ, a3LODQ,
                a0a1LODQ, a1a0LODQ, a1a2LODQ, a2a1LODQ
            }}, ","),
            ref_alt1_alt2_alt3[0].second->note
        }}, ":")
    }}, "\t") + "\n";
    out_string += bcfline;
    return GL4raw[0].second - MAX3(GL4raw[1].second, GL4raw[2].second, GL4raw[3].second);
}

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

int
fill_by_symbol(bcfrec::BcfFormat & fmt,
        // const std::vector<std::pair<AlignmentSymbol, unsigned int>> & symb_index_vec,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet12, 
        unsigned int refpos, 
        const AlignmentSymbol symbol, 
        const std::string & refstring, 
        unsigned int refstring_offset, 
        const std::vector<std::pair<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>>> & mutform2count4vec_bq,
        std::set<size_t> indices_bq,
        const std::vector<std::pair<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>>> & mutform2count4vec_fq,
        std::set<size_t> indices_fq,
        unsigned int minABQ, // = 25
        unsigned int minMQ1, 
        unsigned int maxMQ,
        unsigned int phred_max_sscs,
        unsigned int phred_max_dscs,
        bool use_deduplicated_reads, 
        bool use_only_deduplicated_reads,
        bool is_rescued, 
        const bool prev_is_tumor, 
        const std::string & repeatunit, 
        const unsigned int repeatnum, 
        const auto & tki,
        const double any_mul_contam_frac, // addition params
        const double t2n_mul_contam_frac,
        const double t2n_add_contam_frac,
        const double t2n_add_contam_transfrac, // additional params
        unsigned int min_edge_dist,
        unsigned int central_readlen,
        unsigned int baq_per_aligned_base,
        double powlaw_exponent,
        const auto & bq_ins_2bdepths,
        const auto & bq_del_2bdepths,
        // ,
        // const auto & bq_indel_adjmax_depths
        const bool somaticGT,
        const bool is_ref_bias_aware,
        // const unsigned int phred_umi_dimret_qual, 
        const double phred_umi_dimret_mult,
        
        // const bool is_proton,
        const double illumina_BQ_sqr_coef,
        const unsigned int phred_varcall_err_per_map_err_per_base,
        const auto phred_snv_to_indel_ratio,
        const bool is_proton,
        
        const double powlaw_anyvar_base, // = 90.0;
        const bool tUseHD1, // = false;
        const double phred_triallelic_indel, // = 30.0;
        const double powlaw_dscs_inc,
        const double powlaw_sscs_inc,
        // const auto & indelstrings,
        // const auto & bad0a_indelstring_tkiidx_vec,
        const auto   indelbdepth,
        const auto & indelstring,
        const int specialflag) {
    
    fmt.note = symbol2CountCoverageSet12.additional_note.getByPos(refpos).at(symbol);
    uint64_t bq_qsum_sqrMQ_tot = 0; 
    for (unsigned int strand = 0; strand < 2; strand++) {
        
        fmt.bAltBQ[strand] = (symbol2CountCoverageSet12.bq_qual_p1sum.at(strand).getByPos(refpos).getSymbolCount(symbol));
        // fmt.bAltBQ[strand] = div_by_20(symbol2CountCoverageSet12.bq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(symbol));
        if (use_deduplicated_reads) {
            fmt.cAltBQ[strand] = (symbol2CountCoverageSet12.fq_qual_p1sum.at(strand).getByPos(refpos).getSymbolCount(symbol));
            //fmt.cAltBQ2[strand] = div_by_20(symbol2CountCoverageSet12.fq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(symbol));
        } else {
            fmt.cAltBQ[strand] = fmt.bAltBQ[strand];
            //fmt.cAltBQ2[strand] = div_by_20(symbol2CountCoverageSet12.bq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(symbol));
        }
        fmt.cAltHD[strand] = symbol2CountCoverageSet12.fq_hiqual_dep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.aDB[strand]  = symbol2CountCoverageSet12.du_bias_dedup.at(strand).getByPos(refpos).getSymbolCount(symbol);

        fmt.bPTL[strand] = symbol2CountCoverageSet12.bq_amax_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.bPTR[strand] = symbol2CountCoverageSet12.bq_amax_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.bPBL[strand] = symbol2CountCoverageSet12.bq_bias_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bPBR[strand] = symbol2CountCoverageSet12.bq_bias_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.bMMT[strand] = symbol2CountCoverageSet12.bq_amax_nvars.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.bMMB[strand] = symbol2CountCoverageSet12.bq_bias_nvars.at(strand).getByPos(refpos).getSymbolCount(symbol);
        
        fmt.bSDL[strand] = symbol2CountCoverageSet12.bq_bsum_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bSDR[strand] = symbol2CountCoverageSet12.bq_bsum_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bSB1[strand] = symbol2CountCoverageSet12.bq_bias_1stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bSBR[strand] = symbol2CountCoverageSet12.bq_bias_2stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bNSB[strand] = symbol2CountCoverageSet12.bq_n_seq_bases.at(strand).getByPos(refpos).getSymbolCount(symbol);
        
        fmt.bAD1[strand] = symbol2CountCoverageSet12.bq_tsum_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // total allele depth
        fmt.bAD2[strand] = symbol2CountCoverageSet12.bq_pass_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth
        fmt.bQT2[strand] = symbol2CountCoverageSet12.bq_pass_thres.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  threshold
        fmt.bAD3[strand] = symbol2CountCoverageSet12.bq_vars_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth
        fmt.bADB[strand] = symbol2CountCoverageSet12.bq_vars_badep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        if (MAX_IMBA_DEP == fmt.bADB[strand]) { fmt.bADB[strand] = (uint32_t)(UINT32_MAX); }
        fmt.bQT3[strand] = symbol2CountCoverageSet12.bq_vars_thres.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  threshold
        fmt.bVQ3[strand] = symbol2CountCoverageSet12.bq_vars_vqual.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth 
        
        for (unsigned int seqdir = 0; seqdir < 2; seqdir++) {
            fmt.aAD[strand*2+seqdir] = symbol2CountCoverageSet12.bq_dirs_count.at(strand*2+seqdir).getByPos(refpos).getSymbolCount(symbol);
            fmt.aNMAD[strand*2+seqdir] = symbol2CountCoverageSet12.bq_nms_count.at(strand*2+seqdir).getByPos(refpos).getSymbolCount(symbol);
            fmt.aBQAD[strand*2+seqdir] = symbol2CountCoverageSet12.bq_dirs_bqsum.at(strand*2+seqdir).getByPos(refpos).getSymbolCount(symbol);
        }
        fmt.aBAQADR[1] = symbol2CountCoverageSet12.bq_baq_sum.at(0).getByPos(refpos).getSymbolCount(symbol);
        
        assert(fmt.aB.size() == symbol2CountCoverageSet12.bq_bias_sedir.size());
        for (unsigned int idx = 0; idx < symbol2CountCoverageSet12.bq_bias_sedir.size(); idx++) { 
            fmt.aB[idx] = symbol2CountCoverageSet12.bq_bias_sedir.at(idx).getByPos(refpos).getSymbolCount(symbol);
        }
        // double bq_qsum_rawMQ = (double)symbol2CountCoverageSet12.bq_qsum_rawMQ.at(strand).getByPos(refpos).getSymbolCount(symbol);
        double bq_qsum_sqrMQ = (double)symbol2CountCoverageSet12.bq_qsum_sqrMQ.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bMQ1[strand] = sqrt(bq_qsum_sqrMQ / (DBL_MIN + (double)fmt.bAD1[strand])); // total allele depth
        // fmt.bMQ2[strand] = bq_qsum_sqrMQ / (DBL_MIN + bq_qsum_rawMQ);
        bq_qsum_sqrMQ_tot += bq_qsum_sqrMQ;
        
        // double bq_qual_p1sum = (double)symbol2CountCoverageSet12.bq_qual_p1sum.at(strand).getByPos(refpos).getSymbolCount(symbol);
        double bq_qual_p2sum = (double)symbol2CountCoverageSet12.bq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bBQ1[strand] = sqrt(bq_qual_p2sum / (DBL_MIN + (double)fmt.bAD1[strand]));
        //fmt.bBQ2[strand] = bq_qual_p2sum / (DBL_MIN + bq_qual_p1sum);
        
        // fmt.bEDAD[strand] = symbol2CountCoverageSet12.bq_n_edit_ops.at(strand).getByPos(refpos).getSymbolCount(symbol);
        
        //double fq_qual_p1sum = (double)symbol2CountCoverageSet12.fq_qual_p1sum.at(strand).getByPos(refpos).getSymbolCount(symbol);
        //double fq_qual_p2sum = (double)symbol2CountCoverageSet12.fq_qual_p2sum.at(strand).getByPos(refpos).getSymbolCount(symbol);
        //fmt.cCQ1[strand] = sqrt(fq_qual_p2sum / (DBL_MIN + (double)fmt.cAD1[strand]));
        //fmt.cCQ2[strand] = fq_qual_p2sum / (DBL_MIN + fq_qual_p1sum); 
        
        fmt.bADLQ[strand] = symbol2CountCoverageSet12.bq_tsum_LQdep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        
        // family quality
        fmt.cPTL[strand] = symbol2CountCoverageSet12.fq_amax_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cPTR[strand] = symbol2CountCoverageSet12.fq_amax_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cPBL[strand] = symbol2CountCoverageSet12.fq_bias_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cPBR[strand] = symbol2CountCoverageSet12.fq_bias_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cMMT[strand] = symbol2CountCoverageSet12.fq_amax_nvars.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cMMB[strand] = symbol2CountCoverageSet12.fq_bias_nvars.at(strand).getByPos(refpos).getSymbolCount(symbol);

        fmt.cSDL[strand] = symbol2CountCoverageSet12.fq_bsum_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cSDR[strand] = symbol2CountCoverageSet12.fq_bsum_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cSB1[strand] = symbol2CountCoverageSet12.fq_bias_1stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cSBR[strand] = symbol2CountCoverageSet12.fq_bias_2stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cNSB[strand] = symbol2CountCoverageSet12.fq_n_seq_bases.at(strand).getByPos(refpos).getSymbolCount(symbol);
        
        fmt.cAD1[strand] = symbol2CountCoverageSet12.fq_tsum_depth.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cAD2[strand] = symbol2CountCoverageSet12.fq_pass_depth.at(strand).getByPos(refpos).getSymbolCount(symbol);  
        fmt.cQT2[strand] = symbol2CountCoverageSet12.fq_pass_thres.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cAD3[strand] = symbol2CountCoverageSet12.fq_vars_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth
        fmt.cADB[strand] = symbol2CountCoverageSet12.fq_vars_badep.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  threshold
        if (MAX_IMBA_DEP == fmt.cADB[strand]) { fmt.cADB[strand] = (uint32_t)(UINT32_MAX); }
        fmt.cQT3[strand] = symbol2CountCoverageSet12.fq_vars_thres.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cVQ3[strand] = symbol2CountCoverageSet12.fq_vars_vqual.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth 
        
        fmt.cMajor[strand] = symbol2CountCoverageSet12.major_amplicon.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cMinor[strand] = symbol2CountCoverageSet12.minor_amplicon.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cADTT[strand] = symbol2CountCoverageSet12.fam_total_dep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cADT1[strand] = symbol2CountCoverageSet12.fam_size1_dep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        auto fmt_cADTN = symbol2CountCoverageSet12.fam_nocon_dep.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cADTC[strand] = fmt.cADTT[strand] - fmt.cADT1[strand] - fmt_cADTN;
        // fmt.gapNum[strand] = 0;
        if ((0 < fmt.bAD1[strand]) && (isSymbolIns(symbol) || isSymbolDel(symbol))) {
            // auto cADdiff_cADtotal = fill_by_indel_info(fmt, symbol2CountCoverageSet12, strand, refpos, symbol, refstring, repeatunit, repeatnum);
            //fmt.gapcADD[strand] = cADdiff_cADtotal[0]; // diff
            //fmt.gapcADT[strand] = cADdiff_cADtotal[1];
        }
    }
    
    for (unsigned int region = 0; region < symbol2CountCoverageSet12.bq_regs_count.size(); region++) {
        fmt.aPBAD[region] = symbol2CountCoverageSet12.bq_regs_count.at(region).getByPos(refpos).getSymbolCount(symbol);
    }
    
    fmt.bHap = mutform2count4map_to_phase(mutform2count4vec_bq, indices_bq);
    fmt.cHap = mutform2count4map_to_phase(mutform2count4vec_fq, indices_fq);
    
    fmt.dAD1 = symbol2CountCoverageSet12.duplex_tsum_depth.getByPos(refpos).getSymbolCount(symbol);
    fmt.dAD3 = symbol2CountCoverageSet12.duplex_pass_depth.getByPos(refpos).getSymbolCount(symbol);
    
    unsigned int regionpos = refpos - refstring_offset;
    const std::string vcfref = refstring.substr(regionpos, 1);
    const std::string vcfalt = std::string(SYMBOL_TO_DESC_ARR[symbol]);
    
    bool is_novar = (symbol == LINK_M || (isSymbolSubstitution(symbol) && vcfref == vcfalt));
    
    fmt.bDP = fmt.bDP1[0] + fmt.bDP1[1];
    fmt.bADR = {{ fmt.bRD1[0] + fmt.bRD1[1], fmt.bAD1[0] + fmt.bAD1[1] }}; 
    
    fmt.cDP = fmt.cDPTT[0] + fmt.cDPTT[1];
    fmt.cADR = {{ fmt.cRDTT[0] + fmt.cRDTT[1], fmt.cADTT[0] + fmt.cADTT[1] }}; 
    
    auto fmtAD = SIGN2UNSIGN(0);
    if (use_deduplicated_reads) {
        fmt.DP = fmt.cDP;
        fmtAD = fmt.cADR[1];
        fmt.FA = ((double)fmt.cADR[1]) / ((double)fmt.cDP + DBL_EPSILON);
        fmt.FR = ((double)fmt.cADR[0]) / ((double)fmt.cDP + DBL_EPSILON);
    } else {
        fmt.DP = fmt.bDP;
        fmtAD = fmt.bADR[1];
        fmt.FA = ((double)fmt.bADR[1]) / ((double)fmt.bDP + DBL_EPSILON);
        fmt.FR = ((double)fmt.bADR[0]) / ((double)fmt.bDP + DBL_EPSILON);
    }
    
    fmt.DPHQ = fmt.bDP - (fmt.bDPLQ[0] + fmt.bDPLQ[1]);
    fmt.ADHQ = fmt.bADR[1] - (fmt.bADLQ[0] + fmt.bADLQ[1]);
    assert(fmt.DPHQ >= 0);
    assert(fmt.ADHQ >= 0);
    double dsBQ1 = sqrt((mathsquare(fmt.bBQ1[0]) * (double)fmt.bAD1[0] + mathsquare(fmt.bBQ1[1]) * (double)fmt.bAD1[1]) / (double)(DBL_MIN + fmt.bAD1[0] + fmt.bAD1[1])); 
    fmt.BQ = (unsigned int)dsBQ1;
    fmt.MQ = (unsigned int)sqrt((double)bq_qsum_sqrMQ_tot / (DBL_MIN + (double)(fmt.bAD1[0] + fmt.bAD1[1])));
    
    /* Schematics of the alignment hierarchy of reads for del
     * (012345RR89abcd) is the reference and R is the repeated sequence.
     *  012345
     *   12345R  reject, not realignable
     *    2345R8  reject, partially re-alignable, clipped
     *     345-R89 accept
     *      45-R89a accept
     *       5-R89ab accept
     *         R89abc reject, partially re-alignable, clipped (not counted for depth)
     *          89abcd reject, not realignable (not counted for depth)
     * (012345RRR89abcd) is the reference and R is the repeated sequence.
     *  012345
     *   12345R reject, not realignable
     *    2345RR reject, not realignable
     *     345RR8 reject, partially re-alignable, clipped
     *      45-RR89 accept
     *       5-RR89a accept
     *         RR89ab reject, partially re-alignable, clipped (not counted for depth)
     *          R89abc reject, not realignable (not counted for depth)
     *           89abcd reject, not realignable (not counted for depth)
     */
    /* Schematics of the alignment hierarchy of reads for ins
     * (012345-R89abcd) is the reference and R is the repeated sequence.
     *  012345
     *   12345R reject, not realignable 
     *    2345RR reject, partially realignable, clipped
     *     345RR8 reject, fully realignable, clipped
     *      45RR89 accept
     *       5RR89a accept
     *        RR89ab reject, fully realignable, clipped
     *         R89abc reject, partially realignable, clipped (not counted for depth)
     *          89abcd reject, not realignable (not counted for depth)
     */
    double refmul = 1.0;
    double altmul = 1.0;
    bool isInDel = (isSymbolIns(symbol) || isSymbolDel(symbol));
    if (fmtAD > 0 || is_rescued) {
        assert(fmt.FA >= 0);
        unsigned int ref_bias = 0; // 2; // 4 * 2; // this is rather empirical
        if (isInDel) {
            uint64_t totsize_cnt = 0;
            uint64_t totsize_sum = 0;
            for (unsigned int k = 0; k < fmt.gapSeq.size(); k++) {
                totsize_cnt += fmt.gapcAD1[k]; 
                totsize_sum += fmt.gapSeq[k].size() * fmt.gapcAD1[k];
            }
            auto indel_bias_len = ((totsize_sum * 100UL) / (totsize_cnt * 100UL + 1UL));
            auto context_bias = (repeatunit.size() * MAX(1UL, repeatnum));
            // 1 = match score, 4 = mismatch-penalty, 5 = softclip-penalty, 6 = gap-open-penalty
            if (isSymbolIns(symbol)) {
                ref_bias = (indel_bias_len * 3) + (6 * 2) + (4 * 2) + context_bias;
            } else {
                ref_bias = (indel_bias_len * 2) + (6 * 2) + (4 * 2) + context_bias;
            }
        }
        fmt.RefBias = ref_bias;
        if (!is_ref_bias_aware) { ref_bias = 0; }
        unsigned int readlen = MAX(central_readlen, 30U); // 30 is the min alignment score
        altmul = (double)(readlen - MIN(readlen - 30, ref_bias)) / (double)readlen; // 50.0 / (double)(ref_bias + 50);
        refmul = 2.0 - altmul;
        auto ref_mul_contam_frac = ((double)MIN(ref_bias, readlen/2) / (double)(readlen));
        fill_TN_germline(
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
    fmt.HQ[0] = 0; 
    fmt.HQ[1] = 0;
    
    std::array<RepNumCluster<int>, RCC_NUM> rep_num_clusters;
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
    }
    
    indel_fill_rep_num_clusters(rep_num_clusters, 
             fmt.cADR[0], repeatunit, iseq2cnt_vec, dlen2cnt_vec);
    
    // std::string indelstring = indel_get_majority(fmt, prev_is_tumor, tki); 
    for (unsigned int i = 0; i < RCC_NUM; i++) {
        fmt.RCC[i*RCC_NFS  ] = (int)rep_num_clusters[i].mode;
        fmt.RCC[i*RCC_NFS+1] = (int)rep_num_clusters[i].cnt2m; // indel of two units
        fmt.RCC[i*RCC_NFS+2] = (int)rep_num_clusters[i].cnt1m;
        fmt.RCC[i*RCC_NFS+3] = (int)rep_num_clusters[i].cnt0;
        fmt.RCC[i*RCC_NFS+4] = (int)rep_num_clusters[i].cnt1p;
        fmt.RCC[i*RCC_NFS+5] = (int)rep_num_clusters[i].cnt2p;
    } // fmt.cFR;
    
    unsigned int refpo1 = (isInDel ? (MAX(refpos, 1 + bq_ins_2bdepths.at(0).getIncluBegPosition()) - 1) : refpos);
    unsigned int refpo2 = MIN(refpos + 1, bq_ins_2bdepths.at(0).getExcluEndPosition() - 1);
    // Implicitly assume that there is a symmetry between refpos and refpo2
    fmt.gapbNRD = {
        MAX(bq_ins_2bdepths.at(0).getByPos(refpo1), bq_ins_2bdepths.at(0).getByPos(refpo2)), 
        MAX(bq_ins_2bdepths.at(1).getByPos(refpo1), bq_ins_2bdepths.at(1).getByPos(refpo2)), 
        MAX(bq_del_2bdepths.at(0).getByPos(refpo1), bq_del_2bdepths.at(0).getByPos(refpo2)), 
        MAX(bq_del_2bdepths.at(1).getByPos(refpo1), bq_del_2bdepths.at(1).getByPos(refpo2)) 
    };
    // fmt.gapbNNRD = { bq_indel_adjmax_depths.at(0).getByPos(refpos), bq_indel_adjmax_depths.at(1).getByPos(refpos) };
    
    fmt.VType = SYMBOL_TO_DESC_ARR[symbol];
    double lowestVAQ = prob2phred(1 / (double)(fmt.bAD1[0] + fmt.bAD1[1] + 1)) * ((fmt.bAD1[0] + fmt.bAD1[1]) / (fmt.bDP1[0] + fmt.bDP1[1] + DBL_MIN)) / (double)2;
    double stdVAQs[2] = {0, 0};
    std::array<unsigned int, 2> weightedQT3s = {{0, 0}};
    for (size_t i = 0; i < 2; i++) {
        auto minAD1 = 1; 
        auto gapAD1 = 0;
        if (use_deduplicated_reads) {
            if (use_only_deduplicated_reads) {
                minAD1 = 0;
                gapAD1 = 1;
            } else {
                minAD1 = MIN(fmt.bAD1[i], fmt.cAD1[i]); // prevent symbol with suffix = _N
                gapAD1 = MAX(fmt.bAD1[i], fmt.cAD1[i]) - minAD1;
                if (fmt.bVQ3[i] < fmt.cVQ3[i]) {
                    gapAD1 *= 3;
                }
            }
        }
        double currVAQ = (fmt.bVQ3[i] * minAD1 + fmt.cVQ3[i] * gapAD1) / (double)(minAD1 + gapAD1 + DBL_MIN); // prevent div by zero
        weightedQT3s[i] = (fmt.bQT3[i] * minAD1 + fmt.cQT3[i] * gapAD1) / (double)(minAD1 + gapAD1 + DBL_MIN);
        stdVAQs[i] = currVAQ;
        
        /*
        if (fmt.bSDL[i] < fmt.bAD1[i] * min_edge_dist && fmt.cSDL[i] < fmt.cAD1[i] * min_edge_dist) {
            double edgeVAQ = ((double)baq_per_aligned_base) * MAX((double)fmt.bSDL[i] / (double)(fmt.bAD1[i] + DBL_MIN), (double)fmt.cSDL[i] / (double)(fmt.cAD1[i] + DBL_MIN));
            stdVAQs[i] = MIN(stdVAQs[i], edgeVAQ);
        }
        if (fmt.bSDR[i] < fmt.bAD1[i] * min_edge_dist && fmt.cSDR[i] < fmt.cAD1[i] * min_edge_dist) {
            double edgeVAQ = ((double)baq_per_aligned_base) * MAX((double)fmt.bSDR[i] / (double)(fmt.bAD1[i] + DBL_MIN), (double)fmt.cSDR[i] / (double)(fmt.cAD1[i] + DBL_MIN));
            stdVAQs[i] = MIN(stdVAQs[i], edgeVAQ);
        }
        */
        /*
        if ((fmt.bBQ1[i] < minABQ)) {
            // the following line of code is not theoretically sound
            // stdVAQs[i] = MIN(stdVAQs[i], fmt.bBQ1[i]);
        }
        if ((unsigned int)fmt.bMQ1[i] < minMQ1) {
            // the following line of code is not theoretically sound
            // stdVAQs[i] = MIN(stdVAQs[i], MIN((unsigned int)(fmt.bMQ1[i] * 2), maxMQ)); // 60 is max MAPQ of bwa
        }
        */
        fmt.cVAQ1[i] = currVAQ;
    }
    // Ideally (no read supports a variant at its border, there is no mismatch in any read other than at the variant site) the variable below has a value of 4.
    // Even more ideally (there is no mismatch in the entire contig other than at the variant site) the variable below increases logarithmically as a function of variant depth.
    const double contig_to_frag_len_ratio = (double)2;
    // double vaqBQcap = (double)FLT_MAX; // ((fmt.BQ < minABQ) ? ((double)fmt.BQ) : ((double)FLT_MAX));
    // sqrt(bq_qual_p2sum / (DBL_MIN + (double)fmt.bAD1[strand]));
    
    // double vaqMQcap = ((fmt.MQ < minMQ1) ? (fmt.MQ * contig_to_frag_len_ratio) : ((double)FLT_MAX));
    //double minVAQ = MIN(stdVAQs[0], stdVAQs[1]);
    //double stdVAQ = MAX(stdVAQs[0], stdVAQs[1]);
    
    double weightsum = MIN((double)(weightedQT3s[0] + weightedQT3s[1]), phred_max_dscs);
    double doubleVAQfw = stdVAQs[0] + stdVAQs[1] * MIN(1.0, (weightsum - weightedQT3s[0]) / (weightedQT3s[0] + DBL_MIN));
    double doubleVAQrv = stdVAQs[1] + stdVAQs[0] * MIN(1.0, (weightsum - weightedQT3s[1]) / (weightedQT3s[1] + DBL_MIN));
    fmt.cVAQ2 = {(float)doubleVAQfw, (float)doubleVAQrv};
    
    //double doubleVAQ_multnorm =(double)(1 + fmt.gapcADD[0] + fmt.gapcADD[1]) / (double)(1 + fmt.gapcADT[0] + fmt.gapcADT[1]);
    double doubleVAQ = MAX(doubleVAQfw, doubleVAQrv);
    //double doubleVAQ_norm = doubleVAQ * doubleVAQ_multnorm;
    // double doubleVAQ = stdVAQ + (minVAQ * (phred_max_dscs - phred_max_sscs) / (double)phred_max_sscs);
    double duplexVAQ = (double)fmt.dAD3 * (double)(phred_max_dscs - phred_max_sscs) - (double)(fmt.dAD1 - fmt.dAD3); // h01_to
    duplexVAQ = MIN(duplexVAQ, 200); // Similar to many other upper bounds, the 200 here has no theoretical foundation.
    double duprate = 1.0 - MIN(1.0, (double)(fmt.cDP + 1) / (double)(fmt.bDP + 1));
    double dimret_coef = phred_umi_dimret_mult * duprate + (1.0 - duprate);
    
    
    // const double pl_cap = 90.0;
    // const bool tUseHD1 = false;
    // const double phred_triallelic_indel = 30.0;
    
    //double rawVQ1 = dimret_coef * MIN(vaqBQcap, MAX(lowestVAQ, doubleVAQ + duplexVAQ));       // / 1.5;
    //double rawVQ2 = dimret_coef * MIN(vaqBQcap, MAX(lowestVAQ, doubleVAQ_norm + duplexVAQ));
    // std::string indelstring;
    double somatic_var_pq = 5.0;
    if (LINK_M == symbol) {
        somatic_var_pq = 0;
    } else if (isInDel) {
        /*
        if (is_rescued) {
            bad0a_indelstring_tkiidx_vec.
        } else {
            //const auto indelstrings = indel_get_majority(fmt,// prev_is_tumor, tki, false,
            //    NULL, refpos, symbol, false);
            indelstring = indelstrings.at(sivi).second;
        }
        */
        somatic_var_pq = (int)MIN(indel_phred(8.0, indelstring.size(), repeatunit.size(), repeatnum), 24) + (5-3) - (int)phred_snv_to_indel_ratio;
    } 
    double rawVQ1 = dimret_coef * MAX(lowestVAQ, doubleVAQ + duplexVAQ) / (tUseHD1 ? 1.0 : sqrt(altmul)) + somatic_var_pq;       // / 1.5;
    // double rawVQ2 = dimret_coef * MAX(lowestVAQ, doubleVAQ_norm + duplexVAQ) / sqrt(altmul) + somatic_var_pq;
    
    // intuition for the somatic_var_pq formula: InDel candidate with higher baseline-noise frequency is simply more likely to be a true mutation too
    // assuming that in-vivo biological error and in-vitro technical error are positively correlated with each other. 
    // https://doi.org/10.1007%2Fs00239-010-9377-4
    
    unsigned int pcap_baq = sqrt(fmt.aBAQADR[1] * (double)THE_BAQ_MAX / MAX(SUMVEC(fmt.aAD), 1));
    // const double pcap_bcq = (isInDel ? 200.0 : (mathsquare(fmt.BQ) * illumina_BQ_sqr_coef - 10.0/log(10.0) * log((double)MAX(100, fmt.aB[0]) / 100.0))); // based on heuristics
    //const unsigned int fwTotBQ = (fmt.aBQAD[0] + fmt.aBQAD[2]);
    //const unsigned int rvTotBQ = (fmt.aBQAD[1] + fmt.aBQAD[3]);
    //const unsigned int d2TotBQ = MAX(fwTotBQ, rvTotBQ) + MIN(fwTotBQ, rvTotBQ) * 5/2;
    //const auto pcap_bcq = (unsigned int)(isInDel ? 200.0 : (d2TotBQ / ((double)SUMVEC(fmt.aAD) + DBL_EPSILON) * illumina_BQ_sqr_coef)); // based on heuristics
    
    const unsigned int fwTotAD = (fmt.aAD[0] + fmt.aAD[2]);
    const unsigned int rvTotAD = (fmt.aAD[1] + fmt.aAD[3]);
    const unsigned int fwTotOD = (fmt.aDP[0] + fmt.aDP[2]) - fwTotAD;
    const unsigned int rvTotOD = (fmt.aDP[1] + fmt.aDP[3]) - rvTotAD;
    double coefBQ = 0;
    if (fwTotAD < rvTotAD) {
        coefBQ = 2.0 + 2.0 * (fwTotAD + 0.5) / (double)(fwTotAD + rvTotAD + 2) - (fwTotOD + 0.5) / (double)(fwTotOD + rvTotOD + fwTotAD + 2); 
    } else {
        coefBQ = 2.0 + 2.0 * (rvTotAD + 0.5) / (double)(fwTotAD + rvTotAD + 2) - (rvTotOD + 0.5) / (double)(fwTotOD + rvTotOD + rvTotAD + 2);
    }
    // const auto pcap_bcq = (unsigned int)((isInDel || is_proton) ? 200.0 : (d2TotBQ / ((double)SUMVEC(fmt.aAD) + DBL_EPSILON) * illumina_BQ_sqr_coef)); 
    const auto pcap_bcq = (unsigned int)((isInDel || is_proton) ? 200.0 : (unsigned int)(coefBQ * mathsquare(fmt.BQ) / 20.0));
    
    const auto pcap_tmq1 = (unsigned int)MIN(maxMQ, fmt.MQ) + phred_varcall_err_per_map_err_per_base; // bases on heuristics
    const auto pcap_tmq2 = (unsigned int)MAX(3.0, fmt.MQ - 10.0/log(10.0) * log((double)(fmt.DP + 1) / (double)(fmt.DP * fmt.FA + 0.5))) * ((double)fmt.DP * fmt.FA); // readjustment by MQ
    
    // aln quality, base quality, read count,
    double normBAQAD = fmt.aBAQADR[1] / altmul;
    double normBAQOD = (fmt.aBAQDP - fmt.aBAQADR[1]) / refmul;
    const double mul_per_baq = (isInDel ? INDEL_MUL_PER_BAQ : SNV_MUL_PER_BAQ);
    double alsFA2 = (normBAQAD * mul_per_baq + 0.5) / (normBAQAD * mul_per_baq + normBAQOD + 1.0); // laplace smoothing
    double blsFA1 = (fmt.bADR[1] / altmul + 0.5) / ((fmt.bDP - fmt.bADR[1]) / refmul + fmt.bADR[1] / altmul + 1.0);
    double clsFA1 = (fmt.cADR[1] / altmul + 0.5) / ((fmt.cDP - fmt.cADR[1]) / refmul + fmt.cADR[1] / altmul + 1.0);
    const double powlaw_qual = powlaw_exponent * 10.0 / log(10.0) * log(MIN3(alsFA2, blsFA1, clsFA1)) + powlaw_anyvar_base;
    
    const double sampling_qual = 40.0 * pow(0.5, MAX((double)fmt.bADR[1], ((double)(fmt.bDP + 1)) / ((double)(fmt.cDP + 1)) - 1.0));
    double indel_ic = 0.0; // applies to SNV too
    double indel_p2 = 0;
    if (isInDel && !tUseHD1) {
        int n_str_units = (isSymbolIns(symbol) ? 1 : (-1)) * (int)(indelstring.size() / repeatunit.size());
        const double symbol_to_allele_frac = 1.0 - pow((isSymbolIns(symbol) ? 0.9 : (isSymbolDel(symbol) ? 0.95 : 1.0)), indelstring.size());
        double tAD0a = fmt.cADR[1] * (double)(indelbdepth + 0.5) / (double)(1 + fmt.bADR[1]); // fmt.cADR[1] * ((double)(fmt.gapDP4[2] + 1) / (double)(fmt.gapDP4[0] + 1));
        indel_p2 =  penal_indel_2(tAD0a, n_str_units, fmt.RCC, phred_triallelic_indel);
        indel_ic = 10.0/log(10.0) * log((double)MAX(indelstring.size(), 1U) / (double)(repeatunit.size() * (MAX(1, repeatnum) - 1) + 1));
        // rawVQ1 = MIN(rawVQ1 * (tAD0a + (double)(fmt.cADR[1] - tAD0a) * symbol_to_allele_frac) / ((double)fmt.cADR[1]), rawVQ1 - indel_p2); // ad-hoc adjustment
        rawVQ1 *= (tAD0a + (double)(fmt.cADR[1] - tAD0a) * symbol_to_allele_frac) / ((double)fmt.cADR[1]);
    }
    
    fmt.VQ1.clear();
    fmt.VQ2.clear();
    fmt.VQ3.clear();
    fmt.VQ4.clear();
    fmt.VQ5.clear();
    fmt.VAQ.clear();
    
    clear_push(fmt.VQ1, (int)(pcap_baq + (isInDel ? 10 : 0))); // base alignment quality
    clear_push(fmt.VQ2, (int)pcap_bcq); // base quality
    clear_push(fmt.VQ3, (int)MIN(pcap_tmq1, pcap_tmq2)); // mapping quality
    clear_push(fmt.VQ4, (int)(powlaw_qual + (tUseHD1 ? (fmt.dAD3 > 0 ? powlaw_dscs_inc : powlaw_sscs_inc) : 0) + indel_ic - indel_p2));
    clear_push(fmt.VQ5, (int)(rawVQ1 - sampling_qual));
    auto theVAQ = MIN5(fmt.VQ1[0], fmt.VQ2[0], fmt.VQ3[0], fmt.VQ4[0], fmt.VQ5[0]);
    clear_push(fmt.VAQ, MAX(0, theVAQ));
    fmt.ALODQ = MAX(0, fmt.VAQ[0]);
    //return (int)(fmt.bAD1[0] + fmt.bAD1[1]);
    return 0;
};

#include "version.h"
std::string 
generate_vcf_header(const char *ref_fasta_fname, 
        const char *platform, 
        unsigned int central_readlen,
        const unsigned int minABQ_pcr_snv, 
        const unsigned int minABQ_pcr_indel, 
        const unsigned int minABQ_cap_snv, 
        const unsigned int minABQ_cap_indel, 
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
    ret += std::string("") + "##variantCallerInferredParameters=<" + "inferred_sequencing_platform=" + platform + ",central_readlen=" + std::to_string(central_readlen) + ",minABQs=("
            + std::to_string(minABQ_pcr_snv) + "x" + std::to_string(minABQ_pcr_indel) + "x" +  std::to_string(minABQ_cap_snv) + "x" + std::to_string(minABQ_cap_indel) + ")>\n";
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
#if (1 < N_MODELS)
    for (int i = 0; i < N_MODELS; i++) {
        // ret += std::string(";TQ") + std::to_string(i) + "=" + std::to_string(testquals[i]); 
        ret += "##INFO=<ID=TQ" + std::to_string(i) +",Number=1,Type=Float,Description=\"Variant quality computed by the model " + std::to_string(i) +"\">\n";
    }
#endif
    ret += "##INFO=<ID=SomaticQ,Number=A,Type=Float,Description=\"Somatic quality of the variant, the PHRED-scale probability that this variant is not somatic.\">\n";
    ret += "##INFO=<ID=TLODQ,Number=A,Type=Float,Description=\"Tumor log-of-data-likelihood quality, the PHRED-scale probability that this variant is not of biological origin (i.e., artifactual).\">\n";
    ret += "##INFO=<ID=NLODQ,Number=A,Type=Float,Description=\"Normal log-of-data-likelihood quality, the PHRED-scale probability that this variant is of germline origin.\">\n";
    ret += "##INFO=<ID=QUAL1,Number=A,Type=Float,Description=\"Tumor-only quality without applying universality (deprecated but kept for compatibility with older versions).\">\n";
    // TLODQ1 may be useful in theory, but more data is needed to assess TLODQ1.
    // ret += "##INFO=<ID=TLODQ1,Number=1,TypeFloat,Description=\"Tumor log-of-data-likelihood quality, the PHRED-scale probability that this variant is not of biological origin (i.e., artifactual).\">\n";
    ret += "##INFO=<ID=REFQs,Number=.,Type=Float,Description=\"Non-germline qualities: normal-ALT, tumor-ALT, normal-OTHER, tumor-OTHER\">\n";
    
    ret += "##INFO=<ID=TNQs,Number=.,Type=Float,Description=\"TvsN qualities: baseline, read-transfer contamination by normal (RTCN) baseline, RTCN error, and systematic error.\">\n";
    ret += "##INFO=<ID=TNORQs,Number=.,Type=Float,Description=\"TvsN qualities: ATL-to-NONALT (ANA) odds ratio, defined as tumor ANA divided by normal ANA, using AD and BQ-sum\">\n";
    ret += "##INFO=<ID=TN1Qs,Number=.,Type=Float,Description=\"TvsN qualities: AD power-law, BQ-sum power-law, AD binomial, and BQ-sum binomial\">\n";
    ret += "##INFO=<ID=TN2Qs,Number=.,Type=Float,Description=\"TvsN qualities: TN1Qs for the OTHER covering alleles that are neither ALT nor REF\">\n";
    
    ret += "##INFO=<ID=TAQs,Number=.,Type=Float,Description=\"Tumor-only  qualities: baseline, penalties by adjacent InDels (to only InDel) and nearby InDels (to both SNV and InDel).\">\n";
    ret += "##INFO=<ID=TSQs,Number=.,Type=Float,Description=\"Tumor-only  qualities: normalized power-law, binomial, FA power-law, and coverage-depth adjustment.\">\n";
    ret += "##INFO=<ID=NSQs,Number=.,Type=Float,Description=\"Normal-only qualities: normalized power-law, binomial, FA power-law, and coverage-depth adjustment.\">\n";
    
    // ret += "##INFO=<ID=tDP,Number=1,Type=Integer,Description=\"Tumor-sample DP\">\n";
    ret += "##INFO=<ID=tFA,Number=1,Type=Float,Description=\"Tumor-sample FA (deprecated, equivalent to AD[1] divided by DP in tAD)\">\n";
    // ret += "##INFO=<ID=tFR,Number=1,Type=Float,Description=\"Tumor-sample FR (deprecated, equivalent to AD[0] divided by DP in tAD)\">\n";
    ret += "##INFO=<ID=tADR,Number=R,Type=Integer,Description=\" Tumor-sample RD and AD for tissue (RD: depth of the REF, AD: depth of each ALT)\">\n";
    ret += "##INFO=<ID=tDP,Number=1,Type=Integer,Description=\" Tumor-sample DP for tissue\">\n";
    ret += "##INFO=<ID=nADR,Number=R,Type=Integer,Description=\"Normal-sample RD and AD for tissue\">\n";
    ret += "##INFO=<ID=nDP,Number=1,Type=Integer,Description=\"Normal-sample DP for tissue\">\n";
    ret += "##INFO=<ID=tADCR,Number=R,Type=Integer,Description=\" Tumor-sample RD and AD for ctDNA (RD: depth of the REF, AD: depth of each ALT)\">\n";
    ret += "##INFO=<ID=tDPC,Number=1,Type=Integer,Description=\" Tumor-sample DP for ctDNA\">\n";
    ret += "##INFO=<ID=nADCR,Number=R,Type=Integer,Description=\"Normal-sample RD and AD for ctDNA\">\n";
    ret += "##INFO=<ID=nDPC,Number=1,Type=Integer,Description=\"Normal-sample DP for ctDNA\">\n";
    ret += "##INFO=<ID=tADBQR,Number=R,Type=Integer,Description=\"Tumor-sample (cRefBQ,cAltBQ) or (bRefBQ,bAltBQ), depending on the command-line option\">\n";
    ret += "##INFO=<ID=tDPBQ,Number=1,Type=Integer,Description=\"Tumor-sample cAllBQ or bAllBQ, depending on the command-line option\">\n";
    //ret += "##INFO=<ID=tAltBQ2,Number=1,Type=Integer,Description=\"Tumor-sample cAltBQ2\">\n";
    //ret += "##INFO=<ID=tAllBQ2,Number=1,Type=Integer,Description=\"Tumor-sample cAllBQ2\">\n";
    //ret += "##INFO=<ID=tRefBQ2,Number=1,Type=Integer,Description=\"Tumor-sample cRefBQ2\">\n";
    ret += "##INFO=<ID=tADHDR,Number=R,Type=Integer,Description=\"Tumor-sample (cRefHD,cAltHD) or (bRefHD,bAltHD), depending on the command-line option\">\n";
    ret += "##INFO=<ID=tDPHD,Number=1,Type=Integer,Description=\"Tumor-sample cAllHD or bAllHD, depending on the command-line option\">\n";
    ret += "##INFO=<ID=tFTS,Number=1,Type=String,Description=\"Tumor-sample FTS where the filter strings are separated by amperstand (&)\">\n";
    ret += "##INFO=<ID=tbDP,Number=1,Type=Integer,Description=\"Tumor-sample bDP\">\n";
    ret += "##INFO=<ID=tcHap,Number=1,Type=String,Description=\"Tumor-sample cHap\">\n";
    ret += "##INFO=<ID=tMQ,Number=.,Type=Float,Description=\"Tumor-sample MQ\">\n"; 
    ret += "##INFO=<ID=tEROR,Number=5,Type=Integer,Description=\"Tumor-sample EROR\">\n";
    ret += "##INFO=<ID=tgapDP4,Number=4,Type=Integer,Description=\"Tumor-sample gapDP4\">\n"; 
    ret += "##INFO=<ID=tRCC,Number=" + std::to_string(RCC_NFS*RCC_NUM) + ",Type=Integer,Description=\"Tumor-sample RCC\">\n";
    // ret += "##INFO=<ID=tGLa,Number=3,Type=Integer,Description=\"Tumor-sample GLa\">\n";
    // ret += "##INFO=<ID=tGLb,Number=3,Type=Integer,Description=\"Tumor-sample GLb\">\n";
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

auto
fmtFTSupdate(auto & maxval, std::string & ft, std::vector<unsigned int> & ftv, const char *fkey , const auto fthres, const auto fval) {
    maxval = MAX(maxval, fval);
    if ((unsigned int)fthres <= (unsigned int)fval) {
        ft  += (std::string(fkey) + "&"); // amperstand-separated list of filter strings
        ftv.push_back((unsigned int)fval);
    }
    return fval;
}

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
append_vcf_record(std::string & out_string, 
        std::string & out_string_pass, 
        VcStats & vc_stats,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet, 
        const char *tname, 
        unsigned int refpos, 
        const AlignmentSymbol symbol, 
        bcfrec::BcfFormat & fmtvar, 
        const std::string & refstring,
        const unsigned int extended_inclu_beg_pos, 
        const double vcfqual_thres,
        const bool should_output_all, 
        const bool should_let_all_pass,
        auto & tki, 
        const bool prev_is_tumor,
        unsigned int homref_gt_phred,
        const unsigned int uni_bias_thres,
        const bcf_hdr_t *g_bcf_hdr, 
        const bool is_tumor_format_retrieved,
        const unsigned int highqual_thres,
        double highqual_min_ratio,
        const double any_mul_contam_frac,
        const double t2n_mul_contam_frac,
        const double t2n_add_contam_frac,
        const double t2n_add_contam_transfrac,
        const std::string & repeatunit, 
        const unsigned int repeatnum,
        // const bool is_proton,
        // const unsigned int maxMQ,
        const unsigned int central_readlen,
        const unsigned int phred_triallelic_indel,
        const unsigned int phred_max_sscs,
        const unsigned int phred_max_dscs,
        const double phred_pow_sscs_origin,
        const double phred_pow_dscs_origin,
        const unsigned int vad_thres,
        const bool is_somatic_snv_filtered_by_any_nonref_germline_snv,
        const bool is_somatic_indel_filtered_by_any_nonref_germline_indel,
        // const double illumina_BQ_sqr_coef,
        // const double phred_varcall_err_per_map_err_per_base,
        const double powlaw_exponent,
        const double powlaw_anyvar_base,
        const double syserr_maxqual,
        const double syserr_norm_devqual,
        // const auto phred_umi_dimret_qual,
        const auto phred_umi_dimret_mult,
        const auto bitflag_InDel_penal_t_UMI_n_UMI,
        uint64_t haplo_in_diplo_allele_perc,
        uint64_t diplo_oneside_posbias_perc,
        uint64_t diplo_twoside_posbias_perc,
        uint64_t haplo_oneside_posbias_perc,
        uint64_t haplo_twoside_posbias_perc,
        const double phred_snv_to_indel_ratio,
        const RegionalTandemRepeat & rtr1,
        const RegionalTandemRepeat & rtr2,
        // const auto & indelbdepth,
        const auto & indelstring,
        const unsigned int  specialflag) {
    
    const bcfrec::BcfFormat & fmt = fmtvar; 
    
    assert(repeatunit.size() > 0);
    assert(repeatnum > 0);
    assert(refpos >= extended_inclu_beg_pos);
    assert(refpos - extended_inclu_beg_pos < refstring.size());
    
    const bool is_rescued = (tki.DP > 0);
    if (prev_is_tumor && (!is_rescued)) { return -1; }
    const unsigned int regionpos = refpos - extended_inclu_beg_pos;
    const char *altsymbolname = SYMBOL_TO_DESC_ARR[symbol];
    std::string vcfref;
    std::string vcfalt;
    unsigned int vcfpos;

    // std::string indelstring;
    const bool isInDel = (isSymbolIns(symbol) || isSymbolDel(symbol));
    if (isInDel) {
        vcfpos = refpos; // refpos > 0?
        vcfref = (regionpos > 0 ? refstring.substr(regionpos-1, 1) : "n");
        vcfalt = vcfref;
        //indelstring = indel_get_majority(fmt, prev_is_tumor, tki, 
        //    is_rescued, tname, refpos, symbol); 
        fmtvar.gapDP4 = {0, 0, 0, 0};
        for (size_t si = 0; si < fmt.gapSeq.size(); si++) {
            if (fmt.gapSeq[si] == indelstring) {
                fmtvar.gapDP4[2] += fmt.gapbAD1[si];
                fmtvar.gapDP4[3] += fmt.gapcAD1[si];
            }
            fmtvar.gapDP4[0] += fmt.gapbAD1[si];
            fmtvar.gapDP4[1] += fmt.gapcAD1[si];
        }
        if (indelstring.size() == 0) {
            vcfalt = altsymbolname;
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
        vcfalt = altsymbolname;
    }
    
    // const double indel_ic = ((!isInDel) ? 0.0 : (10.0/log(10.0) * log((double)MAX(indelstring.size(), 1U) / (double)(repeatunit.size() * (repeatnum - 1) + 1))));
    const bool is_novar = (symbol == LINK_M || (isSymbolSubstitution(symbol) && vcfref == vcfalt));
    if (is_novar && (!should_let_all_pass) && (!should_output_all)) {
        //return -1;
    }
    std::string vcffilter;
    if (is_novar) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::noVar]) + ";");
    }
    if (0 < vcffilter.size()) {
        vcffilter.pop_back();
    }
    
    auto aADsum = SUMVEC(fmt.aAD);
    auto aDPsum = SUMVEC(fmt.aDP);
    auto aPBADbb = SUMVEC(fmt.aPBAD);            // both side bias
    auto aPBADlb = fmt.aPBAD[0] + fmt.aPBAD[1]; // left   bias
    auto aPBADrb = fmt.aPBAD[0] + fmt.aPBAD[2]; // right  bias
    auto aPBDPbb = SUMVEC(fmt.aPBDP);            // normalization factor
    auto aPBDPlb = fmt.aPBDP[0] + fmt.aPBDP[1];
    auto aPBDPrb = fmt.aPBDP[0] + fmt.aPBDP[2];
    const bool is_highfrac_var = (100UL * aADsum > haplo_in_diplo_allele_perc * aDPsum);
    fmtvar.FT = "";
    
    bool babs_seqbias = (100UL * (1 + aADsum - aPBADbb) < diplo_twoside_posbias_perc * (1 + aPBADbb));
    bool labs_seqbias = (100UL * (1 + aADsum - aPBADlb) < diplo_oneside_posbias_perc * (1 + aPBADlb));
    bool rabs_seqbias = (100UL * (1 + aADsum - aPBADrb) < diplo_oneside_posbias_perc * (1 + aPBADrb));
    bool brel_seqbias = (100UL * (1 + aADsum - aPBADbb) * (1 + aPBDPbb) < haplo_twoside_posbias_perc * (1 + aPBADbb) * (1 + aDPsum - aPBDPbb));
    bool lrel_seqbias = (100UL * (1 + aADsum - aPBADlb) * (1 + aPBDPlb) < haplo_oneside_posbias_perc * (1 + aPBADlb) * (1 + aDPsum - aPBDPlb));
    bool rrel_seqbias = (100UL * (1 + aADsum - aPBADrb) * (1 + aPBDPrb) < haplo_oneside_posbias_perc * (1 + aPBADrb) * (1 + aDPsum - aPBDPrb));
    
    if (babs_seqbias && is_highfrac_var) { fmtvar.FT += std::string(bcfrec::FILTER_IDS[bcfrec::GPBLR2]) + ";"; }
    if (labs_seqbias && is_highfrac_var) { fmtvar.FT += std::string(bcfrec::FILTER_IDS[bcfrec::GPBL2] ) + ";"; }
    if (rabs_seqbias && is_highfrac_var) { fmtvar.FT += std::string(bcfrec::FILTER_IDS[bcfrec::GPBR2] ) + ";"; }
    if (brel_seqbias) { fmtvar.FT += std::string(bcfrec::FILTER_IDS[bcfrec::GPBLR1]) + ";"; }
    if (lrel_seqbias) { fmtvar.FT += std::string(bcfrec::FILTER_IDS[bcfrec::GPBL1] ) + ";"; }
    if (rrel_seqbias) { fmtvar.FT += std::string(bcfrec::FILTER_IDS[bcfrec::GPBR1] ) + ";"; }
    
    if ("" == fmtvar.FT) {
        fmtvar.FT = "PASS";
    } else {
        fmtvar.FT.pop_back();
    }
    
    // This hard-filtering can be done by bcftools but much more slowly
    int64_t bDP1_0 = (int64_t) fmt.bDP1[0];
    int64_t bDP1_1 = (int64_t) fmt.bDP1[1];
    int64_t cDP1_0 = (int64_t) fmt.cDP1[0];
    int64_t cDP1_1 = (int64_t) fmt.cDP1[1];
    fmtvar.FTS = "";
    fmtvar.FTSV.clear();
    unsigned int maxbias = 100;
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::DB1],    uni_bias_thres, hmean(fmt.aDB [0], bDP1_0, fmt.aDB [1], bDP1_1));
    auto db2  = fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::DB2],    uni_bias_thres, hmean(fmt.aDB [0], cDP1_0, fmt.aDB [1], cDP1_1));
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::MB1],    uni_bias_thres, hmean(fmt.bMMB[0], bDP1_0, fmt.bMMB[1], bDP1_1));
    auto mb2  = fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::MB2],    uni_bias_thres, hmean(fmt.cMMB[0], cDP1_0, fmt.cMMB[1], cDP1_1));
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::PB1L],   uni_bias_thres, hmean(fmt.bPBL[0], bDP1_0, fmt.bPBL[1], bDP1_1));
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::PB1R],   uni_bias_thres, hmean(fmt.bPBR[0], bDP1_0, fmt.bPBR[1], bDP1_1));
    auto pb2l = fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::PB2L],   uni_bias_thres, hmean(fmt.cPBL[0], cDP1_0, fmt.cPBL[1], cDP1_1));
    auto pb2r = fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::PB2R],   uni_bias_thres, hmean(fmt.cPBR[0], cDP1_0, fmt.cPBR[1], cDP1_1));
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::SB1],    uni_bias_thres, hmean(fmt.bSBR[0], bDP1_0, fmt.bSBR[1], bDP1_1));
    auto sb2  = fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::SB2],    uni_bias_thres, hmean(fmt.cSBR[0], cDP1_0, fmt.cSBR[1], cDP1_1));
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::QTD1],   uni_bias_thres, ceil(100*(MAX(fmt.bQT3[0], fmt.bQT3[1]) / (FLT_MIN+(double)MAX(fmt.bQT2[0], fmt.bQT2[1])))));
    // InDel base quality is only correctly computed for InDel ALT and not for non-InDel REF. Thus, the QTD2 bias is not applicable to InDels. TODO: fix this more upstream if possible.
    if (!isInDel) {
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::QTD2],   uni_bias_thres, ceil(100*(MAX(fmt.cQT3[0], fmt.cQT3[1]) / (FLT_MIN+(double)MAX(fmt.cQT2[0], fmt.cQT2[1])))));
    }
    auto fmt_bFA = ((double)fmt.bADR[1]) / (fmt.bDP+1.0);
    auto fmt_cFA = ((double)fmt.cADR[1]) / (fmt.cDP+1.0);
                fmtFTSupdate(maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::DBthis], uni_bias_thres, ceil(100*(fmt_cFA  / fmt_bFA)));
    if ((fmt_cFA < 0.8 || fmt_bFA < 0.8)) {
        fmtFTSupdate(
                maxbias, fmtvar.FTS, fmtvar.FTSV, bcfrec::FILTER_IDS[bcfrec::DBrest], uni_bias_thres, ceil(((1.0-fmt_cFA) / (1.0-fmt_bFA + FLT_MIN)))); 
    }
    if (0 < fmtvar.FTS.size()) {
        fmtvar.FTS.pop_back(); // not passed
    } else {
        fmtvar.FTS += "PASS";
    }
    fmtvar.EROR = {db2, MAX(pb2l, pb2r), sb2, mb2, maxbias};
    
    std::string ref_alt;
    std::string infostring = (prev_is_tumor ? "SOMATIC" : "ANY_VAR");
    
    auto ref_bias = fmt.RefBias;
    unsigned int readlen = MAX(30U, central_readlen);
    const double altmul = (double)(readlen - MIN(readlen - 30, ref_bias)) / (double)readlen;
    const double refmul = 2.0 - altmul;
    
    const auto & nfm = (prev_is_tumor ? fmt : FORMAT_UNCOV);
    if (!prev_is_tumor) {
        ref_alt = vcfref + "\t" + vcfalt;
        tki.FTS = fmt.FTS;
        tki.VAQ = (fmt.VAQ.size() > 0 ? fmt.VAQ[0] : 0);
        
        //static_assert(tki.tDPs.size() == 3);
        //static_assert(tki.nDPs.size() == 3);
        //static_assert(tki.tDPFs.size() == 3);
        //static_assert(tki.nDPFs.size() == 3);
        tki.cADTC = SUM2(fmt.cADTC);
        tki.cRDTC = SUM2(fmt.cRDTC);
        tki.cDPTC = SUM2(fmt.cDPTC);
        tki.dAD3 = fmt.dAD3;
        
        tki.DP = fmt.DP;
        tki.FA = fmt.FA;
        tki.FR = fmt.FR;
        tki.BQ = fmt.BQ;
        tki.MQ = fmt.MQ;
        tki.bDP = fmt.bDP;
        // tki.bFA = fmt.bFA;
        for (size_t i = 0; i < fmt.bAD1.size(); i++) { tki.bAD1[i] = fmt.bAD1[i]; }
        tki.autoBestAltBQ = SUM2(fmt.cAltBQ);
        tki.autoBestAllBQ = SUM2(fmt.cAllBQ);
        tki.autoBestRefBQ = SUM2(fmt.cRefBQ);
        tki.autoBestAltHD = SUM2(fmt.cAltHD);
        tki.autoBestAllHD = SUM2(fmt.cAllHD);
        tki.autoBestRefHD = SUM2(fmt.cRefHD);
        tki.cAltBQ = SUM2(fmt.cAltBQ);
        tki.cAllBQ = SUM2(fmt.cAllBQ);
        tki.cRefBQ = SUM2(fmt.cRefBQ);
        for (size_t i = 0; i < fmt.gapDP4.size(); i++) { tki.gapDP4[i] = fmt.gapDP4[i]; }
        for (size_t i = 0; i < fmt.RCC.size(); i++) { tki.RCC[i] = fmt.RCC[i]; }
        tki.GLa = fmt.GLa;
        tki.GLb = fmt.GLb;
        for (size_t i = 0; i < fmt.EROR.size(); i++) { tki.EROR[i] = fmt.EROR[i]; }
        for (size_t i = 0; i < fmt.gapbNRD.size(); i++) { tki.gapbNRD[i] = fmt.gapbNRD[i]; }
        // for (size_t i = 0; i < fmt.gapbNNRD.size(); i++) { tki.gapbNNRD[i] = fmt.gapbNNRD[i]; }
        for (size_t i = 0; i < fmt.aAD.size(); i++) { tki.aAD[i] = fmt.aAD[i]; }
    } else {
        vcfpos = (tki.ref_alt != "." ? (tki.pos + 1) : vcfpos);
        ref_alt = (tki.ref_alt != "." ? tki.ref_alt : vcfref + "\t" + vcfalt);
    }
    assert(tki.autoBestAllBQ >= tki.autoBestRefBQ + tki.autoBestAltBQ || is_novar || !fprintf(stderr, "%u >= %u + %u failed for %s:%u %s!\n", 
        tki.autoBestAllBQ, tki.autoBestRefBQ, tki.autoBestAltBQ, tname, refpos, SYMBOL_TO_DESC_ARR[symbol]));
    double vcfqual = -(double)FLT_MAX;
    bool keep_var = false;
    {
        const bool tUseHD1 = ((tki.bDP > tki.DP * highqual_min_ratio) && (tki.cADTC       > 0));
        const bool nUseHD1 = ((nfm.bDP > nfm.DP * highqual_min_ratio) && (SUM2(nfm.cADTC) > 0) && tUseHD1); 
        
        const bool tUseHD = (tUseHD1 && !isInDel);
        const bool nUseHD = (nUseHD1 && !isInDel); 
        
        double nAltBQ = SUM2(nfm.cAltBQ);
        double nAllBQ = SUM2(nfm.cAllBQ);
        double nRefBQ = SUM2(nfm.cRefBQ);
        double nAltHD = SUM2(nfm.cAltHD);
        double nAllHD = SUM2(nfm.cAllHD);
        double nRefHD = SUM2(nfm.cRefHD);
        double nAltCD = SUM2(nfm.cADTC);
        double nAllCD = SUM2(nfm.cDPTC);
        
        double nDP0 = (double)((nUseHD) ? (nAllHD) :                      (double)nfm.DP) + DBLFLT_EPS / 2.0;
        double nAD0 = (double)((nUseHD) ? (nAltHD) :             nfm.FA * (double)nfm.DP) + DBLFLT_EPS;
        if (nUseHD1 && nfm.FA > DBLFLT_EPS) {
            nDP0 = MAX3(nDP0, nAD0 / nfm.FA, nAD0 / ((nAltCD + 0.5) / (nAllCD + 1.0)));
        }
        double nRD0 = (double)((nUseHD) ? (nRefHD) :             nfm.FR * (double)nfm.DP) + DBLFLT_EPS / 2.0;
        
        double tAltBQ = tki.autoBestAltBQ;
        double tAllBQ = tki.autoBestAllBQ;
        // double tRefBQ = tki.autoBestRefBQ;
        double tAltHD = tki.autoBestAltHD;
        double tAllHD = tki.autoBestAllHD;
        // double tRefHD = tki.autoBestRefHD;
        double tAltCD = tki.cADTC;
        double tAllCD = tki.cDPTC; 
        
        double tDP0 = (double)((tUseHD) ? (tAllHD) :                      (double)tki.DP) + DBLFLT_EPS;
        double tAD0 = (double)((tUseHD) ? (tAltHD) :             tki.FA * (double)tki.DP) + DBLFLT_EPS / 2.0;
        if (tUseHD1 && tki.FA > DBLFLT_EPS) {
            tDP0 = MAX3(tDP0, tAD0 / tki.FA, tAD0 / ((tAltCD + 0.5) / (tAllCD + 1.0)));
        }
        
        // double tRD0 = (double)((tUseHD) ? (tRefHD) :             tki.FR * (double)tki.DP) + DBLFLT_EPS / 2.0; 
        
        // double nAD0a = nAD0 * ((isInDel && !nUseHD1) ? ((double)(nfm.gapDP4[2] + 1) / (double)(nfm.gapDP4[0] + 1)) : 1.0); // not used ?
        // double tAD0a = tAD0 * ((isInDel && !tUseHD1) ? ((double)(tki.gapDP4[2] + 1) / (double)(tki.gapDP4[0] + 1)) : 1.0);
        
        double nAD1 = (nUseHD ? (highqual_thres * nAltHD) : nAltBQ) + DBLFLT_EPS / 2.0;
        double nDP1 = (nUseHD ? (highqual_thres * nAllHD) : nAllBQ) + DBLFLT_EPS;
        double tAD1 = (tUseHD ? (highqual_thres * tAltHD) : tAltBQ) + DBLFLT_EPS / 2.0;
        double tDP1 = (tUseHD ? (highqual_thres * tAllHD) : tAllBQ) + DBLFLT_EPS;
        double nRD1 = (nUseHD ? (highqual_thres * nRefHD) : nRefBQ) + DBLFLT_EPS / 2.0;
        // double tRD1 = (tUseHD ? (highqual_thres * tRefHD) : tRefBQ) + DBLFLT_EPS / 2.0; 
        
        //const bool is_nonref_snp_excluded = true; // false for treating tri-allelic and tetra-allelic SNV site as potentially somatic
        //const bool is_nonref_indel_excluded = true; // false for treating tri-allelic and tetra-allelic InDel site as potentially somatic
        const bool is_nonref_germline_excluded = (isInDel ? is_somatic_indel_filtered_by_any_nonref_germline_indel : is_somatic_snv_filtered_by_any_nonref_germline_snv);
        
        // the genotype likelihoods here are special in that they can somehow normalized for the purpose of computing nonref probabilities
        // HOWEVER, no prior can be given to the raw genotype likelihoods.
        int32_t nonalt_qual = nfm.GLa[0] - MAX(nfm.GLa[1], nfm.GLa[2]); // somatic_var_pq ?
        int32_t excalt_qual = nfm.GLb[0] - MAX(nfm.GLb[1], nfm.GLb[2]); // somatic_var_pq ?
        
        // https://www.researchgate.net/figure/Distribution-of-the-Variant-Allele-Fraction-VAF-of-somatic-mutations-in-one-sample-of_fig7_278912701
        // https://onlinelibrary.wiley.com/doi/full/10.1002/humu.23674
        // http://mathworld.wolfram.com/ParetoDistribution.html
        const int median_prior_like = (int)round(log(1.0 / pow(0.5, 1.0 / powlaw_exponent)) / log(pow(10.0, 0.1)) * powlaw_exponent); // normalized with respect to median
        int32_t nonalt_tu_q = MAX(median_prior_like - (int)(MAX(tki.GLa[1], tki.GLa[2])), 0);
        int32_t excalt_tu_q = MAX(median_prior_like - (int)(MAX(tki.GLb[1], tki.GLb[2])), 0);
        
        const double tE0 = 1.0;
        const double nE0 = 1.0;
        const double tnE0 = 1.0;
        const double tnDP0ratio = (double)(tDP0 +           1) / (double)(nDP0 + 1);
        const double tnFA0 =      (double)(tAD0 + nAD0 + tnE0) / (double)(tDP0 + nDP0 + tnE0 * 2.0);
        const double tAD0pc0 = 0.5        * tnDP0ratio;
        const double tDP0pc0 = 0.5        * tnDP0ratio / tnFA0;
        const double nAD0pc0 = 0.5       ;
        const double nDP0pc0 = 0.5        / tnFA0;
        
        //const double tE1  = ((10.0/log(10.0) * log((double)(tDP0 + 2))));
        //const double nE1  = ((10.0/log(10.0) * log((double)(nDP0 + 2))));
        const double tnE1 = ((10.0/log(10.0) * log((double)(tDP0 + nDP0 + 2))));
        //const double tnDP1ratio = (double)(tDP1 + tE1        ) / (double)(nDP1 + nE1);
        const double tnFA1 =      (double)(tAD1 + nAD1 + tnE1) / (double)(tDP1 + nDP1 + tnE1 * 2.0);
        const double tAD1pc0 = 0.5 * tnE1 * tnDP0ratio       ;
        const double tDP1pc0 = 0.5 * tnE1 * tnDP0ratio / tnFA1;
        const double nAD1pc0 = 0.5 * tnE1;
        const double nDP1pc0 = 0.5 * tnE1 / tnFA1;
        
        const double pl_exponent = powlaw_exponent;
        
        // const double pl_exponent_t = powlaw_exponent;
        // const double pl_exponent_n = powlaw_exponent;
        
        // const double symbol_to_allele_frac = 1.0 - pow((isSymbolIns(symbol) ? 0.9 : (isSymbolDel(symbol) ? 0.95 : 1.0)), indelstring.size());
        
        // This is the penalty induced by EROR bias. 
        // So far I do not have enough data to verify this power-law. 
        // Also, the theory to support this power-law is rather unclear.
        // Hence, a simple hard-filter with rejection for EROR>=180 is used for now.
        // const double eror_penal = (pl_exponent) * 10.0/log(10.0) * log((double)(MAX(100, maxbias) - 100) / 50.0 + 1.0); 
        // QUESTION: 
        // 1. What if the duplication rate is in-between the one of UMI and the one of non-UMI? 
        // 2. Does EROR bias follows a power-law distribution too? A lot of true positive UMI variants are required to see the power-law if it does.
        //    If EROR bias does follow a power-law, then what is the dimension of EROR bias (4 or 3) given that the dimension of quality is equal to 3? My intuition says it's 3.
        /*
        const double pcap_tmax = powlaw_anyvar_base + (tUseHD1
                ? ((tki.dAD3 > 0)
                    ? ((double)(phred_max_dscs - phred_pow_dscs_origin))
                    : ((double)(phred_max_sscs - phred_pow_sscs_origin) 
                      * (double)(MAX(tki.DP, tki.bDP) - tki.DP) / (double)(tki.bDP + 1)
                      ))
                : 0.0);
        // prob[mapping-error] * prob[false-positive-variant-per-base-position] / num-alts-per-base-positon
        const double pcap_nmax = powlaw_anyvar_base + (nUseHD1
                ? ((nfm.dAD3 > 0)
                    ? ((double)(phred_max_dscs - phred_pow_dscs_origin))
                    : ((double)(phred_max_sscs - phred_pow_sscs_origin) 
                    * (double)(MAX(nfm.DP, nfm.bDP) - nfm.DP) / (double)(nfm.bDP + 1)
                    ))
                : 0.0);
        */
        // double t_indel_penal = 0.0;
        // const double pcap_bcq = ((isInDel || is_proton) ? 200.0 : (mathsquare(tki.BQ) * illumina_BQ_sqr_coef)); // based on heuristics
        // const double pcap_tmq = MIN((double)maxMQ, tki.MQ) + phred_varcall_err_per_map_err_per_base; // bases on heuristics
        // const double pcap_tmq2 = MAX(0.0, tki.MQ - 10.0/log(10.0) * log((double)(tDP0 + tDP0pc0) / (double)(tAD0+tAD0pc0))) * (double)tAD0; // readjustment by MQ
        // const double tn_t_altmul = (tUseHD1 ? 1.0 : altmul);
        // const double tn_t_refmul = (tUseHD1 ? 1.0 : refmul);
        // const double tn_tpo1q = 10.0 / log(10.0) * log((double)(tAD0 / tn_t_altmul + tE0) / ((tDP0 - tAD0) / tn_t_refmul + tAD0 / tn_t_altmul + 2.0 * tE0)) * (pl_exponent_t) + (pcap_tmax);
        // const double tn_tsamq = 40.0 * pow(0.5, MAX((double)tAD0, ((double)(tki.bDP + 1)) / ((double)(tki.DP + 1)) - 1.0));
        const double tn_tra1q = ((double)tki.VAQ); // non-intuitive prior like
        // double _tn_tpo2q = tn_tpo1q; //MIN4(tn_tpo1q, pcap_tmq, pcap_bcq, pcap_tmq2);
        // double _tn_tra2q = MAX(0.0, tn_tra1q /*tn_t_altmul*/ - tn_tsamq); //  + somatic_var_pq); // tumor  BQ bias is stronger as raw variant (biased to high BQ) is called   from each base
        
        /*
        if (isInDel) {
            if (!tUseHD1) {
                int n_str_units = (isSymbolIns(symbol) ? 1 : (-1)) * (int)(indelstring.size() / repeatunit.size());
                t_indel_penal = penal_indel_2(tAD0a, n_str_units, tki.RCC, phred_triallelic_indel);
                _tn_tra2q *= ((double)tAD0a + (double)(tAD0 - tAD0a) * symbol_to_allele_frac) / ((double)tAD0); // ad-hoc adjustment
            }
        }
        
        double t_penal_by_nearby_indel = 0.0; 
        unsigned int indelmax = MAX(tki.gapbNRD[0] + tki.gapbNRD[1], tki.gapbNRD[2] + tki.gapbNRD[3]);
        if ((SUM2(tki.bAD1) < (int)indelmax) && (indelmax > 0)) {
            uint64_t max_pb2 = tki.EROR[1];
            double bAD1frac = (double)(indelmax - SUM2(tki.bAD1)) / (double)indelmax;
            double sqr_pb2 = bAD1frac * (double)(max_pb2 * MIN(200UL, max_pb2));
            t_penal_by_nearby_indel = 10.0/log(10.0) *log(MAX(1.0, sqr_pb2 / 1e4)) * pl_exponent;
        }
        if ((isInDel) && !tUseHD1) {
            _tn_tpo2q -= MAX(t_indel_penal, t_penal_by_nearby_indel);
        }
        if ((!isInDel) && !tUseHD1) {
            _tn_tpo2q -= t_penal_by_nearby_indel;
        }
        */
        // const double tn_tpowq = _tn_tpo2q;
        // const double tn_trawq = _tn_tra2q;
        // QUESTION: why BQ-bias generations are discrete? Because of the noise with each observation of BQ? why indel generations are discrete?
        
        // const double pcap_nbq = ((isInDel || is_proton) ? 200.0 : mathsquare(nfm.BQ) * illumina_BQ_sqr_coef); // based on heuristics
        // const double pcap_nmq = MIN(maxMQ, nfm.MQ) * phred_varcall_err_per_map_err_per_base; // based on heuristics
        // const double pcap_nmq2 = (nfm.MQ - 10.0/log(10.0) * log((double)(nDP0 + nDP0pc0) / (double)(nAD0+nAD0pc0))) * (double)nAD0; // readjsutment by MQ
        // const double tn_n_altmul = (nUseHD1 ? 1.0 : altmul);
        // const double tn_n_refmul = (nUseHD1 ? 1.0 : refmul);
        // const double tn_npo1q = 10.0 / log(10.0) * log((double)(nAD0 / tn_n_altmul + nE0) / ((nDP0 - tAD0) / tn_n_refmul + tAD0 / tn_n_altmul + 2.0 * nE0)) * (pl_exponent_n) + (pcap_nmax);
        //const double tn_nsamq = 40.0 * pow(0.5, MAX((double)nAD0, ((double)(nfm.bDP + 1)) / ((double)(nfm.DP + 1)) - 1.0));
        //const double tn_nra1q = (double)(nfm.VAQ.size() > 0 ? (fmt.VAQ[0]) : 0);
        // double _tn_npo2q = tn_npo1q; // MIN4(tn_npo1q, pcap_nmq, pcap_nbq, pcap_nmq2);
        // double _tn_nra2q = MAX(0.0, tn_nra1q /*tn_n_altmul*/ - tn_nsamq/4.0); // + somatic_var_pq); // normal BQ bias is weaker   as res variant (biased to high BQ) is filtered from each variant
        // if (isInDel) {
        //    if (!nUseHD1) {
        //        _tn_nra2q *= ((double)nAD0a + (double)(nAD0 - nAD0a) * symbol_to_allele_frac) / ((double)nAD0);
        //    }
        // }
        // const double tn_npowq = _tn_npo2q;
        // const double tn_nrawq = _tn_nra2q;
        
        double symfrac = 1.0; // degenerated
        
        const double t2n_rawq0 = ((true || nDP0 <= tDP0) // CHECK: check if the symmetry makes sense. If it does then enable the else part
            ? calc_binom_10log10_likeratio((double)(tDP0 - tAD0) / (double)tDP0, (double)(nDP0 - nAD0),                (double)(nAD0)                      )
            : calc_binom_10log10_likeratio((double)       (nAD0) / (double)nDP0, (double)(tAD0),                       (double)(tDP0 - tAD0)               )); 
        const double t2n_rawq1 = ((true || nDP0 <= tDP0) // CHECK: check if the symmetry makes sense. If it does then enable the else part
            ? calc_binom_10log10_likeratio(        (tDP1 - tAD1) /         tDP1,         (nDP1 - nAD1) / nDP1 * nDP0,          (nAD1 / nDP1) * nDP0)
            : calc_binom_10log10_likeratio(               (nAD1) /         nDP1,         (tAD1 / tDP1) * tDP0,                 (tDP1 - tAD1) / tDP1 * tDP0)); 
        const double t2n_rawq = (isInDel ? t2n_rawq0 : t2n_rawq1 / MAX(1.0, 0.5 + 0.5*tki.EROR[3]/100.0));
        
        auto n3D0 = nDP0 - MIN(nDP0 - (DBLFLT_EPS / 2.0), nAD0 + nRD0);
        auto n3D1 = nDP1 - MIN(nDP1 - (DBLFLT_EPS / 2.0), nAD1 + nRD1);
        const double t3n_rawq0 = ((true || nDP0 <= tDP0) // CHECK: check if the symmetry makes sense. If it does then enable the else part
            ? calc_binom_10log10_likeratio((double)(tDP0 - tAD0) / (double)tDP0, (double)(nDP0 - n3D0),                (double)(n3D0)                      )
            : calc_binom_10log10_likeratio((double)       (n3D0) / (double)nDP0, (double)(tAD0),                       (double)(tDP0 - tAD0)               )); 
        const double t3n_rawq1 = ((true || nDP0 <= tDP0) // CHECK: check if the symmetry makes sense. If it does then enable the else part
            ? calc_binom_10log10_likeratio(        (tDP1 - tAD1) /         tDP1,         (nDP1 - n3D1) / nDP1 * nDP0,          (n3D1 / nDP1) * nDP0)
            : calc_binom_10log10_likeratio(               (n3D1) /         nDP1,         (tAD1 / tDP1) * tDP0,                 (tDP1 - tAD1) / tDP1 * tDP0)); 
        const double t3n_rawq = (isInDel ? t3n_rawq0 : t3n_rawq1 / MAX(1.0, 0.5 + 0.5*tki.EROR[3]/100.0));
        
        const double t2n_or0 = ((double)(tAD0 + symfrac * tAD0pc0) / (double)(tDP0 - tAD0 + symfrac * (tDP0pc0 - tAD0pc0))) 
                             / ((double)(nAD0 +           nAD0pc0) / (double)(nDP0 - nAD0 +           (nDP0pc0 - nAD0pc0)));
        const double t2n_or1 = ((double)(tAD1 + symfrac * tAD1pc0) / (double)(tDP1 - tAD1 + symfrac * (tDP1pc0 - tAD1pc0))) 
                             / ((double)(nAD1 +           nAD1pc0) / (double)(nDP1 - nAD1 +           (nDP1pc0 - nAD1pc0))); 
        
        const double t2n_po0q0 = 10.0/log(10.0) * (1.0+log(symfrac*symfrac)/10.0) * pl_exponent * log(t2n_or0 /symfrac);
        const double t2n_po0q1 = 10.0/log(10.0) * (1.0+log(symfrac*symfrac)/10.0) * pl_exponent * log(t2n_or1 /symfrac);
        const double t2n_po2q0 = MIN(MAX(-syserr_maxqual, t2n_po0q0), syserr_maxqual + 10);
        const double t2n_po2q1 = MIN(MAX(-syserr_maxqual, t2n_po0q1), syserr_maxqual + 10);
        const double t2n_powq0 = MIN(MAX(-syserr_maxqual, t2n_po0q0), syserr_maxqual);
        const double t2n_powq1 = MIN(MAX(-syserr_maxqual, t2n_po0q1), syserr_maxqual);
        const double t2n_powq = (isInDel ? t2n_powq0 : t2n_powq1);
        const double t2t_powq = syserr_maxqual;
        
        const double t3n_or0 = ((double)(tAD0 + symfrac * tAD0pc0) / (double)(tDP0 - tAD0 + symfrac * (tDP0pc0 - tAD0pc0))) 
                             / ((double)(n3D0 +           nAD0pc0) / (double)(nDP0 - n3D0 +           (nDP0pc0 - nAD0pc0)));
        const double t3n_or1 = ((double)(tAD1 + symfrac * tAD1pc0) / (double)(tDP1 - tAD1 + symfrac * (tDP1pc0 - tAD1pc0))) 
                             / ((double)(n3D1 +           nAD1pc0) / (double)(nDP1 - n3D1 +           (nDP1pc0 - nAD1pc0))); 
        const double t3n_po0q0 = 10.0/log(10.0) * (1.0+log(symfrac*symfrac)/10.0) * pl_exponent * log(t3n_or0 /symfrac);
        const double t3n_po0q1 = 10.0/log(10.0) * (1.0+log(symfrac*symfrac)/10.0) * pl_exponent * log(t3n_or1 /symfrac);
        const double t3n_po2q0 = MIN(MAX(-syserr_maxqual, t3n_po0q0), syserr_maxqual + 10);
        const double t3n_po2q1 = MIN(MAX(-syserr_maxqual, t3n_po0q1), syserr_maxqual + 10);
        
        double t2n_contam_q = MIN(calc_binom_10log10_likeratio(t2n_add_contam_transfrac, 
                (isSymbolSubstitution(symbol) ? ((SUM2(nfm.cAltBQ) + DBL_MIN) / (double)(SUM2(nfm.cAllBQ) + 2.0 * DBL_MIN)) : nfm.FA) * (double)nfm.DP,
                (isSymbolSubstitution(symbol) ? ((    (tki.cAltBQ) + DBL_MIN) / (double)(    (tki.cAllBQ) + 2.0 * DBL_MIN)) : tki.FA) * (double)tki.DP), 200.0);
        // double t2n_limq = 200.0; 
        
        /*
        if (isInDel && prev_is_tumor) {
            // This is an ad-hoc adjustment to lower TLODQ when there is contamination, and to further lower TLODQ if the variant has sequencing position bias (aka IGV position bias)
            // TODO: rewrite this piece of code with the sequencing position bias redefined in terms of EROR instead of simple hard filter
            // 3.0 is the tolerance for errors in contamination estimation
            double t2n_add_coontam_transfrac_low = t2n_add_contam_transfrac / 3.0;
            double t2n_nonalt_frac = ((1.0 - nfm.FA) * (double)nfm.DP + nE0) / ((1.0 - tki.FA) * (double)tki.DP + tE0);
            double max_tn_ratio = (1.0 - t2n_add_coontam_transfrac_low) / t2n_add_coontam_transfrac_low * MIN(2.0, t2n_nonalt_frac);
            t2n_limq = 10.0/log(10.0) * log(MAX(max_tn_ratio, 1.0)) - ((LINK_D1 == symbol) ? 6 : 0); 
            uint64_t edge_perc = 100UL * (tki.aPBAD[0] + MAX(tki.aPBAD[1], tki.aPBAD[2]));
            if (edge_perc >= 75UL * SUMVEC(tki.aAD)) {
                t2n_limq -= MIN(MAX(0, fmt.VAQ[0]), 60.0);
            }
        }
        */
#if 0   // This if branch of macro is disabled
        // The following line works well in practice but has no theory supporting it yet.
        double t2n_syserr_q0 = (isInDel ? 0.0 : MIN(MAX(0.0, MIN(tn_npowq, tn_nrawq)) * MIN(1.0, 4.0 / mathsquare(t2n_or1 + 1.0)), SYS_QMAX)); // 50.0
        // double t2n_syserr_q = (isInDel ? 0.0 : MIN(MAX(0.0, MIN3(tn_npowq, tn_nrawq, 45.0) - 20.0 * mathsquare(MAX(0.0, t2n_or1 - 1.0))), 2.0*SYS_QMAX)); // 50.0
#else  
        double t2n_reward_q0 = max_min01_sub02(CENTER(t2n_rawq, t2t_powq), CENTER(t2n_rawq, t2n_powq), t2n_contam_q);
        double t2n_syserr_q0 = MIN3(syserr_maxqual, fmt.VAQ[0] - syserr_norm_devqual * mathsquare(MAX(0.0, t2n_or1 - 1.0)),
                ((double)nAD0 + 1.0) * 10.0 / log(10.0) * log((double)nDP0 + 1.0));
#endif
        // double t2n_reward_q0 = t2n_finq; // MIN(t2n_finq, t2n_limq); 
        const bool use_reward = is_bitflag_checked(bitflag_InDel_penal_t_UMI_n_UMI, isInDel, false, tUseHD1, nUseHD1);
        double t2n_reward_q = (use_reward ? MAX(0.0, t2n_reward_q0) : 0.0);
        const bool use_penalt = is_bitflag_checked(bitflag_InDel_penal_t_UMI_n_UMI, isInDel, true,  tUseHD1, nUseHD1);
        double t2n_syserr_q = (use_penalt ? MAX(0.0, t2n_syserr_q0) : 0.0);
        
        double t_base_q = (double)tki.VAQ;
        // double t_base_q = MIN(tn_trawq, tn_tpowq + (double)indel_ic);
        // fmtvar.ALODQ = MAX(0, t_base_q);
        //if (is_novar && (!should_let_all_pass) && (!should_output_all)) {
        //    return -1;
        //}
        // double t2n_infodist = calc_binom_10log10_likeratio(tAD0 / (tDP0 + DBL_EPSILON), nAD0, (nDP0 + DBL_EPSILON));
        double tlodq = t_base_q + t2n_reward_q - MIN(t2n_contam_q, t2n_syserr_q);
                
        //double n2t_red_qual = MIN(tn_npowq, tn_nrawq + (double)indel_ic) * MIN(1.0, n2t_or1) * MIN(1.0, n2t_or1); // / (t2n_or1 * t2n_or1);
        //double n2t_orr_qual = MIN(tn_npowq, tn_nrawq + (double)indel_ic) * MIN(1.0, n2t_or1);
        //double n2t_or2_qual = MIN(tn_npowq, tn_nrawq + (double)indel_ic) * MIN(1.0, n2t_or1/2.0);
        
        const int32_t noisy_germ_phred = (isInDel ? 5 : 0); // 5; // likelihood that tumor is alt1/alt2 hetero and normal is ref/alt1 hetero, normalized to zero for SNPs.
        
        std::array<double, N_MODELS> testquals = {{0}};
        unsigned int tqi = 0;

        double a_no_alt_qual = add01_between_min01_max01(nonalt_qual, nonalt_tu_q)
                - MIN(t2n_contam_q, t2n_syserr_q)
                + MAX(0, CENTER(t2n_rawq, (isInDel ? t2n_po2q0 : t2n_po2q1) - 10)) // ad-hoc: the 10 is kind of arbitrary.
                ;
        double a_ex_alt_qual = add01_between_min01_max01(excalt_qual, excalt_tu_q) + noisy_germ_phred
                + MAX(0, CENTER(t3n_rawq, (isInDel ? t3n_po2q0 : t3n_po2q1) - 10)) 
                ;
        const int32_t a_nogerm_q = homref_gt_phred + (is_nonref_germline_excluded ? 
                MIN(a_no_alt_qual, a_ex_alt_qual) : a_no_alt_qual);
        
        testquals[tqi++] = MIN(tlodq, (double)a_nogerm_q); // - 5.0;
        
        // testquals[tqi++] = MIN(tn_trawq - tn_nrawq + 0       , tn_tpowq - MAX(0.0, tn_npowq - tvn_or_q) + tvn_powq);
        // testquals[tqi++] = MIN(tn_trawq - tn_nrawq + tvn_rawq, tn_tpowq - MAX(0.0, tn_npowq - tvn_or_q) + tvn_powq); 
        // testquals[tqi++] = MIN(tn_trawq, tn_tpowq) + MIN(tvn_rawq, tvn_powq) - MAX(0.0 , MIN(tn_nrawq, tn_npowq) - tvn_or_q); // 6
        // testquals[tqi++] = MIN(tn_trawq, tn_tpowq)               + tvn_powq  - MIN(MIN(tn_nrawq, tn_npowq           ), contam_phred); // 20
        
        // FIXME: probabilities of germline polymorphism and somatic mutation at a locus are both strongly correlated with the STR pattern for InDels
        //        here we assumed that the likelihood of germline polymorphism is proportional to the likelihood of somatic mutation.
        //        however, in practice the two likelihoods may be different from each other.
        const int MODEL_SEP_1 = 1;
        vcfqual = calc_non_negative(prev_is_tumor ? ((double)testquals[0]) : ((double)tlodq));
        if (prev_is_tumor) {
            unsigned int median_intq = (unsigned int)MIN(MAX(0, (int)vcfqual), VCFQUAL_NUM_BINS - 1);
            vc_stats.vcfqual_to_count[median_intq].nvars+= 1;
            vc_stats.vcfqual_to_count[median_intq].tuDP += tDP0;
            vc_stats.vcfqual_to_count[median_intq].tuAD += tAD0;
            vc_stats.vcfqual_to_count[median_intq].noDP += nDP0;
            vc_stats.vcfqual_to_count[median_intq].noAD += nAD0;
        }
        
        keep_var = ((vcfqual >= vcfqual_thres || (tki.DP * tki.FA) >= vad_thres) && !is_novar); 
        if ((!keep_var) && (!should_output_all) && (!should_let_all_pass)) {
            return -2;
        }
        
#if (1 < N_MODELS)
        for (int i = 0; i < N_MODELS; i++) {
            infostring += std::string(";TQ") + std::to_string(i) + "=" + std::to_string(testquals[i]); 
        }
#endif
        infostring += std::string(";SomaticQ=") + std::to_string(testquals[0]);
        infostring += std::string(";TLODQ=")    + std::to_string(tlodq);
        infostring += std::string(";NLODQ=")    + std::to_string(a_nogerm_q);
        // infostring += std::string(";QUAL1=")    + std::to_string(calc_non_negative(tki.VAQ));
        //unsigned int tlodq1 = (maxbias < uni_bias_thres ?(10* tlodq) : (10*tlodq*10 - 10*60));
        //infostring += std::string(";TLODQ1=")  + std::to_string((int)(tlodq1))
        infostring += std::string(";REFQs=") + string_join(std::array<std::string, 4>({
                std::to_string(nonalt_qual), std::to_string(nonalt_tu_q),
                std::to_string(excalt_qual), std::to_string(excalt_tu_q) 
        }));
        
        infostring += std::string(";TNQs=") + string_join(std::array<std::string, 4>({
                std::to_string((int)t2n_reward_q0),
                // std::to_string(t2n_limq),
                std::to_string((int)t2n_syserr_q0),
                std::to_string((int)t2n_contam_q)
        }));
        infostring += std::string(";TNORQs=") + string_join(std::array<std::string, 2+2*0>({
                std::to_string(t2n_or0),     std::to_string(t2n_or1) // ,
                // std::to_string(n2t_or0),     std::to_string(n2t_or1)
        }));
        infostring += std::string(";TN1Qs=")  + string_join(std::array<std::string, 4>({
                std::to_string(t2n_po0q0), std::to_string(t2n_rawq0), 
                std::to_string(t2n_po0q1), std::to_string(t2n_rawq1)      
        }));
        infostring += std::string(";TN2Qs=") + string_join(std::array<std::string, 4>({
                std::to_string(t3n_rawq0)  , std::to_string(t3n_rawq1),
                std::to_string(t3n_po0q0)  , std::to_string(t3n_po0q1),
        }));
        /*
        infostring += std::string(";TAQs=") + string_join(std::array<std::string, 3>({
            std::to_string(t_base_q),
            // std::to_string(t_indel_penal), std::to_string(t_penal_by_nearby_indel),
        }));
        infostring += std::string(";TSQs=") + string_join(std::array<std::string, 4>({
                std::to_string(tn_tpowq)  , std::to_string(_tn_tra2q),
                // std::to_string(tn_tpo1q)  , 
                std::to_string(tn_tra1q), std::to_string(tn_tsamq), // std::to_string(pcap_tmax)
        }));
        infostring += std::string(";NSQs=") + string_join(std::array<std::string, 4>({
                std::to_string(tn_npowq)  , std::to_string(_tn_nra2q),
                // std::to_string(tn_npo1q)  , 
                std::to_string(tn_nra1q), std::to_string(tn_nsamq),
        }));
        */ 
        infostring += std::string(";tFA=") + std::to_string(tki.FA);
        // infostring += std::string(";tFR=") + std::to_string(tki.FR);
                
            infostring += std::string(";tADR=") + string_join(std::array<std::string, 2>({
                std::to_string((int)(tki.DP * tki.FR + 0.5)),
                std::to_string((int)(tki.DP * tki.FA + 0.5))}));
            infostring += std::string(";tDP=") + std::to_string(tki.DP);
        if (prev_is_tumor) {
            infostring += std::string(";nADR=") + string_join(std::array<std::string, 2>({
                std::to_string((int)(fmt.DP * fmt.FR + 0.5)),
                std::to_string((int)(fmt.DP * fmt.FA + 0.5))}));
            infostring += std::string(";nDP=") + std::to_string(fmt.DP);
        }
            infostring += std::string(";tADCR=") + string_join(std::array<std::string, 2>({
                std::to_string(tki.cRDTC),
                std::to_string(tki.cADTC)}));
            infostring += std::string(";tDPC=") + std::to_string(tki.cDPTC);
        if (prev_is_tumor) {
            infostring += std::string(";nADCR=") + string_join(std::array<std::string, 2>({
                std::to_string(SUM2(fmt.cRDTC)),
                std::to_string(SUM2(fmt.cADTC))}));
            infostring += std::string(";nDPC=") + std::to_string(SUM2(fmt.cDPTC));
        }
        infostring += std::string(";tADBQR=") + string_join(std::array<std::string, 2>({
            std::to_string(tki.autoBestRefBQ),
            std::to_string(tki.autoBestAltBQ)}));
        infostring += std::string(";tDPBQ=") + std::to_string(tki.autoBestAllBQ);
        //infostring += std::string(";tAltBQ2=") + std::to_string(tki.cAltBQ2);
        //infostring += std::string(";tAllBQ2=") + std::to_string(tki.cAllBQ2);
        //infostring += std::string(";tRefBQ2=") + std::to_string(tki.cRefBQ2);
        infostring += std::string(";tADHDR=") + string_join(std::array<std::string, 2>({
            std::to_string(tki.autoBestRefHD),
            std::to_string(tki.autoBestAltHD)}));
        infostring += std::string(";tDPHD=") + std::to_string(tki.autoBestAllHD);
        infostring += std::string(";tFTS=") + tki.FTS;
        // infostring += std::string(";tcHap=") + tki.cHap; // tcHap is linked with tFTS
        infostring += std::string(";tbDP=") + std::to_string(tki.bDP);
        infostring += std::string(";tMQ=") + std::to_string(tki.MQ);
        infostring += std::string(";tEROR=") + other_join(tki.EROR);
        infostring += std::string(";tgapDP4=") + other_join(tki.gapDP4);
        infostring += std::string(";tRCC=") + other_join(tki.RCC);
    } 
    infostring += std::string(";RU=") + repeatunit + ";RC=" + std::to_string(repeatnum);
    auto rtr1_chrpos = ((0 == rtr1.tracklen) ? 0 : (extended_inclu_beg_pos + rtr1.begpos));
    auto rtr2_chrpos = ((0 == rtr2.tracklen) ? 0 : (extended_inclu_beg_pos + rtr2.begpos));
    infostring += std::string(";R3X2=") + other_join(std::array<unsigned int, 6>{{rtr1_chrpos, rtr1.tracklen, rtr1.unitlen, rtr2_chrpos, rtr2.tracklen, rtr2.unitlen}}, ",");
    if (0 == vcffilter.size()) {
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
            vcffilter += ("PASS");
        }
    }
    
    if (0 < vcffilter.size() && ';' == vcffilter[vcffilter.size()-1]) {
        vcffilter.pop_back();
    }
    
    if (keep_var || should_let_all_pass) {
        out_string_pass += 
                std::string(tname) + "\t" + std::to_string(vcfpos) + "\t.\t" + ref_alt + "\t" 
                + std::to_string(vcfqual) + "\t" + vcffilter + "\t" + infostring + "\t" + bcfrec::FORMAT_STR_PER_REC + "\t";
        bcfrec::streamAppendBcfFormat(out_string_pass, fmt);
        out_string_pass += ((g_bcf_hdr != NULL && is_tumor_format_retrieved) ? bcf1_to_string(g_bcf_hdr, tki.bcf1_record) : std::string("")) + "\n";
    }
    
    if (should_output_all) {
        out_string += 
                std::string(tname) + "\t" + std::to_string(vcfpos) + "\t.\t" + ref_alt + "\t"
                + std::to_string(vcfqual) + "\t" + vcffilter + "\t" + infostring + "\t" + bcfrec::FORMAT_STR_PER_REC + "\t";
        bcfrec::streamAppendBcfFormat(out_string, fmt);
        out_string += ((g_bcf_hdr != NULL && is_tumor_format_retrieved) ? bcf1_to_string(g_bcf_hdr, tki.bcf1_record) : std::string("")) + "\n";
    }

    return 0;
}

#endif
