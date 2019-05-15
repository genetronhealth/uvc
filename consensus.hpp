#ifndef consensus_INCLUDED
#define consensus_INCLUDED

#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "htslib/sam.h"
#include "htslib/hts.h"

#include "common.h"
#include "conversion.hpp"
#include "logging.hpp"
#include "utils.hpp"

#include "bcf_formats.step1.c"

// #define bam_phredi(b, i) (phred2bucket(bam_get_qual((b))[(i)]))
#define bam_phredi(b, i) (bam_get_qual((b))[(i)])

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
    BASE_NN, //  = 5, // ambigous base after collapsing different reads
    //BASE_P = 6; // padded in deleted sequence
    LINK_M,  //   = 6, // absence of any gap
    LINK_D3P,// = 7, // deletion of length 3 or plus
    LINK_D2, //  = 8,  // deletion of length 2
    LINK_D1, //  = 9,
    LINK_I3P, //  = 10, // insertion of length 1 // where the inserted sequence is not a repeat
    LINK_I2, //  = 11, 
    LINK_I1,// = 12, 
    LINK_NN, //  = 13, // ambiguous link between bases
    // LINK_P = 13; // padded in deleted sequence
    END_ALIGNMENT_SYMBOLS,
};

struct _CharToSymbol {
    AlignmentSymbol data[256];
    _CharToSymbol() {
        for (int i = 0; i < 256; i++) {
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

const char* SYMBOL_TO_DESC_ARR[] = {
    [BASE_A] = "A", [BASE_C] = "C", [BASE_G] = "G", [BASE_T] = "T", [BASE_N] = "N",
    [BASE_NN] = "<BN>", 
    [LINK_M] = "<LR>", 
    [LINK_D3P] = "<LD3P>", [LINK_D2] = "<LD2>", [LINK_D1] = "<LD1>",
    [LINK_I3P] = "<LI3P>", [LINK_I2] = "<LI2>", [LINK_I1] = "<LI1>",
    [LINK_NN] = "<LN>",
};

bool areSymbolsMutated(AlignmentSymbol ref, AlignmentSymbol alt) {
    if (alt <= BASE_NN) {
        return ref != alt && ref < BASE_N && alt < BASE_N;
    } else {
        return alt != LINK_M && alt != LINK_NN;
    }
};

constexpr bool isSymbolIns(const AlignmentSymbol symbol) {
    if (LINK_I3P == symbol || LINK_I2 == symbol || LINK_I1 == symbol) {
        return true;
    } else {
        return false;
    }
}

constexpr bool isSymbolDel(const AlignmentSymbol symbol) {
     if (LINK_D3P == symbol || LINK_D2 == symbol || LINK_D1 == symbol) {
        return true;
    } else {
        return false;
    }  
}

const AlignmentSymbol insLenToSymbol(unsigned int len) {
    assert(len > 0);
    return 1 == len ? LINK_I1 : (2 == len ? LINK_I2 : LINK_I3P);
}

const AlignmentSymbol delLenToSymbol(unsigned int len) {
    assert(len > 0);
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

bool isSymbolSubstitution(AlignmentSymbol symbol) {
    return (SYMBOL_TYPE_TO_INCLU_BEG[BASE_SYMBOL] <= symbol && symbol <= SYMBOL_TYPE_TO_INCLU_END[BASE_SYMBOL]);
}

template <class T>
class TDistribution {
protected:
    std::array<T, NUM_ALIGNMENT_SYMBOLS> symbol2data; // [NUM_ALIGNMENT_SYMBOLS];
};

typedef \
uint32_t \
molcount_t;

// edge distance bucket
#define NUM_EDBUCKS (11+1) // (11*12/2+1)
// #define EDBUCK_SIZE 4

const std::array<unsigned int, (NUM_EDBUCKS*(NUM_EDBUCKS-1)/2+1+1)> 
DIST_TO_EDBUCK = {0, 
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
};
const std::array<unsigned int, NUM_EDBUCKS>
EDBUCK_TO_DIST = {
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
};

unsigned int pos2edbuck(unsigned int pos) {
    //return MIN(pos / EDBUCK_SIZE, NUM_EDBUCKS - 1);
    return DIST_TO_EDBUCK[MIN(pos, DIST_TO_EDBUCK.size()-1)];
}

unsigned int edbuck2pos(unsigned int edbuck) {
    assert(EDBUCK_TO_DIST.size() > edbuck);
    // return edbuck * EDBUCK_SIZE;
    return EDBUCK_TO_DIST[edbuck];
}

typedef std::array<molcount_t, NUM_BUCKETS> Bucket2Count;
typedef std::array<molcount_t, NUM_EDBUCKS> Bucket2CountEdgeDist;

const int _print_Bucket2CountEdgeDist(const Bucket2CountEdgeDist & arg) {
    for (size_t i = 0; i < NUM_EDBUCKS; i++) {
        LOG(logINFO) << arg.at(i) << "\t";
    }
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
        assert(sizeof(TB2C) * NUM_ALIGNMENT_SYMBOLS == sizeof(this->symbol2data) || !fprintf(stderr, "%d * %d != %d\n", sizeof(TB2C), NUM_ALIGNMENT_SYMBOLS, sizeof(this->symbol2data)));
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

template <class TInteger>
class GenericSymbol2Count : TDistribution<TInteger> {
public:
    TInteger
    getSymbolCount(AlignmentSymbol symbol) const {
        return this->symbol2data[symbol];
    };
    
    template <ValueType TUpdateType = SYMBOL_COUNT_SUM>
    int 
    incSymbolCount(const AlignmentSymbol symbol, const TInteger increment) {
        static_assert(BASE_QUALITY_MAX == TUpdateType || SYMBOL_COUNT_SUM == TUpdateType);
        if (SYMBOL_COUNT_SUM == TUpdateType) {
            this->symbol2data[symbol] += increment;
        } else {
            this->symbol2data[symbol] = MAX(this->symbol2data[symbol], increment);
        }
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
            assert(false);
        }
    };

    int
    _fillConsensusCounts(
            AlignmentSymbol & count_argmax, unsigned int & count_max, unsigned int & count_sum,
            AlignmentSymbol incluBeg, AlignmentSymbol incluEnd) const {
        assert (incluBeg <= incluEnd);

        count_argmax = incluEnd; // AlignmentSymbol(NUM_ALIGNMENT_SYMBOLS); // TODO:investigate further
        count_max = 0;
        count_sum = 0;
        for (AlignmentSymbol symb = incluBeg; symb <= incluEnd; symb = AlignmentSymbol(((unsigned int)symb) + 1)) {
            if (count_max < this->symbol2data[symb]) {
                count_argmax = symb;
                count_max = this->symbol2data[symb];
            }
            count_sum += this->symbol2data[symb];
        }
        
        assert (incluBeg <= count_argmax && count_argmax <= incluEnd || !fprintf(stderr, "The value %d is not between %d and %d", count_argmax, incluBeg, incluEnd));
        return 0;
    };
    
    const int
    fillConsensusCounts(
            AlignmentSymbol & count_argmax, unsigned int & count_max, unsigned int & count_sum,
            const SymbolType symbolType) const {
        if (symbolType == BASE_SYMBOL) {
            return this->_fillConsensusCounts(count_argmax, count_max, count_sum, BASE_A, BASE_NN);
        } else if (symbolType == LINK_SYMBOL) {
            return this->_fillConsensusCounts(count_argmax, count_max, count_sum, LINK_M, LINK_NN);
        } else {
            assert(false);
        }
    };
    
    template<ValueType T_SymbolCountType>
    const AlignmentSymbol
    _updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count,
            const SymbolType symbolType, const AlignmentSymbol ambig_pos, unsigned int incvalue2) {
        AlignmentSymbol argmax_count = END_ALIGNMENT_SYMBOLS; AlignmentSymbol(0);
        unsigned int max_count = 0;
        unsigned int sum_count = 0;
        thatSymbol2Count.fillConsensusCounts(argmax_count, max_count, sum_count, symbolType);
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

    template<ValueType T_SymbolCountType>
    std::array<AlignmentSymbol, 2>
    updateByConsensus(const GenericSymbol2Count<TInteger> & thatSymbol2Count, unsigned int incvalue = 1) {
        AlignmentSymbol baseSymb = this->_updateByConsensus<T_SymbolCountType>(thatSymbol2Count, BASE_SYMBOL, BASE_NN, incvalue);
        AlignmentSymbol linkSymb = this->_updateByConsensus<T_SymbolCountType>(thatSymbol2Count, LINK_SYMBOL, LINK_NN, incvalue);
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
            other.fillConsensusCounts(consalpha, countalpha, totalalpha, symbolType);
            if (countalpha >= (thres.getSymbolCount(consalpha)) && countalpha > 0) {
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
    
    const uint32_t incluBegPosition; // end_pos = incluBegPosition + idx2symbol2data
     
protected:
    std::vector<T> idx2symbol2data;
    std::map<uint32_t, std::map<uint32_t   , uint32_t>> pos2dlen2data;
    std::map<uint32_t, std::map<std::string, uint32_t>> pos2iseq2data;
    const size_t _extern2intern4pos(size_t extern_ref_pos) {
        assert(extern_ref_pos >= incluBegPosition);
        assert(extern_ref_pos <  incluBegPosition + idx2symbol2data.size() || !fprintf(stderr, "%d is not within (%d - %d)", extern_ref_pos, incluBegPosition, incluBegPosition + idx2symbol2data.size()));
        return extern_ref_pos - incluBegPosition;
    };
public:
    const char rchars_rpos_to_rchar(const std::string & refstring, unsigned int extern_ref_pos) {
        auto internpos = _extern2intern4pos(extern_ref_pos);
        assert(internpos < refstring.size());
        return refstring.at(internpos);
    };
    
    const uint32_t tid;
    
    CoveredRegion(uint32_t tid, unsigned int beg, unsigned int end): tid(tid), incluBegPosition(beg)  {
        assert (beg < end || !fprintf(stderr, "assertion %d < %d failed!\n", beg, end));
        this->idx2symbol2data = std::vector<T>(end-beg); // TODO: see if + 1 is needed ehre
    };
    
    T &
    getRefByPos(const unsigned int pos, const bam1_t *bam = NULL) {
        assert(pos >= this->incluBegPosition || !fprintf(stderr, "%d >= %d failed for qname %s !!\n", pos, this->incluBegPosition, (NULL != bam ? bam_get_qname(bam) : "?")));
        unsigned int pos2 = pos - this->incluBegPosition;
        assert(pos2 < idx2symbol2data.size() || !fprintf(stderr, "%d  < %d failed for qname %s !!\n", pos, this->incluBegPosition + idx2symbol2data.size(), (NULL != bam ? bam_get_qname(bam) : "?")));
        return this->idx2symbol2data[pos2];
    };
    
    const T &
    getByPos(const unsigned int pos, const bam1_t *bam = NULL) const {
        assert(pos >= this->incluBegPosition || !fprintf(stderr, "%d >= %d failed for qname %s !\n", pos, this->incluBegPosition, (NULL != bam ? bam_get_qname(bam) : "?")));
        unsigned int pos2 = pos - this->incluBegPosition;
        assert(pos2 < idx2symbol2data.size() || !fprintf(stderr, "%d  < %d failed for qname %s !\n", pos, this->incluBegPosition + idx2symbol2data.size(), (NULL != bam ? bam_get_qname(bam) : "?")));
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
    getPosToDlenToData() const {
        return pos2dlen2data;
    };
    const std::map<uint32_t, std::map<std::string, uint32_t>> & 
    getPosToIseqToData() const {
        return pos2iseq2data;
    };
    
    std::map<uint32_t, std::map<uint32_t   , uint32_t>> & 
    getRefPosToDlenToData() {
        return pos2dlen2data;
    };
    std::map<uint32_t, std::map<std::string, uint32_t>> & 
    getRefPosToIseqToData() {
        return pos2iseq2data;
    };


    template <LinkType TLinkType>
    const auto &
    getPosToIndelToData() const {
        static_assert(INS_LINK == TLinkType || DEL_LINK == TLinkType);
        if (TLinkType == INS_LINK) {
            return pos2iseq2data;
        } else {
            return pos2dlen2data;
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
        assert (tid == UINT32_MAX || aln->core.tid == tid);
        tid = aln->core.tid;
        inc_beg = MIN(inc_beg, aln->core.pos);
        exc_end = MAX(exc_end, bam_endpos(aln)) + 1; // accounts for insertion at the end 
    }
    assert (tid != UINT32_MAX);
    assert (inc_beg < exc_end);
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
};

template <class TSymbol2Count>
class GenericSymbol2CountCoverage : public CoveredRegion<TSymbol2Count> {
public:
    GenericSymbol2CountCoverage() : CoveredRegion<TSymbol2Count>(0, 0, 1) { /* assert(false); */ };
    GenericSymbol2CountCoverage(auto tid, auto beg, auto end) : CoveredRegion<TSymbol2Count>(tid, beg, end) {}
    
    void
    assertUpdateIsLegal(const GenericSymbol2CountCoverage<TSymbol2Count> & other) const {
        assert(this->tid == other.tid);
        assert(this->getIncluBegPosition() <= other.getIncluBegPosition() || !fprintf(stderr, "%d <= %d failed!", this->getIncluBegPosition(), other.getIncluBegPosition()));
        assert(this->getExcluEndPosition() >= other.getExcluEndPosition() || !fprintf(stderr, "%d >= %d failed!", this->getExcluEndPosition(), other.getExcluEndPosition())); 
    }
    
    std::vector<unsigned int> computeZeroBasedPosToInsLenVec(unsigned int & totInsLen) {
        std::vector<unsigned int> ret(this->getExcluEndPosition() - this->getIncluBegPosition(), 0);
        for (auto & refPosToIseqToData : this->getRefPosToIseqToData()) {
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
                AlignmentSymbol consymbol = this->getRefByPos(epos).updateByRepresentative<TIsIncVariable>(other.getByPos(epos), incvalue);
                if (update_pos2indel2count) {
                    if (isSymbolIns(consymbol)) {
                        posToIndelToCount_updateByRepresentative<TIsIncVariable>(this->pos2iseq2data, other.getPosToIseqToData(), epos, incvalue);
                    } else if (isSymbolDel(consymbol)) {
                        posToIndelToCount_updateByRepresentative<TIsIncVariable>(this->pos2dlen2data, other.getPosToDlenToData(), epos, incvalue);
                    }
                }
            }
        }
    }

    // mainly for merging R1 and R2 into one read
    template<ValueType T_ConsensusType> void
    updateByConsensus(const GenericSymbol2CountCoverage<TSymbol2Count> &other, unsigned int incvalue = 1,
            const bool update_pos2indel2count = true, const bool update_idx2symbol2data = true) {
        this->assertUpdateIsLegal(other);
        
        if (update_idx2symbol2data) {
            for (size_t epos = other.getIncluBegPosition(); epos < other.getExcluEndPosition(); epos++) {
                const std::array<AlignmentSymbol, NUM_SYMBOL_TYPES> consymbols = this->getRefByPos(epos).updateByConsensus<T_ConsensusType>(other.getByPos(epos), incvalue);
                if (update_pos2indel2count) {
                    if (isSymbolIns(consymbols[1])) {
                        posToIndelToCount_updateByConsensus(this->pos2iseq2data, other.getPosToIseqToData(), epos, incvalue);
                    } else if (isSymbolDel(consymbols[1])) {
                        posToIndelToCount_updateByConsensus(this->pos2dlen2data, other.getPosToDlenToData(), epos, incvalue);
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
                    posToIndelToCount_updateByConsensus(this->pos2iseq2data, other.getPosToIseqToData(), epos, incvalue);
                } else if (isSymbolDel(consymbols[1])) {
                    posToIndelToCount_updateByConsensus(this->pos2dlen2data, other.getPosToDlenToData(), epos, incvalue);
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
    incIns(const unsigned int epos, const std::string & iseq, const uint32_t incvalue = 1) {
        assert (incvalue > 0); 
        assert (iseq.size() > 0);
        size_t ipos = epos; //  _extern2intern4pos(epos);
        posToIndelToCount_inc(this->pos2iseq2data, ipos, iseq, incvalue);
    };

    void // GenericSymbol2CountCoverage<TSymbol2Count>::
    incDel(const unsigned int epos, const uint32_t dlen, const uint32_t incvalue = 1) {
        assert (incvalue > 0);
        assert (dlen > 0);
        size_t ipos = epos; // _extern2intern4pos(epos);
        posToIndelToCount_inc(this->pos2dlen2data, ipos, dlen, incvalue);
    };
    
    template<ValueType TUpdateType>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByAln(const bam1_t *const b, const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred, uint32_t primerlen = 0) {
        static_assert(BASE_QUALITY_MAX == TUpdateType || SYMBOL_COUNT_SUM == TUpdateType);
        assert(this->tid == b->core.tid);
        assert(this->getIncluBegPosition() <= b->core.pos   || !fprintf(stderr, "%d <= %d failed", this->getIncluBegPosition(), b->core.pos));
        assert(this->getExcluEndPosition() >= bam_endpos(b) || !fprintf(stderr, "%d >= %d failed", this->getExcluEndPosition(), bam_endpos(b)));
        
        unsigned int qpos = 0;
        unsigned int rpos = b->core.pos;
        const uint32_t n_cigar = b->core.n_cigar;
        const uint32_t *cigar =  bam_get_cigar(b);
        const uint8_t *bseq = bam_get_seq(b);
        unsigned int incvalue = 1;
        for (unsigned int i = 0; i < n_cigar; i++) {
            uint32_t c = cigar[i];
            unsigned int cigar_op = bam_cigar_op(c);
            unsigned int cigar_oplen = bam_cigar_oplen(c);
            if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
                for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                    assert(rpos >= b->core.pos && rpos < bam_endpos(b) 
                            || !fprintf(stderr, "Bam line with QNAME %s has rpos that is not within the range (%d - %d)", bam_get_qname(b), b->core.pos, bam_endpos(b)));
                    if (i2 > 0) {
                        if (TUpdateType == BASE_QUALITY_MAX) {
                            incvalue = MIN(bam_phredi(b, qpos-1), bam_phredi(b, qpos)) + symbolType2addPhred[LINK_SYMBOL];
                        }
                        this->inc<TUpdateType>(rpos, LINK_M, incvalue, b);
                    }
                    unsigned int base4bit = bam_seqi(bseq, qpos);
                    unsigned int base3bit = seq_nt16_int[base4bit];
                    if (TUpdateType == BASE_QUALITY_MAX) {
                        incvalue = bam_phredi(b, qpos) + symbolType2addPhred[BASE_SYMBOL];
                    }
                    this->inc<TUpdateType>(rpos, AlignmentSymbol(base3bit), incvalue, b);
                    rpos += 1;
                    qpos += 1;
                }
            } else if (cigar_op == BAM_CINS) {
                if  (TUpdateType == BASE_QUALITY_MAX) {
                    if (0 == qpos || qpos + cigar_oplen >= b->core.l_qseq) {
                        LOG(logWARNING) << "Query " << bam_get_qname(b) << " has insertion of legnth " << cigar_oplen << " at " << qpos
                                << " which is not exclusively between 0 and " << b->core.l_qseq << " aligned to tid " << b->core.tid << " and position " << rpos;
                        incvalue = (0 != qpos ? bam_phredi(b, qpos-1) : ((qpos + cigar_oplen < b->core.l_qseq) ? bam_phredi(b, qpos + cigar_oplen) : 1)) 
                                + symbolType2addPhred[LINK_SYMBOL];
                    } else {
                        incvalue = MIN(bam_phredi(b, qpos-1), bam_phredi(b, qpos + cigar_oplen)) + symbolType2addPhred[LINK_SYMBOL];
                    }
                }
                this->inc<TUpdateType>(rpos, insLenToSymbol(cigar_oplen), incvalue, b);
                std::string iseq;
                iseq.reserve(cigar_oplen);
                unsigned int incvalue2 = 128;
                for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                    unsigned int base4bit = bam_seqi(bseq, qpos+i2);
                    const char base8bit = seq_nt16_str[base4bit];
                    iseq.push_back(base8bit);
                    if (TUpdateType == BASE_QUALITY_MAX) {
                        incvalue2 = MIN(incvalue2, bam_seqi(bseq, qpos+i2));
                    }
                }
                if (TUpdateType == BASE_QUALITY_MAX) {
                    incvalue = incvalue2 + symbolType2addPhred[LINK_SYMBOL];
                }
                this->incIns(rpos, iseq, incvalue);
                qpos += cigar_oplen;
            } else if (cigar_op == BAM_CDEL) {
                if (TUpdateType == BASE_QUALITY_MAX) {
                    incvalue = MIN(bam_phredi(b, qpos), bam_phredi(b, qpos+1)) + symbolType2addPhred[LINK_SYMBOL];
                }
                this->inc<TUpdateType>(rpos, delLenToSymbol(cigar_oplen), incvalue, b);
                this->incDel(rpos, cigar_oplen, incvalue);
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
    }

    template<ValueType TUpdateType>
    int // GenericSymbol2CountCoverage<TSymbol2Count>::
    updateByRead1Aln(std::vector<bam1_t *> aln_vec, const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred) {
        for (bam1_t *aln : aln_vec) {
            this->updateByAln<TUpdateType>(aln, symbolType2addPhred);
        }
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
    
    //std::array<Symbol2CountCoverage, 2> bq_imba_depth;
    std::array<Symbol2CountCoverage, 2> bq_qual_phsum;
    std::array<Symbol2CountCoverage, 2> du_bias_dedup;  
    
    std::array<Symbol2CountCoverage, 2> bq_amax_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> bq_amax_rdist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> bq_bias_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> bq_bias_rdist; // edge distance to end
    
    std::array<Symbol2CountCoverage, 2> bq_bsum_ldist;
    std::array<Symbol2CountCoverage, 2> bq_bsum_rdist;
    std::array<Symbol2CountCoverage, 2> bq_bias_1stra;
    std::array<Symbol2CountCoverage, 2> bq_bias_2stra;
 
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

    std::array<Symbol2CountCoverage, 2> fq_qual_phsum;
    std::array<Symbol2CountCoverage, 2> fq_amax_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> fq_amax_rdist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> fq_bias_ldist; // edge distance to end
    std::array<Symbol2CountCoverage, 2> fq_bias_rdist; // edge distance to end
    
    std::array<Symbol2CountCoverage, 2> fq_bsum_ldist;
    std::array<Symbol2CountCoverage, 2> fq_bsum_rdist;
    std::array<Symbol2CountCoverage, 2> fq_bias_1stra;
    std::array<Symbol2CountCoverage, 2> fq_bias_2stra;

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
    
    Symbol2CountCoverageString additional_note;
    
    Symbol2CountCoverageSet(unsigned int t, unsigned int beg, unsigned int end):
        tid(t), incluBegPosition(beg), excluEndPosition(end)
        //, bq_imba_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_qual_phsum({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , du_bias_dedup({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , bq_amax_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_amax_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , bq_bsum_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bsum_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_1stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_bias_2stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})

        , bq_tsum_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_pass_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_pass_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_badep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , bq_vars_vqual({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
       
        // , fq_imba_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_qual_phsum({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_amax_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_amax_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , fq_bsum_ldist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bsum_rdist({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_1stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_bias_2stra({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , fq_tsum_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_pass_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_pass_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_depth({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_badep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_thres({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fq_vars_vqual({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        
        , major_amplicon({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , minor_amplicon({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fam_total_dep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fam_size1_dep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})
        , fam_nocon_dep({Symbol2CountCoverage(t, beg, end), Symbol2CountCoverage(t, beg, end)})

        , duplex_pass_depth(Symbol2CountCoverage(t, beg, end))
        , duplex_tsum_depth(Symbol2CountCoverage(t, beg, end))
        , dedup_ampDistr({Symbol2Bucket2CountCoverage(t, beg, end), Symbol2Bucket2CountCoverage(t, beg, end)})
        , pb_dist_lpart({Symbol2Bucket2CountCoverageEdgeDist(t, beg, end), Symbol2Bucket2CountCoverageEdgeDist(t, beg, end)})
        , pb_dist_rpart({Symbol2Bucket2CountCoverageEdgeDist(t, beg, end), Symbol2Bucket2CountCoverageEdgeDist(t, beg, end)})
        
        , additional_note(Symbol2CountCoverageString(t, beg, end))
    {
        assert(beg < end);
        assert(this->bq_tsum_depth[0].getIncluBegPosition() == beg);
        assert(this->bq_tsum_depth[1].getIncluBegPosition() == beg);
        Symbol2CountCoverage cov(t, beg, end);
    };
   
    template <bool TIsFilterStrong=false>
    int getbest( // auto & qual_psum, 
            auto & max_pqual, auto & best_phred, auto & best_count,
            const auto & ampDistrByPos, const double symbolTypeSum, const AlignmentSymbol symbol, const unsigned int bias_adjusted_mincount, 
            const unsigned int phred_max, const double homogeneity = 0) const {
        max_pqual = 0;
        best_phred = 0;
        best_count = 0;
        auto tot_count = 0;
        for (unsigned int rev_buc_idx = 0; rev_buc_idx < NUM_BUCKETS; rev_buc_idx++) {
            auto bucket = NUM_BUCKETS - 1 - rev_buc_idx;
            auto count = ampDistrByPos.getSymbolBucketCount(symbol, bucket);
            tot_count += count;
            auto phred = MIN(bucket2phred(bucket), phred_max);
            auto tot_pqual = 0;
            assert(tot_count <= symbolTypeSum || !fprintf(stderr, "%d <= %f failed for symbol %d and bucket %d !!!\n", tot_count, symbolTypeSum, symbol, bucket));
            if (0 < count) {
                if (TIsFilterStrong) {
                    if (tot_count - count <= (bias_adjusted_mincount)) {
                        tot_pqual = h01_to_phredlike<false>(phred2prob(phred), 1 + DBL_EPSILON, MIN(tot_count, bias_adjusted_mincount), symbolTypeSum);
                    }
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
    };
    
    std::pair<unsigned int, unsigned int>
    adabias(const auto & t0v, const auto & t1v) {
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
        unsigned int prev_biasfact100 = 0;
        unsigned int pre2_biasfact100 = 0;
        for (size_t i = 0; i < usize - 1; i++) {
            cur0 += (double)t0v[i];
            cur1 += (double)t1v[i];
            unsigned int curr_biasfact100 = any4_to_biasfact100(sum0 - cur0, cur0, sum1 - cur1, cur1);
            // LOG(logINFO) << "At " << i << " : " << cur0 << " , " << sum0 - cur0 << " , " << cur1 << " , " << sum1 - cur1;
            unsigned int norm_biasfact100 = MIN(curr_biasfact100, MIN(prev_biasfact100, pre2_biasfact100));
            if (norm_biasfact100 > max_biasfact100) {
                max_biasfact100 = norm_biasfact100;
                argmax = i+1;
            }
            pre2_biasfact100 = prev_biasfact100;
            prev_biasfact100 = curr_biasfact100;
        }
        return std::make_pair(edbuck2pos(argmax), max_biasfact100);
    }
    
    template<bool TUsePrev = false>
    double 
    adafilter(
            auto & du_bias_dedup,
            const auto & bq_qual_phsum, const auto & bq_tsum_depth,
            auto & amax_ldist, auto & amax_rdist, auto & bias_ldist, auto & bias_rdist, 
            const auto & bsum_ldist, const auto & bsum_rdist, auto & bias_1stra, auto & bias_2stra,
            const auto & curr_tsum_depth, 
            auto & pass_thres, auto & pass_depth,
            auto & vars_thres, auto & vars_depth, auto & vars_badep, auto & vars_vqual,
            auto & dedup_ampDistr, const auto & prev_tsum_depth,
            auto & pb_dist_lpart, auto & pb_dist_rpart, auto & additional_note, bool should_add_note, unsigned int phred_max) {
        
        assert(dedup_ampDistr.at(0).getIncluBegPosition() == dedup_ampDistr.at(1).getIncluBegPosition());
        assert(dedup_ampDistr.at(0).getExcluEndPosition() == dedup_ampDistr.at(1).getExcluEndPosition());
        for (unsigned int strand = 0; strand < 2; strand++) {
            for (auto pos = dedup_ampDistr.at(strand).getIncluBegPosition(); pos < dedup_ampDistr.at(strand).getExcluEndPosition(); pos++) {
                for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+(unsigned int)symbolType)) {
                    // prepare duplication bias
                    const auto prev_depth_typesum = (TUsePrev ? prev_tsum_depth[strand].getByPos(pos).sumBySymbolType(symbolType) : 0);
                    const auto curr_depth_typesum = curr_tsum_depth[strand].getByPos(pos).sumBySymbolType(symbolType); 
                    
                    // prepare positional bias
                    Bucket2CountEdgeDist vsum_pb_dist_lpart = pb_dist_lpart[strand].getByPos(pos).vectorsumBySymbolType(symbolType);
                    Bucket2CountEdgeDist vsum_pb_dist_rpart = pb_dist_rpart[strand].getByPos(pos).vectorsumBySymbolType(symbolType);
                    
                    // prepare strand bias
                    auto typebsum_uqual_v0 = bq_qual_phsum[1-strand].getByPos(pos).sumBySymbolType(symbolType);
                    auto typetsum_depth_v0 = bq_tsum_depth[1-strand].getByPos(pos).sumBySymbolType(symbolType);
                    double typesum_uqual_v0_avg = typebsum_uqual_v0/ (double)(typetsum_depth_v0 + DBL_MIN);
                    
                    auto typebsum_ldist_v0 = bsum_ldist[1-strand].getByPos(pos).sumBySymbolType(symbolType);
                    auto typebsum_rdist_v0 = bsum_rdist[1-strand].getByPos(pos).sumBySymbolType(symbolType); 
                    
                    const auto dp0 = curr_tsum_depth.at(1-strand).getByPos(pos).sumBySymbolType(symbolType);
                    const auto dp1 = curr_tsum_depth.at(0+strand).getByPos(pos).sumBySymbolType(symbolType);
                    
                    for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType];
                            symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType];
                            symbol = AlignmentSymbol(1+((unsigned int)symbol))) {
                        auto curr_depth_symbsum = curr_tsum_depth[strand].getByPos(pos).getSymbolCount(symbol);
                        unsigned int max_imba_depth = curr_depth_symbsum + 1; // magic number meaning no limit on imba depth
if (curr_depth_symbsum * 4 <= curr_depth_typesum * 3 && curr_depth_symbsum > 0 && SYMBOL_TYPE_TO_AMBIG[symbolType] != symbol) {
                        // compute duplication bias
                        double dup_imba = 1;
                        if (TUsePrev) {
                            auto prev_depth_symbsum = prev_tsum_depth[strand].getByPos(pos).getSymbolCount(symbol);
                            auto db100 = any4_to_biasfact100(curr_depth_symbsum, curr_depth_typesum, prev_depth_symbsum, prev_depth_typesum, true);
                            du_bias_dedup[strand].getRefByPos(pos).incSymbolCount(symbol, db100);
                            dup_imba = biasfact100_to_imba(db100);
                        }
                        
                        // compute positional bias
                        auto pb_ldist_pair = adabias(vsum_pb_dist_lpart, pb_dist_lpart[strand].getByPos(pos).getSymbolCounts(symbol));
                        amax_ldist[strand].getRefByPos(pos).incSymbolCount(symbol, pb_ldist_pair.first);
                        bias_ldist[strand].getRefByPos(pos).incSymbolCount(symbol, pb_ldist_pair.second);
                        
                        auto pb_rdist_pair = adabias(vsum_pb_dist_rpart, pb_dist_rpart[strand].getByPos(pos).getSymbolCounts(symbol));
                        amax_rdist[strand].getRefByPos(pos).incSymbolCount(symbol, pb_rdist_pair.first);
                        bias_rdist[strand].getRefByPos(pos).incSymbolCount(symbol, pb_rdist_pair.second);
                        
                        auto pb_ldist_imba = biasfact100_to_imba(bias_ldist[strand].getRefByPos(pos).getSymbolCount(symbol));
                        auto pb_rdist_imba = biasfact100_to_imba(bias_rdist[strand].getRefByPos(pos).getSymbolCount(symbol));
                        
                        if (should_add_note) {
                            unsigned int allcurr, altcurr, allrest, altrest;
                            allcurr = altcurr = allrest = altrest = 0;
                            for (int i = 0; i < vsum_pb_dist_lpart.size(); i++) {
                                allrest += vsum_pb_dist_lpart[i];
                                altrest += pb_dist_lpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                            }
                            this->additional_note.getRefByPos(pos).at(symbol) += "//(";
                            for (int i = 0; i < vsum_pb_dist_lpart.size(); i++) {
                                allcurr += vsum_pb_dist_lpart[i];
                                altcurr += pb_dist_lpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                                this->additional_note.getRefByPos(pos).at(symbol) += std::to_string(i) + "(" + std::to_string(edbuck2pos(i)) + "/" 
                                        + std::to_string(allrest-allcurr) + "/" + std::to_string(allcurr) + "/" + std::to_string(altrest-altcurr) + "/" + std::to_string(altcurr)  + ")";
                                        //+ std::to_string(altcurr) + "/" + std::to_string(altrest-altcurr) + "/" + std::to_string(allcurr) + "/" + std::to_string(allrest-allcurr)  + ")";
                            }

                            allcurr = altcurr = allrest = altrest = 0;
                            for (int i = 0; i < vsum_pb_dist_rpart.size(); i++) {
                                allrest += vsum_pb_dist_rpart[i];
                                altrest += pb_dist_rpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                            }
                            for (int i = 0; i < vsum_pb_dist_rpart.size(); i++) {
                                allcurr += vsum_pb_dist_rpart[i];
                                altcurr += pb_dist_rpart[strand].getByPos(pos).getSymbolCounts(symbol)[i];
                                this->additional_note.getRefByPos(pos).at(symbol) += std::to_string(i) + "(" + std::to_string(edbuck2pos(i)) + "/" 
                                        + std::to_string(allrest-allcurr) + "/" + std::to_string(allcurr) + "/" + std::to_string(altrest-altcurr) + "/" + std::to_string(altcurr)  + ")";
                                        // + std::to_string(altcurr) + "/" + std::to_string(altrest-altcurr) + "/" + std::to_string(allcurr) + "/" + std::to_string(allrest-allcurr)  + ")";
                            }
                            this->additional_note.getRefByPos(pos).at(symbol) += ")//";
                        }

                        // compute strand bias
                        auto symbval_uqual_v1 = bq_qual_phsum[0+strand].getByPos(pos).getSymbolCount(symbol);
                        auto symbval_count_v1 = bq_tsum_depth[0+strand].getByPos(pos).getSymbolCount(symbol);
                        double symbval_uqual_v1_avg = symbval_uqual_v1 / (double)(symbval_count_v1 + DBL_MIN);
                        auto uqual_avg_imba = pow((double)10, (symbval_uqual_v1_avg - typesum_uqual_v0_avg) / (double)10);
                        
                        auto ad0 = curr_tsum_depth.at(1-strand).getByPos(pos).getSymbolCount(symbol);
                        auto ad1 = curr_tsum_depth.at(0+strand).getByPos(pos).getSymbolCount(symbol); 
                        
                        auto bsum_ldist_v1 = bsum_ldist[0+strand].getByPos(pos).getSymbolCount(symbol);
                        auto bsum_rdist_v1 = bsum_rdist[0+strand].getByPos(pos).getSymbolCount(symbol);
                        auto bsum_dist_imba0 = ((bsum_ldist_v1 + 1) / (double)(ad1 + 1)) / ((typebsum_rdist_v0 + 1) / (double)(dp0 + 1));
                        auto bsum_dist_imba1 = ((bsum_rdist_v1 + 1) / (double)(ad1 + 1)) / ((typebsum_ldist_v0 + 1) / (double)(dp0 + 1));
                        
                        auto sb100raw = any4_to_biasfact100((double)dp0, (double)dp1, (double)ad0, (double)ad1, false , 2);
                        bias_1stra[strand].getRefByPos(pos).incSymbolCount(symbol, sb100raw);
                        assert(bsum_dist_imba0 + bsum_dist_imba1 > 0 || !fprintf(stderr, "%g + %g > 0 failed! (will encounter division by zero)\n", bsum_dist_imba0, bsum_dist_imba1));
                        
                        auto sb100fin = (unsigned int)(sb100raw / MAX(uqual_avg_imba, MAX(bsum_dist_imba0, bsum_dist_imba1)));
                        bias_2stra[strand].getRefByPos(pos).incSymbolCount(symbol, sb100fin);
                        auto str_imba = biasfact100_to_imba(sb100fin);
                        
                        max_imba_depth = (unsigned int)(curr_depth_symbsum / MAX(dup_imba, MAX(MAX(pb_ldist_imba, pb_rdist_imba), str_imba)) + 0.5);
                        if (should_add_note) {
                            this->additional_note.getRefByPos(pos).at(symbol) += "//(" +
                                    std::to_string(uqual_avg_imba) + "/" + 
                                    std::to_string(bsum_dist_imba0) + "/" + 
                                    std::to_string(bsum_dist_imba1) + "/" + 
                                    ")//";
                        }
}
                        // find best cutoff from families
                        double max_pqual = 0;
                        unsigned int best_phred = 0;
                        unsigned int best_count = 0;
                        if (curr_depth_symbsum > 0) {
                            getbest<false>( //qual_phsum_val, 
                                    max_pqual, best_phred, best_count,
                                    dedup_ampDistr[strand].getRefByPos(pos), curr_depth_typesum, symbol, max_imba_depth, phred_max);
                        } else {
                            best_phred = 0;
                            best_count = 0;
                        }
                        pass_thres[strand].getRefByPos(pos).incSymbolCount(symbol, best_phred);
                        pass_depth[strand].getRefByPos(pos).incSymbolCount(symbol, best_count);
                        
                        if (curr_depth_symbsum > 0) {
                            getbest<true> ( //qual_phsum_val,
                                    max_pqual, best_phred, best_count,
                                    dedup_ampDistr[strand].getRefByPos(pos), curr_depth_typesum, symbol, max_imba_depth, phred_max);
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
            }
        }
    };

    int updateByAlns3UsingBQ(
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map,
            const std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> & alns3, 
            const std::basic_string<AlignmentSymbol> & region_symbolvec,
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred,
            bool should_add_note,
            unsigned int phred_max) {
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                for (const auto & alns1 : alns2) {
                    uint32_t tid2, beg2, end2;
                    fillTidBegEndFromAlns1(tid2, beg2, end2, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid, beg2, end2);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln<BASE_QUALITY_MAX>(alns1, symbolType2addPhred);
                    std::basic_string<std::pair<unsigned int, AlignmentSymbol>> pos_symbol_string;
                    unsigned int ldist_inc = 0;
                    unsigned int rdist_inc = 0;
                    std::vector<unsigned int> posToInsertLen = read_ampBQerr_fragWithR1R2.computeZeroBasedPosToInsLenVec(rdist_inc); // extend
                    for (auto epos = read_ampBQerr_fragWithR1R2.getIncluBegPosition(); epos < read_ampBQerr_fragWithR1R2.getExcluEndPosition(); epos++) {
                        unsigned int ldist = 1 + epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition();
                        unsigned int rdist = read_ampBQerr_fragWithR1R2.getExcluEndPosition() - epos;
                        for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+((unsigned int)symbolType))) {
                            AlignmentSymbol con_symbol;
                            unsigned int con_count, tot_count;
                            read_ampBQerr_fragWithR1R2.getByPos(epos).fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                            assert (con_count * 2 >= tot_count);
                            if (0 == tot_count) { continue; }
                            unsigned int phredlike = (con_count * 2 - tot_count); 
                            
                            this->bq_tsum_depth [strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                            unsigned int edge_baq = MIN(ldist, rdist) * 4;
                            unsigned int overallq = MIN(edge_baq, phredlike);
                            this->bq_qual_phsum [strand].getRefByPos(epos).incSymbolCount(con_symbol, overallq);
                            unsigned int pbucket = phred2bucket(overallq);
                            assert (pbucket < NUM_BUCKETS || !fprintf(stderr, "%d < %d failed at position %d and con_symbol %d symboltype %d plusbucket %d\n", 
                                    pbucket,  NUM_BUCKETS, epos, con_symbol, symbolType, symbolType2addPhred[symbolType]));
                            if (isSymbolIns(con_symbol)) {
                                posToIndelToCount_updateByConsensus(this->bq_tsum_depth[strand].getRefPosToIseqToData(), read_ampBQerr_fragWithR1R2.getPosToIseqToData(), epos, 1);
                            }
                            if (isSymbolDel(con_symbol)) {
                                posToIndelToCount_updateByConsensus(this->bq_tsum_depth[strand].getRefPosToDlenToData(), read_ampBQerr_fragWithR1R2.getPosToDlenToData(), epos, 1);
                            }
                            AlignmentSymbol refsymbol = region_symbolvec[epos-this->dedup_ampDistr.at(strand).getIncluBegPosition()]; 
                            if (areSymbolsMutated(refsymbol, con_symbol)) {
                                pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
                            }
                            this->dedup_ampDistr[strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pbucket, 1);
                            ldist_inc += posToInsertLen[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()];
                            this->pb_dist_lpart [strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pos2edbuck(ldist + ldist_inc), 1);
                            this->pb_dist_rpart [strand].getRefByPos(epos).incSymbolBucketCount(con_symbol, pos2edbuck(rdist + rdist_inc), 1);
                            this->bq_bsum_ldist [strand].getRefByPos(epos).incSymbolCount(con_symbol, ldist + ldist_inc);
                            this->bq_bsum_rdist [strand].getRefByPos(epos).incSymbolCount(con_symbol, rdist + rdist_inc);
                            rdist_inc -= posToInsertLen[epos - read_ampBQerr_fragWithR1R2.getIncluBegPosition()];
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
                this->bq_qual_phsum, this->bq_tsum_depth, //  this->bq_qual_phsum, // this->bq_imba_depth,
                this->bq_amax_ldist, this->bq_amax_rdist, this->bq_bias_ldist, this->bq_bias_rdist,
                this->bq_bsum_ldist, this->bq_bsum_rdist, this->bq_bias_1stra, this->bq_bias_2stra,
                this->bq_tsum_depth, this->bq_pass_thres, this->bq_pass_depth,
                this->bq_vars_thres, this->bq_vars_depth, this->bq_vars_badep, this->bq_vars_vqual,
                this->dedup_ampDistr,this->bq_tsum_depth,
                this->pb_dist_lpart, this->pb_dist_rpart, this->additional_note, should_add_note, phred_max);
        return 0;
    }
    
    int
    updateByAlns3UsingFQ(
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, 
            std::array<unsigned int, 2>> & mutform2count4map,
            const auto & alns3, const std::basic_string<AlignmentSymbol> & region_symbolvec,
            const std::array<unsigned int, NUM_SYMBOL_TYPES> & symbolType2addPhred, 
            bool should_add_note, unsigned int phred_max) {
        unsigned int niters = 0;
        for (const auto & alns2pair2dflag : alns3) {
            const auto & alns2pair = alns2pair2dflag.first;
            niters++;
            bool log_alns2 = niters > 1000 && ispowerof2(niters);
            if (log_alns2) {
                LOG(logINFO) << "Processing the " << niters << "/" << alns3.size() << " th duplex-aware families which has " << alns2pair[0].size() << " and " << alns2pair[1].size() << " members.\n" ;
            }
            assert (alns2pair[0].size() != 0 || alns2pair[1].size() != 0);
            for (unsigned int strand = 0; strand < 2; strand++) {
                const auto & alns2 = alns2pair[strand];
                if (alns2.size() == 0) { continue; }
                uint32_t tid2, beg2, end2;
                fillTidBegEndFromAlns2(tid2, beg2, end2, alns2);
                Symbol2CountCoverage read_family_con_ampl(tid2, beg2, end2); 
                Symbol2CountCoverage read_family_amplicon(tid2, beg2, end2); 
                for (const auto & alns1 : alns2) {
                    uint32_t tid1, beg1, end1;
                    fillTidBegEndFromAlns1(tid1, beg1, end1, alns1);
                    Symbol2CountCoverage read_ampBQerr_fragWithR1R2(tid1, beg1, end1);
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln<BASE_QUALITY_MAX>(alns1, symbolType2addPhred);
                    read_family_con_ampl.updateByConsensus<SYMBOL_COUNT_SUM>(read_ampBQerr_fragWithR1R2);
                    int updateresult = read_family_amplicon.updateByFiltering(read_ampBQerr_fragWithR1R2, this->bq_pass_thres[strand], 1, true, strand);
                    if (log_alns2) {
                        LOG(logDEBUG) << "num-updates = " << updateresult;
                    }
                }
                for (size_t epos = read_family_amplicon.getIncluBegPosition(); epos < read_family_amplicon.getExcluEndPosition(); epos++) {
                    const auto & con_ampl_symbol2count = read_family_amplicon.getByPos(epos);
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
                        // we can use EM algorithm to improve this estimator, but it is probably overkill
                        if (1 >= con_count) { continue; } // 0 means no coverage, 1 means no error correction, 2 means low quality family if the symbols disagree with each other
                        
                        for (AlignmentSymbol 
                                symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; 
                                symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; 
                                symbol = AlignmentSymbol(1+(unsigned int)symbol)) {
                            if (con_symbol != symbol || con_count * 2 <= tot_count) {
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
            bool log_alns2 = niters > 1000 && ispowerof2(niters);
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
                    read_ampBQerr_fragWithR1R2.updateByRead1Aln<BASE_QUALITY_MAX>(aln_vec, symbolType2addPhred);
                    // read_family_amplicon.updateByConsensus<SYMBOL_COUNT_SUM>(read_ampBQerr_fragWithR1R2);
                    int updateresult = read_family_amplicon.updateByFiltering(read_ampBQerr_fragWithR1R2, this->bq_pass_thres[strand], 1, true, strand);
                    if (log_alns2) {
                        LOG(logDEBUG) << "num-updates = " << updateresult << " (second time)";
                    }
                }
                if ((2 == alns2pair2dflag.second) && alns2pair[0].size() > 0 && alns2pair[1].size() > 0) { // is duplex
                    read_duplex_amplicon.updateByConsensus<SYMBOL_COUNT_SUM>(read_family_amplicon);
                }
                std::basic_string<std::pair<unsigned int, AlignmentSymbol>> pos_symbol_string; // this->dedup_ampDistr
                unsigned int ldist_inc = 0;
                unsigned int rdist_inc = 0;
                std::vector<unsigned int> posToInsertLen = read_family_amplicon.computeZeroBasedPosToInsLenVec(rdist_inc); // extend
                for (size_t epos = read_family_amplicon.getIncluBegPosition(); epos < read_family_amplicon.getExcluEndPosition(); epos++) {
                    unsigned int ldist = 1 + epos - read_family_amplicon.getIncluBegPosition();
                    unsigned int rdist = read_family_amplicon.getExcluEndPosition() - epos;
                    for (SymbolType symbolType = SymbolType(0); symbolType < NUM_SYMBOL_TYPES; symbolType = SymbolType(1+(unsigned int)symbolType)) {
                        AlignmentSymbol con_symbol; // = END_ALIGNMENT_SYMBOLS;
                        unsigned int con_count, tot_count;
                        read_family_amplicon.getRefByPos(epos).fillConsensusCounts(con_symbol, con_count, tot_count, symbolType);
                        if (0 == tot_count) { continue; }
                                                
                        unsigned int majorcount = this->major_amplicon[strand].getByPos(epos).getSymbolCount(con_symbol);
                        unsigned int minorcount = this->minor_amplicon[strand].getByPos(epos).getSymbolCount(con_symbol);
                        auto con_bq_pass_thres = this->bq_pass_thres[strand].getByPos(epos).getSymbolCount(con_symbol);
                        double con_bq_pass_prob = phred2prob(con_bq_pass_thres) * (1 - DBL_EPSILON); // < 1
                        assert (con_bq_pass_prob >= pow(10, ((double)-51)/10) 
                                || !fprintf(stderr, "%f >= phred51 failed at position %d and symbol %d!\n", con_bq_pass_prob, epos, con_symbol));
                        double phredlike = h01_to_phredlike<true>(minorcount + 1, majorcount + minorcount + (1.0 / con_bq_pass_prob), con_count, tot_count);
                        // no base quality stuff
                        
                        this->fq_tsum_depth [strand].getRefByPos(epos).incSymbolCount(con_symbol, 1);
                        unsigned int edge_baq = MIN(ldist, rdist) * 4;
                        unsigned int overallq = MIN(edge_baq, phredlike);
                        this->fq_qual_phsum [strand].getRefByPos(epos).incSymbolCount(con_symbol, overallq);
                        unsigned int pbucket = phred2bucket(overallq);
                        assert (pbucket < NUM_BUCKETS || !fprintf(stderr, "%d < %d failed at position %d and con_symbol %d symboltype %d plusbucket %d\n", 
                                pbucket,  NUM_BUCKETS, epos, con_symbol, symbolType, symbolType2addPhred[symbolType]));
                        if (isSymbolIns(con_symbol)) {
                            posToIndelToCount_updateByConsensus(this->fq_tsum_depth[strand].getRefPosToIseqToData(), read_family_amplicon.getPosToIseqToData(), epos, 1);
                        }
                        if (isSymbolDel(con_symbol)) {
                            posToIndelToCount_updateByConsensus(this->fq_tsum_depth[strand].getRefPosToDlenToData(), read_family_amplicon.getPosToDlenToData(), epos, 1);
                        }
                        AlignmentSymbol refsymbol = region_symbolvec[epos - this->dedup_ampDistr.at(strand).getIncluBegPosition()]; 
                        if (areSymbolsMutated(refsymbol, con_symbol)) {
                            pos_symbol_string.push_back(std::make_pair(epos, con_symbol));
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
                this->bq_qual_phsum, this->bq_tsum_depth, //  this->fq_qual_phsum, // this->fq_imba_depth,
                this->fq_amax_ldist, this->fq_amax_rdist, this->fq_bias_ldist, this->fq_bias_rdist, 
                this->fq_bsum_ldist, this->fq_bsum_rdist, this->fq_bias_1stra, this->fq_bias_2stra,
                this->fq_tsum_depth, this->fq_pass_thres, this->fq_pass_depth,
                this->fq_vars_thres, this->fq_vars_depth, this->fq_vars_badep, this->fq_vars_vqual,
                this->dedup_ampDistr,this->bq_tsum_depth,
                this->pb_dist_lpart, this->pb_dist_rpart, this->additional_note, should_add_note, phred_max);
    };
    
    std::basic_string<AlignmentSymbol> string2symbolseq(const std::string & instring) {
        std::basic_string<AlignmentSymbol> ret;
        ret.reserve(instring.size());
        for (size_t i = 0; i < instring.size(); i++) {
            ret[i] = CHAR_TO_SYMBOL.data[instring[i]];
        }
        return ret;
    };
    
    int 
    updateHapMap(std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map, const auto & tsum_depth, unsigned int max_ploidy = 4) {
        for (auto it = mutform2count4map.begin(); it != mutform2count4map.end();) {
            std::basic_string<std::pair<unsigned int, AlignmentSymbol>> mutform = it->first;
            auto counts = it->second;
            std::array<unsigned int, 2> minAD = {UINT_MAX, UINT_MAX};
            for (unsigned int strand = 0; strand < 2; strand++) {
                for (std::pair<unsigned int, AlignmentSymbol> simplemut : mutform) {
                    unsigned int ad = tsum_depth.at(strand).getByPos(simplemut.first).getSymbolCount(simplemut.second);
                    minAD[strand] = MIN(minAD[strand], ad);
                }
            }
            if (counts[0] * max_ploidy <= minAD[0] && counts[1] * max_ploidy <= minAD[1]) {
                mutform2count4map.erase(it++);
            } else {
                it++;
            }
        }
    };
    
    int updateByRegion3Aln(
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map_bq,
            std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> & mutform2count4map_fq,
            const std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> & alns3, 
            const std::string & refstring,
            unsigned int bq_phred_added_misma, unsigned int bq_phred_added_indel, bool should_add_note, unsigned int phred_max_sscs, unsigned int phred_max_dscs) {
        const std::array<unsigned int, NUM_SYMBOL_TYPES> symbolType2addPhred = {bq_phred_added_misma, bq_phred_added_indel};
        std::basic_string<AlignmentSymbol> ref_symbol_string = string2symbolseq(refstring);
        updateByAlns3UsingBQ(mutform2count4map_bq, alns3, ref_symbol_string, symbolType2addPhred, should_add_note, phred_max_sscs); // base qualities
        updateHapMap(mutform2count4map_bq, this->bq_tsum_depth);
        updateByAlns3UsingFQ(mutform2count4map_fq, alns3, ref_symbol_string, symbolType2addPhred, should_add_note, phred_max_sscs); // family qualities
        updateHapMap(mutform2count4map_fq, this->fq_tsum_depth);
    };
};

std::array<unsigned int, 2>
BcfFormat_init(bcfrec::BcfFormat & fmt, const Symbol2CountCoverageSet & symbolDistrSets12, unsigned int refpos, SymbolType symbolType) {
    for (unsigned int strand = 0; strand < 2; strand++) {
        fmt.aAllBQ[strand] = symbolDistrSets12.bq_qual_phsum.at(strand).getByPos(refpos).sumBySymbolType(symbolType); 
        fmt.bDP1[strand] = symbolDistrSets12.bq_tsum_depth.at(strand).getByPos(refpos).sumBySymbolType(symbolType);
        fmt.cDP1[strand] = symbolDistrSets12.fq_tsum_depth.at(strand).getByPos(refpos).sumBySymbolType(symbolType);
        fmt.cDPTT[strand] = symbolDistrSets12.fam_total_dep.at(strand).getByPos(refpos).sumBySymbolType(symbolType);
    }
    fmt.gapSeq.clear();
    fmt.gapbAD1.clear();
    fmt.gapcAD1.clear();
    fmt.gapNum[0] = 0;
    fmt.gapNum[1] = 0;

    fmt.dDP1 = symbolDistrSets12.duplex_tsum_depth.getByPos(refpos).sumBySymbolType(symbolType);
    return {fmt.bDP1[0] + fmt.bDP1[1], fmt.cDP1[0] + fmt.cDP1[1]};
};

#define INDEL_ID 1
#include "instcode.hpp"
#undef INDEL_ID
#define INDEL_ID 2
#include "instcode.hpp"
void
fillByIndelInfo(bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet, 
        const unsigned int strand, const unsigned int refpos, const AlignmentSymbol symbol, 
        const std::string & refstring) {
    assert(isSymbolIns(symbol) || isSymbolDel(symbol));
    if (isSymbolIns(symbol)) {
        fillByIndelInfo2_1(fmt, symbol2CountCoverageSet, strand, refpos, symbol,
                symbol2CountCoverageSet.bq_tsum_depth.at(strand).getPosToIseqToData(),
                symbol2CountCoverageSet.fq_tsum_depth.at(strand).getPosToIseqToData(),
                refstring);
    } else {
        fillByIndelInfo2_2(fmt, symbol2CountCoverageSet, strand, refpos, symbol,
                symbol2CountCoverageSet.bq_tsum_depth.at(strand).getPosToDlenToData(),
                symbol2CountCoverageSet.fq_tsum_depth.at(strand).getPosToDlenToData(),
                refstring);
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
                unsigned int mutpos;
                if (SYMBOL_TYPE_TO_INCLU_BEG[BASE_SYMBOL] <= pos2symbol4pair.first && pos2symbol4pair.first <= SYMBOL_TYPE_TO_INCLU_END[BASE_SYMBOL]) {
                    mutpos = pos2symbol4pair.first + 1;
                } else {
                    mutpos = pos2symbol4pair.first;
                }
                phase_string += std::string("(") + std::to_string(mutpos) + "&" + SYMBOL_TO_DESC_ARR[symbol] + ")";
            }
            phase_string += std::string("&") + std::to_string(counts[0]) + "&" + std::to_string(counts[1]) + ")";
        }
    }
    return phase_string;
}

int
fillBySymbol(bcfrec::BcfFormat & fmt, const Symbol2CountCoverageSet & symbol2CountCoverageSet12, 
        unsigned int refpos, const AlignmentSymbol symbol, const std::string & refstring, unsigned int refstring_offset, 
        const std::vector<std::pair<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>>> & mutform2count4vec_bq,
        std::set<size_t> indices_bq,
        const std::vector<std::pair<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>>> & mutform2count4vec_fq,
        std::set<size_t> indices_fq,
        unsigned int minABQ, // = 25
        unsigned int phred_max_sscs,
        unsigned int phred_max_dscs) {
    fmt.note = symbol2CountCoverageSet12.additional_note.getByPos(refpos).at(symbol);
    for (unsigned int strand = 0; strand < 2; strand++) {
        
        fmt.aAltBQ[strand] = symbol2CountCoverageSet12.bq_qual_phsum.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.aDB[strand]  = symbol2CountCoverageSet12.du_bias_dedup.at(strand).getByPos(refpos).getSymbolCount(symbol);

        fmt.bPTL[strand] = symbol2CountCoverageSet12.bq_amax_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.bPTR[strand] = symbol2CountCoverageSet12.bq_amax_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.bPBL[strand] = symbol2CountCoverageSet12.bq_bias_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bPBR[strand] = symbol2CountCoverageSet12.bq_bias_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        
        fmt.bSDL[strand] = symbol2CountCoverageSet12.bq_bsum_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bSDR[strand] = symbol2CountCoverageSet12.bq_bsum_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bSB1[strand] = symbol2CountCoverageSet12.bq_bias_1stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bSBR[strand] = symbol2CountCoverageSet12.bq_bias_2stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
        
        fmt.bAD1[strand] = symbol2CountCoverageSet12.bq_tsum_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // total allele depth
        fmt.bAD2[strand] = symbol2CountCoverageSet12.bq_pass_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth
        fmt.bQT2[strand] = symbol2CountCoverageSet12.bq_pass_thres.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  threshold
        fmt.bAD3[strand] = symbol2CountCoverageSet12.bq_vars_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth
        fmt.bADB[strand] = symbol2CountCoverageSet12.bq_vars_badep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.bQT3[strand] = symbol2CountCoverageSet12.bq_vars_thres.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  threshold
        fmt.bVQ3[strand] = symbol2CountCoverageSet12.bq_vars_vqual.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth 
        
        // family quality
        fmt.cPTL[strand] = symbol2CountCoverageSet12.fq_amax_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cPTR[strand] = symbol2CountCoverageSet12.fq_amax_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cPBL[strand] = symbol2CountCoverageSet12.fq_bias_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cPBR[strand] = symbol2CountCoverageSet12.fq_bias_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol); 

        fmt.cSDL[strand] = symbol2CountCoverageSet12.fq_bsum_ldist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cSDR[strand] = symbol2CountCoverageSet12.fq_bsum_rdist.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cSB1[strand] = symbol2CountCoverageSet12.fq_bias_1stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cSBR[strand] = symbol2CountCoverageSet12.fq_bias_2stra.at(strand).getByPos(refpos).getSymbolCount(symbol);
 
        fmt.cAD1[strand] = symbol2CountCoverageSet12.fq_tsum_depth.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cAD2[strand] = symbol2CountCoverageSet12.fq_pass_depth.at(strand).getByPos(refpos).getSymbolCount(symbol);  
        fmt.cQT2[strand] = symbol2CountCoverageSet12.fq_pass_thres.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cAD3[strand] = symbol2CountCoverageSet12.fq_vars_depth.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth
        fmt.cADB[strand] = symbol2CountCoverageSet12.fq_vars_badep.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  threshold
        fmt.cQT3[strand] = symbol2CountCoverageSet12.fq_vars_thres.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        fmt.cVQ3[strand] = symbol2CountCoverageSet12.fq_vars_vqual.at(strand).getByPos(refpos).getSymbolCount(symbol); // pass  allele depth 
        
        fmt.cMajor[strand] = symbol2CountCoverageSet12.major_amplicon.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cMinor[strand] = symbol2CountCoverageSet12.minor_amplicon.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cADTT[strand] = symbol2CountCoverageSet12.fam_total_dep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cADT1[strand] = symbol2CountCoverageSet12.fam_size1_dep.at(strand).getByPos(refpos).getSymbolCount(symbol);
        fmt.cADTN[strand] = symbol2CountCoverageSet12.fam_nocon_dep.at(strand).getByPos(refpos).getSymbolCount(symbol); 
        
        fmt.gapNum[strand] = 0;
        if ((0 < fmt.bAD1[strand]) && (isSymbolIns(symbol) || isSymbolDel(symbol))) {
            fillByIndelInfo(fmt, symbol2CountCoverageSet12, strand, refpos, symbol, refstring);
        }
    }
    
    fmt.bHap = mutform2count4map_to_phase(mutform2count4vec_bq, indices_bq);
    fmt.cHap = mutform2count4map_to_phase(mutform2count4vec_fq, indices_fq);
    
    fmt.dAD1 = symbol2CountCoverageSet12.duplex_tsum_depth.getByPos(refpos).getSymbolCount(symbol);
    fmt.dAD3 = symbol2CountCoverageSet12.duplex_pass_depth.getByPos(refpos).getSymbolCount(symbol);
    
    unsigned int regionpos = refpos - refstring_offset;
    const std::string vcfref = refstring.substr(regionpos, 1);
    const std::string vcfalt = std::string(SYMBOL_TO_DESC_ARR[symbol]);
    bool is_novar = (symbol == LINK_M || (isSymbolSubstitution(symbol) && vcfref == vcfalt));
 
    fmt.DP = fmt.cDPTT[0] + fmt.cDPTT[1];
    auto fmtAD = fmt.cADTT[0] + fmt.cADTT[1];
    fmt.FA = (double)(fmtAD) / (double)(fmt.DP);
    
    fmt.bDP = fmt.bDP1[0] + fmt.bDP1[1];
    auto fmtbAD = fmt.cADTT[0] + fmt.cADTT[1];
    fmt.bFA = (double)(fmtbAD) / (double)(fmt.bDP);
    
    if (fmtAD > 0) {
        assert(fmt.FA >= 0);
        if (fmt.FA > (0.8 - DBL_EPSILON)) {
            fmt.GT = (is_novar ? "0|0" : "1|1");
            fmt.GQ = (unsigned int)calc_phred10_likeratio(0.5,  fmtAD, fmt.DP - fmtAD); // homo, so assume hetero is the alternative
        } else if (fmt.FA < (0.2 + DBL_EPSILON)) {
            fmt.GT = (is_novar ? "1|1" : "0|0");
            fmt.GQ = (unsigned int)calc_phred10_likeratio(0.5,  fmtAD, fmt.DP - fmtAD); // homo, so assume hetero is the alternative
        } else {
            fmt.GT = "0|1";
            fmt.GQ = (unsigned int)calc_phred10_likeratio(0.1,  fmtAD, fmt.DP - fmtAD); // hetero, so assume homo is the alternative
        }
    } else {
        fmt.GT = ".|.";
        fmt.GQ = 0;
    }
    // fmt.GQ = (); // 0;
    fmt.HQ[0] = 0; 
    fmt.HQ[1] = 0;
    
    fmt.VType = SYMBOL_TO_DESC_ARR[symbol];
    double lowestVAQ = prob2phred(1 / (double)(fmt.bAD1[0] + fmt.bAD1[1] + 1)) * ((fmt.bAD1[0] + fmt.bAD1[1]) / (fmt.bDP1[0] + fmt.bDP1[1] + DBL_MIN)) / (double)2;
    double maxVAQs[2] = {0, 0};
    std::array<unsigned int, 2> weightedQT3s = {0, 0};
    for (size_t i = 0; i < 2; i++) {
        auto minAD1 = MIN(fmt.bAD1[i], fmt.cAD1[i]); // prevent symbol with suffix = _N
        auto maxAD1 = MAX(fmt.bAD1[i], fmt.cAD1[i]);
        auto gapAD1 = maxAD1 - minAD1;
        if (fmt.bVQ3[i] < fmt.cVQ3[i]) {
            gapAD1 *= 2;
        }
        double currVAQ = (fmt.bVQ3[i] * minAD1 + fmt.cVQ3[i] * gapAD1) / (double)(minAD1 + gapAD1 + DBL_MIN); // prevent div by zero
        weightedQT3s[i] = (fmt.bQT3[i] * minAD1 + fmt.cQT3[i] * gapAD1) / (double)(minAD1 + gapAD1 + DBL_MIN);
        if (fmt.bQT2[i] >= minABQ) {
            maxVAQs[i] = MAX(maxVAQs[i], currVAQ);
        } else {
            maxVAQs[i] = MAX(maxVAQs[i], MIN(currVAQ, minABQ));
        }
    }
    //double minVAQ = MIN(maxVAQs[0], maxVAQs[1]);
    //double maxVAQ = MAX(maxVAQs[0], maxVAQs[1]);
    
    double weightsum = MIN((double)(weightedQT3s[0] + weightedQT3s[1]), phred_max_dscs);
    double doubleVAQfw = maxVAQs[0] + maxVAQs[1] * MIN(1.0, (weightsum - weightedQT3s[0]) / (weightedQT3s[0] + DBL_EPSILON));
    double doubleVAQrv = maxVAQs[1] + maxVAQs[0] * MIN(1.0, (weightsum - weightedQT3s[1]) / (weightedQT3s[1] + DBL_EPSILON));
    
    double doubleVAQ = MAX(doubleVAQfw, doubleVAQrv);
    // double doubleVAQ = maxVAQ + (minVAQ * (phred_max_dscs - phred_max_sscs) / (double)phred_max_sscs);
    double duplexVAQ = (double)fmt.dAD3 * (double)(phred_max_dscs - phred_max_sscs) - (double)(fmt.dAD1 - fmt.dAD3); // h01_to
    fmt.VAQ = MAX(lowestVAQ, doubleVAQ + duplexVAQ);
    return (int)(fmt.bAD1[0] + fmt.bAD1[1]);
};

#include "version.h"
std::string 
generateVcfHeader(const char *ref_fasta_fname, const char *platform, const unsigned int minABQ, 
        unsigned int argc,      const char *const *argv,  
        unsigned int n_targets, const char *const *target_name, const uint32_t *target_len,
        const char *const sampleName) {
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
        ret += std::string("") + std::string(argv[i]) + "\t";
    }
    ret += "\n";
    ret += std::string("") + "##variantCallerInferredParameters=<" + "platform=" + platform + ",minABQ=" + std::to_string(minABQ) + ">\n";
    ret += std::string("") + "##reference=" + ref_fasta_fname + "\n";
    for (size_t i = 0; i < n_targets; i++) {
        ret += std::string("") + "##contig=<ID=" + target_name[i] + ",length=" + std::to_string(target_len[i]) + ">\n";
    }
    ret += std::string("") + "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location, where POS (start position) is one-based inclusive.\">\n";
    
    for (unsigned int i = 0; i < bcfrec::FILTER_NUM; i++) {
        ret += std::string("") + bcfrec::FILTER_LINES[i] + "\n";
    }
    
    for (unsigned int i = 0; i < bcfrec::FORMAT_NUM; i++) {
        ret += std::string("") + bcfrec::FORMAT_LINES[i] + "\n";
    }
    ret += std::string("") + "##FORMAT=<ID=gbDP,Number=1,Type=Integer,Description=\"Minimum duped   fragment depths in the genomic block for SNV and InDel\">\n";
    ret += std::string("") + "##FORMAT=<ID=gcDP,Number=1,Type=Integer,Description=\"Minimum deduped fragment depths in the genomic block for SNV and InDel\">\n";
    ret += std::string("") + "##FORMAT=<ID=gSTS,Number=2,Type=Integer,Description=\"Variant types for start and end positions, where 0 means SNV and 1 means InDel.\">\n";
    ret += std::string("") + "##FORMAT=<ID=gBEG,Number=1,Type=Integer,Description=\"Begin position of the genomic block (one-based inclusive)\">\n";
    ret += std::string("") + "##FORMAT=<ID=gEND,Number=1,Type=Integer,Description=\"End position of the genomic block (one-based inclusive)\">\n";
    ret += std::string("") + "##phasing=partial\n";
    ret += std::string("") + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sampleName + "\n";
    return ret;
}

void
appendVcfRecord(std::string & out_string, std::string & out_string_pass, const Symbol2CountCoverageSet & symbol2CountCoverageSet, 
        const char *tname, unsigned int refpos, 
        const AlignmentSymbol symbol, const bcfrec::BcfFormat & fmt, 
        const std::string & refstring,
        const unsigned int extended_inclu_beg_pos, 
        const double vcfqual_thres,
        const bool should_output_all
        ) {
    
    assert(refpos >= extended_inclu_beg_pos);
    assert(refpos - extended_inclu_beg_pos < refstring.size());
    
    unsigned int editdist = 1;
    const unsigned int regionpos = refpos - extended_inclu_beg_pos;
    const char *altsymbolname = SYMBOL_TO_DESC_ARR[symbol];
    std::string vcfref;
    std::string vcfalt;
    unsigned int vcfpos;
    
    if (isSymbolIns(symbol) || isSymbolDel(symbol)) {
        vcfpos = refpos; // refpos > 0?
        vcfref = (regionpos > 0 ? refstring.substr(regionpos-1, 1) : "n");
        vcfalt = vcfref;
        std::string indelstring;
        if (fmt.gapNum[0] <= 0 && fmt.gapNum[1] <= 0) {
            std::cerr << "Invalid indel detected : " << tname << ", " << refpos << ", " << SYMBOL_TO_DESC_ARR[symbol] << std::endl;
            std::string msg;
            bcfrec::streamAppendBcfFormat(msg, fmt);
            std::cerr << msg << "\n";
            assert(false);
        } else if (fmt.gapNum[0] <= 0 || fmt.gapNum[1] <= 0) {
            indelstring = fmt.gapSeq[0];
        } else {
            indelstring = (fmt.gapbAD1[0] > fmt.gapbAD1[fmt.gapNum[0]] ? fmt.gapSeq[0] : fmt.gapSeq.at(fmt.gapNum[0]));
        }
        editdist = MAX(1, indelstring.size());
        if (indelstring.size() == 0) {
            vcfalt = altsymbolname;
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
    
    const float vcfqual = fmt.VAQ;
    
    bool is_novar = (symbol == LINK_M || (isSymbolSubstitution(symbol) && vcfref == vcfalt));
    std::string vcffilter;
    if (is_novar) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::noVar]) + ";");
    }
    if (vcfqual < 50) {
        vcffilter += (std::string(bcfrec::FILTER_IDS[bcfrec::q50]) + ";");
    }
    if (0 < vcffilter.size()) {
        vcffilter.pop_back();
    }
    
    if (0 == vcffilter.size()) {
        vcffilter += ("PASS");
    }

    if (!is_novar && vcfqual >= vcfqual_thres) {
        out_string_pass += std::string(tname) + "\t" + std::to_string(vcfpos) + "\t.\t" + vcfref + "\t" + vcfalt + "\t" 
            + std::to_string(vcfqual) + "\t" + vcffilter + "\t.\t" + bcfrec::FORMAT_STR_PER_REC + "\t";
        bcfrec::streamAppendBcfFormat(out_string_pass, fmt);
        out_string_pass += "\n";
    }
    
    if (should_output_all) {
        out_string += std::string(tname) + "\t" + std::to_string(vcfpos) + "\t.\t" + vcfref + "\t" + vcfalt + "\t" 
                + std::to_string(vcfqual) + "\t" + vcffilter + "\t.\t" + bcfrec::FORMAT_STR_PER_REC + "\t";
        bcfrec::streamAppendBcfFormat(out_string, fmt);
        out_string += "\n";
    }
}

#endif
