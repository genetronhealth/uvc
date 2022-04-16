#ifndef main_consensus_hpp_INCLUDED
#define main_consensus_hpp_INCLUDED

#include "main_conversion.hpp"
#include "common.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <vector>


enum ConsensusBlockCigarType {
    CONSENSUS_BLOCK_CINS,
    CONSENSUS_BLOCK_CSOFT_CLIP,
    CONSENSUS_BLOCK_END,
};
#define NUM_CONSENSUS_BLOCK_CIGAR_TYPES 2
static_assert(NUM_CONSENSUS_BLOCK_CIGAR_TYPES == CONSENSUS_BLOCK_END);
#define ALL_CONSENSUS_BLOCK_CIGAR_TYPES (std::array<ConsensusBlockCigarType, NUM_CONSENSUS_BLOCK_CIGAR_TYPES> {{ CONSENSUS_BLOCK_CINS, CONSENSUS_BLOCK_CSOFT_CLIP }})

// typedef std::pair<char, int8_t> FastqBaseAndQual;
typedef std::basic_string<std::pair<char, int8_t>> FastqRecord;
typedef std::array<uvc1_qual_t, BASE_NN> BaseToCount;
typedef std::vector<BaseToCount> ConsensusBlock;

void
reverseAndComplement(FastqRecord & fqrec) {
    std::reverse(fqrec.begin(), fqrec.end());
    for (auto & baseBQ :fqrec) {
        baseBQ.second = STATIC_REV_COMPLEMENT.data[baseBQ.second];
    }
}

const FastqRecord
consensusBlockToSeqQual(const ConsensusBlock & cb1) {
    FastqRecord ret;
    // ConsensusBlock & cb1 ; // = pos2conblock4it.first->second;
    for (size_t inspos = 0; inspos < cb1.size(); inspos++) {
        AlignmentSymbol conbase = BASE_NN;
        uvc1_readnum_t concount = 0;
        uvc1_readnum_t totcount = 0;
        for (const AlignmentSymbol posbase : SYMBOL_TYPE_TO_NON_NN_SYMBOLS[BASE_SYMBOL]) {
            if (cb1[inspos][posbase] > concount) {
                conbase = posbase;
                concount = cb1[inspos][posbase];
            }
            totcount += cb1[inspos][posbase];
        }
        const char *desc = SYMBOL_TO_DESC_ARR[conbase];
        assert (strlen(desc) == 1);
        ret.push_back(std::make_pair(desc[0], non_neg_minus(concount * 2, totcount)));
    }
    return ret;
}

struct ConsensusBlockSet {
    std::map<uvc1_refgpos_t, ConsensusBlock> pos2conblock;
    void
    incByPosSeqQual(uvc1_readpos_t pos, const std::string & seq, const auto & qual) {
        assert(seq.size() == qual.size());
        auto pos2conblock4it = this->pos2conblock.insert(std::make_pair(pos, ConsensusBlock()));
        ConsensusBlock & cb2 = pos2conblock4it.first->second;
        while (cb2.size() < seq.size()) {
            cb2.push_back(std::array<uvc1_qual_t, BASE_NN> {{ 0 }});
        }
        for (size_t inspos = 0; inspos < seq.size(); inspos++) {
            AlignmentSymbol posbase = CHAR_TO_SYMBOL.data[seq[inspos]];
            cb2[inspos][posbase] += qual[inspos];
        }
    };
    
    FastqRecord
    returnSeqQualVec(uvc1_readpos_t pos) {
        FastqRecord ret;
        auto pos2conblock4it = this->pos2conblock.find(pos);
        if (pos2conblock4it != this->pos2conblock.end()) {
            return consensusBlockToSeqQual(pos2conblock4it->second);
        } else {
            return FastqRecord();
        }
    };

    void
    incByConsensus(const ConsensusBlockSet & cbset, uvc1_readnum_t incvalue = 1) {
        for (const auto & pos2conblock4it1 : cbset.pos2conblock) {
            uvc1_refgpos_t pos = pos2conblock4it1.first;
            auto pos2conblock4it = this->pos2conblock.insert(std::make_pair(pos, ConsensusBlock()));
            const ConsensusBlock & cb1 = cbset.pos2conblock.at(pos);
            ConsensusBlock & cb2 = pos2conblock4it.first->second;
            while (cb2.size() < cb1.size()) {
                cb2.push_back(std::array<uvc1_qual_t, BASE_NN> {{ 0 }});
            }
            for (size_t inspos = 0; inspos < cb1.size(); inspos++) {
                AlignmentSymbol conbase = BASE_NN;
                uvc1_readnum_t concount = 0;
                for (const AlignmentSymbol posbase : SYMBOL_TYPE_TO_NON_NN_SYMBOLS[BASE_SYMBOL]) {
                    if (cb1[inspos][posbase] > concount) {
                        conbase = posbase;
                        concount = cb1[inspos][posbase];
                    }
                }
                cb2[inspos][conbase] += incvalue;
            }
        }
    };
    
    void
    incByMajorMinusMinor(const ConsensusBlockSet & cbset) {
        for (const auto & pos2conblock4it1 : cbset.pos2conblock) {
            uvc1_refgpos_t pos = pos2conblock4it1.first;
            auto pos2conblock4it = this->pos2conblock.insert(std::make_pair(pos, ConsensusBlock()));
            const ConsensusBlock & cb1 = cbset.pos2conblock.at(pos);
            ConsensusBlock & cb2 = pos2conblock4it.first->second;
            while (cb2.size() < cb1.size()) {
                cb2.push_back(std::array<uvc1_qual_t, BASE_NN> {{ 0 }});
            }
            for (size_t inspos = 0; inspos < cb1.size(); inspos++) {
                AlignmentSymbol conbase = BASE_NN;
                uvc1_readnum_t concount = 0;
                uvc1_readnum_t totcount = 0;
                for (const AlignmentSymbol posbase : SYMBOL_TYPE_TO_NON_NN_SYMBOLS[BASE_SYMBOL]) {
                    if (cb1[inspos][posbase] > concount) {
                        conbase = posbase;
                        concount = cb1[inspos][posbase];
                    }
                    totcount += cb1[inspos][posbase];
                }
                cb2[inspos][conbase] += non_neg_minus(concount * 2, totcount);
            }
        }
    };
    
    /*
    void
    incBySummation(const ConsensusBlockSet & cbset) {
        for (const auto & pos2conblock4it1 : cbset.pos2conblock) {
            uvc1_refgpos_t pos = pos2conblock4it1.first;
            auto pos2conblock4it = this->pos2conblock.insert(std::make_pair(pos, ConsensusBlock()));
            const ConsensusBlock & cb1 = cbset.pos2conblock.at(pos);
            ConsensusBlock & cb2 = pos2conblock4it.first->second;
            while (cb2.size() < cb1.size()) {
                cb2.push_back(std::array<uvc1_qual_t, BASE_NN> {{ 0 }});
            }
            for (size_t inspos = 0; inspos < cb1.size(); inspos++) {
                for (size_t posbase = 0; posbase < BASE_NN; posbase++) {
                    cb2[inspos][posbase] += cb1[inspos][posbase];
                }
            }
        }
    };
    */
};

#endif
