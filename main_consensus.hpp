#ifndef main_consensus_hpp_INCLUDED
#define main_consensus_hpp_INCLUDED

#include "main_conversion.hpp"
#include "common.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <vector>


enum ConsensusBlockCigarType {
    CONSENSUS_BLOCK_CSOFT_CLIP_FIXED_LEFT_TO_VAR_RIGHT,
    CONSENSUS_BLOCK_CINS,
    CONSENSUS_BLOCK_CSOFT_CLIP_FIXED_RIGHT_TO_VAR_LEFT,
    CONSENSUS_BLOCK_END,
};
#define NUM_CONSENSUS_BLOCK_CIGAR_TYPES 3
static_assert(NUM_CONSENSUS_BLOCK_CIGAR_TYPES == CONSENSUS_BLOCK_END);
#define ALL_CONSENSUS_BLOCK_CIGAR_TYPES (std::array<ConsensusBlockCigarType, NUM_CONSENSUS_BLOCK_CIGAR_TYPES> \
    {{ CONSENSUS_BLOCK_CSOFT_CLIP_FIXED_LEFT_TO_VAR_RIGHT, CONSENSUS_BLOCK_CINS, CONSENSUS_BLOCK_CSOFT_CLIP_FIXED_RIGHT_TO_VAR_LEFT }})

bool 
is_ConsensusBlockCigarType_right2left(ConsensusBlockCigarType cigartype) { 
    return (CONSENSUS_BLOCK_CSOFT_CLIP_FIXED_RIGHT_TO_VAR_LEFT == cigartype); 
}; 

// typedef std::pair<char, int8_t> FastqBaseAndQual;
typedef std::basic_string<std::pair<char, int8_t>> FastqRecord;
typedef std::array<uvc1_qual_t, BASE_NN> BaseToCount;
typedef std::vector<BaseToCount> ConsensusBlock;

void
reverseAndComplement(FastqRecord & fqrec) {
    std::reverse(fqrec.begin(), fqrec.end());
    for (auto & baseBQ :fqrec) {
        baseBQ.first = STATIC_REV_COMPLEMENT.data[(size_t)baseBQ.first];
    }
}

ConsensusBlock
ConsensusBlock_trim(const ConsensusBlock & conblock, uvc1_qual_t percDP_thres = 20, uvc1_readpos_t n_consec_positions_thres = 3) {
    uvc1_qual_t maxDP = 0;
    for (const auto base2cnt : conblock) {
        uvc1_qual_t thisDP = 0;
        for (const AlignmentSymbol posbase : SYMBOL_TYPE_TO_NON_NN_SYMBOLS[BASE_SYMBOL]) {
            thisDP += base2cnt[posbase];
        }
        UPDATE_MAX(maxDP, thisDP);
    }
    ConsensusBlock ret;
    uvc1_refgpos_t prev_pos = 0;
    uvc1_refgpos_t curr_pos = 0;
    uvc1_refgpos_t n_consec_positions = 0;
    for (const auto base2cnt : conblock) {
        curr_pos++;
        uvc1_qual_t thisDP = 0;
        for (const AlignmentSymbol posbase : SYMBOL_TYPE_TO_NON_NN_SYMBOLS[BASE_SYMBOL]) {
            thisDP += base2cnt[posbase];
        }
        if (thisDP * 100 >= maxDP * percDP_thres) {
            if (prev_pos + 1 == curr_pos) { n_consec_positions++ ;}
            else { n_consec_positions = 1; }
            if (n_consec_positions >= n_consec_positions_thres) {
                for (uvc1_refgpos_t i = 1; i < n_consec_positions; i++) {
                    ret.pop_back();
                }
                break; 
            }
        }
        prev_pos = curr_pos;
        ret.push_back(base2cnt);
    }
    return ret;
}

const FastqRecord
consensusBlockToSeqQual(const ConsensusBlock & cb1, bool is_right2left, uvc1_readnum_t n_frag_supports) {
    FastqRecord ret;
    // ConsensusBlock & cb1 ; // = pos2conblock4it.first->second;
    for (size_t inspos1 = 0; inspos1 < cb1.size(); inspos1++) {
        size_t inspos = (is_right2left ? (cb1.size() - inspos1 - 1) : inspos1);
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
        ret.push_back(std::make_pair(desc[0], non_neg_minus(concount * 2, totcount) / n_frag_supports));
    }
    return ret;
}

struct ConsensusBlockSet {
    std::map<uvc1_refgpos_t, ConsensusBlock> pos2conblock;
    void setIsRightToLeft(bool is_r2l) { is_right2left = is_r2l; }
    
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
    returnSeqQualVec(uvc1_readpos_t pos, uvc1_readnum_t n_frag_supports, uvc1_qual_t percDP = 20, uvc1_refgpos_t n_consec_positions = 3) {
        FastqRecord ret;
        auto pos2conblock4it = this->pos2conblock.find(pos);
        if (pos2conblock4it != this->pos2conblock.end()) {
            return consensusBlockToSeqQual(ConsensusBlock_trim(pos2conblock4it->second, percDP), is_right2left, n_frag_supports);
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
private: 
    bool is_right2left = false;
};

#endif
