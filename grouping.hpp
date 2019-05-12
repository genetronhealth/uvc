#ifndef hts_parser_hpp_INCLUDED
#define hts_parser_hpp_INCLUDED

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "htslib/sam.h"

#include "common.h"
#include "logging.hpp"

#define logDEBUGx1 logDEBUG // logINFO

// position of 5' is the starting position, but position of 3' is unreliable without mate info.

bool ispowof2(auto num) {
    return (num & (num-1)) == 0;
}

auto min(auto a, auto b) {
    return a < b ? a : b;
}

auto max(auto a, auto b) {
    return a > b ? a : b;
}

int 
bed_fname_to_contigs(
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> & tid_beg_end_e2e_vec,
        const std::string & bed_fname, const bam_hdr_t *bam_hdr) {
    std::map<std::string, unsigned int> tname_to_tid;
    for (size_t i = 0; i < bam_hdr->n_targets; i++) {
        tname_to_tid[bam_hdr->target_name[i]] = i;
    }
    std::ifstream bedfile(bed_fname);
    while (bedfile.good()) {
        std::string line;
        getline(bedfile, line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream linestream(line);
        std::string tname;
        unsigned int tbeg;
        unsigned int tend;
        linestream >> tname;
        linestream >> tbeg;
        linestream >> tend;
        assert (tbeg < tend || std::cerr 
                << "The bedfile " << bed_fname << " does not have its end after its start at: " << tname << "\t" << tbeg << "\t" << tend);
        assert (tname_to_tid.find(tname) != tname_to_tid.end() || std::cerr 
                << "The reference template name " << tname << " from the bedfile " << bed_fname << " is not in the input sam file");
        bool end2end = false;
        std::string token;
        while (linestream.good()) {
            linestream >> token;
            if (token.find("end-to-end") != std::string::npos) {
                end2end = true;
            }
        }
        tid_beg_end_e2e_vec.push_back(std::make_tuple(tname_to_tid[tname], tbeg, tend, end2end, 100*1000)); // assume 100*1000 reads fall into this region
    }
}

int 
sam_fname_to_contigs(
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> & tid_beg_end_e2e_vec,
        std::vector<std::tuple<std::string, unsigned int>> & tid_to_tname_tlen_tuple_vec,
        const std::string & input_bam_fname, const std::string & bed_fname) {
    tid_beg_end_e2e_vec.clear();
    tid_to_tname_tlen_tuple_vec.clear();
    samFile *sam_infile = sam_open(input_bam_fname.c_str(), "r"); // AlignmentFile(samfname, "rb")
    bam_hdr_t * samheader = sam_hdr_read(sam_infile);
    tid_to_tname_tlen_tuple_vec.reserve(samheader->n_targets);
    
    for (unsigned int tid = 0; tid < samheader->n_targets; tid++) {
        tid_to_tname_tlen_tuple_vec.push_back(std::make_tuple(std::string(samheader->target_name[tid]), samheader->target_len[tid]));
    }
    int ret;
    if (bed_fname != NOT_PROVIDED) {
        ret = bed_fname_to_contigs(tid_beg_end_e2e_vec, bed_fname, samheader);    
    } else {
        ret = 0;
        unsigned int tid = UINT_MAX;
        unsigned int tbeg = -1;
        unsigned int tend = -1;
        uint64_t nreads = 0;
        bam1_t *alnrecord = bam_init1();
        while (sam_read1(sam_infile, samheader, alnrecord) >= 0) {
            ret++;
            if (BAM_FUNMAP & alnrecord->core.flag) {
                continue;
            }
            bool is_uncov = (alnrecord->core.tid != tid || alnrecord->core.pos > tend);
            uint64_t n_overlap_positions = (1 + tend - min(tend, alnrecord->core.pos));
            uint64_t npositions = (tend - min(tbeg, tend));
            bool has_many_positions = (npositions * npositions > n_overlap_positions * (1024UL*1024UL*16UL));
            bool has_many_reads = (nreads * nreads > n_overlap_positions * (1024UL*1024UL*1024UL*1024UL));
            if (is_uncov || has_many_positions || has_many_reads) {
                if (tid != UINT_MAX) {
                    tid_beg_end_e2e_vec.push_back(std::make_tuple(tid, tbeg, tend, false, nreads));
                }
                tid = alnrecord->core.tid;
                if (is_uncov) {
                    tbeg = alnrecord->core.pos;
                    tend = bam_endpos(alnrecord);
                } else {
                    tbeg = tend;
                    tend = max(tbeg, bam_endpos(alnrecord)) + 1;
                }
                nreads = 0;
            }
            tend = max(tend, bam_endpos(alnrecord));
            nreads += 1;
        }
        if (tid != UINT_MAX) {
            tid_beg_end_e2e_vec.push_back(std::make_tuple(tid, tbeg, tend, false, nreads));
        }
        bam_destroy1(alnrecord);
    }
    bam_hdr_destroy(samheader);
    sam_close(sam_infile);
    return ret; 
}

/*
enum {
    POS_BOTH_5PS_ARE_RELIABLE = 0;
    POS_ONE_5P = 1;
    POS_TWO_5P = 2;
    //POS_ONLY1_5P_UNPAIRED_READ = 1;
    //POS_ONLY1_5P_NOT_PROPERLY_PAIRED_READ = 2;
    POS_3P_OR_5P_IS_OUT_OF_RANGE = 4;
    POS_OTHER_5P_DOMINATES_THIS_5P = 8; // the other end of the paired read is very high in frequency compared with this end
    POS_OTHER_5P_DOMINATES_EACH_5P = 16; // the other end of the paired read is very high in frequency overall
} PositionType;
*/

const unsigned int ARRPOS_MARGIN = 500;
const int8_t ARRPOS_OUTER_RANGE = 10;
const int8_t ARRPOS_INNER_RANGE = 3;

/**
 * Schematics of coverage
 *        -
 *       ---   -
 *    - ----- ---
 *   ------------
 *   OOOIIICIIIOO
 * where:
 * O denotes outer region (Non-center positions just around the inner regions)
 * I denotes inner region (positions that the center covers)
 * C denotes center (position at which a maximal number of reads start or end)
 * - denotes the start or end of any read passing the threshold
 */

enum FilterReason {
    NOT_FILTERED,
    NOT_MAPPED,
    NOT_PRIMARY_ALN,
    LOW_MAPQ,
    LOW_ALN_LEN,
    OUT_OF_RANGE,
    NOT_END_TO_END,
    NUM_FILTER_REASONS
};

enum FilterReason
fill_isrc_isr2_beg_end_with_aln(bool & isrc, bool & isr2, uint32_t & tBeg, uint32_t & tEnd, unsigned int &num_seqs,
        const bam1_t *aln, uint32_t fetch_tbeg, uint32_t fetch_tend,
        unsigned int min_mapq, unsigned int min_alnlen, bool end2end) {
    num_seqs = 0;
    if (aln->core.flag & 0x4) {
        return NOT_MAPPED; // unmapped
    }
    if ((aln->core.flag & 0x900) != 0) {
        return NOT_PRIMARY_ALN; // unmapped
    }

    if (aln->core.qual < min_mapq) { 
        return LOW_MAPQ; // mapq too low 
    }
    if ((bam_endpos(aln) - aln->core.pos) < min_alnlen) {
        return LOW_ALN_LEN; // alignment length too short
    }

    isrc = ((aln->core.flag & 0x10) == 0x10);
    isr2 = ((aln->core.flag & 0x80) == 0x80 && (aln->core.flag & 0x1) == 0x1);
    const uint32_t begpos = aln->core.pos;
    const uint32_t endpos = bam_endpos(aln) - 1;
    int ret = 0;
    if (((aln->core.flag & 0x1) == 0) || ((aln->core.flag & 0x2) == 0) || (aln->core.flag & 0x8) ||  aln->core.isize == 0) {
        tBeg = (isrc ? endpos : begpos);
        tEnd = (isrc ? begpos : endpos);
        num_seqs = 1;
    } else {
        auto tBegP1 = min(begpos, aln->core.mpos);
        auto tEndP1 = tBegP1 + abs(aln->core.isize) - 1;
        tBeg = (isrc ? tEndP1 : tBegP1); 
        tEnd = (isrc ? tBegP1 : tEndP1);
        num_seqs = 2;
    }
    auto tOrdBeg = min(tBeg, tEnd);
    auto tOrdEnd = max(tBeg, tEnd);
    if (tOrdBeg + (ARRPOS_MARGIN - ARRPOS_OUTER_RANGE) <= fetch_tbeg || fetch_tend - 1 + (ARRPOS_MARGIN - ARRPOS_OUTER_RANGE) <= tOrdEnd) {
        return OUT_OF_RANGE; // is out of range
    }
    if (end2end && !(tOrdBeg <= fetch_tbeg && tOrdEnd >= fetch_tend)) {
        return NOT_END_TO_END; // did not span the entire genomic region
    }
    return NOT_FILTERED;
}

unsigned int unsigned_diff(unsigned int a, unsigned int b) {
    return (a > b ? a - b : b - a);
}

int 
poscounter_to_pos2pcenter(
              std::vector<unsigned int> & pos_to_center_pos,
        const std::vector<unsigned int> & pos_to_count, 
        const unsigned int center_mult=5) {
    
    for (size_t locov_pos = ARRPOS_INNER_RANGE; locov_pos < pos_to_count.size() - ARRPOS_INNER_RANGE; locov_pos++) {
        unsigned int locov_count = pos_to_count[locov_pos];
        pos_to_center_pos[locov_pos] = locov_pos; // identity mapping by default
        unsigned int max_count = locov_count;
        // check if inner_pos is attracted by outer position
        for (int hicov_pos = locov_pos - ARRPOS_INNER_RANGE; hicov_pos < locov_pos + ARRPOS_INNER_RANGE + 1; hicov_pos++) {
            unsigned int hicov_count = pos_to_count[hicov_pos];
            if (hicov_count > max_count && hicov_count > locov_count * pow(center_mult, unsigned_diff(locov_pos, hicov_pos))) {
                pos_to_center_pos[locov_pos] = hicov_pos;
                max_count = hicov_count;
            }
        }
    }
    return 0;
}

struct _RevComplement {
    char data[128];
    char table16[16];
    _RevComplement() {
        for (int i = 0; i < 128; i++) {
            data[i] = (char)i; 
        }
        data['A'] = 'T';
        data['T'] = 'A'; 
        data['C'] = 'G'; 
        data['G'] = 'C'; 
        data['a'] = 't';
        data['t'] = 'a'; 
        data['c'] = 'g'; 
        data['g'] = 'c'; 
        for (int i = 0; i < 16; i++) {
            table16[i] = (char)i;
        }
        table16[1] = 8/1;
        table16[2] = 8/2;
        table16[4] = 8/4;
        table16[8] = 8/8;
    }
};

const _RevComplement THE_REV_COMPLEMENT;

template <class TContainer>
int revcompln(TContainer & str, size_t len) {
    for (size_t i = 0; i < len / 2; i++) {
        autoswap(str[i], str[len-1-i]);
    }
    for (size_t i = 0; i < len; i++) {
        str[i] = THE_REV_COMPLEMENT.data[str[i]];    
    }
}

template <class TContainer>
int revcompl(TContainer & str) {
    return revcompln(str, str.size());
}

//// one-way converion of data into hash values

// https://en.wikipedia.org/wiki/Universal_hashing#Hashing_strings
template <class T> uint64_t 
strnhash(const T *str, size_t n) {
    uint64_t ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * 31UL + ((uint64_t)str[i]);
    }
    return ret;
}

template <class T> uint64_t 
strnhash_rc(const T *str, size_t n) {
    uint64_t ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * 31UL + THE_REV_COMPLEMENT.data[((uint64_t)str[n-i-1])];
    }
    return ret;
}

template<class T> uint64_t 
strhash(const T *str) {
    return strnhash(str, SIZE_MAX);
}

uint64_t 
hash2hash(uint64_t hash1, uint64_t hash2) {
    // return hash1 ^ hash2; 
    return hash1 * ((1UL<<(31UL)) - 1UL) + hash2;
}

int 
clean_fill_strand_umi_readset(
        std::vector<std::array<std::vector<std::vector<bam1_t *>>, 2>> &umi_strand_readset) {
    for (auto & strand_readset : umi_strand_readset) {
        for (unsigned int strand = 0; strand < 2; strand++) {
            for (auto & read : strand_readset[strand]) {
                for (bam1_t *aln : read) {
                    bam_destroy1(aln);
                }
            }
        }
    }
}

int 
fill_strand_umi_readset_with_strand_to_umi_to_reads(
        std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> &umi_strand_readset,
        // std::array<std::vector<std::vector<std::vector<bam1_t *>>>, 2> &strand_umi_readset,
        std::map<uint64_t, std::pair<std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>, int>> &umi_to_strand_to_reads
        // std::array<std::map<uint64_t, std::map<uint64_t, std::vector<bam1_t *>>>, 2> &strand_to_umi_to_reads
        ) {
#if 1
    for (auto & umi_to_strand_to_reads_element : umi_to_strand_to_reads) {
        const auto umi = umi_to_strand_to_reads_element.first;
        const auto strand_to_reads = umi_to_strand_to_reads_element.second.first;
        const auto dflag = umi_to_strand_to_reads_element.second.second;
        umi_strand_readset.push_back(std::make_pair(std::array<std::vector<std::vector<bam1_t *>>, 2>(), dflag));
        for (unsigned int strand = 0; strand < 2; strand++) {
            // umi_strand_readset.back()[strand] = std::vector<std::vector<bam1_t *>>();
            for (auto read : strand_to_reads[strand]) {
                const std::vector<bam1_t *> alns = read.second;
                umi_strand_readset.back().first[strand].push_back(std::vector<bam1_t *>());
                for (auto aln : alns) {
                    umi_strand_readset.back().first[strand].back().push_back(aln);
                }
            }
        }
    }
#else // this is the old data structure that is not duplex aware
    for (unsigned int strand = 0; strand < 2; strand++) {
        for (auto umi_to_reads : strand_to_umi_to_reads.at(strand)) {
            strand_umi_readset.at(strand).push_back(std::vector<std::vector<bam1_t *>>());
            for (auto reads : umi_to_reads.second) {
                strand_umi_readset.at(strand).back().push_back(std::vector<bam1_t*>());
                for (auto read : reads.second) {
                    strand_umi_readset.at(strand).back().back().push_back(read);
                }
            }
        }
    }
#endif
    return 0;
};

template <bool is_rc>
uint64_t 
bam2umihash(int & is_umi_found, const bam1_t *aln, const std::vector<uint8_t> & UMI_STRUCT, const int max_begin_diff_umi2read = 5) {
    LOG(logDEBUGx1) << "Going over " << UMI_STRUCT.size() << " bases in the pattern";

    auto *bamseq = bam_get_seq(aln);
    
    for (int i = 0; i < max_begin_diff_umi2read; i++) {
        int patpos = 0;
        uint64_t umihash = 0;
        for (int j = i; j < aln->core.l_qseq && patpos < UMI_STRUCT.size(); j++) {
            // int offset = i + j;
            char int4base;
            if (is_rc) {
                char int4base2 = bam_seqi(bamseq, aln->core.l_qseq - 1 - j);
                int4base = THE_REV_COMPLEMENT.table16[int4base2];
            } else {
                int4base = bam_seqi(bamseq, j);
            }
            if (UMI_STRUCT[patpos] == int4base || 15 == UMI_STRUCT[patpos]) {
                if (0xF == UMI_STRUCT[patpos]) {
                    umihash = umihash * 16 + int4base;
                }
                patpos++;
            } else {
                LOG(logDEBUGx1) << "Misma at query position " << j << " (" << (int)int4base << ") and pattern position " << patpos << " (" << (int)UMI_STRUCT[patpos] << ") for read " << bam_get_qname(aln);
                break;
            }
        }
        if (UMI_STRUCT.size() == patpos) {
            is_umi_found++;
            LOG(logDEBUGx1) << "UMI-is-found: " << patpos << " / " << UMI_STRUCT.size() << " with flag " << is_umi_found << " and hash value " << umihash;
            return umihash;
        } else {
            LOG(logDEBUGx1) << "Fraction of bases in UMI that are found: " << patpos << " / " << UMI_STRUCT.size() << " ";
        }
    }
    return 0;
};

int
bamfname_to_strand_to_familyuid_to_reads(
        std::map<uint64_t, std::pair<std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>, int>> &umi_to_strand_to_reads,
        unsigned int & extended_inclu_beg_pos, unsigned int & extended_exclu_end_pos,
        const std::string input_bam_fname, ErrorCorrectionType input_mode, 
        unsigned int tid, unsigned int fetch_tbeg, unsigned int fetch_tend, 
        bool end2end, unsigned int min_mapq, unsigned int min_alnlen, 
        unsigned int contig_ordinal,
        const std::string UMI_STRUCT_STRING, 
        const hts_idx_t * hts_idx
        ) {
    assert (fetch_tend > fetch_tbeg);
    
    std::vector<uint8_t> umi_struct_string16;
    for (auto ch : UMI_STRUCT_STRING) {
        umi_struct_string16.push_back(seq_nt16_table[ch]);
    }
    for (auto base : umi_struct_string16) {
        LOG(logDEBUGx1) << "Base " << (int)base;
    }
    extended_inclu_beg_pos = INT32_MAX;
    extended_exclu_end_pos = 0;
    
    unsigned int alignmentpassed, lenpassed, cigarpassed, umipassed, pcrpassed;
    alignmentpassed = lenpassed = cigarpassed = umipassed = pcrpassed = 0;
   
    samFile *sam_infile = sam_open(input_bam_fname.c_str(), "r"); // AlignmentFile(input_bam_fname, "rb");
    LOG(logINFO) << "Start processing the chunk " << tid << ":" << fetch_tbeg << "-" << fetch_tend << " (contig no " << contig_ordinal << ")";
    
    unsigned int fetch_size = fetch_tend - fetch_tbeg + (ARRPOS_MARGIN + ARRPOS_OUTER_RANGE) * 2; // arrposIsSpecial
    
    std::vector<unsigned int> inicount(fetch_size, 0);
    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_beg_count = { inicount, inicount, inicount, inicount };
    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_end_count = { inicount, inicount, inicount, inicount };
    
    // hts_idx_t * hts_idx = sam_index_load2(sam_infile, input_bam_fname.c_str(), NULL);
    hts_itr_t * hts_itr;
    bam1_t *aln = bam_init1();
    
    std::array<unsigned int, NUM_FILTER_REASONS> fillcode_to_num_alns;
    unsigned int num_pass_alns = 0;
    unsigned int num_iter_alns = 0;
    hts_itr = sam_itr_queryi(hts_idx, tid, fetch_tbeg, fetch_tend);
    while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) {
        bool isrc = false;
        bool isr2 = false;
        unsigned int tBeg = 0;
        unsigned int tEnd = 0;
        unsigned int num_seqs = 0;
        FilterReason filterReason = fill_isrc_isr2_beg_end_with_aln(isrc, isr2, tBeg, tEnd, num_seqs, 
                aln, fetch_tbeg, fetch_tend, min_alnlen, min_mapq, end2end);
        if (NOT_FILTERED == filterReason) {
            isrc_isr2_to_beg_count[isrc * 2 + isr2][tBeg + ARRPOS_MARGIN - fetch_tbeg] += 1;
            isrc_isr2_to_end_count[isrc * 2 + isr2][tEnd + ARRPOS_MARGIN - fetch_tbeg] += 1;
            num_pass_alns += 1;
            extended_inclu_beg_pos = min(extended_inclu_beg_pos, aln->core.pos);
            extended_exclu_end_pos = max(extended_exclu_end_pos, bam_endpos(aln));
        }
        fillcode_to_num_alns[filterReason]++;
        num_iter_alns += 1;
    }
    sam_itr_destroy(hts_itr);
    LOG(logINFO) << "Iteration 2 begins!";

    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_beg2bcenter = { inicount, inicount, inicount, inicount };
    for (unsigned int isrc_isr2 = 0; isrc_isr2 < 4; isrc_isr2++) {
        auto beg_to_count = isrc_isr2_to_beg_count[isrc_isr2];
        poscounter_to_pos2pcenter(isrc_isr2_to_beg2bcenter[isrc_isr2], beg_to_count);
    }
    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_end2ecenter = { inicount, inicount, inicount, inicount };
    for (unsigned int isrc_isr2 = 0; isrc_isr2 < 4; isrc_isr2++) {
        auto end_to_count = isrc_isr2_to_end_count[isrc_isr2];
        poscounter_to_pos2pcenter(isrc_isr2_to_end2ecenter[isrc_isr2], end_to_count);
    }
    
    unsigned int beg_peak_max = 0;
    for (auto beg_count : isrc_isr2_to_beg_count) {
        for (auto countval : beg_count) {
            beg_peak_max = max(beg_peak_max, countval);
        }
    }
    
    LOG(logINFO) << "Iteration 3 begins!";
    size_t alnidx = 0;
    hts_itr = sam_itr_queryi(hts_idx, tid, fetch_tbeg, fetch_tend);
    while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) {
        bool isrc = false;
        bool isr2 = false;
        unsigned int tBeg = 0;
        unsigned int tEnd = 0;
        unsigned int num_seqs = 0;
        FilterReason filterReason = fill_isrc_isr2_beg_end_with_aln(isrc, isr2, tBeg, tEnd, num_seqs, 
                aln, fetch_tbeg, fetch_tend, min_alnlen, min_mapq, end2end);
        if (NOT_FILTERED != filterReason) {
            continue;
        }
        //LOG(logINFO) << "Iteration 3.1 begins!";
        const char *qname = bam_get_qname(aln);
        const uint64_t qhash = strhash(qname);
        const char *umi_beg1 = strchr(qname,   '#');
        const char *umi_beg = ((NULL != umi_beg1) ? (umi_beg1 + 1) : (qname + aln->core.l_qname));
        const char *umi_end1 = strchr(umi_beg, '#');
        const char *umi_end = ((NULL != umi_end1) ? (umi_end1    ) : (qname + aln->core.l_qname)); 
       
        int is_umi_found = (umi_beg + 1 < umi_end); // UMI has at least one letter
        int is_duplex_found = 0;
        uint64_t umihash = 0;
        if (is_umi_found) {
            size_t umi_len = umi_end - umi_beg;
            size_t umi_half = (umi_end - umi_beg - 1) / 2;
            if ((umi_len % 2 == 1 ) && ( '+' == umi_beg[umi_half])) {
                uint64_t umihash_part1 = strnhash(umi_beg               , umi_half); // alpha
                uint64_t umihash_part2 = strnhash(umi_beg + umi_half + 1, umi_half); // beta
                umihash = ((isrc ^ isr2) ? hash2hash(umihash_part1, umihash_part2) : hash2hash(umihash_part2, umihash_part1));
                is_duplex_found++;
            } else {
                umihash = strnhash(umi_beg, umi_end-umi_beg); // + strnhash_rc(umi_beg, umi_end-umi_beg);
            }
        } else if ((aln->core.flag & 0x1) == 0 && umi_struct_string16.size() > 0) { // should be proton
            umihash = bam2umihash<false>(is_umi_found, aln, umi_struct_string16);
            if (!is_umi_found) {
                umihash = bam2umihash<true>(is_umi_found, aln, umi_struct_string16); 
            }
        }
        unsigned int isrc_isr2 = isrc * 2 + isr2;
        unsigned int beg2 = isrc_isr2_to_beg2bcenter[isrc_isr2][tBeg + ARRPOS_MARGIN - fetch_tbeg];
        // beg2 = (beg2 != 0 ? beg2 :tBeg);
        unsigned int end2 = isrc_isr2_to_end2ecenter[isrc_isr2][tEnd + ARRPOS_MARGIN - fetch_tbeg];
        // end2 = (end2 != 0 ? end2 :tEnd);
        unsigned int beg2count = isrc_isr2_to_beg_count[isrc_isr2][beg2];
        unsigned int end2count = isrc_isr2_to_end_count[isrc_isr2][end2];
        
        //LOG(logINFO) << "Iteration 3.2 begins!";
        unsigned int beg2surrcount = 0;
        for (int8_t i = -(int)ARRPOS_OUTER_RANGE; i < (int)ARRPOS_OUTER_RANGE+1; i++) {
            if (i > ARRPOS_INNER_RANGE || i < -ARRPOS_INNER_RANGE) {
                assert(i+(int)beg2 < fetch_size || !fprintf(stderr, "beg2 index %d + %d = %d is too big!", i, beg2, i+(int)beg2));
                unsigned int beg_count = isrc_isr2_to_beg_count.at(isrc_isr2).at(i+(int)beg2);
                beg2surrcount = max(beg2surrcount, beg_count);
            }
        }
        unsigned int end2surrcount = 0;
        for (int8_t i = -(int)ARRPOS_OUTER_RANGE; i < (int)ARRPOS_OUTER_RANGE+1; i++) {
            if (i > ARRPOS_INNER_RANGE && i < -ARRPOS_INNER_RANGE) {
                assert(i+(int)end2 < fetch_size || !fprintf(stderr, "end2 index %d + %d = %d is too big!", i, end2, i+(int)end2));
                unsigned int end_count = isrc_isr2_to_end_count.at(isrc_isr2).at(i+(int)end2);
                end2surrcount = max(end2surrcount, end_count);
            }
        }
        //LOG(logINFO) << "Iteration 3.3 begins!";
        double begfrac = (double)(beg2count) / (double)(beg2surrcount + 1);
        double endfrac = (double)(end2count) / (double)(end2surrcount + 1);
        double peakimba;
        double peakfrac;
        uint64_t umilabel;
        if ((CORRECTION_NONE == input_mode) || (CORRECTION_BASEQUAL == input_mode) || (CORRECTION_SINGLETON == input_mode)) {
            peakimba = 0;
            peakfrac = 0;
            umilabel = qhash;
        } else if (!is_umi_found) { // no UMI
            peakimba = 8;
            peakfrac = 8;
            if (begfrac > 16 || endfrac > 16) { // evidence for PCR
                pcrpassed += 1;
                umilabel = qhash;
            } else {
                umilabel = 0;
            }
        } else {
            peakimba = 4;
            peakfrac = 4;
            umilabel = umihash; 
        }
        
        //LOG(logINFO) << "Iteration 3.4 begins!";
        int begIsVariable = 0;
        int endIsVariable = 0;
        if (end2count  * peakimba < beg2count) {
            endIsVariable = 1; // special flag indicating no end as lots of reads begin at the same position but end at different positions
        } else if (beg2count * peakimba < end2count) {
            begIsVariable = 1;
        } else if (beg2count > end2count && begfrac > peakfrac) {
            endIsVariable = 2; // special flag indicating no end as the begin position attracts lots of reads
        } else if (beg2count < end2count && endfrac > peakfrac) {
            begIsVariable = 2; 
        } else {
            // do nothing
        }
        //LOG(logINFO) << "Iteration 3.5 begins!";
        
        assert(!(begIsVariable && endIsVariable));
        uint64_t umi3hash;
        if (begIsVariable) {
            umi3hash = hash2hash(umilabel, hash2hash(0, end2+1)); 
        } else if (endIsVariable) {
            umi3hash = hash2hash(umilabel, hash2hash(beg2+1, 0)); 
        } else {
            auto min2 = min(beg2, end2);
            auto max2 = max(beg2, end2);
            umi3hash = hash2hash(umilabel, hash2hash(min2+1, max2+1));
        }
        // umi3hash = umilabel;
        
        //LOG(logINFO) << "Iteration 3.6 begins!";
        unsigned int strand = isrc ^ isr2;
        int dflag = (is_duplex_found ? 2 : (is_umi_found ? 1 : 0));
        umi_to_strand_to_reads.insert(std::make_pair(umi3hash, std::make_pair(std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>(), dflag)));
        umi_to_strand_to_reads[umi3hash].first[strand].insert(std::make_pair(qhash, std::vector<bam1_t *>()));
        umi_to_strand_to_reads[umi3hash].first[strand][qhash].push_back(bam_dup1(aln));
        
        if ((ispowof2(alnidx+1) || ispowof2(num_pass_alns - alnidx)) && (beg_peak_max > 1000 || ispowof2(contig_ordinal+1))) {
            LOG(logINFO) << "  readname = " << qname << " ; "
                    << "alnidx/num_pass_alns = " << alnidx << " / " << num_pass_alns << " ; "
                    << "isrc = " << isrc << " ; "
                    << "isr2 = " << isr2 << " ; " 
                    << "original-range = " << tBeg << " to " << tEnd << " ; " 
                    << "adjusted-range = " << beg2 << " to " << end2 << " ; "
                    << "adjusted-counts = " << beg2count << " to " << end2count << " ; " 
                    << "adjusted-surrounding-counts = " << beg2surrcount << " to " << end2surrcount << " ; " 
                    << "molecular-barcode = " << (is_umi_found ? umihash : 0) << " ; "
                    << "begIsVariable = " << begIsVariable << " ; " 
                    << "endIsVariable = " << endIsVariable << " ; "
                    << "umi3hash = " << umi3hash << " ; "
                    << "strand = " << strand << " ; "
                    << "qhash = " << qhash;
        }
        alnidx += 1;
        //LOG(logINFO) << "Iteration 3.7 begins!";

    }
    sam_itr_destroy(hts_itr);
    
    bam_destroy1(aln);
    // sam_index_destroy(hts_idx); // TODO: reuse index
    // hts_idx_destroy(hts_idx);
    sam_close(sam_infile);
    LOG(logINFO) << "Final step of dedupping!" << std::endl;
    return (int)num_pass_alns;
}
#endif
