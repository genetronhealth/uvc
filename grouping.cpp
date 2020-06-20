#include "grouping.hpp"
#include "logging.hpp"

// position of 5' is the starting position, but position of 3' is unreliable without mate info.
const unsigned int ARRPOS_MARGIN = 1200;
const int8_t ARRPOS_OUTER_RANGE = 10;
const int8_t ARRPOS_INNER_RANGE = 3;

bool 
ispowof2(auto num) {
    return (num & (num-1)) == 0;
}

auto 
min(auto a, auto b) {
    return a < b ? a : b;
}

auto 
max(auto a, auto b) {
    return a > b ? a : b;
}

int 
bed_fname_to_contigs(
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> & tid_beg_end_e2e_vec,
        const std::string & bed_fname, 
        const bam_hdr_t *bam_hdr) {
    
    std::map<std::string, unsigned int> tname_to_tid;
    for (size_t i = 0; i < SIGN2UNSIGN(bam_hdr->n_targets); i++) {
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
        unsigned int nreads = 2 * (tend - tbeg);
        while (linestream.good()) {
            linestream >> token;
            if (token.find("end-to-end") != std::string::npos) {
                end2end = true;
            }
            if (token.find("NumberOfReadsInThisInterval") != std::string::npos) {
                linestream >> nreads;
            }
        }
        tid_beg_end_e2e_vec.push_back(std::make_tuple(tname_to_tid[tname], tbeg, tend, end2end, nreads)); // assume 100*1000 reads fall into this region
    }
    return 0;
}

int 
SamIter::iternext(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> & tid_beg_end_e2e_vec) {
    int ret = 0;
    unsigned int nreads_tot = 0;
    unsigned int region_tot = 0;
    unsigned int nreads_max = 0;
    unsigned int region_max = 0;
    
    if (NOT_PROVIDED != region_bed_fname) {
        for (; this->_bedregion_idx < this->_tid_beg_end_e2e_vec.size(); this->_bedregion_idx++) {
            ret++;
            const auto & bedreg = (this->_tid_beg_end_e2e_vec[this->_bedregion_idx]);
            tid_beg_end_e2e_vec.push_back(bedreg);
            nreads_tot += std::get<4>(bedreg);
            region_tot += std::get<2>(bedreg) - std::get<1>(bedreg);
            if (nreads_tot > (2000*1000 * nthreads) || region_tot > (1000*1000 * nthreads)) {
                this->_bedregion_idx++;
                return ret;
            }
        }
        return ret;
    }
    while (    (NULL == sam_idx && (sam_read1(this->sam_infile, this->samheader, alnrecord) >= 0))
            || (NULL != sam_idx && (sam_itr_next(this->sam_infile, this->sam_itr, alnrecord) >= 0))) {
        ret++;
        if (BAM_FUNMAP & alnrecord->core.flag) {
            continue;
        }
        bool is_uncov = (SIGN2UNSIGN(alnrecord->core.tid) != tid || SIGN2UNSIGN(alnrecord->core.pos) > tend);
        if (UINT_MAX == endingpos) {
            uint64_t n_overlap_positions = min(SIGN2UNSIGN(64), (SIGN2UNSIGN(16) + tend - min(tend, SIGN2UNSIGN(alnrecord->core.pos))));
            uint64_t npositions = (tend - min(tbeg, tend));
            bool has_many_positions = npositions > n_overlap_positions * (1024);
            bool has_many_reads = nreads > n_overlap_positions * (1024 * 2);
            if (has_many_positions || has_many_reads) {
                endingpos = SIGN2UNSIGN(max(bam_endpos(alnrecord), 
                        min(alnrecord->core.pos, alnrecord->core.mpos) + min(abs(alnrecord->core.isize), (int)ARRPOS_MARGIN)) 
                        + (int)(ARRPOS_OUTER_RANGE * 2));
            }
        }
        next_nreads += (SIGN2UNSIGN(bam_endpos(alnrecord)) > endingpos ? 1 : 0);
        if (is_uncov || endingpos < SIGN2UNSIGN(alnrecord->core.pos)) {
            region_max = max(region_max, tend - tbeg);
            nreads_max = max(nreads_max, nreads);
            region_tot += tend - tbeg;
            nreads_tot += nreads;
            auto prev_nreads = next_nreads;
            if (tid != UINT_MAX) {
                tid_beg_end_e2e_vec.push_back(std::make_tuple(tid, tbeg, tend, false, nreads));
                endingpos = UINT_MAX;
                next_nreads = 0;
            }
            tid = alnrecord->core.tid;
            if (is_uncov) {
                tbeg = SIGN2UNSIGN(alnrecord->core.pos);
                tend = SIGN2UNSIGN(bam_endpos(alnrecord));
            } else {
                tbeg = tend;
                tend = max(tbeg, SIGN2UNSIGN(bam_endpos(alnrecord))) + 1;
            }
            nreads = prev_nreads;
            nreads += 1;
            if (nreads_tot > (2000*1000 * nthreads) || region_tot > (1000*1000 * nthreads)) {
                return ret;
            }
        } else {
            tend = max(tend, SIGN2UNSIGN(bam_endpos(alnrecord)));
            nreads += 1;

        }
    }
    if (tid != UINT_MAX) {
        tid_beg_end_e2e_vec.push_back(std::make_tuple(tid, tbeg, tend, false, nreads));
    }
    return ret;
}

int
samfname_to_tid_to_tname_tseq_tup_vec(
        std::vector<std::tuple<std::string, unsigned int>> & tid_to_tname_tseqlen_tuple_vec, 
        const std::string & bam_input_fname) {
    
    tid_to_tname_tseqlen_tuple_vec.clear();
    samFile *sam_infile = sam_open(bam_input_fname.c_str(), "r");
    bam_hdr_t * samheader = sam_hdr_read(sam_infile);
    tid_to_tname_tseqlen_tuple_vec.reserve(samheader->n_targets);
    for (int tid = 0; tid < samheader->n_targets; tid++) {
        tid_to_tname_tseqlen_tuple_vec.push_back(std::make_tuple(std::string(samheader->target_name[tid]), samheader->target_len[tid]));
    }
    bam_hdr_destroy(samheader);
    sam_close(sam_infile);
    return 0;
}

int 
sam_fname_to_contigs(
        std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> & tid_beg_end_e2e_vec,
        std::vector<std::tuple<std::string, unsigned int>> & tid_to_tname_tlen_tuple_vec,
        const std::string & input_bam_fname, 
        const std::string & bed_fname) {
    
    tid_beg_end_e2e_vec.clear();
    tid_to_tname_tlen_tuple_vec.clear();
    samFile *sam_infile = sam_open(input_bam_fname.c_str(), "r");
    bam_hdr_t * samheader = sam_hdr_read(sam_infile);
    tid_to_tname_tlen_tuple_vec.reserve(samheader->n_targets);
    
    for (int tid = 0; tid < samheader->n_targets; tid++) {
        tid_to_tname_tlen_tuple_vec.push_back(std::make_tuple(std::string(samheader->target_name[tid]), samheader->target_len[tid]));
    }
    int ret;
    if (bed_fname != NOT_PROVIDED) {
        ret = bed_fname_to_contigs(tid_beg_end_e2e_vec, bed_fname, samheader);    
    } else {
        ret = 0;
        unsigned int endingpos = UINT_MAX;
        unsigned int tid = UINT_MAX;
        unsigned int tbeg = -1;
        unsigned int tend = -1;
        uint64_t nreads = 0;
        uint64_t next_nreads = 0;
        bam1_t *alnrecord = bam_init1();
        while (sam_read1(sam_infile, samheader, alnrecord) >= 0) {
            ret++;
            if (BAM_FUNMAP & alnrecord->core.flag) {
                continue;
            }
            bool is_uncov = (SIGN2UNSIGN(alnrecord->core.tid) != tid || SIGN2UNSIGN(alnrecord->core.pos) > tend);
            if (UINT_MAX == endingpos) {
                uint64_t n_overlap_positions = min(SIGN2UNSIGN(64), (SIGN2UNSIGN(16) + tend - min(tend, SIGN2UNSIGN(alnrecord->core.pos))));
                uint64_t npositions = (tend - min(tbeg, tend));
                bool has_many_positions = npositions > n_overlap_positions * (1024);
                bool has_many_reads = nreads > n_overlap_positions * (1024 * 2);
                if (has_many_positions || has_many_reads) {
                    endingpos = SIGN2UNSIGN(max(bam_endpos(alnrecord), 
                            min(alnrecord->core.pos, alnrecord->core.mpos) + min(abs(alnrecord->core.isize), (int)ARRPOS_MARGIN)))
                            + (ARRPOS_OUTER_RANGE * 2);
                }
            }
            next_nreads += (SIGN2UNSIGN(bam_endpos(alnrecord)) > endingpos ? 1 : 0);
            if (is_uncov || (endingpos < SIGN2UNSIGN(alnrecord->core.pos))) {
                auto prev_nreads = next_nreads;
                if (tid != UINT_MAX) {
                    tid_beg_end_e2e_vec.push_back(std::make_tuple(tid, tbeg, tend, false, nreads));
                    endingpos = UINT_MAX;
                    next_nreads = 0;
                }
                tid = SIGN2UNSIGN(alnrecord->core.tid);
                if (is_uncov) {
                    tbeg = alnrecord->core.pos;
                    tend = bam_endpos(alnrecord);
                } else {
                    tbeg = tend;
                    tend = max(tbeg, SIGN2UNSIGN(bam_endpos(alnrecord))) + SIGN2UNSIGN(1);
                }
                nreads = prev_nreads;
            }
            tend = max(tend, SIGN2UNSIGN(bam_endpos(alnrecord)));
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
        const bam1_t *aln, const uint32_t fetch_tbeg, const uint32_t fetch_tend,
        const unsigned int min_mapq, const unsigned int min_alnlen, const bool end2end, const bool is_pair_end_merge_enabled) {
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
    if (SIGN2UNSIGN(bam_endpos(aln) - aln->core.pos) < min_alnlen) {
        return LOW_ALN_LEN; // alignment length too short
    }

    isrc = ((aln->core.flag & 0x10) == 0x10);
    isr2 = ((aln->core.flag & 0x80) == 0x80 && (aln->core.flag & 0x1) == 0x1);
    if (!is_pair_end_merge_enabled) { isr2 = false; }
    const uint32_t begpos = aln->core.pos;
    const uint32_t endpos = bam_endpos(aln) - 1;
    if ((!is_pair_end_merge_enabled) 
            || ((aln->core.flag & 0x1) == 0) 
            // || ((aln->core.flag & 0x2) == 0) // having this line causes problems to SRR2556939_chr3_178936090_178936092
            || (aln->core.flag & 0x8) 
            || (0 == aln->core.isize) 
            || (((unsigned int)abs(aln->core.isize)) >= (ARRPOS_MARGIN))) {
        tBeg = (isrc ? endpos : begpos);
        tEnd = (isrc ? begpos : endpos);
        num_seqs = 1;
    } else {
        auto tBegP1 = min(begpos, SIGN2UNSIGN(aln->core.mpos));
        auto tEndP1 = tBegP1 + abs(aln->core.isize) - 1;
        bool strand = (isrc ^ isr2);
        tBeg = (strand ? tEndP1 : tBegP1); 
        tEnd = (strand ? tBegP1 : tEndP1);
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

unsigned int 
unsigned_diff(unsigned int a, unsigned int b) {
    return (a > b ? a - b : b - a);
}

int 
poscounter_to_pos2pcenter(
              std::vector<unsigned int> & pos_to_center_pos,
        const std::vector<unsigned int> & pos_to_count, 
        const unsigned int dedup_center_mult) {
    
    for (size_t locov_pos = ARRPOS_INNER_RANGE; locov_pos < pos_to_count.size() - ARRPOS_INNER_RANGE; locov_pos++) {
        unsigned int locov_count = pos_to_count[locov_pos];
        pos_to_center_pos[locov_pos] = locov_pos; // identity mapping by default
        unsigned int max_count = locov_count;
        // check if inner_pos is attracted by outer position
        for (size_t hicov_pos = locov_pos - ARRPOS_INNER_RANGE; hicov_pos < locov_pos + ARRPOS_INNER_RANGE + 1; hicov_pos++) {
            unsigned int hicov_count = pos_to_count[hicov_pos];
            if (hicov_count > max_count && hicov_count > locov_count * pow(dedup_center_mult, unsigned_diff(locov_pos, hicov_pos))) {
                pos_to_center_pos[locov_pos] = hicov_pos;
                max_count = hicov_count;
            }
        }
    }
    return 0;
}

// one-way converion of data into hash values
// https://en.wikipedia.org/wiki/Universal_hashing#Hashing_strings
template <class T> 
uint64_t 
strnhash(const T *str, size_t n) {
    uint64_t ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * 31UL + ((uint64_t)str[i]);
    }
    return ret;
}

template <class T> 
uint64_t 
strnhash_rc(const T *str, size_t n) {
    uint64_t ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * 31UL + THE_REV_COMPLEMENT.data[((uint64_t)str[n-i-(size_t)1])];
    }
    return ret;
}

template<class T> 
uint64_t 
strhash(const T *str) {
    return strnhash(str, SIZE_MAX);
}

uint64_t 
hash2hash(uint64_t hash1, uint64_t hash2) {
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
    return 0;
}

int 
apply_baq(bam1_t *aln, const unsigned int baq_per_aligned_base) {
    if (aln->core.l_qseq < (10*2)) { return -1; }
    for (unsigned int i = 0; i < 10; i++) { 
        auto j = aln->core.l_qseq - 1 - i;
        (bam_get_qual(aln))[i] = min((bam_get_qual(aln))[i], baq_per_aligned_base * (i+1));
        (bam_get_qual(aln))[j] = min((bam_get_qual(aln))[j], baq_per_aligned_base * (i+1));
    }
    return 0;
}

int 
fill_strand_umi_readset_with_strand_to_umi_to_reads(
        std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> &umi_strand_readset,
        std::map<uint64_t, std::pair<std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>, int>> &umi_to_strand_to_reads,
        unsigned int baq_per_aligned_base
        ) {
    for (auto & umi_to_strand_to_reads_element : umi_to_strand_to_reads) {
        const auto strand_to_reads = umi_to_strand_to_reads_element.second.first;
        const auto dflag = umi_to_strand_to_reads_element.second.second;
        umi_strand_readset.push_back(std::make_pair(std::array<std::vector<std::vector<bam1_t *>>, 2>(), dflag));
        for (unsigned int strand = 0; strand < 2; strand++) {
            for (auto read : strand_to_reads[strand]) {
                const std::vector<bam1_t *> alns = read.second;
                umi_strand_readset.back().first[strand].push_back(std::vector<bam1_t *>());
                for (auto aln : alns) {
                    apply_baq(aln, baq_per_aligned_base);
                    umi_strand_readset.back().first[strand].back().push_back(aln);
                }
            }
        }
    }
    return 0;
};

template <bool is_rc>
uint64_t 
bam2umihash(int & is_umi_found, const bam1_t *aln, const std::vector<uint8_t> & UMI_STRUCT, const int max_begin_diff_umi2read = 5) {
    LOG(logDEBUGx1) << "Going over " << UMI_STRUCT.size() << " bases in the pattern";

    auto *bamseq = bam_get_seq(aln);
    
    for (int i = 0; i < max_begin_diff_umi2read; i++) {
        size_t patpos = 0;
        uint64_t umihash = 0;
        for (int j = i; j < aln->core.l_qseq && patpos < UMI_STRUCT.size(); j++) {
            char int4base;
            if (is_rc) {
                char int4base2 = bam_seqi(bamseq, aln->core.l_qseq - 1 - j);
                int4base = THE_REV_COMPLEMENT.table16[(int8_t)int4base2];
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

std::array<unsigned int, 3>
bamfname_to_strand_to_familyuid_to_reads(
        std::map<uint64_t, std::pair<std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>, int>> &umi_to_strand_to_reads,
        unsigned int & extended_inclu_beg_pos, 
        unsigned int & extended_exclu_end_pos,
        const std::string input_bam_fname, 
        unsigned int tid, 
        unsigned int fetch_tbeg, 
        unsigned int fetch_tend, 
        bool end2end, 
        unsigned int min_mapq, 
        unsigned int min_alnlen, 
        unsigned int regionbatch_ordinal, 
        unsigned int regionbatch_tot_num,
        const std::string UMI_STRUCT_STRING, 
        const hts_idx_t * hts_idx,
        const bool is_molecule_tag_enabled,
        const bool is_pair_end_merge_enabled,
        bool disable_duplex,
        size_t thread_id,
        unsigned int dedup_center_mult,
        unsigned int dedup_amplicon_count_to_surrcount_frac,
        unsigned int dedup_yes_umi_2ends_peak_frac,
        unsigned int dedup_non_umi_2ends_peak_frac,
        unsigned int dedup_yes_umi_2ends_flat_perc,
        unsigned int dedup_non_umi_2ends_flat_perc,
        bool always_log,
        unsigned int specialflag) {
    assert (fetch_tend > fetch_tbeg);
    
    const bool should_log = (ispowof2(regionbatch_ordinal+1) || ispowof2(regionbatch_tot_num - regionbatch_ordinal));
    std::vector<uint8_t> umi_struct_string16;
    for (auto ch : UMI_STRUCT_STRING) {
        umi_struct_string16.push_back(seq_nt16_table[(int8_t)ch]);
    }
    for (auto base : umi_struct_string16) {
        LOG(logDEBUGx1) << "Base " << (int)base;
    }
    extended_inclu_beg_pos = INT32_MAX;
    extended_exclu_end_pos = 0;
    
    unsigned int pcrpassed, umi_pcrpassed;
    pcrpassed = umi_pcrpassed = 0;
   
    samFile *sam_infile = sam_open(input_bam_fname.c_str(), "r");
    if (should_log) {
        LOG(logINFO) << "Thread " << thread_id << " started dedupping the chunk tid" << tid << ":" << fetch_tbeg << "-" << fetch_tend 
                << " (region no " << regionbatch_ordinal << "/" << regionbatch_tot_num << " in this batch)";
    }
    unsigned int fetch_size = fetch_tend - fetch_tbeg + (ARRPOS_MARGIN + ARRPOS_OUTER_RANGE) * 2;
    
    std::vector<unsigned int> inicount(fetch_size, 0);
    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_beg_count = {{ inicount, inicount, inicount, inicount }};
    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_end_count = {{ inicount, inicount, inicount, inicount }};
    
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
                aln, fetch_tbeg, fetch_tend, min_alnlen, min_mapq, end2end, is_pair_end_merge_enabled);
        if (!is_pair_end_merge_enabled) { assert(!isr2); }
        if (NOT_FILTERED == filterReason) {
            isrc_isr2_to_beg_count[isrc * 2 + isr2][tBeg + ARRPOS_MARGIN - fetch_tbeg] += 1;
            isrc_isr2_to_end_count[isrc * 2 + isr2][tEnd + ARRPOS_MARGIN - fetch_tbeg] += 1;
            num_pass_alns += 1;
            extended_inclu_beg_pos = min(extended_inclu_beg_pos, SIGN2UNSIGN(aln->core.pos));
            extended_exclu_end_pos = max(extended_exclu_end_pos, SIGN2UNSIGN(bam_endpos(aln)));
        }
        fillcode_to_num_alns[filterReason]++;
        num_iter_alns += 1;
    }
    sam_itr_destroy(hts_itr);
    
    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_beg2bcenter = {{ inicount, inicount, inicount, inicount }};
    for (unsigned int isrc_isr2 = 0; isrc_isr2 < 4; isrc_isr2++) {
        auto beg_to_count = isrc_isr2_to_beg_count[isrc_isr2];
        poscounter_to_pos2pcenter(isrc_isr2_to_beg2bcenter[isrc_isr2], beg_to_count, dedup_center_mult);
    }
    std::array<std::vector<unsigned int>, 4> isrc_isr2_to_end2ecenter = {{ inicount, inicount, inicount, inicount }};
    for (unsigned int isrc_isr2 = 0; isrc_isr2 < 4; isrc_isr2++) {
        auto end_to_count = isrc_isr2_to_end_count[isrc_isr2];
        poscounter_to_pos2pcenter(isrc_isr2_to_end2ecenter[isrc_isr2], end_to_count, dedup_center_mult);
    }
    
    unsigned int beg_peak_max = 0;
    for (auto beg_count : isrc_isr2_to_beg_count) {
        for (auto countval : beg_count) {
            beg_peak_max = max(beg_peak_max, countval);
        }
    }
    
    size_t alnidx = 0;
    hts_itr = sam_itr_queryi(hts_idx, tid, fetch_tbeg, fetch_tend);
    while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) {
        bool isrc = false;
        bool isr2 = false;
        unsigned int tBeg = 0;
        unsigned int tEnd = 0;
        unsigned int num_seqs = 0;
        FilterReason filterReason = fill_isrc_isr2_beg_end_with_aln(isrc, isr2, tBeg, tEnd, num_seqs, 
                aln, fetch_tbeg, fetch_tend, min_alnlen, min_mapq, end2end, is_pair_end_merge_enabled);
        if (!is_pair_end_merge_enabled) { assert(!isr2); }
        if (NOT_FILTERED != filterReason) {
            continue;
        }
        const char *qname = bam_get_qname(aln);
        const uint64_t qhash = strhash(qname);
        const char *umi_beg1 = strchr(qname,   '#');
        const char *umi_beg = ((NULL != umi_beg1) ? (umi_beg1 + 1) : (qname + aln->core.l_qname));
        const char *umi_end1 = strchr(umi_beg, '#');
        const char *umi_end = ((NULL != umi_end1) ? (umi_end1    ) : (qname + aln->core.l_qname)); 
       
        int is_umi_found = ((umi_beg + 1 < umi_end) && is_molecule_tag_enabled); // UMI has at least one letter
        int is_duplex_found = 0;
        uint64_t umihash = 0;
        if (is_umi_found) {
            size_t umi_len = umi_end - umi_beg;
            size_t umi_half = (umi_end - umi_beg - 1) / 2;
            if ((umi_len % 2 == 1 ) && ( '+' == umi_beg[umi_half]) && (!disable_duplex)) {
                uint64_t umihash_part1 = strnhash(umi_beg               , umi_half); // alpha
                uint64_t umihash_part2 = strnhash(umi_beg + umi_half + 1, umi_half); // beta
                umihash = ((isrc ^ isr2) ? hash2hash(umihash_part1, umihash_part2) : hash2hash(umihash_part2, umihash_part1));
                is_duplex_found++;
            } else {
                umihash = strnhash(umi_beg, umi_end-umi_beg);
            }
        } else if ((aln->core.flag & 0x1) == 0 && umi_struct_string16.size() > 0) { // should be proton
            umihash = bam2umihash<false>(is_umi_found, aln, umi_struct_string16);
            if (!is_umi_found) {
                umihash = bam2umihash<true>(is_umi_found, aln, umi_struct_string16); 
            }
        }
        unsigned int isrc_isr2 = isrc * 2 + isr2;
        unsigned int beg2 = isrc_isr2_to_beg2bcenter[isrc_isr2][tBeg + ARRPOS_MARGIN - fetch_tbeg];
        unsigned int end2 = isrc_isr2_to_end2ecenter[isrc_isr2][tEnd + ARRPOS_MARGIN - fetch_tbeg];
        unsigned int beg2count = isrc_isr2_to_beg_count[isrc_isr2][beg2];
        unsigned int end2count = isrc_isr2_to_end_count[isrc_isr2][end2];
        
        unsigned int beg2surrcount = 0;
        for (int8_t i = -(int)ARRPOS_OUTER_RANGE; i < (int)ARRPOS_OUTER_RANGE+1; i++) {
            if (i > ARRPOS_INNER_RANGE || i < -ARRPOS_INNER_RANGE) {
                assert(i+(int)beg2 < (int)fetch_size || !fprintf(stderr, "beg2 index %d + %d = %d is too big!", i, (int)beg2, i+(int)beg2));
                unsigned int beg_count = isrc_isr2_to_beg_count.at(isrc_isr2).at(i+(int)beg2);
                beg2surrcount = max(beg2surrcount, beg_count);
            }
        }
        unsigned int end2surrcount = 0;
        for (int8_t i = -(int)ARRPOS_OUTER_RANGE; i < (int)ARRPOS_OUTER_RANGE+1; i++) {
            if (i > ARRPOS_INNER_RANGE && i < -ARRPOS_INNER_RANGE) {
                assert(i+(int)end2 < (int)fetch_size || !fprintf(stderr, "end2 index %d + %d = %d is too big!", i, (int)end2, i+(int)end2));
                unsigned int end_count = isrc_isr2_to_end_count.at(isrc_isr2).at(i+(int)end2);
                end2surrcount = max(end2surrcount, end_count);
            }
        }
        double begfrac = (double)(beg2count) / (double)(beg2surrcount + 1);
        double endfrac = (double)(end2count) / (double)(end2surrcount + 1);
        
        const bool is_assay_amplicon = (begfrac > dedup_amplicon_count_to_surrcount_frac || endfrac > dedup_amplicon_count_to_surrcount_frac);
        const uint64_t umilabel = (is_umi_found ? umihash : (is_assay_amplicon ? qhash : 0));
        pcrpassed += is_assay_amplicon;
        
        unsigned int beg_peak_flag = 0;
        unsigned int end_peak_flag= 0;
        const uint64_t peakimba = (is_umi_found ? dedup_yes_umi_2ends_peak_frac : dedup_non_umi_2ends_peak_frac);
        const uint64_t flatperc = (is_umi_found ? dedup_yes_umi_2ends_flat_perc : dedup_non_umi_2ends_flat_perc);
        if (beg2count > (beg2surrcount + 1) * peakimba) {
            beg_peak_flag |= 0x1;
        }
        if (beg2count * flatperc > (end2count + 1) * (uint64_t)100UL) {
            beg_peak_flag |= 0x2;
        }
        if (end2count > (end2surrcount + 1) * peakimba) {
            end_peak_flag |= 0x1;
        }
        if (end2count * flatperc > (beg2count + 1) * (uint64_t)100UL) {
            end_peak_flag |= 0x2;
        }
        
        uint64_t umi3hash;
        if        ((0x3 == end_peak_flag)) {
            umi3hash = hash2hash(umilabel, hash2hash(0, end2+1)); 
        } else if ((0x3 == beg_peak_flag)) {
            umi3hash = hash2hash(umilabel, hash2hash(beg2+1, 0)); 
        } else {
            auto min2 = min(beg2, end2);
            auto max2 = max(beg2, end2);
            umi3hash = hash2hash(umilabel, hash2hash(min2+1, max2+1));
        }
        
        unsigned int strand = (isrc ^ isr2);
        int dflag = (is_duplex_found ? 2 : (is_umi_found ? 1 : 0));
        umi_to_strand_to_reads.insert(std::make_pair(umi3hash, std::make_pair(std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>(), dflag)));
        umi_to_strand_to_reads[umi3hash].first[strand].insert(std::make_pair(qhash, std::vector<bam1_t *>()));
        umi_to_strand_to_reads[umi3hash].first[strand][qhash].push_back(bam_dup1(aln));
        
        const bool should_log_read = (ispowof2(alnidx+1) || ispowof2(num_pass_alns - alnidx));
        if (!is_pair_end_merge_enabled) { assert(!isr2); }
        if ((should_log_read && (beg_peak_max >= 2000 || should_log)) || always_log) {
            LOG(logINFO) << "Thread " << thread_id << " ; readname = " << qname << " ; " 
                    << "tid = " << aln->core.tid << " ; "
                    << "alnidx/num_pass_alns = " << alnidx << "/" << num_pass_alns << " ; "
                    << "isrc = " << isrc << " ; "
                    << "isr2 = " << isr2 << " ; "
                    << "num_seqs = " << num_seqs << " ; "
                    << "original-range = " << tBeg << " to " << tEnd << " ; " 
                    << "adjusted-range = " << beg2 << " to " << end2 << " ; "
                    << "adjusted-counts = " << beg2count << " to " << end2count << " ; " 
                    << "adjusted-surrounding-counts = " << beg2surrcount << " to " << end2surrcount << " ; " 
                    << "molecular-barcode = " << (is_umi_found ? umihash : 0) << " ; "
                    << "beg_peak_flag = " << beg_peak_flag << " ; " 
                    << "end_peak_flag = " << end_peak_flag << " ; "
                    << "umi3hash = " << umi3hash << " ; "
                    << "strand = " << strand << " ; "
                    << "qhash = " << qhash;
        }
        alnidx += 1;
    }
    sam_itr_destroy(hts_itr);
    
    bam_destroy1(aln);
    // sam_index_destroy(hts_idx); // TODO: reuse index if possible
    // hts_idx_destroy(hts_idx);
    sam_close(sam_infile);
    if (should_log) { LOG(logINFO) << "Thread " << thread_id << " finished dedupping."; }
    return std::array<unsigned int, 3>({num_pass_alns, pcrpassed, umi_pcrpassed});
}

