#include "grouping.hpp"
#include "logging.hpp"

//#define MAX_NUM_REF_BASES (1000*1000)
//#define MAX_NUM_READS (2000*1000)

// at 150*16 average sequencing depth, the two below amount of bytes are approx equal to each other.
#define NUM_BYTES_PER_GENOMIC_POS ((int64_t)(1024*8)) // estimated
#define NUM_BYTES_PER_READ ((int64_t)(512)) // estimated

#define UPDATE_MIN(a, b) ((a) = MIN((a), (b)));
// x1 and y1 are inclusive whereas x2 and y2 are exclusive
#define ARE_INTERVALS_OVERLAPPING(x1, x2, y1, y2) (!(((x2) <= (y1)) || (((y2) <= (x1)))))

// position of 5' is the starting position, but position of 3' is unreliable without mate info.
const uvc1_readpos_t ARRPOS_MARGIN = MAX_INSERT_SIZE;
const uvc1_readpos_t ARRPOS_OUTER_RANGE = 10;
const uvc1_readpos_t ARRPOS_INNER_RANGE = 3;

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

int 
bed_fname_to_contigs(
        std::vector<bedline_t> & tid_beg_end_e2e_vec,
        const std::string & bed_fname, 
        const bam_hdr_t *bam_hdr) {
    
    std::map<std::string, uvc1_refgpos_t> tname_to_tid;
    for (uvc1_refgpos_t i = 0; i < bam_hdr->n_targets; i++) {
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
        uvc1_refgpos_t tbeg;
        uvc1_refgpos_t tend;
        linestream >> tname;
        linestream >> tbeg;
        linestream >> tend;
        if (!(tbeg < tend)) {
            std::cerr << "The bedfile " << bed_fname << " does not have its end after its start at: " << tname << "\t" << tbeg << "\t" << tend;
            exit (16);
        }
        if (tname_to_tid.find(tname) == tname_to_tid.end()) {
            std::cerr << "The reference template name " << tname << " from the bedfile " << bed_fname << " is not in the input sam file";
            exit (17);
        }
        bool end2end = false;
        std::string token;
        uvc1_readnum_t nreads = 2 * (tend - tbeg);
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

int64_t
SamIter::iternext(std::vector<bedline_t> &tid_beg_end_e2e_vec, const CommandLineArgs &paramset) {
    uvc1_readnum_big_t total_n_reads = 0;
    uvc1_refgpos_t total_n_ref_positions = 0;
    
    if (NOT_PROVIDED != region_bed_fname) {
        for (; this->_bedregion_idx < this->_tid_beg_end_e2e_vec.size(); this->_bedregion_idx++) {
            //ret++;
            const auto & bedreg = (this->_tid_beg_end_e2e_vec[this->_bedregion_idx]);
            tid_beg_end_e2e_vec.push_back(bedreg);
            const auto bed_tid = std::get<0>(bedreg);
            const auto bed_beg = std::get<1>(bedreg);
            const auto bed_end = std::get<2>(bedreg);
            int64_t region_n_reads = bed_in_avg_sequencing_DP; // Please note that left-over reads from the previoous iteration are ignored
            if (-1 == bed_in_avg_sequencing_DP) {
                region_n_reads = 0;
                while (    (NULL == sam_idx && (sam_read1(this->sam_infile, this->samheader, alnrecord) >= 0))
                        || (NULL != sam_idx && (sam_itr_next(this->sam_infile, this->sam_itr, alnrecord) >= 0))) {
                    if ((bed_tid == alnrecord->core.tid) && 
                            ARE_INTERVALS_OVERLAPPING(bed_beg, bed_end, alnrecord->core.pos, bam_endpos(alnrecord))) {
                        region_n_reads++;
                    } else if ((bed_tid < alnrecord->core.tid) || ((bed_tid == alnrecord->core.tid) && (bed_end <= alnrecord->core.pos))) {
                        break;
                    }
                }
            }
            total_n_reads += region_n_reads;
            total_n_ref_positions += bed_end - bed_beg;
            const size_t tot_n_bytes_used_by_reads = INT64MUL(total_n_reads, NUM_BYTES_PER_READ);
            const size_t tot_n_bytes_used_by_rposs = INT64MUL(total_n_ref_positions, 256*4); // estimated from the htslib specs of VCF
            const bool tot_has_too_much_mem = ((tot_n_bytes_used_by_reads + tot_n_bytes_used_by_rposs) 
                        > (((uint64_t) 1024*1024) * this->mem_per_thread * this->nthreads));
            if (tot_has_too_much_mem) {
                this->_bedregion_idx++;
                return total_n_reads;
            }
        }
    } else {
        
        uvc1_refgpos_t block_tid = this->last_it_tid;
        uvc1_refgpos_t block_beg = this->last_it_beg;
        uvc1_refgpos_t block_running_end = this->last_it_end;

        uvc1_readnum_big_t region_n_reads = 0;
        std::set<std::string> visited_qnames;
        
        int sam_read_ret = -1;
        do {
            sam_read_ret = ((NULL != sam_idx) ? (sam_itr_next(this->sam_infile, this->sam_itr, alnrecord))
                : (sam_read1(this->sam_infile, this->samheader, alnrecord)));
            if ((-1 == this->last_it_tid) && (sam_read_ret < 0)) {
                LOG(logWARNING) << "Encountered error while iterating over the first BAM record in the file " << this->input_bam_fname;
                break;
            }
            if (BAM_FUNMAP & alnrecord->core.flag) { continue; }
            NORM_INSERT_SIZE(alnrecord);
            const auto curr_tid = alnrecord->core.tid;
            const auto curr_beg = alnrecord->core.pos;
            const auto curr_end = bam_endpos(alnrecord);
            const size_t n_bytes_used_by_reads = INT64MUL(total_n_reads, NUM_BYTES_PER_READ);
            const size_t n_bytes_used_by_poss = INT64MUL(block_running_end - block_beg, NUM_BYTES_PER_GENOMIC_POS);
            const bool is_template_changed = (curr_tid != block_tid); 
            const bool is_far_jumped = ((curr_tid == block_tid) && (block_running_end + MAX_INSERT_SIZE < curr_beg));
            const bool has_too_much_mem = ((n_bytes_used_by_reads + n_bytes_used_by_poss) > (((uint64_t)1024*1024) / NUM_WORKING_UNITS_PER_THREAD) * this->mem_per_thread);
            
            if (0 == (total_n_reads % (1024*1024))) {
                LOG(logDEBUG4) << "ReadName=" << bam_get_qname(alnrecord) 
                        << " TID=" << (alnrecord->core.tid)
                        << " POS=" << (alnrecord->core.pos)
                        << " is_template_changed=" << is_template_changed 
                        << " is_far_jumped=" << is_far_jumped 
                        << " has_too_much_mem=" 
                        << has_too_much_mem 
                        << " sam_read_ret=" << sam_read_ret
                        << " total_n_reads=" << total_n_reads
                        << " approx total_n_ref_bases=" << (block_running_end - block_beg);
            }

            if (is_template_changed || is_far_jumped || has_too_much_mem || (sam_read_ret < 0)) {
                // flush to output due to ref-genome segmentation
                const int64_t div = 1; // Please note that MGVCF_REGION_MAX_SIZE will be used later so that div is set to one here.
                int64_t block_norm_end = MIN((((block_running_end + div - 1) / div) * div), (uvc1_refgpos_t)(this->samheader->target_len[block_tid]));
                const bool is_1st_read = (-1 == block_tid);                    
                const bool is_min_DP_failed_1_ = ((UNSIGN2SIGN(visited_qnames.size()) < paramset.min_altdp_thres) && is_far_jumped);
                const bool is_min_DP_failed = ((NOT_PROVIDED == paramset.vcf_tumor_fname) && is_min_DP_failed_1_ && (!paramset.should_output_all));
                const bool is_block_zero_sized =  (block_beg >= block_norm_end); 
                if (!is_1st_read && !is_min_DP_failed && !is_block_zero_sized) {
                    const auto bed_region = std::make_tuple(block_tid, block_beg, block_norm_end, false, region_n_reads);
                    tid_beg_end_e2e_vec.push_back(bed_region);
                    total_n_ref_positions += (block_norm_end - block_beg);
                    total_n_reads += region_n_reads;
                    region_n_reads = 0;
                    visited_qnames.clear();
                }
                block_tid = curr_tid;
                const auto new_block_beg = MAX(block_beg, (curr_beg / div) * div); // skip over non-covered bases
                block_beg = (is_template_changed ? curr_beg : MAX(new_block_beg, block_norm_end));
                const size_t tot_n_bytes_used_by_reads = INT64MUL(total_n_reads, NUM_BYTES_PER_READ);
                const size_t tot_n_bytes_used_by_rposs = INT64MUL(total_n_ref_positions, 256*4); // estimated from the htslib specs of VCF
                const bool tot_has_too_much_mem = ((tot_n_bytes_used_by_reads + tot_n_bytes_used_by_rposs) 
                        > (((uint64_t) 1024*1024) * this->mem_per_thread * this->nthreads));
                if (tot_has_too_much_mem) {
                    this->last_it_tid = block_tid;
                    this->last_it_beg = block_beg;
                    this->last_it_end = MAX(block_beg, block_norm_end);
                    return (total_n_reads);
                }
            }
            // can check for stronger condition
            if (UNSIGN2SIGN(visited_qnames.size()) <= paramset.min_altdp_thres) {
                visited_qnames.insert(bam_get_qname(alnrecord));
            }
            block_running_end = (is_template_changed ? curr_end : MAX(block_running_end, curr_end));
            region_n_reads++;
            
        } while (sam_read_ret >= 0);    
    }
    return total_n_reads;
}

int
samfname_to_tid_to_tname_tseq_tup_vec(
        std::vector<std::tuple<std::string, uvc1_refgpos_t>> & tid_to_tname_tseqlen_tuple_vec, 
        const std::string & bam_input_fname) {
    
    tid_to_tname_tseqlen_tuple_vec.clear();
    samFile *sam_infile = sam_open(bam_input_fname.c_str(), "r");
    bam_hdr_t * samheader = sam_hdr_read(sam_infile);
    tid_to_tname_tseqlen_tuple_vec.reserve(samheader->n_targets);
    for (uvc1_refgpos_t tid = 0; tid < UNSIGN2SIGN(samheader->n_targets); tid++) {
        tid_to_tname_tseqlen_tuple_vec.push_back(std::make_tuple(std::string(samheader->target_name[tid]), samheader->target_len[tid]));
    }
    bam_hdr_destroy(samheader);
    sam_close(sam_infile);
    return 0;
}

enum FilterReason {
    NOT_FILTERED,
    NOT_MAPPED,
    NOT_PRIMARY_ALN,
    LOW_MAPQ,
    LOW_ALN_LEN,
    LOW_ISIZE,
    HIGH_ISIZE,
    ZERO_ISIZE,
    OUT_OF_RANGE,
    NOT_END_TO_END,
    NUM_FILTER_REASONS
};

template <class T>
enum FilterReason
fill_isrc_isr2_beg_end_with_aln(bool & isrc, bool & isr2, uvc1_refgpos_t & tBeg, uvc1_refgpos_t & tEnd, T &num_seqs,
        const bam1_t *aln, const uvc1_refgpos_t fetch_tbeg, const uvc1_refgpos_t fetch_tend,
        const uvc1_qual_t min_mapqual, 
        const uvc1_readpos_t min_aln_len, 
        const uvc1_readpos_t min_isize, 
        const uvc1_readpos_t max_isize, 
        const bool is_zero_isize_discarded,
        const bool end2end, const bool is_pair_end_merge_enabled) {
    num_seqs = 0;
    if (aln->core.flag & 0x4) {
        return NOT_MAPPED;
    }
    if ((aln->core.flag & 0x900) != 0) {
        return NOT_PRIMARY_ALN;
    }

    if (aln->core.qual < min_mapqual) { 
        return LOW_MAPQ;
    }
    if (SIGN2UNSIGN(bam_endpos(aln) - aln->core.pos) < min_aln_len) {
        return LOW_ALN_LEN;
    }
    if (0 == (aln->core.isize)) {
        if (is_zero_isize_discarded) {
            return ZERO_ISIZE;
        }
    } else {
        if (abs(aln->core.isize) < min_isize) {
            return LOW_ISIZE;
        }
        if (abs(aln->core.isize) > max_isize) {
            return HIGH_ISIZE;
        }
    }
    
    isrc = ((aln->core.flag & 0x10) == 0x10);
    isr2 = ((aln->core.flag & 0x80) == 0x80 && (aln->core.flag & 0x1) == 0x1);
    if (!is_pair_end_merge_enabled) { isr2 = false; }
    const auto begpos = aln->core.pos;
    const auto endpos = bam_endpos(aln) - 1;
    if ((!is_pair_end_merge_enabled) 
            || ((aln->core.flag & 0x1) == 0) 
            // || ((aln->core.flag & 0x2) == 0) // having this line causes problems to SRR2556939_chr3_178936090_178936092
            || (aln->core.flag & 0x8) 
            || (0 == aln->core.isize) 
            || ((abs(aln->core.isize)) >= (ARRPOS_MARGIN))) {
        tBeg = (isrc ? endpos : begpos);
        tEnd = (isrc ? begpos : endpos);
        num_seqs = 1;
    } else {
        auto tBegP1 = MIN(begpos, SIGN2UNSIGN(aln->core.mpos));
        auto tEndP1 = tBegP1 + abs(aln->core.isize) - 1;
        bool strand = (isrc ^ isr2);
        tBeg = (strand ? tEndP1 : tBegP1); 
        tEnd = (strand ? tBegP1 : tEndP1);
        num_seqs = 2;
    }
    auto tOrdBeg = MIN(tBeg, tEnd);
    auto tOrdEnd = MAX(tBeg, tEnd);
    if (tOrdBeg + (ARRPOS_MARGIN - ARRPOS_OUTER_RANGE) <= fetch_tbeg || fetch_tend - 1 + (ARRPOS_MARGIN - ARRPOS_OUTER_RANGE) <= tOrdEnd) {
        return OUT_OF_RANGE;
    }
    if (end2end && !(tOrdBeg <= fetch_tbeg && tOrdEnd >= fetch_tend)) {
        return NOT_END_TO_END;
    }
    return NOT_FILTERED;
}

uvc1_unsigned_int_t 
unsigned_diff(uvc1_unsigned_int_t a, uvc1_unsigned_int_t b) {
    return (a > b ? a - b : b - a);
}

int 
poscounter_to_pos2pcenter(
              std::vector<uvc1_refgpos_t> & pos_to_center_pos,
        const std::vector<uvc1_readnum_t> & pos_to_count, 
        const double dedup_center_mult) {
    
    for (uvc1_refgpos_t locov_pos = ARRPOS_INNER_RANGE; locov_pos < UNSIGN2SIGN(pos_to_count.size()) - ARRPOS_INNER_RANGE; locov_pos++) {
        auto locov_count = pos_to_count[locov_pos];
        pos_to_center_pos[locov_pos] = locov_pos;
        auto max_count = locov_count;
        // check if inner_pos is attracted by outer position
        for (auto hicov_pos = locov_pos - ARRPOS_INNER_RANGE; hicov_pos < locov_pos + ARRPOS_INNER_RANGE + 1; hicov_pos++) {
            auto hicov_count = pos_to_count[hicov_pos];
            if ((hicov_count > max_count) && ((hicov_count + 1) > (locov_count + 1) * pow(dedup_center_mult, unsigned_diff(locov_pos, hicov_pos)))) {
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
auto
strnhash(const T *str, size_t n) {
    uvc1_hash_t  ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * 31UL + ((uvc1_hash_t)str[i]);
    }
    return ret;
}

template <class T> 
auto
strnhash_rc(const T *str, size_t n) {
    uvc1_hash_t ret = 0;
    for (size_t i = 0; i < n && str[i]; i++) {
        ret = ret * 31UL + THE_REV_COMPLEMENT.data[((uvc1_hash_t)str[n-i-(size_t)1])];
    }
    return ret;
}

template<class T> 
uvc1_hash_t 
strhash(const T *str) {
    return strnhash(str, SIZE_MAX);
}

uvc1_hash_t 
hash2hash(uvc1_hash_t hash1, uvc1_hash_t hash2) {
    return hash1 * ((1UL<<(31UL)) - 1UL) + hash2;
}

int 
clean_fill_strand_umi_readset(
        std::vector<std::array<std::vector<std::vector<bam1_t *>>, 2>> &umi_strand_readset) {
    for (auto & strand_readset : umi_strand_readset) {
        for (int strand = 0; strand < 2; strand++) {
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
apply_bq_err_correction3(bam1_t *aln, const uvc1_qual_t assay_sequencing_BQ_max, const uvc1_qual_t assay_sequencing_BQ_inc) {
    if ((0 == aln->core.l_qseq) || (aln->core.flag & 0x4)) { return -1; }
    
    for (uvc1_readpos_t i = 0; i < aln->core.l_qseq; i++) {
        uvc1_qual_t bq = bam_get_qual(aln)[i];
        bam_get_qual(aln)[i] = MIN(bq + assay_sequencing_BQ_inc, assay_sequencing_BQ_max);
    }
    
    const auto cigar = bam_get_cigar(aln);
    const int strand = ((aln->core.flag & 0x10) ? 1 : 0);
    uvc1_readpos_t inclu_beg_poss[2] = {0, aln->core.l_qseq - 1};
    uvc1_readpos_t exclu_end_poss[2] = {aln->core.l_qseq, 0 - 1};
    uvc1_readpos_t end_clip_len = 0;
    if (aln->core.n_cigar > 0) {
        auto cigar_1elem = cigar[0];
        if (bam_cigar_op(cigar_1elem) == BAM_CSOFT_CLIP) {
            if (0 == strand) {
                inclu_beg_poss[0] += bam_cigar_oplen(cigar_1elem);
            } else {
                exclu_end_poss[1] += bam_cigar_oplen(cigar_1elem);
                end_clip_len = bam_cigar_oplen(cigar_1elem);
            }
        }
        cigar_1elem = cigar[aln->core.n_cigar-1];
        if (bam_cigar_op(cigar_1elem) == BAM_CSOFT_CLIP) {
            if (1 == strand) {
                inclu_beg_poss[1] -= bam_cigar_oplen(cigar_1elem);
            } else {
                exclu_end_poss[0] -= bam_cigar_oplen(cigar_1elem);
                end_clip_len = bam_cigar_oplen(cigar_1elem);
            }
        }
    }
    
    const uvc1_refgpos_t pos_incs[2] = {1, -1};
    {
        uint8_t prev_b = 0;
        uvc1_unsigned_int_t distinct_cnt = 0;
        int termpos = exclu_end_poss[strand] - pos_incs[strand];
        for (; termpos != inclu_beg_poss[strand] - pos_incs[strand]; termpos -= pos_incs[strand]) {
            uint8_t b = bam_seqi(bam_get_seq(aln), termpos);
            auto q = bam_get_qual(aln)[termpos];
            if (b != prev_b && q >= 20) {
                prev_b = b;
                distinct_cnt += 1;
                if (2 == distinct_cnt) { break; }
            }
        }
        uvc1_readpos_t homopol_tracklen = abs(termpos - (exclu_end_poss[strand] - pos_incs[strand]));
        uvc1_qual_t tail_penal = (end_clip_len >= 20 ? 1 : 0)
                + (homopol_tracklen >= 15 ? 2 : (homopol_tracklen >= 10 ? 1 : 0));
        if (tail_penal > 0) {
            const bool is_in_log_reg = (aln->core.tid == 0 && aln->core.pos < 9509431 && aln->core.pos > 9509400);
            if (is_in_log_reg) {
                LOG(logINFO) << "tail_penal = " << tail_penal << " for read " << bam_get_qname(aln);
            }
            for (uvc1_refgpos_t pos = exclu_end_poss[strand] - pos_incs[strand]; pos != (inclu_beg_poss[strand] - pos_incs[strand]) && pos != termpos; pos -= pos_incs[strand]) {
                const uvc1_qual_t q = bam_get_qual(aln)[pos];
                bam_get_qual(aln)[pos] = MAX(bam_get_qual(aln)[pos], tail_penal + 1) - tail_penal;
                if (is_in_log_reg) {
                    LOG(logINFO) << "\tQuality adjustment at pos " << pos << " : " << q << " -> " << (int)bam_get_qual(aln)[pos]; 
                }
            }
        }
    }
    {
        uvc1_refgpos_t homopol_len = 0;
        uint8_t prev_b = 0;
        // https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-451
        for (uvc1_refgpos_t pos = inclu_beg_poss[strand]; pos != exclu_end_poss[strand]; pos += pos_incs[strand]) {
            const uint8_t b = bam_seqi(bam_get_seq(aln), pos);
            if (b == prev_b) {
                homopol_len++;
                if (homopol_len >= 4 && b == seq_nt16_table['G']) {
                    bam_get_qual(aln)[pos] = MAX(bam_get_qual(aln)[pos], 1 + 1) - 1;
                }
            } else {
                prev_b = b;
                homopol_len = 1;
            }
        }
    }
    return 0;
}

int 
fill_strand_umi_readset_with_strand_to_umi_to_reads(
        std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, MolecularBarcode>> &umi_strand_readset,
        std::map<uvc1_hash_t, std::pair<std::array<std::map<uvc1_hash_t, std::vector<bam1_t *>>, 2>, MolecularBarcode>> &umi_to_strand_to_reads,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    for (auto & umi_to_strand_to_reads_element : umi_to_strand_to_reads) {
        const auto strand_to_reads = umi_to_strand_to_reads_element.second.first;
        const auto dflag = umi_to_strand_to_reads_element.second.second;
        umi_strand_readset.push_back(std::make_pair(std::array<std::vector<std::vector<bam1_t *>>, 2>(), dflag));
        for (int strand = 0; strand < 2; strand++) {
            for (auto read : strand_to_reads[strand]) {
                const std::vector<bam1_t *> alns = read.second;
                umi_strand_readset.back().first[strand].push_back(std::vector<bam1_t *>());
                for (auto aln : alns) {
                    apply_bq_err_correction3(aln, paramset.assay_sequencing_BQ_max, paramset.assay_sequencing_BQ_inc);
                    umi_strand_readset.back().first[strand].back().push_back(aln);
                }
            }
        }
    }
    return 0;
};

template <bool is_rc>
uvc1_hash_t 
bam2umihash(int & is_umi_found, const bam1_t *aln, const std::vector<uint8_t> & UMI_STRUCT, const int max_begin_diff_umi2read = 5) {
    LOG(logDEBUGx1) << "Going over " << UMI_STRUCT.size() << " bases in the pattern";

    auto *bamseq = bam_get_seq(aln);
    
    for (int i = 0; i < max_begin_diff_umi2read; i++) {
        size_t patpos = 0;
        uvc1_hash_t umihash = 0;
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

std::array<uvc1_readnum_t, 3>
bamfname_to_strand_to_familyuid_to_reads(
        std::map<uvc1_hash_t, std::pair<std::array<std::map<uvc1_hash_t, std::vector<bam1_t *>>, 2>, MolecularBarcode>> &umi_to_strand_to_reads,
        uvc1_refgpos_t & extended_inclu_beg_pos, 
        uvc1_refgpos_t & extended_exclu_end_pos,
        uvc1_refgpos_t tid, 
        uvc1_refgpos_t fetch_tbeg, 
        uvc1_refgpos_t fetch_tend, 
        bool end2end, 
        size_t regionbatch_ordinal, 
        size_t regionbatch_tot_num,
        const std::string UMI_STRUCT_STRING, 
        // const hts_idx_t * hts_idx,
        size_t thread_id,
        const std::vector<const bam1_t*> & bam_list,
        const CommandLineArgs & paramset,
        const uvc1_flag_t specialflag IGNORE_UNUSED_PARAM) {
    assert (fetch_tend > fetch_tbeg);
    
    const bool is_pair_end_merge_enabled = (PAIR_END_MERGE_NO != paramset.pair_end_merge);

    const bool should_log = (ispowerof2(regionbatch_ordinal+1) || ispowerof2(regionbatch_tot_num - regionbatch_ordinal));
    std::vector<uint8_t> umi_struct_string16;
    for (auto ch : UMI_STRUCT_STRING) {
        umi_struct_string16.push_back(seq_nt16_table[(int8_t)ch]);
    }
    for (auto base : umi_struct_string16) {
        LOG(logDEBUGx1) << "Base " << (int)base;
    }
    extended_inclu_beg_pos = INT32_MAX;
    extended_exclu_end_pos = 0;
    
    uvc1_readnum_t pcrpassed, umi_pcrpassed;
    pcrpassed = umi_pcrpassed = 0;
   
    //samFile *sam_infile = sam_open(paramset.bam_input_fname.c_str(), "r");
    if (should_log) {
        LOG(logINFO) << "Thread " << thread_id << " started dedupping the chunk tid" << tid << ":" << fetch_tbeg << "-" << fetch_tend 
                << " (region no " << regionbatch_ordinal << "/" << regionbatch_tot_num << " in this batch)";
    }
    uvc1_refgpos_t fetch_size = fetch_tend - fetch_tbeg + (ARRPOS_MARGIN + ARRPOS_OUTER_RANGE) * 2;
    
    std::vector<uvc1_readnum_t> inicount(fetch_size, 0);
    std::array<std::vector<uvc1_readnum_t>, 4> isrc_isr2_to_beg_count = {{ inicount, inicount, inicount, inicount }};
    std::array<std::vector<uvc1_readnum_t>, 4> isrc_isr2_to_end_count = {{ inicount, inicount, inicount, inicount }};
    std::vector<uvc1_readnum_big_t> inicount64(fetch_size + 1, 0);
    std::array<std::vector<uvc1_readnum_big_t>, 4> isrc_isr2_to_border_count_prefixsum = {{ inicount64, inicount64, inicount64, inicount64 }};;
    
    //hts_itr_t * hts_itr;
    // bam1_t *aln = bam_init1();
    
    std::array<uvc1_readnum_t, NUM_FILTER_REASONS> fillcode_to_num_alns;
    uvc1_readnum_t num_pass_alns = 0;
    uvc1_readnum_t num_iter_alns = 0;
    // hts_itr = sam_itr_queryi(hts_idx, tid, fetch_tbeg, fetch_tend);
    // while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) 
    for (const bam1_t *aln : bam_list){
        bool isrc = false;
        bool isr2 = false;
        uvc1_refgpos_t tBeg = 0;
        uvc1_refgpos_t tEnd = 0;
        uvc1_unsigned_int_t num_seqs = 0;
        FilterReason filterReason = fill_isrc_isr2_beg_end_with_aln(isrc, isr2, tBeg, tEnd, num_seqs, 
                aln, fetch_tbeg, fetch_tend, 
                paramset.kept_aln_min_aln_len, 
                paramset.kept_aln_min_mapqual,
                paramset.kept_aln_min_isize,
                paramset.kept_aln_max_isize,
                paramset.kept_aln_is_zero_isize_discarded,
                end2end, is_pair_end_merge_enabled);
        if (!is_pair_end_merge_enabled) { assert(!isr2); }
        if (NOT_FILTERED == filterReason) {
            isrc_isr2_to_beg_count[isrc * 2 + isr2][tBeg + ARRPOS_MARGIN - fetch_tbeg] += 1;
            isrc_isr2_to_end_count[isrc * 2 + isr2][tEnd + ARRPOS_MARGIN - fetch_tbeg] += 1;
            num_pass_alns += 1;
            extended_inclu_beg_pos = MIN(extended_inclu_beg_pos, SIGN2UNSIGN(aln->core.pos));
            extended_exclu_end_pos = MAX(extended_exclu_end_pos, SIGN2UNSIGN(bam_endpos(aln)));
        }
        fillcode_to_num_alns[filterReason]++;
        num_iter_alns += 1;
    }
    // sam_itr_destroy(hts_itr);
    for (size_t isrc_isr2 = 0; isrc_isr2 < 4; isrc_isr2++) {
        uvc1_readnum_big_t beg_prefixsum = 0;
        uvc1_readnum_big_t end_prefixsum = 0;
        isrc_isr2_to_border_count_prefixsum[isrc_isr2][0] = (beg_prefixsum + end_prefixsum);
        for (size_t i = 0; i < isrc_isr2_to_border_count_prefixsum[isrc_isr2].size() - 1; i++) {
            beg_prefixsum += isrc_isr2_to_beg_count[isrc_isr2][i];
            end_prefixsum += isrc_isr2_to_end_count[isrc_isr2][i];
            isrc_isr2_to_border_count_prefixsum[isrc_isr2][i+1] = (beg_prefixsum + end_prefixsum);
        }
    }
    std::array<std::vector<uvc1_refgpos_t>, 4> isrc_isr2_to_beg2bcenter = {{ inicount, inicount, inicount, inicount }};
    for (size_t isrc_isr2 = 0; isrc_isr2 < 4; isrc_isr2++) {
        auto beg_to_count = isrc_isr2_to_beg_count[isrc_isr2];
        poscounter_to_pos2pcenter(isrc_isr2_to_beg2bcenter[isrc_isr2], beg_to_count, paramset.dedup_center_mult);
    }
    std::array<std::vector<uvc1_refgpos_t>, 4> isrc_isr2_to_end2ecenter = {{ inicount, inicount, inicount, inicount }};
    for (size_t isrc_isr2 = 0; isrc_isr2 < 4; isrc_isr2++) {
        auto end_to_count = isrc_isr2_to_end_count[isrc_isr2];
        poscounter_to_pos2pcenter(isrc_isr2_to_end2ecenter[isrc_isr2], end_to_count, paramset.dedup_center_mult);
    }
    
    uvc1_readnum_t beg_peak_max = 0;
    for (auto beg_count : isrc_isr2_to_beg_count) {
        for (auto countval : beg_count) {
            beg_peak_max = MAX(beg_peak_max, countval);
        }
    }
    
    size_t alnidx = 0;
    // hts_itr = sam_itr_queryi(hts_idx, tid, fetch_tbeg, fetch_tend);
    // while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) 
    for (const bam1_t * aln : bam_list) {
        bool isrc = false;
        bool isr2 = false;
        uvc1_refgpos_t tBeg = 0;
        uvc1_refgpos_t tEnd = 0;
        uvc1_unsigned_int_t num_seqs = 0;
        FilterReason filterReason = fill_isrc_isr2_beg_end_with_aln(isrc, isr2, tBeg, tEnd, num_seqs, 
                aln, fetch_tbeg, fetch_tend, 
                paramset.kept_aln_min_aln_len, 
                paramset.kept_aln_min_mapqual,
                paramset.kept_aln_min_isize,
                paramset.kept_aln_max_isize,
                paramset.kept_aln_is_zero_isize_discarded,
                end2end, is_pair_end_merge_enabled);
        if (!is_pair_end_merge_enabled) { assert(!isr2); }
        if (NOT_FILTERED != filterReason) {
            continue;
        }
        const char *qname = bam_get_qname(aln);
        const uvc1_hash_t qname_hash = strhash(qname);
        const size_t qname_len = strlen(qname);
        const char *umi_beg1 = strchr(qname,   '#');
        const char *umi_beg = ((NULL != umi_beg1) ? (umi_beg1 + 1) : (qname + qname_len));
        const char *umi_end1 = strchr(umi_beg, '#');
        const char *umi_end = ((NULL != umi_end1) ? (umi_end1    ) : (qname + qname_len)); 
       
        int is_umi_found = ((umi_beg + 1 < umi_end) && (MOLECULE_TAG_NONE != paramset.molecule_tag)); // UMI has at least one letter
        int is_duplex_found = 0;
        uvc1_hash_t umihash = 0;
        const size_t umi_len = umi_end - umi_beg;
        if (is_umi_found) {
            size_t umi_half = (umi_end - umi_beg - 1) / 2;
            if ((umi_len % 2 == 1 ) && ( '+' == umi_beg[umi_half]) && (!paramset.disable_duplex)) {
                uvc1_hash_t umihash_part1 = strnhash(umi_beg               , umi_half); // alpha
                uvc1_hash_t umihash_part2 = strnhash(umi_beg + umi_half + 1, umi_half); // beta
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
        size_t isrc_isr2 = isrc * 2 + isr2;
        uvc1_refgpos_t beg1 = tBeg + ARRPOS_MARGIN - fetch_tbeg;
        uvc1_refgpos_t end1 = tEnd + ARRPOS_MARGIN - fetch_tbeg;
        uvc1_refgpos_t beg2 = isrc_isr2_to_beg2bcenter[isrc_isr2][beg1];
        uvc1_refgpos_t end2 = isrc_isr2_to_end2ecenter[isrc_isr2][end1];
        uvc1_readnum_t beg2count = isrc_isr2_to_beg_count[isrc_isr2][beg2];
        uvc1_readnum_t end2count = isrc_isr2_to_end_count[isrc_isr2][end2];
        const auto insert2posL = MIN(beg2, end2);
        const auto insert2posR = MAX(beg2, end2);
        const auto insert_cov_totDP = isrc_isr2_to_border_count_prefixsum[isrc_isr2][insert2posR] - isrc_isr2_to_border_count_prefixsum[isrc_isr2][insert2posL];
        
        /*
        uvc1_readnum_t beg2surrcount = 0;
        for (auto i = -ARRPOS_OUTER_RANGE; i < ARRPOS_OUTER_RANGE + 1; i++) {
            if (i > ARRPOS_INNER_RANGE || i < -ARRPOS_INNER_RANGE) {
                assert(i+beg2 < fetch_size || !fprintf(stderr, "beg2 index %d + %d = %d is too big!", i, beg2, i+beg2));
                uvc1_readnum_t beg_count = isrc_isr2_to_beg_count.at(isrc_isr2).at(i + beg2);
                beg2surrcount = MAX(beg2surrcount, beg_count);
            }
        }
        uvc1_readnum_t end2surrcount = 0;
        for (auto i = -ARRPOS_OUTER_RANGE; i < ARRPOS_OUTER_RANGE + 1; i++) {
            if (i > ARRPOS_INNER_RANGE && i < -ARRPOS_INNER_RANGE) {
                assert(i + end2 < fetch_size || !fprintf(stderr, "end2 index %d + %d = %d is too big!", i, end2, i + end2));
                uvc1_readnum_t end_count = isrc_isr2_to_end_count.at(isrc_isr2).at(i + end2);
                end2surrcount = MAX(end2surrcount, end_count);
            }
        }
        */
        const uvc1_readnum_big_t tot_ins_cov_border_DP = 
                isrc_isr2_to_border_count_prefixsum[isrc_isr2][insert2posR] 
              - isrc_isr2_to_border_count_prefixsum[isrc_isr2][insert2posL];
        // in the denominator we have a) -2 to take out two positions at beg2 and end2 and b) +2 to add pseudocount, and these two +-2 cancel out each other.
        double begratio = (double)(beg2count * (insert2posR - insert2posL) + 1) / (double)((tot_ins_cov_border_DP - beg2count) + (insert2posR - insert2posL));
        double endratio = (double)(end2count * (insert2posR - insert2posL) + 1) / (double)((tot_ins_cov_border_DP - end2count) + (insert2posR - insert2posL));
        const bool is_beg_amplicon = (begratio > paramset.dedup_amplicon_border_to_insert_cov_weak_avgDP_ratio);
        const bool is_end_amplicon = (endratio > paramset.dedup_amplicon_border_to_insert_cov_weak_avgDP_ratio);
        const bool is_beg_strong_amplicon = (begratio > paramset.dedup_amplicon_border_to_insert_cov_strong_avgDP_ratio);
        const bool is_end_strong_amplicon = (endratio > paramset.dedup_amplicon_border_to_insert_cov_strong_avgDP_ratio);

    /*
        double begfrac = (double)(beg2count + 1) / (double)(beg2surrcount + 2);
        double endfrac = (double)(end2count + 1) / (double)(end2surrcount + 2);
        const bool is_beg_amplicon = (begfrac > paramset.dedup_amplicon_count_to_surrcount_ratio_twosided);
        const bool is_end_amplicon = (endfrac > paramset.dedup_amplicon_count_to_surrcount_ratio_twosided);
        const bool is_beg_strong_amplicon = (begfrac > paramset.dedup_amplicon_count_to_surrcount_ratio);
        const bool is_end_strong_amplicon = (endfrac > paramset.dedup_amplicon_count_to_surrcount_ratio);
        
                double is_insert_amplicon_1 = (MIN(begratio, endratio) > paramset.dedup_amplicon_border_to_insert_cov_avgDP_ratio_of_min);
        double is_insert_amplicon_2 = (MAX(begratio, endratio) > paramset.dedup_amplicon_border_to_insert_cov_avgDP_ratio_of_max);

        */
        const bool is_assay_amplicon = (is_beg_strong_amplicon || is_end_strong_amplicon
                || (is_beg_amplicon && is_end_amplicon));
        pcrpassed += is_assay_amplicon;
        
        // beg end qname UMI = 1 2 4 8
        // IonTorrent amplicon without UMI: beg + end + qname
        // IonTorrent capture  without UMI: beg + end
        // IonTorrent amplicon with    UMI: beg + UMI
        // IonTorrent capture  with    UMI: beg + UMI
        // Illumina   amplicon without UMI: beg + end + qname
        // Illumina   capture  without UMI: beg + end
        // Illumina   amplicon with    UMI: beg + end + UMI
        // Illumina   capture  with    UMI: beg + end + UMI
        //
        // For Illumina amplicon with UMI:
        //     if beg * frac > end, then: beg + UMI
        //     if end * frac > beg, then: end + UMI
        
        uvc1_flag_t dedup_idflag = 0x0;
        if (paramset.dedup_flag != 0) {
            dedup_idflag = paramset.dedup_flag;
        } else if ((SEQUENCING_PLATFORM_IONTORRENT == paramset.inferred_sequencing_platform)) { // is_proton
            if (is_umi_found) { dedup_idflag = 0x9; }
            else if ((is_beg_amplicon && is_end_amplicon) || is_beg_strong_amplicon) {
                dedup_idflag = 0x7;
            } else {
                dedup_idflag = 0x3;
            }
        } else {
            if (is_umi_found) {
                if (is_beg_strong_amplicon && is_end_amplicon 
                        && beg2count > end2count * paramset.dedup_amplicon_end2end_ratio) {
                    dedup_idflag = 0x9;
                } else if (is_end_strong_amplicon && is_beg_amplicon 
                        && end2count > beg2count * paramset.dedup_amplicon_end2end_ratio) {
                    dedup_idflag = 0xA;
                } else {
                    dedup_idflag = 0xB;
                }
            } else if (is_beg_amplicon && is_end_amplicon) {
                dedup_idflag = 0x7;
            } else {
                dedup_idflag = 0x3;
            }
        }
        
        uvc1_hash_t molecule_hash = 0;
        if (0x3 == (0x3 & dedup_idflag)) {
            auto min2 = MIN(beg2, end2);
            auto max2 = MAX(beg2, end2);
            molecule_hash = hash2hash(molecule_hash, hash2hash(min2, max2)); 
        } else if (0x1 & dedup_idflag) {
            molecule_hash = hash2hash(molecule_hash, beg2+1); 
        } else if (0x2 & dedup_idflag) {
            molecule_hash = hash2hash(molecule_hash, end2+1); 
        }
        if (0x4 & dedup_idflag) {
            molecule_hash = hash2hash(molecule_hash, qname_hash);
        }
        if (0x8 & dedup_idflag) {
            molecule_hash = hash2hash(molecule_hash, umihash); 
        }
        
        int strand = (isrc ^ isr2);
        MolecularBarcode mb;
        mb.umistring = (is_umi_found ? std::string(umi_beg, umi_len) : "");
        mb.duplexflag = (is_umi_found ? 0x1 : 0) + (is_duplex_found ? 0x2 : 0) + (is_assay_amplicon ? 0x4 : 0);
        
        umi_to_strand_to_reads.insert(std::make_pair(molecule_hash, std::make_pair(std::array<std::map<uvc1_hash_t, std::vector<bam1_t *>>, 2>(), mb)));
        umi_to_strand_to_reads[molecule_hash].first[strand].insert(std::make_pair(qname_hash, std::vector<bam1_t *>()));
        umi_to_strand_to_reads[molecule_hash].first[strand][qname_hash].push_back(bam_dup1(aln));
        
        const bool should_log_read = (ispowerof2(alnidx + 1) || ispowerof2(num_pass_alns - alnidx));
        if (!is_pair_end_merge_enabled) { assert(!isr2); }
        if ((should_log_read && (beg_peak_max >= 2000 || should_log)) || paramset.always_log) {
            LOG(logINFO) << "thread_id = " << thread_id << " ; "
                    << "readname = " << qname << " ; "
                    << "alnidx = " << alnidx << " ; "
                    << "num_pass_alns = " << num_pass_alns << " ; "
                    << "isrc = " << isrc << " ; "
                    << "isr2 = " << isr2 << " ; "
                    << "strand = " << strand << " ; "
                    << "num_seqs = " << num_seqs << " ; "
                    << "dedup_idflag = " << dedup_idflag << " ; "
                    << "is_assay_amplicon = " << is_assay_amplicon << " ; "
                    << "tid = " << aln->core.tid << " ; "
                    << "fastaseq_range = " << tBeg << "," << tEnd << " ; "
                    << "original_range = " << beg1 << "," << end1 << " ; "
                    << "adjusted_rdiff = " << (beg2 - beg1) << "," << (end2 - end1) << " ; "
                    << "adjusted_count = " << beg2count << "," << end2count << " ; " 
                    //<< "adjusted_surrounding_counts = " << beg2surrcount << "," << end2surrcount << " ; " 
                    << "insert_cov_totDP =  " << insert_cov_totDP << " from-" << beg2 << "-to-" << end2 << " ; " 
                    << "barcode_umihash = " << (is_umi_found ? umihash : 0) << " ; "
                    << "molecule_hash = " << molecule_hash << " ; "
                    << "qname_hash = " << qname_hash << " ; "
                    << "dflag = " << mb.duplexflag << " ; "
                    << "UMIstring = " << umi_beg << " ; "
                    << "UMIsize = " << umi_len << " ; "
                    << "num_qname_from_molecule_so_far = " << umi_to_strand_to_reads[molecule_hash].first[strand].size() << " ; ";
        }
        alnidx += 1;
    }
    //sam_itr_destroy(hts_itr);
    
    //bam_destroy1(aln);
    //sam_close(sam_infile);
    if (should_log) { LOG(logINFO) << "Thread " << thread_id << " finished dedupping."; }
    return std::array<uvc1_readnum_t, 3>({num_pass_alns, pcrpassed, umi_pcrpassed});
}

