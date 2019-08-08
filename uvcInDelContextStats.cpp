// This code is used to analyze the context information of InDels, e.g. for ATCATCATC, RU=ATC and RC=3
// This code is actually not used in calling somatic variants but may be used to adjust the parameters to uvc

#include <array>
#include <map>
#include <string>
#include <vector>

#include <assert.h>
#include <limits.h>
#include <string.h>

#include "htslib/faidx.h"
#include "htslib/sam.h"

struct _CharToSymbol {
std::array<uint8_t, 128> data;
    _CharToSymbol() {
        for (int i = 0; i < 128; i++) {
            data[i] = 0;
        }
        data['A'] = data['a'] = 0x1;
        data['C'] = data['c'] = 0x2;
        data['G'] = data['g'] = 0x4;
        data['T'] = data['t'] = 0x8;
        data['I'] = data['i'] = 0x0;
        data['-'] = data['_'] = 0x0;
    }
};

const auto CHAR_TO_SYMBOL = _CharToSymbol();

int
updateInDelDepthsByAln(std::vector<unsigned int> & region_var_dp, std::vector<unsigned int> & region_ref_dp,
        const bam1_t *const b, const unsigned int region_offset, const char* region_str, unsigned int region_str_len) {
    if (b->core.flag & 0x4) { return - 1; }
    if (0 == (b->core.pos % (1024*64))) {
        fprintf(stderr, "rpos=%d, region_offset=%d, region_var_dp.size()=%d, region_str_len=%d\n", b->core.pos, region_offset, region_var_dp.size(), region_str_len);
    }
    assert(region_var_dp.size() == region_ref_dp.size());
    unsigned int qpos = 0;
    unsigned int rpos = b->core.pos;
    const uint32_t n_cigar = b->core.n_cigar;
    const uint32_t *cigar =  bam_get_cigar(b);
    const uint8_t *bseq = bam_get_seq(b);
    for (unsigned int i = 0; i < n_cigar; i++) {
        uint32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        unsigned int cigar_oplen = bam_cigar_oplen(c);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            for (unsigned int i2 = 0; i2 < cigar_oplen; i2++) {
                if (i2 > 0) {
                    region_ref_dp[rpos-region_offset] += 1;
                }
                qpos += 1;
                rpos += 1;
            }
        } else if (cigar_op == BAM_CINS) {
            region_var_dp[rpos-region_offset] += 1;
            qpos += cigar_oplen;
        } else if (cigar_op == BAM_CDEL) {
            region_var_dp[rpos-region_offset] += 1;
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

int indelpos_to_context(
        std::string & repeatunit, unsigned int & max_repeatnum,
        const std::string & refstring, unsigned int refpos) {
    max_repeatnum = 0;
    if (refpos >= refstring.size()) {
        repeatunit = "";
        return -1;
    }
    unsigned int repeatsize_at_max_repeatnum = 0;
    for (unsigned int repeatsize = 1; repeatsize < 6*2; repeatsize++) {
        unsigned int qidx = refpos;
        while (qidx + repeatsize < refstring.size() && refstring[qidx] == refstring[qidx+repeatsize]) {
            qidx++;
        }
        unsigned int repeatnum = (qidx - refpos) / repeatsize + 1;
        if (repeatnum > max_repeatnum) {
            max_repeatnum = repeatnum;
            repeatsize_at_max_repeatnum = repeatsize;
        }
    }
    repeatunit = refstring.substr(refpos, repeatsize_at_max_repeatnum);
}

struct IndelContext {
    std::string repeatunit;
    unsigned int repeatsize;
    int lendelta;
    IndelContext(std::string ru, unsigned int rs, int ld) {
        repeatunit = ru;
        repeatsize = rs;
        lendelta = ld;
    }
    bool operator<(const IndelContext & other) const {
        if (repeatsize != other.repeatsize) { return repeatsize < other.repeatsize; }
        if (lendelta != other.lendelta) { return lendelta < other.lendelta; }
        if (repeatunit != other.repeatunit) { return repeatunit < other.repeatunit; }
        return 0;
    }
};

int
map_inc(std::map<IndelContext, unsigned int> & indelmap,
        const std::string & repeatunit, const unsigned int repeatsize, const int lendelta) {
    IndelContext indelcontext(repeatunit, repeatsize, lendelta);
    indelmap.insert(std::make_pair(indelcontext, 0));
    indelmap[indelcontext]++;
}

int
map_print(const char *prefix, std::map<IndelContext, unsigned int> & indelmap) {
    for (auto context2count: indelmap) {
        auto context = context2count.first;
        auto count = context2count.second;
        printf("%s\t%10s\t%d\t%d\t%d\n", prefix, context.repeatunit.c_str(), context.repeatsize, context.lendelta, count);
    }
}

int
updateInDelStringsByAln(
        std::map<IndelContext, unsigned int> & errors_map,
        std::map<IndelContext, unsigned int> & homozy_map,
        std::map<IndelContext, unsigned int> & hetero_map,
        std::map<IndelContext, unsigned int> & lowdep_map,
        const std::vector<unsigned int> & region_var_dp, const std::vector<unsigned int> & region_ref_dp,
        const bam1_t *const b, const unsigned int region_offset, const std::string & region_string, unsigned int region_str_len) {
    if (b->core.flag & 0x4) { return - 1; }
    if (0 == (b->core.pos % (1024*64))) {
        fprintf(stderr, "rpos=%d, region_offset=%d, region_var_dp.size()=%d, region_str_len=%d\n", b->core.pos, region_offset, region_var_dp.size(), region_str_len);
    }
    assert(region_var_dp.size() == region_ref_dp.size());
    unsigned int qpos = 0;
    unsigned int rpos = b->core.pos;
    const uint32_t n_cigar = b->core.n_cigar;
    const uint32_t *cigar =  bam_get_cigar(b);
    const uint8_t *bseq = bam_get_seq(b);
    for (unsigned int i = 0; i < n_cigar; i++) {
        uint32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        unsigned int cigar_oplen = bam_cigar_oplen(c);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            qpos += cigar_oplen;
            rpos += cigar_oplen;
        } else if (cigar_op == BAM_CINS || cigar_op == BAM_CDEL) {
            unsigned int mapidx = (cigar_op == BAM_CINS ? 0 : 1);
            int lendelta = (cigar_op == BAM_CINS ? ((int)cigar_oplen) : (-int(cigar_oplen)));
            std::string repeatunit;
            unsigned int repeatsize = 0;
            indelpos_to_context(repeatunit, repeatsize, region_string, rpos-region_offset);
            unsigned int ref_dp = region_ref_dp.at(rpos-region_offset); 
            unsigned int var_dp = region_var_dp.at(rpos-region_offset); 
            
            if (ref_dp + var_dp >= 25) {
                if        (var_dp * 4 <= ref_dp) {
                    map_inc(errors_map, repeatunit, repeatsize, lendelta);
                } else if (ref_dp * 4 <= var_dp) {
                    map_inc(hetero_map, repeatunit, repeatsize, lendelta);
                } else {
                    map_inc(homozy_map, repeatunit, repeatsize, lendelta);
                }
            } else {
                map_inc(lowdep_map, repeatunit, repeatsize, lendelta);
            }
            if (cigar_op == BAM_CINS) {
                qpos += cigar_oplen;
            } else {
                rpos += cigar_oplen;
            }
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

int
gen_bed_region(const char *tname, std::vector<unsigned int> & region_var_dp, std::vector<unsigned int> & region_ref_dp) {
    assert(region_var_dp.size() == region_ref_dp.size());
    for (unsigned i = 0; i < region_var_dp.size(); i++) {
        if (region_var_dp[i] >= 4 && region_var_dp[i] * (200-1) > region_ref_dp[i]) {
            printf("%s\t%d\t%d\t%d/%d\n", tname, (i > 2 ? (i-2) : 0), (i+2)+1, region_var_dp[i], region_ref_dp[i]);
        }
    }
}

int main(int argc, char **argv) {
    std::map<IndelContext, unsigned int> errors_map;
    std::map<IndelContext, unsigned int> homozy_map;
    std::map<IndelContext, unsigned int> hetero_map;
    std::map<IndelContext, unsigned int> lowdep_map;

    faidx_t *fai = fai_load(argv[1]);
    
    samFile *sam_infile = sam_open(argv[2], "r");
    bam_hdr_t * sam_header = sam_hdr_read(sam_infile);
    hts_idx_t * hts_idx = sam_index_load2(sam_infile, argv[2], NULL);
    
    unsigned int tid_beg = 0;
    unsigned int tid_end = UINT_MAX;
    if (argc == 5) {
        tid_beg = atoi(argv[3]);
        tid_end = atoi(argv[4]);
    }
    if (tid_end > sam_header->n_targets) {
        tid_end = sam_header->n_targets;
    }
    for (unsigned int tid = tid_beg; tid < tid_end; tid++) {
        unsigned int seqlen = sam_header->target_len[tid];
        std::vector<unsigned int> region_var_dp(seqlen, 0);
        std::vector<unsigned int> region_ref_dp(seqlen, 0);
        
        const char *tname = faidx_iseq(fai, tid);
        assert(0 == strcmp(tname, sam_header->target_name[tid]));
        int regionlen;
        char *fetchedseq = faidx_fetch_seq(fai, tname, 0, seqlen - 1, &regionlen);
        assert(regionlen == seqlen);
        std::string region_string(fetchedseq, regionlen); 
        
        for (unsigned int i = 0; i < regionlen; i++) {
            fetchedseq[i] = CHAR_TO_SYMBOL.data[fetchedseq[i]]; 
        }
        // std::string region_string = std::string(fetchedseq);
        
        bam1_t *aln = bam_init1();
        hts_itr_t * hts_itr = sam_itr_queryi(hts_idx, tid, 0, seqlen);
        while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) {
            updateInDelDepthsByAln(region_var_dp, region_ref_dp,
                    aln, 0, fetchedseq, regionlen);
        }
        hts_itr_destroy(hts_itr);
        
                
        hts_itr = sam_itr_queryi(hts_idx, tid, 0, seqlen);
        while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) {
            updateInDelStringsByAln(
                    errors_map, homozy_map, hetero_map, lowdep_map,
                    region_var_dp, region_ref_dp,
                    aln, 0, region_string, regionlen);
        }
        bam_destroy1(aln);
        free(fetchedseq);
#if 0
        for (unsigned int k = 0; k < region_string.size(); k++) {
            printf("%c", region_string[k]);
            if (0 == (k % 128)) {
                printf("\n");
            }
        }
        printf("\n");
#endif
        // printf("%s\n", region_string.c_str());
        // gen_bed_region(tname, region_var_dp, region_ref_dp);
    }
    
    hts_idx_destroy(hts_idx);
    bam_hdr_destroy(sam_header);
    sam_close(sam_infile);
    
    map_print("errors", errors_map);
    map_print("homozy", homozy_map);
    map_print("hetero", hetero_map);
    map_print("lowdep", lowdep_map);
}


