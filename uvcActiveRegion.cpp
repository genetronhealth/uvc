// This code is used for generating a list of positions that may be of variant candidates in BED format. 
// This code is not actually used for calling somatic variants because the speedup gained from using this bed file is too low.
// However, this code can still be useful in some scenarios (such as preparing a bed file for calling high-confidence variants in low-coverage bam file).
// Thus, this source file is kept.

#include <array>
#include <string>
#include <vector>

#include "assert.h"
#include "string.h"

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
updateDepthsByAln(std::vector<unsigned int> & region_var_dp, std::vector<unsigned int> & region_ref_dp,
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
                unsigned int base4bit = bam_seqi(bseq, qpos);
                if (region_str[rpos-region_offset] != base4bit && region_str[rpos-region_offset] && bam_get_qual(b)[qpos] > 20) {
                   region_var_dp[rpos-region_offset] += 1;
                } else {
                   region_ref_dp[rpos-region_offset] += 1;
                }
                qpos += 1;
                rpos += 1;
            }
        } else if (cigar_op == BAM_CINS) {
            region_var_dp[rpos-region_offset-1] += 1;
            region_var_dp[rpos-region_offset-0] += 1;
            qpos += cigar_oplen;
        } else if (cigar_op == BAM_CDEL) {
            region_var_dp[rpos-region_offset-1] += 1;
            region_var_dp[rpos-region_offset+cigar_oplen] += 1; 
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
    faidx_t *fai = fai_load(argv[1]);
    
    samFile *sam_infile = sam_open(argv[2], "r");
    bam_hdr_t * sam_header = sam_hdr_read(sam_infile);
    hts_idx_t * hts_idx = sam_index_load2(sam_infile, argv[2], NULL);
    
    for (unsigned int tid = 0; tid < sam_header->n_targets; tid++) {
        unsigned int seqlen = sam_header->target_len[tid];
        std::vector<unsigned int> region_var_dp(seqlen, 0);
        std::vector<unsigned int> region_ref_dp(seqlen, 0);
        
        const char *tname = faidx_iseq(fai, tid);
        assert(0 == strcmp(tname, sam_header->target_name[tid]));
        int regionlen;
        char *fetchedseq = faidx_fetch_seq(fai, tname, 0, seqlen - 1, &regionlen);
        assert(regionlen == seqlen);
        for (unsigned int i = 0; i < regionlen; i++) {
            fetchedseq[i] = CHAR_TO_SYMBOL.data[fetchedseq[i]]; 
        }
        // std::string region_string = std::string(fetchedseq);
        
        hts_itr_t * hts_itr = sam_itr_queryi(hts_idx, tid, 0, seqlen);
        bam1_t *aln = bam_init1();
        while (sam_itr_next(sam_infile, hts_itr, aln) >= 0) {
            updateDepthsByAln(region_var_dp, region_ref_dp,
                    aln, 0, fetchedseq, regionlen);
        }
        bam_destroy1(aln);
        hts_itr_destroy(hts_itr);
        free(fetchedseq);
        
        gen_bed_region(tname, region_var_dp, region_ref_dp);
    }
    
    hts_idx_destroy(hts_idx);
    bam_hdr_destroy(sam_header);
    sam_close(sam_infile);
}


