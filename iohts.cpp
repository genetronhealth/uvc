#include "iohts.hpp"
#include "common.hpp"

#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

#include <vector>

bool BedLine::is_valid() {
    return ((tid >= 0) || (tname.size() > 0)) && (beg_pos < end_pos);
}

std::vector<bam1_t *>
load_bam_records(
        const char * bam_fname,
        const hts_idx_t * hts_idx, 
        const uvc1_refgpos_t query_tid, 
        const uvc1_refgpos_t query_beg, 
        const uvc1_refgpos_t query_end) {
    
    std::vector<bam1_t *> ret;
    samFile *sam_infile = sam_open(bam_fname, "r");
    bam1_t *aln = bam_init1();
    hts_itr_t *hts_itr = sam_itr_queryi(hts_idx, query_tid, query_beg, query_end);
    int itr_result = 0;
    while ((itr_result = sam_itr_next(sam_infile, hts_itr, aln)) >= 0) { 
        ret.push_back(bam_dup1(aln));
    }
    sam_itr_destroy(hts_itr);
    bam_destroy1(aln);
    return ret;
}

