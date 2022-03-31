#ifndef IS_IOHTS_INCLUDED
#define IS_IOHTS_INCLUDED

#include "common.hpp"

#include "htslib/sam.h"
#include "htslib/hts.h"

#include <string>
#include <vector>

#define BED_END_TO_END_BIT 0x1

struct BedLine {
    std::string tname; // can be set to empty if tid != -1
    uvc1_refgpos_t tid; // can be set to -1 if !tname.isEmpty()
    uvc1_refgpos_t beg_pos;
    uvc1_refgpos_t end_pos;
    uvc1_flag_t region_flag;
    uvc1_readnum_big_t n_reads;
    
    BedLine(
            uvc1_refgpos_t a_tid,
            uvc1_refgpos_t a_beg_pos,
            uvc1_refgpos_t a_end_pos,
            const uvc1_flag_t a_region_flag,
            const uvc1_readnum_big_t a_n_reads) {
        this->tid = a_tid;
        this->beg_pos = a_beg_pos;
        this->end_pos = a_end_pos;
        this->region_flag = a_region_flag;
        this->n_reads = a_n_reads;
    }; 
    bool is_valid();
};

std::vector<bam1_t *>
load_bam_records(
        int & sam_itr_queryi_ret,
        samFile *samfile,
        const hts_idx_t * hts_idx,
        const uvc1_refgpos_t query_tid,
        const uvc1_refgpos_t query_beg, 
        const uvc1_refgpos_t query_end);

#endif

