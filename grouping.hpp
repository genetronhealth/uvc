#ifndef grouping_hpp_INCLUDED
#define grouping_hpp_INCLUDED

#include "CmdLineArgs.hpp"
#include "common.hpp"

#include "htslib/sam.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define logDEBUGx1 logDEBUG // logINFO

typedef std::tuple<uvc1_refgpos_t, uvc1_refgpos_t, uvc1_refgpos_t, bool, uvc1_readnum_t> bedline_t;

int 
bed_fname_to_contigs(
        std::vector<bedline_t> & tid_beg_end_e2e_vec,
        const std::string & bed_fname, const bam_hdr_t *bam_hdr);

struct SamIter {
    const std::string input_bam_fname;
    const std::string & tier1_target_region; 
    const std::string region_bed_fname;
    const size_t nthreads;
    samFile *sam_infile = NULL;
    bam_hdr_t *samheader = NULL;
    hts_idx_t *sam_idx = NULL; 
    hts_itr_t *sam_itr = NULL;
    
    uvc1_refgpos_t endingpos = INT32_MAX;
    uvc1_refgpos_t tid = -1;
    uvc1_refgpos_t tbeg = INT32_MAX;
    uvc1_refgpos_t tend = INT32_MAX;
    uvc1_refgpos_t prev_tbeg = 0;
    uvc1_refgpos_t prev_tend = 0;
    uvc1_readnum_big_t nreads = 0;
    uvc1_readnum_big_t next_nreads = 0;
    bam1_t *alnrecord = bam_init1();
    
    std::vector<bedline_t> _tid_beg_end_e2e_vec;
    size_t _bedregion_idx = 0;
    
    SamIter(const std::string & in_bam_fname, 
            const std::string & tier1_target_reg, 
            const std::string & reg_bed_fname, 
            const uvc1_unsigned_int_t nt): 
            input_bam_fname(in_bam_fname), 
            tier1_target_region(tier1_target_reg), 
            region_bed_fname(reg_bed_fname), 
            nthreads(nt) {
        this->sam_infile = sam_open(input_bam_fname.c_str(), "r");
        if (NULL == this->sam_infile) {
            fprintf(stderr, "Failed to open the file %s!", input_bam_fname.c_str());
            abort();
        }
        this->samheader = sam_hdr_read(sam_infile);
        if (NULL == this->samheader) {
            fprintf(stderr, "Failed to read the header of the file %s!", input_bam_fname.c_str());
            abort();
        }
        if (NOT_PROVIDED != this->tier1_target_region) {
            this->sam_idx = sam_index_load(this->sam_infile, input_bam_fname.c_str());
            if (NULL == this->sam_idx) {
                fprintf(stderr, "Failed to load the index for the file %s!", input_bam_fname.c_str());
                abort();
            }
            this->sam_itr = sam_itr_querys(this->sam_idx, this->samheader, this->tier1_target_region.c_str());
            if (NULL == this->sam_itr) {
                fprintf(stderr, "Failed to load the region %s in the indexed file %s!", tier1_target_region.c_str(), input_bam_fname.c_str());
                abort();
            }
        }
        
        if (NOT_PROVIDED != this->region_bed_fname) {
            bed_fname_to_contigs(this->_tid_beg_end_e2e_vec, this->region_bed_fname, this->samheader); 
        }
    }
    ~SamIter() {
        bam_destroy1(alnrecord);
        if (NULL != sam_itr) { sam_itr_destroy(sam_itr); }
        if (NULL != sam_idx) { hts_idx_destroy(sam_idx); }
        bam_hdr_destroy(samheader);
        sam_close(sam_infile);
    }
    
    int 
    iternext(std::vector<bedline_t> & tid_beg_end_e2e_vec);
};

int
samfname_to_tid_to_tname_tseq_tup_vec(
        std::vector<std::tuple<std::string, uvc1_refgpos_t>> & tid_to_tname_tseqlen_tuple_vec, 
        const std::string & bam_input_fname);

int 
clean_fill_strand_umi_readset(
        std::vector<std::array<std::vector<std::vector<bam1_t *>>, 2>> &umi_strand_readset);

int 
fill_strand_umi_readset_with_strand_to_umi_to_reads(
        std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, uvc1_flag_t>> &umi_strand_readset,
        std::map<uvc1_hash_t, std::pair<std::array<std::map<uvc1_hash_t, std::vector<bam1_t *>>, 2>, uvc1_flag_t>> &umi_to_strand_to_reads,
        const CommandLineArgs & paramset,
        uvc1_flag_t specialflag);

std::array<uvc1_readnum_t, 3>
bamfname_to_strand_to_familyuid_to_reads(
        std::map<uvc1_hash_t, std::pair<std::array<std::map<uvc1_hash_t, std::vector<bam1_t *>>, 2>, uvc1_flag_t>> &umi_to_strand_to_reads,
        uvc1_refgpos_t & extended_inclu_beg_pos,
        uvc1_refgpos_t & extended_exclu_end_pos,
        uvc1_refgpos_t tid, 
        uvc1_refgpos_t fetch_tbeg, 
        uvc1_refgpos_t fetch_tend, 
        bool end2end, 
        size_t regionbatch_ordinal, 
        size_t regionbatch_tot_num,
        const std::string UMI_STRUCT_STRING, 
        const hts_idx_t * hts_idx,
        size_t thread_id,
        const CommandLineArgs & paramset,
        uvc1_flag_t specialflag);
#endif

