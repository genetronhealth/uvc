// This is a simple SNV caller designed to estimate the level of contamination from the tumor into the matched normal

#include <assert.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "htslib/sam.h"
//#include "htslib/faidx.h"
//#include "htslib/kstring.h"
//#include "htslib/khash.h"

#define MAX_ISIZE 1000
#define POS_BASE_ARRSIZE (MAX_ISIZE*100)

#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))

#define bam_phredi(b, i) (bam_get_qual((b))[(i)])

typedef struct {
    uint32_t data[4];
} base_distr_t;

const char SEQ_ACGT[4] = {'A', 'C', 'G', 'T'}; 

int 
vcf_hdr_print(int argc, char **argv) {
    time_t rawtime;
    time(&rawtime);
    char timestring[80];
    strftime(timestring, 80, "%F %T", localtime(&rawtime));
    printf("##fileformat=VCFv4.2\n");
    printf("##fileDate=%s\n", timestring);
    printf("##variantCallerCommand=");
    for (int i = 0; i < argc; i++) {
        printf(" %s", argv[i]);
    }
    printf("\n");
    // printf("##variantCallerVersion=%s\n");
    printf("##INFO=<ID=TAD,Number=1,Type=Integer,Description=\"Tumor allele depth\">\n");
    printf("##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Tumor total depth\">\n");
    printf("##INFO=<ID=NAD,Number=1,Type=Integer,Description=\"Normal allele depth\">\n");
    printf("##INFO=<ID=NDP,Number=1,Type=Integer,Description=\"Normal total depth\">\n");
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}

int
posidx_to_basedistr_update(base_distr_t *posidx_to_basedistr, const size_t offset, const bam1_t *b) {
    if (b->core.isize > MAX_ISIZE || b->core.qual < 50) { // insert size of at most MAX_ISIZE and MQ of at least 50
        return -1;
    }
    unsigned int qpos = 0;
    unsigned int rpos = b->core.pos;
    const uint32_t n_cigar = b->core.n_cigar;
    const uint32_t *cigar =  bam_get_cigar(b);
    const uint8_t *bseq = bam_get_seq(b);
    const uint8_t *bqual = bam_get_qual(b);
    unsigned int qsum = 0;
    for (int i = 0; i < b->core.l_qseq; i++) {
        qsum += bqual[i];
    }
    unsigned int qavg = qsum / b->core.l_qseq;
    for (unsigned int i = 0; i < n_cigar; i++) {
        uint32_t c = cigar[i];
        unsigned int cigar_op = bam_cigar_op(c);
        unsigned int cigar_oplen = bam_cigar_oplen(c);
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            for (int j = 0; j < cigar_oplen; j++) {
                unsigned int base4bit = bam_seqi(bseq, qpos);
                unsigned int base3bit = seq_nt16_int[base4bit];
                unsigned int qual8bit = bam_phredi(b, qpos);
                if (base3bit < 4 && qual8bit + 7 >= qavg) { // BQ of at least 30 ?
                    posidx_to_basedistr[rpos-offset].data[base3bit]++;
                }
                qpos += 1;
                rpos += 1;
            }
        } else if (cigar_op == BAM_CINS) {
            qpos += cigar_oplen;
        } else if (cigar_op == BAM_CDEL) {
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
            return -3;
        } else {
            return -4;
        }
    }
}

int
pos_base_distr_print(uint64_t *t_ad_sum_ptr, uint64_t *t_dp_sum_ptr, uint64_t *n_ad_sum_ptr, uint64_t *n_dp_sum_ptr,
        const char *tname, unsigned int beg, unsigned int end, const base_distr_t *t_base_distr_arr, const base_distr_t *n_base_distr_arr) {
    for (int i = 0; i < end-beg; i++) {
        unsigned int t_dp = 0;
        unsigned int n_dp = 0;
        for (int j = 0; j < 4; j++) {
            t_dp += t_base_distr_arr[i].data[j];
            n_dp += n_base_distr_arr[i].data[j];
        }
        if (t_dp < 30 || n_dp < 30) { continue; } // at least 30 depth of read coverage
        for (int j = 0; j < 4; j++) {
            unsigned int t_ad = t_base_distr_arr[i].data[j];
            unsigned int n_ad = n_base_distr_arr[i].data[j];
            if (t_ad < 5) { continue; } // at least 5 reads supporting the variant
            if (t_ad * 20 < t_dp) { continue; } // at least 5% allele fraction
            if (t_ad < n_ad * 3) { continue; } // at least 3 times more variant-supporting reads in the tumor than in the normal
            //if (t_ad * MAX(n_dp - n_ad, n_dp / 2) * 2 < 
            //    n_ad * MAX(t_dp - t_ad, t_dp / 2) * 5) { continue; } // at least 2.5 times more variant-supporting reads in the tumor than in the normal normalized by allele fraction
            if (t_ad * n_dp * 2 < n_ad * t_dp * 5) { continue; }
            printf("%s\t%d\t.\t.\t%c\tTAD=%d;TDP=%d;NAD=%d;NDP=%d\n", tname, beg+i+1, SEQ_ACGT[j], t_ad, t_dp, n_ad, n_dp); 
            *t_ad_sum_ptr += t_ad;
            *t_dp_sum_ptr += MAX(t_dp - t_ad, t_dp / 2);
            *n_ad_sum_ptr += n_ad;
            *n_dp_sum_ptr += MAX(n_dp - n_ad, n_dp / 2); 
        }
    }
}

int main(int argc,char** argv) {
    fprintf(stderr, "Usage is: %s tumor-bam normal-bam\n", argv[0]);
    if (argc < 3) return -1;
    
    samFile *t_in = sam_open(argv[1], "r");
    if (NULL == t_in) return -1;
    bam_hdr_t *t_hdr = sam_hdr_read(t_in);
    if (NULL == t_hdr) return -1;     
    bam1_t *t_b = bam_init1();
    base_distr_t t_base_distr_arr[POS_BASE_ARRSIZE+MAX_ISIZE];
    memset(t_base_distr_arr, 0, sizeof(t_base_distr_arr));
    int t_sam_read_ret = sam_read1(t_in, t_hdr, t_b);

    samFile *n_in = sam_open(argv[2], "r");
    if (NULL == n_in) return -1;
    bam_hdr_t *n_hdr = sam_hdr_read(n_in);
    if (NULL == n_hdr) return -1;
    bam1_t *n_b = bam_init1();
    base_distr_t n_base_distr_arr[POS_BASE_ARRSIZE+MAX_ISIZE];
    memset(n_base_distr_arr, 0, sizeof(n_base_distr_arr));
    int n_sam_read_ret = sam_read1(n_in, n_hdr, n_b); 
    
    uint64_t t_ad_sum = 0;
    uint64_t t_dp_sum = 0; 
    uint64_t n_ad_sum = 0;
    uint64_t n_dp_sum = 0;

    assert (t_hdr->n_targets == n_hdr->n_targets);
    size_t n_targets = t_hdr->n_targets;
    int tid = 0;
    vcf_hdr_print(argc, argv);
    for (int tid = 0; tid < n_targets; tid++) {
        assert (t_hdr->target_len[tid] == n_hdr->target_len[tid]);  
        size_t target_len = t_hdr->target_len[tid];  
        
        for (int pos_beg = 0; pos_beg < target_len; pos_beg += POS_BASE_ARRSIZE) {
            size_t pos_end = MIN(target_len, pos_beg + POS_BASE_ARRSIZE);
            if (pos_beg / POS_BASE_ARRSIZE % 100 == 0) {
                fprintf(stderr, "Processing %s:%d-%d\n", t_hdr->target_name[tid], pos_beg, pos_end);
            }
            int t_num_updates = 0;
            while (t_sam_read_ret >= 0 && t_b->core.tid <= tid && t_b->core.pos < pos_end) {
                posidx_to_basedistr_update(t_base_distr_arr, pos_beg, t_b);
                t_num_updates++;
                t_sam_read_ret = sam_read1(t_in, t_hdr, t_b);
            }
            int n_num_updates = 0;
            while (n_sam_read_ret >= 0 && n_b->core.tid <= tid && n_b->core.pos < pos_end) {
                posidx_to_basedistr_update(n_base_distr_arr, pos_beg, n_b);
                n_num_updates++;
                n_sam_read_ret = sam_read1(n_in, n_hdr, n_b);
            }
            if (t_num_updates && n_num_updates) {
                pos_base_distr_print(&t_ad_sum, &t_dp_sum, &n_ad_sum, &n_dp_sum,
                        t_hdr->target_name[tid], pos_beg, pos_end, t_base_distr_arr, n_base_distr_arr);
            }
            if (t_num_updates || n_num_updates) {
                for (int i = 0; i < MAX_ISIZE; i++) {
                    for (int j = 0; j < 4; j++) {
                        t_base_distr_arr[i].data[j] = t_base_distr_arr[POS_BASE_ARRSIZE+i].data[j];
                        n_base_distr_arr[i].data[j] = n_base_distr_arr[POS_BASE_ARRSIZE+i].data[j];
                    }
                }
                size_t t_tlen = t_hdr->target_len[t_b->core.tid];
                size_t n_tlen = n_hdr->target_len[n_b->core.tid];
                size_t restlen = MIN(MAX(t_tlen, n_tlen), POS_BASE_ARRSIZE);
                memset(&(t_base_distr_arr[MAX_ISIZE]), 0, sizeof(base_distr_t) * restlen);
                memset(&(n_base_distr_arr[MAX_ISIZE]), 0, sizeof(base_distr_t) * restlen);
            }
        }
    }
    
    double addfrac = (double)n_ad_sum / (double)(n_ad_sum + t_ad_sum + DBL_EPSILON);
    fprintf(stderr, "##(TADsum,TDPsum,NADsum,NDPsum)=(%d,%d,%d,%d)\n", t_ad_sum, t_dp_sum, n_ad_sum, n_dp_sum);
    fprintf(stderr, "##estimatedContaminationFromTumorToNormal_additiveFraction=%f (recommended=%f)\n", addfrac, floor(addfrac*175.0) / 100.0);
    fprintf(stderr, "##estimatedContaminationFromTumorToNormal_multiplicativeFraction=%f\n", (
            double)((double)n_ad_sum * (double)t_dp_sum) / (double)(((double)t_ad_sum * (double)n_dp_sum) + ((double)n_ad_sum * (double)t_dp_sum)));
    
    bam_destroy1(t_b);
    bam_hdr_destroy(t_hdr);
    sam_close(t_in);
    
    bam_destroy1(n_b);
    bam_hdr_destroy(n_hdr);
    sam_close(n_in);
    
    return 0;
}

