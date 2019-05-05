
/*
 * https://www.snip2code.com/Snippet/589504/answer-to-https---www-biostars-org-p-151
 * answer to https://www.biostars.org/p/151053/  htslib C API: getting reads from region. 
 * Compilation  gcc -Isamtools  -Ihtslib -Lsamtools -Lhtslib biostar151053.c htslib/libhts.a samtools/libbam.a  -lz -lpthread
 * by Pierre Lindenbaum  @ Pierre Lindenbaum
 * */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "samtools.h"

int main(int argc,char** argv)
{
    hts_itr_t *iter=NULL;
    hts_idx_t *idx=NULL;
    samFile *in = NULL;
    bam1_t *b= NULL;
    bam_hdr_t *header = NULL;
    if(argc!=3) return -1;
    in = sam_open(argv[1], "r");
    if(in==NULL) return -1;
    if ((header = sam_hdr_read(in)) == 0) return -1;
    idx = sam_index_load(in,  argv[1]);
    if(idx==NULL) return -1;
    iter  = sam_itr_querys(idx, header, argv[2]); 
    if(iter==NULL) return -1;
    b = bam_init1();
    while ( sam_itr_next(in, iter, b) >= 0) 
    {
        fputs("DO STUFF\n",stdout); 
    }
    hts_itr_destroy(iter);
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
    return 0;
}

