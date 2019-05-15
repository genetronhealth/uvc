#include <htslib/vcf.h>
#include <float.h>
#include <math.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

const unsigned int LOD_THRES = 30;

void xfree(void *ptr) {
    if (ptr != NULL) { free(ptr); }
}

int main(int argc, char **argv) {
    
    htsFile *infile = hts_open(argv[1], "r");
    bcf_hdr_t *inhdr = bcf_hdr_read(infile);
    htsFile *outfile = hts_open(argv[2], "wz");
    bcf_hdr_write(outfile, inhdr);
    bcf1_t *bcfrec = bcf_init();
    
    int32_t *gEND = NULL;
    int32_t *GQ = NULL;
    float *FA = NULL; 
    float *VAQ = NULL;
    int32_t *gt_arr = NULL; 

    unsigned int nsmpl = bcf_hdr_nsamples(inhdr);
    assert (nsmpl == 2);
    int32_t prev_normal_gEND = 0;
    int32_t curr_normal_gEND = 0;
    int32_t curr_normal_rid = -1;
    int32_t curr_normal_has_germline_var = 0;
    int32_t curr_normal_GQ = 0;
    
    while (bcf_read(infile, inhdr, bcfrec) >= 0) {
        // int ndst = 0;
        bcf_unpack(bcfrec, BCF_UN_ALL);
        //fprintf(stderr, "Processing rid %d pos %d\n", bcfrec->rid, bcfrec->pos);
        if (0 == strcmp("<NON_REF>", bcfrec->d.allele[1])) {
            
            //fprintf(stderr, "Processing NON_REF rid %d pos %d\n", bcfrec->rid, bcfrec->pos);
            int mem_gEND = 0; 
            int32_t num_gEND = bcf_get_format_int32(inhdr, bcfrec, "gEND", &gEND, &mem_gEND);
            //fprintf(stderr, "Processing NON_REF num_gEND rid %d pos %d , nvals = %d\n", bcfrec->rid, bcfrec->pos, num_gEND);
            assert (num_gEND == nsmpl);
            int32_t normal_gEND = gEND[1];
            if (normal_gEND != bcf_int32_missing) {
                
                //fprintf(stderr, "Processing NON_REF GQ rid %d pos %d with normal_gEND = %d \n", bcfrec->rid, bcfrec->pos, normal_gEND); 
                int mem_GQ = 0; 
                int32_t num_GQ = bcf_get_format_int32(inhdr, bcfrec, "GQ", &GQ, &mem_GQ);
                curr_normal_GQ = GQ[1]; //100; // (int)GQ[1];
                
                //fprintf(stderr, "Processing NON_REF gEND rid %d pos %d\n", bcfrec->rid, bcfrec->pos);
                prev_normal_gEND = curr_normal_gEND;
                curr_normal_gEND = normal_gEND;
                curr_normal_rid = bcfrec->rid;
                
                int ngt_arr = 0;
                int num_GT = bcf_get_genotypes(inhdr, bcfrec, &gt_arr, &ngt_arr);
                assert (num_GT == nsmpl * 2);
                int max_ploidy = num_GT/nsmpl; // diploid 
                curr_normal_has_germline_var = 0;
                for (int j = 0; j < max_ploidy; j++) {
                    int allele = bcf_gt_allele(gt_arr[max_ploidy+j]);
                    curr_normal_has_germline_var += allele;
                }
                // assert (0);
            }
        } else {
            //fprintf(stderr, "Processing FA rid %d pos %d\n", bcfrec->rid, bcfrec->pos);
            int mem_FA = 0;
            int num_FA = bcf_get_format_float(inhdr, bcfrec, "FA", &FA, &mem_FA);
            assert (num_FA == nsmpl || !fprintf(stderr, "%d != %d", num_FA, nsmpl));
            float  tumor_FA = ((bcf_float_missing == FA[0] || isnan(FA[0])) ? 0 : FA[0]) + DBL_EPSILON; // FA[0];
            float normal_FA = ((bcf_float_missing == FA[1] || isnan(FA[1])) ? 0 : FA[1]) + DBL_EPSILON;
            assert( tumor_FA >= 0 || !fprintf(stderr, " tumor_FA %f >= 0 failed at rid %d and position %d!\n",  tumor_FA, bcfrec->rid, bcfrec->pos));
            assert(normal_FA >= 0 || !fprintf(stderr, "normal_FA %f >= 0 failed at rid %d and position %d!\n", normal_FA, bcfrec->rid, bcfrec->pos));
            //fprintf(stderr, "Processing VAQ rid %d pos %d\n", bcfrec->rid, bcfrec->pos);
            int mem_VAQ = 0;
            int num_VAQ = bcf_get_format_float(inhdr, bcfrec, "VAQ", &VAQ, &mem_VAQ);
            assert (num_VAQ == nsmpl);
            float  tumor_VAQ = ((bcf_float_missing == VAQ[0] || isnan(VAQ[0])) ? 0 : VAQ[0]) + DBL_EPSILON;
            float normal_VAQ = ((bcf_float_missing == VAQ[1] || isnan(VAQ[1])) ? 0 : VAQ[1]) + DBL_EPSILON;
            assert ( tumor_VAQ >= 0);
            assert (normal_VAQ >= 0);
            //fprintf(stderr, "Processing LOD rid %d pos %d\n", bcfrec->rid, bcfrec->pos);
            float tlod = (tumor_VAQ - normal_VAQ) * MAX(tumor_FA - normal_FA, 0) / (tumor_FA + normal_FA + DBL_EPSILON);
            float nlod = (((!curr_normal_has_germline_var) && bcfrec->pos <= curr_normal_gEND && bcfrec->rid == curr_normal_rid) ? curr_normal_GQ : 0) + 30;
            
            bcfrec->qual = MIN(tlod, nlod);
            if (bcfrec->qual >= LOD_THRES) {
                bcf_write(outfile, inhdr, bcfrec);
            }
        }
    }
    
    xfree(gEND);
    xfree(GQ);
    xfree(FA); 
    xfree(VAQ);
    xfree(gt_arr);
    
    bcf_destroy(bcfrec);
    bcf_hdr_destroy(inhdr);
    bcf_close(outfile);
    bcf_close(infile);
}

// p=subprocess.Popen(['bcftools', 'merge', '--force-samples', '-m', 'none', sys.argv[1], sys.argv[2]], stdout=subprocess.PIPE)
   /* 
    if record.REF != record.ALT[0]:
        normalBGrecord.POS = record.POS
        normalBGrecord.REF = record.REF[0:1]
        if prevNormalBGrecord != normalBGrecord: vcf_writer.write_record(normalBGrecord)
        vcf_writer.write_record(record)
     */   
