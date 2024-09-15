
#include "main.hpp"

#include "CmdLineArgs.hpp"
#include "common.hpp"
#include "grouping.hpp"
#include "iohts.hpp"
#include "version.h"

#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/synced_bcf_reader.h"

#include <chrono>
#include <ctime>
#include <future>
#include <thread>
#include <tuple>

#if !defined(USE_STDLIB_THREAD)
#include "omp.h"
#endif

void 
xfree(void *ptr) {
    if (NULL != ptr) { free(ptr); }
}

template <bool TIsInputCompressed>
int 
clearstring(BGZF * bgzip_file, const std::string & outstring_allp, bool is_output_to_stdout = false, bool flush = !TIsInputCompressed) {
    if (is_output_to_stdout) {
        std::cout << outstring_allp;
        return outstring_allp.size();
    }
    if (NULL == bgzip_file) { return -1; }
    int ret = 0;
    if (TIsInputCompressed) {
        ret = bgzf_raw_write(bgzip_file, outstring_allp.c_str(), outstring_allp.size());
        LOG(logINFO) << "Written " << ret << " bytes of compressed data from " << outstring_allp.size()  << " bytes of compressed data.";
    } else {
        ret = bgzf_write(bgzip_file, outstring_allp.c_str(), outstring_allp.size());
        LOG(logINFO) << "Written " << ret << " bytes of compressed data from " << outstring_allp.size()  << " bytes of raw data.";
    }
    if (flush) { 
        int flushret = bgzf_flush(bgzip_file); 
        if (flushret != 0) {
            return flushret;
        }
    }
    return ret;
};

std::string 
load_refstring(const faidx_t *ref_faidx, uvc1_refgpos_t tid, uvc1_refgpos_t incbeg, uvc1_refgpos_t excend) {
    assertUVC(incbeg < excend);
    if (NULL == ref_faidx) {
        return std::string(excend - incbeg, 'n');
    }
    const char *tname = faidx_iseq(ref_faidx, tid);
    int regionlen;
    char *fetchedseq = faidx_fetch_seq(ref_faidx, tname, incbeg, excend - 1, &regionlen);
    assertUVC (regionlen == (excend - incbeg) || !fprintf(stderr, "%d == %u - %u failed for %s:%d-%d\n", regionlen, excend, incbeg, tname, incbeg, excend));
    std::string ret(fetchedseq);
    for (size_t i = 0; i < ret.size(); i++) {
        ret[i] = toupper(ret[i]);
    }
    free(fetchedseq);
    return ret;
};

template <class TKey, class TVal>
std::vector<std::pair<TKey, TVal>>
map2vector(const std::map<TKey, TVal> & key2val4map) {
    std::vector<std::pair<TKey, TVal>> ret;
    for (auto kv : key2val4map) {
        ret.push_back(kv);
    }
    return ret;
};

std::map<std::pair<uvc1_refgpos_t, AlignmentSymbol>, std::set<size_t>>
mutform2count4vec_to_simplemut2indices(
        std::vector<HapLink> mutform2count4vec) {
    
    std::map<std::pair<uvc1_refgpos_t, AlignmentSymbol>, std::set<size_t>> simplemut2indices;
    for (size_t i = 0; i < mutform2count4vec.size(); i++) {
        auto counts = mutform2count4vec[i].fr_cnts; 
        if (counts[0] + counts[1] < 2) { continue; }
        std::basic_string<std::pair<uvc1_refgpos_t, AlignmentSymbol>> mutset = mutform2count4vec[i].pos_symb_string;
        for (auto simplemut : mutset) {
            simplemut2indices.insert(std::make_pair(simplemut, std::set<size_t>()));
            simplemut2indices[simplemut].insert(i);
        }
    }
    return simplemut2indices;
};

int 
bgzip_string(std::string & compressed_outstring, const std::string & uncompressed_outstring) {
    char *compressed_outstr = (char*)malloc(uncompressed_outstring.size() * sizeof(char));
    if (NULL == compressed_outstr) {
        fprintf(stderr, "The library function malloc failed at line %d in file %s!\n", __LINE__, __FILE__);
        exit(-1);
    }
    size_t compressed_capacity = uncompressed_outstring.size();
    size_t compressed_totlen = 0;
    size_t uncompress_totlen = 0;
    do {
        if (compressed_totlen + BGZF_BLOCK_SIZE >= compressed_capacity) {
            char *compressed_outstr_tmp = (char*)realloc(compressed_outstr, (compressed_capacity * 2 + BGZF_BLOCK_SIZE) * sizeof(char));
            if (NULL == compressed_outstr_tmp) {
                fprintf(stderr, "The library function realloc failed at line %d in file %s!\n", __LINE__, __FILE__);
                exit(-2);
            }
            compressed_outstr = compressed_outstr_tmp;
            compressed_capacity = compressed_capacity * 2 + BGZF_BLOCK_SIZE;
        }
        size_t block_len = MIN(uncompressed_outstring.size() - uncompress_totlen, (size_t)BGZF_BLOCK_SIZE);
        if (0 == block_len) { break; }
        size_t compressed_len = 0;
        bgzf_compress(compressed_outstr + compressed_totlen, &compressed_len, uncompressed_outstring.c_str() + uncompress_totlen, block_len, 5);
        uncompress_totlen += block_len;
        compressed_totlen += compressed_len;
    } while (uncompress_totlen < uncompressed_outstring.size());
    
    compressed_outstring += std::string(compressed_outstr, compressed_totlen);
    free(compressed_outstr);
    return 0;
}

struct BatchArg {
    std::array<std::string, NUM_FQLIKE_CON_OUT_FILES> outstring3fastq; // outstring_fastq;
    std::string outstring_allp;
    std::string outstring_pass;
    int thread_id;
    samFile *samfile;
    hts_idx_t *hts_idx;
    faidx_t *ref_faidx;
    bcf_hdr_t *bcf_hdr;
    bcf_srs_t *sr;
    
    BedLine prev_bedline;
    BedLine bedline;
    std::tuple<std::string, uint32_t> tname_tseqlen_tuple;
    size_t regionbatch_ordinal;
    size_t regionbatch_tot_num;

    const CommandLineArgs paramset;
    const std::string UMI_STRUCT_STRING;
    const bool is_vcf_out_pass_to_stdout;
};

std::string 
als_to_string(const char *const* const allele, uint32_t n_allele) {
    std::string ret;
    ret.reserve(n_allele*2);
    for (uint32_t i = 0; i < n_allele; i++) {
        if (0 == i) {
            ret += allele[i];
        } else if (1 == i) {
            ret += std::string("\t") + allele[i];
        } else {
            ret += std::string(",") + allele[i];
        }
    }
    return ret;
}

#define assertAlways(x) assert((x));

#define BCF_GET_FORMAT_INT32_WITH_CHECK(v) \
        ndst_val = 0; \
        valsize = bcf_get_format_int32(bcf_hdr, line, v, &bcfints, &ndst_val); \
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for %s and line %ld!\n", ndst_val, valsize, v, line->pos));

#define BCF_GET_FORMAT_INT32B_WITH_CHECK(k,v) \
        ndst_val = 0; \
        valsize = bcf_get_format_int32(bcf_hdr, line, v, &bcfints, &ndst_val); \
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for %s and line %ld!\n", ndst_val, valsize, v, line->pos)); \
        (k) = (bcfints[0] + bcfints[1]);

template <class T1, class T2, class T3>
const std::map<std::tuple<uvc1_refgpos_t, uvc1_refgpos_t, AlignmentSymbol>, std::vector<TumorKeyInfo>>
rescue_variants_from_vcf(
        const T1 & bedlines,
        const T2 & tid_to_tname_tlen_tuple_vec,
        const std::string & vcf_tumor_fname, 
        const T3 *bcf_hdr,
        const bool is_tumor_format_retrieved) {
    std::map<std::tuple<uvc1_refgpos_t, uvc1_refgpos_t, AlignmentSymbol>, std::vector<TumorKeyInfo>> ret;
    if (ISNT_PROVIDED(vcf_tumor_fname)) {
        return ret;
    }
    std::string regionstring;
    std::string targetstring;
    int64_t prev_tid = -1;
    int64_t prev_beg = -1;
    int64_t prev_end = -1;
    for (const auto & bedline : bedlines) {
        if (prev_tid != bedline.tid) {
            if (prev_tid != -1) {
                if (0 < regionstring.size()) { regionstring += std::string(","); }
                regionstring += (std::get<0>(tid_to_tname_tlen_tuple_vec[prev_tid]) + ":" + std::to_string(prev_beg+1) + "-" + std::to_string(prev_end+1-1));
            }
            prev_tid = bedline.tid;
            prev_beg = bedline.beg_pos;
        }

        const uvc1_refgpos_t tid = bedline.tid; // std::get<0>(tid_beg_end_e2e);
        const uvc1_refgpos_t rpos_inclu_beg = bedline.beg_pos; // std::get<1>(tid_beg_end_e2e);
        const uvc1_refgpos_t rpos_exclu_end = bedline.end_pos; // std::get<2>(tid_beg_end_e2e);
        
        const auto & tname_tseqlen_tuple = tid_to_tname_tlen_tuple_vec[tid];
        if (0 < targetstring.size()) { targetstring += std::string(","); }
        std::string targetspec = (std::get<0>(tname_tseqlen_tuple) + ":" + std::to_string(rpos_inclu_beg+1) + "-" + std::to_string(rpos_exclu_end+1-1));
        targetstring += targetspec;
        
        prev_end = bedline.end_pos;        
    }
    
    if (prev_tid != -1) {
        if (0 < regionstring.size()) { regionstring += std::string(","); }
        regionstring += (std::get<0>(tid_to_tname_tlen_tuple_vec[prev_tid]) + ":" + std::to_string(prev_beg+1) + "-" + std::to_string(prev_end+1-1));
    }
    
    { 
        LOG(logINFO) << "BCF region is " << regionstring;
        LOG(logINFO) << "BCF target is " << targetstring;

    }
    if (targetstring.size() == 0) {
        return ret;
    }
    bcf_srs_t *const sr = bcf_sr_init();
    if (NULL == sr) {
        LOG(logCRITICAL) << "Failed to initialize bcf sr";
        exit(-6);
    }
    
    bcf_sr_set_regions(sr, regionstring.c_str(), false);
    bcf_sr_set_targets(sr, targetstring.c_str(), false, 0);

    int sr_set_opt_retval = bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    if (sr_set_opt_retval < 0) {
        LOG(logCRITICAL) << "The validation of the index file for the tumor vcf " << vcf_tumor_fname << " failed with return code " << sr_set_opt_retval;
        exit(-8);
    }
    int sr_add_reader_retval = bcf_sr_add_reader(sr, vcf_tumor_fname.c_str());
    if (sr_add_reader_retval != 1) {
        LOG(logCRITICAL) << "Failed to synchronize-read the tumor vcf " << vcf_tumor_fname << " with return code " << sr_add_reader_retval;
        exit(-9);
    }
    
    int valsize = 0; 
    int ndst_val = 0;
    char *bcfstring = NULL;
    float *bcffloats = NULL;
    int32_t *bcfints = NULL;
    
    while (bcf_sr_next_line(sr)) {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        bcf_unpack(line, BCF_UN_ALL);
        
        // skip over all symbolic alleles except MGVCF_SYMBOL and ADDITIONAL_INDEL_CANDIDATE_SYMBOL
        bool should_continue = false;
        for (uint32_t i = 1; i < line->n_allele; i++) {
            if ('<' == line->d.allele[i][0] && ((strcmp("<NON_REF>", line->d.allele[i]) && strcmp("<ADDITIONAL_INDEL_CANDIDATE>", line->d.allele[i])) || !is_tumor_format_retrieved)) {
                should_continue = true;
            }
        }

        if (should_continue) { continue; }
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "VTI", &bcfints, &ndst_val);
        if (valsize <= 0) { continue; }
        assertUVC((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for VTI and line %ld!\n", ndst_val, valsize, line->pos));
        assertUVC((2 == line->n_allele) || !fprintf(stderr, "Bcf line %ld has %d alleles!\n", line->pos, line->n_allele));
        const AlignmentSymbol symbol = AlignmentSymbol(bcfints[1]);
        
        auto symbolpos = ((isSymbolSubstitution(symbol) || MGVCF_SYMBOL == symbol || ADDITIONAL_INDEL_CANDIDATE_SYMBOL == symbol) ? (line->pos) : (line->pos + 1));
        TumorKeyInfo tki;
        tki.VTI = bcfints[1];

if (MGVCF_SYMBOL == symbol) {
        LOG(logINFO) << "gVCFblock with pos " << symbolpos << " was retrieved";
}
if (ADDITIONAL_INDEL_CANDIDATE_SYMBOL == symbol) {
        LOG(logINFO) << "ADDITIONAL_INDEL_CANDIDATE symbol with pos " << symbolpos << " was retrieved";
}

if (MGVCF_SYMBOL != symbol && ADDITIONAL_INDEL_CANDIDATE_SYMBOL != symbol) {

        BCF_GET_FORMAT_INT32B_WITH_CHECK(tki.BDP, "BDPb");
        /*
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "BDPf", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for BDPf and line %ld!\n", ndst_val, valsize, line->pos));
        tki.BDP = bcfints[0];
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "BDPr", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for BDPr and line %ld!\n", ndst_val, valsize, line->pos));
        tki.BDP += bcfints[0];
        */
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "bDPf", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for bDPf and line %ld!\n", ndst_val, valsize, line->pos));
        tki.bDP = bcfints[1];
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "bDPr", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for bDPr and line %ld!\n", ndst_val, valsize, line->pos));
        tki.bDP += bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "CDP1x", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for CDP1x and line %ld!\n", ndst_val, valsize, line->pos));
        tki.CDP1x = bcfints[0];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cDP1x", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cDP1x and line %ld!\n", ndst_val, valsize, line->pos));
        tki.cDP1x = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cVQ1", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cVQ1 and line %ld!\n", ndst_val, valsize, line->pos));
        tki.cVQ1 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cPCQ1", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cPCQ1 and line %ld!\n", ndst_val, valsize, line->pos));
        tki.cPCQ1 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "CDP2x", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for CDP2x and line %ld!\n", ndst_val, valsize, line->pos));
        tki.CDP2x = bcfints[0];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cDP2x", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cDP2x and line %ld!\n", ndst_val, valsize, line->pos));
        tki.cDP2x = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cVQ2", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cVQ2 and line %ld!\n", ndst_val, valsize, line->pos));
        tki.cVQ2 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cPCQ2", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cPCQ2 and line %ld!\n", ndst_val, valsize, line->pos));
        tki.cPCQ2 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "bNMQ", &bcfints, &ndst_val);
        assertAlways((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for bNMQ and line %ld!\n", ndst_val, valsize, line->pos));
        tki.bNMQ = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "vHGQ", &bcfints, &ndst_val);
        assertAlways((1 == ndst_val && 1 == valsize) || !fprintf(stderr, "1 == %d && 1 == %d failed for vHGQ and line %ld!\n", ndst_val, valsize, line->pos));
        tki.vHGQ = bcfints[0];
        
        // extra code for backward compatibility
        
        BCF_GET_FORMAT_INT32B_WITH_CHECK(tki.tDP, "CDP1b");
        BCF_GET_FORMAT_INT32_WITH_CHECK("cDP1f");
        tki.tADR = {{ bcfints[0], bcfints[1] }};
        BCF_GET_FORMAT_INT32_WITH_CHECK("cDP1r");
        for (int i = 0; i < 2; i++) { tki.tADR[i] += bcfints[i]; }
        
        BCF_GET_FORMAT_INT32B_WITH_CHECK(tki.tDPC, "CDP2b");
        BCF_GET_FORMAT_INT32_WITH_CHECK("cDP2f");
        tki.tADCR = {{ bcfints[0], bcfints[1] }};
        BCF_GET_FORMAT_INT32_WITH_CHECK("cDP2r");
        for (int i = 0; i < 2; i++) { tki.tADCR[i] += bcfints[i]; }
        
}
        tki.pos = line->pos;
        tki.ref_alt = als_to_string(line->d.allele, line->n_allele);
        if (is_tumor_format_retrieved) {
            tki.bcf1_record = bcf_dup(line);
        }
        
        ndst_val = 0;
        valsize = bcf_get_format_char(bcf_hdr, line, "_C2XP", &bcfstring, &ndst_val);
        tki.enable_tier2_consensus_format_tags = (valsize > 0);
        
        const auto retkey = std::make_tuple(line->rid, symbolpos, symbol);
        ret.insert(std::make_pair(retkey, std::vector<TumorKeyInfo>()));
        ret[retkey].push_back(tki);
    }
    xfree(bcfstring);
    xfree(bcfints);
    xfree(bcffloats);
    bcf_sr_destroy(sr);
    return ret;
}

template <bool TIsAnyTandemRepeat = false>
CoveredRegion<uvc1_qual_big_t> 
region_repeatvec_to_baq_offsetarr(
        const std::vector<RegionalTandemRepeat> & region_repeatvec, 
        uvc1_refgpos_t tid,
        uvc1_refgpos_t extended_inclu_beg_pos, 
        uvc1_refgpos_t extended_exclu_end_pos,
        const CommandLineArgs & paramset) {
    auto ret = CoveredRegion<uvc1_qual_big_t>(tid, extended_inclu_beg_pos, extended_exclu_end_pos);    
    uvc1_qual_big_t baq_prefixsum = 0;
    for (auto i = extended_inclu_beg_pos; i < extended_exclu_end_pos; i++) {
        auto rtr_idx = i - extended_inclu_beg_pos;
        const auto & rtr = region_repeatvec[rtr_idx];
        const auto tracklen2 = (TIsAnyTandemRepeat ? rtr.anyTR_tracklen : rtr.tracklen);
        assertUVC (rtr.begpos <= rtr_idx);
        assertUVC (rtr.unitlen > 0);
        assertUVC (tracklen2 >= rtr.unitlen);
        if (tracklen2 / rtr.unitlen >= 3 || (tracklen2 / rtr.unitlen >= 2 && tracklen2 >= (uvc1_readpos_t)round(paramset.indel_polymerase_size))) {
            baq_prefixsum += (paramset.indel_str_phred_per_region * 10) / (tracklen2) + 1;
            ret.getRefByPos(i) = baq_prefixsum;
        } else {
            baq_prefixsum += paramset.indel_nonSTR_phred_per_base * 10;
            ret.getRefByPos(i) = baq_prefixsum;
        }
    }
    for (auto i = extended_inclu_beg_pos; i < extended_exclu_end_pos; i++) {
        ret.getRefByPos(i) /= 10;
    }
    return ret;
}

bool 
are_depths_diff(uvc1_readnum_t currDP, uvc1_readnum_t prevDP, uvc1_readnum_t mul_perc_ratio, uvc1_readnum_t add_num_ratio) {
    const uvc1_readnum_t minDP = MIN(currDP, prevDP);
    const uvc1_readnum_t maxDP = MAX(currDP, prevDP);
    if (minDP * mul_perc_ratio >= maxDP * 100) {
        return false;
    }
    if (minDP + add_num_ratio >= maxDP) {
        return false;
    }
    return true;
}

void
umi_strand_readset_uvc_destroy(auto & umi_strand_readset) {
    for (auto strand_readset : umi_strand_readset) {
        for (int strand = 0; strand < 2; strand++) {
            auto readset = strand_readset.first[strand]; 
            for (auto read : readset) {
                for (bam1_t *b : read) {
                    bam_destroy1(b);
                } 
            }
        }
    }
}

template <class T>
int 
process_batch(
        std::string & uncompressed_vcf_string,
        std::array<std::string, NUM_FQLIKE_CON_OUT_FILES> & uncompressed_3fastq_string,
        BatchArg & arg,
        const T & tid_pos_symb_to_tkis) {
    
    const faidx_t *const ref_faidx = arg.ref_faidx;
    
    const bcf_hdr_t *const bcf_hdr = arg.bcf_hdr;
    const CommandLineArgs & paramset = arg.paramset;
    const std::string UMI_STRUCT_STRING = arg.UMI_STRUCT_STRING;
    const BedLine bedline = arg.bedline;
    const std::tuple<std::string, uvc1_refgpos_t> tname_tseqlen_tuple = arg.tname_tseqlen_tuple;
    const auto regionbatch_ordinal = arg.regionbatch_ordinal;
    const auto regionbatch_tot_num = arg.regionbatch_tot_num;
    const auto thread_id = arg.thread_id;
    
    bool is_loginfo_enabled = (ispowerof2(regionbatch_ordinal + 1) || ispowerof2(regionbatch_tot_num - regionbatch_ordinal) || paramset.always_log);
    std::string raw_out_string;
    std::string raw_out_string_pass;
    
    auto tid = bedline.tid;
    
    const auto incluBegPosition = bedline.beg_pos;
    const auto excluEndPosition = bedline.end_pos;
    bool end2end = (bedline.region_flag & BED_END_TO_END_BIT); 
    
    std::map<MolecularBarcode, std::pair<std::array<std::map<uvc1_hash_t, std::vector<bam1_t *>>, 2>, MolecularBarcode>> umi_to_strand_to_reads;
    uvc1_refgpos_t bam_inclu_beg_pos, bam_exclu_end_pos; 
    std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, MolecularBarcode>> umi_strand_readset;

    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts bamfname_to_strand_to_familyuid_to_reads with pair_end_merge = " << paramset.pair_end_merge; }
    std::array<uvc1_readnum_big_t, 3> passed_pcrpassed_umipassed = bamfname_to_strand_to_familyuid_to_reads(
            umi_to_strand_to_reads,
            bam_inclu_beg_pos,
            bam_exclu_end_pos,
            tid,
            incluBegPosition, 
            excluEndPosition,
            end2end,
            regionbatch_ordinal,
            regionbatch_tot_num,
            UMI_STRUCT_STRING,
            arg.samfile,
            arg.hts_idx,
            thread_id,
            paramset,
            0);
    const auto num_passed_reads = passed_pcrpassed_umipassed[0]; // -1 means that min read depth is not satisfied. 
    const auto num_pcrpassed_reads = passed_pcrpassed_umipassed[1];
    const bool is_by_capture = ((num_pcrpassed_reads) * 2 <= num_passed_reads);
    const AssayType inferred_assay_type = ((ASSAY_TYPE_AUTO == paramset.assay_type) ? (is_by_capture ? ASSAY_TYPE_CAPTURE : ASSAY_TYPE_AMPLICON) : (paramset.assay_type));
    
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts converting umi_to_strand_to_reads with is_by_capture = " << is_by_capture << "  " ;}
    fill_strand_umi_readset_with_strand_to_umi_to_reads(
            umi_strand_readset, 
            umi_to_strand_to_reads, 
            paramset,
            0); // bigger region
    
    if ((0 == num_passed_reads) || (-1 == num_passed_reads)) { 
        umi_strand_readset_uvc_destroy(umi_strand_readset);
        return -1; 
    };
    const uvc1_qual_t minABQ_snv = ((ASSAY_TYPE_AMPLICON == inferred_assay_type) ? paramset.syserr_minABQ_pcr_snv : paramset.syserr_minABQ_cap_snv);
    const uvc1_qual_t minABQ_indel = ((ASSAY_TYPE_AMPLICON == inferred_assay_type) ? paramset.syserr_minABQ_pcr_indel : paramset.syserr_minABQ_cap_indel);
    
    const uvc1_refgpos_t rpos_inclu_beg = MAX(incluBegPosition, bam_inclu_beg_pos);
    const uvc1_refgpos_t rpos_exclu_end = MIN(excluEndPosition, bam_exclu_end_pos);
    const uvc1_refgpos_t extended_inclu_beg_pos = MAX(0, non_neg_minus(MIN(incluBegPosition, bam_inclu_beg_pos), MAX_STR_N_BASES));
    const uvc1_refgpos_t extended_exclu_end_pos = MIN(std::get<1>(tname_tseqlen_tuple), MAX(excluEndPosition, bam_exclu_end_pos) + MAX_STR_N_BASES);
    
    const auto tkis_beg = tid_pos_symb_to_tkis.lower_bound(std::make_tuple(tid, extended_inclu_beg_pos    , AlignmentSymbol(0)));
    const auto tkis_end = tid_pos_symb_to_tkis.upper_bound(std::make_tuple(tid, extended_exclu_end_pos + 1, AlignmentSymbol(0)));
    std::vector<bool> extended_posidx_to_is_rescued(extended_exclu_end_pos - extended_inclu_beg_pos + 1, false);
    int num_rescued = 0;
    for (auto tkis_it = tkis_beg; tkis_it != tkis_end; tkis_it++) {
        auto symbolpos = std::get<1>(tkis_it->first);
        extended_posidx_to_is_rescued[symbolpos - extended_inclu_beg_pos] = true;
        num_rescued++;
        if (is_loginfo_enabled) {
            // The short deletion at 22:17946835 in NA12878-NA24385 mixture is overwhelmed by the false positve long del spanning the true positive short del.
            // Due to the VCF specs and left alignment as best practice, the symmetry between left and right alignments are broken. 
            // Thus, this variant is true positive if using left alignemnt but false positive if using right alignment. 
            // TODO: modifications to the VCF specs or somatic-variant definitions, which is beyond the scope of this program.
            LOG(logDEBUG4) << "Thread " << thread_id << " iterated over symbolpos " << symbolpos << " and symbol " << std::get<2>(tkis_it->first) << " as a rescued var";
        }
    } // bigger region
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " deals with " << num_rescued << " tumor-sample variants in region " << extended_inclu_beg_pos << " to " << extended_exclu_end_pos + 1 ;}
    
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts constructing symbolToCountCoverageSet12 with " << extended_inclu_beg_pos << (" , ") << extended_exclu_end_pos; }
    // + 1 accounts for insertion at the end of the region, this should happen rarely for only re-aligned reads at around once per one billion base pairs
    if (is_loginfo_enabled) { LOG(logINFO)<< "Thread " << thread_id << " starts updateByRegion3Aln with " << umi_strand_readset.size() << " families"; }
    const std::string refstring = load_refstring(ref_faidx, tid, extended_inclu_beg_pos, extended_exclu_end_pos);
    std::vector<RegionalTandemRepeat> region_repeatvec = refstring2repeatvec(
            refstring, 
            paramset.indel_str_repeatsize_max,
            paramset.indel_vntr_repeatsize_max, // https://en.wikipedia.org/wiki/Variable_number_tandem_repeat
            paramset.indel_BQ_max,
            paramset.indel_polymerase_slip_rate,
            paramset.indel_del_to_ins_err_ratio,
            0);
    const auto & baq_offsetarr = region_repeatvec_to_baq_offsetarr(region_repeatvec, tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1, paramset);
    const auto & baq_offsetarr2 = region_repeatvec_to_baq_offsetarr<true>(region_repeatvec, tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1, paramset);

    assertUVC(((baq_offsetarr.getExcluEndPosition() - baq_offsetarr.getIncluBegPosition()) == UNSIGN2SIGN(region_repeatvec.size())) 
            || !fprintf(stderr, "%d - %d == %lu failed (baq == repeat in size)!\n", baq_offsetarr.getExcluEndPosition(), baq_offsetarr.getIncluBegPosition(), region_repeatvec.size()));

    // begin of smaller regions
    Symbol2CountCoverageSet symbolToCountCoverageSet12(tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1); 
    
    std::vector<HapLink> mutform2count4vec_bq;
    std::vector<HapLink> mutform2count4vec_fq;
    std::vector<HapLink> mutform2count4vec_f2q;
    std::array<std::string, NUM_FQLIKE_CON_OUT_FILES> fqdata3;
    symbolToCountCoverageSet12.updateByRegion3Aln(
            fqdata3,
            mutform2count4vec_bq,
            mutform2count4vec_fq,
            mutform2count4vec_f2q,
            umi_strand_readset,
            
            refstring,
            region_repeatvec,
            baq_offsetarr,
            baq_offsetarr2,
            
            arg.prev_bedline,
            arg.bedline,
            
            paramset,
            0);

    std::string buf_out_string_pass;

if (paramset.inferred_is_vcf_generated) {

    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts analyzing phasing info"; }
    auto simplemut2indices_bq = mutform2count4vec_to_simplemut2indices(mutform2count4vec_bq);
    auto simplemut2indices_fq = mutform2count4vec_to_simplemut2indices(mutform2count4vec_fq);
    auto simplemut2indices_f2q = mutform2count4vec_to_simplemut2indices(mutform2count4vec_f2q);

    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " starts generating block gzipped vcf"; }
    
    const std::set<size_t> empty_size_t_set;
    uvc1_readpos_t prev_tracklen = 0;
    uvc1_readpos_t curr_tracklen = 0;

    for (uvc1_refgpos_t zerobased_pos = rpos_inclu_beg; zerobased_pos <= rpos_exclu_end; zerobased_pos++, prev_tracklen = curr_tracklen) {
        std::string repeatunit;
        uvc1_readpos_t repeatnum = 0;
        
        uvc1_rp_diff_t rridx = zerobased_pos - extended_inclu_beg_pos;
        indelpos_to_context(repeatunit, repeatnum, refstring, rridx, paramset.indel_str_repeatsize_max);
        curr_tracklen = repeatnum * UNSIGN2SIGN(repeatunit.size());
        
        const std::array<AlignmentSymbol, 2> symboltype_to_refsymbol = {{
                ((UNSIGN2SIGN(refstring.size()) == ((zerobased_pos - 1) - extended_inclu_beg_pos) || (-1 == ((zerobased_pos - 1) - extended_inclu_beg_pos))) 
                ? BASE_NN : CHAR_TO_SYMBOL.data.at(refstring.at((zerobased_pos - 1) - extended_inclu_beg_pos))),
                LINK_M, 
        }};
        const auto refsize = UNSIGN2SIGN(refstring.size());
        const auto refidx = zerobased_pos - extended_inclu_beg_pos;
        const AlignmentSymbol prev_base1 = ((refidx     >= 2      ) ? CHAR_TO_SYMBOL.data.at(refstring.at((refidx - 2))) : BASE_NN);
        const AlignmentSymbol prev_base2 = ((refidx     >= 3      ) ? CHAR_TO_SYMBOL.data.at(refstring.at((refidx - 3))) : BASE_NN);
        const AlignmentSymbol next_base1 = ((refidx     <  refsize) ? CHAR_TO_SYMBOL.data.at(refstring.at((refidx)))     : BASE_NN);
        const AlignmentSymbol next_base2 = ((refidx + 1 <  refsize) ? CHAR_TO_SYMBOL.data.at(refstring.at((refidx + 1))) : BASE_NN);
        
        std::array<bcfrec::BcfFormat, 2> st_to_init_fmt = {{bcfrec::BcfFormat(), bcfrec::BcfFormat()}};
        
        std::array<std::vector<std::tuple<bcfrec::BcfFormat, TumorKeyInfo>>, NUM_SYMBOL_TYPES> st_to_fmt_tki_tup_vec;
        uvc1_readnum_t ins_bdepth = 0;
        uvc1_readnum_t del_bdepth = 0;
        uvc1_readnum_t ins_cdepth = 0;
        uvc1_readnum_t del_cdepth = 0;
        
        uvc1_readnum_t ins1_bdepth = 0;
        uvc1_readnum_t del1_bdepth = 0;
        uvc1_readnum_t ins1_cdepth = 0;
        uvc1_readnum_t del1_cdepth = 0;
        
        for (const SymbolType symboltype : SYMBOL_TYPE_ARR)
        {
            if (zerobased_pos == rpos_inclu_beg && BASE_SYMBOL == symboltype) { continue; } 
            const uvc1_refgpos_t refpos = (BASE_SYMBOL == symboltype ? (zerobased_pos - 1) : zerobased_pos);
            
            TumorKeyInfo THE_DUMMY_TUMOR_KEY_INFO;
            const AlignmentSymbol refsymbol = symboltype_to_refsymbol[symboltype];
            std::array<uvc1_readnum_t, 2> bDPcDP = BcfFormat_symboltype_init(
                    st_to_init_fmt[symboltype], 
                    symbolToCountCoverageSet12, 
                    refpos, 
                    symboltype, 
                    refsymbol, 
                    0);
            if ((paramset.outvar_flag & OUTVAR_MGVCF) && ((((refpos) % MGVCF_REGION_MAX_SIZE) == 0) 
                    || (refpos == incluBegPosition)) && (SYMBOL_TYPE_ARR[0] == symboltype)) {
                const auto & frag_format_depth_sets = symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets;
                const auto & fam_format_depth_sets = symbolToCountCoverageSet12.symbol_to_fam_format_depth_sets_2strand;
                
                const uvc1_qual_t init_refQ = (INT_MAX / 2 + 1);
                uvc1_readnum_t prev_tot_bdepth = 0;
                uvc1_readnum_t prev_tot_cdepth = 0;
                uvc1_readnum_t prev_tot_cdep12 = 0;
                uvc1_qual_t prev_refQ = init_refQ;
                std::vector<uvc1_rp_diff_t> pos_stype_BDP_CDP_refQ_1dvec;
                const auto rp2end = MIN(refpos + MGVCF_REGION_MAX_SIZE + 1, symbolToCountCoverageSet12.getUnifiedExcluEndPosition());
                for (uvc1_refgpos_t rp2 = refpos; rp2 < rp2end; rp2++) {
                    for (SymbolType stype : SYMBOL_TYPES_IN_VCF_ORDER) {
                        const uvc1_rp_diff_t refstring_offset = rp2 - extended_inclu_beg_pos;
                        if (refstring_offset > UNSIGN2SIGN(refstring.size())) { 
                            fprintf(stderr, "The refstring offset %d at tid %d pos %d is invalid!\n\n", refstring_offset, tid, refpos); 
                            abort(); 
                        }
                        const AlignmentSymbol base_m = ((refstring_offset < UNSIGN2SIGN(refstring.size())) 
                            ? CHAR_TO_SYMBOL.data[refstring.substr(refstring_offset, 1)[0]] 
                            : BASE_N);
                        const auto refsymbol = ((BASE_SYMBOL == stype) ? (base_m) : (LINK_M));
                        uvc1_readnum_t curr_tot_bdepth = 0;
                        uvc1_readnum_t curr_tot_cdepth = 0;
                        uvc1_readnum_t curr_tot_cdep12 = 0;
                        for (int strand = 0; strand < 2; strand++) {
                            curr_tot_bdepth += formatSumBySymbolType(frag_format_depth_sets[strand].getByPos(rp2), stype, FRAG_bDP);
                            curr_tot_cdepth += formatSumBySymbolType(fam_format_depth_sets[strand].getByPos(rp2), stype, FAM_cDP1);
                            curr_tot_cdep12 += formatSumBySymbolType(fam_format_depth_sets[strand].getByPos(rp2), stype, FAM_cDP12);
                        }
                        const auto ref_cdepth = 
                                fam_format_depth_sets[0].getByPos(rp2)[refsymbol][FAM_cDP12] + 
                                fam_format_depth_sets[1].getByPos(rp2)[refsymbol][FAM_cDP12];
                        const auto nonref_cdepth = curr_tot_cdep12 - ref_cdepth;
                        
                        // contam, lower  nonref-FA -> close to zero qual
                        const auto ref_like_binom = -calc_binom_10log10_likeratio(paramset.contam_any_mul_frac, nonref_cdepth + 0.5, curr_tot_cdepth + 1.0);
                        const auto ref_like_powlaw = -MAX(0, paramset.powlaw_exponent * (10/log(10)) 
                                * logit2((nonref_cdepth + 0.5) / (curr_tot_cdepth + 1.0), paramset.contam_any_mul_frac));
                        
                        // hetero, higher nonref-FA (lower ref-FA) -> close to zero qual
                        const auto nonref_like_binom = -calc_binom_10log10_likeratio(paramset.germ_hetero_FA, ref_cdepth + 0.5, curr_tot_cdepth + 1.0);
                        const auto nonref_like_powlaw = -MAX(0, paramset.powlaw_exponent * (10/log(10)) 
                                * logit2((ref_cdepth + 0.5) / (curr_tot_cdepth + 1.0), paramset.germ_hetero_FA));
                        
                        const uvc1_qual_t curr_refQ = paramset.germ_phred_hetero_snp 
                                + (uvc1_qual_t)round(MAX(ref_like_binom, ref_like_powlaw) 
                                - (uvc1_qual_t)round(MAX(nonref_like_binom, nonref_like_powlaw)));
                        if ((init_refQ == prev_refQ) || (abs(curr_refQ - prev_refQ) > 10)
                                || are_depths_diff(curr_tot_bdepth, prev_tot_bdepth, 100 + 30, 3)
                                || are_depths_diff(curr_tot_cdepth, prev_tot_cdepth, 100 + 30, 3)
                                || are_depths_diff(curr_tot_cdep12, prev_tot_cdep12, 100 + 30, 3)) {
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(rp2 + ((BASE_SYMBOL == stype) ? (1) : (0)));
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(1+(int32_t)stype);
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(INT32_MIN);
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(curr_tot_bdepth);
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(curr_tot_cdepth);
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(curr_tot_cdep12);
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(curr_refQ);
                            pos_stype_BDP_CDP_refQ_1dvec.push_back(INT32_MIN);
                            prev_tot_bdepth = curr_tot_bdepth;
                            prev_tot_cdepth = curr_tot_cdepth;
                            prev_tot_cdep12 = curr_tot_cdep12;
                            prev_refQ = curr_refQ;
                        }
                    }
                }
                const auto vcfREF = refstring.substr(refpos - extended_inclu_beg_pos, 1);
                const AlignmentSymbol match_refsymbol = CHAR_TO_SYMBOL.data[vcfREF[0]];
                const std::string gvcf_blockline = string_join(std::vector<std::string>{{
                    std::get<0>(tname_tseqlen_tuple), 
                    std::to_string(refpos + 1),
                    std::string("."), 
                    vcfREF,
                    std::string("<NON_REF>"),
                    std::string("."), 
                    std::string("."), 
                    std::string("MGVCF_BLOCK"),
                    std::string("GT:VTI:POS_VT_BDP_CDP_HomRefQ"),
                    std::string(".") + ":" + std::to_string(match_refsymbol) + "," + std::to_string(MGVCF_SYMBOL) + ":" 
                            + int32t_join(pos_stype_BDP_CDP_refQ_1dvec, ",") + "," + std::to_string(rp2end),
                }}, "\t");
                
                std::string tumor_gvcf_format = "";
                if (paramset.is_tumor_format_retrieved && IS_PROVIDED(paramset.vcf_tumor_fname)) { 
                    const auto tkis_it = tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, MGVCF_SYMBOL));
                    if (tkis_it != tid_pos_symb_to_tkis.end()) {
                        const auto & tkis = tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, MGVCF_SYMBOL))->second;
                        if (tkis.size() == 1) {
                            tumor_gvcf_format = bcf1_to_string(bcf_hdr, tkis[0].bcf1_record);
                            LOG(logDEBUG4) << "gVCFblock at " << refpos << " is indeed found, tumor_gvcf_format == " << tumor_gvcf_format; 
                        } else {
                            tumor_gvcf_format = std::string("\t.:.,.:-1");
                            LOG(logDEBUG4) << "gVCFblock at " << refpos << " is not found, tkis.size() == " << tkis.size();
                        }
                    } else {
                        tumor_gvcf_format = std::string("\t.:.,.:.");
                        LOG(logDEBUG4) << "gVCFblock at " << refpos << " is not found at all.";
                    }
                }
                buf_out_string_pass += gvcf_blockline + tumor_gvcf_format + "\n";
            }
            
            const auto aCDP = symbolToCountCoverageSet12.seg_format_prep_sets.getByPos(refpos).segprep_a_near_long_clip_dp;
            const auto ADP = symbolToCountCoverageSet12.seg_format_prep_sets.getByPos(refpos).segprep_a_dp;
            
            const bool is_in_long_track = (curr_tracklen > MAX(paramset.microadjust_alignment_tracklen_min - 1, prev_tracklen));
            const bool is_in_clip_region = ((aCDP >= paramset.microadjust_alignment_clip_min_count) 
                    && (aCDP >= ADP * (paramset.microadjust_alignment_clip_min_frac - DBL_EPSILON)));
            if ((OUTVAR_ADDITIONAL_INDEL_CANDIDATE & paramset.outvar_flag)
                    && (SYMBOL_TYPE_ARR[0] == symboltype)
                    && (is_in_long_track || is_in_clip_region)
                    && (ADP >= 2 * paramset.microadjust_alignment_clip_min_count)) {
                const auto vcfREF = refstring.substr(refpos - extended_inclu_beg_pos, 1);
                const AlignmentSymbol match_refsymbol = CHAR_TO_SYMBOL.data[vcfREF[0]];
                const std::string vcfline = string_join(std::vector<std::string>{{
                    std::get<0>(tname_tseqlen_tuple), // chrom
                    std::to_string(refpos + 1), // pos
                    std::string("."), // id
                    vcfREF, // ref
                    SYMBOL_TO_DESC_ARR[ADDITIONAL_INDEL_CANDIDATE_SYMBOL], // alt
                    std::string("."), // qual
                    std::string("."), // filter
                    (std::string("ADDITIONAL_INDEL_CANDIDATE;RU=") + repeatunit + ";RC=" + std::to_string(repeatnum)), // info
                    std::string("GT:VTI:clipDP"), // format
                    std::string(".") + ":" + std::to_string(match_refsymbol) + "," + std::to_string(ADDITIONAL_INDEL_CANDIDATE_SYMBOL) 
                            + ":" + std::to_string(ADP) + "," + std::to_string(aCDP) // format values
                }}, "\t");
                std::string tumor_format = "";
                if (paramset.is_tumor_format_retrieved && IS_PROVIDED(paramset.vcf_tumor_fname)) { 
                    const auto tkis_it = tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, ADDITIONAL_INDEL_CANDIDATE_SYMBOL));
                    if (tkis_it != tid_pos_symb_to_tkis.end()) {
                        const auto & tkis = tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, ADDITIONAL_INDEL_CANDIDATE_SYMBOL))->second;
                        if (tkis.size() == 1) {
                            tumor_format = bcf1_to_string(bcf_hdr, tkis[0].bcf1_record);
                        } else {
                            tumor_format = std::string("\t.:-1,-1:-1,-1");
                        }
                    } else {
                        tumor_format = std::string("\t.:.,.:.,.");
                    }
                }
                buf_out_string_pass += vcfline + tumor_format + "\n";
            }
            
            const auto ref_bdepth = 
                        symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[0].getByPos(refpos)[refsymbol][FRAG_bDP]
                      + symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[1].getByPos(refpos)[refsymbol][FRAG_bDP];
            
            for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symboltype])
            {
                const bool is_pos_rescued = (IS_PROVIDED(paramset.vcf_tumor_fname) && (extended_posidx_to_is_rescued[refpos - extended_inclu_beg_pos]));
                const bool is_var_rescued = (is_pos_rescued && (tid_pos_symb_to_tkis.end() != tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, symbol)))); 
                const auto bdepth = 
                        symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[0].getByPos(refpos)[symbol][FRAG_bDP]
                      + symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[1].getByPos(refpos)[symbol][FRAG_bDP];
                const auto cdepth = 
                        MAX(symbolToCountCoverageSet12.symbol_to_fam_format_depth_sets_2strand[0].getByPos(refpos)[symbol][FAM_cDP1],
                            symbolToCountCoverageSet12.symbol_to_fam_format_depth_sets_2strand[0].getByPos(refpos)[symbol][FAM_cDP12])
                      + MAX(symbolToCountCoverageSet12.symbol_to_fam_format_depth_sets_2strand[1].getByPos(refpos)[symbol][FAM_cDP1],
                            symbolToCountCoverageSet12.symbol_to_fam_format_depth_sets_2strand[1].getByPos(refpos)[symbol][FAM_cDP12]);
                if (isSymbolIns(symbol)) {
                    ins_bdepth += bdepth;
                    ins_cdepth += cdepth;
                    if (LINK_I1 == symbol) {
                        ins1_bdepth += bdepth;
                        ins1_cdepth += cdepth;
                    }
                } else if (isSymbolDel(symbol)) {
                    del_bdepth += bdepth;
                    del_cdepth += cdepth;
                    if (LINK_D1 == symbol) {
                        del1_bdepth += bdepth;
                        del1_cdepth += cdepth;
                    }
                }
                if ((ISNT_PROVIDED(paramset.vcf_tumor_fname))
                        &&    (((refsymbol != symbol) && (bdepth < paramset.min_altdp_thres))
                            || ((refsymbol == symbol) && (bDPcDP[0] - ref_bdepth < paramset.min_altdp_thres)))
                        && (!paramset.should_output_all)) {
                    continue;
                }
                
                if (IS_PROVIDED(paramset.vcf_tumor_fname) && (!is_pos_rescued)) {
                    continue;
                }
                const auto simplemut = std::make_pair(refpos, symbol);
                const auto indices_bq = (simplemut2indices_bq.find(simplemut) != simplemut2indices_bq.end() ? simplemut2indices_bq[simplemut] : empty_size_t_set); 
                const auto indices_fq = (simplemut2indices_fq.find(simplemut) != simplemut2indices_fq.end() ? simplemut2indices_fq[simplemut] : empty_size_t_set);
                const auto indices_f2q = (simplemut2indices_f2q.find(simplemut) != simplemut2indices_f2q.end() ? simplemut2indices_f2q[simplemut] : empty_size_t_set);
                std::vector<TumorKeyInfo> tkis;
                if (is_var_rescued) {
                    tkis = tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, symbol))->second;
                }
                
                auto fmt = st_to_init_fmt[symboltype];
                int tkiidx = 0;
                std::vector<std::tuple<uvc1_readnum_t, uvc1_readnum_t, std::string, TumorKeyInfo>> bcad0a_indelstring_tki_vec;
                if (isSymbolIns(symbol) || isSymbolDel(symbol)) {
                    for (int strand = 0; strand < 2; strand++) {
                        if (0 < symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[strand].getByPos(refpos)[symbol][FRAG_bDP]) {
                            fill_by_indel_info(
                                    fmt, 
                                    symbolToCountCoverageSet12,
                                    strand,
                                    refpos,
                                    symbol,
                                    refstring,
                                    0);
                        }
                    }
                    if (is_var_rescued) {
                        for (const auto & tki : tkis) {
                            const auto tabpos = tki.ref_alt.find("\t");
                            const auto & vcfref = tki.ref_alt.substr(0, tabpos);
                            const auto & vcfalt = tki.ref_alt.substr(tabpos + 1);
                            std::string indelstring;
                            if (vcfref.size() > vcfalt.size()) {
                                indelstring = vcfref.substr(vcfalt.size());
                            } else {
                                assertUVC (vcfref.size() < vcfalt.size());
                                indelstring = vcfalt.substr(vcfref.size());
                            }
                            bcad0a_indelstring_tki_vec.push_back(std::make_tuple(bdepth, cdepth, indelstring, tki));
                            tkiidx++;
                        }
                    } else {
                        const auto bcad0a_arr_indelstring_pairvec = indel_get_majority(
                                fmt,
                                std::get<0>(tname_tseqlen_tuple).c_str(), 
                                refpos, 
                                symbol, 
                                (IS_PROVIDED(paramset.vcf_tumor_fname)),
                                0);
                        for (const auto & bcad0a_arr_indelstring_pair : bcad0a_arr_indelstring_pairvec) {
                            bcad0a_indelstring_tki_vec.push_back(std::make_tuple(
                                    bcad0a_arr_indelstring_pair.first[0], 
                                    bcad0a_arr_indelstring_pair.first[1], 
                                    bcad0a_arr_indelstring_pair.second, 
                                    THE_DUMMY_TUMOR_KEY_INFO));
                        }
                    }
                } else if (is_var_rescued) {
                    for (const auto & tki : tkis) {
                        bcad0a_indelstring_tki_vec.push_back(std::make_tuple(bdepth, cdepth, std::string(""), tki));
                    }
                } else {
                    bcad0a_indelstring_tki_vec.push_back(std::make_tuple(bdepth, cdepth, std::string(""), THE_DUMMY_TUMOR_KEY_INFO));
                }
                
                for (const auto & bcad0a_indelstring_tki : bcad0a_indelstring_tki_vec) {
                    const auto & indelstring = std::get<2>(bcad0a_indelstring_tki);
                    const auto & tki = std::get<3>(bcad0a_indelstring_tki);
                    const bool is_homopol_1bp = (prev_base1 == refsymbol && next_base1 == refsymbol);
                    const bool is_homopol_2bp = (prev_base2 == refsymbol && next_base2 == refsymbol);
                    BcfFormat_symbol_init(
                            fmt,
                            symbolToCountCoverageSet12,
                            refpos,
                            symbol,
                            mutform2count4vec_bq,     
                            indices_bq,
                            mutform2count4vec_fq, 
                            indices_fq,
                            mutform2count4vec_f2q, 
                            indices_f2q,

                            std::get<0>(bcad0a_indelstring_tki),
                            std::get<1>(bcad0a_indelstring_tki),
                            indelstring,
                            (isSymbolSubstitution(symbol) 
                                    ? non_neg_minus(minABQ_snv, (is_homopol_1bp ? (is_homopol_2bp ? 20 : 10) : 0))
                                    : minABQ_indel),
                            paramset,
                            0);
                    BcfFormat_symbol_calc_DPv(
                            fmt,
                            region_repeatvec[MAX(refpos - extended_inclu_beg_pos, 3) - 3],
                            region_repeatvec[MIN(refpos - extended_inclu_beg_pos + 3, UNSIGN2SIGN(region_repeatvec.size()) - 1)],
                            ((is_var_rescued && (tki.VTI == LAST(fmt.VTI))) ? ((double)(tki.cDP1x + 1) / (double)(tki.CDP1x + 2)) : -1.0), /* tpfa */
                            tki,
                            refsymbol,
                            symbolToCountCoverageSet12,
                            refpos,
                            paramset,
                            0);
                    st_to_fmt_tki_tup_vec[symboltype].push_back(std::make_tuple(fmt, std::get<3>(bcad0a_indelstring_tki)));
                }
            }
        } // end of iterations within symboltype
        const size_t string_pass_old_size = buf_out_string_pass.size();
        auto st_to_nlodq_fmtptr1_fmtptr2_tup = std::array<std::tuple<uvc1_qual_t, bcfrec::BcfFormat*, bcfrec::BcfFormat*>, NUM_SYMBOL_TYPES>();
        std::array<int32_t, NUM_SYMBOL_TYPES> curr_vAC = {{ 0 }};
        for (SymbolType symboltype : SYMBOL_TYPE_ARR) {
            if (zerobased_pos == rpos_inclu_beg && BASE_SYMBOL == symboltype) { continue; } 
            const uvc1_refgpos_t refpos = (BASE_SYMBOL == symboltype ? (zerobased_pos - 1) : zerobased_pos);
            const AlignmentSymbol refsymbol = symboltype_to_refsymbol[symboltype];
            
            auto & fmt_tki_tup_vec = st_to_fmt_tki_tup_vec[symboltype];
            
            if (fmt_tki_tup_vec.size() == 0) { continue; }
                BcfFormat_symbol_sum_DPv(fmt_tki_tup_vec);
                
                std::vector<std::tuple<uvc1_qual_t, uvc1_qual_t, uvc1_qual_t, AlignmentSymbol, std::string>> maxVQ_VQ1_VQ2_symbol_indelstr_tup_vec;
                for (auto & fmt_tki_tup : fmt_tki_tup_vec) 
                {
                    const auto VTI = LAST(std::get<0>(fmt_tki_tup).VTI);
                    assertUVC((VTI >= 0 && VTI <= END_ALIGNMENT_SYMBOLS) || !fprintf(stderr, "VTI %d at pos %d is invalid!\n", VTI, refpos));
                    AlignmentSymbol symbol = (VTI >= 0 ? AlignmentSymbol(VTI) : END_ALIGNMENT_SYMBOLS);
                    auto & fmt = std::get<0>(fmt_tki_tup);
                    const auto & tki = std::get<1>(fmt_tki_tup);
                    BcfFormat_symbol_calc_qual(
                            fmt,
                            
                            ins_cdepth,
                            del_cdepth,
                            
                            ins1_cdepth,
                            del1_cdepth,
                            
                            repeatunit,
                            repeatnum,
                            
                            IS_PROVIDED(paramset.vcf_tumor_fname),
                            region_repeatvec.at(MAX(refpos - extended_inclu_beg_pos, 3) - 3),
                            region_repeatvec.at(MIN(refpos - extended_inclu_beg_pos + 3, UNSIGN2SIGN(region_repeatvec.size()) - 1)),
                            tid,
                            refpos,
                            refsymbol,
                            ((IS_PROVIDED(paramset.vcf_tumor_fname) && (tki.VTI == LAST(fmt.VTI))) ? 
                                ((double)(tki.bDP + 0.5) / (double)(tki.BDP + 1.0)) : -1.0),
                            symbolToCountCoverageSet12,
                            paramset,
                            0);
                    if (refsymbol != symbol) {
                        uvc1_qual_t cVQ1 = LAST(std::get<0>(fmt_tki_tup).cVQ1);
                        uvc1_qual_t cVQ2 = LAST(std::get<0>(fmt_tki_tup).cVQ2);
                        maxVQ_VQ1_VQ2_symbol_indelstr_tup_vec.push_back(std::make_tuple(MAX(cVQ1, cVQ2), cVQ1, cVQ2, symbol, LAST(std::get<0>(fmt_tki_tup).gapSa)));
                        const uvc1_qual_t germ_phred_het3al = ((BASE_SYMBOL == symboltype) ? paramset.germ_phred_het3al_snp : paramset.germ_phred_het3al_indel);
                        if (MAX(cVQ1, cVQ2) >= germ_phred_het3al) {
                            curr_vAC[symboltype] += 1;
                        }
                    }
                }
                std::sort(maxVQ_VQ1_VQ2_symbol_indelstr_tup_vec.rbegin(), maxVQ_VQ1_VQ2_symbol_indelstr_tup_vec.rend());
                for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                    std::get<0>(fmt_tki_tup).cVQ1M = {{ -999 }};
                    std::get<0>(fmt_tki_tup).cVQ2M = {{ -999 }};
                    std::get<0>(fmt_tki_tup).cVQAM = {{ SYMBOL_TO_DESC_ARR[END_ALIGNMENT_SYMBOLS] }};
                    std::get<0>(fmt_tki_tup).cVQSM = {{ "" }};
                    size_t tup_vec_idx = 0;
                    for (const auto & maxVQ_VQ1_VQ2_symbol_indelstr_tup : maxVQ_VQ1_VQ2_symbol_indelstr_tup_vec) {
                        std::get<0>(fmt_tki_tup).cVQ1M[tup_vec_idx] = std::get<1>(maxVQ_VQ1_VQ2_symbol_indelstr_tup);
                        std::get<0>(fmt_tki_tup).cVQ2M[tup_vec_idx] = std::get<2>(maxVQ_VQ1_VQ2_symbol_indelstr_tup);
                        std::get<0>(fmt_tki_tup).cVQAM[tup_vec_idx] = SYMBOL_TO_DESC_ARR[std::get<3>(maxVQ_VQ1_VQ2_symbol_indelstr_tup)];
                        std::get<0>(fmt_tki_tup).cVQSM[tup_vec_idx] = std::get<4>(maxVQ_VQ1_VQ2_symbol_indelstr_tup);
                        tup_vec_idx += 1;
                        if (std::get<0>(fmt_tki_tup).cVQ1M.size() == tup_vec_idx) { break; }
                    }
                }
                auto reffmt = st_to_init_fmt[symboltype];
                auto & init_fmt = st_to_init_fmt[symboltype];
                bool is_ref_found = false;
                for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                    if (refsymbol == (AlignmentSymbol)(LAST(std::get<0>(fmt_tki_tup).VTI))) {  
                        reffmt = std::get<0>(fmt_tki_tup); 
                        is_ref_found = true; 
                    }
                }
                if(!is_ref_found) {
                    fprintf(stderr, "The position %s:%d with symboltype %d has no REF allele!\n",
                            std::get<0>(tname_tseqlen_tuple).c_str(), refpos, symboltype);
                    abort();
                }
                for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                    streamFrontPushBcfFormatR(std::get<0>(fmt_tki_tup), reffmt);
                }
                                    
                std::vector<std::pair<AlignmentSymbol, bcfrec::BcfFormat*>> symbol_format_vec;
                for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                    auto & fmt = std::get<0>(fmt_tki_tup);
                    auto symbol = (AlignmentSymbol)(LAST(fmt.VTI));
                    if (paramset.should_add_note) {
                        fmt.note += std::string("/symb/") + std::to_string(symbol) + 
                            std::string("/gVQ1/CONTQ/") + std::to_string(LAST(fmt.gVQ1)) + "/" + std::to_string(LAST(fmt.CONTQ))  + "//";
                    }
                    if (symbol != BASE_NN) {
                        symbol_format_vec.push_back(std::make_pair(symbol, &fmt));
                    }
                }
                clear_push(init_fmt.VTI, (int32_t)END_ALIGNMENT_SYMBOLS); 
                clear_push(init_fmt.gVQ1, 0); // can use a very negative number to force out all homref alleles
                clear_push(init_fmt.CONTQ, 0);
                clear_push(init_fmt.cDP12f, 0);
                clear_push(init_fmt.cDP12r, 0);
                clear_push(init_fmt.cDP0a, 0);
                clear_push(init_fmt.cDP1v, 50);
                while (symbol_format_vec.size() <= 4) {
                    symbol_format_vec.push_back(std::make_pair(END_ALIGNMENT_SYMBOLS, &init_fmt));
                }
                auto nlodq_fmtptr1_fmtptr2_tup = output_germline(
                        buf_out_string_pass,
                        refsymbol,
                        symbol_format_vec,
                        std::get<0>(tname_tseqlen_tuple).c_str(),
                        refstring,
                        refpos,
                        extended_inclu_beg_pos,
                        paramset,
                        0);
                st_to_nlodq_fmtptr1_fmtptr2_tup[symboltype] = nlodq_fmtptr1_fmtptr2_tup;
                for (auto & fmt_tki_tup : st_to_fmt_tki_tup_vec[symboltype]) {
                    auto & tmpref = std::get<0>(fmt_tki_tup).vNLODQ[symboltype];
                    tmpref = std::get<0>(nlodq_fmtptr1_fmtptr2_tup);
                }
        } // end of iterations within symboltype
        const bool is_germline_var_generated = (buf_out_string_pass.size() > string_pass_old_size);
        for (const auto symboltype : SYMBOL_TYPE_ARR) {
            if (zerobased_pos == rpos_inclu_beg && BASE_SYMBOL == symboltype) { continue; } 
            const uvc1_refgpos_t refpos = (BASE_SYMBOL == symboltype ? (zerobased_pos - 1) : zerobased_pos);
            const auto refsymbol = symboltype_to_refsymbol[symboltype];
            
            auto & fmt_tki_tup_vec = st_to_fmt_tki_tup_vec[symboltype];
            const auto & nlodq_fmtptr1_fmtptr2_tup = st_to_nlodq_fmtptr1_fmtptr2_tup[symboltype];
            
            uvc1_qual_t nlodq = ((BASE_SYMBOL == symboltype) ? paramset.germ_phred_hetero_snp : paramset.germ_phred_hetero_indel);
            AlignmentSymbol argmin_nlodq_symbol = END_ALIGNMENT_SYMBOLS; 
            for (auto & fmt_tki_tup : fmt_tki_tup_vec) 
            {
                auto & fmt = std::get<0>(fmt_tki_tup);
                assertUVC(fmt.vAC.size() == curr_vAC.size());
                for (size_t i = 0; i < fmt.vAC.size(); i++) { fmt.vAC[i] = curr_vAC[i]; }
                const auto & symbol = (AlignmentSymbol)(LAST(fmt.VTI));
                auto & tki = std::get<1>(fmt_tki_tup);
                const bool will_generate_out = (ISNT_PROVIDED(paramset.vcf_tumor_fname)
                        ? (paramset.outvar_flag & OUTVAR_ANY)
                        : (tki.ref_alt.size() > 0 && (paramset.outvar_flag & OUTVAR_SOMATIC)));
                const bool is_out_blocked = (
                        ((BASE_NN == symbol) && !(OUTVAR_BASE_NN & paramset.outvar_flag)) 
                     || ((LINK_NN == symbol) && !(OUTVAR_LINK_NN & paramset.outvar_flag)));
                if (will_generate_out && !is_out_blocked) {
                    fmt.GT = "./1";
                    const auto germ_phred = (isSymbolSubstitution(symbol) ? paramset.germ_phred_hetero_snp : paramset.germ_phred_hetero_indel);
                    const auto nlodq_singlesite = std::get<0>(nlodq_fmtptr1_fmtptr2_tup);
                    // the 3 accounts for mapping error and copy-number variation
                    const auto nlodq_singlesample = nlodq_singlesite - 3 + germ_phred;
                    
                    assertUVC ((ISNT_PROVIDED(paramset.vcf_tumor_fname)) == (0 == tki.ref_alt.size()));
                    if (IS_PROVIDED(paramset.vcf_tumor_fname)) {
                        uvc1_qual_t nlodq_inc = 999;
                        const auto fmtptrs = std::vector<bcfrec::BcfFormat*> {{ std::get<1>(nlodq_fmtptr1_fmtptr2_tup), std::get<2>(nlodq_fmtptr1_fmtptr2_tup) }} ;
                        for (const auto *fmtptr : fmtptrs) {
                            const AlignmentSymbol normsymbol = AlignmentSymbol(LAST(fmtptr->VTI));
                            auto bgerr_norm_max_ad = collectget(fmtptr->cDP1x, 1, 50); 
                            double tAD = (tki.cDP1x + 1*50) / 100.0;
                            double tDP = (tki.CDP1x + 2*50) / 100.0;
                            double nAD = (bgerr_norm_max_ad + 1*50) / 100.0;
                            double nDP = (fmtptr->CDP1x[0] + 2*50) / 100.0;
                            double bjpfrac = ((tAD) / (tDP)) / ((nAD) / (nDP));
                            uvc1_qual_t binom_b10log10like = calc_binom_10log10_likeratio((tDP - tAD) / (tDP), nDP - nAD, nAD);
                            uvc1_qual_t powlaw_b10log10like = (paramset.powlaw_exponent * 10 / log(10) * log(bjpfrac));
                            uvc1_qual_t phred_het3al_chance_inc_snp = 2 * paramset.germ_phred_hetero_snp - paramset.germ_phred_het3al_snp;
                            uvc1_qual_t phred_het3al_chance_inc_indel = 2 * paramset.germ_phred_hetero_indel - paramset.germ_phred_het3al_indel;
                            uvc1_qual_t triallele_inc = ((normsymbol != symbol) ? 
                                    (isSymbolSubstitution(symbol) ? phred_het3al_chance_inc_snp : phred_het3al_chance_inc_indel) : 0);
                            uvc1_qual_t triallele_thr = 3; // 8;
                            const auto new_nlodq_inc = BETWEEN(MIN(binom_b10log10like, powlaw_b10log10like), -triallele_thr, paramset.powlaw_anyvar_base) 
                                    + triallele_inc;
                            if (nlodq_inc > new_nlodq_inc) {
                                nlodq_inc = new_nlodq_inc;
                                argmin_nlodq_symbol = normsymbol;
                                if (paramset.should_add_note) {
                                    fmt.note += std::string("/nlodqDeltaIs/") 
                                        + other_join(std::array<double, 5> {{(double)tAD, (double)tDP, (double)nAD, (double)nDP, (double)nlodq_inc }}, "/") + "#" 
                                        + std::to_string(std::get<0>(nlodq_fmtptr1_fmtptr2_tup)) + "#"
                                        + SYMBOL_TO_DESC_ARR[FIRST(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->VTI)] + "//"
                                        + SYMBOL_TO_DESC_ARR[LAST(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->VTI)] + "//"
                                        + SYMBOL_TO_DESC_ARR[LAST(std::get<2>(nlodq_fmtptr1_fmtptr2_tup)->VTI)] + "//"
                                        + other_join(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->gVQ1, "/") + "//"
                                        + other_join(std::get<2>(nlodq_fmtptr1_fmtptr2_tup)->gVQ1, "/") + "//"
                                        + other_join(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->CONTQ, "/") + "//"
                                        + other_join(std::get<2>(nlodq_fmtptr1_fmtptr2_tup)->CONTQ, "/") + "//";
                                }
                            }
                        }
                        const auto totBDP = SUMPAIR(fmt.BDPb); // (fmt.BDPf[0] + fmt.BDPr[0]);
                        const auto n_norm_alts = (totBDP - (FIRST(fmt.bDPf) + FIRST(fmt.bDPr))) + (LAST(fmt.bDPf) + LAST(fmt.bDPr));
                        //nlodq = MAX(nlodq_singlesite + nlodq_inc, tki.vHGQ + MIN(3, totBDP - n_norm_alts * (int)round(0.5 / paramset.contam_any_mul_frac)));
                        nlodq = MAX(MAX(nlodq_singlesite, germ_phred + nlodq_inc), tki.vHGQ + MIN(3, totBDP - n_norm_alts * (int)round(0.5 / paramset.contam_any_mul_frac)));
                    } else {
                        nlodq = nlodq_singlesample;
                    }
                    fmt.vHGQ = nlodq_singlesample;
                    append_vcf_record(
                            buf_out_string_pass,
                            std::get<0>(tname_tseqlen_tuple).c_str(),
                            refpos,
                            extended_inclu_beg_pos,
                            refstring,
                            region_repeatvec,
                            repeatunit,
                            repeatnum,
                            refsymbol,
                            symbol,
                            fmt,
                            tki,
                            nlodq,
                            argmin_nlodq_symbol,
                            (paramset.should_output_all || is_germline_var_generated),
                            bcf_hdr,
                            baq_offsetarr,
                            paramset,
                            0);
                }
            } // fmt_tki_tup_vec
        } // end of SYMBOL_TYPE_ARR
    } // zerobased_pos
}
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " starts destroying bam records"; }
    umi_strand_readset_uvc_destroy(umi_strand_readset);
    /*
    if (!is_vcf_out_pass_to_stdout) {
        bgzip_string(outstring_pass, buf_out_string_pass);
    } else {
        outstring_pass += buf_out_string_pass;
    }
    */
if (paramset.inferred_is_vcf_generated) {
    uncompressed_vcf_string += buf_out_string_pass;
}
    if (paramset.fam_consensus_out_fastq.size() > 0) {
        for (size_t i = 0; i < NUM_FQLIKE_CON_OUT_FILES; i++) { 
            uncompressed_3fastq_string[i] += fqdata3[i];
        }
    }
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " is done with current task"; }
    return 0;
};


BGZF *bgzip_open_wrap1(const std::string &filename) {
    BGZF *fp = bgzf_open(filename.c_str(), "w");
    if (NULL == fp) {
        LOG(logERROR) << "Unable to open the bgzip file " << filename << " for writinng!";
        exit(-9);
    }
    return fp;
}

int gzip_close_wrap1(BGZF *fp, const std::string &filename) {
    if (NULL != fp) {
        int closeresult = bgzf_close(fp);
        if (closeresult != 0) {
            LOG(logERROR) << "Unable to close the bgzip file " << filename << " for writing!";
        }
        return closeresult;
    } else {
        return 0;
    }
}

static_assert((sizeof(size_t) > 4), "Error: 32-bit architectures are not supported!");
static_assert((sizeof(void*) > 4), "Error: 32-bit memory systems are not supported!");
int 
main(int argc, char **argv) {
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    const char *UMI_STRUCT = getenv("ONE_STEP_UMI_STRUCT");
    const std::string UMI_STRUCT_STRING = ((UMI_STRUCT != NULL && strlen(UMI_STRUCT) > 0) ? std::string(UMI_STRUCT) : std::string(""));
    CommandLineArgs paramset;
    int parsing_result_flag = -1;
    int parsing_result_ret = paramset.initFromArgCV(parsing_result_flag, argc, argv);
    if (paramset.bam_input_fname.compare(OPT_ONLY_PRINT_VCF_HEADER) == 0) {
        std::string header_outstring = generate_vcf_header(
            argc, 
            argv, 
            0, 
            NULL, 
            NULL,
            OPT_ONLY_PRINT_VCF_HEADER,
            paramset);
        std::cout << header_outstring << std::endl;
        return 0;
    }
    if (paramset.bam_input_fname.compare(OPT_ONLY_PRINT_DEBUG_DETAIL) == 0) {
        LOG(logINFO) << "Use about " << SIZE_PER_GENOMIC_POS << " bytes per genomic position";
        return 0;
    }
    if (parsing_result_ret || parsing_result_flag) {
        return parsing_result_ret;
    }
    LOG(logINFO) << "Program " << argv[0] << " version " << VERSION_DETAIL;
    LOG(logINFO) << "<GIT_DIFF_FULL_DISPLAY_MSG>"; 
    LOG(logINFO) << GIT_DIFF_FULL;
    LOG(logINFO) << "</GIT_DIFF_FULL_DISPLAY_MSG>";
    
    std::vector<std::tuple<std::string, uvc1_refgpos_t>> tid_to_tname_tseqlen_tuple_vec;
    samfname_to_tid_to_tname_tseq_tup_vec(tid_to_tname_tseqlen_tuple_vec, paramset.bam_input_fname);
    
    const int nthreads = paramset.max_cpu_num;
    bool is_vcf_out_pass_empty_string = (std::string("") == paramset.vcf_out_pass_fname);
    bool is_vcf_out_pass_to_stdout = (std::string("-") == paramset.vcf_out_pass_fname);
    BGZF *fp_pass = ((!is_vcf_out_pass_empty_string && !is_vcf_out_pass_to_stdout) ? bgzip_open_wrap1(paramset.vcf_out_pass_fname) : NULL);
    std::array<BGZF*, NUM_FQLIKE_CON_OUT_FILES> fastq_fps;
    std::array<std::string, NUM_FQLIKE_CON_OUT_FILES> fastq_filenames = {{ "" }};
    for (size_t i = 0; i < NUM_FQLIKE_CON_OUT_FILES; i++) {
        fastq_filenames[i] = paramset.fam_consensus_out_fastq + FASTQ_LIKE_SUFFIXES[i];
        fastq_fps[i] = ((paramset.fam_consensus_out_fastq.size() > 0) ? (bgzip_open_wrap1(fastq_filenames[i])) :  NULL);
    }
    // Commented out for now due to lack of good documentation for these bgzf APIs. Can investigate later.
    /*
    if (paramset.vcf_output_fname.size() != 0 && paramset.vcf_output_fname != "-") {
        bgzf_index_build_init(fp_allp);
    }
    */
    // bgzf_mt(fp_allp, nthreads, 128);
    // samFile *sam_infile = sam_open(paramset.bam_input_fname.c_str(), "r");
    
    std::ofstream bed_out;
    if (IS_PROVIDED(paramset.bed_out_fname)) {
        bed_out.open(paramset.bed_out_fname, std::ios::out);
    }

#if defined(USE_STDLIB_THREAD)
    const size_t nidxs = nthreads * 2 + 1;
#else
    const size_t nidxs = nthreads;
#endif
    
    bcf_hdr_t *g_bcf_hdr = NULL;
    const char *g_sample = NULL;
    if (IS_PROVIDED(paramset.vcf_tumor_fname)) {
        htsFile *infile = hts_open(paramset.vcf_tumor_fname.c_str(), "r");
        g_bcf_hdr = bcf_hdr_read(infile);
        g_sample = "";
        if (bcf_hdr_nsamples(g_bcf_hdr) > 0) {
            g_sample = g_bcf_hdr->samples[0];
        };
        bcf_close(infile);
    }
    std::vector<hts_idx_t*> sam_idxs(nidxs, NULL);
    std::vector<samFile*> samfiles(nidxs, NULL);
    std::vector<faidx_t*> ref_faidxs(nidxs, NULL);
    std::vector<bcf_srs_t*> srs(nidxs, NULL);
    for (size_t i = 0; i < nidxs; i++) {
        samfiles[i] = sam_open(paramset.bam_input_fname.c_str(), "r");
        if (NULL == samfiles[i]) {
            LOG(logCRITICAL) << "Failed to load BAM file " << paramset.bam_input_fname << " for thread with ID = " << i;
            exit(-3);
        }
        sam_idxs[i] = sam_index_load(samfiles[i], paramset.bam_input_fname.c_str());
        if (NULL == sam_idxs[i]) {
            LOG(logCRITICAL) << "Failed to load BAM index " << paramset.bam_input_fname << " for thread with ID = " << i;
            exit(-4);
        }
        if (paramset.fasta_ref_fname.size() > 0) {
            ref_faidxs[i] = fai_load(paramset.fasta_ref_fname.c_str());
            if (NULL == ref_faidxs[i]) {
                LOG(logCRITICAL) << "Failed to load reference index for file " << paramset.fasta_ref_fname << " for thread with ID = " << i;
                exit(-5);
            }
        }
    }
    
    bam_hdr_t * samheader = sam_hdr_read(samfiles[0]);
    std::string header_outstring = generate_vcf_header(
            argc, 
            argv, 
            samheader->n_targets, 
            samheader->target_name, 
            samheader->target_len,
            g_sample, 
            paramset);
    clearstring<false>(fp_pass, header_outstring, is_vcf_out_pass_to_stdout);

    std::vector<BedLine> bedlines1;
    std::vector<BedLine> bedlines2;
    
    std::map<std::tuple<uvc1_refgpos_t, uvc1_refgpos_t, AlignmentSymbol>, std::vector<TumorKeyInfo>> tid_pos_symb_to_tkis1; 
    std::map<std::tuple<uvc1_refgpos_t, uvc1_refgpos_t, AlignmentSymbol>, std::vector<TumorKeyInfo>> tid_pos_symb_to_tkis2; 
    BedLine prev_bedline_tmp = BedLine(-1, 0, 0, 0, 0);
    SamIter samIter(paramset);
    int64_t n_sam_iters = 0;
    uvc1_flag_t iter_ret_flag;
    int64_t iter_nreads = samIter.iternext(iter_ret_flag, bedlines1, 0);
    
    LOG(logINFO) << "PreProcessed " << iter_nreads << " reads in tier-1-region no " << (n_sam_iters);
    // rescue_variants_from_vcf
    tid_pos_symb_to_tkis1 = rescue_variants_from_vcf(bedlines1, tid_to_tname_tseqlen_tuple_vec, paramset.vcf_tumor_fname, g_bcf_hdr, paramset.is_tumor_format_retrieved);
    LOG(logINFO) << "Rescued/retrieved " << tid_pos_symb_to_tkis1.size() << " variants in tier-1-region no " << (n_sam_iters);
    while (iter_nreads > 0) {
        n_sam_iters++;
        auto read_bam_thread = std::async([&bedlines2, &tid_pos_symb_to_tkis2, &samIter, &iter_nreads, &iter_ret_flag, &n_sam_iters, &paramset, &tid_to_tname_tseqlen_tuple_vec, g_bcf_hdr]() {
            bedlines2.clear();
            iter_nreads = samIter.iternext(iter_ret_flag, bedlines2, 0);
            LOG(logINFO) << "PreProcessed " << iter_nreads << " reads in tier-1-region no " << (n_sam_iters);
            
            tid_pos_symb_to_tkis2 = rescue_variants_from_vcf(bedlines2, tid_to_tname_tseqlen_tuple_vec, paramset.vcf_tumor_fname, g_bcf_hdr, paramset.is_tumor_format_retrieved);
            LOG(logINFO) << "Rescued/retrieved " << tid_pos_symb_to_tkis2.size() << " variants in tier-1-region no " << (n_sam_iters);
        });
        // read_bam_thread.start();
        const auto & bedlines = bedlines1;
        
        const size_t allridx = 0;
        const size_t incvalue = bedlines.size();
        
        size_t nreads = 0;
        size_t npositions = 0;
        for (size_t j = 0; j < incvalue; j++) {
            auto region_idx = allridx + j;
            nreads += bedlines[region_idx].n_reads;
            npositions += bedlines[region_idx].end_pos - bedlines[region_idx].beg_pos; 
        }
        
        assertUVC(incvalue > 0);
        
        // distribute inputs as evenly as possible
#if defined(USE_STDLIB_THREAD)
        const size_t UNDERLOAD_RATIO = 1;
#else
        const size_t UNDERLOAD_RATIO = NUM_WORKING_UNITS_PER_THREAD;
#endif
        uvc1_refgpos_t last_tid = ((bedlines.size() > 0) ? (bedlines[0].tid) : -1);
        uvc1_readnum_big_t curr_nreads = 0;
        uvc1_refgpos_t curr_npositions = 0;
        size_t curr_zerobased_region_idx = 0;
        std::vector<std::pair<size_t, size_t>> beg_end_pair_vec;
        for (size_t j = 0; j < incvalue; j++) {
            auto region_idx = allridx + j;
            const auto curr_tid = bedlines[region_idx].tid;
            curr_nreads += bedlines[region_idx].n_reads;
            curr_npositions +=(bedlines[region_idx].end_pos) - (bedlines[region_idx].beg_pos);
            if ((j > 0) && ((last_tid != curr_tid)
                    || (curr_nreads * nthreads * UNDERLOAD_RATIO > nreads)
                    || (curr_npositions * nthreads * UNDERLOAD_RATIO > npositions))) {
                beg_end_pair_vec.push_back(std::make_pair(curr_zerobased_region_idx, j));
                curr_zerobased_region_idx = j;
                last_tid = curr_tid;
                curr_nreads = 0;
                curr_npositions = 0;
            }
        }
        beg_end_pair_vec.push_back(std::make_pair(curr_zerobased_region_idx, incvalue));

        uvc1_readnum_big_t t1_tot_n_reads = 0;
        uvc1_refgpos_big_t t1_tot_n_bases = 0;

        std::string bedstring = "";
        LOG(logINFO) << "Start-bam-iteration-" << (n_sam_iters-1) << ": will-process-the-following-tier1-regions (all indices are zero-based): ";
        for (size_t t2_idx = 0; t2_idx < beg_end_pair_vec.size(); t2_idx++) {
            size_t beg = beg_end_pair_vec[t2_idx].first;
            size_t end = beg_end_pair_vec[t2_idx].second;
            uvc1_readnum_big_t t2_tot_n_reads = 0;
            uvc1_readnum_big_t t2_tot_n_bases = 0;
            for (size_t t3_idx = beg; t3_idx < end; t3_idx++) {
                const auto & bedline = bedlines[t3_idx];
                const auto n_bases = bedline.end_pos - bedline.beg_pos;
                bedstring += (std::get<0>(tid_to_tname_tseqlen_tuple_vec[bedline.tid ])
                          + "\t" + std::to_string(bedline.beg_pos)
                          + "\t" + std::to_string(bedline.end_pos)
                          + "\tBedLineFlag\t" + std::to_string(bedline.region_flag)
                          + "\tNumberOfReadsInThisInterval\t" + std::to_string(bedline.n_reads)
                          + "\tNumberOfRefBasesInThisInterval\t" + std::to_string(n_bases)
                          + "\tTier1regionIndex\t" + std::to_string(n_sam_iters - 1)
                          + "\tTier2regionIndex\t" + std::to_string(t2_idx)
                          + "\tTier3regionIndex\t" + std::to_string(t3_idx)
                          + "\n");
                t2_tot_n_reads += bedline.n_reads;
                t2_tot_n_bases += n_bases;
            }
            LOG(logINFO) << "End-of-tier-2-region-no-" << t2_idx << " tot_n_reads=" << t2_tot_n_reads << " tot_n_ref_bases=" << t2_tot_n_bases;
            t1_tot_n_reads += t2_tot_n_reads;
            t1_tot_n_bases += t2_tot_n_bases;
        }
        assertUVC(((size_t)t1_tot_n_bases) == npositions || !fprintf(stderr, "%lu ==%lu failed!", ((size_t)t1_tot_n_bases), npositions));
        LOG(logINFO) << "End-of-tier-1-region-no-" << (n_sam_iters-1) << " tot_n_reads=" << t1_tot_n_reads << " tot_n_ref_bases=" << t1_tot_n_bases;
        if (bed_out.is_open()) { 
            bed_out << bedstring; 
        }
        LOG(logINFO) << "The " << (n_sam_iters-1) << "-th tier-1 BED region is as follows " 
                << "(each tier-1 region is parsed as one unit, each tier-2 region is processed by one thread, "
                << "and each tier-3 region is one unit of to call variants within):\n" << bedstring;
        
        LOG(logINFO) << "Start processing the chunks from " << allridx << " to " << allridx + incvalue
                << " which contains approximately " << nreads << " reads and " << npositions << " positions divided into " 
                << beg_end_pair_vec.size() << " sub-chunks";
        
#if defined(USE_STDLIB_THREAD)
        if ( nidxs <= beg_end_pair_vec.size()) {abort();}
        std::vector<std::thread> threads; 
        threads.reserve(beg_end_pair_vec.size());
#endif
        std::vector<BatchArg> batchargs;
        batchargs.reserve(beg_end_pair_vec.size());
        for (size_t beg_end_pair_idx = 0; beg_end_pair_idx < beg_end_pair_vec.size(); beg_end_pair_idx++) {
            struct BatchArg a = {
                    outstring3fastq : (std::array<std::string, NUM_FQLIKE_CON_OUT_FILES> {{ std::string("") }}),
                    outstring_allp : "",
                    outstring_pass : "",
                    thread_id : 0,
                    samfile: NULL,
                    hts_idx : NULL, 
                    ref_faidx : NULL,
                    bcf_hdr : g_bcf_hdr,
                    sr : NULL,
                    
                    prev_bedline: prev_bedline_tmp,
                    bedline: BedLine(-1, 0, 0, 0, 0), // bedlines.at(0),
                    tname_tseqlen_tuple : tid_to_tname_tseqlen_tuple_vec.at(0),
                    regionbatch_ordinal : 0,
                    regionbatch_tot_num : 0,

                    paramset : paramset, 
                    UMI_STRUCT_STRING : UMI_STRUCT_STRING,
                    is_vcf_out_pass_to_stdout : is_vcf_out_pass_to_stdout,
            };
            batchargs.push_back(a);
        }
        size_t beg_end_pair_size = beg_end_pair_vec.size();

#if defined(_OPENMP) && !defined(USE_STDLIB_THREAD)
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
#endif
        for (size_t beg_end_pair_idx = 0; beg_end_pair_idx < beg_end_pair_size; beg_end_pair_idx++) {
            
#if defined(_OPENMP) && !defined(USE_STDLIB_THREAD)
            size_t thread_id = omp_get_thread_num();
#elif defined(USE_STDLIB_THREAD)
            size_t thread_id = beg_end_pair_idx;
#else
            size_t thread_id = 0;
#endif
            auto & batcharg = batchargs[beg_end_pair_idx];
            batcharg.thread_id = thread_id;
            batcharg.samfile = samfiles[thread_id];
            batcharg.hts_idx = sam_idxs[thread_id];
            batcharg.ref_faidx = ref_faidxs[thread_id];
            batcharg.sr = srs[thread_id];

            std::pair<size_t, size_t> beg_end_pair = beg_end_pair_vec[beg_end_pair_idx];
#if defined(USE_STDLIB_THREAD)
            std::thread athread([
                        &batcharg, allridx, beg_end_pair, beg_end_pair_idx, &bedlines, &tid_to_tname_tseqlen_tuple_vec, &tid_pos_symb_to_tkis1
                        ]() {
#endif
                    LOG(logINFO) << "Thread " << batcharg.thread_id << " will process the sub-chunk " << beg_end_pair_idx << " which ranges from " 
                            << beg_end_pair.first << " to " << beg_end_pair.second;
                    assertUVC (beg_end_pair.first < beg_end_pair.second);
                    std::string uncompressed_vcf_string;
                    std::array<std::string, NUM_FQLIKE_CON_OUT_FILES> uncompressed_3fastq_string;
                    for (size_t j = beg_end_pair.first; j < beg_end_pair.second; j++) {
                        batcharg.regionbatch_ordinal = j;
                        batcharg.regionbatch_tot_num = beg_end_pair.second;
                        assertUVC (((size_t)(allridx + j) < bedlines.size())
                                || !fprintf(stderr, "%lu + %lu < %lu failed!\n", allridx, j, bedlines.size()));
                        if (allridx + j > 0) {
                            batcharg.prev_bedline = bedlines.at(allridx + j - 1);
                        }
                        batcharg.bedline = bedlines.at(allridx + j);
                        assertUVC (((size_t)(batcharg.bedline.tid)) < tid_to_tname_tseqlen_tuple_vec.size() 
                                || !fprintf(stderr, "%lu < %lu failed!\n", (size_t)(batcharg.bedline.tid), tid_to_tname_tseqlen_tuple_vec.size()));
                        batcharg.tname_tseqlen_tuple = tid_to_tname_tseqlen_tuple_vec.at((batcharg.bedline.tid));                        
                        process_batch(uncompressed_vcf_string, uncompressed_3fastq_string, batcharg, tid_pos_symb_to_tkis1);
                    }
                    if (batcharg.is_vcf_out_pass_to_stdout) {
                        batcharg.outstring_pass += uncompressed_vcf_string;
                    } else {
                        bgzip_string(batcharg.outstring_pass, uncompressed_vcf_string);
                    }
                    for (size_t i = 0; i < NUM_FQLIKE_CON_OUT_FILES; i++) {
                        bgzip_string(batcharg.outstring3fastq[i], uncompressed_3fastq_string[i]);
                    }
#if defined(USE_STDLIB_THREAD)
            });
            threads.push_back(std::move(athread));
#endif
            
        }
#if defined(USE_STDLIB_THREAD)
        for (auto & t : threads) {
            t.join();
        }
#endif
        for (size_t beg_end_pair_idx = 0; beg_end_pair_idx < beg_end_pair_vec.size(); beg_end_pair_idx++) {
            if (batchargs[beg_end_pair_idx].outstring_pass.size() > 0) {
                // empty string means end of file
                clearstring<true>(fp_pass, batchargs[beg_end_pair_idx].outstring_pass, is_vcf_out_pass_to_stdout);
            }
            for (size_t i = 0; i < fastq_fps.size(); i++) { 
                if (batchargs[beg_end_pair_idx].outstring3fastq[i].size() > 0) {
                    clearstring<true>(fastq_fps[i], batchargs[beg_end_pair_idx].outstring3fastq[i]); 
                }
            } 
        }
        LOG(logINFO) << "Start-waiting-for read_bam_thread with tier-1-region no " << (n_sam_iters);
        // read_bam_thread.join(); // end this iteration
        read_bam_thread.get();
        LOG(logINFO) << "End-waiting-for read_bam_thread with tier-1-region no " << (n_sam_iters);
        for (auto tid_pos_symb_to_tkis1_pair: tid_pos_symb_to_tkis1) {
            for (auto tki : tid_pos_symb_to_tkis1_pair.second) {
                if (NULL != tki.bcf1_record) {
                    bcf_destroy(tki.bcf1_record); 
                }
            }
        }
        prev_bedline_tmp = (bedlines1.size() ? LAST(bedlines1) : prev_bedline_tmp);
        autoswap(bedlines1, bedlines2);
        autoswap(tid_pos_symb_to_tkis1, tid_pos_symb_to_tkis2);
    }
    
    clearstring<true>(fp_pass, std::string(""), is_vcf_out_pass_to_stdout); // write end of file
    for (auto fastq_fp : fastq_fps) { clearstring<true>(fastq_fp, std::string("")); } // empty string means end of file

    bam_hdr_destroy(samheader);
    if (NULL != g_bcf_hdr) {
        bcf_hdr_destroy(g_bcf_hdr);
    }
    for (size_t i = 0; i < nidxs; i++) { 
        if (NULL != srs[i]) {
            bcf_sr_destroy(srs[i]);
        }
        if (NULL != ref_faidxs[i]) { 
            fai_destroy(ref_faidxs[i]); 
        }
        if (NULL != sam_idxs[i]) {
            hts_idx_destroy(sam_idxs[i]);
        }
        if (NULL != samfiles[i]) {
            sam_close(samfiles[i]);
        }
    }
    // bgzf_flush is internally called by bgzf_close
    gzip_close_wrap1(fp_pass, paramset.vcf_out_pass_fname);
    for (size_t i = 0; i < fastq_fps.size(); i++) { gzip_close_wrap1(fastq_fps[i], fastq_filenames[i]); }
    
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
 
    std::cerr << std::fixed << std::setprecision(2) << "CPU time used: "
              << 1.0 * (c_end-c_start) / CLOCKS_PER_SEC << " seconds\n"
              << "Wall clock time passed: "
              << std::chrono::duration<double>(t_end-t_start).count()
              << " seconds\n";
    return 0;
}

