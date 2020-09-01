
#include "main.hpp"

#include "CmdLineArgs.hpp"
#include "common.hpp"
#include "grouping.hpp"
#include "version.h"

#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/synced_bcf_reader.h"

#include <chrono>
#include <ctime>
#include <thread>
#include <tuple>

#if !defined(USE_STDLIB_THREAD)
#include "omp.h"
#endif

const unsigned int G_BLOCK_SIZE = 1000;

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
load_refstring(const faidx_t *ref_faidx, unsigned int tid, unsigned int incbeg, unsigned int excend) {
    assert(incbeg < excend);
    if (NULL == ref_faidx) {
        return std::string(excend - incbeg, 'n');
    }
    const char *tname = faidx_iseq(ref_faidx, tid);
    int regionlen;
    char *fetchedseq = faidx_fetch_seq(ref_faidx, tname, incbeg, excend - 1, &regionlen);
    assert (regionlen == (int)(excend - incbeg) || !fprintf(stderr, "%d == %u - %u failed", regionlen, excend, incbeg));
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

std::map<std::pair<unsigned int, AlignmentSymbol>, std::set<size_t>>
mutform2count4vec_to_simplemut2indices(std::vector<std::pair<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>>> mutform2count4vec) {
    std::map<std::pair<unsigned int, AlignmentSymbol>, std::set<size_t>> simplemut2indices;
    for (size_t i = 0; i < mutform2count4vec.size(); i++) {
        auto counts = mutform2count4vec[i].second; 
        if (counts[0] + counts[1] < 2) { continue; }
        std::basic_string<std::pair<unsigned int, AlignmentSymbol>> mutset = mutform2count4vec[i].first;
        for (auto simplemut : mutset) {
            simplemut2indices.insert(std::make_pair(simplemut, std::set<size_t>()));
            simplemut2indices[simplemut].insert(i);
            
            //for (std::pair<unsigned int, AlignmentSymbol> simplemut : mutforms) {
            //    simplemut2indices.insert(std::make_pair(simplemut, std::set<size_t>()));
            //    simplemut2indices[simplemut].insert(i);
            //}
        }
    }
    return simplemut2indices;
};

/*
std::map<unsigned int, std::set<size_t>> simplemut2indices
simplemut2indices4vec_to_pos2idx( mutform2count4vec) {
    std::map<unsigned int, std::set<size_t>> simplemut2indices;
    for (size_t i = 0; i < mutform2count4vec_bq.size(); i++) {
        for (auto mutform : mutform2count4vec_bq[i].first) {
            for (auto pos_symb : mutform) {
                unsigned int pos = pos_symb.first;
                if (pos_symb.second > 1) {
                    simplemut2indices.insert(std::make_pair(pos, std::set<size_t>()));
                    simplemut2indices[pos].insert(i);
                }
            }
        }
    }
}
*/

int 
bgzip_string(std::string & compressed_outstring, const std::string & uncompressed_outstring) {
    char *compressed_outstr = (char*)malloc(uncompressed_outstring.size() * sizeof(char));
    if (NULL == compressed_outstr) {
        fprintf(stderr, "The library function malloc failed at line %d in file %s !\n", __LINE__, __FILE__);
        exit(-1);
    }
    size_t compressed_capacity = uncompressed_outstring.size();
    size_t compressed_totlen = 0;
    size_t uncompress_totlen = 0;
    do {
        if (compressed_totlen + BGZF_BLOCK_SIZE >= compressed_capacity) {
            char *compressed_outstr_tmp = (char*)realloc(compressed_outstr, (compressed_capacity * 2 + BGZF_BLOCK_SIZE) * sizeof(char));
            if (NULL == compressed_outstr_tmp) {
                fprintf(stderr, "The library function realloc failed at line %d in file %s !\n", __LINE__, __FILE__);
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
        // fprintf(stderr, "Compressed %d out of %d raw data into %d\n", uncompressed_totlen, uncompressed_outstring.size(), compressed_totlen);
    } while (uncompress_totlen < uncompressed_outstring.size());
    compressed_outstring += std::string(compressed_outstr, compressed_totlen);
    free(compressed_outstr);
    return 0;
}

struct BatchArg {
    std::string outstring_allp;
    std::string outstring_pass;
    // VcStats vc_stats;
    unsigned int thread_id;
    hts_idx_t *hts_idx;
    faidx_t *ref_faidx;
    bcf_hdr_t *bcf_hdr;
    bcf_srs_t *sr;
    
    std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int> tid_beg_end_e2e_tuple;
    std::tuple<std::string, unsigned int> tname_tseqlen_tuple;
    unsigned int region_ordinal; // deprecated
    unsigned int region_tot_num; // deprecated
    unsigned int regionbatch_ordinal;
    unsigned int regionbatch_tot_num;

    const CommandLineArgs paramset;
    const std::string UMI_STRUCT_STRING;
    const bool is_vcf_out_pass_to_stdout;
    const bool is_vcf_out_empty_string;
};

unsigned int
gen_fq_tsum_depths(const auto & fq_tsum_depth, unsigned int refpos) {
    std::array<unsigned int, 3> fq_tsum_depths = {{0, 0, 0}};
    for (unsigned int strand = 0; strand < 2; strand++) {
        fq_tsum_depths[0] += fq_tsum_depth.at(strand).getByPos(refpos+0).sumBySymbolType(LINK_SYMBOL);
        fq_tsum_depths[1] += fq_tsum_depth.at(strand).getByPos(refpos+0).sumBySymbolType(BASE_SYMBOL);
        fq_tsum_depths[2] += fq_tsum_depth.at(strand).getByPos(refpos+1).sumBySymbolType(LINK_SYMBOL);
    }
    return MIN(MIN(fq_tsum_depths[0], fq_tsum_depths[2]), fq_tsum_depths[1]);
}

std::vector<unsigned int>
gen_dp100(const auto & fq_tsum_depth, unsigned int inclu_beg, unsigned int exclu_end) {
    assert (inclu_beg <= exclu_end);
    std::vector<unsigned int> dp100;
    dp100.reserve(exclu_end - inclu_beg);
    //for (unsigned int rpos2 = refpos; rpos2 < MIN(refpos+100, rpos_exclu_end); rpos2++) {
    for (unsigned int rpos2 = inclu_beg; rpos2 < exclu_end; rpos2++) {
        unsigned int fq_tsum_depth_g = gen_fq_tsum_depths(fq_tsum_depth, rpos2);
        dp100.push_back(fq_tsum_depth_g);
    }
    return dp100;
}

std::string 
dp100_to_string(const std::vector<unsigned int> & bg1dp100, const std::vector<unsigned int> & bg2dp100, 
        const std::string & chromosome, unsigned int refpos, bool is_rev) {
    assert (bg1dp100.size() == bg2dp100.size());
    std::string ret = chromosome + "\t" + std::to_string(refpos) + 
            "\t.\tN\t<DP100" + (is_rev ? "RV" : "FW") + ">\t.\t.\t.\tGT:bgNPOS:bg1DPS:bg2DPS\t.:" + std::to_string(bg1dp100.size()) + ":";
    for (auto dp : bg1dp100) {
        ret += std::to_string(dp) + ",";
    }
    ret += "-1:";
    for (auto dp : bg2dp100) {
        ret += std::to_string(dp) + ",";
    }
    ret += "-1\n";
    return ret;
}

std::string
genomic_reg_info_to_string(const std::string & chromosome, 
        unsigned int incluBeg, SymbolType stypeBeg,
        unsigned int incluEnd, SymbolType stypeEnd,
        const unsigned int gbDPmin, const unsigned int gcDPmin,
        const std::string &gfGTmm2, const unsigned int gfGQmin,
        const std::string & refstring, unsigned int refstring_offset) {
    unsigned int begpos = incluBeg; // = (stypeBeg == BASE_SYMBOL ? (incluBeg+1) : incluBeg);
    unsigned int endpos = incluEnd; // = (stypeEnd == BASE_SYMBOL ? (excluEnd+1) : (excluEnd));
    unsigned int refstring_idx = begpos - refstring_offset;
    const std::string begchar = (refstring_idx > 0 ? refstring.substr(refstring_idx - 1, 1) : "n");
    std::string ret = chromosome + "\t" + std::to_string(begpos)
            + "\t.\t" + begchar + "\t<NON_REF" + ">\t.\t.\t.\tGT:GQ:gbDP:gcDP:gSTS:gBEG:gEND\t"
            +               (gfGTmm2) + ":"
            + std::to_string(gfGQmin) + ":"
            + std::to_string(gbDPmin) + ":"
            + std::to_string(gcDPmin) + ":"
            + std::to_string(stypeBeg) + "," + std::to_string(stypeEnd) + ":"
            + std::to_string(begpos) + ":"
            + std::to_string(endpos) + "\n";
    return ret;
}

const bool
is_sig_higher_than(auto a, auto b, unsigned int mfact, unsigned int afact) {
    return (a * 100 > b * (100 + mfact)) && (a > b + afact);
}

const bool
is_sig_out(auto a, auto minval, auto maxval, unsigned int mfact, unsigned int afact) {
    return is_sig_higher_than(a, minval, mfact, afact) || is_sig_higher_than(maxval, a, mfact, afact);
}

std::string 
als_to_string(const char *const* const allele, unsigned int n_allele) {
    std::string ret;
    ret.reserve(n_allele*2);
    for (unsigned int i = 0; i < n_allele; i++) {
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

const std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, std::vector<TumorKeyInfo>>
rescue_variants_from_vcf(const auto & tid_beg_end_e2e_vec, const auto & tid_to_tname_tlen_tuple_vec, const std::string & vcf_tumor_fname, const auto *bcf_hdr, 
        const bool is_tumor_format_retrieved) {
    std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, std::vector<TumorKeyInfo>> ret;
    if (NOT_PROVIDED == vcf_tumor_fname) {
        return ret;
    }
    std::string regionstring;
    for (const auto & tid_beg_end_e2e : tid_beg_end_e2e_vec) {
        const unsigned int tid = std::get<0>(tid_beg_end_e2e);
        const unsigned int rpos_inclu_beg = std::get<1>(tid_beg_end_e2e);
        const unsigned int rpos_exclu_end = std::get<2>(tid_beg_end_e2e);
        
        const auto & tname_tseqlen_tuple = tid_to_tname_tlen_tuple_vec[tid];
        if (0 < regionstring.size()) { regionstring += std::string(","); }
        regionstring += (std::get<0>(tname_tseqlen_tuple) + ":" + std::to_string(rpos_inclu_beg+1) + "-" + std::to_string(rpos_exclu_end+1-1)); 
    }
    
    { 
        LOG(logINFO) << "Region is " << regionstring;
    }
    if (regionstring.size() == 0) {
        return ret;
    }
    bcf_srs_t *const sr = bcf_sr_init();
    if (NULL == sr) {
        LOG(logCRITICAL) << "Failed to initialize bcf sr";
        exit(-6);
    }
    
    // int retflag = // bcf_sr_set_opt(srs[i], BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF); 
    bcf_sr_set_regions(sr, regionstring.c_str(), false);
    /*int sr_set_opt_retval =*/ 
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    int sr_add_reader_retval = bcf_sr_add_reader(sr, vcf_tumor_fname.c_str());
    if (sr_add_reader_retval != 1) {
        LOG(logCRITICAL) << "Failed to synchronize-read the tumor vcf " << vcf_tumor_fname << " with return code " << sr_add_reader_retval;
        exit(-7);
    }
    
    int valsize = 0; 
    int ndst_val = 0;
    char *bcfstring = NULL;
    float *bcffloats = NULL;
    int32_t *bcfints = NULL;
    
    while (bcf_sr_next_line(sr)) {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        bcf_unpack(line, BCF_UN_ALL);
        
        // skip over all symbolic alleles
        bool should_continue = false;
        for (unsigned int i = 1; i < line->n_allele; i++) {
            if ('<' == line->d.allele[i][0]) {
                should_continue = true;
            }
        }
        if (should_continue) { continue; }
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "VTI", &bcfints, &ndst_val);
        if (valsize <= 0) { continue; }
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for VTI and line %d!\n", ndst_val, valsize, line->pos));
        assert((2 == line->n_allele) || !fprintf(stderr, "Bcf line %d has %d alleles!\n", line->pos, line->n_allele));
        // std::string desc(bcfstring);
        // LOG(logINFO) << "Trying to retrieve the symbol " << desc << " at pos " << line->pos << " valsize = " << valsize << " ndst_val = " << ndst_val;
        const AlignmentSymbol symbol = AlignmentSymbol(bcfints[1]); //DESC_TO_SYMBOL_MAP.at(desc);
        
        auto symbolpos = (isSymbolSubstitution(symbol) ? (line->pos) : (line->pos + 1));
        /*
        if (symbolpos < extended_inclu_beg_pos && isSymbolDel(symbol)) {
            continue;
        }
        auto intpos = symbolpos - extended_inclu_beg_pos;
        if (!(intpos < extended_posidx_to_symbol_to_tkinfo.size())) {
            LOG(logINFO) << "Warning!!! Trying to retrieve the symbol " << desc << " at pos " << line->pos << " valsize = " << valsize << " ndst_val = " << ndst_val;
            fprintf(stderr, "%d < %d failed with regionstring %s!\n", intpos, extended_posidx_to_symbol_to_tkinfo.size(), regionstring.c_str());
            abort();
        }
        */
        TumorKeyInfo tki;
        tki.VTI = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "BDPf", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for BDPf and line %d!\n", ndst_val, valsize, line->pos));
        tki.BDP = bcfints[0];
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "BDPr", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for BDPr and line %d!\n", ndst_val, valsize, line->pos));
        tki.BDP += bcfints[0];

        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "bDPf", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for bDPf and line %d!\n", ndst_val, valsize, line->pos));
        tki.bDP = bcfints[1];
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "bDPr", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for bDPr and line %d!\n", ndst_val, valsize, line->pos));
        tki.bDP += bcfints[1];

        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "CDP1f", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for CDP1f and line %d!\n", ndst_val, valsize, line->pos));
        tki.CDP1 = bcfints[0];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "CDP1r", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for CDP1r and line %d!\n", ndst_val, valsize, line->pos));
        tki.CDP1 += bcfints[0];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cDP1f", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cDP1f and line %d!\n", ndst_val, valsize, line->pos));
        tki.cDP1 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cDP1r", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cDP1r and line %d!\n", ndst_val, valsize, line->pos));
        tki.cDP1 += bcfints[1];
        

        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "CDP1x", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for CDP1x and line %d!\n", ndst_val, valsize, line->pos));
        tki.CDP1x = bcfints[0];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cDP1x", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cDP1x and line %d!\n", ndst_val, valsize, line->pos));
        tki.cDP1x = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cVQ1", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cVQ1 and line %d!\n", ndst_val, valsize, line->pos));
        tki.cVQ1 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cPCQ1", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cPCQ1 and line %d!\n", ndst_val, valsize, line->pos));
        tki.cPCQ1 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "CDP2x", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for CDP2x and line %d!\n", ndst_val, valsize, line->pos));
        tki.CDP2x = bcfints[0];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cDP2x", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cDP2x and line %d!\n", ndst_val, valsize, line->pos));
        tki.cDP2x = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cVQ2", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cVQ2 and line %d!\n", ndst_val, valsize, line->pos));
        tki.cVQ2 = bcfints[1];
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line, "cPCQ2", &bcfints, &ndst_val);
        assert((2 == ndst_val && 2 == valsize) || !fprintf(stderr, "2 == %d && 2 == %d failed for cPCQ2 and line %d!\n", ndst_val, valsize, line->pos));
        tki.cPCQ2 = bcfints[1];

        tki.pos = line->pos;
        tki.ref_alt = als_to_string(line->d.allele, line->n_allele);
        if (is_tumor_format_retrieved) {
            tki.bcf1_record = bcf_dup(line);
        }
        const auto retkey = std::make_tuple(line->rid, symbolpos, symbol);            
        ret.insert(std::make_pair(retkey, std::vector<TumorKeyInfo>()));
        ret[retkey].push_back(tki);
        // ret.insert(std::make_pair(retkey, tki));
    }
    xfree(bcfstring);
    xfree(bcfints);
    xfree(bcffloats);
    bcf_sr_destroy(sr);
    return ret;
}

int 
process_batch(BatchArg & arg, const auto & tid_pos_symb_to_tkis) {
    
    std::string & outstring_allp = arg.outstring_allp;
    std::string & outstring_pass = arg.outstring_pass;
    const hts_idx_t *const hts_idx = arg.hts_idx;
    const faidx_t *const ref_faidx = arg.ref_faidx;
    
    const bcf_hdr_t *const bcf_hdr = arg.bcf_hdr;
    const CommandLineArgs & paramset = arg.paramset;
    const std::string UMI_STRUCT_STRING = arg.UMI_STRUCT_STRING;
    const std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int> tid_beg_end_e2e_tuple = arg.tid_beg_end_e2e_tuple;
    const std::tuple<std::string, unsigned int> tname_tseqlen_tuple = arg.tname_tseqlen_tuple;
    const unsigned int regionbatch_ordinal = arg.regionbatch_ordinal;
    const unsigned int regionbatch_tot_num = arg.regionbatch_tot_num;
    const unsigned int thread_id = arg.thread_id;
    //  const bool should_let_all_pass = paramset.should_let_all_pass;
    const bool is_vcf_out_pass_to_stdout = arg.is_vcf_out_pass_to_stdout;
    
    bool is_loginfo_enabled = (ispowerof2(regionbatch_ordinal + 1) || ispowerof2(regionbatch_tot_num - regionbatch_ordinal));
    std::string raw_out_string;
    std::string raw_out_string_pass;
    // faidx_t *ref_faidx = (fasta_ref_fname.size() > 0 ? fai_load(fasta_ref_fname.c_str()) : NULL);
    
    auto tid = std::get<0>(tid_beg_end_e2e_tuple);
    auto incluBegPosition = std::get<1>(tid_beg_end_e2e_tuple);
    auto excluEndPosition = std::get<2>(tid_beg_end_e2e_tuple);
    auto end2end = std::get<3>(tid_beg_end_e2e_tuple);
    // auto nreads = std::get<4>(tid_beg_end_e2e_tuple);
    
    std::map<uint64_t, std::pair<std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>, int>> umi_to_strand_to_reads;
    unsigned int bam_inclu_beg_pos, bam_exclu_end_pos; 
    std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> umi_strand_readset;

    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts bamfname_to_strand_to_familyuid_to_reads with pair_end_merge = " << paramset.pair_end_merge; }
    std::array<unsigned int, 3> passed_pcrpassed_umipassed = bamfname_to_strand_to_familyuid_to_reads(
            umi_to_strand_to_reads,
            bam_inclu_beg_pos,
            bam_exclu_end_pos,
            paramset.bam_input_fname,
            tid,
            incluBegPosition, 
            excluEndPosition,
            end2end,
            paramset.min_mapqual,
            paramset.min_aln_len,
            regionbatch_ordinal,
            regionbatch_tot_num,
            UMI_STRUCT_STRING,
            hts_idx,
            ASSAY_TYPE_CAPTURE != paramset.assay_type,
            PAIR_END_MERGE_NO != paramset.pair_end_merge,
            paramset.disable_duplex,
            thread_id,
            paramset.dedup_center_mult,
            paramset.dedup_amplicon_count_to_surrcount_ratio,
            paramset.dedup_amplicon_end2end_ratio,
            paramset.always_log,
            (SEQUENCING_PLATFORM_IONTORRENT == paramset.sequencing_platform),
            paramset.dedup_flag,
            0);
    unsigned int num_passed_reads = passed_pcrpassed_umipassed[0];
    unsigned int num_pcrpassed_reads = passed_pcrpassed_umipassed[1];
    unsigned int num_umipassed_reads = passed_pcrpassed_umipassed[2];
    bool is_umi_barcoded = (num_umipassed_reads * 2 > num_passed_reads);
    bool is_by_capture = ((num_pcrpassed_reads) * 2 <= num_passed_reads);
    const AssayType inferred_assay_type = ((ASSAY_TYPE_AUTO == paramset.assay_type) ? (is_by_capture ? ASSAY_TYPE_CAPTURE : ASSAY_TYPE_AMPLICON) : (paramset.assay_type));
    
    if (0 == num_passed_reads) { return -1; };
    unsigned int minABQ_snv = ((ASSAY_TYPE_AMPLICON == inferred_assay_type) ? paramset.minABQ_pcr_snv : paramset.minABQ_cap_snv);
    unsigned int minABQ_indel = ((ASSAY_TYPE_AMPLICON == inferred_assay_type) ? paramset.minABQ_pcr_indel : paramset.minABQ_cap_indel);
    
    const unsigned int rpos_inclu_beg = MAX(incluBegPosition, bam_inclu_beg_pos);
    const unsigned int rpos_exclu_end = MIN(excluEndPosition, bam_exclu_end_pos); 
    const unsigned int extended_inclu_beg_pos = (unsigned int)MAX(0, MIN(((int)incluBegPosition) - 100, (int)bam_inclu_beg_pos));
    const unsigned int extended_exclu_end_pos = MIN(std::get<1>(tname_tseqlen_tuple), MAX(excluEndPosition + 100, bam_exclu_end_pos));
    
    //std::vector<std::array<TumorKeyInfo, NUM_ALIGNMENT_SYMBOLS>> extended_posidx_to_symbol_to_tkinfo(extended_exclu_end_pos - extended_inclu_beg_pos + 1);
    const auto tkis_beg = tid_pos_symb_to_tkis.lower_bound(std::make_tuple(tid, extended_inclu_beg_pos    , AlignmentSymbol(0)));
    const auto tkis_end = tid_pos_symb_to_tkis.upper_bound(std::make_tuple(tid, extended_exclu_end_pos + 1, AlignmentSymbol(0)));
    std::vector<bool> extended_posidx_to_is_rescued(extended_exclu_end_pos - extended_inclu_beg_pos + 1, false);
    unsigned int num_rescued = 0;
    for (auto tkis_it = tkis_beg; tkis_it != tkis_end; tkis_it++) {
        auto symbolpos = std::get<1>(tkis_it->first);
        extended_posidx_to_is_rescued[symbolpos - extended_inclu_beg_pos] = true;
        num_rescued++;
        if (is_loginfo_enabled) {
            // NOTE: the true positive short del at 22:17946835 in NA12878-NA24385 mixture is overwhelmed by the false positve long del spanning the true positive short del.
            // However, so far manual check with limited experience cannot confirm that the true positive is indeed a true positive.
            // TODO: have to see more examples of this case and adjust code accordingly.
            LOG(logDEBUG4) << "Thread " << thread_id << " iterated over symbolpos " << symbolpos << " and symbol " << std::get<2>(tkis_it->first) << " as a rescued var";
        }
    }
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " deals with " << num_rescued << " tumor-sample variants in region " << extended_inclu_beg_pos << " to " << extended_exclu_end_pos + 1 ;}
    
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts converting umi_to_strand_to_reads with is_by_capture = " << is_by_capture << "  " ;}
    fill_strand_umi_readset_with_strand_to_umi_to_reads(umi_strand_readset, umi_to_strand_to_reads, paramset.baq_per_aligned_base);
    
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts constructing symbolToCountCoverageSet12 with " << extended_inclu_beg_pos << (" , ") << extended_exclu_end_pos; }
    // + 1 accounts for insertion at the end of the region, this should happen rarely for only re-aligned reads at around once per one billion base pairs
    Symbol2CountCoverageSet symbolToCountCoverageSet12(tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1); 
    if (is_loginfo_enabled) { LOG(logINFO)<< "Thread " << thread_id << " starts updateByRegion3Aln with " << umi_strand_readset.size() << " families"; }
    std::string refstring = load_refstring(ref_faidx, tid, extended_inclu_beg_pos, extended_exclu_end_pos);
    std::vector<RegionalTandemRepeat> region_repeatvec = refstring2repeatvec(refstring);
    // repeatvec_LOG(region_repeatvec, extended_inclu_beg_pos);

    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> adjcount_x_rpos_x_misma_vec;
    std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> mutform2count4map_bq;
    std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> mutform2count4map_fq;
    const PhredMutationTable sscs_mut_table(
                paramset.phred_max_sscs_transition_CG_TA, 
                paramset.phred_max_sscs_transition_TA_CG, 
                paramset.phred_max_sscs_transversion_any,
                paramset.phred_max_sscs_indel_open,
                paramset.phred_max_sscs_indel_ext);
    symbolToCountCoverageSet12.updateByRegion3Aln(
            mutform2count4map_bq, 
            mutform2count4map_fq,
            umi_strand_readset, 
            refstring,
            region_repeatvec,
            paramset.bq_phred_added_misma, 
            paramset.bq_phred_added_indel, 
            paramset.should_add_note, 
            paramset.phred_max_frag_indel_basemax, 
            sscs_mut_table,
            !paramset.disable_dup_read_merge, 
            SEQUENCING_PLATFORM_IONTORRENT == paramset.sequencing_platform,
            ((ASSAY_TYPE_AMPLICON == inferred_assay_type) ? paramset.primerlen : 0),
            paramset.phred_frac_indel_error_before_barcode_labeling,
            paramset.baq_per_aligned_base,
            0);
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts analyzing phasing info"; }
    auto mutform2count4vec_bq = map2vector(mutform2count4map_bq);
    auto simplemut2indices_bq = mutform2count4vec_to_simplemut2indices(mutform2count4vec_bq);
    auto mutform2count4vec_fq = map2vector(mutform2count4map_fq);
    auto simplemut2indices_fq = mutform2count4vec_to_simplemut2indices(mutform2count4vec_fq);
    
#if 0
    LOG(logINFO) << "Thread " << thread_id << " is dealing with " << mutform2count4vec_bq.size() << " phasing units.";
    if (mutform2count4map_bq.size() < 20) {
        for (const auto & mutform2count : mutform2count4vec_bq) {
            LOG(logINFO) << "Muatation " << (mutform2count.second[0]) << "," 
                    << (mutform2count.second[1]);
            for (const auto & pos_symb : mutform2count.first) {
                LOG(logINFO) << (pos_symb.first) << "," 
                        << (pos_symb.second);
            }
        }
    }
#endif
    
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " starts generating block gzipped vcf"; }
    
    std::string buf_out_string_pass;
    const std::set<size_t> empty_size_t_set;
    
    // This code can result in false negative InDels, but it looks like these false negative InDels can be avoided with proper (but potentially very complex) normalization
    // Therefore, this code and all other code related to this code are commented out for now.
    /*
    std::array<CoveredRegion<uint32_t>, 2> bq_indel_adjmax_depths = {
        CoveredRegion<uint32_t>(tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1),
        G
        CoveredRegion<uint32_t>(tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1)
    };
    for (unsigned int refpos = extended_inclu_beg_pos; refpos <= extended_exclu_end_pos; refpos++) {
        auto totins = bq_ins_2bdepths[0].getByPos(refpos) + bq_ins_2bdepths[1].getByPos(refpos);
        if (totins > 0) {
            for (unsigned int p2 = MAX(extended_inclu_beg_pos + 20, refpos) - 20; p2 <= MIN(refpos + 20, extended_exclu_end_pos); p2++) {
                bq_indel_adjmax_depths[0].getRefByPos(p2) = MAX(bq_indel_adjmax_depths[0].getByPos(p2), totins);
            }
        }
        auto totdel = bq_del_2bdepths[0].getByPos(refpos) + bq_del_2bdepths[1].getByPos(refpos);
        if (totdel > 0) {
            for (unsigned int p2 = MAX(extended_inclu_beg_pos + 20, refpos) - 20; p2 <= MIN(refpos + 20, extended_exclu_end_pos); p2++) {
                bq_indel_adjmax_depths[1].getRefByPos(p2) = MAX(bq_indel_adjmax_depths[1].getByPos(p2), totdel);
            }
        }
    }
    */
    
    for (unsigned int refpos = rpos_inclu_beg; refpos <= rpos_exclu_end; refpos++) {
        std::string repeatunit;
        unsigned int repeatnum = 0;
        
        unsigned int rridx = refpos - extended_inclu_beg_pos;
        indelpos_to_context(repeatunit, repeatnum, refstring, rridx);
        const std::array<SymbolType, 2> stype_to_immediate_prev = {{LINK_SYMBOL, BASE_SYMBOL}};
        
        std::array<unsigned int, 2> stype2refqual = {{0, 0}};
        for (SymbolType symbolType : SYMBOL_TYPES_IN_VCF_ORDER) {
            TumorKeyInfo THE_DUMMY_TUMOR_KEY_INFO;
            bcfrec::BcfFormat init_fmt;
            const AlignmentSymbol refsymbol = (LINK_SYMBOL == symbolType ? LINK_M : (
                    refstring.size() == (refpos - extended_inclu_beg_pos) ? BASE_NN : CHAR_TO_SYMBOL.data.at(refstring.at(refpos - extended_inclu_beg_pos))));
            std::array<unsigned int, 2> bDPcDP = BcfFormat_symboltype_init(
                    init_fmt, 
                    symbolToCountCoverageSet12, 
                    refpos, 
                    symbolType, 
                    !paramset.disable_dup_read_merge, 
                    refsymbol, 
                    0);
            unsigned int var_bdepth = bDPcDP[0];
            unsigned int var_cdepth = bDPcDP[1];
            
            unsigned int ins_bdepth = 0;
            unsigned int del_bdepth = 0;
            unsigned int ins_cdepth = 0;
            unsigned int del_cdepth = 0;
            
            unsigned int ins1_bdepth = 0;
            unsigned int del1_bdepth = 0;
            unsigned int ins1_cdepth = 0;
            unsigned int del1_cdepth = 0;
            
            /*
            std::array<unsigned int, 2> refsymbol_bDPcDP = BcfFormat_symbol_init(
                    init_fmt,
                    symbolToCountCoverageSet12,
                    refpos,
                    refsymbol,
                    0,
                    mutform2count4vec_bq,             
                    empty_size_t_set,
                    mutform2count4vec_fq,
                    empty_size_t_set,
                    0);
            BcfFormat_symbol_calc_qual(
                    init_fmt,
                    refsymbol_bDPcDP[0],
                    refsymbol_bDPcDP[1],
                    "",
                    repeatunit,
                    repeatnum,
                    0,
                    0);
            */
            if (rpos_exclu_end != refpos) {
                std::vector<std::tuple<bcfrec::BcfFormat, TumorKeyInfo>> fmt_tki_tup_vec;
                
                for (AlignmentSymbol symbol : SYMBOL_TYPE_TO_SYMBOLS[symbolType]) 
                {
                    // unsigned int phred_max_sscs_all = sscs_mut_table.toPhredErrRate(refsymbol, symbol);
                    const bool is_pos_rescued = (NOT_PROVIDED != paramset.vcf_tumor_fname && (extended_posidx_to_is_rescued[refpos - extended_inclu_beg_pos]));
                    const bool is_var_rescued = (is_pos_rescued && (tid_pos_symb_to_tkis.end() != tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, symbol)))); 
                    const auto bdepth = 
                            symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[0].getByPos(refpos)[symbol][FRAG_bDP]
                          + symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[1].getByPos(refpos)[symbol][FRAG_bDP];
                    const auto cdepth = 
                            symbolToCountCoverageSet12.symbol_to_fam_format_depth_sets_2strand[0].getByPos(refpos)[symbol][FAM_cDP1]
                          + symbolToCountCoverageSet12.symbol_to_fam_format_depth_sets_2strand[1].getByPos(refpos)[symbol][FAM_cDP1];
                    if (refsymbol == symbol) {
                        var_bdepth = bDPcDP[0] - bdepth;
                        var_cdepth = bDPcDP[1] - cdepth;
                    }
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
                    if ((NOT_PROVIDED == paramset.vcf_tumor_fname)
                            && (bdepth < paramset.min_altdp_thres) 
                            && (refsymbol != symbol)
                            && (!paramset.should_output_all)) {
                        continue;
                    }
                    
                    if (NOT_PROVIDED != paramset.vcf_tumor_fname
                            && (!is_pos_rescued
                                //  || ((!is_var_rescued) && refsymbol != symbol)
                                )) {
                        continue;
                    }
                    const auto simplemut = std::make_pair(refpos, symbol);
                    auto indices_bq = (simplemut2indices_bq.find(simplemut) != simplemut2indices_bq.end() ? simplemut2indices_bq[simplemut] : empty_size_t_set); 
                    auto indices_fq = (simplemut2indices_fq.find(simplemut) != simplemut2indices_fq.end() ? simplemut2indices_fq[simplemut] : empty_size_t_set);
                    std::vector<TumorKeyInfo> tkis;
                    if (is_var_rescued) {
                        tkis = tid_pos_symb_to_tkis.find(std::make_tuple(tid, refpos, symbol))->second;
                    }
                    
                    auto fmt = init_fmt;
                    int tkiidx = 0;
                    std::vector<std::tuple<unsigned int, unsigned int, std::string, TumorKeyInfo>> bcad0a_indelstring_tki_vec;
                    if (isSymbolIns(symbol) || isSymbolDel(symbol)) {
                        for (size_t strand = 0; strand < 2; strand++) {
                            if (0 < symbolToCountCoverageSet12.symbol_to_frag_format_depth_sets[strand].getByPos(refpos)[symbol][FRAG_bDP]) {
                                fill_by_indel_info(
                                        fmt, 
                                        symbolToCountCoverageSet12,
                                        strand,
                                        refpos,
                                        symbol,
                                        refstring,
                                        repeatunit,
                                        repeatnum
                                );
                            }
                        }
                        if (is_var_rescued) {
                            for (const auto & tki : tkis) {
                                unsigned int tabpos = tki.ref_alt.find("\t");
                                const auto & vcfref = tki.ref_alt.substr(0, tabpos);
                                const auto & vcfalt = tki.ref_alt.substr(tabpos + 1);
                                std::string indelstring;
                                if (vcfref.size() > vcfalt.size()) {
                                    indelstring = vcfref.substr(vcfalt.size());
                                } else {
                                    assert (vcfref.size() < vcfalt.size());
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
                                    true);
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
                        std::array<unsigned int, 2> altsymbol_bDPcDP = BcfFormat_symbol_init(
                                fmt,
                                symbolToCountCoverageSet12,
                                refpos,
                                symbol,
                                mutform2count4vec_bq,     
                                indices_bq,
                                mutform2count4vec_fq, 
                                indices_fq,
                                std::get<0>(bcad0a_indelstring_tki),
                                std::get<1>(bcad0a_indelstring_tki),
                                indelstring,
                                (isSymbolSubstitution(symbol) ? minABQ_snv : minABQ_indel),
                                0);
                        BcfFormat_symbol_calc_DPv(
                                fmt,
                                refpos,
                                region_repeatvec.at(MAX(refpos - extended_inclu_beg_pos, 3) - 3),
                                region_repeatvec.at(MIN(refpos - extended_inclu_beg_pos + 3, region_repeatvec.size() - 3)),
                                ((is_var_rescued && (tki.VTI == LAST(fmt.VTI))) ? ((double)(tki.cDP1x + 1) / (double)(tki.CDP1x + 2)) : -1.0), /* tpfa */
                                (ASSAY_TYPE_AMPLICON == inferred_assay_type),
                                0);
                        fmt_tki_tup_vec.push_back(std::make_tuple(fmt, std::get<3>(bcad0a_indelstring_tki)));
                    }
                }
                
                if (rpos_exclu_end != refpos && fmt_tki_tup_vec.size() > 0) {
                    BcfFormat_symbol_sum_DPv(fmt_tki_tup_vec);
                    
                    //int max_cVQ1 = -999;
                    //int max_cVQ2 = -999;
                    //AlignmentSymbol argmax_cVQ1 = END_ALIGNMENT_SYMBOLS;
                    //AlignmentSymbol argmax_cVQ2 = END_ALIGNMENT_SYMBOLS;
                    //std::string strmax_cVQ1 = "";
                    //std::string strmax_cVQ2 = "";
                    
                    std::vector<std::tuple<int, int, int, AlignmentSymbol, std::string>> maxVQ_VQ1_VQ2_symbol_indelstr_tup_vec;
                    for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                        AlignmentSymbol symbol = AlignmentSymbol(LAST(std::get<0>(fmt_tki_tup).VTI));
                        unsigned int phred_max_sscs_all = sscs_mut_table.toPhredErrRate(refsymbol, symbol);
                        auto & fmt = std::get<0>(fmt_tki_tup);
                        const auto & tki = std::get<1>(fmt_tki_tup);
                        BcfFormat_symbol_calc_qual(
                                fmt,
                                var_bdepth,
                                var_cdepth,
                                
                                ins_bdepth,
                                ins_cdepth,
                                del_bdepth,
                                del_cdepth,
                                
                                ins1_bdepth,
                                ins1_cdepth,
                                del1_bdepth,
                                del1_cdepth,
                                
                                repeatunit,
                                repeatnum,
                                90,
                                3,
                                phred_max_sscs_all,
                                (phred_max_sscs_all - (isSymbolSubstitution(symbol) ? paramset.phred_pow_sscs_origin : paramset.phred_pow_sscs_indel_origin)),
                                paramset.phred_varcall_err_per_map_err_per_base,
                                NOT_PROVIDED != paramset.vcf_tumor_fname,
                                // (is_mut_transition(refsymbol, AlignmentSymbol(LAST(std::get<0>(fmt_tki_tup).VTI))) ? 2 : 0),
                                region_repeatvec.at(MAX(refpos - extended_inclu_beg_pos, 3) - 3),
                                region_repeatvec.at(MIN(refpos - extended_inclu_beg_pos + 3, region_repeatvec.size() - 3)),
                                refpos,
                                tid,
                                ((NOT_PROVIDED != paramset.vcf_tumor_fname && (tki.VTI == LAST(fmt.VTI))) ? 
                                    ((double)(tki.cDP1x + 1) / (double)(tki.CDP1x + 2)) : -1.0), /* tpfa */
                                0);
                        if (refsymbol != symbol) {
                            int cVQ1 = LAST(std::get<0>(fmt_tki_tup).cVQ1);
                            int cVQ2 = LAST(std::get<0>(fmt_tki_tup).cVQ2);
                            maxVQ_VQ1_VQ2_symbol_indelstr_tup_vec.push_back(std::make_tuple(MAX(cVQ1, cVQ2), cVQ1, cVQ2, symbol, LAST(std::get<0>(fmt_tki_tup).gapSa)));
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
                    auto reffmt = init_fmt;
                    bool is_ref_found = false;
                    for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                        if (refsymbol == (AlignmentSymbol)(LAST(std::get<0>(fmt_tki_tup).VTI))) {  
                            reffmt = std::get<0>(fmt_tki_tup); 
                            is_ref_found = true; 
                        }
                    }
                    assert(is_ref_found || !fprintf(stderr, "The position %s:%d with symbolType %d has no ref!\n", 
                            std::get<0>(tname_tseqlen_tuple).c_str(), refpos, symbolType));
                    for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                        streamFrontPushBcfFormatR(std::get<0>(fmt_tki_tup), reffmt);
                    }
                    
                    int nlodq = 31;
                    std::string note1;
                    std::vector<std::pair<AlignmentSymbol, bcfrec::BcfFormat*>> symbol_format_vec;
                    for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                        auto & fmt = std::get<0>(fmt_tki_tup);
                        auto symbol = (AlignmentSymbol)(LAST(fmt.VTI)); // std::get<1>(fmtinfo);
                        note1 += std::string("/symb/") + std::to_string(symbol) + 
                                std::string("/gVQ1/CONTQ/") + std::to_string(LAST(fmt.gVQ1)) + "/" + std::to_string(LAST(fmt.CONTQ))  + "//";
                        if (symbol != BASE_NN && symbol != LINK_NN) {
                            symbol_format_vec.push_back(std::make_pair(symbol, &fmt));
                        }
                    }
                    clear_push(init_fmt.gVQ1, 0); // can be a very negative number to force out all homref alleles
                    clear_push(init_fmt.CONTQ, 0);
                    clear_push(init_fmt.cDP1f, 0);
                    clear_push(init_fmt.cDP1r, 0);
                    clear_push(init_fmt.cDP1a, 0);
                    clear_push(init_fmt.cDP1v, 50);
                    while (symbol_format_vec.size() < 4) {
                        symbol_format_vec.push_back(std::make_pair(SYMBOL_TYPE_TO_NO_VAR_SYMBOL[symbolType], &init_fmt));
                    }
                    const size_t string_pass_old_size = buf_out_string_pass.size();
                    auto nlodq_fmtptr1_fmtptr2_tup = output_germline(
                            buf_out_string_pass,
                            refsymbol,
                            symbol_format_vec,
                            std::get<0>(tname_tseqlen_tuple).c_str(), // tname,
                            refstring,
                            refpos,
                            extended_inclu_beg_pos, // regionpos,
                            paramset.central_readlen,
                            paramset.outvar_flag, // OUTVAR_GERMLINE
                            paramset.should_output_all_germline,
                            NOT_PROVIDED != paramset.vcf_tumor_fname,
                            0);
                    const bool is_germline_var_generated = (buf_out_string_pass.size() > string_pass_old_size);

                    for (auto & fmt_tki_tup : fmt_tki_tup_vec) {
                        auto & fmt = std::get<0>(fmt_tki_tup);
                        const auto & symbol = (AlignmentSymbol)(LAST(fmt.VTI)); // std::get<1>(fmtinfo);
                        auto & tki = std::get<1>(fmt_tki_tup);
                        const bool will_generate_out = (NOT_PROVIDED == paramset.vcf_tumor_fname 
                                ? (paramset.outvar_flag & OUTVAR_ANY)
                                : (tki.ref_alt.size() > 0 && (paramset.outvar_flag & OUTVAR_SOMATIC)));
                        if (will_generate_out) {
                            fmt.GT = TT_HETERO[0];
                            nlodq = std::get<0>(nlodq_fmtptr1_fmtptr2_tup);
                            assert ((NOT_PROVIDED == paramset.vcf_tumor_fname) == (0 == tki.ref_alt.size()));
                            if (NOT_PROVIDED != paramset.vcf_tumor_fname) {
                                int nlodq_inc = 999;
                                const auto fmtptrs = std::vector<bcfrec::BcfFormat*> {{ std::get<1>(nlodq_fmtptr1_fmtptr2_tup), std::get<2>(nlodq_fmtptr1_fmtptr2_tup) }} ;
                                for (const auto *fmtptr : fmtptrs) {
                                    const AlignmentSymbol normsymbol = AlignmentSymbol(LAST(fmtptr->VTI));
                                    auto bgerr_norm_max_ad = collectget(fmtptr->cDP1x, 1, 50); 
                                    //auto bgerr_norm_max_ad = MAX(
                                    //        collectget(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->cDP1x, 1, 50),
                                    //        collectget(std::get<2>(nlodq_fmtptr1_fmtptr2_tup)->cDP1x, 1, 50));
                                    double tAD = (tki.cDP1x + 1) / 100.0;
                                    double tDP = (tki.CDP1x + 2) / 100.0;
                                    double nAD = (bgerr_norm_max_ad + 1) / 100.0;
                                    double nDP = (fmtptr->CDP1x[0] + 2) / 100.0;
                                    double bjpfrac = ((tAD) / (tDP)) / ((nAD) / (nDP));
                                    int binom_b10log10like = (int)calc_binom_10log10_likeratio((tDP - tAD) / (tDP), nDP - nAD, nAD);
                                    int powlaw_b10log10like = (int)(3 * 10 / log(10) * log(bjpfrac));
                                    int triallele_inc = ((normsymbol != symbol) ? (isSymbolSubstitution(symbol) ? 25 : 6) : 0);
                                    const auto new_nlodq_inc = BETWEEN(MIN(binom_b10log10like, powlaw_b10log10like) - 8, 0, 80) + triallele_inc;
                                    if (nlodq_inc > new_nlodq_inc) {
                                        nlodq_inc = new_nlodq_inc;
                                        fmt.note += std::string("/nlodqDeltaIs/") 
                                        + other_join(std::array<double, 4> {{ tAD, tDP, nAD, nDP }}, "/") + "#" 
                                        + std::to_string(std::get<0>(nlodq_fmtptr1_fmtptr2_tup)) + "#"
                                        + SYMBOL_TO_DESC_ARR[FIRST(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->VTI)] + "//"
                                        + SYMBOL_TO_DESC_ARR[LAST(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->VTI)] + "//"
                                        + SYMBOL_TO_DESC_ARR[LAST(std::get<2>(nlodq_fmtptr1_fmtptr2_tup)->VTI)] + "//"
                                        + other_join(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->gVQ1, "/") + "//"
                                        + other_join(std::get<2>(nlodq_fmtptr1_fmtptr2_tup)->gVQ1, "/") + "//"
                                        + other_join(std::get<1>(nlodq_fmtptr1_fmtptr2_tup)->CONTQ, "/") + "//"
                                        + other_join(std::get<2>(nlodq_fmtptr1_fmtptr2_tup)->CONTQ, "/") + "//"
                                        ;
                                    }
                                }
                                nlodq += nlodq_inc; // BETWEEN(MIN(binom_b10log10like, powlaw_b10log10like), 0, 80);
                                fmt.note += note1;
                            }
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
                                    // std::get<3>(fmtinfo), // indelstring,
                                    fmt,
                                    tki,
                                    paramset.is_tumor_format_retrieved,
                                    nlodq,
                                    paramset.vqual,
                                    paramset.vdp,
                                    paramset.vad,
                                    (paramset.should_output_all || is_germline_var_generated),
                                    bcf_hdr,
                                    0);
                        }
                    }
                }
            }
        }    
    }
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " starts destroying bam records"; }
    for (auto strand_readset : umi_strand_readset) {
        for (unsigned int strand = 0; strand < 2; strand++) {
            auto readset = strand_readset.first[strand]; 
            for (auto read : readset) {
                for (bam1_t *b : read) {
                    bam_destroy1(b);
                } 
            }
        }
    }
    if (!is_vcf_out_pass_to_stdout) {
        bgzip_string(outstring_pass, buf_out_string_pass);
    } else {
        outstring_pass += buf_out_string_pass;
    }
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " is done with current task"; }
    return 0;
};

int 
main(int argc, char **argv) {
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    const char *UMI_STRUCT = getenv("ONE_STEP_UMI_STRUCT");
    const std::string UMI_STRUCT_STRING = ((UMI_STRUCT != NULL && strlen(UMI_STRUCT) > 0) ? std::string(UMI_STRUCT) : std::string(""));
    CommandLineArgs paramset;
    int parsing_result_flag = -1;
    SequencingPlatform inferred_sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    int parsing_result_ret = paramset.initFromArgCV(parsing_result_flag, inferred_sequencing_platform, argc, argv);
    if (parsing_result_ret || parsing_result_flag) {
        return parsing_result_ret; 
    }
    LOG(logINFO) << "Program " << argv[0] << " version " << VERSION_DETAIL;
    LOG(logINFO) << "<GIT_DIFF_FULL_DISPLAY_MSG>"; 
    LOG(logINFO) << GIT_DIFF_FULL;
    LOG(logINFO) << "</GIT_DIFF_FULL_DISPLAY_MSG>";
    
    std::vector<std::tuple<std::string, unsigned int>> tid_to_tname_tseqlen_tuple_vec;
    samfname_to_tid_to_tname_tseq_tup_vec(tid_to_tname_tseqlen_tuple_vec, paramset.bam_input_fname);
    
    const unsigned int nthreads = paramset.max_cpu_num;
    bool is_vcf_out_empty_string = (std::string("") == paramset.vcf_output_fname);
    BGZF *fp_allp = NULL;
    if (!is_vcf_out_empty_string) { 
        fp_allp = bgzf_open(paramset.vcf_output_fname.c_str(), "w");
        if (NULL == fp_allp) {
            LOG(logERROR) << "Unable to open the bgzip file " << paramset.vcf_output_fname;
            exit(-8);
        }
    }
    bool is_vcf_out_pass_empty_string = (std::string("") == paramset.vcf_out_pass_fname);
    bool is_vcf_out_pass_to_stdout = (std::string("-") == paramset.vcf_out_pass_fname);
    BGZF *fp_pass = NULL;
    if (!is_vcf_out_pass_empty_string && !is_vcf_out_pass_to_stdout) {  
        fp_pass = bgzf_open(paramset.vcf_out_pass_fname.c_str(), "w");
        if (NULL == fp_pass) {
            LOG(logERROR) << "Unable to open the bgzip file " << paramset.vcf_out_pass_fname;
            exit(-9);
        }
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
    if (NOT_PROVIDED != paramset.bed_out_fname) {
        bed_out.open(paramset.bed_out_fname, std::ios::out);
    }

#if defined(USE_STDLIB_THREAD)
    const unsigned int nidxs = nthreads * 2 + 1;
#else
    const unsigned int nidxs = nthreads;
#endif
    
    bcf_hdr_t *g_bcf_hdr = NULL;
    const char *g_sample = NULL;
    if (NOT_PROVIDED != paramset.vcf_tumor_fname) {
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
    
    // VcStats all_vc_stats;

    bam_hdr_t * samheader = sam_hdr_read(samfiles[0]);
    std::string header_outstring = generate_vcf_header(
            paramset.fasta_ref_fname.c_str(), 
            SEQUENCING_PLATFORM_TO_DESC.at(inferred_sequencing_platform).c_str(), 
            paramset.central_readlen, 
            //paramset.minABQ_pcr_snv, 
            //paramset.minABQ_pcr_indel, 
            //paramset.minABQ_cap_snv, 
            //paramset.minABQ_cap_indel, 
            argc, 
            argv, 
            samheader->n_targets, 
            samheader->target_name, 
            samheader->target_len,
            paramset.sample_name.c_str(), 
            g_sample, 
            paramset.is_tumor_format_retrieved);
    clearstring<false>(fp_allp, header_outstring);
    clearstring<false>(fp_pass, header_outstring, is_vcf_out_pass_to_stdout);

    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> tid_beg_end_e2e_tuple_vec1;
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> tid_beg_end_e2e_tuple_vec2;
    std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, std::vector<TumorKeyInfo>> tid_pos_symb_to_tkis1; 
    std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, std::vector<TumorKeyInfo>> tid_pos_symb_to_tkis2; 
    SamIter samIter(paramset.bam_input_fname, paramset.tier1_target_region, paramset.bed_region_fname, nthreads); 
    unsigned int n_sam_iters = 0;    
    int iter_nreads = samIter.iternext(tid_beg_end_e2e_tuple_vec1);
    LOG(logINFO) << "PreProcessed " << iter_nreads << " reads in super-contig no " << (n_sam_iters);
    // rescue_variants_from_vcf
    tid_pos_symb_to_tkis1 = rescue_variants_from_vcf(tid_beg_end_e2e_tuple_vec1, tid_to_tname_tseqlen_tuple_vec, paramset.vcf_tumor_fname, g_bcf_hdr, paramset.is_tumor_format_retrieved);
    LOG(logINFO) << "Rescued/retrieved " << tid_pos_symb_to_tkis1.size() << " variants in super-contig no " << (n_sam_iters);
    while (iter_nreads > 0) {
        n_sam_iters++;
        std::thread read_bam_thread([&tid_beg_end_e2e_tuple_vec2, &tid_pos_symb_to_tkis2, &samIter, &iter_nreads, &n_sam_iters, &paramset, &tid_to_tname_tseqlen_tuple_vec, g_bcf_hdr]() {
            tid_beg_end_e2e_tuple_vec2.clear();
            iter_nreads = samIter.iternext(tid_beg_end_e2e_tuple_vec2);
            LOG(logINFO) << "PreProcessed " << iter_nreads << " reads in super-contig no " << (n_sam_iters);
            
            tid_pos_symb_to_tkis2 = rescue_variants_from_vcf(tid_beg_end_e2e_tuple_vec2, tid_to_tname_tseqlen_tuple_vec, paramset.vcf_tumor_fname, g_bcf_hdr, paramset.is_tumor_format_retrieved);
            LOG(logINFO) << "Rescued/retrieved " << tid_pos_symb_to_tkis2.size() << " variants in super-contig no " << (n_sam_iters);
        });
        const auto & tid_beg_end_e2e_tuple_vec = tid_beg_end_e2e_tuple_vec1; 
        const std::string bedstring_header = std::string("The BED-genomic-region is as follows (") + std::to_string(tid_beg_end_e2e_tuple_vec.size()) 
                + " chunks) for super-contig no " + std::to_string(n_sam_iters-1) + "\n";
        std::string bedstring = "";
        for (const auto & tid_beg_end_e2e_tuple : tid_beg_end_e2e_tuple_vec) {
            bedstring += (std::get<0>(tid_to_tname_tseqlen_tuple_vec[std::get<0>(tid_beg_end_e2e_tuple)]) + "\t"
                  + std::to_string(std::get<1>(tid_beg_end_e2e_tuple)) + "\t"
                  + std::to_string(std::get<2>(tid_beg_end_e2e_tuple)) + "\t"
                  + std::to_string(std::get<3>(tid_beg_end_e2e_tuple)) + "\t"
                  + "NumberOfReadsInThisInterval\t"
                  + std::to_string(std::get<4>(tid_beg_end_e2e_tuple)) + "\t" 
                  + "\n");
        }
        LOG(logINFO) << bedstring_header << bedstring;
        if (bed_out.is_open()) { bed_out << bedstring; }
        const unsigned int allridx = 0;  
        const unsigned int incvalue = tid_beg_end_e2e_tuple_vec.size();
        
        unsigned int nreads = 0;
        unsigned int npositions = 0;
        for (unsigned int j = 0; j < incvalue; j++) {
            auto region_idx = allridx + j;
            nreads += std::get<4>(tid_beg_end_e2e_tuple_vec[region_idx]);
            npositions += std::get<2>(tid_beg_end_e2e_tuple_vec[region_idx]) - std::get<1>(tid_beg_end_e2e_tuple_vec[region_idx]); 
        }
        
        assert(incvalue > 0);
        
        // distribute inputs as evenly as possible
#if defined(USE_STDLIB_THREAD)
        const unsigned int UNDERLOAD_RATIO = 1;
#else
        const unsigned int UNDERLOAD_RATIO = 4;
#endif
        unsigned int curr_nreads = 0;
        unsigned int curr_npositions = 0;
        unsigned int curr_zerobased_region_idx = 0;
        std::vector<std::pair<unsigned int, unsigned int>> beg_end_pair_vec;
        for (unsigned int j = 0; j < incvalue; j++) {
            auto region_idx = allridx + j;
            curr_nreads += std::get<4>(tid_beg_end_e2e_tuple_vec[region_idx]);
            curr_npositions += std::get<2>(tid_beg_end_e2e_tuple_vec[region_idx]) - std::get<1>(tid_beg_end_e2e_tuple_vec[region_idx]); 
            if (curr_nreads * nthreads * UNDERLOAD_RATIO > nreads || curr_npositions * nthreads * UNDERLOAD_RATIO > npositions || (j == incvalue - 1)) {
                beg_end_pair_vec.push_back(std::make_pair(curr_zerobased_region_idx, j+1));
                curr_nreads = 0;
                curr_npositions = 0;
                curr_zerobased_region_idx = j+1;
            }
        }
        
        LOG(logINFO) << "Will process the chunks from " << allridx << " to " << allridx + incvalue
                << " which contains approximately " << nreads << " reads and " << npositions << " positions divided into " 
                << beg_end_pair_vec.size() << " sub-chunks";
        
#if defined(USE_STDLIB_THREAD)
        if ( nidxs <= beg_end_pair_vec.size()) {abort();}
        std::vector<std::thread> threads; 
        threads.reserve(beg_end_pair_vec.size());
#endif
        std::vector<BatchArg> batchargs;
        batchargs.reserve(beg_end_pair_vec.size());
        for (unsigned int beg_end_pair_idx = 0; beg_end_pair_idx < beg_end_pair_vec.size(); beg_end_pair_idx++) {
            struct BatchArg a = {
                    outstring_allp : "",
                    outstring_pass : "",
                    // vc_stats : VcStats(),
                    thread_id : 0,
                    hts_idx : NULL, 
                    ref_faidx : NULL,
                    bcf_hdr : g_bcf_hdr,
                    sr : NULL,
                    
                    tid_beg_end_e2e_tuple : tid_beg_end_e2e_tuple_vec.at(0),
                    tname_tseqlen_tuple : tid_to_tname_tseqlen_tuple_vec.at(0),
                    region_ordinal : n_sam_iters,
                    region_tot_num : (unsigned int)(INT32_MAX - 1),
                    regionbatch_ordinal : 0,
                    regionbatch_tot_num : 0,

                    paramset : paramset, 
                    UMI_STRUCT_STRING : UMI_STRUCT_STRING,
                    is_vcf_out_pass_to_stdout : is_vcf_out_pass_to_stdout,
                    is_vcf_out_empty_string : is_vcf_out_empty_string
            };
            batchargs.push_back(a);
        }
        unsigned int beg_end_pair_size = beg_end_pair_vec.size();

#if defined(_OPENMP) && !defined(USE_STDLIB_THREAD)
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
#endif
        for (unsigned int beg_end_pair_idx = 0; beg_end_pair_idx < beg_end_pair_size; beg_end_pair_idx++) {
            
#if defined(_OPENMP) && !defined(USE_STDLIB_THREAD)
            size_t thread_id = omp_get_thread_num();
#elif defined(USE_STDLIB_THREAD)
            size_t thread_id = beg_end_pair_idx;
#else
            size_t thread_id = 0;
#endif
            auto & batcharg = batchargs[beg_end_pair_idx];
            batcharg.thread_id = thread_id;
            batcharg.hts_idx = sam_idxs[thread_id];
            batcharg.ref_faidx = ref_faidxs[thread_id];
            batcharg.sr = srs[thread_id];

            std::pair<unsigned int, unsigned int> beg_end_pair = beg_end_pair_vec[beg_end_pair_idx];
#if defined(USE_STDLIB_THREAD)
            std::thread athread([
                        &batcharg, allridx, beg_end_pair, beg_end_pair_idx, &tid_beg_end_e2e_tuple_vec, &tid_to_tname_tseqlen_tuple_vec, &tid_pos_symb_to_tkis1
                        ]() {
#endif
                    LOG(logINFO) << "Thread " << batcharg.thread_id << " will process the sub-chunk " << beg_end_pair_idx << " which ranges from " 
                            << beg_end_pair.first << " to " << beg_end_pair.second;
                    
                    for (unsigned int j = beg_end_pair.first; j < beg_end_pair.second; j++) {
                        batcharg.regionbatch_ordinal = j;
                        batcharg.regionbatch_tot_num = beg_end_pair.second;
                        batcharg.tid_beg_end_e2e_tuple = tid_beg_end_e2e_tuple_vec.at(allridx + j);
                        batcharg.tname_tseqlen_tuple = tid_to_tname_tseqlen_tuple_vec.at(std::get<0>(batcharg.tid_beg_end_e2e_tuple));
                        process_batch(batcharg, tid_pos_symb_to_tkis1);
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
        for (unsigned int beg_end_pair_idx = 0; beg_end_pair_idx < beg_end_pair_vec.size(); beg_end_pair_idx++) {
            if (batchargs[beg_end_pair_idx].outstring_allp.size() > 0) {
                clearstring<true>(fp_allp, batchargs[beg_end_pair_idx].outstring_allp); // empty string means end of file
            }
            if (batchargs[beg_end_pair_idx].outstring_pass.size() > 0) {
                clearstring<true>(fp_pass, batchargs[beg_end_pair_idx].outstring_pass); // empty string means end of file
            }
            // all_vc_stats.update(batchargs[beg_end_pair_idx].vc_stats);
        }
        read_bam_thread.join(); // end this iter
        for (auto tid_pos_symb_to_tkis1_pair: tid_pos_symb_to_tkis1) {
            for (auto tki : tid_pos_symb_to_tkis1_pair.second) {
                if (NULL != tki.bcf1_record) {
                    bcf_destroy(tki.bcf1_record); 
                }
            }
        }
        autoswap(tid_beg_end_e2e_tuple_vec1, tid_beg_end_e2e_tuple_vec2);
        autoswap(tid_pos_symb_to_tkis1, tid_pos_symb_to_tkis2);
    }
    if (NOT_PROVIDED != paramset.vc_stats_fname) {
        // std::ofstream vc_stats_ofstream(paramset.vc_stats_fname.c_str());
        // all_vc_stats.write_tsv(vc_stats_ofstream);
    } else {
        // all_vc_stats.write_tsv(std::cerr);
    }
    clearstring<true>(fp_allp, std::string("")); // write end of file
    clearstring<true>(fp_pass, std::string(""), is_vcf_out_pass_to_stdout);
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
    if (fp_allp) {
        int closeresult = bgzf_close(fp_allp);
        if (closeresult != 0) {
            LOG(logERROR) << "Unable to close the bgzip file " << paramset.vcf_output_fname;
        }
    }
    if (fp_pass) {
        int closeresult = bgzf_close(fp_pass);
        if (closeresult != 0) {
            LOG(logERROR) << "Unable to close the bgzip file " << paramset.vcf_output_fname;
        }
    }
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
 
    std::cerr << std::fixed << std::setprecision(2) << "CPU time used: "
              << 1.0 * (c_end-c_start) / CLOCKS_PER_SEC << " seconds\n"
              << "Wall clock time passed: "
              << std::chrono::duration<double>(t_end-t_start).count()
              << " seconds\n";
    return 0;
}

