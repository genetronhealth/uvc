#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/synced_bcf_reader.h"
#include "CmdLineArgs.hpp"
#include "consensus.hpp"
#include "grouping.hpp"
#include "version.h"

#include <chrono>
#include <ctime>
#include <thread>
#if defined(USE_STDLIB_THREAD)
#else
#include "omp.h"
#endif

#include "common.h"

const unsigned int G_BLOCK_SIZE = 1000;

void xfree(void *ptr) {
    if (NULL != ptr) { free(ptr); }
}

template <bool TIsInputCompressed>
int clearstring(BGZF * bgzip_file, const std::string & outstring_allp, bool is_output_to_stdout = false, bool flush = !TIsInputCompressed) {
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

std::string load_refstring(const faidx_t *ref_faidx, unsigned int tid, unsigned int incbeg, unsigned int excend) {
    assert(incbeg < excend);
    if (NULL == ref_faidx) {
        return std::string(excend - incbeg, 'n');
    }
    const char *tname = faidx_iseq(ref_faidx, tid);
    int regionlen;
    char *fetchedseq = faidx_fetch_seq(ref_faidx, tname, incbeg, excend - 1, &regionlen);
    assert (regionlen == (excend - incbeg) || !fprintf(stderr, "%d == %d - %d failed", regionlen, excend, incbeg));
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
        fprintf(stderr, "The library function malloc failed at line %s in file %s !\n", __LINE__, __FILE__);
        exit(-1);
    }
    size_t compressed_capacity = uncompressed_outstring.size();
    size_t compressed_totlen = 0;
    size_t uncompress_totlen = 0;
    do {
        if (compressed_totlen + BGZF_BLOCK_SIZE >= compressed_capacity) {
            char *compressed_outstr_tmp = (char*)realloc(compressed_outstr, (compressed_capacity * 2 + BGZF_BLOCK_SIZE) * sizeof(char));
            if (NULL == compressed_outstr_tmp) {
                fprintf(stderr, "The library function realloc failed at line %s in file %s !\n", __LINE__, __FILE__);
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
}

struct BatchArg {
    std::string outstring_allp;
    std::string outstring_pass;
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
    std::array<unsigned int, 3> fq_tsum_depths = {0, 0, 0};
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
genomicRegionInfoToString(const std::string & chromosome, 
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

struct TumorKeyInfo {
    std::string ref_alt;
    int32_t pos = 0;
    float VAQ = 0;
    int32_t DP = 0;
    float FA = 0;
    float FR = 0;
    int32_t bDP = 0;
    float bFA = 0;
    int32_t AutoBestAllBQ = 0;
    int32_t AutoBestAltBQ = 0;
    int32_t AutoBestRefBQ = 0;
    bcf1_t *bcf1_record = NULL;
    /*
    ~TumorKeyInfo() {
        if (bcf1_record != NULL) {
            // this line must be elsewhere apparently due to the subtle differences between C and C++.
            // bcf_destroy(bcf1_record);
        }
    }
    */
};

std::string als_to_string(const char *const* const allele, unsigned int m_allele) {
    std::string ret;
    ret.reserve(m_allele*2);
    for (unsigned int i = 0; i < m_allele; i++) {
        if (0 == i) {
            ret += allele[i];
        } else if (1 == i) {
            ret += std::string("\t") + allele[i];
        } else {
            ret += std::string(",") + allele[i];
        }
    }
    /*
    ret.reserve(m_als+1);
    unsigned int n_alleles = 0;
    for (unsigned int i = 0; i < m_als; i++) {
        if ('\0' == als[i]) {
            if (0 == n_alleles) {
                ret += "\t";
            } else if ((i+1) != m_als) {
                ret += ",";
            } else {
                // pass
            }
            n_alleles++;
        } else {
            ret += als[i];
        }
    }
    */
    return ret;
}

const std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, TumorKeyInfo>
rescue_variants_from_vcf(const auto & tid_beg_end_e2e_vec, const auto & tid_to_tname_tlen_tuple_vec, const std::string & vcf_tumor_fname, const auto *bcf_hdr, 
        const bool is_tumor_format_retrieved) {
    std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, TumorKeyInfo> ret;
    if (vcf_tumor_fname.size() == 0) {
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
    
    int retflag = bcf_sr_set_regions(sr, regionstring.c_str(), false);
    // bcf_sr_set_opt(srs[i], BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF); 
    int sr_set_opt_retval = bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    int sr_add_reader_retval = bcf_sr_add_reader(sr, vcf_tumor_fname.c_str());
    if (sr_add_reader_retval != 1) {
        LOG(logCRITICAL) << "Failed to synchronize-read the tumor vcf " << vcf_tumor_fname << " with return code " << sr_add_reader_retval;
        exit(-7);
    }
    
    int valsize = 0; 
    int ndst_val = 0;
    char *tVType = NULL;
    float *tVAQ = NULL;
    int32_t *tDP = NULL;
    float *tFA = NULL;
    float *tFR = NULL;
    int32_t *tAutoBestAllBQ = NULL;
    int32_t *tAutoBestAltBQ = NULL;
    int32_t *tAutoBestRefBQ = NULL;
    
    while (bcf_sr_next_line(sr)) {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        ndst_val = 0;
        valsize = bcf_get_format_char(bcf_hdr, line, "VType", &tVType, &ndst_val);
        if (valsize <= 0) { continue; }
        assert(ndst_val == valsize);
        std::string desc(tVType);
        // LOG(logINFO) << "Trying to retrieve the symbol " << desc << " at pos " << line->pos << " valsize = " << valsize << " ndst_val = " << ndst_val;
        AlignmentSymbol symbol = DESC_TO_SYMBOL_MAP.at(desc);
        
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
        ndst_val = 0;
        valsize = bcf_get_format_float(bcf_hdr, line, "VAQ", &tVAQ, &ndst_val);
        assert(ndst_val == valsize && valsize > 0 || !fprintf(stderr, "%d == %d && %d > 0 failed!", ndst_val, valsize, valsize));
        ndst_val = 0;
        valsize = bcf_get_format_float(bcf_hdr, line,  "FA", &tFA,  &ndst_val);
        assert(ndst_val == valsize && valsize > 0 || !fprintf(stderr, "%d == %d && %d > 0 failed!", ndst_val, valsize, valsize));
        ndst_val = 0;
        valsize = bcf_get_format_float(bcf_hdr, line,  "FR", &tFR,  &ndst_val);
        assert(ndst_val == valsize && valsize > 0 || !fprintf(stderr, "%d == %d && %d > 0 failed!", ndst_val, valsize, valsize)); 
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line,  "DP", &tDP,  &ndst_val);
        assert(ndst_val == valsize && valsize > 0 || !fprintf(stderr, "%d == %d && %d > 0 failed!", ndst_val, valsize, valsize));
        
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line,  "cAllBQ", &tAutoBestAllBQ, &ndst_val);
        assert(2 == ndst_val && 2 == valsize || !fprintf(stderr, "2 == %d == %d failed!", ndst_val, valsize));
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line,  "cAltBQ", &tAutoBestAltBQ, &ndst_val);
        assert(2 == ndst_val && 2 == valsize || !fprintf(stderr, "2 == %d == %d failed!", ndst_val, valsize));
        ndst_val = 0;
        valsize = bcf_get_format_int32(bcf_hdr, line,  "cRefBQ", &tAutoBestRefBQ, &ndst_val);
        assert(2 == ndst_val && 2 == valsize || !fprintf(stderr, "2 == %d == %d failed!", ndst_val, valsize));
          
        TumorKeyInfo tki;
        tki.VAQ = tVAQ[0];
        tki.DP = tDP[0];
        tki.FA = tFA[0];
        tki.FR = tFR[0];
        tki.AutoBestAllBQ = tAutoBestAllBQ[0] + tAutoBestAllBQ[1];
        tki.AutoBestAltBQ = tAutoBestAltBQ[0] + tAutoBestAltBQ[1];
        tki.AutoBestRefBQ = tAutoBestRefBQ[0] + tAutoBestRefBQ[1];
        tki.pos = line->pos;
        bcf_unpack(line, BCF_UN_STR);
        tki.ref_alt = als_to_string(line->d.allele, line->d.m_allele);
        if (is_tumor_format_retrieved) {
            tki.bcf1_record = bcf_dup(line);
            // bcf_unpack(tki.bcf1_record, BCF_UN_STR);
        }
        const auto retkey = std::make_tuple(line->rid, symbolpos, symbol);            
        ret.insert(std::make_pair(retkey, tki));
    }
    xfree(tVType);
    xfree(tVAQ);
    xfree(tDP);
    xfree(tFA);
    xfree(tFR);
    xfree(tAutoBestAllBQ);
    xfree(tAutoBestAltBQ);
    xfree(tAutoBestRefBQ);
    bcf_sr_destroy(sr);
    return ret;
}

int 
process_batch(BatchArg & arg, const auto & tid_pos_symb_to_tki) {
    
    std::string & outstring_allp = arg.outstring_allp;
    std::string & outstring_pass = arg.outstring_pass;
    const hts_idx_t *const hts_idx = arg.hts_idx;
    const faidx_t *const ref_faidx = arg.ref_faidx;
    
    const bcf_hdr_t *const bcf_hdr = arg.bcf_hdr;
    // bcf_srs_t *const sr = arg.sr;
    
    const CommandLineArgs & paramset = arg.paramset;
    const std::string UMI_STRUCT_STRING = arg.UMI_STRUCT_STRING;
    const std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int> tid_beg_end_e2e_tuple = arg.tid_beg_end_e2e_tuple;
    const std::tuple<std::string, unsigned int> tname_tseqlen_tuple = arg.tname_tseqlen_tuple;
    //const unsigned int region_ordinal = arg.region_ordinal;
    //const unsigned int region_tot_num = arg.region_tot_num;
    const unsigned int regionbatch_ordinal = arg.regionbatch_ordinal;
    const unsigned int regionbatch_tot_num = arg.regionbatch_tot_num;
    const unsigned int thread_id = arg.thread_id;
    const bool should_output_all = !arg.is_vcf_out_empty_string;
    const bool should_let_all_pass = paramset.should_let_all_pass;
    const bool is_vcf_out_pass_to_stdout = arg.is_vcf_out_pass_to_stdout;
    
    bool is_loginfo_enabled = (ispowerof2(regionbatch_ordinal + 1) || ispowerof2(regionbatch_tot_num - regionbatch_ordinal));
    std::string raw_out_string;
    std::string raw_out_string_pass;
    // faidx_t *ref_faidx = (fasta_ref_fname.size() > 0 ? fai_load(fasta_ref_fname.c_str()) : NULL);
    
    auto tid = std::get<0>(tid_beg_end_e2e_tuple);
    auto incluBegPosition = std::get<1>(tid_beg_end_e2e_tuple);
    auto excluEndPosition = std::get<2>(tid_beg_end_e2e_tuple);
    auto end2end = std::get<3>(tid_beg_end_e2e_tuple);
    auto nreads = std::get<4>(tid_beg_end_e2e_tuple);
    
    std::map<uint64_t, std::pair<std::array<std::map<uint64_t, std::vector<bam1_t *>>, 2>, int>> umi_to_strand_to_reads;
    unsigned int extended_inclu_beg_pos, extended_exclu_end_pos; 
    std::vector<std::pair<std::array<std::vector<std::vector<bam1_t *>>, 2>, int>> umi_strand_readset;

    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts bamfname_to_strand_to_familyuid_to_reads with pair_end_merge = " << paramset.pair_end_merge; }
    std::array<unsigned int, 3> passed_pcrpassed_umipassed = bamfname_to_strand_to_familyuid_to_reads(umi_to_strand_to_reads, 
            extended_inclu_beg_pos, extended_exclu_end_pos,
            paramset.bam_input_fname,
            tid, incluBegPosition, excluEndPosition,
            end2end, paramset.min_mapqual, paramset.min_aln_len,
            regionbatch_ordinal, regionbatch_tot_num, UMI_STRUCT_STRING, hts_idx, 
            ASSAY_TYPE_CAPTURE != paramset.assay_type, PAIR_END_MERGE_NO != paramset.pair_end_merge, paramset.disable_duplex, thread_id);
    
    unsigned int num_passed_reads = passed_pcrpassed_umipassed[0];
    unsigned int num_pcrpassed_reads = passed_pcrpassed_umipassed[1];
    unsigned int num_umipassed_reads = passed_pcrpassed_umipassed[2];
    bool is_by_capture = ((num_pcrpassed_reads) * 2 <= num_passed_reads);
    AssayType inferred_assay_type = ((ASSAY_TYPE_AUTO == paramset.assay_type) ? (is_by_capture ? ASSAY_TYPE_CAPTURE : ASSAY_TYPE_AMPLICON) : (paramset.assay_type));
    
    if (0 == num_passed_reads) { return -1; };
    unsigned int minABQ_snv = ((ASSAY_TYPE_CAPTURE != inferred_assay_type) ? paramset.minABQ_pcr_snv : paramset.minABQ_cap_snv);
    unsigned int minABQ_indel = ((ASSAY_TYPE_CAPTURE != inferred_assay_type) ? paramset.minABQ_pcr_indel : paramset.minABQ_cap_indel);
    
    const unsigned int rpos_inclu_beg = MAX(incluBegPosition, extended_inclu_beg_pos);
    const unsigned int rpos_exclu_end = MIN(excluEndPosition, extended_exclu_end_pos); 

    //std::vector<std::array<TumorKeyInfo, NUM_ALIGNMENT_SYMBOLS>> extended_posidx_to_symbol_to_tkinfo(extended_exclu_end_pos - extended_inclu_beg_pos + 1);
    const auto tki_beg = tid_pos_symb_to_tki.lower_bound(std::make_tuple(tid, extended_inclu_beg_pos    , AlignmentSymbol(0)));
    const auto tki_end = tid_pos_symb_to_tki.upper_bound(std::make_tuple(tid, extended_exclu_end_pos + 1, AlignmentSymbol(0)));
    std::vector<bool> extended_posidx_to_is_rescued(extended_exclu_end_pos - extended_inclu_beg_pos + 1, false);
    unsigned int num_rescued = 0;
    for (auto tki_it = tki_beg; tki_it != tki_end; tki_it++) {
        auto symbolpos = std::get<1>(tki_it->first);
        extended_posidx_to_is_rescued[symbolpos - extended_inclu_beg_pos] = true;
        num_rescued++;
        if (is_loginfo_enabled) {
            // NOTE: the true positive short del at 22:17946835 in NA12878-NA24385 mixture is overwhelmed by the false positve long del spanning the true positive short del.
            // However, so far manual check with limited experience cannot confirm that the true positive is indeed a true positive.
            // TODO: have to see more examples of this case and adjust code accordingly.
            LOG(logDEBUG4) << "Thread " << thread_id << " iterated over symbolpos " << symbolpos << " and symbol " << std::get<2>(tki_it->first) << " as a rescued var";
        }
    }
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " deals with " << num_rescued << " tumor-sample variants in region " << extended_inclu_beg_pos << " to " << extended_exclu_end_pos + 1 ;}
    // auto tki_it = tki_beg; 
    //if (paramset.vcf_tumor_fname.size() != 0) {
        // do not check tumor vcf here.
    //}
    
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts converting umi_to_strand_to_reads with is_by_capture = " << is_by_capture << "  " ;}
    fill_strand_umi_readset_with_strand_to_umi_to_reads(umi_strand_readset, umi_to_strand_to_reads);
    
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts constructing symbolToCountCoverageSet12 with " << extended_inclu_beg_pos << (" , ") << extended_exclu_end_pos; }
    //std::array<Symbol2CountCoverageSet, 2> symbolToCountCoverageSet12 = {
    //    Symbol2CountCoverageSet(tid, extended_inclu_beg_pos, extended_exclu_end_pos),
    //    Symbol2CountCoverageSet(tid, extended_inclu_beg_pos, extended_exclu_end_pos)
    //};
    // + 1 accounts for insertion at the end of the region, this should happen RARELY (like 1e-20 chance)
    Symbol2CountCoverageSet symbolToCountCoverageSet12(tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1); 
    if (is_loginfo_enabled) { LOG(logINFO)<< "Thread " << thread_id << " starts updateByRegion3Aln with " << umi_strand_readset.size() << " families"; }
    std::string refstring = load_refstring(ref_faidx, tid, extended_inclu_beg_pos, extended_exclu_end_pos);
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> adjcount_x_rpos_x_misma_vec;
    // unsigned int maxvalue;
    std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> mutform2count4map_bq;
    std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> mutform2count4map_fq;
    const PhredMutationTable sscs_mut_table(
                paramset.phred_max_sscs_transition_CG_TA, 
                paramset.phred_max_sscs_transition_TA_CG, 
                paramset.phred_max_sscs_transversion_any,
                paramset.phred_max_sscs_indel_open,
                paramset.phred_max_sscs_indel_ext);
    symbolToCountCoverageSet12.updateByRegion3Aln(
            mutform2count4map_bq, mutform2count4map_fq,
            // adjcount_x_rpos_x_misma_vec, // maxvalue, 
            umi_strand_readset, refstring, 
            paramset.bq_phred_added_misma, paramset.bq_phred_added_indel, paramset.should_add_note, 
            paramset.phred_max_frag_indel_ext, paramset.phred_max_frag_indel_basemax,  
            sscs_mut_table,
            // paramset.phred_max_dscs, 
            minABQ_snv, // minABQ_indel,
            // ErrorCorrectionType(paramset.seq_data_type),
            (is_by_capture ? paramset.ess_georatio_dedup_cap : paramset.ess_georatio_dedup_pcr), paramset.ess_georatio_duped_pcr,
            !paramset.disable_dup_read_merge, 
            is_loginfo_enabled, thread_id, paramset.fixedthresBQ, paramset.nogap_phred);
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id << " starts analyzing phasing info"; }
    auto mutform2count4vec_bq = map2vector(mutform2count4map_bq);
    auto simplemut2indices_bq = mutform2count4vec_to_simplemut2indices(mutform2count4vec_bq);
    auto mutform2count4vec_fq = map2vector(mutform2count4map_fq);
    auto simplemut2indices_fq = mutform2count4vec_to_simplemut2indices(mutform2count4vec_fq);
    
    /*
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
    */
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " starts generating block gzipped vcf"; }
    
    std::string buf_out_string;
    std::string buf_out_string_pass;
    const std::set<size_t> empty_size_t_set;
    const unsigned int capDP = 10*1000*1000;
    
    unsigned int gbDPmin = capDP;
    unsigned int gbDPmax = 0;
    unsigned int gcDPmin = capDP;
    unsigned int gcDPmax = 0;
    std::string gfGTmm2 = "./.";
    unsigned int gfGQmin = capDP;
    unsigned int gfGQmax = 0;
    
    unsigned int prevPosition = rpos_inclu_beg;
    SymbolType prevSymbolType = NUM_SYMBOL_TYPES;
    for (unsigned int refpos = rpos_inclu_beg; refpos <= rpos_exclu_end; refpos++) {        
        // while (tki_it != tki_end && std::get<1>(*tki_it) < refpos) { tki_it++; }
        const std::array<SymbolType, 2> allSymbolTypes = {LINK_SYMBOL, BASE_SYMBOL};
        const std::array<SymbolType, 2> stype_to_immediate_prev = {LINK_SYMBOL, BASE_SYMBOL};
        for (unsigned int stidx = 0; stidx < 2; stidx++) {
            const SymbolType symbolType = allSymbolTypes[stidx];
            bcfrec::BcfFormat init_fmt;
            const AlignmentSymbol refsymbol = (LINK_SYMBOL == symbolType ? LINK_M : (
                    refstring.size() == (refpos - extended_inclu_beg_pos) ? BASE_NN : CHAR_TO_SYMBOL.data.at(refstring.at(refpos - extended_inclu_beg_pos))));
            std::array<unsigned int, 2> bDPcDP = BcfFormat_init(init_fmt, symbolToCountCoverageSet12, refpos, symbolType, !paramset.disable_dup_read_merge, refsymbol);
            AlignmentSymbol most_confident_symbol = END_ALIGNMENT_SYMBOLS;
            float most_confident_qual = 0;
            std::string most_confident_GT = "./.";
            float most_confident_GQ = 0;
            std::vector<bcfrec::BcfFormat> fmts(SYMBOL_TYPE_TO_INCLU_END[symbolType] - SYMBOL_TYPE_TO_INCLU_BEG[symbolType] + 1, init_fmt);
            if (rpos_exclu_end != refpos && bDPcDP[0] >= paramset.min_depth_thres) {
                for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; symbol = AlignmentSymbol(1+(unsigned int)symbol)) {
                    const auto simplemut = std::make_pair(refpos, symbol);
                    auto indices_bq = (simplemut2indices_bq.find(simplemut) != simplemut2indices_bq.end() ? simplemut2indices_bq[simplemut] : empty_size_t_set); 
                    auto indices_fq = (simplemut2indices_fq.find(simplemut) != simplemut2indices_fq.end() ? simplemut2indices_fq[simplemut] : empty_size_t_set);
                    const bool is_rescued = (
                            extended_posidx_to_is_rescued[refpos - extended_inclu_beg_pos] &&
                            (tid_pos_symb_to_tki.end() != tid_pos_symb_to_tki.find(std::make_tuple(tid, refpos, symbol)))); 
                            //(extended_posidx_to_symbol_to_tkinfo[refpos-extended_inclu_beg_pos][symbol].DP > 0); 
                    unsigned int phred_max_sscs = sscs_mut_table.to_phred_rate(refsymbol, symbol);
                    int altdepth = fillBySymbol(fmts[symbol - SYMBOL_TYPE_TO_INCLU_BEG[symbolType]], symbolToCountCoverageSet12, 
                            refpos, symbol, refstring, extended_inclu_beg_pos, mutform2count4vec_bq, indices_bq, mutform2count4vec_fq, indices_fq, 
                            ((BASE_SYMBOL == symbolType) ? minABQ_snv : minABQ_indel),
                            paramset.minMQ1, paramset.maxMQ,
                            phred_max_sscs, paramset.phred_dscs_minus_sscs + phred_max_sscs,
                            // ErrorCorrectionType(paramset.seq_data_type), 
                            !paramset.disable_dup_read_merge, !paramset.enable_dup_read_vqual,
                            is_rescued);
                }
                for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; symbol = AlignmentSymbol(1+(unsigned int)symbol)) {
                    float vaq = fmts[symbol - SYMBOL_TYPE_TO_INCLU_BEG[symbolType]].VAQ;
                    if (vaq >= most_confident_qual) {
                        most_confident_symbol = symbol;
                        most_confident_qual = vaq;
                        auto GTval = fmts[symbol - SYMBOL_TYPE_TO_INCLU_BEG[symbolType]].GT;
                        auto GQval = fmts[symbol - SYMBOL_TYPE_TO_INCLU_BEG[symbolType]].GQ;
                        most_confident_GT = GTval;
                        most_confident_GQ = GQval;
                    }
                }
            }
            if ((!paramset.is_tumor_format_retrieved) && (rpos_exclu_end != refpos || allSymbolTypes[0] == symbolType)) {
                auto bDPval = bDPcDP[0];
                auto cDPval = bDPcDP[1];
                auto fGTmm2 = (most_confident_GQ >= 10 ? most_confident_GT : "./.");
                auto fGQval = most_confident_GQ;
                bool gbDPhasWideRange = (is_sig_out(bDPval, gbDPmin, gbDPmax, 130,  3) || (0 == gbDPmax && bDPval > gbDPmax) || (0 == bDPval && gbDPmin > bDPval));
                bool gcDPhasWideRange = (is_sig_out(cDPval, gcDPmin, gcDPmax, 130,  3));
                bool gfGQhasWideRange = (is_sig_out(fGQval, gfGQmin, gfGQmax, 130, 10) || (std::string(fGTmm2) != gfGTmm2));
                if (gbDPhasWideRange || gcDPhasWideRange || gfGQhasWideRange ||
                        (refpos - prevPosition >= G_BLOCK_SIZE) || (refpos == rpos_exclu_end)) {
                    auto iendPosition = (LINK_SYMBOL == symbolType ? (refpos - 1) : refpos);
                    auto iendSymbolType = stype_to_immediate_prev[symbolType];
                    
                    unsigned int begpos = (BASE_SYMBOL == prevSymbolType ? (prevPosition+1) : prevPosition);
                    unsigned int endpos = (BASE_SYMBOL == iendSymbolType ? (iendPosition+1) : iendPosition);
                    
                    std::string genomicInfoString = ((0 == gbDPmax) ? "" : genomicRegionInfoToString(
                            std::get<0>(tname_tseqlen_tuple),
                            begpos, prevSymbolType,
                            endpos, iendSymbolType,
                            gbDPmin, gcDPmin,
                            gfGTmm2, gfGQmin,
                            refstring, extended_inclu_beg_pos));
                    raw_out_string += genomicInfoString + buf_out_string;
                    raw_out_string_pass += genomicInfoString + buf_out_string_pass;
                    buf_out_string.clear();
                    buf_out_string_pass.clear();
                    prevPosition = refpos;
                    prevSymbolType = symbolType;
                    gbDPmin = bDPval;
                    gbDPmax = bDPval;
                    gcDPmin = cDPval;
                    gcDPmax = cDPval;
                    gfGTmm2 = fGTmm2;
                    gfGQmin = fGQval;
                    gfGQmax = fGQval;
                } else {
                    UPDATE_MIN(gbDPmin, bDPval);
                    UPDATE_MAX(gbDPmax, bDPval);
                    UPDATE_MIN(gcDPmin, cDPval);
                    UPDATE_MAX(gcDPmax, cDPval);
                    UPDATE_MIN(gfGQmin, fGQval);
                    UPDATE_MAX(gfGQmax, fGQval);
                }
            }
            if (rpos_exclu_end != refpos && bDPcDP[0] >= paramset.min_depth_thres) {
                for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; symbol = AlignmentSymbol(1+(unsigned int)symbol)) {
                    auto & fmt = fmts[symbol - SYMBOL_TYPE_TO_INCLU_BEG[symbolType]];
                    const bool pass_thres = ((fmt.bAD1[0] + fmt.bAD1[1]) >= paramset.min_altdp_thres); //  (paramset.min_alt  );
                    const bool is_rescued = (
                            //(tki_it != tki_end) 
                            //&& (std::get<1>(*tki_it) == refpos) 
                            extended_posidx_to_is_rescued[refpos - extended_inclu_beg_pos] &&
                            (tid_pos_symb_to_tki.end() != tid_pos_symb_to_tki.find(std::make_tuple(tid, refpos, symbol)))); 
                            //(extended_posidx_to_symbol_to_tkinfo[refpos-extended_inclu_beg_pos][symbol].DP > 0); 
                    TumorKeyInfo tki;
                    if (is_rescued) {
                        tki = tid_pos_symb_to_tki.find(std::make_tuple(tid, refpos, symbol))->second; 
                    }
                    if (is_rescued || pass_thres) {
                        fmt.CType = SYMBOL_TO_DESC_ARR[most_confident_symbol];
                        fmt.CAQ = most_confident_qual;
                        appendVcfRecord(buf_out_string, buf_out_string_pass, symbolToCountCoverageSet12,
                                std::get<0>(tname_tseqlen_tuple).c_str(), refpos, symbol, fmt,
                                refstring, extended_inclu_beg_pos, paramset.vqual, should_output_all, should_let_all_pass,
                                tki, paramset.vcf_tumor_fname.size() != 0, paramset.phred_germline_polymorphism, 
                                paramset.nonref_to_alt_frac_snv, paramset.nonref_to_alt_frac_indel,
                                paramset.tnq_mult_snv, paramset.tnq_mult_indel
                                , paramset.mai_tier_qual // = 40;
                                , paramset.mai_tier_abq // = 40;
                                , paramset.str_tier_qual // = 50;
                                , paramset.str_tier_len // = 16;
                                , paramset.uni_bias_thres // = 180
                                , bcf_hdr);
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
    bgzip_string(outstring_allp, raw_out_string);
    if (!is_vcf_out_pass_to_stdout) {
        bgzip_string(outstring_pass, raw_out_string_pass);
    } else {
        outstring_pass += raw_out_string_pass;
    }
    if (is_loginfo_enabled) { LOG(logINFO) << "Thread " << thread_id  << " is done with current task"; }
};

int main(int argc, char **argv) {
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    const char *UMI_STRUCT = getenv("ONE_STEP_UMI_STRUCT");
    const std::string UMI_STRUCT_STRING = (UMI_STRUCT != NULL ? std::string(UMI_STRUCT) : std::string(""));
    CommandLineArgs paramset;
    int parsing_result_flag = -1;
    SequencingPlatform inferred_sequencing_platform = SEQUENCING_PLATFORM_AUTO;
    int parsing_result_ret = paramset.initFromArgCV(parsing_result_flag, inferred_sequencing_platform, argc, argv);
    if (parsing_result_ret || parsing_result_flag) {
        return parsing_result_ret; 
    }
    LOG(logINFO) << "Program " << argv[0] << " version " << VERSION_DETAIL;
    std::vector<std::tuple<std::string, unsigned int>> tid_to_tname_tseqlen_tuple_vec;
    /*
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> tid_beg_end_e2e_tuple_vec;
    LOG(logINFO) << "step1: sam_fname_to_contigs"; 
    sam_fname_to_contigs(tid_beg_end_e2e_tuple_vec, tid_to_tname_tseqlen_tuple_vec, paramset.bam_input_fname, paramset.bed_region_fname);
    LOG(logINFO) << "The BED region is as follows (" << tid_beg_end_e2e_tuple_vec.size() << " chunks):";
    for (auto tid_beg_end_e2e_tuple : tid_beg_end_e2e_tuple_vec) {
        std::cerr << std::get<0>(tid_to_tname_tseqlen_tuple_vec[std::get<0>(tid_beg_end_e2e_tuple)]) << "\t" 
                  << std::get<1>(tid_beg_end_e2e_tuple) << "\t"
                  << std::get<2>(tid_beg_end_e2e_tuple) << "\t"
                  << std::get<3>(tid_beg_end_e2e_tuple) << "\t"
                  << std::get<4>(tid_beg_end_e2e_tuple) << "\n";
    }
    */
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
    /*
    if (paramset.vcf_output_fname.size() != 0 && paramset.vcf_output_fname != "-") {
        bgzf_index_build_init(fp_allp);
    }
    */
    // bgzf_mt(fp_allp, nthreads, 128);
    // samFile *sam_infile = sam_open(paramset.bam_input_fname.c_str(), "r");

#if defined(USE_STDLIB_THREAD)
    const unsigned int nidxs = nthreads * 2 + 1;
#else
    const unsigned int nidxs = nthreads;
#endif
    
    bcf_hdr_t *g_bcf_hdr = NULL;
    const char *g_sample = NULL;
    if (paramset.vcf_tumor_fname.size() > 0) {
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
        sam_idxs[i] = sam_index_load2(samfiles[i], paramset.bam_input_fname.c_str(), NULL);
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
        if (paramset.vcf_tumor_fname.size() > 0) {
            /*
            srs[i] = bcf_sr_init();
            if (NULL == srs[i]) {
                LOG(logCRITICAL) << "Failed to initialize bcf sr for thread with ID = " << i;
                exit(-6);
            }
            // bcf_sr_set_opt(srs[i], BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
            
            int sr_set_opt_retval = bcf_sr_set_opt(srs[i], BCF_SR_REQUIRE_IDX);
            int sr_add_reader_retval = bcf_sr_add_reader(srs[i], paramset.vcf_tumor_fname.c_str());
            if (sr_add_reader_retval != 1) {
                LOG(logCRITICAL) << "Failed to synchronize-read the tumor vcf " << paramset.vcf_tumor_fname << " for thread with ID = " << i;
                exit(-7);
            }
            */
        }
    }
    
    bam_hdr_t * samheader = sam_hdr_read(samfiles[0]);
    std::string header_outstring = generateVcfHeader(paramset.fasta_ref_fname.c_str(), SEQUENCING_PLATFORM_TO_DESC.at(inferred_sequencing_platform).c_str(), 
            paramset.minABQ_pcr_snv, paramset.minABQ_pcr_indel, paramset.minABQ_cap_snv, paramset.minABQ_cap_indel, argc, argv, 
            samheader->n_targets, samheader->target_name, samheader->target_len,
            paramset.sample_name.c_str(), g_sample);
    clearstring<false>(fp_allp, header_outstring);
    clearstring<false>(fp_pass, header_outstring, is_vcf_out_pass_to_stdout);
    
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> tid_beg_end_e2e_tuple_vec1;
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> tid_beg_end_e2e_tuple_vec2;
    std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, TumorKeyInfo> tid_pos_symb_to_tki1; 
    std::map<std::tuple<unsigned int, unsigned int, AlignmentSymbol>, TumorKeyInfo> tid_pos_symb_to_tki2; 
    SamIter samIter(paramset.bam_input_fname, paramset.bed_region_fname, nthreads); 
    unsigned int n_sam_iters = 0;    
    int iter_nreads = samIter.iternext(tid_beg_end_e2e_tuple_vec1);
    LOG(logINFO) << "PreProcessed " << iter_nreads << " reads in super-contig no " << (n_sam_iters);
    // rescue_variants_from_vcf
    tid_pos_symb_to_tki1 = rescue_variants_from_vcf(tid_beg_end_e2e_tuple_vec1, tid_to_tname_tseqlen_tuple_vec, paramset.vcf_tumor_fname, g_bcf_hdr, paramset.is_tumor_format_retrieved);
    LOG(logINFO) << "Rescued/retrieved " << tid_pos_symb_to_tki1.size() << " variants in super-contig no " << (n_sam_iters);
    while (iter_nreads > 0) {
        n_sam_iters++;
        std::thread read_bam_thread([&tid_beg_end_e2e_tuple_vec2, &tid_pos_symb_to_tki2, &samIter, &iter_nreads, &n_sam_iters, &paramset, &tid_to_tname_tseqlen_tuple_vec, g_bcf_hdr]() {
            tid_beg_end_e2e_tuple_vec2.clear();
            iter_nreads = samIter.iternext(tid_beg_end_e2e_tuple_vec2);
            LOG(logINFO) << "PreProcessed " << iter_nreads << " reads in super-contig no " << (n_sam_iters);
            
            tid_pos_symb_to_tki2 = rescue_variants_from_vcf(tid_beg_end_e2e_tuple_vec2, tid_to_tname_tseqlen_tuple_vec, paramset.vcf_tumor_fname, g_bcf_hdr, paramset.is_tumor_format_retrieved);
            LOG(logINFO) << "Rescued/retrieved " << tid_pos_symb_to_tki2.size() << " variants in super-contig no " << (n_sam_iters);
        });
        const auto & tid_beg_end_e2e_tuple_vec = tid_beg_end_e2e_tuple_vec1; 
        std::string bedstring = std::string("The BED-genomic-region is as follows (") + std::to_string(tid_beg_end_e2e_tuple_vec.size()) 
                + " chunks) for super-contig no " + std::to_string(n_sam_iters-1) + "\n";
        for (const auto & tid_beg_end_e2e_tuple : tid_beg_end_e2e_tuple_vec) {
            bedstring += (std::get<0>(tid_to_tname_tseqlen_tuple_vec[std::get<0>(tid_beg_end_e2e_tuple)]) + "\t"
                  + std::to_string(std::get<1>(tid_beg_end_e2e_tuple)) + "\t"
                  + std::to_string(std::get<2>(tid_beg_end_e2e_tuple)) + "\t"
                  + std::to_string(std::get<3>(tid_beg_end_e2e_tuple)) + "\t"
                  + "NumberOfReadsInThisInterval\t"
                  + std::to_string(std::get<4>(tid_beg_end_e2e_tuple)) + "\t" 
                  + "\n");
        }
        LOG(logINFO) << bedstring;
        
        const unsigned int allridx = 0;  
        const unsigned int incvalue = tid_beg_end_e2e_tuple_vec.size();
        
        unsigned int nreads = 0;
        unsigned int npositions = 0;
        for (unsigned int j = 0; j < incvalue; j++) {
            auto region_idx = allridx + j;
            nreads += std::get<4>(tid_beg_end_e2e_tuple_vec[region_idx]);
            npositions += std::get<2>(tid_beg_end_e2e_tuple_vec[region_idx]) - std::get<1>(tid_beg_end_e2e_tuple_vec[region_idx]); 
        }

        /*
    unsigned int incvalue = 0;
    for (unsigned int allridx = 0; allridx < tid_beg_end_e2e_tuple_vec.size(); allridx += incvalue) {
        incvalue = 0;
        unsigned int nreads = 0;
        unsigned int npositions = 0;
        while (allridx + incvalue < tid_beg_end_e2e_tuple_vec.size() && nreads < (2000*1000 * nthreads)  && npositions < (1000*1000 * nthreads)) {
            nreads += std::get<4>(tid_beg_end_e2e_tuple_vec[allridx + incvalue]);
            npositions += std::get<2>(tid_beg_end_e2e_tuple_vec[allridx + incvalue]) - std::get<1>(tid_beg_end_e2e_tuple_vec[allridx + incvalue]); 
            incvalue++;
        }
        */
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
            // for (unsigned int j = 0; j < incvalue; j++) {

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
                        &batcharg, allridx, beg_end_pair, beg_end_pair_idx, &tid_beg_end_e2e_tuple_vec, &tid_to_tname_tseqlen_tuple_vec, &tid_pos_symb_to_tki1
                        //beg_end_pair_idx, allridx, is_vcf_out_pass_to_stdout, &beg_end_pair_vec, 
                        //&outstrings, &outstrings_pass, paramset, &tid_beg_end_e2e_tuple_vec, tid_to_tname_tseqlen_tuple_vec, UMI_STRUCT, 
                        //sam_idxs, samfiles, ref_faidxs
                        ]() {
#endif
                    LOG(logINFO) << "Thread " << batcharg.thread_id << " will process the sub-chunk " << beg_end_pair_idx << " which ranges from " 
                            << beg_end_pair.first << " to " << beg_end_pair.second;
                    
                    //auto *this_sam_idx = sam_index_load2(samfiles[allridx], paramset.bam_input_fname.c_str(), NULL);
                    for (unsigned int j = beg_end_pair.first; j < beg_end_pair.second; j++) {
                        batcharg.regionbatch_ordinal = j;
                        batcharg.regionbatch_tot_num = beg_end_pair.second;
                        batcharg.tid_beg_end_e2e_tuple = tid_beg_end_e2e_tuple_vec.at(allridx + j);
                        batcharg.tname_tseqlen_tuple = tid_to_tname_tseqlen_tuple_vec.at(std::get<0>(batcharg.tid_beg_end_e2e_tuple));
                        process_batch(batcharg, tid_pos_symb_to_tki1);
                    }
#if defined(USE_STDLIB_THREAD)
                    //hts_idx_destroy(this_sam_idx);
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
        }
        read_bam_thread.join(); // end this iter
        for (auto tid_pos_symb_to_tki1_pair: tid_pos_symb_to_tki1) {
            if (NULL != tid_pos_symb_to_tki1_pair.second.bcf1_record) {
                bcf_destroy(tid_pos_symb_to_tki1_pair.second.bcf1_record); 
            }
        }
        autoswap(tid_beg_end_e2e_tuple_vec1, tid_beg_end_e2e_tuple_vec2);
        autoswap(tid_pos_symb_to_tki1, tid_pos_symb_to_tki2);
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
    // sam_close(sam_infile);
    // bgzf_flush(fp_allp);
    /*
    if (paramset.vcf_output_fname.size() != 0 && paramset.vcf_output_fname != "-") {
        int index_dump_exitcode = bgzf_index_dump(fp_allp, paramset.vcf_output_fname.c_str(), ".csi");
        if (0 != index_dump_exitcode) {
            fprintf(stderr, "bgzf_index_dump to file %s.csi returned exit code %d\n",  paramset.vcf_output_fname.c_str(), index_dump_exitcode);
        }
    }
    */
    if (fp_allp) {
        // bgzf_flush(fp_allp);
        int closeresult = bgzf_close(fp_allp);
        if (closeresult != 0) {
            LOG(logERROR) << "Unable to close the bgzip file " << paramset.vcf_output_fname;
        }
    }
    if (fp_pass) {
        // bgzf_flush(fp_pass);
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
}

