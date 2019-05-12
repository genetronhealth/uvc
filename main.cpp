#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "CmdLineArgs.hpp"
#include "consensus.hpp"
#include "grouping.hpp"
#include "version.h"

#if defined(USE_STDLIB_THREAD)
#include <thread>
#else
#include "omp.h"
#endif

const unsigned int G_BLOCK_SIZE = 1000;

/*
int _old_clearstring(BGZF * bgzip_file, const std::string & outstring_allp, bool is_output_to_stdout = false, bool flush = true) {
    if (is_output_to_stdout) {
        std::cout << outstring_allp;
        return outstring_allp.size();
    }
    if (NULL == bgzip_file) { return -1; }
    int ret = bgzf_write(bgzip_file, outstring_allp.c_str(), outstring_allp.size());
    LOG(logINFO) << "Written " << ret << " bytes of compressed data from " << outstring_allp.size()  << " bytes of raw data.";
    if (flush) { 
        int flushret = bgzf_flush(bgzip_file); 
        if (flushret != 0) {
            return flushret;
        }
    }
    return ret;
};
*/

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
    
    std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int> tid_beg_end_e2e_tuple;
    std::tuple<std::string, unsigned int> tname_tseqlen_tuple;
    unsigned int region_ordinal;
    unsigned int regionbatch_ordinal;

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
genomicRegionInfoToString(const std::string & chromosome, unsigned int incluBeg, unsigned int excluEnd,
        const std::array<unsigned int, 2> & gbDPmin, const std::array<unsigned int, 2> & gcDPmin, 
        const std::array<std::string,  2> & gfGTmm2, const std::array<unsigned int, 2> & gfGQmin
        //const std::array<float, 2> gCAQamin
        ) {
    std::string ret = chromosome + "\t" + std::to_string(incluBeg) 
            + "\t.\tN\t<NON_REF" + ">\t.\t.\t.\tGT:gbDP:gcDP:" /*+ gCAQ */ + "gGT:gGQ" + ":gEND\t.:" 
            + std::to_string(gbDPmin[0]) + "," + std::to_string(gbDPmin[1]) + ":" 
            + std::to_string(gcDPmin[0]) + "," + std::to_string(gcDPmin[1]) + ":" 
            +               (gfGTmm2[0]) + "," +               (gfGTmm2[1]) + ":" 
            + std::to_string(gfGQmin[0]) + "," + std::to_string(gfGQmin[1]) + ":" 
            // + std::to_string(gCAQamin[0]) + "," + std::to_string(gCAQamin[1]) + ":" 
            + std::to_string(excluEnd) + "\n";
    return ret;
}

const bool
is_sig_higher_than(auto a, auto b, unsigned int mfact = 30, unsigned int afact = 3) {
    return (a * 100 > b * (100 + mfact)) && (a > b + afact);
}

const bool
is_sig_out(auto a, auto minval, auto maxval, unsigned int mfact = 30, unsigned int afact = 3) {
    return is_sig_higher_than(a, minval, mfact, afact) || is_sig_higher_than(maxval, a, mfact, afact);
}

int 
process_batch(BatchArg & arg) {
    
    std::string & outstring_allp = arg.outstring_allp;
    std::string & outstring_pass = arg.outstring_pass;
    const hts_idx_t *const hts_idx = arg.hts_idx;
    const faidx_t *const ref_faidx = arg.ref_faidx;
    const CommandLineArgs & paramset = arg.paramset;
    const std::string UMI_STRUCT_STRING = arg.UMI_STRUCT_STRING;
    const std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int> tid_beg_end_e2e_tuple = arg.tid_beg_end_e2e_tuple;
    const std::tuple<std::string, unsigned int> tname_tseqlen_tuple = arg.tname_tseqlen_tuple;
    const unsigned int region_ordinal = arg.region_ordinal;
    const unsigned int regionbatch_ordinal = arg.regionbatch_ordinal;
    const unsigned int thread_id = arg.thread_id;
    const bool should_output_all = !arg.is_vcf_out_empty_string;
    const bool is_vcf_out_pass_to_stdout = arg.is_vcf_out_pass_to_stdout;
    
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

    LOG(logINFO) << "Thread " << thread_id << " starts bamfname_to_strand_to_familyuid_to_reads";
    int num_pass_reads = bamfname_to_strand_to_familyuid_to_reads(umi_to_strand_to_reads, 
            extended_inclu_beg_pos, extended_exclu_end_pos,
            paramset.bam_input_fname, ErrorCorrectionType(paramset.seq_data_type),
            tid, incluBegPosition, excluEndPosition,
            end2end, paramset.min_mapqual, paramset.min_aln_len,
            region_ordinal, UMI_STRUCT_STRING, hts_idx);
    
    if (0 == num_pass_reads) { return -1; };
    
    LOG(logINFO) << "Thread " << thread_id << " starts converting umi_to_strand_to_reads";
    fill_strand_umi_readset_with_strand_to_umi_to_reads(umi_strand_readset, umi_to_strand_to_reads);
    
    LOG(logINFO) << "Thread " << thread_id << " starts constructing symbolToCountCoverageSet12 with " << extended_inclu_beg_pos << (" , ") << extended_exclu_end_pos;
    //std::array<Symbol2CountCoverageSet, 2> symbolToCountCoverageSet12 = {
    //    Symbol2CountCoverageSet(tid, extended_inclu_beg_pos, extended_exclu_end_pos),
    //    Symbol2CountCoverageSet(tid, extended_inclu_beg_pos, extended_exclu_end_pos)
    //};
    // + 1 accounts for insertion at the end of the region, this should happen RARELY (like 1e-20 chance)
    Symbol2CountCoverageSet symbolToCountCoverageSet12(tid, extended_inclu_beg_pos, extended_exclu_end_pos + 1); 
    LOG(logINFO)<< "Thread " << thread_id << " starts updateByRegion3Aln with " << umi_strand_readset.size() << " families";
    std::string refstring = load_refstring(ref_faidx, tid, extended_inclu_beg_pos, extended_exclu_end_pos);
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> adjcount_x_rpos_x_misma_vec;
    // unsigned int maxvalue;
    std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> mutform2count4map_bq;
    std::map<std::basic_string<std::pair<unsigned int, AlignmentSymbol>>, std::array<unsigned int, 2>> mutform2count4map_fq;
    symbolToCountCoverageSet12.updateByRegion3Aln(
            mutform2count4map_bq, mutform2count4map_fq,
            // adjcount_x_rpos_x_misma_vec, // maxvalue, 
            umi_strand_readset, refstring, 
            paramset.bq_phred_added_misma, paramset.bq_phred_added_indel, paramset.should_add_note, paramset.phred_max_sscs, paramset.phred_max_dscs);
    LOG(logINFO) << "Thread " << thread_id << " starts analyzing phasing info";
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
    LOG(logINFO) << "Thread " << thread_id  << " starts generating block gzipped vcf";
    
    std::string buf_out_string;
    std::string buf_out_string_pass;
    const std::set<size_t> empty_size_t_set;
    const unsigned int rpos_inclu_beg = MAX(incluBegPosition, extended_inclu_beg_pos);
    const unsigned int rpos_exclu_end = MIN(excluEndPosition, extended_exclu_end_pos);
    
    const unsigned int capDP = 10*1000*1000;
    // const float capCAQ = 1e6;
    
    std::array<unsigned int, 2> gbDPmin = {capDP, capDP};
    std::array<unsigned int, 2> gbDPmax = {0, 0};
    std::array<unsigned int, 2> gcDPmin = {capDP, capDP};
    std::array<unsigned int, 2> gcDPmax = {0, 0};
    
    std::array<std::string,  2> gfGTmm2 = {""   , ""   };
    std::array<unsigned int, 2> gfGQmin = {capDP, capDP};
    std::array<unsigned int, 2> gfGQmax = {0, 0};

    //std::array<unsigned int, 2> gCAQmin = {capDP, capDP};
    //std::array<unsigned int, 2> gCAQmax = {0, 0};
    
    std::array<unsigned int, 2> bDPval = {capDP, capDP};
    std::array<unsigned int, 2> cDPval = {capDP, capDP};
    
    std::array<std::string,  2> fGTmm2 = {""   , ""   };
    std::array<unsigned int, 2> fGQval = {capDP, capDP};

    //std::array<unsigned int, 2> CAQval = {capDP, capDP};
    
    unsigned int prevPosition = rpos_inclu_beg;
    for (unsigned int refpos = rpos_inclu_beg; refpos <= rpos_exclu_end; refpos++) {
        
if (rpos_exclu_end != refpos) {
        const std::array<SymbolType, 2> allSymbolTypes = {LINK_SYMBOL, BASE_SYMBOL};
        for (unsigned int stidx = 0; stidx < 2; stidx++) {
            const SymbolType symbolType = allSymbolTypes[stidx];
            bcfrec::BcfFormat init_fmt;
            std::array<unsigned int, 2> bDPcDP = BcfFormat_init(init_fmt, symbolToCountCoverageSet12, refpos, symbolType);
            AlignmentSymbol most_confident_symbol = END_ALIGNMENT_SYMBOLS;
            float most_confident_qual = 0;
            std::string most_confident_GT = "";
            float most_confident_GQ = 0;
            if (bDPcDP[0] >= paramset.min_depth_thres) {
                std::vector<bcfrec::BcfFormat> fmts(SYMBOL_TYPE_TO_INCLU_END[symbolType] - SYMBOL_TYPE_TO_INCLU_BEG[symbolType] + 1, init_fmt);
                for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; symbol = AlignmentSymbol(1+(unsigned int)symbol)) {
                    const auto simplemut = std::make_pair(refpos, symbol);
                    auto indices_bq = (simplemut2indices_bq.find(simplemut) != simplemut2indices_bq.end() ? simplemut2indices_bq[simplemut] : empty_size_t_set); 
                    auto indices_fq = (simplemut2indices_fq.find(simplemut) != simplemut2indices_fq.end() ? simplemut2indices_fq[simplemut] : empty_size_t_set);
                    
                    int altdepth = fillBySymbol(fmts[symbol - SYMBOL_TYPE_TO_INCLU_BEG[symbolType]], symbolToCountCoverageSet12, 
                            refpos, symbol, refstring, mutform2count4vec_bq, indices_bq, mutform2count4vec_fq, indices_fq, 
                            paramset.minABQ, paramset.phred_max_sscs, paramset.phred_max_dscs);
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
                for (AlignmentSymbol symbol = SYMBOL_TYPE_TO_INCLU_BEG[symbolType]; symbol <= SYMBOL_TYPE_TO_INCLU_END[symbolType]; symbol = AlignmentSymbol(1+(unsigned int)symbol)) {
                    auto & fmt = fmts[symbol - SYMBOL_TYPE_TO_INCLU_BEG[symbolType]];
                    if ((fmt.bAD1[0] + fmt.bAD1[1]) < paramset.min_altdp_thres) {
                        continue;
                    }
                    fmt.CType = SYMBOL_TO_DESC_ARR[most_confident_symbol];
                    fmt.CAQ = most_confident_qual;
                    appendVcfRecord(buf_out_string, buf_out_string_pass, symbolToCountCoverageSet12,
                            std::get<0>(tname_tseqlen_tuple).c_str(), refpos, symbol, fmt,
                            refstring, extended_inclu_beg_pos, paramset.vqual, should_output_all);
                }
            }
            bDPval[symbolType] = bDPcDP[0];
            cDPval[symbolType] = bDPcDP[1];
            fGTmm2[symbolType] = (most_confident_GQ >= 10 ? most_confident_GT : ".|.");
            fGQval[symbolType] = most_confident_GQ;
            // CAQval[symbolType] = most_confident_qual;
        }
}
        bool gbDPhasWideRange0 = (is_sig_out(bDPval[0], gbDPmin[0], gbDPmax[0], 130, 3) || 0 == gbDPmin[0]);
        bool gbDPhasWideRange1 = (is_sig_out(bDPval[1], gbDPmin[1], gbDPmax[1], 130, 3) || 0 == gbDPmin[1]);
        bool gcDPhasWideRange0 = (is_sig_out(cDPval[0], gcDPmin[0], gcDPmax[0], 130, 3) || 0 == gcDPmin[0]);
        bool gcDPhasWideRange1 = (is_sig_out(cDPval[1], gcDPmin[1], gcDPmax[1], 130, 3) || 0 == gcDPmin[1]);
        bool gfGQhasWideRange0 = (is_sig_out(fGQval[0], gfGQmin[0], gfGQmax[0], 130, 10) || std::string(fGTmm2[0]) != gfGTmm2[0]);
        bool gfGQhasWideRange1 = (is_sig_out(fGQval[1], gfGQmin[1], gfGQmax[1], 130, 10) || std::string(fGTmm2[1]) != gfGTmm2[1]);
        //bool gCAQhasWideRange0 = (is_sig_out(CAQval[0], gCAQmin[0], gCAQmax[1], 130, 10));
        //bool gCAQhasWideRange1 = (is_sig_out(CAQval[0], gCAQmin[0], gCAQmax[1], 130, 10));
        if ((   gbDPhasWideRange0 || gbDPhasWideRange1 ||
                gcDPhasWideRange0 || gcDPhasWideRange1 ||
                gfGQhasWideRange0 || gfGQhasWideRange1 ||
                // gCAQhasWideRange0 || gCAQhasWideRange1 ||
                (refpos - prevPosition >= G_BLOCK_SIZE) ||
                (refpos == rpos_exclu_end)) 
                // && refpos != rpos_inclu_beg
                ) {
            std::string genomicInfoString = ((0 == gbDPmax[0] || 0 == gbDPmax[1]) ? "" : genomicRegionInfoToString(
                    std::get<0>(tname_tseqlen_tuple),
                    prevPosition,
                    refpos,
                    gbDPmin,
                    gcDPmin,
                    // gCAQmin
                    gfGTmm2,
                    gfGQmin
                    ));
            raw_out_string += buf_out_string;
            raw_out_string_pass += genomicInfoString + buf_out_string_pass;
            buf_out_string.clear();
            buf_out_string_pass.clear();
            prevPosition = refpos;
            gbDPmin = bDPval;
            gbDPmax = bDPval;
            gcDPmin = cDPval;
            gcDPmax = cDPval;
            //gCAQmin = CAQval;
            //gCAQmax = CAQval;
            gfGTmm2 = fGTmm2;
            gfGQmin = fGQval;
            gfGQmax = fGQval;
        } else {
            UPDATE_MIN2(gbDPmin, bDPval);
            UPDATE_MAX2(gbDPmax, bDPval);
            UPDATE_MIN2(gcDPmin, cDPval);
            UPDATE_MAX2(gcDPmax, cDPval);
            //UPDATE_MIN2(gCAQmin, CAQval);
            //UPDATE_MAX2(gCAQmax, CAQval);
            UPDATE_MIN2(gfGQmin, fGQval);
            UPDATE_MAX2(gfGQmax, fGQval);
        }
    }
    LOG(logINFO) << "Thread " << thread_id  << " starts destroying bam records"; 
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
    LOG(logINFO) << "Thread " << thread_id  << " is done with current task";
};

int main(int argc, char **argv) {
    const char *UMI_STRUCT = getenv("ONE_STEP_UMI_STRUCT");
    const std::string UMI_STRUCT_STRING = (UMI_STRUCT != NULL ? std::string(UMI_STRUCT) : std::string(""));
    CommandLineArgs paramset;
    int parsing_result_flag = -1;
    int parsing_result_ret = paramset.initFromArgCV(parsing_result_flag, argc, argv);
    if (parsing_result_ret || parsing_result_flag) {
        return parsing_result_ret; 
    }
    LOG(logINFO) << "Program " << argv[0] << " version " << VERSION_DETAIL;
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, bool, unsigned int>> tid_beg_end_e2e_tuple_vec;
    std::vector<std::tuple<std::string, unsigned int>> tid_to_tname_tseqlen_tuple_vec;
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
    std::vector<hts_idx_t*> sam_idxs(nidxs, NULL);
    std::vector<samFile*> samfiles(nidxs, NULL);
    std::vector<faidx_t*> ref_faidxs(nidxs, NULL);
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
    }
    
    bam_hdr_t * samheader = sam_hdr_read(samfiles[0]);
    std::string header_outstring = generateVcfHeader(paramset.fasta_ref_fname.c_str(), argc, argv, 
            samheader->n_targets, samheader->target_name, samheader->target_len, 
            paramset.sample_name.c_str());
    clearstring<false>(fp_allp, header_outstring);
    clearstring<false>(fp_pass, header_outstring, is_vcf_out_pass_to_stdout);
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
                    
                    tid_beg_end_e2e_tuple : tid_beg_end_e2e_tuple_vec.at(0),
                    tname_tseqlen_tuple : tid_to_tname_tseqlen_tuple_vec.at(0),
                    region_ordinal : 0,
                    regionbatch_ordinal : 0,
                    
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

            std::pair<unsigned int, unsigned int> beg_end_pair = beg_end_pair_vec[beg_end_pair_idx];
#if defined(USE_STDLIB_THREAD)
            std::thread athread([
                        &batcharg, allridx, beg_end_pair, beg_end_pair_idx, &tid_beg_end_e2e_tuple_vec, &tid_to_tname_tseqlen_tuple_vec
                        //beg_end_pair_idx, allridx, is_vcf_out_pass_to_stdout, &beg_end_pair_vec, 
                        //&outstrings, &outstrings_pass, paramset, &tid_beg_end_e2e_tuple_vec, tid_to_tname_tseqlen_tuple_vec, UMI_STRUCT, 
                        //sam_idxs, samfiles, ref_faidxs
                        ]() {
#endif
                    LOG(logINFO) << "Thread " << batcharg.thread_id << " will process the sub-chunk " << beg_end_pair_idx << " which ranges from " 
                            << beg_end_pair.first << " to " << beg_end_pair.second;
                    
                    //auto *this_sam_idx = sam_index_load2(samfiles[allridx], paramset.bam_input_fname.c_str(), NULL);
                    for (unsigned int j = beg_end_pair.first; j < beg_end_pair.second; j++) {
                        batcharg.region_ordinal = allridx + j;
                        batcharg.regionbatch_ordinal = j;
                        batcharg.tid_beg_end_e2e_tuple = tid_beg_end_e2e_tuple_vec.at(allridx + j);
                        batcharg.tname_tseqlen_tuple = tid_to_tname_tseqlen_tuple_vec.at(std::get<0>(batcharg.tid_beg_end_e2e_tuple));
                        process_batch(batcharg
                                //outstrings[j], outstrings_pass[j], paramset, tid_beg_end_e2e_tuple_vec.at(i+j), tid_to_tname_tseqlen_tuple_vec, 
                                //i+j, j, UMI_STRUCT, this_sam_idx, //sam_idxs[thread_id], 
                                //ref_faidxs[thread_id], thread_id, is_vcf_out_pass_to_stdout
                        );
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
    }
    clearstring<true>(fp_allp, std::string("")); // write end of file
    clearstring<true>(fp_pass, std::string(""), is_vcf_out_pass_to_stdout);
    bam_hdr_destroy(samheader);
    for (size_t i = 0; i < nidxs; i++) { 
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
}

