
#undef INDELTYPE
#if INDEL_ID == 1
#define INDELTYPE std::string
#else
#define INDELTYPE uint32_t
#endif

// this code is instantiated multiple times, with INDELTYPE as instantiation parameter.
// template <AlignmentSymbol link1, AlignmentSymbol link2, AlignmentSymbol link3p> // , LinkType linkType>
void   
#if INDEL_ID == 1
fillByIndelInfo2_1
#else
fillByIndelInfo2_2
#endif
(bcfrec::BcfFormat & fmt,
        const Symbol2CountCoverageSet & symbol2CountCoverageSet, 
        const unsigned int strand, const unsigned int refpos, const AlignmentSymbol symbol,
        
        //const std::map<std::uint32_t, std::map<INDELTYPE, uint32_t>> & rawdu_amplicon_pos2indel2data,
        //const std::map<std::uint32_t, std::map<INDELTYPE, uint32_t>> & rawdu_ampBQerr_pos2indel2data,
        //const std::map<std::uint32_t, std::map<INDELTYPE, uint32_t>> & dedup_amplicon_pos2indel2data,
        //const std::map<std::uint32_t, std::map<INDELTYPE, uint32_t>> & size1_amplicon_pos2indel2data, 
        const std::map<std::uint32_t, std::map<INDELTYPE, uint32_t>> & bq_tsum_depth, //  _ amplicon_pos2indel2data,
        const std::map<std::uint32_t, std::map<INDELTYPE, uint32_t>> & fq_tsum_depth, // size1_amplicon_pos2indel2data, 
        const std::string & refchars) {
    
    assert(isSymbolIns(symbol) || isSymbolDel(symbol));
    AlignmentSymbol link1  = (isSymbolIns(symbol) ? LINK_I1  : LINK_D1 ); // , AlignmentSymbol link2, AlignmentSymbol link3p;
    AlignmentSymbol link2  = (isSymbolIns(symbol) ? LINK_I2  : LINK_D2 ); 
    AlignmentSymbol link3p = (isSymbolIns(symbol) ? LINK_I3P : LINK_D3P);
    
    assert (link1 == symbol || link2 == symbol || link3p == symbol || 
            !fprintf(stderr, "Symbol %d does not match any of {%s, %s, %s}", 
            SYMBOL_TO_DESC_ARR[symbol], SYMBOL_TO_DESC_ARR[link1], SYMBOL_TO_DESC_ARR[link2], SYMBOL_TO_DESC_ARR[link3p]));
    // std::vector<std::tuple<uint32_t, uint32_t, std::string>> rawdu2_dedup_size1_mutform_tuples; 
    
    std::vector<std::tuple<uint32_t, uint32_t, std::string>> bqfq_depth_mutform_tuples; // bqfq_tsum_tuples;
    // assert(rawdu_amplicon_pos2indel2data.find(refpos) != rawdu_amplicon_pos2indel2data.end());
    assert(bq_tsum_depth.find(refpos) != bq_tsum_depth.end());
    // assert(fq_tsum_depth.find(refpos) != fq_tsum_depth.end()); // may not find because of filtering
    //assert(_pos2indel2data.find(refpos) != rawdu_ampBQerr_pos2indel2data.end());
    //assert(rawdu_ampBQerr_pos2indel2data.find(refpos) != rawdu_ampBQerr_pos2indel2data.end());
    //assert(rawdu_ampBQerr_pos2indel2data.find(refpos) != rawdu_ampBQerr_pos2indel2data.end());
    
    // assert(rawdu_amplicon_pos2indel2data.at(refpos).size() == rawdu_ampBQerr_pos2indel2data.at(refpos).size());
    // for (auto indel2data4it = rawdu_amplicon_pos2indel2data.at(refpos).begin(); indel2data4it != rawdu_amplicon_pos2indel2data.at(refpos).end(); indel2data4it++) {
    for (auto indel2data4 : bq_tsum_depth.at(refpos)) { // = bq_tsum_depth.at(refpos).begin(); indel2data4it != bq_tsum_depth.at(refpos).end(); indel2data4it++) {
        const auto indel = indel2data4.first;
#if INDEL_ID == 1
        //if (INS_LINK == linkType) {
        const std::string indelstring = indel2data4.first;
#else
        //} else {
        const std::string indelstring = refchars.substr(refpos - symbol2CountCoverageSet.duplex_tsum_depth.getIncluBegPosition(), indel2data4.first); 
        //}
#endif
        //const uint32_t rawdu_amplicon_data = indel2data4.second;
        assert (indelstring.size() > 0);
        //assert (rawdu_amplicon_data > 0);
        //assert (rawdu_ampBQerr_pos2indel2data.at(refpos).find(indel) != rawdu_ampBQerr_pos2indel2data.at(refpos).end());
        
        const uint32_t bqdata = posToIndelToData_get(bq_tsum_depth, refpos, indel);
        const uint32_t fqdata = posToIndelToData_get(fq_tsum_depth, refpos, indel);
        assert(bqdata > 0);
        //if (link1 == symbol && indelstring.size() != 1 || link2 == symbol && indelstring.size() != 2 || link3p == symbol && indelstring.size() < 3) {
        //    continue;
        //}
        if (link1 == symbol && indelstring.size() == 1 || link2 == symbol && indelstring.size() == 2 || link3p == symbol && indelstring.size() >= 3) {
            // rawdu2_dedup_size1_mutform_tuples.push_back(
            //        std::make_tuple(rawdu_amplicon_data, rawdu_ampBQerr_data, dedup_amplicon_data, size1_amplicon_data, indelstring)
            // );
            bqfq_depth_mutform_tuples.push_back(std::make_tuple(bqdata, fqdata, indelstring));
        }
    }
    unsigned int gapT1AD1sum = 0;
    unsigned int gapT2AD1sum = 0;
    // std::sort(rawdu2_dedup_size1_mutform_tuples.rbegin(), rawdu2_dedup_size1_mutform_tuples.rend());
    // bq_tsum_depth.push_back(std::make_tuple(bq_tsum_depth, fq_tsum_depth, indelstring));
    std::sort(bqfq_depth_mutform_tuples.rbegin(), bqfq_depth_mutform_tuples.rend());
    // size_t ituple = 0;
    fmt.gapNum[strand] = bqfq_depth_mutform_tuples.size() ; // rawdu2_dedup_size1_mutform_tuples.size();
    //unsigned int gapAD2sum = 0;
    //unsigned int gapADrsum = 0;
    //for (auto rawdu_dedup_size1_mutform_tuple : rawdu2_dedup_size1_mutform_tuples) {
    for (auto bqfq_depth_mutform : bqfq_depth_mutform_tuples) {
        assert(std::get<2>(bqfq_depth_mutform).size() > 0);
        fmt.gapSeq.push_back(std::get<2>(bqfq_depth_mutform));
        fmt.gapT1AD1.push_back(std::get<0>(bqfq_depth_mutform));
        fmt.gapT2AD1.push_back(std::get<1>(bqfq_depth_mutform));
        //fmt.gapSeq.push_back(std::get<4>(rawdu_dedup_size1_mutform_tuple)); // this is a std::string
        //fmt.gapAD2.push_back(std::get<2>(rawdu_dedup_size1_mutform_tuple));
        //fmt.gapAD4.push_back(std::get<3>(rawdu_dedup_size1_mutform_tuple));
        //fmt.gapADe.push_back(std::get<1>(rawdu_dedup_size1_mutform_tuple));
        //fmt.gapADr.push_back(std::get<0>(rawdu_dedup_size1_mutform_tuple));
        //gapAD2sum += fmt.gapAD2.back();
        //gapADrsum += fmt.gapADr.back();
        gapT1AD1sum += fmt.gapT1AD1.back();
        gapT2AD1sum += fmt.gapT2AD1.back();
        // ituple++; if (2 == ituple) { break; }
    }
    if (fmt.gapT1AD1[strand] != gapT1AD1sum) {
        std::string msg = std::to_string(strand) + "\t" + std::to_string(refpos) + "\t" + std::to_string(symbol);
        bcfrec::streamAppendBcfFormat(msg, fmt);
        std::cerr << msg << "\n";
    }
    // "/4+16" is a probabilistic check
    /*
    assert(fmt.AD2[strand] >= gapAD2sum && (fmt.AD2[strand] <= gapAD2sum * 5 / 4 + 16) ||
            !(std::cerr << fmt.AD2[strand] << " >= | <=5/4+16 " << gapAD2sum
            << " failed for AD2 and gapAD2sum for strand " << strand << " at position " << refpos << " for symbol " << SYMBOL_TO_DESC_ARR[symbol] 
            << " and gapNum " << fmt.gapNum[strand] << " or equiv " << rawdu2_dedup_size1_mutform_tuples.size() << std::endl));
    assert(fmt.ADr[strand] >= gapADrsum && (fmt.ADr[strand] <= gapADrsum * 5 / 4 + 16) ||
            !(std::cerr << fmt.ADr[strand] << " >= | <=5/4+16 " << gapADrsum
            << " failed for ADr and gapADrsum for strand " << strand << " at position " << refpos << " for symbol " << SYMBOL_TO_DESC_ARR[symbol] << " and gapNum " << fmt.gapNum[strand] << std::endl));
    */
};

