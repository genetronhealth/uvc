#include <assert.h>
#include <stdint.h>
#include <map>

#include <iostream> // for debugging

// T=uint32_t for deletion and T=string for insertion

template <class T> 
uint32_t 
posToIndelToData_get(const std::map<uint32_t, std::map<T, uint32_t>> & pos2indel2data, const uint32_t refpos, const T & indel, uint32_t default1 = 0, uint32_t default2 = 0) {
    if (pos2indel2data.find(refpos) == pos2indel2data.end()) { return default1; }
    if (pos2indel2data.at(refpos).find(indel) == pos2indel2data.at(refpos).end()) { return default2; }
    return pos2indel2data.at(refpos).at(indel);
}

template <class T> 
void
posToIndelToCount_inc(std::map<uint32_t, std::map<T, uint32_t>> & pos2indel2count, uint32_t pos, const T indel, uint32_t incvalue = 1) {
    assert(incvalue > 0);
    auto pos2indel2count4it = pos2indel2count.insert(std::make_pair(pos, std::map<T, uint32_t>()));
    auto indel2count4it = pos2indel2count4it.first->second.insert(std::make_pair(indel, 0));
    indel2count4it.first->second += incvalue;
    assert (posToIndelToData_get(pos2indel2count, pos, indel) > 0); 
}

template <class T> 
T
posToIndelToCount_updateByConsensus(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src, uint32_t epos, uint32_t incvalue = 1) {
    const auto pos2indel2count4it = src.find(epos);
    assert(pos2indel2count4it != src.end());
    
    auto indel2count = pos2indel2count4it->second;
    // if (indel2count.size() > 1) { return T(); }
    
    const T & src_indel = (indel2count.size() > 1 ? T() : indel2count.begin()->first);
    // const uint32_t src_count = indel2count.begin()->second;
    assert (indel2count.begin()->second > 0);
    posToIndelToCount_inc<T>(dst, epos, src_indel, incvalue);
    return src_indel;
    //return (int)src_count;
}

template <bool TIsIncVariable, class T> 
int
posToIndelToCount_updateByRepresentative(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src, uint32_t epos, uint32_t incvalue = 1) {
    const auto pos2indel2count4it = src.find(epos);
    assert(pos2indel2count4it != src.end());
    
    T max_indel;
    uint32_t max_count = 0;
    uint32_t sum_count = 0;
    for (auto indel2count : pos2indel2count4it->second) {
        auto indel = indel2count.first;
        auto count = indel2count.second;
        if (count > max_count) {
            max_indel = indel;
            max_count = count;
        }
        if (TIsIncVariable) {sum_count += count; }
    }
    assert(0 < max_count);
    if (TIsIncVariable) {
        posToIndelToCount_inc<T>(dst, epos, max_indel, sum_count); 
    } else {
        posToIndelToCount_inc<T>(dst, epos, max_indel, incvalue);
    }
    return (int)max_count;
}

template <class T> 
void
posToIndelToCount_updateByConsensus_old(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src, uint32_t incvalue = 1) {
    for (auto src_pos2indel2count4it : src) {
        auto src_pos = src_pos2indel2count4it.first;
        auto src_indel2count = src_pos2indel2count4it.second;
        if (src_indel2count.size() == 1) {
            assert (src_indel2count.begin()->second >= incvalue 
                    || !fprintf(stderr, "Condition %d >= %d failed!\n", src_indel2count.begin()->second, incvalue)
                    || !(std::cerr << "value is " << src_indel2count.begin()->first << std::endl));
            posToIndelToCount_inc<T>(dst, src_pos, src_indel2count.begin()->first, incvalue); 
        } else if (src_indel2count.size() > 1) {
            // auto insresult2 = insresult->insert(std::make_pair(T(), 0)); insresult2-second++; // do nothing
        } else { assert(false); }
    }
}

template <class T> 
void
posToIndelToCount_updateBySummation(std::map<uint32_t, std::map<T, uint32_t>> & dst, const std::map<uint32_t, std::map<T, uint32_t>> & src) {
    for (auto src_pos2indel2count4it : src) {
        auto src_pos = src_pos2indel2count4it.first;
        auto src_indel2count = src_pos2indel2count4it.second;
        for (auto src_indel2count4it : src_indel2count) {
            assert(src_indel2count4it.second > 0 || !(
                std::cerr << src_indel2count4it.second << " > 0 failed for the key " << src_indel2count4it.first << " at position " << src_pos << std::endl
            ));
            posToIndelToCount_inc<T>(dst, src_pos, src_indel2count4it.first, src_indel2count4it.second);
        }
    }
}

void
posToIndelToCount_DlenToDseq(std::map<uint32_t, std::map<std::string, uint32_t>> & dst, const std::map<uint32_t, std::map<uint32_t, uint32_t>> & src,
        const std::string & refchars, uint32_t incluBegPos) {
    for (auto src_pos2indel2count4it : src) {
        uint32_t src_pos = src_pos2indel2count4it.first;
        auto src_del2count = src_pos2indel2count4it.second;
        for (auto src_del2count4it : src_del2count) {
            std::string dseq = refchars.substr(src_pos - incluBegPos, src_del2count4it.first);
            posToIndelToCount_inc(dst, src_pos, dseq, src_del2count4it.second);
        }
    }
}

