#include "MolecularID.hpp"

#include "Hash.hpp"

uvc1_hash_t MolecularBarcode::calcHash() const {
        uvc1_hash_t ret = 0;
        ret += hash2hash(ret, hash2hash(beg_tidpos_pair.first, beg_tidpos_pair.second));
        ret += hash2hash(ret, hash2hash(end_tidpos_pair.first, end_tidpos_pair.second));
        ret += hash2hash(ret, strhash(qnamestring.c_str()));
        ret += hash2hash(ret, strhash(umistring.c_str()));
        ret += hash2hash(ret, duplexflag);
        ret += hash2hash(ret, dedup_idflag);
        return ret;
}

