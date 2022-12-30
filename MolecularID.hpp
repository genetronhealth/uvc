#ifndef MolecularID_hpp_INCLUDED
#define MolecularID_hpp_INCLUDED

#include "common.hpp"

#include <string>

struct MolecularBarcode {

    std::pair<uvc1_refgpos_t, uvc1_refgpos_t> beg_tidpos_pair;
    std::pair<uvc1_refgpos_t, uvc1_refgpos_t> end_tidpos_pair;
    std::string qnamestring = "";
    std::string umistring = "";

    uvc1_flag_t duplexflag = 0x0;
    uvc1_flag_t dedup_idflag = 0x0;

    uvc1_hash_t hashvalue;

    MolecularBarcode
    createKey() const {
        MolecularBarcode mb;

        mb.beg_tidpos_pair = std::make_pair(-1, -1);
        mb.end_tidpos_pair = std::make_pair(-1, -1);
        if (0x3 == (0x3 & dedup_idflag)) {
            auto min2 = MIN(this->beg_tidpos_pair, this->end_tidpos_pair);
            auto max2 = MAX(this->beg_tidpos_pair, this->end_tidpos_pair);
            mb.beg_tidpos_pair = min2;
            mb.end_tidpos_pair = max2;
        } else if (0x1 & dedup_idflag) {
            mb.beg_tidpos_pair = this->beg_tidpos_pair;
        } else if (0x2 & dedup_idflag) {
            mb.end_tidpos_pair = this->end_tidpos_pair;
        }

        if (0x4 & dedup_idflag) {
            mb.qnamestring = this->qnamestring;
        } else {
            mb.qnamestring = "";
        }
        if (0x8 & dedup_idflag) {
            mb.umistring = this->umistring;
        } else {
            mb.umistring = "";
        }

        mb.duplexflag = this->duplexflag;
        mb.dedup_idflag = this->dedup_idflag;
        return mb;
    }
    bool
    operator<(const MolecularBarcode & that) const {
        bool isdiff, isless;
        compare_diff_less(isdiff, isless, this->beg_tidpos_pair, that.beg_tidpos_pair);
        if (isdiff) { return isless; }
        compare_diff_less(isdiff, isless, this->end_tidpos_pair, that.end_tidpos_pair);
        if (isdiff) { return isless; }
        compare_diff_less(isdiff, isless, this->qnamestring, that.qnamestring);
        if (isdiff) { return isless; }
        compare_diff_less(isdiff, isless, this->umistring, that.umistring);
        if (isdiff) { return isless; }
        compare_diff_less(isdiff, isless, this->duplexflag, that.duplexflag);
        if (isdiff) { return isless; }
        compare_diff_less(isdiff, isless, this->dedup_idflag, that.dedup_idflag);
        if (isdiff) { return isless; }
        return (this->hashvalue < that.hashvalue);
    }
    uvc1_hash_t calcHash() const;
};

#endif
