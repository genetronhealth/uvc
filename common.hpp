#ifndef common_hpp_INCLUDED
#define common_hpp_INCLUDED

#include <array>
#include <string>
#include <vector>

#define NUM_BUCKETS (64+8)
#define SIGN2UNSIGN(a) ((unsigned int)(a))

#define VCFQUAL_NUM_BINS 120 // max vcfqual is at 105, so vcfqual that is higher than the max is reserved
#define RCC_NUM 2
#define RCC_NFS 6

enum BiasType {
    BIAS_FRAG_DUP = 1,
    BIAS_FRAG_POS = 2,
    BIAS_FRAG_STR = 4,
    BIAS_FRAG_MIS = 8,
    BIAS_SSEG_POS = 32,
    BIAS_SSEG_STR = 64,
    BIAS_SSEG_END = 128,
};

const std::vector<std::string> BIAS_TYPE_TO_MSG = {
    "Fragment duplication bias with bitmask " + std::to_string(BIAS_FRAG_DUP),
    "Fragment position bias with bitmask " + std::to_string(BIAS_FRAG_POS),
    "Fragment strand bias with bitmask " + std::to_string(BIAS_FRAG_STR),
    "Fragment mismatch bias with bitmas " + std::to_string(BIAS_FRAG_MIS),
    "Sequencing-segment position bias with bitmask " + std::to_string(BIAS_SSEG_POS),
    "Sequencing-segment strand bias with bitmask "  + std::to_string(BIAS_SSEG_STR),
    "Sequencing-segment two-sided position bias with bitmask "  + std::to_string(BIAS_SSEG_END)
};

const char *const NOT_PROVIDED = "";

enum AssayType {
    ASSAY_TYPE_AUTO,
    ASSAY_TYPE_CAPTURE,
    ASSAY_TYPE_AMPLICON,
};

const std::vector<std::string> ASSAY_TYPE_TO_MSG = {
    [ASSAY_TYPE_AUTO] = "Automatically infer assay type from the data",
    [ASSAY_TYPE_CAPTURE] = "Data is generatd from a capture-based assay",
    [ASSAY_TYPE_AMPLICON] = "Data is generated from an amplicon-based assay", 
};

enum MoleculeTag {
    MOLECULE_TAG_AUTO,
    MOLECULE_TAG_NONE,
    MOLECULE_TAG_BARCODING,
    MOLECULE_TAG_DUPLEX,
};

const std::vector<std::string> MOLECULE_TAG_TO_MSG = {
    [MOLECULE_TAG_AUTO] = "Automatically infer assay type from the data",
    [MOLECULE_TAG_NONE] = "Molecule is not tagged",
    [MOLECULE_TAG_BARCODING] = "Molecule is tagged with an unique molecular identifer on one strand as in SAFE-SEQ", 
    [MOLECULE_TAG_DUPLEX] = "Molecule is tagged with a duplex tag",
};

enum SequencingPlatform {
    SEQUENCING_PLATFORM_AUTO,
    SEQUENCING_PLATFORM_ILLUMINA,
    SEQUENCING_PLATFORM_IONTORRENT,
    SEQUENCING_PLATFORM_OTHER,
};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_MSG = {
    [SEQUENCING_PLATFORM_AUTO] = "Automatically infer assay type from the data",
    [SEQUENCING_PLATFORM_ILLUMINA] = "Illumina sequencing platform (compatible with BGI)",
    [SEQUENCING_PLATFORM_IONTORRENT] = "IonTorrent sequencing platform by Life Technologies",
    [SEQUENCING_PLATFORM_OTHER] = "Other sequencing platform (for example, Nanopore)",
};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_DESC = {
    [SEQUENCING_PLATFORM_AUTO] = "AUTO",
    [SEQUENCING_PLATFORM_ILLUMINA] = "ILLUMINA/BGI",
    [SEQUENCING_PLATFORM_IONTORRENT] = "IONTORRENT/LifeTechnologies",
    [SEQUENCING_PLATFORM_OTHER] = "OTHER",
};

enum PairEndMerge {
    PAIR_END_MERGE_YES,
    PAIR_END_MERGE_NO,
};

const std::vector<std::string> PAIR_END_MERGE_TO_MSG = {
    [PAIR_END_MERGE_YES] = "paired-end sequenced segments are merged",
    [PAIR_END_MERGE_NO]  = "paired-end sequenced segments are not merged",
};


const std::array<std::string, 2> GT_HETERO = {{"0/1", "0/."}};
const std::array<std::string, 2> GT_HOMREF = {{"0/0", "0/0"}};
const std::array<std::string, 2> GT_HOMALT = {{"1/1", "./."}};

const std::array<std::string, 2> TT_HETERO = {{"./1", "1/."}};
const std::array<std::string, 2> TT_HOMREF = {{"./1", "1/."}};
const std::array<std::string, 2> TT_HOMALT = {{"1/1", "1/."}};

struct TnDP4 {
    unsigned int nvars= 0;
    unsigned int tuAD = 0;
    unsigned int noAD = 0;
    unsigned int tuDP = 0;
    unsigned int noDP = 0;
    TnDP4() {}
};

struct VcStats {
    std::array<TnDP4, VCFQUAL_NUM_BINS> vcfqual_to_count = {{TnDP4()}};
    
    VcStats() {}
    
    void 
    update(const VcStats & other) {
        for (unsigned int i = 0; i < VCFQUAL_NUM_BINS; i++) { 
            this->vcfqual_to_count[i].nvars+= other.vcfqual_to_count[i].nvars;
            this->vcfqual_to_count[i].tuAD += other.vcfqual_to_count[i].tuAD;
            this->vcfqual_to_count[i].noAD += other.vcfqual_to_count[i].noAD;
            this->vcfqual_to_count[i].tuDP += other.vcfqual_to_count[i].tuDP;
            this->vcfqual_to_count[i].noDP += other.vcfqual_to_count[i].noDP;
        }
    };
    
    int
    write_tsv(auto & ostream) {
        ostream <<"##Mass-function-aka-density-histogram\n";
        ostream <<"#Variant-quality\tnumber-of-calls\ttumor-AD\tnormal-AD\ttumor-DP\tnormal-DP\n";
        for (unsigned int i = 0; i < VCFQUAL_NUM_BINS; i++) { 
            ostream << i << "\t"; 
            ostream << this->vcfqual_to_count[i].nvars<< "\t";
            ostream << this->vcfqual_to_count[i].tuAD << "\t";
            ostream << this->vcfqual_to_count[i].noAD << "\t";
            ostream << this->vcfqual_to_count[i].tuDP << "\t";
            ostream << this->vcfqual_to_count[i].noDP << "\n";
        }
        unsigned int tot_nvars = 0;
        unsigned int tot_tuAD = 0;
        unsigned int tot_noAD = 0;
        unsigned int tot_tuDP = 0;
        unsigned int tot_noDP = 0;
        bool vcfqual_is_too_low = false;
        ostream <<"##Complementary-cumulative-function\n";
        ostream <<"#variant-quality\tnumber-of-calls\ttumor-AD\tnormal-AD\ttumor-DP\tnormal-DP"
                << "\testimated-additive-contamination-fraction\testimated-multiplicative-contamination-fraction"
                << "\tTUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION_STATUS_FOR_VCF_QUAL\n";
        for (int i = VCFQUAL_NUM_BINS - 1; i >= 0; i--) { 
            tot_nvars += this->vcfqual_to_count[i].nvars;
            tot_tuAD += this->vcfqual_to_count[i].tuAD;
            tot_noAD += this->vcfqual_to_count[i].noAD;
            tot_tuDP += this->vcfqual_to_count[i].tuDP;
            tot_noDP += this->vcfqual_to_count[i].noDP;
            
            ostream << i << "\t"; 
            ostream << tot_nvars<< "\t";
            ostream << tot_tuAD << "\t";
            ostream << tot_noAD << "\t";
            ostream << tot_tuDP << "\t";
            ostream << tot_noDP << "\t";

            ostream << (0.03 + (double)tot_noAD)                  / (1.0 + (double)tot_tuAD                 ) << "\t"; 
            ostream << (0.03 + (double)tot_noAD*(double)tot_tuDP) / (1.0 + (double)tot_tuAD*(double)tot_noDP) << "\t"; 
            
            if ((!vcfqual_is_too_low) && tot_nvars >= 25 && tot_noDP >= 25*100 && tot_tuDP >= 25*100) {
                ostream << "VCF_QUAL_EQUALS_THE_THRESHOLD_FOR_TUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION" << "\n";  
                vcfqual_is_too_low = true;
            } else if (vcfqual_is_too_low) {
                ostream << "VCF_QUAL_IS_BELOW_THE_THRESHOLD_FOR_TUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION" << "\n"; 
            } else {
                ostream << "VCF_QUAL_IS_ABOVE_THE_THRESHOLD_FOR_TUMOR_TO_NORMAL_CONTAMINATION_ESTIMATION" << "\n";  
            }
        }
        return 0;
    }
};

// region_pos32_unitlen8_repeatnum16_qual8_vec
struct RegionTandemRepeat {
    uint32_t begpos = 0;
    uint16_t tracklen = 0;
    uint8_t unitlen = 0;
    // uint8_t edgeBAQ;
};

// I tried to use the following correction types in the past. 
// However, correction types were very confusing to the end user and always missed some edge cases.
// Hence, I used assay types, sequencing platforms, UMI-awareness, etc. instead of correction types.
// If anyone thinks that correction types can be reworked into something intuitive to the end user, please let me know.
/*
enum ErrorCorrectionType {
    CORRECTION_AUTO,
    CORRECTION_NONE,
    CORRECTION_BASEQUAL,
    CORRECTION_SINGLETON,
    CORRECTION_DUPLICATE,
    CORRECTION_BARCODE,
    CORRECTION_DUPLEX,
    END_ERROR_CORRECTION_TYPES
};
*/

bool ispowerof2(auto num) {
    return (num & (num-1)) == 0;
}

#endif

