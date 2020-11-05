#include "common.hpp"

const std::vector<std::string> BIAS_TYPE_TO_MSG = {
    "Fragment duplication bias with bitmask " + std::to_string(BIAS_FRAG_DUP),
    "Fragment position bias with bitmask " + std::to_string(BIAS_FRAG_POS),
    "Fragment strand bias with bitmask " + std::to_string(BIAS_FRAG_STR),
    "Fragment mismatch bias with bitmas " + std::to_string(BIAS_FRAG_MIS),
    "Sequencing-segment position bias with bitmask " + std::to_string(BIAS_SSEG_POS),
    "Sequencing-segment strand bias with bitmask "  + std::to_string(BIAS_SSEG_STR),
    "Sequencing-segment two-sided position bias with bitmask "  + std::to_string(BIAS_SSEG_END)
};

const std::vector<std::string> ASSAY_TYPE_TO_MSG = {
    [ASSAY_TYPE_AUTO] = "Automatically infer assay type from the data",
    [ASSAY_TYPE_CAPTURE] = "Data is generatd from a capture-based assay",
    [ASSAY_TYPE_AMPLICON] = "Data is generated from an amplicon-based assay",
};

const std::vector<std::string> MOLECULE_TAG_TO_MSG = {
    [MOLECULE_TAG_AUTO] = "Automatically infer assay type from the data",
    [MOLECULE_TAG_NONE] = "Molecule is not tagged",
    [MOLECULE_TAG_BARCODING] = "Molecule is tagged with an unique molecular identifer on one strand as in SAFE-SEQ",
    [MOLECULE_TAG_DUPLEX] = "Molecule is tagged with a duplex tag",
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





