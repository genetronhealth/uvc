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
    [ASSAY_TYPE_AUTO] = "Assay type of each molecule fragment will be automatically inferred from the data",
    [ASSAY_TYPE_CAPTURE] = "Data is generatd from a capture-based assay with selection by probe hybridization",
    [ASSAY_TYPE_AMPLICON] = "Data is generated from an amplicon-based assay with targeted amplification by PCR",
};

const std::vector<std::string> MOLECULE_TAG_TO_MSG = {
    [MOLECULE_TAG_AUTO] = "Molecule tag of each molecule fragment will be automatically inferred from the data",
    [MOLECULE_TAG_NONE] = "Molecule is not tagged",
    [MOLECULE_TAG_BARCODING] = "Molecule is tagged with a unique molecular identifer (UMI) on one strand as in Safe-SeqS",
    [MOLECULE_TAG_DUPLEX] = "Molecule is tagged with a duplex UMI",
};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_MSG = {
    [SEQUENCING_PLATFORM_AUTO] = "Unknown sequencing platform that will be automatically inferred from the data",
    [SEQUENCING_PLATFORM_ILLUMINA] = "Illumina sequencing platform (compatible with BGI and MGI)",
    [SEQUENCING_PLATFORM_IONTORRENT] = "IonTorrent sequencing platform by Life Technologies and ThermoFisher",
    [SEQUENCING_PLATFORM_OTHER] = "Other sequencing platform (for example, Nanopore)",
};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_NAME = {
    [SEQUENCING_PLATFORM_AUTO] = "AUTO",
    [SEQUENCING_PLATFORM_ILLUMINA] = PLAT_ILLUMINA_LIKE,
    [SEQUENCING_PLATFORM_IONTORRENT] = PLAT_ION_LIKE,
    [SEQUENCING_PLATFORM_OTHER] = "OtherSequencingPlatform",
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





