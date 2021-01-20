#include "common.hpp"

const std::vector<std::string> ASSAY_TYPE_TO_MSG = {
    [ASSAY_TYPE_AUTO] = std::string("Assay type of each molecule fragment will be automatically inferred from the data"),
    [ASSAY_TYPE_CAPTURE] = std::string("Data is generatd from a capture-based assay with selection by probe hybridization"),
    [ASSAY_TYPE_AMPLICON] = std::string("Data is generated from an amplicon-based assay with targeted amplification by PCR"),
};

const std::vector<std::string> MOLECULE_TAG_TO_MSG = {
    [MOLECULE_TAG_AUTO] = std::string("Molecule tag of each molecule fragment will be automatically inferred from the data"),
    [MOLECULE_TAG_NONE] = std::string("Molecule is not tagged"),
    [MOLECULE_TAG_BARCODING] = std::string("Molecule is tagged with a unique molecular identifer (UMI) on one strand as in Safe-SeqS"),
    [MOLECULE_TAG_DUPLEX] = std::string("Molecule is tagged with a duplex UMI"),
};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_MSG = {
    [SEQUENCING_PLATFORM_AUTO] = std::string("Unknown sequencing platform that will be automatically inferred from the data"),
    [SEQUENCING_PLATFORM_ILLUMINA] = std::string("Illumina sequencing platform (compatible with BGI and MGI)"),
    [SEQUENCING_PLATFORM_IONTORRENT] = std::string("IonTorrent sequencing platform by Life Technologies and ThermoFisher"),
    [SEQUENCING_PLATFORM_OTHER] = std::string("Other sequencing platform (for example, Nanopore)"),
};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_NAME = {
    [SEQUENCING_PLATFORM_AUTO] = std::string("AUTO"),
    [SEQUENCING_PLATFORM_ILLUMINA] = std::string(PLAT_ILLUMINA_LIKE),
    [SEQUENCING_PLATFORM_IONTORRENT] = std::string(PLAT_ION_LIKE),
    [SEQUENCING_PLATFORM_OTHER] = std::string("OtherSequencingPlatform"),
};

const std::vector<std::string> PAIR_END_MERGE_TO_MSG = {
    [PAIR_END_MERGE_YES] = std::string("paired-end sequenced segments are merged"),
    [PAIR_END_MERGE_NO]  = std::string("paired-end sequenced segments are not merged"),
};
