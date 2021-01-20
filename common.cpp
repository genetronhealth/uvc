#include "common.hpp"

const std::vector<std::string> ASSAY_TYPE_TO_MSG = {{
    [ASSAY_TYPE_AUTO] = "Assay type of each molecule fragment will be automatically inferred from the data",
    [ASSAY_TYPE_CAPTURE] = "Data is generatd from a capture-based assay with selection by probe hybridization",
    [ASSAY_TYPE_AMPLICON] = "Data is generated from an amplicon-based assay with targeted amplification by PCR",
}};

const std::vector<std::string> MOLECULE_TAG_TO_MSG = {{
    [MOLECULE_TAG_AUTO] = "Molecule tag of each molecule fragment will be automatically inferred from the data",
    [MOLECULE_TAG_NONE] = "Molecule is not tagged",
    [MOLECULE_TAG_BARCODING] = "Molecule is tagged with a unique molecular identifer (UMI) on one strand as in Safe-SeqS",
    [MOLECULE_TAG_DUPLEX] = "Molecule is tagged with a duplex UMI",
}};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_MSG = {{
    [SEQUENCING_PLATFORM_AUTO] = "Unknown sequencing platform that will be automatically inferred from the data",
    [SEQUENCING_PLATFORM_ILLUMINA] = "Illumina sequencing platform (compatible with BGI and MGI)",
    [SEQUENCING_PLATFORM_IONTORRENT] = "IonTorrent sequencing platform by Life Technologies and ThermoFisher",
    [SEQUENCING_PLATFORM_OTHER] = "Other sequencing platform (for example, Nanopore)",
}};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_NAME = {{
    [SEQUENCING_PLATFORM_AUTO] = "AUTO",
    [SEQUENCING_PLATFORM_ILLUMINA] = PLAT_ILLUMINA_LIKE,
    [SEQUENCING_PLATFORM_IONTORRENT] = PLAT_ION_LIKE,
    [SEQUENCING_PLATFORM_OTHER] = "OtherSequencingPlatform",
}};

const std::vector<std::string> PAIR_END_MERGE_TO_MSG = {{
    [PAIR_END_MERGE_YES] = "paired-end sequenced segments are merged",
    [PAIR_END_MERGE_NO]  = "paired-end sequenced segments are not merged",
}};

