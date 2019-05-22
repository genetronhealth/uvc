#ifndef common_INCLUDED
#define common_INCLUDED
#define NUM_BUCKETS (64)

#include <string>
#include <vector>

const char *const NOT_PROVIDED = "";

enum AssayType {
    ASSAY_TYPE_AUTO,
    ASSAY_TYPE_CAPTURE,
    ASSAY_TYPE_AMPLICON,
};

const std::vector<std::string> ASSAY_TYPE_TO_MSG = {
    [ASSAY_TYPE_AUTO] = "Automatically infer assay type from the data (根据数据自动判断试验类型)",
    [ASSAY_TYPE_CAPTURE] = "Data is generatd from a capture-based assay (捕获试验)",
    [ASSAY_TYPE_AMPLICON] = "Data is generated from an amplicon-based assay  (扩增子试验)" 
};

enum MoleculeTag {
    MOLECULE_TAG_AUTO,
    MOLECULE_TAG_NONE,
    MOLECULE_TAG_BARCODING,
    MOLECULE_TAG_DUPLEX,
};

const std::vector<std::string> MOLECULE_TAG_TO_MSG = {
    [MOLECULE_TAG_AUTO] = "Automatically infer assay type from the data (根据数据自动判断试验类型)",
    [MOLECULE_TAG_NONE] = "Molecule is not tagged (没有分子标签)",
    [MOLECULE_TAG_BARCODING] = "Molecule is tagged with an unique molecular identifer on one strand as in SAFE-SEQ (有单链分子标签)", 
    [MOLECULE_TAG_DUPLEX] = "Molecule is tagged with a duplex tag (有双链分子标签)",
};

enum SequencingPlatform {
    SEQUENCING_PLATFORM_AUTO,
    SEQUENCING_PLATFORM_ILLUMINA,
    SEQUENCING_PLATFORM_IONTORRENT,
    SEQUENCING_PLATFORM_OTHER,
};

const std::vector<std::string> SEQUENCING_PLATFORM_TO_MSG = {
    [SEQUENCING_PLATFORM_AUTO] = "Automatically infer assay type from the data (根据数据自动判断试验类型)",
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
    [PAIR_END_MERGE_YES] = "paired-end sequenced segments are merged (合并R1和R2)",
    [PAIR_END_MERGE_NO]  = "paired-end sequenced segments are not merged (不合并R1和R2)",
};

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

const char *const CORRECTION_TYPE_TO_MSG[] = {
    [CORRECTION_AUTO] = "Automatically infer from the input BAM file, it is recommended to use this mode (自动识别去重找错机制，建议使用此模式)",
    [CORRECTION_NONE] = "Do not attempt to find and/or correct errors, consider R1 and R2 as two reads (不去重不找错)",
    [CORRECTION_BASEQUAL] = "Find and/or correct errors by using base quality (用碱基质量找错)",
    [CORRECTION_SINGLETON] = "PCR without UMI which is experimental (还在测试阶段的方法)",
    [CORRECTION_DUPLICATE] = "Find and/or correct errors by dedupping with the density distribution of start/end coordinates of mapped reads (用起始终止的位置去重找错)",
    [CORRECTION_BARCODE] = "Find and/or correct errors by dedupping with UMI (unique molecular identifier) signle-strand barcode and the denstiy distribution of start/end coordinates of mapped reads "
            "(用分子签和起始终止位置的分布去重找错)",
    [CORRECTION_DUPLEX] = "Find and/or correct errors by dedupping with duplex barcode and the denstiy distribution of start/end coordinates of mapped reads (用duplex标签和起始终止位置的分布去重找错)",
};
*/

bool ispowerof2(auto num) {
    return (num & (num-1)) == 0;
}

#endif

