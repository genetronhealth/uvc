#ifndef common_INCLUDED
#define common_INCLUDED
#define NUM_BUCKETS (64)

const char *const NOT_PROVIDED = "";

enum LibPrepType {
    LIB_PREP_AUTO,
    LIB_PREP_CAPTURE,
    LIB_PREP_PCR,
};

enum SequencingType {
    SEQ_TYPE_AUTO,
    SEQ_TYPE_NORMAL,
    SEQ_TYPE_BARCODE,
    SEQ_TYPE_DUPLEX,
};

enum SequencingPlatform {
    SEQ_PLAT_AUTO,
    SEQ_PLAT_ILLUMINA,
    SEQ_PLAT_IONTORRENT,
};

enum ReadEnd {
    READ_END_AUTO,
    READ_END_SINGLE,
    READ_END_PAIR,
};

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

bool ispowerof2(auto num) {
    return (num & (num-1)) == 0;
}

#endif

