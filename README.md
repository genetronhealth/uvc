UVC is a very accurate and reasonably fast somatic variant caller.
The executable uvc1 in the bin direcotry takes one BAM file as input and generates one block-gzipped VCF file as output.
The script uvcTN.sh in the bin directory takes two BAM files corresponding to tumor and normal as input and generate two block-gzipped VCF files (tumor-variant VCF and normal-filtered VCF) as output.

# How to install

UVC requires BASH 4.0+ (4.0 is the minimum version required) and a compiler that supports the C++14 standard.
The Makefile in this directory compiles with g++, but the Makefile can be easily modified to use another compiler instead of g++ (for example, clang).
To install from scratch, please run: (./install-dependencies.sh && make clean && make all -j4 && make deploy). 
Please note that ./install-dependencies.sh requires bzip2 to decompress the downloaded files with the (.tar.bz2) extension.
UVC depends on git 2.12+, htslib 1.6+ and bcftools 1.6+ (lower versions of htslib and bcftools may also work, but are not tested).
If these two dependencies were already installed, then install-dependencies.sh may not be necessary.
For trouble-shooting with the installation of htslib and bcftools, please check their official repositories at https://github.com/samtools/htslib and https://github.com/samtools/bcftools.
More specifically, if any error message containing "error while loading shared libraries" pops up, please use the command (./configure --disable-plugins --disable-bz2 --disable-lzma --disable-libcurl --disable-s3 --disable-largefile) to build the corresponding htslib required by UVC first, then build UVC.
Although not required, it is highly recommmended that bcftools is installed at a system location (a location that can be found in the PATH environment variable).
The UVC binary uses multi-threading efficiently for up to 16 threads. 
After reaching 16 threads, adding more threads no longer significantly reduces wall-clock runtime.
However, more efficient speed-up can still be gained by runing with GNU parallel or qsub to use one job per chromosome.

In total, the installation should take about 5 minutes.

# How to use

The script uvcTN.sh in the bin directory is used for analyzing tumor-normal pairs.
Run uvcTN.sh without any command-line argument will display its usage help.
The usage help for uvcTN.sh refers to the executable uvc1, which performs the actual variant calling.
The executable uvc1 can perform each of the following tasks:
 1. tumor-only variant call to generate a tumor-only bgzipped vcf with vcf.gz as file extension.
 2. filtering of tumor variants in tumor-only bgzipped vcf with its matched normal.
The script uvcTN.sh simply wraps around the binary executable uvc1.

For UMI (unique molecular identifier, a.k.a. molecular barcode) to be detected, the read name (QNAME) in the input BAM file should be in the format of originalName#UMI.
For example, the UMI-labeled read name can be
 1. "H5G5ABBCC:4:1209:10114:63736#ACGTAACCA" (ACGTAACCA is the single-strand barcode) or 
 2. "H5G5ABBCC:1:3010:10412:33669#AGTA+TGGT" (AGTA+TGGT is the duplex barcode).
The auxiliary tool debarcode can be used for appending UMI sequences into read names.
Running debarcode without any command-line argument will display its usage help.

It is recommended to manually check the following outlier variant candidates if very high precision is required:
 1. for non-UMI data, variant candidates with FORMAT/FTS consisting of three or more filter strings (typically less than 1.5% of all variants).
 2. for UMI data, variant candidates with FORMAT/FTS consisting of one or more filter strings (typically less than 6% of all variants).

If manual check is still too labor-intensive, then it is recommended to keep such outlier variant candidate if the candidate
 1. is at a hotspot (for example, if the candidate shows high-frequency occurence in the COSMIC database) and
 2. does not show germline risk (such as low-frequency occurence or absence in dbSNP).

Otherwise, it is recommended to reject such variant candidate.

# What to report if a runtime error arises

In fact, uvc1 and some other executables all generated the same output given the same input. Their differences are as follows.
 1. uvc1: the release version that runs the fastest with multi-threading. 
    Because of extreme runtime optimization, this program probably will not generate any useful error message when an error arises.
 2. uvc.mt.out: the debug version that runs with multi-threading. 
    If the release version encounters any problem. please use this debug version to generate some useful error message. 
    The error message can then be used by the software maintainer for debugging and/or testing.
 3. uvc.st.out: the debug version that runs without multi-threading. 
    If the debug version with multi-threading encounters any problem, please use this debug version to generate some useful error message.
    The error message can then be used by the software maintainer for debugging and/or testing.
 4. uvc.cppt.out: similar to uvc1 except that uvc.cppt.out uses c++14 thread instead of OpenMP for multi-threading.

All bug reports, feature requests, and ideas for improvement are welcome (although not all of them may be addressed in time)!

# Other things

The environment variable ONE_STEP_UMI_STRUCT has special meaning to UVC.
Please make sure that ONE_STEP_UMI_STRUCT is either not set (e.g., by using the unset shell command) or set to the empty string before running UVC.
The python script extract-barcodes.py is obsolete and is replaced by debarcode.
Compared with extract-barcodes.py, debarcode generates equivalent output but consumes only 40% of its runtime.
The outputs of these two programs may not be the same in compressed space but are exactly the same in decompressed space.
The script bin/uvcnorm.sh can be used for normalizing variants.
By default, the normalization generates one SNV record per position and one InDel record per position.
The script bin/uvcSurrogateAlign.sh is still under development and should be be used.


For more information, please check the wiki.

