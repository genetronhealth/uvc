* UVC

UVC is a very accurate and reasonably fast somatic variant caller.
The executable uvc1 in the bin directory takes one BAM file as input and generates one block-gzipped VCF file as output.
The script uvcTN.sh in the bin directory takes two BAM files corresponding to tumor and normal as input and generate two block-gzipped VCF files (tumor-variant VCF and normal-filtered VCF) as output.

# how to install
1. install libhts library
see https://github.com/samtools/htslib

2. compile the uvc source codes with g++ (supports the C++14 standard) and cmake (>= 3.4)
```
git clone https://github.com/Schaudge/uvc.git 
cd uvc && mkdir build && cd build
cmake.. && make
```
We does not recommend compile (and install) the software from Makefile by the make command currently, unless you known what you mean!

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


# Other things

The environment variable ONE_STEP_UMI_STRUCT has special meaning to UVC.
Please make sure that ONE_STEP_UMI_STRUCT is either not set (e.g., by using the unset shell command) or set to the empty string before running UVC.
The python script extract-barcodes.py is obsolete and is replaced by debarcode.
Compared with extract-barcodes.py, debarcode generates equivalent output but consumes only 40% of its runtime.
The outputs of these two programs may not be the same in compressed space but are exactly the same in decompressed space.
The script bin/uvcnorm.sh can be used for normalizing variants.
By default, the normalization generates one SNV record per position and one InDel record per position.

