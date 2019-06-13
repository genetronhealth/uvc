#!/usr/bin/env python

import subprocess
import sys

def recall_to_precision_curve(binseq):
    rp = []
    npositives = sum(binseq)
    nposs = 0
    for i, label in enumerate(binseq):
        assert label == 0 or label == 1
        nposs += label
        precision = nposs # / float(i + 1)
        recall = i + 1 # nposs # / npositives
        rp.append((recall, precision))
    return sorted(rp)

def accumulate_quals(file, quals, get_label, selection):
    for line in file:
        line = line.strip()
        if line.startswith('#'): 
            if line.startswith('##fileformat=VCF'):
                fileformat = 'vcf'
            continue
        if line.startswith('Chrom'): 
            fileformat = 'position.xls'
            continue
        tokens = line.split('\t')
        #sys.stderr.write(line + '\n')
        start = tokens[1]
        end = tokens[2]
        ref = tokens[3]
        alt = tokens[4]
        if fileformat == 'vcf':
            qual = float(tokens[5]) if tokens[5] != '.' else (0 if len(ref) <= 9 else float('inf')) #0
        elif fileformat == 'position.xls':
            qual = float(tokens[6])
        if ref != alt and (
                (ref.lower() in 'a c g t n' and alt.lower() in 'a c g t n') 
                or (alt.lower().startswith('del') or alt == '-') 
                or len(ref) >= 4 and len(alt) == 1):
            if ('all' in selection
                    or 'false' in selection and 'PredictionStatus=FalsePositive' in line 
                    or 'true'  in selection and ('PredictionStatus=TruePositive' in line or 'VariantIsRefStdPositive' in line)
                    or 'hd734' in selection and 'PredictionStatus=HD734Positive' in line): 
                quals.append((qual, get_label, line))

def vcfs_to_quals(vcf_files, get_label = 0, selection = 'all'):
    quals = []
    if vcf_files: 
        for vcf_file in vcf_files:
            fileformat = 'none'
            with smartopen(vcf_file) as file:
                accumulate_quals(file, quals, get_label, selection)
    else:
        accumulate_quals(sys.stdin, quals, get_label, selection)
    return quals

def neg_pos_vcfs_to_recall_to_precision(negvcfs, posvcfs, negselection):
    quals0 = vcfs_to_quals(negvcfs, 0, negselection)
    quals1 = vcfs_to_quals(posvcfs, 1, 'hd734,true' + (',all' if 'pos-all' in negselection else ''))
    quals = sorted(quals0 + quals1, key = lambda x:-x[0])
    qs = [q[1] for q in quals]
    return zip(recall_to_precision_curve(qs), [q[0] for q in quals], [q[2] for q in quals])
    
def smartopen(fname):
    if fname.endswith('.bcf') or fname.endswith('.vcf.gz') or fname.endswith('.vcf'):
        cmd = ['/biocluster/data/bioexec/software/bcftools-1.3.1/bcftools', 'view', 
            '-i', ( sys.argv[4] if len(sys.argv) > 4 else ''), # '"DP>=1000"', 
            fname]
        sys.stderr.write('{}\n'.format(cmd))
        return subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
    else:
        return open(fname)

def infname_to_vcfs(infname):
    ret = []
    with smartopen(infname) as infile:
        for i, line in enumerate(infile):
            if i == 0 and (line.startswith('##fileformat=VCF') or line.startswith('Chrom\t')):
                return [infname] 
            ret.append(line.strip())
    return ret

def main():
    if len(sys.argv) > 3:
        negselection = sys.argv[3]
    else:
        negselection = 'all'
    if len(sys.argv) > 2: 
        posvcfs = infname_to_vcfs(sys.argv[2])
    else:
        posvcfs = []
    if len(sys.argv) > 1: 
        negvcfs = infname_to_vcfs(sys.argv[1])
    else:
        negvcfs = []
    #print(posvcfs)
    rps = neg_pos_vcfs_to_recall_to_precision(negvcfs, posvcfs, negselection)
    visited_strand_muts = set([])
    print('num_retrieved\tnum_relevant\tsort_by_value\tnum_dedup_retrived\toriginal_dataline')
    for (r, p), q, line in rps:
        tokens = line.split('\t')
        chrom = tokens[0]
        start = tokens[1]
        end = tokens[2]
        ref = tokens[3].upper()
        alt = tokens[4].upper()
        print('{}\t{}\t{}\t{}\t{}'.format(r, p, q, len(visited_strand_muts), line))
        visited_strand_muts.add((chrom, start, ref, alt))

if __name__ == '__main__':
    main()
 
