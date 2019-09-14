#!/usr/bin/env python
import sys
import gzip

def iter_fq(fqfile):
    hdr, seq, comment, qual = ('', '', '', '')
    for i, line in enumerate(fqfile):
        if i % 4 == 0:
            assert line.startswith('@')
            hdr = line[1:] # SRR3478158.1565
            #hdrtokens = hdr.split()[0].split('.')
            #if len(hdrtokens) > 2 and hdrtokens[2] in ['1', '2']:
            #    hdrtokens[2] = '12'
            #    hdr = '.'.join(hdrtokens)
        elif i % 4 == 1:
            seq = line
        elif i % 4 == 2:
            assert line.startswith('+') 
            comment = line[1:]
        else:
            qual = line
            yield hdr, seq, comment, qual

def change_hdr(hdr, seq1, umi_beg1, umi_end1, seq2, umi_beg2, umi_end2):
    hdr2 = hdr.rstrip().split()
    if (len(hdr2) > 1 and len(hdr2[1].split(':')[0]) >= 3):
        hdr3 = hdr2[1]
    else:
        hdr3 = hdr2[0]
        hdrtokens = hdr3.split('.')
        if len(hdrtokens) > 2 and hdrtokens[2] in ['1', '2']:
           hdrtokens[2] = '12'
           hdr3 = '.'.join(hdrtokens)
    r1tag = (seq1[umi_beg1:umi_end1] if len(seq1) > umi_end1 else ''.join(['N' for _ in range(umi_end1-umi_beg1)]))
    r2tag = (seq2[umi_beg2:umi_end2] if len(seq2) > umi_end2 else ''.join(['N' for _ in range(umi_end2-umi_beg2)]))
    hdr4 = hdr3 + '#' + ((r1tag + '+' + r2tag) if (r1tag != '' and r2tag != '') else (r1tag + r2tag))
    return hdr4 + '\n'

def proc(r1infile, r2infile, r1outfile, r2outfile, incluBeg1 = 0, excluEnd1 = 11, incluBeg2 = 0, excluEnd2 = 11):
    g1, g2 = iter_fq(r1infile), iter_fq(r2infile)
    rec1 = next(g1, -1)
    rec2 = next(g2, -1)
    while rec1 != -1 and rec2 != -1:
        #hdr1 = change_hdr(rec1[0], rec1[1])
        #hdr2 = change_hdr(rec2[0], rec1[1])
        hdr1 = change_hdr(rec1[0], rec1[1], incluBeg1, excluEnd1, rec2[1], incluBeg2, excluEnd2)
        hdr2 = change_hdr(rec2[0], rec1[1], incluBeg1, excluEnd1, rec2[1], incluBeg2, excluEnd2)
        r1outfile.write('@{}{}+{}{}'.format(hdr1, rec1[1], rec1[2], rec1[3]))
        r2outfile.write('@{}{}+{}{}'.format(hdr2, rec2[1], rec2[2], rec2[3])) 
        rec1 = next(g1, -1)
        rec2 = next(g2, -1)
    assert rec1 == -1 and rec2 == -1    
 
r1in = sys.argv[1]
r2in = sys.argv[2]

r1out = sys.argv[3]
r2out = sys.argv[4]

assert 4 == len(set([r1in, r2in, r1out, r2out]))

if len(sys.argv) > 6:
    incluBeg = int(sys.argv[5])
    excluEnd = int(sys.argv[6])
else:
    incluBeg = 0
    excluEnd = 23

if len(sys.argv) > 8:
    incluBeg2 = int(sys.argv[5])
    excluEnd2 = int(sys.argv[6])
else:
    incluBeg2 = 0
    excluEnd2 = 0

with gzip.open(r1in) as r1infile, gzip.open(r2in) as r2infile, gzip.open(r1out, 'wb', compresslevel=1) as r1outfile, gzip.open(r2out, 'wb', compresslevel=1) as r2outfile:
    proc(r1infile, r2infile, r1outfile, r2outfile, incluBeg, excluEnd, incluBeg2, excluEnd2)


