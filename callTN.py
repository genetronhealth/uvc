import subprocess
import sys
import vcf

def get0(d, k, v=0):
    #print('Trying to get value with key {} from the dictionary {}'.format(k, d))
    if k in d and d[k]: 
        return d[k]
    else: 
        return v

LOD_THRES = 50

normalGT = ''
normalGQ = 0
normalGE = 0
normalBGrecord = None

p=subprocess.Popen(['bcftools', 'merge', '--force-samples', '-m', 'none', sys.argv[1], sys.argv[2]], stdout=subprocess.PIPE)
vcf_reader = vcf.Reader(p.stdout, 'r')
vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
for record in vcf_reader:
    assert record.FORMAT
    assert 2 == len(record.samples)
    tumor = record.samples[0].data.__dict__
    normal = record.samples[1].data.__dict__
    if (record.ALT and 1 == len(record.ALT) and str(record.ALT[0]) == '<NON_REF>' and
        'gEND' in normal and normal['gEND']):
        normalGT = normal['GT']
        normalGQ = normal['GQ']
        normalGE = normal['gEND']
        normalBGrecord = record
        #print('ALT is {}'.format(record.ALT[0]))
    
    tumorFA =  get0(tumor, 'FA')  + 1e-7
    normalFA = get0(normal, 'FA') + 1e-7
    tlod = (get0(tumor, 'VAQ') - get0(normal, 'VAQ')) * max((tumorFA - normalFA, 0)) / (tumorFA + normalFA)
    if tlod < LOD_THRES: continue
    #print('The following record passed TLOD threshold: {}'.format(record))
    nlod = (0 if (record.POS > normalGE or normalGT in ['0|1', '1|1', '0/1', '1/1']) else normalGQ) + 30
    if nlod < LOD_THRES: continue
    #print('The following record passed NLOD threshold: {}'.format(record)) 
    record.QUAL = min((tlod, nlod))
    if record.REF != record.ALT[0]:
        vcf_writer.write_record(normalBGrecord)
        vcf_writer.write_record(record)
        
