
from __future__ import print_function
import sys
import os
import pysam
from pysam import VCF
import numpy
import csv
import cPickle as pickle


vcf_file='mainset_November2015_chr22.vcf.gz'
thresh=[1,5,10,15,20,25,30,50,100]
print('#chrom','pos','mean','median', ' '.join([str(t) for t in thresh]))
vcf=pysam.VariantFile(vcf_file)
for v in vcf:
    chrom=v.chrom
    pos=v.pos
    dp=[int(v.samples[s]['DP'] or 0) for s in v.samples]
    n=len(dp)
    mu='%.2f' % numpy.mean(dp)
    med='%.2f' % numpy.median(dp)
    print(chrom, pos, mu, med, ' '.join(['%.4f'%(float(len(filter(lambda x: x>t, dp)))/float(n)) for t in thresh]))




