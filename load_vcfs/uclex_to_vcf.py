#! /sw/arch/bin/python
#
# *************************************************************
#
# $Source: $
# $Revision: $                                                                 
# $State: $                                                                     
# $Date: $                                                      
# $Author: $  
#
# $Log: $
#
#
# *************************************************************

from __future__ import print_function
import argparse
import os
import os.path
import sys
#import tabix
import pysam
from pysam import VCF
import csv
import numpy
import cPickle as pickle
# python data tables!
import pandas
from collections import defaultdict

#BASEDIR='~pontikos'
#BASEDIR='/Users/pontikos/'

usage_example = ' '
parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, epilog = usage_example) 
#compulsory arguments
parser.add_argument('--infile', dest='infile', help = "chromosome on which to do imputation", required=False, default=None)
#parser.add_argument('--phenotypes', dest='phenotypes', help = "file to which to write genotypes", default=None)
args = parser.parse_args()

#vcf_file='/goon2/scratch2/vyp-scratch2/vincent/GATK/mainset_November2015/mainset_November2015_chr22.vcf.gz'
vcf_file=args.infile

uclex_phenotypes='/goon2/scratch2/vyp-scratch2/UCLexInfo/uclex-samples.csv'
#uclex_phenotypes=args.phenotypes

#os.path.join(BASEDIR,'UCLexInfo','uclex-samples.csv')

#d=pandas.read_csv(file(os.path.join(BASEDIR,'UCLexInfo','uclex-samples.csv'),'r'),delimiter=',')
POPS=pandas.read_csv('~rmhanpo/uclex_browser/pop.csv')
POPS.index=POPS.phenotype

samplePheno=dict()
phenoSample=dict()
for l in csv.DictReader(file(uclex_phenotypes,'r'),delimiter=','):
    # "sample","owner","phenotype","ethnicity","sequencing.platform","capture","first.integrated"
    samplePheno[l['sample']]=l['phenotype']
    phenoSample[l['phenotype']]=phenoSample.get(l['phenotype'],[])+[l['sample']]


#tb=tabix.open(vcf)
#records=tb.querys('22:16084640-16957586')

#tb=pysam.TabixFile(vcf)
#records=tb.fetch(region='22:16084640-16957586')

#headers=[h for h in tb.header]
#headers=(headers[len(headers)-1]).strip().split('\t')

H=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

print('##fileformat=VCFv4.1')
print('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|HGNC_ID|BIOTYPE|CANONICAL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|EXON|INTRON|DOMAINS|HGVSc|HGVSp|GMAF|UCLEX2_MAF|AMR_MAF|ASN_MAF|EUR_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF">')
#print('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SIFT|PolyPhen|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|Condel|CAROL|CADD_PHRED|CADD_RAW|GO|ExAC_AF|ExAC_AF_AFR|ExAC_AF_AMR|ExAC_AF_EAS|ExAC_AF_FIN|ExAC_AF_NFE|ExAC_AF_OTH|ExAC_AF_SAS">')
print( '\t'.join(H) )

 
# Adjusted Alt Allele Counts (DP >= 10 & GQ >= 20)
# ##INFO=<ID=AC_Adj,Number=A,Type=Integer,Description="Adjusted Allele Counts">
# AC_Adj <= AC

# Number of Heterozygous Individuals (DP >= 10 & GQ >= 20)
# ##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Adjusted Heterozygous Counts">

# Number of Homozygous Alt Allele Individuals (DP >= 10 & GQ >= 20)
# ##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Adjusted Homozygous Counts">

# For chr1-22:
# sum(AC_Adj) = sum(AC_Het) + 2*sum(AC_Hom) 

vcf=pysam.VariantFile(vcf_file,'r')

#n=len(dp)
#mu='%.2f' % numpy.mean(dp)
#med='%.2f' % numpy.median(dp)
#print(chrom, pos, mu, med, ' '.join(['%.4f'%(float(len(filter(lambda x: x>t, dp)))/float(n)) for t in thresh]))

for v in vcf:
    # we only care about coding variants
    if 'CSQ' not in v.info.keys(): continue
    ac=defaultdict(list)
    genotypes=dict([(s,v.samples[s]['GT'],) for s in v.samples])
    geno=[v.samples[s]['GT'] for s in v.samples]
    for i,alt_allele, in enumerate(v.alts):
        het=geno.count((v.ref,alt_allele))
        hom=geno.count((alt_allele,alt_allele))
        ac['AC_Het'].append(het)
        #','.join(map(lambda x: str(x), [geno.count((v.ref,v.alts[i])) for i in range(0,len(v.alts))]))
        ac['AC_Hom'].append(hom)
        #','.join(map(lambda x: str(x), [geno.count((v.alts[i],v.alts[i])) for i in range(0,len(v.alts))]))
        AC=het+2*hom
        ac['AC'].append(AC)
        #','.join([int(x)+2*int(y) for x, y, in zip(ac['AC_Het'].split(','),ac['AC_Hom'].split(','))])
        AN=2*len(geno)
        ac['AN'].append(AN)
        ac['AF'].append(float(AC)/float(AN))
        if v.chrom == 'X':
            #AC_Hemi is a count Male alt alleles, where each site (excluding PAR) only has one allele.
            #AC_Hom is a count of Female alt alleles
            #AC_Het is a count of Female alt alleles
            #chrX (non-PAR)
            #sum(AC_Adj) = sum(AC_Hemi) + sum(AC_Het) + 2*sum(AC_Hom) 
            #AN_Adj (and all population AN on chrX) = 2*n_Female + n_Male
            #AN_Adj (and all population AN on chrY) = n_Male
            #ac['AC_Hemi'].append(
            #sum(AC_Adj) = sum(AC_Hemi) + sum(AC_Het) + 2*sum(AC_Hom) 
            ac['AC_Adj'].append(AC)
            #AN_Adj (and all population AN on chrX) = 2*n_Female + n_Male
            ac['AN_Adj'].append(AN)
        elif v.chrom == 'Y':
            #For chrY
            #sum(AC_Adj) = sum(AC_Hemi) 
            #sum(AC_Adj) = sum(AC_Hemi) 
            #AN_Adj (and all population AN on chrY) = n_Male
            ac['AC_Adj'].append(AC)
            ac['AN_Adj'].append(AN)
        else:
            ac['AC_Adj'].append(AC)
            ac['AN_Adj'].append(AN)
        for p in phenoSample:
            g=[genotypes[s] for s in phenoSample[p]]
            het=g.count((v.ref,alt_allele))
            hom=g.count((v.alts[i],alt_allele))
            ac['Hom_%s'%POPS.loc[p]['code']].append(hom)
            ac['Het_%s'%POPS.loc[p]['code']].append(het)
            ac['AC_%s'%POPS.loc[p]['code']].append(het+2*hom)
            ac['AN_%s'%POPS.loc[p]['code']].append(2*len(g))
            #print( p, len(g), g.count('./.'), g.count('0/0'), g.count('0/1'), g.count('1/1') )
        ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SIFT|PolyPhen|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|Condel|CAROL|CADD_PHRED|CADD_RAW|GO|ExAC_AF|ExAC_AF_AFR|ExAC_AF_AMR|ExAC_AF_EAS|ExAC_AF_FIN|ExAC_AF_NFE|ExAC_AF_OTH|ExAC_AF_SAS">
    print([','.join(map(lambda x: str(x), ac[k])) for k in ac])
    if len(v.alts)>1: break
#','.join(map(lambda x: str(x), ac
    CSQ=dict(zip("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SIFT|PolyPhen|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|Condel|CAROL|CADD_PHRED|CADD_RAW|GO|ExAC_AF|ExAC_AF_AFR|ExAC_AF_AMR|ExAC_AF_EAS|ExAC_AF_FIN|ExAC_AF_NFE|ExAC_AF_OTH|ExAC_AF_SAS".split('|'), v.info['CSQ'].split('|')))
    ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|HGNC_ID|BIOTYPE|CANONICAL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|EXON|INTRON|DOMAINS|HGVSc|HGVSp|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF">
    INFO = dict(v.info.items())
    CSQ='|'.join([CSQ.get(k,'') for k in ['Allele','Gene','Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','ALLELE_NUM','DISTANCE','STRAND','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','SIFT','PolyPhen','EXON','INTRON','DOMAINS','HGVSc','HGVSp','GMAF','AFR_MAF','AMR_MAF','ASN_MAF','EUR_MAF','AA_MAF','EA_MAF','CLIN_SIG','SOMATIC','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','LoF_info','LoF_flags','LoF_filter','LoF']])
    INFO['CSQ']=CSQ
    INFO=INFO.items()+[(k,ac[k],) for k in ac]
    INFO=';'.join(['%s=%s' % (a, b,) for a, b, in INFO] )
    #H=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    #print('\t'.join( [v[k] for k in H] ))
    FILTER=','.join([v.filter[k].name for k in v.filter]) or '.'
    print('\t'.join(map(lambda x: str(x), [v.chrom,v.pos,v.id or '.',v.ref,','.join(v.alts),v.qual,FILTER,INFO])))
    #print(v.filter.values())
    #print(dir(v.filter))


    
##INFO=<ID=DP_HIST,Number=A,Type=String,Description="Histogram for DP; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">
##INFO=<ID=GQ_HIST,Number=A,Type=String,Description="Histogram for GQ; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">

#d['FORMAT']
#INFO field
#d['InbreedingCoeff']
#d['MQ']
#d['MQ0']
#d['MQRankSum']
#d['NCC']
#d['QD']
#d['ReadPosRankSum']
#d['VQSLOD']



#1   13372   .   G   C   608.91  PASS    AC=3;AC_AFR=0;AC_AMR=0;AC_Adj=2;AC_EAS=0;AC_FIN=0;AC_Het=0;AC_Hom=1;AC_NFE=0;AC_OTH=0;AC_SAS=2;AF=6.998e-05;AN=42870;AN_AFR=770;AN_AMR=134;AN_Adj=8432;AN_EAS=254;AN_FIN=16;AN_NFE=2116;AN_OTH=90;AN_SAS=5052;BaseQRankSum=0.727;ClippingRankSum=1.15;DP=139843;FS=0.000;GQ_MEAN=12.48;GQ_STDDEV=15.18;Het_AFR=0;Het_AMR=0;Het_EAS=0;Het_FIN=0;Het_NFE=0;Het_OTH=0;Het_SAS=0;Hom_AFR=0;Hom_AMR=0;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=1;InbreedingCoeff=-0.0844;MQ=35.72;MQ0=0;MQRankSum=0.727;NCC=60853;QD=23.42;ReadPosRankSum=0.727;VQSLOD=-1.687e+00;culprit=MQ;DP_HIST=14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;GQ_HIST=1012|14971|172|100|3161|259|127|30|8|9|5|16|1162|274|59|45|17|2|3|3,0|0|0|0|0|1|0|0|0|0|0|0|0|1|0|0|0|0|0|0;

#CSQ=C|ENSG00000223972|ENST00000456328|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|620||||||1||1|DDX11L1|HGNC|37102|processed_transcript|YES||||||||3/3|||ENST00000456328.2:n.620G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000450305|Transcript|splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant|412||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene|||||||||5/6|||ENST00000450305.2:n.412G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000515242|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|613||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene|||||||||3/3|||ENST00000515242.2:n.613G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000518655|Transcript|intron_variant&non_coding_transcript_variant|||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene||||||||||2/3||ENST00000518655.2:n.482-31G>C|||||||||||||||||||,C||ENSR00000528767|RegulatoryFeature|regulatory_region_variant|||||||1||||||regulatory_region|||||||||||||||||||||||||||||||

#AC=3
#AC_AFR=0
#AC_AMR=0
#AC_Adj=2
#AC_EAS=0
#AC_FIN=0
#AC_Het=0
#AC_Hom=1
#AC_NFE=0
#AC_OTH=0
#AC_SAS=2
#AF=6.998e-05
#AN=42870
#AN_AFR=770
#AN_AMR=134
#AN_Adj=8432
#AN_EAS=254
#AN_FIN=16
#AN_NFE=2116
#AN_OTH=90
#AN_SAS=5052
#BaseQRankSum=0.727
#ClippingRankSum=1.15
#DP=139843
#FS=0.000
#GQ_MEAN=12.48
#GQ_STDDEV=15.18
#Het_AFR=0 #Het_AMR=0 #Het_EAS=0 #Het_FIN=0 #Het_NFE=0 #Het_OTH=0 #Het_SAS=0
#Hom_AFR=0 #Hom_AMR=0 #Hom_EAS=0 #Hom_FIN=0 #Hom_NFE=0 #Hom_OTH=0 #Hom_SAS=1
#InbreedingCoeff=-0.0844
#MQ=35.72
#MQ0=0
#MQRankSum=0.727
#NCC=60853
#QD=23.42
#ReadPosRankSum=0.727
#VQSLOD=-1.687e+00
#culprit=MQ
#DP_HIST=14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0
#GQ_HIST=1012|14971|172|100|3161|259|127|30|8|9|5|16|1162|274|59|45|17|2|3|3,0|0|0|0|0|1|0|0|0|0|0|0|0|1|0|0|0|0|0|0
#CSQ=C|ENSG00000223972|ENST00000456328|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|620||||||1||1|DDX11L1|HGNC|37102|processed_transcript|YES||||||||3/3|||ENST00000456328.2:n.620G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000450305|Transcript|splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant|412||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene|||||||||5/6|||ENST00000450305.2:n.412G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000515242|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|613||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene|||||||||3/3|||ENST00000515242.2:n.613G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000518655|Transcript|intron_variant&non_coding_transcript_variant|||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene||||||||||2/3||ENST00000518655.2:n.482-31G>C|||||||||||||||||||,C||ENSR00000528767|RegulatoryFeature|regulatory_region_variant|||||||1||||||regulatory_region|||||||||||||||||||||||||||||||





