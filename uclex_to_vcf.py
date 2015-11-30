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
import os.path
import sys
#import tabix
import pysam
import csv

# python data tables!
import pandas

#BASEDIR='~pontikos'
#BASEDIR='/Users/pontikos/'
BASEDIR='../uclex_data/'

#d=pandas.read_csv(file(os.path.join(BASEDIR,'UCLexInfo','uclex-samples.csv'),'r'),delimiter=',')
POPS=pandas.read_csv('pop.csv')
POPS.index=POPS.phenotype


samplePheno=dict()
phenoSample=dict()
for l in csv.DictReader(file(os.path.join(BASEDIR,'UCLexInfo','uclex-samples.csv'),'r'),delimiter=','):
    # "sample","owner","phenotype","ethnicity","sequencing.platform","capture","first.integrated"
    samplePheno[l['sample']]=l['phenotype']
    phenoSample[l['phenotype']]=phenoSample.get(l['phenotype'],[])+[l['sample']]

vcf=os.path.join(BASEDIR,'uclex','mainset_November2015_chr22.vcf.gz') 

#tb=tabix.open(vcf)
#records=tb.querys('22:16084640-16957586')

tb=pysam.TabixFile(vcf)
records=tb.fetch(region='22:16084640-16957586')

headers=[h for h in tb.header]
headers=(headers[len(headers)-1]).strip().split('\t')

H=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

print('##fileformat=VCFv4.1')
print('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|HGNC_ID|BIOTYPE|CANONICAL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|EXON|INTRON|DOMAINS|HGVSc|HGVSp|GMAF|UCLEX2_MAF|AMR_MAF|ASN_MAF|EUR_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF">')
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


for r in records:
    r=r.strip().split('\t')
    d=dict(zip(headers,r))
    if ',' in d['ALT']: continue
    genotypes=dict([ (s, d[s].split(':')[0]) for s in samplePheno ])
    ac=dict()
    geno=[genotypes[s] for s in samplePheno]
    ac['AC_Het']=geno.count('0/1')
    ac['AC_Hom']=geno.count('1/1')
    ac['AC']=ac['AC_Het']+2*ac['AC_Hom']
    #ac['AC_Adj']
    ac['AN']=2*len(geno)
    ac['AF']=float(ac['AC'])/float(ac['AN'])
    for p in phenoSample:
        g=[genotypes[s] for s in phenoSample[p]]
        ac['Hom_%s'%POPS.loc[p]['code']]=g.count('1/1')
        ac['Het_%s'%POPS.loc[p]['code']]=g.count('0/1')
        ac['AC_%s'%POPS.loc[p]['code']]=g.count('0/1')+2*g.count('1/1')
        ac['AN_%s'%POPS.loc[p]['code']]=2*len(g)
        #d.where(d['phenotype']=='unknown')
        #print( p, len(g), g.count('./.'), g.count('0/0'), g.count('0/1'), g.count('1/1') )
    #break
    ac['AC_Adj']=ac['AC']
    ac['AN_Adj']=ac['AN']
    d['INFO']=d['INFO']+';'+';'.join( ['%s=%s'%(k,ac[k],) for k in ac] )
    print('\t'.join( [d[k] for k in H] ))

    





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



#1   13372   .   G   C   608.91  PASS    AC=3;AC_AFR=0;AC_AMR=0;AC_Adj=2;AC_EAS=0;AC_FIN=0;AC_Het=0;AC_Hom=1;AC_NFE=0;AC_OTH=0;AC_SAS=2;AF=6.998e-05;AN=42870;AN_AFR=770;AN_AMR=134;AN_Adj=8432;AN_EAS=254;AN_FIN=16;AN_NFE=2116;AN_OTH=90;AN_SAS=5052;BaseQRankSum=0.727;ClippingRankSum=1.15;DP=139843;FS=0.000;GQ_MEAN=12.48;GQ_STDDEV=15.18;Het_AFR=0;Het_AMR=0;Het_EAS=0;Het_FIN=0;Het_NFE=0;Het_OTH=0;Het_SAS=0;Hom_AFR=0;Hom_AMR=0;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=1;InbreedingCoeff=-0.0844;MQ=35.72;MQ0=0;MQRankSum=0.727;NCC=60853;QD=23.42;ReadPosRankSum=0.727;VQSLOD=-1.687e+00;culprit=MQ;DP_HIST=14728|2455|2120|518|121|499|534|314|111|21|10|2|2|0|0|0|0|0|0|0,1|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;GQ_HIST=1012|14971|172|100|3161|259|127|30|8|9|5|16|1162|274|59|45|17|2|3|3,0|0|0|0|0|1|0|0|0|0|0|0|0|1|0|0|0|0|0|0;CSQ=C|ENSG00000223972|ENST00000456328|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|620||||||1||1|DDX11L1|HGNC|37102|processed_transcript|YES||||||||3/3|||ENST00000456328.2:n.620G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000450305|Transcript|splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant|412||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene|||||||||5/6|||ENST00000450305.2:n.412G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000515242|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|613||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene|||||||||3/3|||ENST00000515242.2:n.613G>C|||||||||||||||||||,C|ENSG00000223972|ENST00000518655|Transcript|intron_variant&non_coding_transcript_variant|||||||1||1|DDX11L1|HGNC|37102|transcribed_unprocessed_pseudogene||||||||||2/3||ENST00000518655.2:n.482-31G>C|||||||||||||||||||,C||ENSR00000528767|RegulatoryFeature|regulatory_region_variant|||||||1||||||regulatory_region|||||||||||||||||||||||||||||||

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





