#! /bin/env python
from __future__ import print_function
import sys
import argparse
import os
import os.path
#import tabix
import pysam
from pysam import VCF
import csv
import cPickle as pickle
# python data tables!
#import pandas
from collections import defaultdict
import pymongo
#from bson import Binary, Code
#from bson.json_util import dumps
from StringIO import StringIO
import json

usage_example = ' '
parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, epilog = usage_example) 
parser.add_argument('--infile', dest='infile', help = "chromosome on which to do imputation", required=True, default=None)
parser.add_argument('--login', dest='login', help = "login for retrieving phenotypes", required=False, default=None)
args = parser.parse_args()

vcf_file=args.infile
vcf=pysam.VariantFile(vcf_file,'r')

# variant database
#client = pymongo.MongoClient(host='phenotips',port=27017)
#variants_db=client['uclex']

# phenotips API to retrieve patient HPO terms
if args.login:
    import phizz
    from phenotips_client import PhenotipsClient
    phenotips=PhenotipsClient(host='phenotips',port='8080',debug=False)

# get the CSQ header
desc=vcf.header.info['CSQ'].record['Description']
CSQ=desc.replace('Consequence annotations from Ensembl VEP. Format:','').replace(' ','').replace('"','').strip().split('|') 

protein_letters = {
            'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
            'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
            'O': 'Pyl', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
            'S': 'Ser', 'T': 'Thr', 'U':'Sec', 'V': 'Val',
            'W': 'Trp', 'Y': 'Tyr', 'X': 'Ter', '*': 'Ter',
            '-':''
}

def phenotypes(variant):
    for s in variant['HOM']+variant['HET']:
        for hpo in phenotips.patient_hpo(eid=s,auth=args.login):
            if hpo.startswith('HP'):
                variant['HPO']=variant.get('HPO',[])+phizz.query_hpo([hpo])
            elif hpo.startswith('MIM'):
                variant['HPO']=variant.get('HPO',[])+phizz.query_disease([hpo])
    #for hpo in variant['HPO']: variant['HPO'+=get_ancestors(hpo['hpo_term'])
    if 'HPO' in variant:
        #variant['HPO']=list(set([h['hpo_term'] for h in variant['HPO']]))
        variant['HPO']=[h['hpo_term'] for h in variant['HPO']]
    else:
        variant['HPO']=[]
    return(variant)


for v in vcf:
    for alt in v.alts:
        variant=dict()
        variant['VARIANT_ID']='%s-%s-%s-%s' % (v.chrom, v.pos,v.ref,alt,)
        variant['pos']=v.pos
        variant['start']=v.start
        variant['stop']=v.stop
        variant['ref']=v.ref
        variant['alt']=alt
        variant['alleles']=v.alleles
        variant['alts']=v.alts
        variant['rlen']=v.rlen
        variant['chrom']=v.chrom
        variant['id']=v.id
        variant['rid']=v.rid
        variant['qual']=v.qual
        #print(v.contig)
        variant['filter']=[k for k in v.filter.keys()]
        variant['format']=dict([(v.format[k].name,v.format[k].id,) for k in v.format.keys()])
        variant['info']=dict(v.info)
        if 'CSQ' in variant['info']:
            # one csq per transcript
            # use transcript id as key
            variant['CSQ']=dict()
            for transcript_csq in variant['info']['CSQ'].split(','):
                transcript_csq=dict(zip(CSQ,transcript_csq.split('|')))
                pp=transcript_csq['Protein_position']
                amino_acids=transcript_csq['Amino_acids']
                if len(amino_acids)==1:
                    aa=protein_letters.get(amino_acids,'')
                    transcript_csq['Amino_acids']="p."+aa+pp+aa
                elif '/' in amino_acids:
                    aa=[ protein_letters.get(x,'') for x in amino_acids.split('/')[0] ]
                    aa2=[ protein_letters.get(x,'') for x in amino_acids.split('/')[1] ]
                    transcript_csq['Amino_acids']='p.%s'%''.join(map(lambda x: ''.join(x),zip(aa,pp.split('-'))))
                else:
                    aa=[ protein_letters.get(x,'') for x in amino_acids ]
                    transcript_csq['Amino_acids']='p.%s'%''.join(map(lambda x: ''.join(x),zip(aa,pp.split('-'))))
                variant['CSQ'][transcript_csq['Feature']]=transcript_csq
            del variant['info']['CSQ']
            variant['SYMBOL']=list(set([variant['CSQ'][trans]['SYMBOL'] for trans in variant['CSQ']]))
            variant['Gene']=list(set([variant['CSQ'][trans]['Gene'] for trans in variant['CSQ']]))
            variant['Transcript']=[variant['CSQ'][trans]['Feature'] for trans in variant['CSQ']]
            #print(variant)
            #print(variant['CSQ']['Amino_acids'].split('/'))
            #print(variant['CSQ']['Allele'])
            # only aa sub of the canonical (sometimes no canonical though)
            variant['Canonical_Transcript']=[t for t in variant['CSQ'] if 'CANONICAL' in variant['CSQ'][t] and variant['CSQ'][t]['CANONICAL']=='YES']
        else:
            variant['CSQ']=[]
        #print(v.samples.keys())
        # this needs to go in a separate table
        variant['MISS']=len([s for s in v.samples if v.samples[s]['GT'].count(None)==2])
        variant['WT']=len([s for s in v.samples if v.samples[s]['GT'].count(variant['ref'])==2])
        # hets / homs collection
        for s in v.samples :
            x={ 'variant_id':variant['VARIANT_ID'], 'sample':s }
            io=StringIO()
            json.dump(x,io)
            if v.samples[s]['GT'].count(variant['alt'])==1:
                print('hets', io.getvalue())
            elif v.samples[s]['GT'].count(variant['alt'])==2:
                print('homs', io.getvalue())
        variant['HET']=[s for s in v.samples if v.samples[s]['GT'].count(variant['alt'])==1]
        variant['HOM']=[s for s in v.samples if v.samples[s]['GT'].count(variant['alt'])==2]
        # lookup the HPO for each sample in which the variant is seen
        if args.login: variant=phenotypes(variant)
        variant['HET']=len(variant['HET'])
        variant['HOM']=len(variant['HOM'])
        #print(v.info.keys())
        #print(variant)
        #variants_db.variants.update_one({'VARIANT_ID':variant['VARIANT_ID']},variant,upsert=True)
        #bson.BSON.encode(variant)
        # flatten info
        if 'info' in variant:
            for k in variant.get('info',[]): variant[k]=variant['info'][k]
            del variant['info']
        if 'SYMBOL' in variant:
            x={ 'variant_id':variant['VARIANT_ID'], 'gene':variant['SYMBOL'][0], 'gene_id':variant['Gene'] }
            io=StringIO()
            json.dump(x,io)
            print('gene', x)
        io=StringIO()
        json.dump(variant,io)
        # variants collection
        print('variants', io.getvalue())




