from views import *
from lookups import *
import rest as annotation
import requests
import primer3
import myvariant
import re
from utils import *
import itertools
import pysam
import csv
#hpo lookup
import phizz
import random
import orm


@app.route('/variant3/<variant_str>')
def variant_page3(variant_str):
    db=get_db()
    variant=db.variants.find_one({'VARIANT_ID':variant_str})
    patients=[p for p in db.patients.find({'external_id':{'$in': variant['HET']+variant['HOM']}})]
    hpo_terms=[p['features'] for p in patients]
    print(hpo_terms)
    print 'Rendering variant: %s' % variant_str
    return render_template( 'test.html', variant=variant)

@app.route('/variant/<variant_str>')
def variant_page(variant_str):
    db = get_db()
    variant=orm.Variant(db=db,variant_id=variant_str)
    # pos, ref, alt = get_minimal_representation(pos, ref, alt)
    #v=load_variant(db,variant_id)
    #xpos = get_xpos(chrom, pos)
    if variant is None:
        variant = {
            'chrom': chrom,
            'pos': pos,
            'xpos': xpos,
            'ref': ref,
            'alt': alt
        }
    consequences = []
    ordered_csqs = []
    # Adds major_consequence
    #base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
    base_coverage = []
    #any_covered = any([x['has_coverage'] for x in base_coverage])
    any_covered = any([x['has_coverage'] for x in base_coverage])
    # check the appropriate sqlite db to get the *expected* number of
    # available bams and *actual* number of available bams for this variant
    print 'Rendering variant: %s' % variant_str
    return render_template(
        'variant.html',
        variant=variant,
        base_coverage=base_coverage,
        consequences=consequences,
        any_covered=any_covered,
        ordered_csqs=ordered_csqs,
        metrics=[]
    )


@app.route('/variant_json/<variant_str>')
def variant_json(variant_str):
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    #mainset_February2016_chrX_filtered.vcf.gz
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    records=tb.fetch(region=region)
    records=[r.split('\t') for r in records]
    for r in records:
        geno=dict(zip(headers, r))
        POS=geno['POS']
        REF=geno['REF']
        print 'POS', POS
        print 'REF', REF
        for i, ALT, in enumerate(geno['ALT'].split(',')):
            print 'ALT', ALT
            # insertion
            if ref=='-' and REF+alt==ALT:
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # deletion
            # replace leftmost
            elif alt=='-' and ALT==REF.replace(ref,''):
                return reponse(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # replace rightmost
            elif alt=='-' and ALT==REF[::-1].replace(ref[::-1], "", 1)[::-1]:
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            # 
            elif alt=='-' and ref==REF and ALT=='*':
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            elif alt==ALT and ref==REF:
                return response(POS=int(POS), REF=REF, ALT=ALT, index=i+1, geno=geno, chrom=chrom, pos=pos)
            continue

@app.route('/set_variant_causal/<individual>/<variant_str>')
def set_variant_causal(individual, variant_str):
    print individual, variant_str
    db=get_db()
    #get_db().patients.update({'patient_id':individual},{'$addToSet':{'causal_variants':variant_str}})
    var=db.variants.find_one({'variant_id':variant_str})
    gene_id=var['genes'][0]
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name_upper']
    print 'GENE_NAME', gene_name
    # update Gene in phenotips
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=individual,auth=auth)
    p['genes']=p.get('genes',[])+[{'gene':gene_name}]
    print conn.update_patient( eid=p['external_id'], auth=auth, patient=p )
    print get_db('patients').patients.update({'external_id':individual},{'$set':p},w=0)
    p=db.patients.find_one({'external_id':individual})
    p['causal_variants']=list(frozenset(p.get('causal_variants',[])+[variant_str]))
    db.patients.update({'external_id':individual},{'$set':{'causal_variants':p['causal_variants']}},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)

@app.route('/unset_variant_causal/<individual>/<variant_str>')
def unset_variant_causal(individual, variant_str):
    print individual, variant_str
    db=get_db()
    p=db.patients.find_one({'external_id':individual})
    if 'causal_variants' in p and not p['causal_variants']: p['causal_variants']=[]
    if variant_str in p.get('causal_variants',[]):
        p['causal_variants']=p['causal_variants'].remove(variant_str)
    db.patients.update({'external_id':individual},{'$set':{'causal_variants':p['causal_variants']}},w=0)
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p2=conn.get_patient(eid=individual,auth=auth)
    p2['genes']=[]
    for var in p['causal_variants']:
        var=db.variants.find_one({'variant_id':var})
        gene_id=var['genes'][0]
        gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name_upper']
        print 'GENE_NAME', gene_name
        p2['genes']=list(frozenset(p2.get('genes',[])+[{'gene':gene_name}]))
    # update Gene in phenotips
    print conn.update_patient( eid=p2['external_id'], auth=auth, patient=p2 )
    print get_db('patients').patients.update({'external_id':individual},{'$set':p2},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)

@app.route('/set_variant_status/<individual>/<variant_str>/<status>')
def set_variant_status(individual, variant_str, status):
    print individual, variant_str, status
    db=get_db()
    #print get_db().patients.update({'patient_id':individual},{'$addToSet':{'variant_status':{variant_str:status}}})
    rare_variants=db.patients.find_one({'external_id':individual},{'rare_variants':1})['rare_variants']
    for rv in rare_variants:
        if rv['variant_id']==variant_str:
            rv['status']=status
    print db.patients.update({'external_id':individual},{'$set':{'rare_variants':rare_variants}})
    return status


