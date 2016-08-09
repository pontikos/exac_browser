from views import *
from lookups import *
import rest as annotation
import requests
import re
from utils import *
import itertools
import pysam
import csv
#hpo lookup
import phizz
import random
import orm

@app.route('/list_genes/')
def list_genes():
    # if gene not ensembl id then translate to
    db=get_db()
    genes=[g for g in db.genes.find(fields={'_id':False})]
    return(jsonify(result=genes))


@app.route('/gene/<gene_id>',methods=['GET'])
def gene_page(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    gene_name=db.genes.find_one({'gene_id':gene_id})['gene_name']
    print(gene_name)
    hpo_string=lookups.get_gene_hpo(get_db('hpo'),gene_name)
    #if gene_id in app.config['GENES_TO_CACHE']:
        #return open(os.path.join(app.config['GENE_CACHE_DIR'], '{}.html'.format(gene_id))).read()
    #else:
    return get_gene_page_content(gene_id,hpo,hpo_string)


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()
def get_gene_page_content(gene_id,hpo_id=None,hpo_string=""):
    db = get_db()
    patients_db=get_db('patients')
    hpo_db=get_db('hpo')
    # scrape exac
    #b=browser.Browser('exac.broadinstitute.org')
    #p=b.get_page('/gene/%s'%gene_id)
    #m = re.compile('window.table_variants\s*=\s*(.*)\s*;')
    #exac_table_variants=dict()
    #if m: exac_table_variants=loads(m.search(p).group(1))
    #gene = lookups.get_gene(db, gene_id)
    gene=orm.Gene(db=db, gene_id=gene_id)
    gene_name=gene.gene_name_upper
    if gene is None: abort(404)
    # gene and hpo term ideally
    cache_key = 't-gene-{}'.format(gene_id)
    t = cache.get(cache_key)
    print 'Rendering %sgene: %s' % ('' if t is None else 'cached ', gene_id)
    if t is not None:
        print 'Rendering gene: %s' % gene_id
        #print(phizz.query_gene(ensembl_id=str(gene_id)))
        return t
    # 1. get variants in gene from db
    #variants_in_gene = [v for v in db.variants.find({'variant_id':{'$in':gene['variant_ids']}},{'_id':0})]
    variants_in_gene=gene.get_variants_in_gene()
    print('variants_in_gene',len(variants_in_gene))
    # which transcript contains most of the variants
    transcript_counter=Counter([t for v in variants_in_gene for t in v.transcripts ])
    print(transcript_counter)
    transcript_with_most_variants=transcript_counter.most_common(1)[0][0]
    print('transcript with most variants',transcript_with_most_variants)
    transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)
    print('transcripts_in_gene',len(transcripts_in_gene))
    # Get some canonical transcript and corresponding info
    transcript_id = gene.canonical_transcript
    #transcript_id = transcript_with_most_variants
    # if none of the variants are on the canonical transcript use the transcript with the most variants on it
    transcript = orm.Transcript(db=db, transcript_id=transcript_id)
    variants_in_transcript = transcript.get_variants_in_transcript()
    #print('variants_in_transcript',len(variants_in_transcript))
    coverage_stats = lookups.get_coverage_for_transcript(db, transcript.xstart - EXON_PADDING, transcript.xstop + EXON_PADDING)
    print('coverage_stats',len(coverage_stats))
    add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
    constraint_info = lookups.get_constraint_for_transcript(db, transcript_id)
    #print('constraint_info',constraint_info)
    # 2. get patients with hpo term and patients without
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    #patient_ids=[p['external_id'] for p in patients]
    #hpo=phizz.query_hpo([hpo_id])[0]
    # samples
    everyone=frozenset(file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().split('\t'))
    #everyone=frozenset([p['external_id'] for p in patients_db.patients.find()])
    #for p in patients_db.patients.find(): if 'external_id' not in p: print(p)
    if hpo_id is None:
        hpo="HP:0000001"
        hpo_name='All'
        cases=frozenset()
        print('num cases',len(cases))
    else:
        #cases
        hpo_name=hpo_db.hpo.find_one({'id':hpo_id})['name'][0]
        print(hpo_name)
        cases=[p['external_id'] for p in lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)]
        cases=frozenset(cases)
        print('num cases',len(cases))
    #controls
    #everyone=frozenset(everyone.tolist())
    controls=everyone-cases
    #everyone=everyone & headers
    #controls=frozenset(everyone) - frozenset(cases)
    print('num controls',len(controls))
    # 3. for each variant tabix,  get counts and chisq
    #chrom, pos, ref, alt = variants_in_gene[0]['variant_id'].split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % gene.chrom,)
    region ='%s:%s-%s' % (str(gene.chrom), str(gene.start), str(gene.stop),)
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip('#').strip().split('\t')
    records=[dict(zip(headers,r.strip().split('\t'))) for r in tb.fetch(region)]
    print(len(records))
    records=dict([('%s-%s-%s-%s' % (r['CHROM'], r['POS'], r['REF'], r['ALT'],),r,) for r in records])
    #for i,_, in enumerate(variants_in_gene):
    for i in []:
        v=variants_in_gene[i]
        variant_str=v.variant_id
        print(variant_str)
        variant_str=str(variant_str).strip().replace('_','-')
        chrom, pos, ref, alt = variant_str.split('-')
        v.pos_coding_noutr= get_xpos(chrom, pos)
        #region=str('%s:%s-%s'%(chrom, pos, int(pos),))
        #records=tb.fetch(region=region)
        #geno=dict(zip(headers, [r.split('\t') for r in records][0]))
        if variant_str not in records:
            v.data["-log10pvalue"]=0
            continue
        geno=records[variant_str]
        #samples=[h for h in geno if geno[h].split(':')[0]=='0/1' or geno[h].split(':')[0]=='1/1']
        geno_cases=[geno[ca].split(':')[0] for ca in cases]
        caco=Counter(geno_cases)
        geno_controls=[geno[co].split(':')[0] for co in controls]
        coco=Counter(geno_controls)
        ca_mut=caco.get('1/1',0)+caco.get('0/1',0)
        ca_wt=caco.get('0/0',0)
        co_mut=coco.get('1/1',0)+coco.get('0/1',0)
        co_wt=coco.get('0/0',0)
        v.co_wt=co_wt
        v.ca_wt=ca_wt
        v.co_mut=co_mut
        v.ca_mut=ca_mut
        if ca_mut==0:
            v.data["-log10pvalue"]=0
        else:
            counts=numpy.array([[ca_mut,ca_wt],[co_mut,co_wt]])
            print(counts)
            stat=chisquare(counts)
            print(stat)
            v.data['-log10pvalue']=-math.log10(stat['p.value'])
            #print(chisquare([ca_mut,ca_wt,co_mut,co_wt]))
            #d=csv.DictReader(file('/data/uclex_files/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
            #headers=file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().replace('#','').split('\t')
            #get EXAC info
        for j,_, in enumerate(exac_table_variants):
            exac_v=exac_table_variants[j]
            if exac_v['variant_id']!=variant_str: continue
            #v.EXAC=exac_v
            v.HGVS=exac_v['HGVS']
            v.HGVSp=exac_v['HGVSp']
            v.HGVSc=exac_v['HGVSc']
            v.major_consequence=exac_v['major_consequence']
            v.exac_allele_freq=exac_v['allele_freq']
            break
        variants_in_gene[i]=v
    # gene hpo
    # dotfile='/slms/UGI/vm_exports/vyp/phenotips/HPO/dot/%s.dot' % gene_name
    # temporary dotfile for test
    dotfile='static/dot/literature_phenotype/%s.dot' % gene_name
    if os.path.isfile(dotfile):
        literature_DOT=file(dotfile,'r').read().replace('\n','\\n')
        # replace single quote
        literature_DOT=re.sub("'", '&#39;', literature_DOT)
        #fontsize=7
        # change fontsize to 7
        #DOT=re.sub(r'fontsize="\d+"', 'fontsize="%d"' % fontsize, DOT)
    else:
        literature_DOT=''
    simreg_DOT=''
    print('all good')
    # this is to set the pos
    for i,_, in enumerate(variants_in_gene):
        variant_id=variants_in_gene[i].variant_id
        for j,_, in enumerate(variants_in_transcript):
            variant_id2=variants_in_transcript[j]['variant_id']
            if variant_id==variant_id2:
                variants_in_transcript[j]['-log10pvalue']=variants_in_gene[i].data['-log10pvalue']
                break
    if not variants_in_transcript: variants_in_transcript=variants_in_gene
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?>",file('templates/variant_row.tmpl','r').read())
    #for v in variants_in_gene: if v.in_exac: print(v.ExAC_freq.keys())
    t=render_template( 'gene.html',
            gene=gene,
            table_headers=table_headers,
            transcript=transcript,
            variants_in_gene=variants_in_gene,
            variants_in_transcript=variants_in_transcript,
            transcripts_in_gene=transcripts_in_gene,
            constraint=constraint_info,
            csq_order=orm.csq_order,
            literature_DOT=hpo_string.replace('\n','\\n'),
            simreg_DOT=simreg_DOT,
            hpo_name=hpo_name,
            coverage_stats=coverage_stats)
    #cache.set(cache_key, t, timeout=1000*60)
    return t


