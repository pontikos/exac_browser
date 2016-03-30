#!/usr/bin/env python2

from scipy.stats import chisquare
import math

from Bio import Entrez
from phenotips_python_client import PhenotipsClient
from phenotips_python_client import browser
from bson.json_util import loads
from mongodb import *
# fizz: hpo lookup
import phizz
import itertools
import json
import os
import pymongo
import pysam
import gzip
from parsing import *
import logging
import lookups
import random
import sys
from utils import *

from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify, send_from_directory
from flask.ext.compress import Compress
from flask.ext.runner import Runner
from flask_errormail import mail_on_500

from flask import Response
from collections import defaultdict, Counter
from werkzeug.contrib.cache import SimpleCache

from multiprocessing import Process
import glob
import sqlite3
import traceback
import time

from functools import wraps
from flask import request, Response

from flask_debugtoolbar import DebugToolbarExtension
from werkzeug.exceptions import default_exceptions, HTTPException

import pandas
import csv
import time
import StringIO

from urlparse import urlparse
import pickle

#import pdb

from flask import Flask, session
from flask.ext.session import Session

# handles live plotting if necessary
import math
import plotly
print plotly.__version__  # version >1.9.4 required
from plotly.graph_objs import Scatter, Layout

# connect to R session
import pyRserve

import numpy


logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().setLevel(logging.INFO)

ADMINISTRATORS = ( 'n.pontikos@ucl.ac.uk',)
app = Flask(__name__)
mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = True
cache = SimpleCache(default_timeout=60*60*24)


REGION_LIMIT = 1E5
EXON_PADDING = 50
# Load default config and override config from an environment variable
#app.config.from_pyfile('uclex.cfg')
app.config.from_pyfile('uclex-old.cfg')

# Check Configuration section for more details
#SESSION_TYPE = 'null'
SESSION_TYPE = 'mongodb'
#SESSION_USE_SIGNER=True
app.config.from_object(__name__)
sess=Session()
sess.init_app(app)


@app.route('/set/<query>')
def set(query):
    value = query
    session['key'] = value
    return value

@app.route('/get/')
def get():
    return session.get('key', 'not set')

def check_auth(username, password):
    """
    This function is called to check if a username / password combination is valid.
    Will try to connect to phenotips instance.
    """
    conn=PhenotipsClient()
    response=conn.get_patient(auth='%s:%s' % (username, password,),number=1)
    if response:
        # setting a session key for pubmedBatch to save result
        session['user'] = username
        session['password'] = password
        return True
    else:
        return False

def authenticate():
    """Sends a 401 response that enables basic auth"""
    return Response( 'Could not verify your access level for that URL.\n' 'You have to login with proper credentials', 401, {'WWW-Authenticate': 'Basic realm="Login Required"'})

def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password): return authenticate()
        return f(*args, **kwargs)
    return decorated

def get_db(dbname=None):
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if dbname is None: dbname=app.config['DB_NAME']
    if not hasattr(g, 'db_conn'):
        g.db_conn=dict()
        g.db_conn[dbname] = connect_db(dbname)
    elif dbname not in g.db_conn:
        g.db_conn[dbname] = connect_db(dbname)
    return g.db_conn[dbname]

def get_R_session():
    if not hasattr(g, 'R_session'): g.R_session=pyRserve.connect()
    return g.R_session

def get_hpo_graph():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if not hasattr(g, 'hpo_graph'):
        from hpo_similarity.ontology import Ontology
        from hpo_similarity.similarity import CalculateSimilarity
        ontology=Ontology(app.config['HPO_OBO'])
        g.hpo_graph=ontology.get_graph()
    return g.hpo_graph


def connect_db(dbname=None):
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    print(client)
    if not dbname: dbname=app.config['DB_NAME']
    print(dbname)
    return client[dbname]


def parse_tabix_file_subset(tabix_filenames, subset_i, subset_n, record_parser):
    """
    Returns a generator of parsed record objects (as returned by record_parser) for the i'th out n subset of records
    across all the given tabix_file(s). The records are split by files and contigs within files, with 1/n of all contigs
    from all files being assigned to this the i'th subset.

    Args:
        tabix_filenames: a list of one or more tabix-indexed files. These will be opened using pysam.Tabixfile
        subset_i: zero-based number
        subset_n: total number of subsets
        record_parser: a function that takes a file-like object and returns a generator of parsed records
    """
    start_time = time.time()
    print(tabix_filenames)
    open_tabix_files = [pysam.Tabixfile(tabix_filename) for tabix_filename in tabix_filenames]
    tabix_file_contig_pairs = [(tabix_file, contig) for tabix_file in open_tabix_files for contig in tabix_file.contigs]
    # get every n'th tabix_file/contig pair
    tabix_file_contig_subset = tabix_file_contig_pairs[subset_i : : subset_n]
    short_filenames = ", ".join(map(os.path.basename, tabix_filenames))
    print(short_filenames)
    num_file_contig_pairs = len(tabix_file_contig_subset)
    print(("Loading subset %(subset_i)s of %(subset_n)s total: %(num_file_contig_pairs)s contigs from %(short_filenames)s") % locals())
    counter = 0
    for tabix_file, contig in tabix_file_contig_subset:
        header_iterator = tabix_file.header
        records_iterator = tabix_file.fetch(contig, 0, 10**9, multiple_iterators=True)
        for parsed_record in record_parser(itertools.chain(header_iterator, records_iterator)):
            counter += 1
            yield parsed_record
            if counter % 100000 == 0:
                seconds_elapsed = int(time.time()-start_time)
                print(("Loaded %(counter)s records from subset %(subset_i)s of %(subset_n)s from %(short_filenames)s " "(%(seconds_elapsed)s seconds)") % locals())
    print("Finished loading subset %(subset_i)s from  %(short_filenames)s (%(counter)s records)" % locals())



def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    autocomplete_strings = []
    for gene in get_db().genes.find():
        autocomplete_strings.append(gene['gene_name'])
        if 'other_names' in gene:
            autocomplete_strings.extend(gene['other_names'])
    f = open(os.path.join(app.config['UCLEX_FILES_DIRECTORY'],'autocomplete_strings.txt'), 'w')
    for s in sorted(autocomplete_strings):
        f.write(s+'\n')
    f.close()
    # create static gene pages for genes in
    if not os.path.exists(app.config['GENE_CACHE_DIR']): os.makedirs(app.config['GENE_CACHE_DIR'])
    # get list of genes ordered by num_variants
    for gene_id in app.config['GENES_TO_CACHE']:
        try:
            page_content = get_gene_page_content(gene_id)
        except Exception as e:
            print e
            continue
        f = open(os.path.join(app.config['GENE_CACHE_DIR'], '{}.html'.format(gene_id)), 'w')
        f.write(page_content)
        f.close()


def precalculate_metrics():
    db = get_db()
    print 'Reading %s variants...' % db.variants.count()
    metrics = defaultdict(list)
    binned_metrics = defaultdict(list)
    progress = 0
    start_time = time.time()
    for variant in db.variants.find(fields=['quality_metrics', 'site_quality', 'allele_num', 'allele_count']):
        for metric, value in variant['quality_metrics'].iteritems():
            metrics[metric].append(float(value))
        qual = float(variant['site_quality'])
        metrics['site_quality'].append(qual)
        if variant['allele_num'] == 0: continue
        if variant['allele_count'] == 1:
            binned_metrics['singleton'].append(qual)
        elif variant['allele_count'] == 2:
            binned_metrics['doubleton'].append(qual)
        else:
            for af in AF_BUCKETS:
                if float(variant['allele_count'])/variant['allele_num'] < af:
                    binned_metrics[af].append(qual)
                    break
        progress += 1
        if not progress % 100000:
            print 'Read %s variants. Took %s seconds' % (progress, int(time.time() - start_time))
    print 'Done reading variants. Dropping metrics database... '
    db.metrics.drop()
    print 'Dropped metrics database. Calculating metrics...'
    for metric in metrics:
        bin_range = None
        data = map(numpy.log, metrics[metric]) if metric == 'DP' else metrics[metric]
        if metric == 'FS':
            bin_range = (0, 20)
        elif metric == 'VQSLOD':
            bin_range = (-20, 20)
        elif metric == 'InbreedingCoeff':
            bin_range = (0, 1)
        if bin_range is not None:
            data = [x if (x > bin_range[0]) else bin_range[0] for x in data]
            data = [x if (x < bin_range[1]) else bin_range[1] for x in data]
        hist = numpy.histogram(data, bins=40, range=bin_range)
        edges = hist[1]
        # mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        lefts = [edges[i] for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': metric,
            'mids': lefts,
            'hist': list(hist[0])
        })
    for metric in binned_metrics:
        hist = numpy.histogram(map(numpy.log, binned_metrics[metric]), bins=40)
        edges = hist[1]
        mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': 'binned_%s' % metric,
            'mids': mids,
            'hist': list(hist[0])
        })
    db.metrics.ensure_index('metric')
    print 'Done pre-calculating metrics!'

# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()


@app.route('/')
@requires_auth
def homepage():
    cache_key = 't-homepage'
    #t = cache.get(cache_key)
    #if t: return t
    db=get_db()
    total_variants=db.variants.count()
    patients_db=get_db('patients')
    total_patients=patients_db.patients.count()
    male_patients=patients_db.patients.find( {'sex':'M'}).count()
    female_patients=patients_db.patients.find( {'sex':'F'}).count()
    unknown_patients=patients_db.patients.find( {'sex':'U'}).count()
    rsession=get_R_session()
    print(rsession.eval('length(variants)'))
    dotfile='static/dot/patients.dot'
    DOT=file(dotfile,'r').read().replace('\n','\\n')
    # replace single quote
    DOT=re.sub("'", '&#39;', DOT)
    #fontsize=7
    # change fontsize to 7
    #DOT=re.sub(r'fontsize="\d+"', 'fontsize="%d"' % fontsize, DOT)
    t = render_template('homepage.html',total_patients=total_patients,total_variants=total_variants,male_patients=male_patients,female_patients=female_patients,unknown_patients=unknown_patients,DOT=DOT)
    #cache.set(cache_key, t)
    return t


@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    if not hasattr(g, 'autocomplete_strings'): g.autocomplete_strings = [s.strip() for s in open(os.path.join(app.config['UCLEX_FILES_DIRECTORY'], 'autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = str(request.args.get('query'))
    #for n in dir(request): print(n, getattr(request,n))
    #print(request.HTTP_REFERER)
    print(request.referrer)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    else:
        referrer=''
    #u.netloc
    print(referrer)
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    print "Searched for %s: %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('{}/gene/{}'.format(referrer,identifier))
    elif datatype == 'transcript':
        return redirect('{}/transcript/{}'.format(referrer,identifier))
    elif datatype == 'variant':
        return redirect('{}/variant/{}'.format(referrer,identifier))
    elif datatype == 'region':
        return redirect('{}/region/{}'.format(referrer,identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('{}/dbsnp/{}'.format(referrer,identifier))
    elif datatype == 'hpo':
        return redirect('{}/hpo/{}'.format(referrer,identifier))
    elif datatype == 'mim':
        return redirect('{}/mim/{}'.format(referrer,identifier))
    elif datatype == 'error':
        return redirect('{}/error/{}'.format(referrer,identifier))
    elif datatype == 'not_found':
        return redirect('{}/not_found/{}'.format(referrer,identifier))
    else:
        raise Exception


@app.route('/patient/<patient_str>')
def get_patient(patient_str):
    pass

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
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    pos = int(pos)
    # pos, ref, alt = get_minimal_representation(pos, ref, alt)
    xpos = get_xpos(chrom, pos)
    variant = lookups.get_variant(db, xpos, ref, alt)
    print(variant)
    if variant is None:
        variant = {
            'chrom': chrom,
            'pos': pos,
            'xpos': xpos,
            'ref': ref,
            'alt': alt
        }
    consequences = None
    ordered_csqs = None
    if 'vep_annotations' in variant:
        variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
        ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
        ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
        consequences = defaultdict(lambda: defaultdict(list))
        for annotation in variant['vep_annotations']:
            annotation['HGVS'] = get_proper_hgvs(annotation)
            consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
    base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
    any_covered = any([x['has_coverage'] for x in base_coverage])
    metrics = lookups.get_metrics(db, variant)
    # check the appropriate sqlite db to get the *expected* number of
    # available bams and *actual* number of available bams for this variant
    sqlite_db_path = os.path.join(
        app.config["READ_VIZ_DIR"],
        "combined_bams",
        chrom,
        "combined_chr%s_%03d.db" % (chrom, pos % 1000))
    print(sqlite_db_path)
    try:
        read_viz_db = sqlite3.connect(sqlite_db_path)
        n_het = read_viz_db.execute("select n_expected_samples, n_available_samples from t " "where chrom=? and pos=? and ref=? and alt=? and het_or_hom=?", (chrom, pos, ref, alt, 'het')).fetchone()
        n_hom = read_viz_db.execute("select n_expected_samples, n_available_samples from t " "where chrom=? and pos=? and ref=? and alt=? and het_or_hom=?", (chrom, pos, ref, alt, 'hom')).fetchone()
        read_viz_db.close()
    except Exception, e:
        logging.debug("Error when accessing sqlite db: %s - %s", sqlite_db_path, e)
        n_het = n_hom = None
    read_viz_dict = {
        'het': {'n_expected': n_het[0] if n_het is not None and n_het[0] is not None else -1, 'n_available': n_het[1] if n_het and n_het[1] else 0,},
        'hom': {'n_expected': n_hom[0] if n_hom is not None and n_hom[0] is not None else -1, 'n_available': n_hom[1] if n_hom and n_hom[1] else 0,},
    }
    for het_or_hom in ('het', 'hom',):
        #read_viz_dict[het_or_hom]['some_samples_missing'] = (read_viz_dict[het_or_hom]['n_expected'] > 0)    and (read_viz_dict[het_or_hom]['n_expected'] - read_viz_dict[het_or_hom]['n_available'] > 0)
        read_viz_dict[het_or_hom]['all_samples_missing'] = (read_viz_dict[het_or_hom]['n_expected'] != 0) and (read_viz_dict[het_or_hom]['n_available'] == 0)
        read_viz_dict[het_or_hom]['readgroups'] = [ '%(chrom)s-%(pos)s-%(ref)s-%(alt)s_%(het_or_hom)s%(i)s' % locals() for i in range(read_viz_dict[het_or_hom]['n_available']) ]   #eg. '1-157768000-G-C_hom1', 
        read_viz_dict[het_or_hom]['urls'] = [ os.path.join('combined_bams', chrom, 'combined_chr%s_%03d.bam' % (chrom, pos % 1000)) for i in range(read_viz_dict[het_or_hom]['n_available']) ]
    print 'Rendering variant: %s' % variant_str
    return render_template(
        'variant.html',
        variant=variant,
        base_coverage=base_coverage,
        consequences=consequences,
        any_covered=any_covered,
        ordered_csqs=ordered_csqs,
        metrics=metrics,
        read_viz=read_viz_dict,
    )


@app.route('/variant2/<variant_str>')
def variant_page2(variant_str):
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    print(region)
    records=tb.fetch(region=region)
    geno=dict(zip(headers, [r.split('\t') for r in records][0]))
    samples=[h for h in geno if geno[h].split(':')[0]=='0/1' or geno[h].split(':')[0]=='1/1']
    d=csv.DictReader(file('/data/uclex_files/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
    #d=csv.DictReader(file('/data/UCLpheno/uclex-hpo.txt','r'),delimiter='\t')
    for r in d:
        #if r['eid'] not in samples:
        if r['sample'] not in samples: continue
        print(r['sample'])
        pheno=r['phenotype']
        print(pheno)
        if pheno.startswith('HP'):
            print(phizz.query_hpo([pheno]))
        elif pheno.startswith('MIM'):
            print(phizz.query_disease([pheno]))
    return('\n\n'.join(samples))


@app.route('/chisqu/<variant_str>')
def chisq(variant_str):
    variant_str=str(variant_str).strip().replace('_','-')
    chrom, pos, ref, alt = variant_str.split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    print(region)
    records=tb.fetch(region=region)
    geno=dict(zip(headers, [r.split('\t') for r in records][0]))
    samples=[h for h in geno if geno[h].split(':')[0]=='0/1' or geno[h].split(':')[0]=='1/1']
    #d=csv.DictReader(file('/data/uclex_files/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
    #headers=file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().replace('#','').split('\t')
    #d=csv.DictReader(file('/data/UCLpheno/uclex-hpo.txt','r'),delimiter='\t')
    for r in d:
        #if r['eid'] not in samples:
        if r['sample'] not in samples: continue
        print(r['sample'])
        pheno=r['phenotype']
        print(pheno)
        if pheno.startswith('HP'):
            print(phizz.query_hpo([pheno]))
        elif pheno.startswith('MIM'):
            print(phizz.query_disease([pheno]))
    return('\n\n'.join(samples))



@app.route('/hpo/<hpo_id>')
def hpo_page(hpo_id):
    patients_db=get_db('patients')
    db=get_db()
    rsession=get_R_session()
    hpo_db=get_db('hpo')
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    #patient_ids=[p['external_id'] for p in patients]
    #hpo=phizz.query_hpo([hpo_id])[0]
    hpo_name=str(rsession.eval('hpo$name["%s"]'%hpo_id).tolist()[0])
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    #r=patients_db.hpo.find_one({'hp_id':hpo_id})
    #if r: external_ids=r['external_ids']
    #else: external_ids=[]
    genes=[lookups.get_gene_by_name(db, r['Gene-Name']) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    print('num genes', len(genes))
    #for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id}): print(r)
    #pmids=[r['pmid'] for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id})]
    pmids=[]
    #hpo_patients=rsession.eval('query.hpo("%s")'%hpo_id,)
    #print('num patients',len(hpo_patients))
    #if type(hpo_patients) is str:
        #hpo_patients=[hpo_patients]
    #else:
        #hpo_patients=hpo_patients.tolist()
    ## only return common variants if there are many individuals
    ##rsession.voidEval('common_variants <- common.variants')
    ## private variants (not seen in others in the cohort)
    ##rsession.voidEval('common_variants <- common.variants')
    #variants=rsession.r.private_variants(hpo_patients)
    #if type(variants) is str:
        #variants=[variants]
    #else:
        #variants=variants.tolist()
    #print('num variants',len(variants),)
    #variants=[db.variants.find_one({'variant_id':v.replace('_','-')}) for v in variants[:100]]
    #[variant for variant in lookups.get_variants_in_gene(db, g['gene_id'])]
       #if variant['major_consequence']!='stop_gained': continue
       #print(variant)
       #break
    #print( lookups.get_variants_in_gene(db, 'CNNM4') )
    #vcf_reader = pysam.VariantFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % '22')
    #for record in vcf_reader:
        #for s in external_ids:
            #r=record.samples[s]
            #if 'GT' in r: print(r['GT'])
    return render_template('phenotype.html',hpo_id=hpo_id,hpo_name=hpo_name,external_ids=[],genes=genes,pmids=pmids,individuals=[],variants=[])

# AJAX
# fetch patients iwth hpo term
@app.route('/fetch_hpo',methods=['GET','POST'])
def fetch_hpo():
    if request.method=='POST':
        hpo_ids=request.form['hpo_ids'].strip().split(',')
    else:
        hpo_ids=request.args.get('hpo_ids').strip().split(',')
    hpo_id=hpo_ids[0]
    print('HPO',hpo_id)
    rsession=get_R_session()
    hpo_names=str(rsession.eval('hpo$name["%s"]'%hpo_id).tolist()[0])
    hpo_patients=rsession.eval('query.hpo("%s")'%hpo_id,)
    print('num patients',len(hpo_patients))
    if type(hpo_patients) is str:
        hpo_patients=[hpo_patients]
    else:
        hpo_patients=hpo_patients.tolist()
    res=jsonify(result=hpo_patients)
    return res

# AJAX
# fetch variants private to patients
# That is variants which are only seen in these patients and no one else.
@app.route('/fetch_private_variants',methods=['GET','POST'])
def fetch_private_variants():
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    rsession=get_R_session()
    variants=rsession.r.private_variants(hpo_patients)
    print('private variants', variants)
    if type(variants) is str:
        variants=[variants]
    else:
        variants=variants.tolist()
    print('num of private variants',len(variants),)
    res=jsonify(result=variants)
    return res

# AJAX
# fetches information from db
@app.route('/fetch_variant',methods=['GET','POST'])
def fetch_variant():
    if request.method=='POST':
        variants=request.form['variants'].strip().split(',')
    else:
        variants=request.args.get('variant').strip().split(',')
    db=get_db()
    print(variants)
    req_len=len(variants)
    variants=[db.variants.find_one({'variant_id':variant_id.replace('_','-')}, fields={'_id': False}) for variant_id in variants]
    ans_len=len(variants)
    print(req_len==ans_len)
    res=jsonify(result=variants)
    return res

@app.route('/mim/<mim_id>')
def mim_page(mim_id):
    db=get_db('patients')
    print(str(mim_id))
    patients=[p for p in db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    patient_ids=[p['external_id'] for p in patients]
    print(phizz.query_disease([hpo_id]))
    print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    return render_template('test.html')

@app.route('/patient/<patient_id>')
def patient_page(patient_id):
    db=get_db()
    patients=[p for p in db.patients.find({'external_id': str(patient_id)})]
    print(patients)
    return None

@app.route('/Exomiser/<path:path>')
@requires_auth
def exomiser_page(path):
    #is this user authorized to see this patient?
    return send_from_directory('Exomiser', path)

@app.route('/example/')
@requires_auth
def example():
    return send_from_directory('templates', 'temp-plot.html')


def get_gene_page_content(gene_id,hpo_id=None):
    db = get_db()
    rsession=get_R_session()
    # scrape exac
    b=browser.Browser('exac.broadinstitute.org')
    p=b.get_page('/gene/%s'%gene_id)
    m = re.compile('window.table_variants\s*=\s*(.*)\s*;')
    exac_table_variants=dict()
    if m: exac_table_variants=loads(m.search(p).group(1))
    gene = lookups.get_gene(db, gene_id)
    gene_name=gene['gene_name']
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
    variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
    print('variants_in_gene',len(variants_in_gene))
    # which transcript contains most of the variants
    transcript_counter=Counter([t for v in variants_in_gene for t in v['transcripts'] ])
    print(transcript_counter)
    transcript_with_most_variants=transcript_counter.most_common(1)[0][0]
    print('transcript with most variants',transcript_with_most_variants)
    transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)
    print('transcripts_in_gene',len(transcripts_in_gene))
    # Get some canonical transcript and corresponding info
    #transcript_id = gene['canonical_transcript']
    transcript_id = transcript_with_most_variants
    # if none of the variants are on the canonical transcript use the transcript with the most variants on it
    transcript = lookups.get_transcript(db, transcript_id)
    variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
    #print('variants_in_transcript',len(variants_in_transcript))
    #coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
    #print('coverage_stats',len(coverage_stats))
    add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
    constraint_info = lookups.get_constraint_for_transcript(db, transcript_id)
    #print('constraint_info',constraint_info)
    # 2. get patients with hpo term and patients without
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    #patient_ids=[p['external_id'] for p in patients]
    #hpo=phizz.query_hpo([hpo_id])[0]
    # samples
    headers=frozenset(file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().split('\t'))
    if hpo_id is None:
        hpo="HP:0000001"
        hpo_name='All'
        cases=frozenset()
        print('num cases',len(cases))
    else:
        #cases
        hpo_name=str(rsession.eval('hpo$name["%s"]'%hpo_id).tolist()[0])
        print(hpo_name)
        cases=rsession.eval('query.hpo("%s")'%hpo_id,)
        cases=frozenset(cases.tolist())
        print('num cases',len(cases))
    #controls
    #everyone=rsession.eval('query.hpo("HP:0000001")')
    #everyone=frozenset(everyone.tolist())
    controls=headers-cases
    #everyone=everyone & headers
    #controls=frozenset(everyone) - frozenset(cases)
    print('num controls',len(controls))
    # 3. for each variant tabix,  get counts and chisq
    rsession.voidEval('chisq_test <- chisq.test')
    chrom, pos, ref, alt = variants_in_gene[0]['variant_id'].split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/chr%s.vcf.gz' % chrom,)
    region ='%s:%s-%s' % (str(gene['chrom']), str(gene['start']), str(gene['stop']),)
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip('#').strip().split('\t')
    records=[dict(zip(headers,r.strip().split('\t'))) for r in tb.fetch(region)]
    print(len(records))
    records=dict([('%s-%s-%s-%s' % (r['CHROM'], r['POS'], r['REF'], r['ALT'],),r,) for r in records])
    for i,_, in enumerate(variants_in_gene):
        v=variants_in_gene[i]
        variant_str=v['variant_id']
        print(variant_str)
        variant_str=str(variant_str).strip().replace('_','-')
        chrom, pos, ref, alt = variant_str.split('-')
        v['pos_coding_noutr']= get_xpos(chrom, pos)
        #region=str('%s:%s-%s'%(chrom, pos, int(pos),))
        #records=tb.fetch(region=region)
        #geno=dict(zip(headers, [r.split('\t') for r in records][0]))
        if variant_str not in records:
            v["-log10pvalue"]=0
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
        v['co_wt']=co_wt
        v['ca_wt']=ca_wt
        v['co_mut']=co_mut
        v['ca_mut']=ca_mut
        if ca_mut==0:
            v["-log10pvalue"]=0
        else:
            counts=numpy.array([[ca_mut,ca_wt],[co_mut,co_wt]])
            print(counts)
            stat=dict(rsession.r.chisq_test(counts).astuples())
            v['-log10pvalue']=-math.log10(stat['p.value'])
            #print(chisquare([ca_mut,ca_wt,co_mut,co_wt]))
            #d=csv.DictReader(file('/data/uclex_files/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
            #headers=file('/slms/UGI/vm_exports/vyp/phenotips/uclex_files/current/headers.txt','r').read().strip().replace('#','').split('\t')
            #get EXAC info
        for j,_, in enumerate(exac_table_variants):
            exac_v=exac_table_variants[j]
            if exac_v['variant_id']!=variant_str: continue
            v['EXAC']=exac_v
            v['HGVS']=exac_v['HGVS']
            v['HGVSp']=exac_v['HGVSp']
            v['HGVSc']=exac_v['HGVSc']
            v['major_consequence']=exac_v['major_consequence']
            v['exac_allele_freq']=exac_v['allele_freq']
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
        variant_id=variants_in_gene[i]['variant_id']
        for j,_, in enumerate(variants_in_transcript):
            variant_id2=variants_in_transcript[j]['variant_id']
            if variant_id==variant_id2:
                variants_in_transcript[j]['-log10pvalue']=variants_in_gene[i]['-log10pvalue']
                break
    t=render_template( 'gene.html',
            gene=gene,
            transcript=transcript,
            variants_in_gene=variants_in_gene,
            variants_in_transcript=variants_in_transcript,
            transcripts_in_gene=transcripts_in_gene,
            constraint=constraint_info,
            csq_order=csq_order,
            literature_DOT=literature_DOT,
            simreg_DOT=simreg_DOT,
            hpo_name=hpo_name, coverage_stats=[])
    #cache.set(cache_key, t, timeout=1000*60)
    return t


@app.route('/gene/<gene_id>',methods=['GET'])
def gene_page(gene_id):
    # if gene not ensembl id then translate to
    db=get_db()
    hpo=request.args.get('hpo')
    if not gene_id.startswith('ENSG'): gene_id = lookups.get_gene_by_name(get_db(), gene_id)['gene_id']
    if gene_id in app.config['GENES_TO_CACHE']:
        return open(os.path.join(app.config['GENE_CACHE_DIR'], '{}.html'.format(gene_id))).read()
    else:
        return get_gene_page_content(gene_id,hpo)


@app.route('/transcript2/<transcript_id>')
def transcript_page2(transcript_id):
    db = get_db()
    try:
        transcript = lookups.get_transcript(db, transcript_id)
        cache_key = 't-transcript-{}'.format(transcript_id)
        t = cache.get(cache_key)
        print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
        if t is None:
            gene = lookups.get_gene(db, transcript['gene_id'])
            gene['transcripts'] = lookups.get_transcripts_in_gene(db, transcript['gene_id'])
            variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
            coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
            add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
            t = render_template(
                'transcript.html',
                transcript=transcript,
                transcript_json=json.dumps(transcript),
                variants_in_transcript=variants_in_transcript,
                variants_in_transcript_json=json.dumps(variants_in_transcript),
                coverage_stats=coverage_stats,
                coverage_stats_json=json.dumps(coverage_stats),
                gene=gene,
                gene_json=json.dumps(gene),
                csq_order=csq_order,
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on transcript:', transcript_id, ';Error=', traceback.format_exc()
        abort(404)



@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    db = get_db()
    try:
        transcript = lookups.get_transcript(db, transcript_id)
        cache_key = 't-transcript-{}'.format(transcript_id)
        t = cache.get(cache_key)
        print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
        if t: return t
        variants=[v for v in db.variants.find({'Transcript':str(transcript_id)})]
        genes=list(set([variants['Gene'] for v in variants]))
        print(genes)
        cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on transcript:', transcript_id, ';Error=', traceback.format_exc()
        abort(404)



@app.route('/region/<region_id>')
def region_page(region_id):
    db = get_db()
    try:
        region = region_id.split('-')
        cache_key = 't-region-{}'.format(region_id)
        t = cache.get(cache_key)
        print 'Rendering %sregion: %s' % ('' if t is None else 'cached ', region_id)
        if t is None:
            chrom = region[0]
            start = None
            stop = None
            if len(region) == 3:
                chrom, start, stop = region
                start = int(start)
                stop = int(stop)
            if start is None or stop - start > REGION_LIMIT or stop < start:
                return render_template(
                    'region.html',
                    genes_in_region=None,
                    variants_in_region=None,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    coverage=None,
                    csq_order=csq_order,
                )
            if start == stop:
                start -= 20
                stop += 20
            genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop)
            variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
            t = render_template(
                'region.html',
                genes_in_region=genes_in_region,
                variants_in_region=variants_in_region,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=coverage_array,
                csq_order=csq_order,
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on region:', region_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/dbsnp/<rsid>')
def dbsnp_page(rsid):
    db = get_db()
    try:
        variants = lookups.get_variants_by_rsid(db, rsid)
        chrom = None
        start = None
        stop = None
        print 'Rendering rsid: %s' % rsid
        return render_template(
            'region.html',
            rsid=rsid,
            variants_in_region=variants,
            chrom=chrom,
            start=start,
            stop=stop,
            coverage=None,
            genes_in_region=None,
            csq_order=csq_order,
        )
    except Exception, e:
        print 'Failed on rsid:', rsid, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/not_found/<query>')
def not_found_page(query):
    return render_template(
        'not_found.html',
        query=query
    )


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    return render_template(
        'error.html',
        query=query
    )


@app.route('/downloads')
def downloads_page():
    return render_template('downloads.html')


@app.route('/about')
def about_page():
    return render_template('about.html')


@app.route('/participants')
def participants_page():
    return render_template('about.html')


@app.route('/terms')
def terms_page():
    return render_template('terms.html')


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


@app.route('/faq')
def faq_page():
    return render_template('faq.html')

@app.route('/samples')
def samples_page():
    samples=pandas.read_csv('/slms/UGI/vm_exports/vyp/phenotips/HPO/hpo.txt')
    return render_template('samples.html',samples=samples.to_html(escape=False))


@app.route('/text')
def text_page():
    db = get_db()
    query = request.args.get('text')
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    if datatype in ['gene', 'transcript']:
        gene = lookups.get_gene(db, identifier)
        link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%(chrom)s%%3A%(start)s-%(stop)s" % gene
        output = '''Searched for %s. Found %s.
%s; Canonical: %s.
%s''' % (query, identifier, gene['full_gene_name'], gene['canonical_transcript'], link)
        output += '' if 'omim_accession' not in gene else '''
In OMIM: %(omim_description)s
http://omim.org/entry/%(omim_accession)s''' % gene
        return output
    elif datatype == 'error' or datatype == 'not_found':
        return "Gene/transcript %s not found" % query
    else:
        return "Search types other than gene transcript not yet supported"


@app.route('/read_viz/<path:path>')
def read_viz_files(path):
    full_path = os.path.abspath(os.path.join(app.config["READ_VIZ_DIR"], path))
    # security check - only files under READ_VIZ_DIR should be accsessible
    if not full_path.startswith(app.config["READ_VIZ_DIR"]):
        return "Invalid path: %s" % path
    logging.info("path: " + full_path)
    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_from_directory(app.config["READ_VIZ_DIR"], path)
    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        logging.error(error_msg)
        return error_msg
    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset
    data = None
    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)
    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))
    logging.info("GET range request: %s-%s %s" % (m.group(1), m.group(2), full_path))
    return rv


@app.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response


### all the mongodb reading/writing code


def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = raw_input('This will drop the database and reload. Are you sure you want to continue? [no] ')
    if not confirm.startswith('y'):
        print('Exiting...')
        sys.exit(1)
    all_procs = []
    for load_function in [load_variants_file, load_dbsnp_file, load_base_coverage, load_gene_models, load_constraint_information]:
        procs = load_function()
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))
    [p.join() for p in all_procs]
    print('Done! Loading MNPs...')
    load_mnps()
    print('Done! Creating cache...')
    #create_cache()
    print('Done!')


def load_base_coverage():
    """ """
    def load_coverage(coverage_files, i, n, db):
        coverage_generator = parse_tabix_file_subset(coverage_files, i, n, get_base_coverage_from_file)
        try:
            db.base_coverage.insert(coverage_generator, w=0)
        except pymongo.errors.InvalidOperation, e:
            print(e)
            # handle error when coverage_generator is empty
            pass  
    db = get_db()
    db.base_coverage.drop()
    print("Dropped db.base_coverage")
    # load coverage first; variant info will depend on coverage
    db.base_coverage.ensure_index('xpos')
    procs = []
    coverage_files = app.config['BASE_COVERAGE_FILES']
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    random.shuffle(app.config['BASE_COVERAGE_FILES'])
    for i in range(num_procs):
        p = Process(target=load_coverage, args=(coverage_files, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs
    #print 'Done loading coverage. Took %s seconds' % int(time.time() - start_time)


def load_variants_file():
    def load_variants(sites_file, i, n, db):
        for f in sites_file:
            print(f)
            variants_generator = parse_tabix_file_subset([f], i, n, get_variants_from_sites_vcf)
            try:
                db.variants.insert(variants_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when variant_generator is empty
    db = get_db()
    db.variants.drop()
    print("Dropped db.variants")
    # grab variants from sites VCF
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')
    sites_vcfs = app.config['SITES_VCFS']
    print(sites_vcfs)
    #if len(sites_vcfs) > 1: raise Exception("More than one sites vcf file found: %s" % sites_vcfs)
    procs = []
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    #pdb.set_trace()
    for i in range(num_procs):
        p = Process(target=load_variants, args=(sites_vcfs, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

    #print 'Done loading variants. Took %s seconds' % int(time.time() - start_time)


def load_constraint_information():
    db = get_db()
    db.constraint.drop()
    print 'Dropped db.constraint.'
    start_time = time.time()
    with gzip.open(app.config['CONSTRAINT_FILE']) as constraint_file:
        for transcript in get_constraint_information(constraint_file):
            db.constraint.insert(transcript, w=0)
    db.constraint.ensure_index('transcript')
    print 'Done loading constraint info. Took %s seconds' % int(time.time() - start_time)


def load_mnps():
    db = get_db()
    start_time = time.time()
    db.variants.ensure_index('has_mnp')
    print 'Done indexing.'
    while db.variants.find_and_modify({'has_mnp' : True}, {'$unset': {'has_mnp': '', 'mnps': ''}}): pass
    print 'Deleted MNP data.'
    with gzip.open(app.config['MNP_FILE']) as mnp_file:
        for mnp in get_mnp_data(mnp_file):
            variant = lookups.get_raw_variant(db, mnp['xpos'], mnp['ref'], mnp['alt'], True)
            db.variants.find_and_modify({'_id': variant['_id']}, {'$set': {'has_mnp': True}, '$push': {'mnps': mnp}}, w=0)
    db.variants.ensure_index('has_mnp')
    print 'Done loading MNP info. Took %s seconds' % int(time.time() - start_time)


def load_gene_models():
    db = get_db()
    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    print 'Dropped db.genes, db.transcripts, and db.exons.'
    start_time = time.time()
    canonical_transcripts = {}
    with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE']) as canonical_transcript_file:
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript
    omim_annotations = {}
    with gzip.open(app.config['OMIM_FILE']) as omim_file:
        for fields in get_omim_associations(omim_file):
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)
    dbnsfp_info = {}
    with gzip.open(app.config['DBNSFP_FILE']) as dbnsfp_file:
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
            dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)
    print 'Done loading metadata. Took %s seconds' % int(time.time() - start_time)
    # grab genes from GTF
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                gene['other_names'] = dbnsfp_info[gene_id][1]
            db.genes.insert(gene, w=0)
    print 'Done loading genes. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name_upper')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    print 'Done indexing gene table. Took %s seconds' % int(time.time() - start_time)
    # and now transcripts
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading transcripts. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    print 'Done indexing transcript table. Took %s seconds' % int(time.time() - start_time)
    # Building up gene definitions
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading exons. Took %s seconds' % int(time.time() - start_time)
    start_time = time.time()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    print 'Done indexing exon table. Took %s seconds' % int(time.time() - start_time)
    return []


def load_dbsnp_file():
    db = get_db()
    def load_dbsnp(dbsnp_file, i, n, db):
        if os.path.isfile(dbsnp_file + ".tbi"):
            dbsnp_record_generator = parse_tabix_file_subset([dbsnp_file], i, n, get_snp_from_dbsnp_file)
            try:
                db.dbsnp.insert(dbsnp_record_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when coverage_generator is empty
        else:
            with gzip.open(dbsnp_file) as f:
                db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)
    db.dbsnp.drop()
    db.dbsnp.ensure_index('rsid')
    db.dbsnp.ensure_index('xpos')
    start_time = time.time()
    dbsnp_file = app.config['DBSNP_FILE']
    print "Loading dbsnp from %s" % dbsnp_file
    if os.path.isfile(dbsnp_file + ".tbi"): num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    else:
        # see if non-tabixed .gz version exists
        if os.path.isfile(dbsnp_file):
            print(("WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
                "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
                "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
                "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
                "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz") % locals())
            num_procs = 1
        else:
            raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())
    procs = []
    for i in range(num_procs):
        p = Process(target=load_dbsnp, args=(dbsnp_file, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs
    #print 'Done loading dbSNP. Took %s seconds' % int(time.time() - start_time)
    #start_time = time.time()
    #db.dbsnp.ensure_index('rsid')
    #print 'Done indexing dbSNP table. Took %s seconds' % int(time.time() - start_time)


"""
Get the most recent common ancestor between two sets of hpo terms.
"""
def mrc_hpo():
    hpo_graph=get_hpo_graph()
    db=get_db()
    for var in db.variants.find():
        hpo_anc=[]
        for eid in list(set(var['HET']+var['HOM'])):
            patient=db.patients.find_one({'external_id':eid})
            if not patient: continue
            if 'features' not in patient: continue
            for f in patient['features']:
                fid=f['id']
                if not fid.startswith('HP'): continue
                hpo_anc.append(set(hpo_graph.get_ancestors(fid)))
        if not hpo_anc: continue
        if 'SYMBOL' not in var: continue
        var['ALL_HPO']=list(set(set.union(*hpo_anc)))
        var['SHARED_HPO']=list(set.intersection(*hpo_anc))
        print(var['VARIANT_ID'],var['SYMBOL'],len(var['HET']+var['HOM']),var['SHARED_HPO'],var['ALL_HPO'])
        db.variants.update({'VARIANT_ID':var['VARIANT_ID']},var,upsert=True)


"""
To get pred_score for pubmedBatch
[D/A].each = 10, [P].each = 5, [C].each = 6, [T/B/N].each = -1
splicing/indel = 1000
"""
def get_pred_score(row):
    pred = 0
    if 'Func' in row and re.search('splic', row['Func']) or 'ExonicFunc' in row and re.search(r'stop|frame|del|insert', row['ExonicFunc']):
            pred = 1000
    else:
        for key in row:
            if re.search('Pred', key):
                if row[key] == 'D' or row[key] == 'A':
                    pred += 10
                elif row[key] == 'P':
                    pred += 5
                elif row[key] == 'C':
                    pred += 6
                elif row[key] == 'T' or row[key] == 'B' or row[key] == 'N':
                    pred -= 1
    return pred

"""
for pubmedBatch
check title and abstract is truely relevant. Assign to both this gene and each ref
"""
def scrutinise():
    pass


@app.route('/pubmedbatch/', methods=['GET', 'POST'])
@requires_auth
def pubmedbatch_main():
    # this is the main page of pubmedBatch
    # It allows a user to create a "folder" to hold results in the mangodb. Post request will handle the 'folder' creation
    user = session.get('user')
    db = get_db('pubmedbatch')
    user_db = db.results.find_one({'user_id':user})
    
    if request.method == 'POST':
        # time to create a folder
        # first check if the folder already exists. Yes? pass. No? create
        folder = request.form['create-folder']
        user_folder =  [d['folder_name'] for d in user_db['folder'] if 'folder_name' in d]
        if folder not in user_folder:
            db.results.update({'user_id': user}, {'$push': {'folder': {'folder_name': folder, 'files': []}}})
    else:
        # get folders. If user not exists in pubmed db, create it
        #print(res)
        if not user_db:
            # A new user dropping by...
            print('A "user" is being made in pubmedDB!')
            db.results.insert({'user_id':user,'folder':[]})
        
    # let's render the template  
    user_db = db.results.find_one({'user_id':user})
    #print user_db
    return render_template( 'pubmedbatch_main.html', 
            user = user,
            folders = [d['folder_name'] for d in user_db['folder'] if 'folder_name' in d])


@app.route('/pubmedbatch/<folder>', methods=['GET', 'POST'])
@requires_auth
def pubmedbatch(folder):
    # This is the main business
    user = session.get('user')
    db = get_db('pubmedbatch')
    
    if request.method == 'POST':
        # post. get form data, return JSON
        ##########################################################
        # get form data
        #########################################################
        column = int(request.form['column'])
        Entrez.email = request.form['email']
        AND_term = request.form['AND']
        OR_term = request.form['OR']
        verbose = request.form.get('verbose','')
        csv_file = request.files['csv_upload']
        known_genes = request.files['known_genes'] or ''
        mask_genes = request.files['mask_genes'] or ''
        #########################################################
        # parse known genes and mask genes
        known_genes = known_genes.read().split() if known_genes else ''
        mask_genes = mask_genes.read().split() if mask_genes else ''
        #########################################################
        # read csv. has to make the csv string looking like a filehandle
        csv_content = csv_file.read()        
        csvreader = csv.reader(StringIO.StringIO(csv_content), delimiter=',', quotechar='"')
        # life time from config
        life = app.config['PUBMEDBATCH_LIFE']
        # number of lines?
        line_num = len(re.findall('\n', csv_content))
        # format terms
        OR = OR_term.split()
        OR.sort()
        AND = AND_term.split()
        AND.sort()
        smashed_OR = ['"' + t + '"' + '[ALL FIELDS]' for t in OR]
        smashed_AND = ['"' + t + '"' + '[Title/Abstract]' for t in AND]
        smashed_OR = ' OR '.join(smashed_OR)
        smashed_AND = ' AND '.join(smashed_AND)
        smashed_term = ' AND (' + smashed_OR + ')'
        if smashed_AND:
            smashed_term += ' AND ' + smashed_AND
        ###########################################################
        # it's time to read the csv file
        ###########################################################
        row_num = -1
        header = [] # header
        output = [] # all the results get pushed here
        for row in csvreader:
            row_num += 1
            # read header
            if not row_num:
                header = row
                # add 2 columns after HUGO
                header[column+1:column+1] = ['ref(pubmedID)', 'pubmed_score']
                # add a pred score at the beginning
                header[0:0] = ['pred_score']
                continue
            # read in real data
            gene_name = row[column]
            if gene_name == 'NA':
                continue
            # get rid of any parentheses and their content
            gene_name = re.sub(r'\([^)]*\)?','',gene_name)

            genes={} # storing 
            print gene_name
            return gene_name
            ####################################################
            # first to see if masked, pass
            # else
            #   db.cache
            #       when in db.cache and up-to-date, get the pubmed result
            #       when in db.cache and out-of-date, add new pubmed result
            #       when not in db.cache, search pubmed and store
            ####################################################
            # db.cache's structure:
            #  {['key': '_'.join([gene_name.upper(), ','.join(OR).lower(), ','.join(AND).lower()]), 
            #   'result': {pubmed_result},
            #   'date': now]}
            if re.search(r'\b' + re.escape(gene_name) + r'\b', mask_genes, re.IGNORECASE):
                # masked. don't search, just return the row
                if verbose:
                    continue
                row[column+1:column+1] = [{total_score: 0, results: ['masked']}, 0]
                # add a placeholder for pred_score
                row[0:0] = 0
                ha = {}
                for col in range(len(header)):
                    ha[header[col]] = row[col]
                ha['pred_score'] = get_pred_score(ha)
                output.append(ha)
            else:
                # not masked
                now = time.mktime(time.localtime()) #get time in seconds
                lag = 0 # use it as a flag of how to search. 0 = search; now-saved['date'] = update; 
                term = '_'.join([gene_name.upper(), ','.join(OR).lower(), ','.join(AND).lower()])
                #print term
                # now check if the result in the db is uptodate
                saved = db.cache.find_one({'key': term})
                if saved:
                    lag = now - saved['date']
                    # has record. let's see if it is out of date
                    if lag  <= life:
                        # not out of date. let's use it
                        genes[gene_name] = saved['result']

                if gene_name not in genes:
                    # need to search
                   # handle = Entrez.einfo()
                   # record = Entrez.read(handle)
                   # handle.close()

                    if lag:
                        lag = lag/3600/24 # convert it to days
                        # need to update
                        search_results = Entrez.read(Entrez.esearch(db='pubmed', term='Retinal dystrophy', reldate=lag, datetype='pdat', usehistory='y'))
                    else:
                        # just search
                        search_results = Entrez.read(Entrez.esearch(db='pubmed', term='Retinal dystrophy', usehistory='y'))
                    # now done the search. let's get results
                    attempt = 1
                    while attempt <= 10:
                        try:
                            handle = Entrez.efetch("pubmed",restart=0,retmax=1,retmode="xml", webenv=search_results['WebEnv'], query_key=search_results['QueryKey'])
                        except HTTPError as err:
                            if 500 <= err.code <= 599:
                                print('Received error from server %s' % err)
                            else:
                                print('Something is wrong while efetch..')
                            print('Attempt %i of 10' % attempt)
                            attempt += 1
                            time.sleep(5)
                                
                    records = Entrez.parse(handle)
                    
                    if records:
                        # got something. let's do some calculation
                        print records
                        return 'ok'
                pass

            return 'ok'
            break
        #handle = Entrez.einfo()
        #record = Entrez.read(handle)
        #handle.close()
        #search_results = Entrez.read(Entrez.esearch(db='pubmed', term='Retinal dystrophy', reldate=365, datetype='pdat', usehistory='y'))
        #count = int(search_results['Count'])
        #print('Found %i results' % count)
        #handle = Entrez.efetch("pubmed",restart=0,retmax=1,retmode="xml", webenv=search_results['WebEnv'], query_key=search_results['QueryKey'])
        #records = Entrez.parse(handle)
        #import pprint
        #pp = pprint.PrettyPrinter(indent=10)
        ##return( '\n'.join([record['MedlineCitation']['Article']['ArticleTitle'] for record in records]) )
        #records=[ r for r in records]
        #for r in records:
        #    x=r['MedlineCitation']['Article']['Abstract']['AbstractText']
        #    #print('\n'.join(map(lambda x: str(x), x)))
        #    for i in x:
        #        print('%s:%s' % (i.attributes['Label'], str(i)))
        ##abstract = '\n'.join(abstract)

    else:
        # get. display page
        # First see if folder exists. if not, return error
        user_folders = db.results.find_one({'user_id': user})['folder']
        folders = [d['folder_name'] for d in user_folders if 'folder_name' in d]
        # get AND an OR field
        AND = app.config['PUBMEDBATCH_AND']
        OR = app.config['PUBMEDBATCH_OR']

        #user_folders = user_folders['folder']
        if folder not in folders:
            return "Error: " + folder + " does not exist!" 
        # get the files in the folder to display
        print(d)
        print(folder)
        print(user_folders)
        for d in [user_folders]:
            if d['folder_name'] == folder:
                files = [e['file_name'] for e in d['files'] if 'file_name' in e]
        return render_template( 'pubmedbatch.html',
                home_pubmedbatch = home_pubmedbatch,
                files = files,
                user = user, 
                folder = folder,
                AND = AND,
                OR = OR)

@app.route('/plot/<gene>')
def plot(gene):
    #db = get_db()
    #var=db.variants.find_one({'VARIANT_ID':'3_8775295_C_T'})
    d=csv.DictReader(file('/slms/UGI/vm_exports/vyp/phenotips/CARDIO/assoc_3.csv','r'),delimiter=',')
    x=[i for i, r, in enumerate(d)]
    d=csv.DictReader(file('/slms/UGI/vm_exports/vyp/phenotips/CARDIO/assoc_3.csv','r'),delimiter=',')
    y=[-math.log10(float(r['HCM.chisq.p'])) for r in d]
    print(x)
    print(y)
    d=csv.DictReader(file('/slms/UGI/vm_exports/vyp/phenotips/CARDIO/assoc_3.csv','r'),delimiter=',')
    #layout = dict( yaxis = dict( type = 'log', tickvals = [ 1.5, 2.53, 5.99999 ]), xaxis = dict( ticktext = [ "green eggs", "& ham", "H2O", "Gorgonzola" ], tickvals = [ 0, 1, 2, 3, 4, 5 ]))
    labels=[r['VARIANT_ID'] for r in d]
    layout = Layout( xaxis = dict( ticktext=labels, tickvals=x ), title="p-value plot" )
    #Layout( title="p-value plot")
    plotly.offline.plot({
        "data": [
                Scatter(
                    x=x,
                    y=y
                    )
                ],
        "layout": layout
        }, filename='genes/%s-pvalues.html' % (gene,), auto_open=False)
    return send_from_directory('genes', '%s-pvalues.html' % gene,)


if __name__ == "__main__":
    # use ssl
    # add some common url. Would be good if can generate the url in real time
    home = ''
    home_pubmedbatch = '/pubmedbatch'
    #from OpenSSL import SSL
    # altnerative
    #context = SSL.Context(SSL.SSLv23_METHOD)
    #context = ssl.SSLContext(ssl.PROTOCOL_SSLv23)
    #context = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
    #context.load_cert_chain( '/home/rmhanpo/phenotips.crt', '/home/rmhanpo/phenotips.key' )
    #context.use_privatekey_file('/home/rmhanpo/ssl/phenotips.cs.ucl.ac.uk.key')
    #context.use_privatekey_file('/home/rmhanpo/ssl/phenotips.cs.ucl.ac.uk.key')
    #context.use_privatekey_file('/home/rmhanpo/phenotips.key')
    #context.use_certificate_file('/home/rmhanpo/phenotips.crt')
    # this is now handled by Apache
    app.run(host='0.0.0.0',port=8000,threaded=True)
    # threaded
    #app.run(threaded=True)
    #app.run(host='127.0.0.1',port=8000, debug = True, ssl_context=context)
    #app.run(host='0.0.0.0', port=8000, ssl_context=context)
    #app.run(host='0.0.0.0', port=8000, debug=True, ssl_context='adhoc')
    #app.run(host='0.0.0.0', port=8000, debug=True)
    #app.run(host='127.0.0.1', port=8000, debug=True)
    #app.run(host='0.0.0.0', port=8000, debug=True, ssl_context=('/home/rmhanpo/phenotips.key', '/home/rmhanpo/phenotips.crt'))
    #toolbar=DebugToolbarExtension(app)
    #runner = Runner(app)  # adds Flask command line options for setting host, port, etc.
    #runner.run()





