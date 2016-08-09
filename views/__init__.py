
#flask import
from flask import Flask, session
from flask.ext.session import Session
from flask import Response
from flask import stream_with_context, request, Response, make_response
from flask import Flask
from flask import request
from flask import send_file
from flask import session
from flask import g
from flask import redirect
from flask import url_for
from flask import abort
from flask import render_template
from flask import flash
from flask import jsonify
from flask import send_from_directory
from flask.ext.compress import Compress
from flask.ext.runner import Runner
from flask_errormail import mail_on_500
from flask_debugtoolbar import DebugToolbarExtension 
import sys
import StringIO
import urllib, base64 
import numpy as np
from jinja2_extensions import *
import md5 
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
from collections import defaultdict, Counter
from collections import OrderedDict
from werkzeug.contrib.cache import SimpleCache 
from multiprocessing import Process
import glob
import sqlite3
import traceback
import time 
from functools import wraps 
from werkzeug.exceptions import default_exceptions, HTTPException 
import pandas
import csv
import time
import StringIO 
from urlparse import urlparse
import pickle 
#import pdb 
# handles live plotting if necessary
import math
import plotly
print plotly.__version__  # version >1.9.4 required
from plotly.graph_objs import Scatter, Layout 
# connect to R session
#import pyRserve 
import numpy
import subprocess
from flask import Flask, render_template, redirect, url_for, request

from load_individual import load_patient 
from Crypto.Cipher import DES
import base64




logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().setLevel(logging.INFO)

#ADMINISTRATORS = ( 'n.pontikos@ucl.ac.uk',)
#app = Flask(__name__)
#mail_on_500(app, ADMINISTRATORS)
#Compress(app)
#app.config['COMPRESS_DEBUG'] = True
##cache = SimpleCache(default_timeout=60*60*24)



app = Flask(__name__)
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



def check_auth(username, password):
    """
    This function is called to check if a username / password combination is valid.
    Will try to connect to phenotips instance.
    """
    print username
    conn=PhenotipsClient()
    response=conn.get_patient(auth='%s:%s' % (username, password,),number=1)
    if response:
        session['password2'] = password
        password=md5.new(password).hexdigest()
        session['user'] = username
        session['password'] = password
        return True
    else: return False
    # check that user name and hash of password exist in database
    db_users=get_db('users')
    # setting a session key for pubmedBatch to save result
    session['password2'] = password
    password=md5.new(password).hexdigest()
    session['user'] = username
    session['password'] = password
    r=db_users.users.find_one({'user':username})
    if r is None:
        return False
    elif md5.new(r['password']).hexdigest() == md5.new(password).hexdigest():
        print('LOGIN', session['user'])
        return True
    else:
        return False


def authenticate():
    """Sends a 401 response that enables basic auth"""
    return Response( 'Could not verify your access level for that URL.\n' 'You have to login with proper credentials', 401, {'WWW-Authenticate': 'Basic realm="Login Required"'})


def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if session:
          if 'user' in session and 'password2' in session and check_auth(session['user'],session['password2']):
             return f(*args, **kwargs)
          else:
             return redirect('https://uclex.cs.ucl.ac.uk/login')
             #return render_template('login.html', error='Invalid Credentials. Please try again.')
        print 'method', request.method
        error=None
        if request.method == 'POST':
          username=request.form['username']
          password=request.form['password']
          if check_auth(username,password):
             return f(*args, **kwargs)
          else:
             # doesn't redirect
             #return render_template('login.html', error='Invalid Credentials. Please try again.')
             #return login()
             return redirect('https://uclex.cs.ucl.ac.uk/login')
    return decorated


# 
@app.route('/login', methods=['GET','POST'])
def login():
    print request.method
    error = None
    print 'login', request.method
    print request.form
    if request.method == 'POST':
       username=request.form['username']
       password=request.form['password']
       if not check_auth(username,password):
          error = 'Invalid Credentials. Please try again.'
       else:
           return redirect('https://uclex.cs.ucl.ac.uk')
    return render_template('login.html', error=error)

# 
@app.route('/logout')
def logout():
    try:
        print session
        del session['user']
        del session['password']
        del session['password2']
        del session
    except NameError:
        return redirect('https://uclex.cs.ucl.ac.uk/login')
    return render_template('login.html', error="You have been logged out")



@app.route('/')
@requires_auth
def homepage():
    cache_key = 't-homepage'
    #t = cache.get(cache_key)
    #if t: return t
    db=get_db()
    total_variants=db.variants.count()
    print('total_variants',total_variants,)
    total_patients=db.patients.count()
    print('total_patients',total_patients,)
    male_patients=db.patients.find( {'sex':'M'}).count()
    print('male_patients',male_patients,)
    female_patients=db.patients.find( {'sex':'F'}).count()
    print('female_patients',female_patients,)
    unknown_patients=db.patients.find( {'sex':'U'}).count()
    dotfile='static/dot/ENSG00000122375_hom_comp.dot'
    DOT=file(dotfile,'r').read().replace('\n','\\n')
    # replace single quote
    DOT=re.sub("'", '&#39;', DOT)
    #fontsize=7
    # change fontsize to 7
    #DOT=re.sub(r'fontsize="\d+"', 'fontsize="%d"' % fontsize, DOT)
    exac_variants=db.variants.find({'in_exac':True}).count()
    print('exac_variants',exac_variants,)
    pass_variants=db.variants.find({'filter':'PASS'}).count()
    print('pass_variants',pass_variants,)
    pass_exac_variants=db.variants.find({'in_exac':True,'filter':'PASS'}).count()
    print('pass_exac_variants',pass_exac_variants,)
    pass_exac_variants=db.variants.find({'in_exac':True,'filter':'PASS'}).count()
    nonexac_variants=db.variants.find({'in_exac':False}).count()
    pass_nonexac_variants=db.variants.find({'in_exac':False,'filter':'PASS'}).count()
    nonpass_variants=(total_variants-pass_variants)
    nonpass_nonexac_variants=nonexac_variants-pass_nonexac_variants
    #labels = 'PASS', 'non-PASS',
    #sizes =[100*pass_variants/float(total_variants),100*(nonpass_variants)/float(total_variants)]
    #print(sizes)
    #colors = ['yellowgreen', 'red']
    #explode = (0.1, 0)
    #plt.figure(figsize=(5,5))
    #plt.margins(1, 1)
    #plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90)
     ## Set aspect ratio to be equal so that pie is drawn as a circle.
    #plt.axis('equal')
    #plt.axis('off')
    #plt.show()
    # word cloud
    #from os import path
    #from wordcloud import WordCloud
    #text = 'HPO HPO HPO HPO all day'
    ## Read the whole text.
    ## take relative word frequencies into account, lower max_font_size
    #wordcloud = WordCloud().generate(text)
    #plt.figure()
    #plt.imshow(wordcloud)
    #plt.axis("off")
    #plt.show()
    #imgdata = StringIO.StringIO()
    #plt.savefig(imgdata, format='svg')
    #imgdata.seek(0)  # rewind the data
    #import urllib
    #image=urllib.quote(base64.b64encode(imgdata.buf))
    #image=imgdata.buf
    #image = '<svg' + image.split('<svg')[1]
    t = render_template('homepage.html',
        total_patients=total_patients,
        male_patients=male_patients,
        female_patients=female_patients,
        unknown_patients=unknown_patients,
        DOT=DOT,
        total_variants=total_variants,
        exac_variants=exac_variants,
        pass_variants=pass_variants,
        pass_exac_variants=pass_exac_variants,
        pass_nonexac_variants=pass_nonexac_variants,
        #image=image.decode('utf8'))
        image="")
    #cache.set(cache_key, t)
    return t



@app.route('/set/<query>')
def set(query):
    value = query
    session['key'] = value
    return value

@app.route('/get/')
def get():
    return session.get('key', 'not set')


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

def response(POS, REF, ALT, index, geno, chrom, pos):
    homozygous_genotype='/'.join([str(index),str(index)])
    heterozygous_genotype='/'.join(['0',str(index)])
    variant=dict()
    variant['pos']=POS
    variant['ref']=REF
    variant['alt']=ALT
    variant['hom_samples']=[h for h in geno if geno[h].split(':')[0]==homozygous_genotype][0:100]
    variant['HOM_COUNT']=len(variant['hom_samples'])
    variant['het_samples']=[h for h in geno if geno[h].split(':')[0]==heterozygous_genotype][0:100]
    variant['HET_COUNT']=len(variant['het_samples'])
    variant['wt_samples']=[h for h in geno if geno[h].split(':')[0]=='0/0'][1:100]
    variant['WT_COUNT']=len([h for h in geno if geno[h].split(':')[0]=='0/0'])
    variant['MISS_COUNT']=len([h for h in geno if geno[h].split(':')[0]=='./.'])
    variant['allele_num']= 2*(variant['HOM_COUNT'] + variant['HET_COUNT']+variant['WT_COUNT'])
    variant['allele_count']=2*variant['HOM_COUNT'] + variant['HET_COUNT']
    #variant['site_quality'] = variant['QUAL']
    #variant['filter'] = variant['FILTER']
    if variant['WT_COUNT']==0:
        variant['allele_freq'] = None
    else:
        variant['allele_freq'] = float(variant['HET_COUNT']+2*variant['HOM_COUNT']) / float(2*variant['WT_COUNT'])
    var2='-'.join([str(chrom),str(pos),variant['ref'],variant['alt']])
    variant['variant_id']=var2
    samples=variant['het_samples']+variant['hom_samples']
    print(samples)
    variant['hpo']=[p for p in get_db('patients').patients.find({'external_id':{'$in':samples}},{'_id':0,'features':1,'external_id':1})]
    return(jsonify(result=variant))



#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

#from uclex_base import *
#from uclex_gene import *

#def get_db(dbname=None):
#    """
#    Opens a new database connection if there is none yet for the
#    current application context.
#    """
#    if dbname is None: dbname=app.config['DB_NAME']
#    if not hasattr(g, 'db_conn'):
#        g.db_conn=dict()
#        g.db_conn[dbname] = connect_db(dbname)
#    elif dbname not in g.db_conn:
#        g.db_conn[dbname] = connect_db(dbname)
#    return g.db_conn[dbname]



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

# AJAX
# Not finished
@app.route('/chisqu/<variant_str>',methods=['GET','POST'])
def chisq(variant_str):
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
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
    res=jsonify(result=hpo_patients)
    return res


def stream_template(template_name, **context):
    app.update_template_context(context)
    t = app.jinja_env.get_template(template_name)
    rv = t.stream(context)
    rv.enable_buffering(5)
    return rv

@app.route('/my-large-page.html')
def render_large_template():
    rows = iter_all_rows()
    return Response(stream_template('the_template.html', rows=rows))



@app.route('/stream')
def streamed_response():
    def generate():
         yield 'Hello '
         yield request.args['name']
         yield '!'
    return Response(stream_with_context(generate()))

def generate_patient_table():
    def get_variants(variant_ids):
        for vid in db.variants.find({'variant_id':{'$in':variant_ids}}):
            yield 

'''
serve the Vincent annotated csv files
'''
@app.route('/download/send_csv', methods=['GET','POST'])
@requires_auth
def download_csv():
    conn=PhenotipsClient()
    p_id = request.args.get('p_id')
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=p_id,auth=auth)
    if not p: return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    folder = request.args.get('folder')
    path = '/slms/UGI/vm_exports/vyp/phenotips/DROPBOX/'
    csv_file = os.path.join(path,folder, p_id + '.csv')
    filename = folder+'_'+p_id+'.csv'
    if not os.path.isfile(csv_file):
        return 'Oops, file not found!'
    return send_file(csv_file,
                     mimetype='text/csv',
                     attachment_filename=filename,
                     as_attachment=True)

def encrypt(s):
    obj=DES.new(session['password'][:8], DES.MODE_ECB)
    s=s+(8-(len(s) % 8))*' '
    s=obj.encrypt(s)
    s=base64.urlsafe_b64encode(s)
    return s

def decrypt(s):
    obj=DES.new(session['password'][:8], DES.MODE_ECB)
    s=base64.urlsafe_b64decode(str(s))
    s=obj.decrypt(s)
    s=s.replace(' ','')
    return s

# shows each individual, 
# all_individuals
@app.route('/individuals')
@requires_auth
def individuals_page():
    page=int(request.args.get('page',0))
    number=int(request.args.get('number',500))
    hpo_db=get_db('hpo')
    def f(p):
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'homozygous_variants_count':1,'compound_hets_count':1, 'rare_variants_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        return p
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    patients=conn.get_patient(auth=auth,start=page*number,number=number).get('patientSummaries',[])
    eids=[p['eid'] for p in patients]
    print(eids)
    patients=get_db('patients').patients.find({'external_id':{'$in':eids}})
    #patients=get_db('patients').patients.find({'external_id':re.compile('^IRDC')},{'pubmedBatch':0})
    individuals=[f(p) for p in patients if 'external_id' in p]
    # family_history":{"consanguinity":true}
    #if session['user']=='demo': for ind in individuals: ind['external_id']=encrypt(ind['external_id'])
    return render_template('individuals_page.html',individuals=individuals)


@app.route('/research_pubmed', methods=['POST'])
def research_pubmed():
    # use new search terms to update the individual-pubmedbatch table
    patient_id = request.form['p_id']
    search_term = request.form['OR']
    # update patient pubmed result status as running (1)
    db=get_db()
    db.patients.update({'external_id':patient_id},{'$set': {'pubmedbatch.status': 1}})
    # do the actual update
    #exit_status = subprocess.call(['python','individual_pubmedBatch.py', '-p', patient_id, '--OR', search_term])
    exit_status=0
    # reset update status to 0
    db.patients.update({'external_id':patient_id},{'$set': {'pubmedbatch.status': 0}})
    return str(exit_status)

@app.route('/hpo')
def hpo_main():
    # HPO summary page
    # major groups, borrowed from phenotips
    major_groups = {'GROWTH PARAMETERS':['HP:0000256','HP:0000252','HP:0000098','HP:0004322','HP:0004324','HP:0004325','HP:0001508','HP:0001528'],'CRANIOFACIAL':['HP:0001363','HP:0000204','HP:0000175','HP:0001999'],'EYE DEFECTS':['HP:0000505','HP:0000481','HP:0000589','HP:0000593','HP:0000518','HP:0000479','HP:0000587','HP:0000568','HP:0000639','HP:0000486','HP:0000601','HP:0000316'],'EAR DEFECTS':['HP:0000407','HP:0000405','HP:0004467','HP:0000384','HP:0000356','HP:0000359'],'CUTANEOUS':['HP:0000953','HP:0001010','HP:0005306','HP:0011276'],'CARDIOVASCULAR':['HP:0001631','HP:0001629','HP:0001674','HP:0001680','HP:0001636','HP:0001638','HP:0011675'],'RESPIRATORY':['HP:0000776','HP:0002088'],'MUSCULOSKELETAL':['HP:0002652','HP:0002659','HP:0009816','HP:0009824','HP:0100490','HP:0001836','HP:0006101','HP:0001770','HP:0100258','HP:0100259','HP:0001180','HP:0001849','HP:0002650','HP:0000925','HP:0001371','HP:0001762'],'GASTROINTESTINAL':['HP:0002032','HP:0002575','HP:0001543','HP:0001539','HP:0002251','HP:0001396','HP:0002910','HP:0001738','HP:0000819'],'GENITOURINARY':['HP:0000107','HP:0000085','HP:0000069','HP:0000795','HP:0000062','HP:0000047','HP:0000028'],'BEHAVIOR, COGNITION AND DEVELOPMENT':['HP:0001263','HP:0010862','HP:0002194','HP:0000750','HP:0001328','HP:0001256','HP:0002342','HP:0010864','HP:0007018','HP:0000717','HP:0000708'],'NEUROLOGICAL':['HP:0001290','HP:0001250','HP:0001251','HP:0001332','HP:0002072','HP:0001257','HP:0010301','HP:0002011']}
    hpo_freq = lookups.get_hpo_size_freq('hpo_freq.tsv')
    return str(hpo_freq)

@app.route('/hpo/<hpo_id>')
def hpo_page(hpo_id):
    patients_db=get_db('patients')
    db=get_db()
    hpo_db=get_db('hpo')
    patients_db=get_db('patients')
    #patients=[p for p in patients_db.patients.find( { 'features': {'$elemMatch':{'id':str(hpo_id)}} } )]
    print(hpo_id)
    if not hpo_id.startswith('HP:'):
        hpo_id=hpo_db.hpo.find_one({'name':hpo_id})['id'][0]
    print(hpo_id)
    hpo_name=hpo_db.hpo.find_one({'id':hpo_id})['name'][0]
    print('HPO ANCESTORS')
    hpo_ancestors=lookups.get_hpo_ancestors(hpo_db,hpo_id)
    print(len(hpo_ancestors))
    print([h['name'] for h in hpo_ancestors])
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HET' : { '$in': patient_ids }})]))
    #print(len([v['VARIANT_ID'] for v in db.variants.find({'HOM' : { '$in': patient_ids }})]))
    #r=patients_db.hpo.find_one({'hp_id':hpo_id})
    #if r: external_ids=r['external_ids']
    #else: external_ids=[]
    genes=[lookups.get_gene_by_name(db, r['Gene-Name']) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    print('num genes', len(genes))
    #for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id}): print(r)
    #pmids=[r['pmid'] for r in hpo_db.hpo_pubmed.find({'hpoid':hpo_id})]
    patients=lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)
    print('num patients', len(patients))
    pmids=[]
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
    return render_template('phenotype.html',hpo_id=hpo_id,hpo_name=hpo_name,individuals=[str(p['external_id']) for p in patients],genes=genes,pmids=pmids,variants=[])

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
    hpo_db=get_db('hpo')
    patients_db=get_db('patients')
    hpo_patients=[p['external_id'] for p in lookups.get_hpo_patients(hpo_db,patients_db,hpo_id)]
    print('num patients',len(hpo_patients))
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
    db=get_db()
    if len(hpo_patients)==1:
        variants=db.variants.find({'PRIVATE_MUT':hpo_patients})
    else:
        #rsession=get_R_session()
        variants=rsession.r.private_variants(hpo_patients)
        #variants=[]
        print('private variants', variants)
        if type(variants) is str:
            variants=[variants]
        else:
            variants=variants.tolist()
    print('num of private variants',len(variants),)
    res=jsonify(result=variants)
    return res

# AJAX
# fetch common variants to patients
# That is variants which are seen in all these patients.
@app.route('/fetch_common_variants',methods=['GET','POST'])
def fetch_common_variants():
    if request.method=='POST':
        hpo_patients=request.form['patients'].strip().split(',')
    else:
        hpo_patients=request.args.get('patients').strip().split(',')
    print('hpo_patients',hpo_patients,)
    #rsession=get_R_session()
    #variants=rsession.r.common_variants(hpo_patients)
    variants=[]
    print('common variants', variants)
    if type(variants) is str:
        variants=[variants]
    else:
        variants=variants.tolist()
    print('num of common variants',len(variants),)
    res=jsonify(result=variants)
    return res


# AJAX
# fetches variant record from db
@app.route('/fetch_variant',methods=['GET','POST'])
def fetch_variant():
    if request.method=='POST':
        variants=request.form['variants'].strip().split(',')
    else:
        variants=request.args.get('variants').strip().split(',')
    db=get_db()
    req_len=len(variants)
    variant_ids=map(lambda x: x.replace('_','-'),variants)
    variants=[v for v in db.variants.find({'variant_id':{'$in':variant_ids}}, fields={'_id': False})]
    ans_len=len(variants)
    print(req_len==ans_len)
    res=jsonify(result=variants)
    return res


# AJAX
# fetches information from db
@app.route('/variant_count',methods=['GET','POST'])
def variant_count():
    if request.method=='POST':
        external_id=request.form['external_id'].strip()
    else:
        external_id=request.args.get('external_id').strip()
    #rsession=get_R_session()
    #res=jsonify(result={'variant_count':rsession.eval('sum(as.logical(variants[["%s"]]))' % external_id) , 'external_id':external_id})
    #return res

# AJAX
# fetches information from db
@app.route('/private_variant_count',methods=['GET','POST'])
def private_variant_count():
    if request.method=='POST':
        external_id=request.form['external_id'].strip()
    else:
        external_id=request.args.get('external_id').strip()
    db=get_db('patients')
    p=db.patients.find_one({'external_id':external_id})
    if 'PRIVATE_MUT' not in p: private_variant_count=0
    else: private_variant_count=len(p['PRIVATE_MUT'])
    res=jsonify(result={'variant_count': private_variant_count, 'external_id':external_id})
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
    db = get_db('exac')
    db.variants.drop()
    print("Dropped db.variants")
    # grab variants from sites VCF
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')
    db.variants.ensure_index('variant_id')
    #sites_vcfs = app.config['SITES_VCFS']
    sites_vcfs=['/slms/UGI/vm_exports/vyp/phenotips/ExAC/0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz']
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


#progressbar
'''
{
    'random_p_id':{
        'total':456,
        'count':123,
        'status':['running','done']
    },
    ...
}
'''
PROGRESS_BAR = {}

'''
initiate a progress instance
arg: total length of genes
return: progress_id
'''

def init_progress_bar(id,length):
    # check id
    if id in PROGRESS_BAR:
        if PROGRESS_BAR[id]['status'] != 'done':
            return 'the id already exists in PROGRESS_BAR'

    # initialise progress_bar
    PROGRESS_BAR[id] = {
            'total': length,
            'count':0,
            'message': '',
            'status':'running'
    }
    return 0

'''
update progress
arg: {
    id: id, 
    message: message,
    step: 1
    }
default step 1
'''

def update_progress_bar(obj):
    # check if id in PROGRESS_BAR
    if not obj['id'] in PROGRESS_BAR:
        return 'ID does not exist in PROGRESS_BAR'

    # update progress
    if not 'step' in obj:
        obj['step'] = 1
    PROGRESS_BAR[obj['id']]['count'] += obj['step']

    PROGRESS_BAR[obj['id']]['message'] = obj['message']
    # done?
    if PROGRESS_BAR[obj['id']]['count'] == PROGRESS_BAR[obj['id']]['total']:
        PROGRESS_BAR[obj['id']]['status'] = 'done'

'''
kill a progress
'''

def kill_progress_bar(key):
    if key in PROGRESS_BAR:
        del PROGRESS_BAR[key]


'''
to check if an iterable is empty
'''
def peek(iterable):
    try:
        first = next(iterable)
    except RuntimeError:
        return None
    except StopIteration:
        return None
    return first, itertools.chain([first], iterable)


'''
find the freaking PID, Title or Abstract no matter what!
'''
def find_item(obj, key):
    if key in obj:
        return obj[key]
    if isinstance(obj, dict):
        for k in obj:
            if isinstance(obj[k], dict):
                item = find_item(obj[k], key)
                if item is not None:
                    return item
            elif isinstance(obj[k], list):
                for i in obj[k]:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item
    elif isinstance(obj, list):
        for k in obj:
            if isinstance(k, dict):
                item = find_item(k, key)
                if item is not None:
                    return item
            elif isinstance(k, list):
                for i in k:
                    if isinstance(i, str):
                        continue
                    item = find_item(i, key)
                    if item is not None:
                        return item


"""
for pubmedBatch
check title and abstract is truely relevant. Assign to both this gene and each ref
"""
def scrutinise(obj):
    print obj['smashed_all']
    if obj['lag']:
        obj['lag'] = obj['lag']/3600/24 # convert it to days
            # need to update
        search_results = Entrez.read(Entrez.esearch(db='pubmed', term=obj['smashed_all'], reldate=obj['lag'], datetype='pdat', usehistory='y'))
    else:
        # just search
        search_results = Entrez.read(Entrez.esearch(db='pubmed',retmax=50, term=obj['smashed_all'], usehistory='y'))
    # now done the search. let's get results
    count = int(search_results["Count"])
    print count
    results = {'results':[], 'total_score':0}
    # get search content
    attempt = 1
    while attempt <= 10:
        try:
            handle = Entrez.efetch("pubmed",
                                   restart=0,
                                   retmax=50,
                                   retmode="xml",
                                   webenv=search_results['WebEnv'],
                                   query_key=search_results['QueryKey']
                                   )
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print('Received error from server %s' % err)
            else:
                print('Something is wrong while efetch..')
            print('Attempt %i of 10' % attempt)
            attempt += 1
            time.sleep(5)
    record = Entrez.parse(handle)
    if peek(record):
        # got something. let's do some calculation
        for r in record:
            # calculate score
            score = 0
            pid = str(find_item(r, 'PMID'))
            abstract_list = find_item(r, 'AbstractText')
            # parse abstract
            abstract = ''
            if abstract_list:
                for a in abstract_list:
                    if hasattr(a, 'attributes') and 'Label' in a.attributes:
                        abstract = abstract + '<b>' + a.attributes['Label'] + ': </b>'
                        abstract = abstract + a + '<br/>'
                    else:
                        abstract = abstract + a
        
            title = find_item(r, 'ArticleTitle')
            if title:
                score = score + len(obj['reg'].findall(title))
            if abstract:
                score = score + len(obj['reg'].findall(abstract))
        
            # add result to genes[gene_name]
            if score:
                results['results'].append({
                    'id': pid,
                    'title': title,
                    'abstract': abstract,
                    'score': score
                })
                results['total_score'] = results['total_score'] + score
    results['results'] = sorted(results['results'], key=lambda k: k['score'], reverse=True)
    return results

def get_pred_score(obj):
    # for the batch_pubmed route.
    # calculate the pred score
    # [D/A].each = 10, [P].each = 5, [C].each = 6, [T/B/N].each = -1. If there is a splicing/insertion/deletion event, the score is set as 1000. Not given is set as 0
    # ref: https://github.com/plagnollab/DNASeq_pipeline/blob/master/GATK_v2/filtering.md
    pred = 0
    if ('Func' in obj and re.search('splic', obj['Func'])) or ('ExonicFunc' in obj and re.search(r'stop|frame|del|insert', obj['ExonicFunc'])):
        pred = 1000;
    else:
        for k in obj:
            if re.search('Pred', k):
                if obj[k] == 'D' or obj[k] == 'A':
                    pred = pred + 10
                elif obj[k] == 'P':
                    pred = pred + 5
                elif obj[k] == 'C':
                    pred = pred + 6
                elif obj[k] == 'T' or obj[k] == 'B' or obj[k] == 'N':
                    pred = pred - 1
                else:
                    pass
    return pred;


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

""" JINJA2 filer """
def highlight(text, list, myclass):
    # wrap list element in text (case insensitive) with <span>
    # note that gene description has to be split by ','
    #  with class to do highlighting
    for l in list:
        # remove (.*), escape +?.*
        l = re.sub(r'\(.*\)', '', l)
        l = re.sub(r'\+','\\+',l)
        l = re.sub(r'\?','\\?',l)
        l = re.sub(r'\.','\\.',l)
        l = re.sub(r'\*','\\*',l)
        l = re.sub(r'\[.*\]','',l)
        l = re.sub(r'\\', '\\\\',l)
        words = l.split(',')
        for w in words:
            # wrap w with brackets to be a catch group
            text = re.sub(r'(\b%s\b)' % w, r'<span class="%s">\1</span>' % myclass, text, flags=re.I)
    return text
jinja2.filters.FILTERS['highlight'] = highlight

def highlight2(text, kw, myclass):
    # wrap list element in text (case insensitive) with <span>
    # note that gene description has to be split by ','
    #  with class to do highlighting
    # remove (.*), escape +?.*
    for w in kw:
        # wrap w with brackets to be a catch group
        text = re.sub(r'(%s)'%w, r'<span class="%s">\1</span>' % myclass, text, flags=re.I)
    return text
jinja2.filters.FILTERS['highlight2'] = highlight2


@app.route('/load_individual/<individual>')
@requires_auth
def load_individual(individual):
    patient = get_db().patients.find_one({'external_id':individual})
    if 'rare_variants' in patient and type(patient['rare_variants']) is list:
        referrer=request.referrer
        if referrer:
            u = urlparse(referrer)
            referrer='%s://%s' % (u.scheme,u.hostname,)
            if u.port: referrer='%s:%s' % (referrer,u.port,)
        else:
            referrer=''
        url=referrer+'/individual/'+individual
        print(url)
        return redirect(url)
    filename='/slms/UGI/vm_exports/vyp/phenotips/DROPBOX/rare_variants/%s.csv' % individual
    if not os.path.isfile(filename): return '%s not found!' % filename 
    auth='%s:%s' % (session['user'],session['password2'],)
    p = Process(target=load_patient, args=(filename,auth))
    p.start()
    return 'Loading %s...' % individual


import views.uclex_gene
import views.uclex_irdc
import views.gene
import views.variant
import views.individual
import views.pubmedbatch
import views.igv



