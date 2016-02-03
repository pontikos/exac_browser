#!/usr/bin/env python2

from phenotips import rest
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

#import pdb

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
app.config.from_pyfile('uclex.cfg')

GENE_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'gene_cache')
GENES_TO_CACHE = {l.strip('\n') for l in open(os.path.join(os.path.dirname(__file__), 'genes_to_cache.txt'))}


def check_auth(username, password):
    """This function is called to check if a username / password combination is valid.  """
    response=rest.get_patient('%s:%s' % (username, password,))
    if response:
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

def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if not hasattr(g, 'db_conn'): g.db_conn = connect_db()
    return g.db_conn


def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    print(client)
    print(app.config['DB_NAME'])
    return client[app.config['DB_NAME']]


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
    f = open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'), 'w')
    for s in sorted(autocomplete_strings):
        f.write(s+'\n')
    f.close()

    # create static gene pages for genes in
    if not os.path.exists(GENE_CACHE_DIR):
        os.makedirs(GENE_CACHE_DIR)

    # get list of genes ordered by num_variants
    for gene_id in GENES_TO_CACHE:
        try:
            page_content = get_gene_page_content(gene_id)
        except Exception as e:
            print e
            continue
        f = open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id)), 'w')
        f.write(page_content)
        f.close()


def precalculate_metrics():
    import numpy
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
    t = cache.get(cache_key)
    if t is None:
        t = render_template('homepage.html')
        cache.set(cache_key, t)
    return t


@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    #for n in dir(request): print(n, getattr(request,n))
    #print(request.HTTP_REFERER)
    if request.referrer:
        referrer=request.referrer
    else:
        referrer=''
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    print "Searched for %s: %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('{}gene/{}'.format(referrer,identifier))
    elif datatype == 'transcript':
        return redirect('{}transcript/{}'.format(referrer,identifier))
    elif datatype == 'variant':
        return redirect('{}variant/{}'.format(referrer,identifier))
    elif datatype == 'region':
        return redirect('{}region/{}'.format(referrer,identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('{}dbsnp/{}'.format(referrer,identifier))
    elif datatype == 'error':
        return redirect('{}error/{}'.format(referrer,identifier))
    elif datatype == 'not_found':
        return redirect('{}not_found/{}'.format(referrer,identifier))
    else:
        raise Exception


@app.route('/variant/<variant_str>')
def variant_page(variant_str):
    db = get_db()
    try:
        chrom, pos, ref, alt = variant_str.split('-')
        pos = int(pos)
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        xpos = get_xpos(chrom, pos)
        variant = lookups.get_variant(db, xpos, ref, alt)
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
    except Exception:
        print 'Failed on variant:', variant_str, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/variant2/<variant_str>')
def variant_page2(variant_str):
    chrom,pos,ref,alt,=str(variant_str.strip()).split('-')
    tb=pysam.TabixFile('/slms/UGI/vm_exports/vyp/phenotips/mainset_January2015/mainset_January2015_chr%s.vcf.gz' % chrom,)
    region=str('%s:%s-%s'%(chrom, pos, int(pos),))
    headers=[h for h in tb.header]
    headers=(headers[len(headers)-1]).strip().split('\t')
    print(region)
    records=tb.fetch(region=region)
    geno=dict(zip(headers, [r.split('\t') for r in records][0]))
    samples=[h for h in geno if geno[h].split(':')[0]=='0/1' or geno[h].split(':')[0]=='1/1']
    d=csv.DictReader(file('/data/uclex_data/UCLexInfo/uclex-samples.csv','r'),delimiter=',')
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


@app.route('/gene/<gene_id>')
def gene_page(gene_id):
    if gene_id in GENES_TO_CACHE:
        return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    else:
        return get_gene_page_content(gene_id)

def get_gene_page_content(gene_id):
    db = get_db()
    try:
        gene = lookups.get_gene(db, gene_id)
        if gene is None:
            abort(404)
        cache_key = 't-gene-{}'.format(gene_id)
        t = cache.get(cache_key)
        print 'Rendering %sgene: %s' % ('' if t is None else 'cached ', gene_id)
        if t is None:
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
            print('variants_in_transcript',len(variants_in_transcript))
            coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
            print('coverage_stats',len(coverage_stats))
            add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
            constraint_info = lookups.get_constraint_for_transcript(db, transcript_id)
            print('constraint_info',constraint_info)
            #print(gene)
            #print(transcript)
            print(variants_in_gene)
            t=render_template(
                'gene.html',
                gene=gene,
                transcript=transcript,
                variants_in_gene=variants_in_gene,
                variants_in_transcript=variants_in_transcript,
                transcripts_in_gene=transcripts_in_gene,
                coverage_stats=coverage_stats,
                constraint=constraint_info,
                csq_order=csq_order,
            )
            print('all good')
            cache.set(cache_key, t, timeout=1000*60)
        print 'Rendering gene: %s' % gene_id
        #print(phizz.query_gene(ensembl_id=str(gene_id)))
        return t
    except Exception, e:
	print(dir(e))
	print(e.args)
	print(e.message)
	print 'Failed on gene:', gene_id, ';Error=', e
        #abort(404)

@app.route('/gene2/<gene_id>')
def gene_page2(gene_id):
    if gene_id in GENES_TO_CACHE:
        return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    else:
        return get_gene_page_content2(gene_id)


def get_gene_page_content2(gene_id):
    db = get_db()
    try:
        gene = lookups.get_gene(db, gene_id)
        if gene is None: abort(404)
        cache_key = 't-gene-{}'.format(gene_id)
        t = cache.get(cache_key)
        print 'Rendering %sgene: %s' % ('' if t is None else 'cached ', gene_id)
        if t is None:
            variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
            print('variants_in_gene',len(variants_in_gene))
            # which transcript contains most of the variants
            transcript_counter=Counter([t for v in variants_in_gene for t in v['transcripts'] ])
            #which of these variant is missense/frameshift
            #what is the HPO enrichment
            for v in variants_in_gene:
                #print(v['vep_annotation']['Consequence'])
                if v['major_consequence']!='stop_gained': continue
                variant_id='%s-%s-%s-%s' % (v['chrom'],v['pos'],v['ref'],v['alt'],)
                print(lookups.get_hpo(variant_id))
                # [u'allele_count', u'pos', u'quality_metrics', u'variant_id', u'alt', u'pop_homs', u'pop_acs', 'category', u'allele_freq', 'major_consequence', u'vep_annotations', 'HGVSc', u'rsid', u'ref', u'xpos', u'site_quality', u'orig_alt_alleles', u'genes', 'HGVSp', u'hom_count', u'chrom', u'xstart', u'allele_num', u'pop_ans', u'filter', 'flags', u'xstop', 'HGVS', u'transcripts', 'CANONICAL']
                # u'chrom', u'xstart', u'allele_num',
                #print(v.keys())
        return 'done'
    except Exception, e:
	print(dir(e))
	print(e.args)
	print(e.message)
	print 'Failed on gene:', gene_id, ';Error=', e



@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
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
    samples=pandas.read_csv('/data/uclex_data/UCLexInfo/uclex-samples.csv')
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

if __name__ == "__main__":
    # use ssl
    #from OpenSSL import SSL
    # altnerative
    import ssl
    #context = SSL.Context(SSL.SSLv23_METHOD)
    #context = ssl.SSLContext(ssl.PROTOCOL_SSLv23)
    #context = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
    #context.load_cert_chain( '/home/rmhanpo/phenotips.crt', '/home/rmhanpo/phenotips.key' )
    #context.use_privatekey_file('/home/rmhanpo/ssl/phenotips.cs.ucl.ac.uk.key')
    #context.use_privatekey_file('/home/rmhanpo/ssl/phenotips.cs.ucl.ac.uk.key')
    #context.use_privatekey_file('/home/rmhanpo/phenotips.key')
    #context.use_certificate_file('/home/rmhanpo/phenotips.crt')
    # this is now handled by Apache
    context=( '/home/rmhanpo/keys/host.crt', '/home/rmhanpo/keys/host.key')
    context=( '/home/rmhanpo/keys/phenotips.crt', '/home/rmhanpo/keys/phenotips.key')
    app.run(host='0.0.0.0',port=8000)
    #app.run(host='127.0.0.1',port=8000, debug = True, ssl_context=context)
    #app.run(host='0.0.0.0', port=8000, ssl_context=context)
    #app.run(host='0.0.0.0', port=8000, debug=True, ssl_context='adhoc')
    #app.run(host='0.0.0.0', port=8000, debug=True)
    #app.run(host='127.0.0.1', port=8000, debug=True)
    #app.run(host='0.0.0.0', port=8000, debug=True, ssl_context=('/home/rmhanpo/phenotips.key', '/home/rmhanpo/phenotips.crt'))
    #toolbar=DebugToolbarExtension(app)
    #runner = Runner(app)  # adds Flask command line options for setting host, port, etc.
    #runner.run()





