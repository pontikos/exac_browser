from views import *
from lookups import *
import requests
import re
from utils import *
import itertools
import pysam
import csv
#hpo lookup
import orm


@app.route('/individual/<individual>')
@requires_auth
def individual_page(individual):
    #if session['user']=='demo': individual=decrypt(str(individual))
    # make sure that individual is accessible by user
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=individual,auth=auth)
    if not p: return 'Sorry you are not permitted to see this patient, please get in touch with us to access this information.'
    db=get_db()
    hpo_db=get_db('hpo')
    patient_db=get_db('patients')
    patient = db.patients.find_one({'external_id':individual})
    patient2 = patient_db.patients.find_one({'external_id':individual})
    if patient2 is None:
        referrer=request.referrer
        if referrer:
            u = urlparse(referrer)
            referrer='%s://%s' % (u.scheme,u.hostname,)
            if u.port: referrer='%s:%s' % (referrer,u.port,)
        else:
            referrer=''
        url=referrer+'/load_individual/'+individual
        print(url)
        return redirect(url)
    patient['report_id']=patient2['report_id']
    patient['features']=patient2.get('features',[])
    patient['sex'] = patient2['sex']
    patient['family_history'] = patient2.get('family_history',[])
    hpo_ids=[f['id'] for f in patient['features'] if f['observed']=='yes']
    # TODO
    # mode of inheritance in hpo terms: HP:0000005
    #print lookups.get_hpo_children(hpo_db, 'HP:0000005')
    patient['global_mode_of_inheritance']=patient2.get('global_mode_of_inheritance',None)
    # minimise it
    hpo_ids = lookups.hpo_minimum_set(hpo_db, hpo_ids)
    hpo_terms = [(i, hpo_db.hpo.find_one({'id':i})['name'][0]) for i in hpo_ids]
    # this has missing HPO ids. see IRDC_batch2_OXF_3001 and #HP:0000593
    hpo_gene=dict()
    for hpo_id,hpo_term, in hpo_terms:
        hpo_gene[hpo_id] = []
        for gene_name in [x['Gene-Name'] for x in hpo_db.ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.find({'HPO-ID':hpo_id},{'Gene-Name':1,'_id':0})]:
            #gene_hpo[gene_name]=gene_hpo.get(gene_name,[])+[{'hpo_id':hpo_id,'hpo_term':hpo_term}]
            hpo_gene[hpo_id]=hpo_gene.get(hpo_id,[])+[gene_name]
    for k in hpo_gene: hpo_gene[k]=list(frozenset(list(hpo_gene[k])))
    print '========'
    print hpo_gene
    print '========'
    # get pubmedbatch table
    pubmedbatch = patient.get('pubmedbatch',{})
    # candidate genes
    patient['genes'] = patient2.get('genes',[])
    # solved genes
    patient['solved'] = patient2.get('solved',[])
    genes = {}
    # is this still updating?
    update_status = pubmedbatch.get('status', 0);
    # get known and retnet genes
    known_genes = open('ret_known_genes.txt', 'r').readline().strip().split()
    RETNET  = json.load(open('retnet.json', 'r'))
    # get combinatorics of features to draw venn diagram
    feature_combo = []
    feature_venn = []
    for i in range(len(hpo_terms)):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    #venn_ind = -1
    for combo in feature_combo:
        # construct features_venn key
        #venn_ind += 1
        dic_key = '","'.join([hpo_terms[i][1] for i in combo])
        dic_key = '"' + dic_key + '"'
        for ind in range(len(combo)):
            if ind == 0:
                x=hpo_terms[combo[ind]][0]
                feature_venn.append({'key': dic_key, 'value':frozenset(hpo_gene.get(x,""))})
            else:
                tem = feature_venn[-1]['value']
                feature_venn[-1]['value'] = feature_venn[-1]['value'] & frozenset(hpo_gene[hpo_terms[combo[ind]][0]])
    gene_info=dict()
    if 'rare_variants' not in patient or type(patient['rare_variants']) is not list:
        referrer=request.referrer
        if referrer:
            u = urlparse(referrer)
            referrer='%s://%s' % (u.scheme,u.hostname,)
            if u.port: referrer='%s:%s' % (referrer,u.port,)
        else:
            referrer=''
        url=referrer+'/load_individual/'+individual
        print(url)
        return redirect(url)
    for v in patient['rare_variants']:
        if 'HUGO' not in v: v['HUGO']=''
        gene=v['HUGO'].upper() 
        gene_info[gene]=dict()
        if gene in known_genes: gene_info[gene]['known']=True
        if gene not in RETNET: continue
        gene_info[gene]['disease'] = RETNET[gene]['disease']
        gene_info[gene]['omim'] = RETNET[gene]['omim']
        gene_info[gene]['mode'] = RETNET[gene]['mode']
    genes['homozygous_variants']=[v.get('HUGO','').upper() for v in patient['homozygous_variants']]
    genes['compound_hets']=[v.get('HUGO','').upper() for v in patient['compound_hets']]
    genes['rare_variants']=[v.get('HUGO','').upper() for v in patient['rare_variants']]
    genes_pubmed=dict()
    for v in patient['rare_variants']:
        hugo=v['HUGO']
        genes_pubmed[hugo]=get_db('pubmedbatch').cache.find_one( {'key':re.compile(hugo+'[_ ].*')} )
    # figure out the order of columns from the variant row
    table_headers=re.findall("<td class='?\"?(.*)-cell'?\"?>",file('templates/variant_row.tmpl','r').read())
    print table_headers
    return render_template('individual.html', 
            external_id = individual,
            patient=patient,
            table_headers=table_headers,
            pubmedbatch=pubmedbatch,
            pubmed_db=get_db('pubmed_cache'),
            features = hpo_terms,
            genes = genes,
            hpo_gene = hpo_gene,
            gene_info=gene_info,
            genes_pubmed = genes_pubmed,
            update_status = update_status,
            feature_venn = feature_venn)


@app.route('/individual_update/<individual>')
def individual_update(individual):
    conn=PhenotipsClient()
    auth='%s:%s' % (session['user'],session['password2'],)
    p=conn.get_patient(eid=individual,auth=auth)
    print 'UPDATE'
    print p
    print get_db('patients').patients.update({'external_id':individual},{'$set':p},w=0)
    if request.referrer:
        referrer=request.referrer
        u = urlparse(referrer)
        referrer='%s://%s' % (u.scheme,u.hostname,)
        if u.port: referrer='%s:%s' % (referrer,u.port,)
    return redirect(referrer+'/individual/'+individual)


