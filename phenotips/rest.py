#!/bin/env python
from __future__ import print_function
import re
import os.path
import sys 
import httplib
import urllib
import urlparse
import re
import socket
import gzip
import StringIO
from StringIO import StringIO
import time
import random
from binascii import b2a_base64, a2b_base64
#from email.utils import fix_eols
import json
import pandas
from browser import Browser

b=Browser('localhost:8080',debug=True,print_requests=True)

# get patient with eid or all patients if not
# specified
def get_patient(auth,eid=None):
    auth=b2a_base64(auth).strip()
    headers={'Authorization':'Basic %s'%auth, 'Accept':'application/json'}
    if not eid:
        p=b.get_page('/rest/patients?number=10000', headers=headers)
        io = StringIO(p)
        try:
            d=json.load( io )
        except:
            return None
    else:
        p=b.get_page('/rest/patients/eid/%s'%eid, headers=headers)
        io = StringIO(p)
        try:
            d=json.load( io )
        except:
            return None
    return d

def patient_exists(auth,eid=None):
    p=get_patient(auth,eid)
    if p is None:
        return False
    else:
        return True

def get_vocabularies(auth):
    auth=b2a_base64(auth).strip()
    # get vocabularies
    headers={'Authorization':'Basic %s'%auth, 'Accept':'application/json; application/xml'}
    p=b.get_page('/rest/vocabularies', headers=headers)
    print(p)
    #io = StringIO(p)
    #d=json.load( io )
    #print(d)
    return p

# create patient
def create_patient(auth, patient):
    headers={'Authorization':'Basic %s'% b2a_base64(auth).strip(), 'Content-Type':'application/json', 'Accept':'application/json'}
    io=StringIO()
    json.dump(patient,io)
    json_patient=io.getvalue()
    p=b.get_page('/rest/patients', headers=headers, post=json_patient)
    print(p)

def update_patient(ID, auth, patient):
    io=StringIO()
    json.dump(patient,io)
    json_patient=io.getvalue()
    if patient_exists(auth,ID):
        print('update')
        print(patient)
        headers={'Authorization':'Basic %s'% b2a_base64(auth).strip(),'Content-Type':'application/json', 'Accept':'application/json'}
        p=b.get_page('/rest/patients/eid/%s'%ID, headers=headers, post=json_patient, special='PUT')
        print(p)
    else:
        print('create')
        print(patient)
        create_patient(auth,patient)

def delete_patient(ID, auth, patient):
    auth=b2a_base64(auth).strip()
    headers={'Authorization':'Basic %s'%auth, 'Content-Type':'application/json', 'Accept':'application/json'}
    io=StringIO()
    json.dump(patient,io)
    json_patient=io.getvalue()
    p=b.get_page('/rest/patients/eid/%s'%ID, headers=headers, post=json_patient, special='DELETE')
    print(p)

def update_phenotips_from_csv(info, owner,password):
    info=pandas.read_csv(info)
    print(info.columns.values)
    for i, r, in info.iterrows():
        if r['owner']!=owner: continue
        patient=dict()
        auth='%s:%s' % (r['owner'],password,)
        #auth=login
        #auth=b2a_base64(auth).strip()
        patient['external_id']=r['sample']
        ethnicity=r['ethnicity']
        gender={'0':'U','1':'M','2':'F'}[str(r['gender'])]
        if isinstance(ethnicity,str):
            patient["ethnicity"]={"maternal_ethnicity":[ethnicity],"paternal_ethnicity":[ethnicity]}
        else:
            patient["ethnicity"]={"maternal_ethnicity":[],"paternal_ethnicity":[]}
        patient["prenatal_perinatal_history"]={}
        patient["prenatal_perinatal_phenotype"]={"prenatal_phenotype":[],"negative_prenatal_phenotype":[]}
        patient['reporter']=r['owner']
        patient['sex']=gender
        patient['solved']='unsolved'
        patient['contact']={ "user_id":r['owner'], "name":r['owner'], "email":r['email'], "institution":'' }
        patient['clinicalStatus']={ "clinicalStatus":"affected" }
        #patient['disorders']=[ { "id":r['phenotype'], 'label':''} ]
        patient['features']=[ { "id":r['phenotype'], 'label':'', 'type':'phenotype', 'observed':'yes' } ]
        #update_patient(ID=r['sample'],auth=auth,patient=patient)
        #delete_patient(ID=r['sample'],auth=auth,patient=patient)
        # if patient exists, update patient, otherwise create patient
        update_patient(patient['external_id'],auth,patient)


def dump_hpo_to_csv(outFile, owner, password):
    auth='%s:%s' % (owner, password,)
    patients=get_patient(auth)['patientSummaries']
    #file(sprintf('uclex_hpo_%d-%d-%d.txt'),)
    hpo_file=open(outFile, 'w+')
    print('eid', 'hpo', 'genes', 'solved', sep='\t',file=hpo_file)
    for p in patients:
        eid=p['eid']
        print(eid)
        patient=get_patient(auth,eid)
        print(patient)
        if 'features' in patient:
            hpo=','.join([f['id'] for f in patient['features']])
        else:
            hpo=''
        if 'genes' in patient:
            genes=','.join([g['gene'] for g in patient['genes']])
        else:
            genes=''
        solved=patient['solved']['status']
        print(eid, hpo, genes, solved, sep='\t',file=hpo_file)


def dump_hpo_to_json(outFile, owner, password):
    auth='%s:%s' % (owner, password,)
    patients=get_patient(auth)['patientSummaries']
    #hpo_file=open(outFile, 'w+')
    for p in patients:
        eid=p['eid']
        print(eid)
        patient=get_patient(auth,eid)
        io=StringIO()
        json.dump(patient,io)
        json_patient=io.getvalue()
        print(json_patient)


def dump_to_mongodb(user, password):
    import pymongo
    client = pymongo.MongoClient(host='localhost', port=27017)
    db=client['exac']
    db.patients.drop()
    db.patients.ensure_index('external_id')
    db.patients.ensure_index('report_id')
    auth='%s:%s' % (user, password,)
    patients=get_patient(auth)['patientSummaries']
    for p in patients:
        eid=p['eid']
        print(eid)
        p=get_patient(auth,eid)
        db.patients.insert(p,w=0)


