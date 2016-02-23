
from __future__ import print_function
import phizz
import hpo_similarity
from hpo_similarity.ontology import Ontology
import sys
import json
from StringIO import StringIO 
import pymongo


mongo_host='localhost'
mongo_port=27017
client = pymongo.MongoClient(host=mongo_host, port=int(mongo_port))
db=client['patients']
db.hpo.drop()
db.hpo.ensure_index('hp_id')



path = '/slms/UGI/vm_exports/vyp/phenotips/uclex_files/hp.obo'
ontology = Ontology(path)
hpo_graph = ontology.get_graph()

#print('\t'.join(['external_id','mode_of_inheritance','hpo_features','hpo_ancestors']))
for l in sys.stdin.readlines():
      io = StringIO(l)
      p=json.load(io)
      external_id=p['external_id']
      mode_of_inheritance= [i['id'] for i in p.get('global_mode_of_inheritance',[])]
      hpo=[]
      for f in p.get('features',[]):
          if f['id'].startswith('MIM'):
              hpo+=[x['hpo_term'] for x in phizz.query_disease([f['id']])]
          else:
              hpo+=[f['id']]
      if hpo:
          hpo_ancestors = list(set(sum([list(hpo_graph.get_ancestors(h)) for h in hpo],[])))
      else:
          hpo_ancestors = []
      #print(external_id,','.join(mode_of_inheritance),','.join(hpo),','.join(hpo_ancestors),sep='\t')
      for h in list(set(hpo+hpo_ancestors)):
          q=db.hpo.find_one({'hp_id':h})
          if q is not None:
              external_ids=q['external_ids']
              external_ids=list(set([external_id]+external_ids))
              r=db.hpo.update({'hp_id':h}, {'$set':{'external_ids':external_ids}},upsert=True)
              print('update',h,external_ids)
          else:
              external_ids=[external_id]
              db.hpo.insert({'hp_id':h,'external_ids':external_ids})
              print('insert',h,external_ids)




