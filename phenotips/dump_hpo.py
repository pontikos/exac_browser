from __future__ import print_function
import sys
import rest
import datetime

now = datetime.datetime.now()

auth=sys.argv[1]

patients=rest.get_patient(auth)['patientSummaries']

#file(sprintf('uclex_hpo_%d-%d-%d.txt'),)
hpo_file=open('uclex-hpo.txt', 'w+')

print('eid', 'hpo', 'genes', 'solved', sep='\t',file=hpo_file)

for p in patients:
    eid=p['eid']
    print(eid)
    patient=rest.get_patient(auth,eid)
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



