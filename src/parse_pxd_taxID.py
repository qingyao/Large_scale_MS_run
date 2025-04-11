## takes as input:
### ../rsc/{paxdb_ver}/species.{string_ver}.txt
### ../{year_month_folder}/pxd_group_fnames.yaml

## creates output:
### ../{year_month_folder}/missing_taxID_name.txt
### ../{year_month_folder}/taxID_perl_mapped.tsv
### ../{year_month_folder}/taxID_mapped_final.yaml
### ../{year_month_folder}/pxdID_taxID.yaml
### ../{year_month_folder}/taxID_name.tsv

import json, os, sys
from collections import defaultdict
import yaml
import subprocess as sbp

with open('../rsc/config.yaml') as f:
    config = yaml.safe_load(f)

paxdb_ver = config['PaxDbVersion']
string_ver = config['StringVersion']

year_month_folder = sys.argv[1]
metadata_folder = config['MetadataFolder']

string_taxIDs = set()
string_taxID_name = {}
with open(f'../rsc/{paxdb_ver}/species.{string_ver}.txt') as f:
    for l in f:
        taxID, _, name = l.split('\t')[:3]
        string_taxIDs.add(taxID)
        string_taxID_name[taxID] = name
        
pxdids = set()
with open(f'../{year_month_folder}/pxd_group_fnames.yaml') as f:
    grouped = yaml.safe_load(f)
pxdids = set(grouped.keys())

pxdID_taxID = defaultdict(list)
taxID_name = {}
original_taxIDs = set()
for j in pxdids:
    with open(f'{metadata_folder}/{j}.json') as f:
        json_data = json.load(f)
        if 'species' in json_data:
            for sp in json_data['species']:
                current_taxID = None
                current_name = None
                for term in sp['terms']:
                    if term['name'].endswith('TaxID'):
                        # print(term['value'])
                        current_taxID = term['value']
                        pxdID_taxID[j].append(current_taxID)
                        original_taxIDs.add(current_taxID)
                    elif term['name'].endswith('scientific name'):
                        # print(term['value'])
                        current_name = term['value']
                if current_taxID and current_name:
                    taxID_name[current_taxID] = current_name

missing_taxIDs = set()
for i in original_taxIDs:
    if i not in string_taxIDs:
        missing_taxIDs.add(i)

with open(f'../{year_month_folder}/missing_taxID_name.txt', 'w') as wf:
    for i in missing_taxIDs:
        print(i, taxID_name[i], sep = '\t', file = wf)
        
### run mapping
sbp.run(f'perl ../rsc/map_taxid.pl ../{year_month_folder} ../{year_month_folder}/missing_taxID_name.txt', shell = True)

### process mapping output
old_new = defaultdict(list)
with open(f'../{year_month_folder}/taxID_perl_mapped.tsv') as f:
    for l in f:
        new, old, score = l.split()[:3]
        old_new[old].append((new,float(score))) 

mapper_final = {}
for old, new_score in old_new.items():
    best = sorted(new_score, key = lambda x: x[1])[-1]
    mapper_final[old] = best[0]

with open(f'../{year_month_folder}/taxID_mapped_final.yaml', 'w') as wf:
    yaml.dump(mapper_final, wf)
    
## update taxID
new_pxdID_taxID = defaultdict(list)
for pxdID, taxIDs in pxdID_taxID.items():
    for taxID in taxIDs:
        if taxID in missing_taxIDs:
            if taxID in mapper_final:
                new_pxdID_taxID[pxdID].append(mapper_final[taxID])
            ## if not, the taxID won't be added
        else:
            new_pxdID_taxID[pxdID].append(taxID)
            
with open(f'../{year_month_folder}/pxdID_taxID.yaml','w') as wf:
    yaml.dump(dict(new_pxdID_taxID), wf)
    
final_taxIDs = (original_taxIDs - missing_taxIDs) | set(mapper_final.keys())
with open(f'../{year_month_folder}/taxID_name.tsv', 'w') as wf:
    
    for taxID, name in taxID_name.items():
        if taxID in final_taxIDs:
            if taxID in missing_taxIDs:
                new_taxID = mapper_final[taxID]
                print(new_taxID, string_taxID_name[new_taxID], sep = '\t', file = wf)
            else:
                print(taxID, name, sep = '\t', file = wf)