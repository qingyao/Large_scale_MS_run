## dataset_id	dataset	taxid	organ_ontology	organ	num_proteins  score	origin_version
import hashlib, os, yaml, sys
import subprocess as sbp
from collections import defaultdict

year_month_folder = sys.argv[1]
data_folder = f'../{year_month_folder}/converted'
with open('../rsc/config.yaml') as f:
    config = yaml.safe_load(f)

paxdb_ver = config['PaxDbVersion']
string_ver = config['StringVersion']

ontology_mapper = {}
## from last version
with open(f'../rsc/{paxdb_ver}/ontology_terms.tsv') as f:
    for l in f:
        onto_id, term = l.strip().split('\t')[:2]
        ontology_mapper[term] = onto_id

## newly added
with open(f'../rsc/{paxdb_ver}/new_ontology_terms.tsv') as f:
    for l in f:
        onto_id, term = l.strip().split('\t')[:2]
        ontology_mapper[term] = onto_id

out_folder = os.path.join(year_month_folder, 'converted')

data_info = defaultdict(dict)
for taxID in os.listdir(data_folder):
    if not taxID.isdigit():
        continue

    for fname in os.listdir(os.path.join(data_folder, taxID)):
        file_basename, ext = os.path.splitext(fname)
        pxdid = file_basename.split('_')[0]
        name = '_'.join(file_basename.split('_')[1:])

        dataset_id = int.from_bytes(hashlib.sha256((paxdb_ver+'-'+file_basename).encode('utf-8')).digest()[:4], 'little')
        data_info[dataset_id]['taxid'] = taxID

        if name.endswith('_si'):
            data_info[dataset_id]['si'] = True
        else:
            data_info[dataset_id]['si'] = False

        if ext == '.abu':
            data_info[dataset_id]['dataset'] = fname
            r = sbp.run(['wc','-l', os.path.join(data_folder, taxID, fname)], capture_output = True)
            line_number = r.stdout.decode('utf-8').split()[0]
            data_info[dataset_id]['num_proteins'] = line_number
        elif ext == '.zscores':
            with open(os.path.join(data_folder, taxID, fname)) as f:
                zscore = f.readline().strip()
                data_info[dataset_id]['zscore'] = zscore

        with open(os.path.join(os.path.pardir, year_month_folder, 'fragpipe_processing_output', pxdid, taxID, 'label_tissue.yaml')) as f:
            label_tissue = yaml.safe_load(f)

        with open(os.path.join(os.path.pardir, year_month_folder, 'fragpipe_processing_output', pxdid, taxID, 'label_name.yaml')) as f:
            label_name = yaml.safe_load(f)

        for l, n in label_name.items():
            if n == name.replace('_si',''):
                label = l
                break
        else:
            label = None

        if label:
            tissue = label_tissue[label].upper()
            data_info[dataset_id]['organ'] = tissue
            data_info[dataset_id]['organ_ontology'] = ontology_mapper[tissue]
            data_info[dataset_id]['exclude'] = False
        else:
            data_info[dataset_id]['exclude'] = True

## by default include dataset with _si
## unless both coverage and score is better.
for dataset_id in list(data_info.keys()):
    if dataset_id in data_info:
        info = data_info[dataset_id]
        if not info['si']:
            ## other_id is processed with _si
            other_id = int.from_bytes(hashlib.sha256((paxdb_ver+'-'+os.path.splitext(info['dataset'])[0]+'_si').encode('utf-8')).digest()[:4], 'little')
            
            si_score = data_info[other_id]['zscore']
            si_num = data_info[other_id]['num_proteins']
            if info['zscore'] > si_score and info['num_proteins'] > si_num:
                del data_info[other_id]
            else:
                del data_info[dataset_id]


with open(f'../{year_month_folder}/dataset_information.tsv', 'w') as wf:
    print('#dataset_id', 'dataset', 'taxid', 'organ_ontology', 'organ', 'num_proteins', 'score', 'origin_version', sep = '\t', file = wf)
    for dataset_id, info in data_info.items():
        if info['exclude']:
            continue
        print(dataset_id, info['dataset'], info['taxid'], info['organ_ontology'], info['organ'], info['num_proteins'], info['zscore'], paxdb_ver, sep = '\t', file = wf)