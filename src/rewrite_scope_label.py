## a manual check before compute abu
## check if data scope is ok.
## input new filename and tissue assignment
## 
## by the month folder
import os, yaml, glob, shutil, requests
import subprocess as sbp
from collections import defaultdict

month = '202502'
paxdb_ver = 'v6.0'
string_ver = 'v12.0'
pxd_info_file = '../shared/all_until_2024Jun03.tsv'
pxd_pmid_map_file = '../shared/PXDID_PMID_doi.map'
pmid_text_file = '../shared/result-allpxd.tsv'
tissue_terms = f'rsc/{paxdb_ver}/ontology_terms.tsv'
taxID_name_file = f'rsc/{paxdb_ver}/species.{string_ver}.txt' 

pxd_pmid = {}
with open(pxd_pmid_map_file) as f:
    for l in f:
        ll = l.split('\t')
        if ll[1]:
            pxd_pmid[ll[0]] = ll[1]
        else:
            pxd_pmid[ll[0]] = ll[2]

# print(len(pxd_pmid))

available_tissues = set()
with open(tissue_terms) as f:
    for l in f:
        available_tissues.add(l.split('\t')[1])

taxID_name = {}
with open(taxID_name_file) as f:
    next(f)
    for l in f:
        taxID, _, name = l.split('\t')[:3]
        taxID_name[taxID] = name
        

root_folder = os.path.join(month, 'fragpipe_processing_output')
for pxdid in os.listdir(root_folder):
    if not pxdid.startswith('PXD'):
        continue
    
    manifest_fps = glob.glob(f'{root_folder}/{pxdid}/*-*.manifest')
    taxIDs = set()
    if len(manifest_fps) > 0:
        for manifest_fp in manifest_fps:
            taxID = os.path.splitext(os.path.basename(manifest_fp))[0].split('-')[1]
            taxIDs.add(taxID)
    
    else:
        print(pxdid, 'No manifest file.')
        continue
    if len(taxIDs) > 1:
        multi_taxID_flag = True
    else:
        multi_taxID_flag = False
    
    
    for taxID in os.listdir(os.path.join(root_folder, pxdid)):
        try:
            int(taxID)
        except ValueError:
            continue
    
        ## pxdid, taxID
        
        group_files = defaultdict(list)
        ## curated manifest has priority
        curated_manifest_fp = os.path.join(root_folder, pxdid, f'{pxdid}-{taxID}-curated.manifest')
        if os.path.isfile(curated_manifest_fp):
            manifest_fp = curated_manifest_fp
        else:
            manifest_fp = os.path.join(root_folder, pxdid, f'{pxdid}-{taxID}.manifest')
        with open(manifest_fp) as f:
            for l in f:
                fp, group = l.split()[:2]
                fn = os.path.basename(fp)
                group_files[group].append(fn)
                
        ## step 0 filter by number of lines in combined_protein.tsv file
        print('step 0: Quantification check')
        comb_prot_file = os.path.join(root_folder, pxdid, taxID, 'combined_protein.tsv')
        if os.path.isfile(comb_prot_file):
            res = sbp.run(['wc', '-l', os.path.join(root_folder, pxdid, taxID, 'combined_protein.tsv')], capture_output=True)
            protein_count = res.stdout.decode().strip().split()[0]
            if int(protein_count) < 500:
                print(pxdid, protein_count, 'proteins: too few\n')
                continue
            else:
                print(pxdid, protein_count, 'proteins\n')
        else:
            print(pxdid, 'No protein results\n')
            continue
            
        print()
        ## step 1 check if data scope is ok and if the species was correct
        print('step 1: Data scope check')
        scope_file = os.path.join(root_folder, pxdid, 'scope.flag')
        res = sbp.run(['grep', pxdid, pxd_info_file], capture_output=True)
        ll = res.stdout.decode().strip().split('\t')
        print('Title:', ll[1])
        print('Keywords:', ll[-1]) 
        pmid = pxd_pmid[pxdid]
        tmp = sbp.run(['grep', pmid, pmid_text_file], capture_output=True)
        if tmp:
            print('Title of paper: ', '\n'.join(tmp.stdout.decode().split('\n')[0].split('\t')[4:]))
        else:
            print('No paper found')
        print()
        
        if os.path.isfile(scope_file):
            with open(scope_file) as f:
                scope = f.read().strip()
                
            if scope == 'False':
                print(pxdid, 'Data scope was checked negative.')
                continue
        else:
            
            scope_txt = input('Does it fit the PaxDb data scope? [Y]es/[N]o (N is no, default yes): ')
            
            if scope_txt == 'N':
                scope = False
                with open(os.path.join(root_folder, pxdid, 'scope.flag'), 'w') as wf:
                    print('False', file = wf)
                continue
            else:
                scope = True
                with open(os.path.join(root_folder, pxdid, 'scope.flag'), 'w') as wf:
                    print('True', file = wf)
            
        print()
        ## step 2 give new names for groups which pass criteria (control/healthy)
        print('step 2: Group inclusion and rename')
        if not os.path.isfile(os.path.join(root_folder, pxdid, taxID, 'label_name.yaml')):
                
            
            ### first print all group_labels to give an overview
            print(pxdid)
            
            if multi_taxID_flag:
                print(taxID, taxID_name[taxID])
            
            print(len(group_files), 'groups:', ','.join(group_files.keys()))
            
            label_name_mapper = {}
            for group_label, fns in group_files.items():
                print(f'Label: {group_label}')
                for fn in fns:
                    print('    '+fn)
                reject_group = input('    The sample group fits PaxDb? yes[Return]/no[type anything] ')
                if reject_group.strip() != '':
                    print('Rejected the group. Next up...\n')
                    continue
                new_label = input('    Accept label [Return] or type new label? ')
                if new_label:
                    label_name_mapper[group_label] = new_label
                else:
                    label_name_mapper[group_label] = group_label
            with open(os.path.join(root_folder, pxdid, taxID, 'label_name.yaml'), 'w') as wf:
                yaml.dump(label_name_mapper, wf)
            
            print()
        ## step 3 give tissue label on groups which pass the criteria (control/healthy in step 2)
        print('step 3: Tissue assignment')
        if os.path.isfile(os.path.join(root_folder, pxdid, taxID, 'label_tissue.yaml')):
            continue
        
        ### if this project is from PRIDE, it may contain additional info
        res=requests.get(f'https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{pxdid}')
        if res.status_code==200:
            try:
                metadata = res.json()
                
                if 'sampleAttributes' in metadata:
                    for samp_attr in metadata['sampleAttributes']:
                        
                        if samp_attr['key']['name'] == 'organism part':
                            
                            print('PRIDE - organism part:', ','.join([i['name'] for i in samp_attr['value']]))
            except:
                pass
        else:
            print('Did not get any tissue info from PRIDE...')
            
        with open(os.path.join(root_folder, pxdid, taxID, 'label_name.yaml')) as f:
            label_name = yaml.safe_load(f)
        print(label_name)
        
        label_tissue_mapper = {}
        new_tissue_pxd = defaultdict(list)
        for group_label in label_name.keys():
            print(f'Label: {group_label}')
            tissue = input('    Whole organism [Return] or type tissue? ').upper()
            if not tissue.strip():
                tissue = 'WHOLE_ORGANISM'
            if not tissue in available_tissues:
                new_tissue_pxd[tissue].append(pxdid + ':' + group_label)
            label_tissue_mapper[group_label] = tissue
            
        with open(os.path.join(root_folder, pxdid, taxID, 'label_tissue.yaml'), 'w') as wf:
            yaml.dump(label_tissue_mapper, wf)
        
        if len(new_tissue_pxd) > 0:
            new_tissue_file =os.path.join(month, 'new_tissue.yaml')
            if os.path.isfile(new_tissue_file):
                with open(new_tissue_file) as f:
                    existing_yaml = yaml.safe_load(f)
                for tissue, id_lists in existing_yaml.items():
                    new_tissue_pxd[tissue].extend(id_lists)
            with open(os.path.join(month, 'new_tissue.yaml'), 'w') as wf:
                yaml.dump(new_tissue_pxd, wf)

        # for i in range(10):
        print()
        print()
        print()
        print()
        print()
        print()
        print()
        print()
        print()
        
        