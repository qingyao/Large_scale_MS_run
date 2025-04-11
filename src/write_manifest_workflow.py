## use the parsed pxd group filenames to write manifest and workflow

import json, os, glob, sys, shutil, yaml
from collections import defaultdict
import subprocess as sbp

def write_manifest(pxdid, meta_data_folder, pxd_group_fnames, dataset_filepath, manifest_dirpath):
    os.makedirs(manifest_dirpath, exist_ok=True)

    group_fpath = defaultdict(list)
    group_count = defaultdict(int)
    for group, fnames in pxd_group_fnames[pxdid].items():
        
        # print(group_label)
        for fn in fnames:
            filepath = f'{dataset_filepath}/{pxdid}/raw/{fn}' 
            ## check file exists
            if os.path.isfile(filepath):
                group_fpath[group].append(filepath)
    with open(f'{manifest_dirpath}/{pxdid}.manifest', 'w') as wf:
        for group_label, fpaths in group_fpath.items():
            for fpath in fpaths:
                print(os.path.abspath(fpath), group_label, '', 'DDA', sep ='\t', file = wf)

def write_workflow(default_workflow, pxdid, pxdid_taxid, workflow_dirpath, year_month_folder, paxdb_ver, string_ver):
    
    if pxdid not in pxdid_taxid: 
        return None
    
    missing_fasta = []
    for i, taxid in enumerate(pxdid_taxid[pxdid]):
        fasta_fpath = f'../rsc/{paxdb_ver}/fasta+decoy/fasta.{string_ver}.{taxid}_decoy.fa'
        if os.path.isfile(fasta_fpath):
            line = 'database.db-path='+fasta_fpath
            with open(f'{workflow_dirpath}/{pxdid}-{taxid}.workflow', 'w') as wf:
                print(line, file = wf)
                with open(default_workflow) as f:
                    for l in f:
                        print(l.strip(), file = wf)
        else:
            with open(f'../{year_month_folder}/logs/prepare_fragpipe_workflow.log', 'a') as wf:
                print(pxdid, 'no available fasta file for', taxid, file = wf)
        
            missing_fasta.append(taxid)
    if len(missing_fasta) == len(pxdid_taxid[pxdid]):
        return None
        
    if i > 0:
        return "need manual check"
    else:
        return "only one taxID, can proceed"
   
    
if __name__ == '__main__':
    with open('../rsc/config.yaml') as f:
        config = yaml.safe_load(f)

    paxdb_ver = config['PaxDbVersion']
    string_ver = config['StringVersion']

    metadata_folder = config['MetadataFolder']
    year_month_folder = sys.argv[1]
    raw_download_dir = f'../{year_month_folder}/raw_downloads'
    processing_output_dir = f'../{year_month_folder}/fragpipe_processing_output'
    
    ## output file
    out_file = f'../{year_month_folder}/pxdid_taxid_to_run_fragpipe.txt'
    if os.path.isfile(out_file):
        os.remove(out_file)
        
    pxdids = set()
    
    with open(f'../{year_month_folder}/pxdids.txt') as f:
        for l in f:
            pxdids.add(l.strip())
            
    with open(f'../{year_month_folder}/pxd_group_fnames.yaml') as f:
        pxd_group_fnames = yaml.safe_load(f)
        
    pxdid_taxid_map_filepath = f'../{year_month_folder}/pxdID_taxID.yaml'
    with open(pxdid_taxid_map_filepath) as f:
        pxdid_taxid = yaml.safe_load(f)
    
    # print(len(pxdids))
    # print(len(set(pxd_group_fnames.keys())& pxdids))
    
    for pxdid in pxdids:
        if len(glob.glob(f'{processing_output_dir}/{pxdid}/*/combined_protein.tsv')) != 0:
            print(pxdid, 'done')
            continue
        if not pxdid in pxd_group_fnames:
            continue
        write_manifest(pxdid, f'{metadata_folder}/meta_data', pxd_group_fnames, raw_download_dir, f'{processing_output_dir}/{pxdid}')
        res = write_workflow('../rsc/lfq-mbr.workflow', pxdid, pxdid_taxid, f'{processing_output_dir}/{pxdid}', year_month_folder, paxdb_ver, string_ver)
        print(pxdid)
        if res:
            print(pxdid)
            manifest_file = f'{processing_output_dir}/{pxdid}/{pxdid}.manifest'
            if res == "need manual check":
                
                if os.path.isfile(manifest_file):
                    for workflow_file in glob.glob(f'{processing_output_dir}/{pxdid}/{pxdid}-*.workflow'):
                        taxid = os.path.splitext(os.path.basename(workflow_file))[0].split('-')[1]
                        if os.path.isfile(f'{processing_output_dir}/{pxdid}/{pxdid}-{taxid}-curated.manifest'):
                            continue
                        shutil.copy(manifest_file, f'{processing_output_dir}/{pxdid}/{pxdid}-{taxid}.manifest')
                        
                        manual_check = input(f'Multiple taxIDs for this project. Please edit the manifest for {pxdid} {taxid} and confirm that the files corresponding to taxonID. [C]onfirm/[S]kip the project/taxID if unsure:')
                        if manual_check == 'C':
                            print(pxdid, 'multiple taxIDs:', taxid)
                            shutil.copy(f'{processing_output_dir}/{pxdid}/{pxdid}-{taxid}.manifest', 
                                        f'{processing_output_dir}/{pxdid}/{pxdid}-{taxid}-curated.manifest')
                            # sbp.run(['./run_fragpipe.sh', pxdid, taxid])
                            with open(out_file, 'a') as wf:
                                print(pxdid, taxid+'-curated', sep = '\t', file = wf)
                        else:
                            with open(f'../{year_month_folder}/logs/prepare_fragpipe_manifest.log', 'a') as wf:
                                print(pxdid, 'no dataset for', taxid, file = wf)
            else:
                
                if os.path.isfile(manifest_file):
                    workflow_file = glob.glob(f'{processing_output_dir}/{pxdid}/{pxdid}-*.workflow')[0]
                    if os.path.isfile(workflow_file):
                        taxid = os.path.splitext(os.path.basename(workflow_file))[0].split('-')[1]
                        shutil.copy(manifest_file, f'{processing_output_dir}/{pxdid}/{pxdid}-{taxid}.manifest')
                        
                        with open(out_file, 'a') as wf:
                            print(pxdid, taxid, sep = '\t', file = wf)
                        # sbp.run(['./run_fragpipe.sh', pxdid, taxid])