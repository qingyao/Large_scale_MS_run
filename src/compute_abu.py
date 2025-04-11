## wrapper to compute abu from fragpipe output
## by the year_month_folder
### Note: fragpipe transforms column name: replace('-','_')
import os
import pandas as pd
from helper import calculate_stringID as cs
import glob, yaml, sys
import subprocess as sbp
from collections import defaultdict
from paxdb import spectral_counting as sc

with open('../rsc/config.yaml') as f:
    config = yaml.safe_load(f)

year_month_folder = sys.argv[1]
if len(sys.argv) > 2 and sys.argv[2] == 'o':
    overwrite = True
else:
    overwrite = False

paxdb_ver = config['PaxDbVersion']
string_ver = config['StringVersion']
fasta_dir = f'../rsc/{paxdb_ver}/fasta'

def write_data_file(data_df, out_column, out_file, vals):

    data_out = pd.DataFrame({'protein ID': data_df.iloc[:,0], out_column: vals/vals.sum() * 10**6})

    data_out= data_out[data_out.iloc[:,1] != 0]
    data_out.dropna(inplace=True)
    data_out.to_csv(out_file, sep = '\t', index = False, float_format = '%.10f', header = None)

    # print(data_df.shape)

def get_label_weights(name, labels, root_folder, pxdid, taxID):
    ### read manifestf
    label_size = {}
    for label in labels:
        curated_fn = os.path.join(root_folder, pxdid, f'{pxdid}-{taxID}-curated.manifest')
        if os.path.isfile(curated_fn):
            manifest_file = curated_fn
        else:
            manifest_file = os.path.join(root_folder, pxdid, f'{pxdid}-{taxID}.manifest')
        res = sbp.run([f'grep {label} {manifest_file} | wc -l'], capture_output = True, shell = True)
        
        if res.returncode == 0:
            label_size[label] = int(res.stdout.decode())
        else:
            label_size[label] = 1 ## equally weighted if manifest file not exist
    return label_size

def compute_data_normalized(label_name, data_df, peptide, out_dir, pxdid, root_folder):
    ## here if multiple group labels are merged into same name
    ## do a weighted average by the number of samples in each group
    data_df['peptide'] = peptide
    data_df.dropna(subset=['peptide'], inplace = True)
    
    name_labels = defaultdict(list)
    for label, name in label_name.items():
        name_labels[name].append(label)    
    for name, labels in name_labels.items():
        out_file = os.path.join(out_dir,f'{pxdid}_{name}.abu')
        if not overwrite and os.path.isfile(out_file):
            continue
        if len(labels) == 1:
            label = labels[0]
            wo_norm =  data_df.loc[:,label.replace('-','_')+' Intensity']
            
        else:
            label_size = get_label_weights(name, labels, root_folder, pxdid, taxID)
            for i, label in enumerate(labels):
                if i == 0:
                    wo_norm = data_df.loc[:,label.replace('-','_')+' Intensity'] * label_size[label]
                else:
                    wo_norm +=  data_df.loc[:,label.replace('-','_')+' Intensity'] * label_size[label]
        
        by_peptide = wo_norm.div(data_df['peptide'], axis = 0)
        
        out_column = name+'_iBAQ_ppm'
        write_data_file(data_df, out_column, out_file, by_peptide)

def write_si_files(label_name, peptide_table, out_folder, root_folder, pxdid, taxID):
    name_labels = defaultdict(list)
    out_files = []
    for label, name in label_name.items():
        name_labels[name].append(label)    
    for name, labels in name_labels.items():
        if len(labels) == 1:
            label = labels[0]
            new_df = peptide_table.loc[:,['Peptide Sequence', label.replace('-','_')+' Intensity']]
        else:
            label_size = get_label_weights(name, labels, root_folder, pxdid, taxID)
            for i, label in enumerate(labels):
                if i == 0:
                    seq_intensity = peptide_table.loc[:,label.replace('-','_')+' Intensity'] * label_size[label]
                else:
                    seq_intensity += peptide_table.loc[:,label.replace('-','_')+' Intensity'] * label_size[label]
            new_df = pd.DataFrame({'Peptide Sequence': peptide_table.iloc[:,0], label.replace('-','_')+' Intensity': seq_intensity})
        new_df = new_df[new_df.iloc[:,1] != 0]
        ## convert everything to int
        min_val = min(new_df.iloc[:,1]) 
        if min_val < 1:
            new_df.iloc[:,1] *= math.ceil(1/min_val)
        new_df.iloc[:,1] = new_df.iloc[:,1].round().astype(int)
        new_df.dropna(inplace=True)
        out_file = os.path.join(out_folder, taxID, f'{pxdid}_{name}.si')
        out_files.append(out_file)
        new_df.to_csv(out_file, sep = '\t', index = False, header = None)
    
    return out_files

root_folder = os.path.join(os.path.pardir, year_month_folder, 'fragpipe_processing_output')
out_folder = os.path.join(os.path.pardir, year_month_folder, 'converted')
os.makedirs(out_folder, exist_ok=True)

for pxdid in os.listdir(root_folder):
    if not pxdid.startswith('PXD'):
        continue
    # if not pxdid == 'PXD000902':
    #     continue
    if len(glob.glob(f'{root_folder}/{pxdid}/*-*.manifest')) > 1:
        multi_taxID_flag = True
    elif len(glob.glob(f'{root_folder}/{pxdid}/*-*.manifest')) == 1:
        multi_taxID_flag = False
    else:
        print(pxdid, 'No manifest file.')
        continue
    
    
    for taxID in os.listdir(os.path.join(root_folder, pxdid)):
        try:
            int(taxID)
        except ValueError:
            continue
        
        ## pxdid, taxID
        print(pxdid, taxID)
        
        label_name_file = os.path.join(root_folder, pxdid, taxID, 'label_name.yaml')
        if not os.path.isfile(label_name_file):
            continue
        with open(label_name_file) as f:
            label_name = yaml.safe_load(f)
        
        
        ### protein_based
        comb_prot_file = os.path.join(root_folder, pxdid, taxID, 'combined_protein.tsv')
        if not os.path.isfile(comb_prot_file):
            continue
        protein_table = pd.read_csv(comb_prot_file, 
                                    sep = '\t', 
                                    usecols=['Protein ID']+[i.replace('-','_')+' Intensity' for i in label_name.keys()]) ### fragpipe converts dash to underscore
        peptide = cs.calculate_peptide(list(protein_table['Protein ID']))
        
        os.makedirs(os.path.join(out_folder, taxID), exist_ok=True)
        compute_data_normalized(label_name, protein_table, peptide, os.path.join(out_folder, taxID), pxdid, root_folder)
        
        
        ### 
        comb_pep_file = os.path.join(root_folder, pxdid, taxID, 'combined_peptide.tsv')
        if not os.path.isfile(comb_pep_file):
            continue
        peptide_table = pd.read_csv(comb_pep_file, 
                                    sep = '\t', 
                                    usecols=['Peptide Sequence']+[i.replace('-','_')+' Intensity' for i in label_name.keys()]) ### fragpipe converts dash to underscore
        
        
        si_files = write_si_files(label_name, peptide_table, out_folder, root_folder, pxdid, taxID)
        
        # compute_data_NSAF
        for si_file in si_files:
            abu_file = si_file.replace('.si', '_si') + '.abu'
            if not os.path.isfile(abu_file) or overwrite:
                sc.calculate_abundance_and_raw_spectral_counts(si_file,
                                                            abu_file, 
                                                            taxID, 
                                                            fasta_dir, 
                                                            string_ver)