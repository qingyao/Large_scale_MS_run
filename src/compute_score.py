# compute score
## by the year_month_folder
import os, yaml, sys
from paxdb import scores

with open('../rsc/config.yaml') as f:
    config = yaml.safe_load(f)

year_month_folder = sys.argv[1]
paxdb_ver = config['PaxDbVersion']
string_ver = config['StringVersion']
links_folder = f'../rsc/{paxdb_ver}/links'
interaction_fn_template = '{}.network_'+ string_ver + '_900.txt'

for r, _, fs in os.walk(os.path.join(os.path.pardir, year_month_folder, 'converted')):
    for f in fs:
        if f.endswith('.abu'):
            taxID = r.split(os.path.sep)[-1]
            out_fn = os.path.splitext(f)[0] + '.zscores'
            out_fp = os.path.join(r,out_fn)
            in_fp = os.path.join(r, f)
            
            if os.path.isfile(out_fp):
                if os.path.getmtime(out_fp) > os.path.getmtime(in_fp):
                    print(f, 'score has been computed.')
                    continue
            
            print(f)
            scores.score_dataset(in_fp, 
                                 out_fp, 
                                 os.path.join(links_folder, interaction_fn_template.format(taxID)))