import os, json, sys, yaml
import subprocess as sbp

with open('../rsc/config.yaml') as f:
    config = yaml.safe_load(f)

metadata_folder = config['MetadataFolder']

year_month_folder = sys.argv[1]
done_pxds = set()
for fn in os.listdir('../done_pxds'):
    if fn.endswith('txt'):
        with open('../done_pxds/' + fn) as f:
            for l in f:
                done_pxds.add(l.strip())

if not os.path.isfile(f'../{year_month_folder}/download_urls.txt'):
    download_urls = []
    with open(f'../{year_month_folder}/pxd_group_fnames.yaml') as f:
        pxd_group_fnames = yaml.safe_load(f)
        
    for pxdid in pxd_group_fnames:
        if pxdid in done_pxds:
            continue
        with open(f'{metadata_folder}/{pxdid}.json') as f:
            
            json_data = json.load(f)
            if 'datasetFiles' in json_data:
                for i in json_data['datasetFiles']:
                    lk = i['value'] 
                    if os.path.splitext(lk)[1].lower() == '.raw':
                        if 'pride-archive' in lk and int(lk.split('/')[4]) > 2020:
                            
                            lk = lk.replace('pride-archive','pride/data/archive')
                        lk = lk.replace('+', '') ## remove the + sign, as in PXD002403
                        lk = lk.replace(' ', '') ## remove space, as in PXD006598 and others
                        filename = os.path.basename(lk)
                        download_urls.append(f'{pxdid}/raw/{filename}\t{lk}') 
                            
    with open(f'../{year_month_folder}/download_urls.txt', 'w') as f:
        for i in download_urls:
            print(i, file = f)
        
os.makedirs(f'../{year_month_folder}/logs/', exist_ok = True)
with open(f'../{year_month_folder}/logs/download.log', 'a') as log_file:
    sbp.run(['./helper/download.sh', f'../{year_month_folder}/raw_downloads', f'../{year_month_folder}/download_urls.txt'], stdout=log_file, stderr=sbp.STDOUT)
