#download_metadata.py

import urllib.request as rq
import os, time, sys, yaml

year_month_folder = sys.argv[1]
pxdid_list_fp = sys.argv[2]

with open('../rsc/config.yaml') as f:
    config = yaml.safe_load(f)

metadata_folder = config['MetadataFolder']

## read file to get pxdids
pxdids = []
with open(pxdid_list_fp) as f:
    for l in f:
        pxdids.append(l.strip())

downloaded = set()            
for i in os.listdir(metadata_folder):
    
    if i.startswith('PXD'):
        downloaded.add(os.path.splitext(i)[0])

os.makedirs(f'{year_month_folder}/logs/', exist_ok = True)
err_f = open(f'../{year_month_folder}/logs/metadata_download_error.log','a')

# PXD_ID = 'PXD020105'
for i, PXD_ID in enumerate(pxdids):
    if PXD_ID in downloaded:
        continue
    print(i)
    with rq.urlopen(f'http://central.proteomexchange.org/cgi/GetDataset?ID={PXD_ID}&outputMode=json') as response:
        res = response.read()
        # data = json.loads(res.decode('UTF-8'))
        try:
            with open(f'{metadata_folder}/{PXD_ID}.json','w') as wf:
                print(res.decode('latin-1'), file = wf)
            
        except Exception as err:
            print(PXD_ID, err, sep = '\t', file = err_f)

        time.sleep(0.1)
