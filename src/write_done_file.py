## write done file
import sys, os

year_month_folder = sys.argv[1]

done_pxd = []
for i in os.listdir(f'../{year_month_folder}/fragpipe_processing_output'):
    if i.startswith('PXD'):
        done_pxd.append(i)
        
with open(f'../done_pxds/{year_month_folder}.txt', 'w') as wf:
    for i in sorted(done_pxd):
        print(i, file = wf)