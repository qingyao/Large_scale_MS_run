## instead of using clustering and extract from group labels
## directly extracted words from file names and grouped by the words
## from longest word if not shared by all then use it to define groups

import numpy as np
import pandas as pd
import json, os, re, yaml, sys
from collections import defaultdict
import inflection
import spacy
# from nltk.corpus import wordnet

nlp = spacy.load("en_core_web_md")
# def is_wordnet_word(token):
#     """Check if the token is a valid word in WordNet"""
#     return bool(wordnet.synsets(token))

def is_spacy_word(token):
    return nlp.vocab[token].has_vector

def sort_fname_by_ext(file_names, file_extensions):
    
    category_string = {ext:set() for ext in file_extensions}
    for ext in file_extensions:
        pn = ''
        for n in file_names:
            if os.path.splitext(n)[1].lower() == ext: ## sometimes it's .RAW
                
                # n_base = os.path.splitext(n)[0]
                
                category_string[ext].add(n)
                
    return category_string

if __name__ == '__main__':
    
    year_month_folder = sys.argv[1]
    pxdid_list_fp = sys.argv[2]
    logfp = f'../{year_month_folder}/logs/metadata_parse.log'
    if os.path.isfile(logfp):
        os.remove(logfp)
        
    with open('../rsc/config.yaml') as f:
        config = yaml.safe_load(f)

    metadata_folder = config['MetadataFolder']
    ## read file to get pxdids
    pxdids = []
    with open(pxdid_list_fp) as f:
        for l in f:
            pxdids.append(l.strip())
    
    pxd_files = {}
    for pxdid in pxdids:
        
        with open(f'{metadata_folder}/{pxdid}.json') as f:
            try:
                json_data = json.load(f)
            except:
                print(pxdid)
                continue
            if 'datasetFiles' in json_data:
                
                file_names = [os.path.basename(i['value']) for i in json_data['datasetFiles']]
                file_extensions = set([os.path.splitext(i['value'])[1] for i in json_data['datasetFiles']])

                fnames = sort_fname_by_ext(file_names, ['.raw'])
                
                pxd_files[pxdid] = fnames['.raw']
                
    # print(pxd_files)
    # pxd_group_differentiating_words = {}
    pxd_groups = {}
    all_tokens = set()
    for pxd, files in pxd_files.items():
        if len(files)>60: ## make it faster in testing
            continue
        tokens = set()
        print(pxd)
        token_to_check = set()
        for fname in files:
            for token in nlp(inflection.underscore(fname.translate({ord(x):'_' for x in '+- .'}).translate({ord(x):'' for x in '0123456789'})).replace('_',' ')):
                t=token.text
                if is_spacy_word(t) and len(t)>2 and t not in ['fraction', 'dia', 'raw'] or t in ['wt','ctrl','ko']:
                    # print(t)
                    tokens.add(t)
                    token_to_check.add(t)

        token_file_count = {}
        for t in token_to_check:
            file_count = sum([1 if t in fname.lower() else 0 for fname in files])
            if file_count < len(files):
                token_file_count[t] = file_count
                tokens.add(t)
        all_tokens.update(tokens)
        tokens = list(tokens) #preserve order
        group_fname = defaultdict(list)
        for fname in files:
            group_component = [] 
            for t in tokens:
                if t in fname.lower():
                    group_component.append(t)
                    
            group_fname['_'.join(group_component)].append(fname)
        pxd_groups[pxd] = group_fname
        # pxd_group_differentiating_words[pxd] = token_file_count
        
    parseable_pxd = set()
    pxd_group_fnames = defaultdict(dict)
    with open(f'../{year_month_folder}/parsed.txt', 'w') as wf:
        for pxd, word_count in pxd_groups.items():
            total_parsable_file_count = 0
            if len(word_count) > 0:
                print(pxd, len(pxd_files[pxd]), file = wf)
                for g, fnames in pxd_groups[pxd].items():
                    if g != '':
                        print(' ', g, len(fnames), file = wf)
                        total_parsable_file_count+= len(fnames)
                        pxd_group_fnames[pxd][g] = []
                        for fname in fnames:
                            pxd_group_fnames[pxd][g].append(fname.replace('+', '').replace(' ', ''))
                            
                # print(sorted([(k,v) for k,v in word_count.items()], key=lambda x: x[1], reverse=True), file = wf)
            if sum([1 if 'TMT' in i else 0 for i in pxd_files[pxd]]) > 1:
                ## exclude processing TMT experiments (no correction matrix)
                with open(logfp, 'a') as logf:
                    print(pxd, 'excluded due to TMT', file = logf)
                
            elif total_parsable_file_count > 0.5* len(pxd_files[pxd]) and len(pxd_groups[pxd]) > 1: ## if more than one group and more than half samples are recognized
                parseable_pxd.add(pxd)
                
            elif len(pxd_groups[pxd]) == 1 and total_parsable_file_count > 0.2* len(pxd_files[pxd]) and  total_parsable_file_count < len(pxd_files[pxd]): # if only one group at least 1/5 of all samples but not all samples. Rarely anything passes this criteria
                with open(logfp, 'a') as logf:
                    print(pxd, 'one group:', list(pxd_group_fnames[pxd].keys())[0], 'with', total_parsable_file_count, 'file out of total', len(pxd_files[pxd]), file = logf)
                parseable_pxd.add(pxd)
    
    new_pxd_group_fnames = {}
    for i in parseable_pxd:
        new_pxd_group_fnames[i]=pxd_group_fnames[i]
                    
    with open(f'../{year_month_folder}/pxd_group_fnames.yaml', 'w') as wf:
        yaml.dump(new_pxd_group_fnames, wf)
        
    with open(f'../{year_month_folder}/parsable_pxd.txt', 'w') as wf:
        for i in parseable_pxd:
            print(i, file = wf)
    
    ## todo, include the token score
    # with open(f'../{year_month_folder}/all_tokens.txt', 'w') as wf:
    #     for t in all_tokens:
    #         print(t, file = wf)