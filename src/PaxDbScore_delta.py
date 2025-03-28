## PaxDbScore

from collections import defaultdict
import sys, random
import numpy as np
from statistics import mean, stdev

# random.seed(42)

abu_file = sys.argv[1]
interaction_file = sys.argv[2]
score_cutoff = 900
N_iter = 500

def read_abundances(abundance_filepath):
    protein_abu = {}
    with open(abundance_filepath) as f:
        for l in f:
            try:
                prot_id, abu = l.strip().split()[:2]
                abu = float(abu)
                if abu == 0:
                    continue
                protein_abu[prot_id] = abu
            except ValueError:
                return None
            
    return protein_abu
            
def get_zscore(real_score, shuffle_list):
    return abs((real_score - mean(shuffle_list)) /stdev(shuffle_list))

# def fisher_yates_shuffle(list_size):
#     max_idx = list_size - 1
#     list_idx = list(range(list_size))
#     for i in range(list_size):
#         new_i = random.randint(i,max_idx)
#         tmp = list_idx[i] 
#         list_idx[i] = list_idx[new_i]
#         list_idx[new_i] = tmp
        
#     return list_idx

def numpy_shuffle(list_size):
    new_abu_index = np.arange(list_size)
    np.random.shuffle(new_abu_index)
    return new_abu_index
    
def read_interactions(interaction_filepath, score_cutoff, protein_abu):
    p1_p2 = defaultdict(set)
    with open(interaction_filepath) as f:
        for l in f:
            pid1, pid2, score = l.strip().split()
            
            if pid1 in protein_abu and pid2 in protein_abu:
                
                if int(score) < score_cutoff:
                    continue
                
                if pid1 not in p1_p2[pid2]: # no inverse redundancy
                    
                    p1_p2[pid1].add(pid2)
                
    return p1_p2

def calculate_median_ratio(protein_abu, p1_p2):
    ratios = []
    for p1, p1_ps in p1_p2.items():
        for p2 in p1_ps:
            ratios.append(protein_abu[p1] / protein_abu[p2])
    
    return np.median(np.absolute(np.log2(ratios)))
        
protein_abu = read_abundances(abu_file)
if not protein_abu:
    sys.exit('abundance file format problem')
p1_p2 = read_interactions(interaction_file, score_cutoff, protein_abu)
   
real = calculate_median_ratio(protein_abu, p1_p2)
shuffle_list = []
abus = list(protein_abu.values())

for i in range(N_iter):
        
    new_abu_index = numpy_shuffle(len(protein_abu))
    protein_abu = {k: abus[new_abu_index[i]] for i,k in enumerate(protein_abu)}
    shuffle_list.append(calculate_median_ratio(protein_abu, p1_p2))
   

# print(real)
# print(max(shuffle_list))
# print(min(shuffle_list))
print(round(get_zscore(real, shuffle_list),2))

