## calculate molecular weight from string ID
## input [stringIDs,]
## output [MWs,]
import os, re
from pathlib import Path

string_version = 'v12.0'
paxdb_version = 'v6.0'

fp = Path(__file__).parent.parent.parent / 'rsc' / 'aa_mw_table.tsv'
aa_mw = {}
with fp.open() as fr:
    for l in fr:
        aa, mw = l.strip().split()
        aa_mw[aa] = int(mw)

def get_species(input_ids):
    return input_ids[0].split('.')[0]

def get_fasta_fn(species_id):
    fn = Path(__file__).parent.parent.parent /'rsc' / paxdb_version / 'fasta' / f'fasta.{string_version}.{species_id}.fa'
    return str(fn)

def get_sequence(input_ids):
    species_id = input_ids[0].split('.')[0]
    fasta_fn = get_fasta_fn(species_id)

    id_seq = {}  
    with open(fasta_fn) as fr:
        for l in fr:
            if l.startswith('>'):
                string_id = l.strip().replace('>', '')
                if string_id in input_ids:
                    wanted = 1
                else:
                    wanted = 0
            elif wanted:
                if string_id in id_seq:
                    id_seq[string_id] += l.strip()
                else:
                    id_seq[string_id] = l.strip()
    
    return id_seq

def calculate_MW(input_ids):
    id_seq = get_sequence(input_ids)
    
    output_mw = []
    for string_id in input_ids:
        
        output_mw.append(sum([aa_mw[i] - 19 if i in aa_mw else 0 for i in id_seq[string_id]] ))

    return output_mw


def calculate_peptide(input_ids, digestion = 'complete'):
    id_seq = get_sequence(input_ids)
    
    output_pep = []
    for string_id in input_ids:
        # print(string_id)
        fragment_detectable = 0
        if string_id in id_seq:
            fragments = [i + 'B' for i in re.split('K|R', id_seq[string_id])] ## cut site added back for the length
            fragments[-1] = fragments[-1][:-1] # remove X
            # print(fragments)
            prev_fragseq = []
            prev_frag = []
            
            if digestion == 'partial':
                for frag in fragments:
                    if len(frag) < 6:
                        
                        prev_frag.append(len(frag))
                        prev_fragseq.append(frag)
                        while sum(prev_frag) > 30:
                            prev_frag.pop(0)
                            prev_fragseq.pop(0)

                        fragment_detectable += len([i for i in range(len(prev_frag)-1) if sum(prev_frag[i:]) >= 6])
                        # print('condition1', frag)
                        # for i in range(len(prev_frag)-1):
                        #     if sum(prev_frag[i:]) >= 6:
                                # print(''.join(prev_fragseq[i:]))


                    elif len(frag) <= 30:
                        fragment_detectable += 1
                        prev_frag.append(len(frag))
                        prev_fragseq.append(frag)
                        while sum(prev_frag) > 30:
                            prev_frag.pop(0)
                            prev_fragseq.pop(0)

                        fragment_detectable += len([i for i in prev_frag[:-1]])
                        # print('condition2', frag)
                        # for i in range(len(prev_frag)):
                            # print(''.join(prev_fragseq[i:]))

                    else:
                        prev_frag = []
                    
                output_pep.append(fragment_detectable)
            
            else: ## complete
                output_pep.append(sum([1 for i in fragments if len(i) >= 6 and len(i) <= 30])+1)

            
        else:
            output_pep.append(None)
            
    return output_pep

if __name__ == '__main__':
    
    # print(calculate_peptide(['9544.ENSMMUP00000001055']))
    # print(calculate_peptide(['9544.ENSMMUP00000001055'], digestion='partial'))
    print(calculate_peptide(['9606.ENSP00000463058', '9606.ENSP00000332604']))
    print(calculate_MW(['100226.SCO3207']))
    # print(calculate_peptide(['9606.ENSP00000451828','9606.ENSP00000258149','9606.ENSP00000269305']))
    # print(calculate_peptide(['9606.ENSP00000451828','9606.ENSP00000258149','9606.ENSP00000269305'], digestion = 'partial'))
    # print(calculate_peptide(['4513.MLOC_48833.1']))
    # print(calculate_peptide(['4513.MLOC_48833.1'], digestion = 'partial'))