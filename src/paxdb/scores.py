#!/usr/bin/python3
#
# Computes z-scores for all datasets. Requires libstatistics-descriptive-perl package installed (ubuntu)!
# 
# For each dataset:
# runs PaxDbScore*.pl 3 times, then takes the median
#
# Input files: 
# 1) abundance file: ie ../9615/GPM_2012_09_Canis_familiaris.sc)
# 2) interaction file: ie ../interactions/protein.links.v9.0.CanisFamiliaris_dog_9615_900.txt
# 
# Output files:
#  * for each abundance file - a new one with zscores_ prepended to the name (ie zscores_file.SC)
# 
####
# Notes:
#
# Some input files had rows with only a single field (protein, no counts).
# To check if there are any use this command:
#   for i in *txt; do echo $i; awk '{ if (NF < 2) {printf ("%s\n", $1);exit;}}' $i; done
#
# To filter these lines out use this command:
#
#   awk '{ if (NF >= 2) printf ("%s\t%s\n", $1, $2);}' troublesome_file.txt > tmp; mv tmp troublesome_file.txt
#
# as a test I compared lines in troublesome_file.txt with 
#   $(cat troublesome_file.txt  | cut -f 2 | sort | grep -v '^$' | wc -l) 
# and looks like it works fine
# 
# 
# Originally written by Gabi in bash
# rewritten in python by Milan Simonovic <milan.simonovic@imls.uzh.ch>
#
import glob

import logging
import os
import subprocess
import sys
from pathlib import Path

def score_dataset(dataset_file, output_file, interactions_file):
    logging.debug('abundance_file: %s', dataset_file)
    scores = compute_scores(dataset_file, interactions_file)
    if scores == None:
        return
    #TODO update google doc!
    write_scores(dataset_file, scores, output_file)


def compute_scores(d, interactions_file):
    # cmd_args = ['perl', '-w', 'PaxDbScore_delta.pl', d, interactions_file]
    cmd_args = ['python3',  str(Path(__file__).parent.parent / 'PaxDbScore_delta.py'), d, interactions_file]
    scores = []
    for i in range(1, 4):  # @UnusedVariable
        
        # cmd_out = subprocess.check_output(cmd_args) ### stderr not captured for logging.
        rs = subprocess.run(cmd_args, capture_output = True)
        if rs.returncode == 0:
            cmd_out = rs.stdout.decode('utf-8')
            scores.append(float(cmd_out.strip()))
        else:
            logging.error('failed to score %s, cmd: %s, error: %s', d, " ".join(cmd_args), rs.stderr.decode('utf-8')) # sys.exc_info()[1]
            return None
    logging.info('score(%s): %s', d, str(scores))
    return scores


def write_scores(datafile, scores, outfilepath):
    if len(scores) != 3:
        logging.warning('%s scores missing, 3 expected: %s', datafile, str(scores))
        return

    with open(outfilepath, 'w') as output:
        # output.write('\n'.join(map(lambda x: str(x), sorted(scores))))
        output.write(str(sorted(scores)[1]) + '\n')  # median


class DatasetSorter():
    def __init__(self):
        pass

    def _get_dataset_score(self, name, scores):
        if name in scores:
            return scores[name]
        logging.warning('no score for %s', name)
        return -sys.maxsize - 1  # just put it at the end

    @staticmethod
    def _read_scores(scores_folder):
        scores = dict()
        for scores_file in glob.glob(scores_folder + "/*.zscores"):
            with open(scores_file) as f:
                line = f.readline().strip()
                # without extension
                scores[os.path.basename(os.path.splitext(scores_file)[0])] = float(line)
        if len(scores) == 0:
            raise Exception('no scores found in ' + scores_folder)
        return scores

    def sort_datasets(self, scores_folder):
        scores = self._read_scores(scores_folder)
        return self.sort(scores)

    def sort(self, scores):
        datasets = scores.keys()
        by_score = sorted(datasets, key=lambda d: self._get_dataset_score(d, scores))
        by_score.reverse()
        return by_score


if __name__ == "__main__":
    
    with open('../../rsc/config.yaml') as f:
        config = yaml.safe_load(f)

    paxdb_version = config['PaxDbVersion']

    s = DatasetSorter()
    for species in os.listdir('../../202502/converted/'):
        if not species.isdigit():
            continue
        d = s.sort_datasets('../../202502/converted/' + species)
        print(d)

    print('done')
