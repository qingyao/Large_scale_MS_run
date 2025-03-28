#!/usr/bin/python

# this is for the spectral counting calculating protein abundance
# and raw spectral count and mapped peptide
#

import shlex
import subprocess
import logging
from pathlib import Path

def calculate_abundance_and_raw_spectral_counts(pepfile, scfile, speid, fasta_dir, fasta_ver):
    cmd = 'python3 ' + str(Path(__file__).parent.parent / 'ComputeAbundanceswithSC.py') + " -s '{1}' -f '{2}/fasta.{3}.{0}.fa'"
    cmd = cmd.format(speid, pepfile, fasta_dir, fasta_ver)
    try:
        cmd_out = subprocess.check_output(shlex.split(cmd))
        with open(scfile, "wb") as ofile:
            ofile.write(cmd_out)
        return scfile
    except Exception as e:
        logging.exception('failed to run %s: %s'.format(cmd, e))
    return None
    