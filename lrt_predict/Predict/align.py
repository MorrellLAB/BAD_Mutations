#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import tempfile
import subprocess
import os

#   And external libraries here
from Bio import SeqIO

#   Import our helper scripts here
from ..General import parse_input
from ..General import set_verbosity
from ..General import check_modules

class PrankAlign:
    def __init__(self, unaligned_sequences, query_sequence, model, verbose):
        self.mainlog = set_verbosity.verbosity('Prank_Align', verbose)
        #   This is file-like object
        self.input_seq = unaligned_sequences
        self.query = query_sequence
        self.model = model
        self.output = None
        return

    #   A function to call the prank alignment
    def prank_align(self):
        #   Get the base directory of the LRT package, based on where this file is
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the prank script
        prank_script = os.path.join(lrt_path, 'Shell_Scripts', 'Prank_Align.sh')
        #   Check for the presence of the prank executable
        prank_path = check_modules.check_executable('prank')
        #   Next create a temporary output file
        prank_out = tempfile.NamedTemporaryFile(mode='w+t', prefix='LRTPredict_PrankAlign_', suffix='_MSA')
        #   This is a bit inefficient, maybe, but we only use the tempfile function to get a name
        prank_out.close()
        self.mainlog.debug('Created temporary file with prefix ' + prank_out.name + ' for holding prank outputs.')
        #   Create the command line
        cmd = ['bash', prank_script, prank_path, self.input_seq.name, self.model, prank_out.name]
        self.mainlog.debug(' '.join(cmd))
        #   Then, we'll execute it
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        return (out, err, prank_out)
