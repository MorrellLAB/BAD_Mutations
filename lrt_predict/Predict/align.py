#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import tempfile
import subprocess
import os
import time

#   And external libraries here
from Bio import SeqIO

#   Import our helper scripts here
from ..General import parse_input
from ..General import set_verbosity
from ..General import check_modules

class PastaAlign:
    def __init__(self, unaligned_sequences, query_sequence, verbose):
        self.mainlog = set_verbosity.verbosity('Pasta_Align', verbose)
        #   This is file-like object
        self.input_seq = unaligned_sequences
        self.query = query_sequence
        self.output = None
        return

    #   A function to call the pasta alignment
    def pasta_align(self):
        #   Get the base directory of the LRT package, based on where this file is
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the pasta script
        pasta_script = os.path.join(lrt_path, 'Shell_Scripts', 'Pasta_Align.sh')
        #   Check for the presence of the pasta executable
        pasta_path = check_modules.check_executable('run_pasta.py')
        #   Next create a temporary output file
        #pasta_out = tempfile.NamedTemporaryFile(mode='w+t', prefix='LRTPredict_PastaAlign_', suffix='_MSA')
        #   Pasta expects a directory for output. We use the system temp directory
        pasta_out = tempfile.gettempdir()
        #   This is a bit inefficient, maybe, but we only use the tempfile function to get a name
        #pasta_out.close()
        #self.mainlog.debug('Created temporary file with prefix ' + pasta_out.name + ' for holding pasta outputs.')
        #   We make a job name from the time in microseconds
        #   This should be good enough...
        pasta_job = 'pastajob_' + '%.6f' % time.time()
        #   Create the command line
        cmd = ['bash', pasta_script, pasta_path, self.input_seq.name, pasta_out, pasta_job]
        self.mainlog.debug(' '.join(cmd))
        #   Then, we'll execute it
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        #   Then, build the output name
        #   The structure of the Pasta output files is
        #       Jobname.marker001.Unaligned_name.aln
        #       Jobname.tre
        aln_out = pasta_out + '/' + pasta_job + '.marker001.' + self.input_seq.name.split('/')[-1] + '.aln'
        tree_out = pasta_out + '/' + pasta_job + '.tre'
        return (out, err, aln_out, tree_out)
