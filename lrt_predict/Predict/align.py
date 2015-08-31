#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import tempfile
import subprocess
import os
import time

#   Import our helper scripts here
from ..General import parse_input
from ..General import set_verbosity
from ..General import check_modules


class PastaAlign:
    def __init__(self, unaligned_sequences, query_sequence, verbose):
        self.mainlog = set_verbosity.verbosity('Pasta_Align', verbose)
        #   This is file-like object
        self.input_seq = unaligned_sequences
        #   This will be populated with sequences for back-translation
        self.input_dict = {}
        self.query = query_sequence
        return

    def prepare_sequences(self):
        """Prepares the CDS sequences for alignment with Pasta. Checks if any
        sequences are not multiples of 3, and appends N if not. Translates
        the nucleotide sequences to amino acids for protein alignment with
        Pasta, then removes the trailing stop codon, if it is present."""
        pass

    def back_translate(self):
        """Back-translates from amino acid to nucleotide, using the original
        input sequences as a guide to avoid ambiguity. Assumes that a non-gap
        character in the amino acid alignment will be faithfully represented
        by a triplet in the source sequence, and will not check identity of
        translated codons."""
        pass

    def pasta_align(self):
        """Align the amino acid sequences with Pasta."""
        #   Get the base directory of the LRT package
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the pasta script
        pasta_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'Pasta_Align.sh')
        #   Check for the presence of the pasta executable
        pasta_path = check_modules.check_executable('run_pasta.py')
        #   Pasta expects a directory for output. We use the system temp dir
        pasta_out = tempfile.gettempdir()
        #   We make a job name from the time in microseconds
        #   This should be good enough...
        pasta_job = 'pastajob_' + '%.6f' % time.time()
        #   Create the command line
        cmd = [
            'bash',
            pasta_script,
            pasta_path,
            self.input_seq.name,
            pasta_out,
            pasta_job]
        self.mainlog.debug(' '.join(cmd))
        #   Then, we'll execute it
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        #   Then, build the output name
        #   The structure of the Pasta output files is
        #       Jobname.marker001.Unaligned_name.aln
        #       Jobname.tre
        aln_out = ''.join([
            pasta_out,
            os.path.sep,
            pasta_job,
            '.marker001.',
            self.input_seq.name.split('/')[-1].replace('.fasta', ''),
            '.aln'])
        tree_out = ''.join([
            pasta_out,
            os.path.sep,
            pasta_job,
            '.tre'])
        return (out, err, aln_out, tree_out)
