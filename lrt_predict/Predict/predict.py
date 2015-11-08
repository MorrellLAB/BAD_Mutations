#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import subprocess
import os

#   Import Biopython library
from Bio import AlignIO
from Bio import SeqIO

#   Import our helper scripts here
from ..General import parse_input
from ..General import set_verbosity
from ..General import check_modules

class LRTPredict(object):
    def __init__(self, nuc_aln, treefile, query, substitutions, verbose):
        self.mainlog = set_verbosity.verbosity('LRT_Prediction', verbose)
        self.nmsa = nuc_aln
        self.phylogenetic = treefile
        self.query = query
        self.substitutions = parse_input.parse_subs(substitutions, self.mainlog)
        self.query_pos = 0
        return

    def get_query_position(self):
        """Get the index of the query sequence in the Pasta alignment."""
        #   Get the name from the query sequence
        qseq = SeqIO.read(self.query, 'fasta')
        #   Read the alignment
        a = AlignIO.read(open(self.nmsa, 'r'), 'fasta')
        self.qname = qseq.id
        #   And step through it, saving the position of the query
        for index, sequence in enumerate(a):
            if sequence.id == qseq.id:
                self.mainlog.debug(qseq.id + ' is at position ' + str(index))
                self.query_pos = index
                break
        #   Return the alignment object, so we don't have to read it up again
        return a

    def get_aligned_positions(self, alignment):
        """Get the position of the query codons in the alignment. Uses the
        query sequence, and counts gap characters to calculate the aligned
        positions to be computed."""
        #   Parse the alignment object and get the list of aligned positions
        #   that need to be predicted.
        self.aligned_pos = []
        real_position = 1
        for index, column in enumerate(alignment[self.query_pos:]):
            if column == '-':
                continue
            else:
                real_position += 1
            if real_position % 3 == 0:
                if real_position / 3 in self.subs:
                    self.aligned_position.append(real_position/3)

    def prepare_hyphy_inputs(self):
        """Prepare the input files for the HYPHY prediction script. Writes the
        paths of the MSA, tree, and substitutions file into a plain text file.
        Also builds the name of the temporary output file to hold the
        predictions."""
        infile = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_HYHPY_In_',
            suffix='.txt'
            )
        #   We write the paths of the MSA, the tree, the positions, and the
        #   query name into the input file.
        infile.write(self.nmsa.name + '\n')
        infile.write(self.phylogenetic + '\n')
        infile.write(self.subs + '\n')
        infile.write(self.qname)
        infile.flush()
        #   And then we create a name for the HYPHY output file
        outfile = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_HYPHY_Out_',
            suffix='.txt'
            )
        self.hyphy_input = infile.name
        self.hyphy_output = outfile.name
        return

    def predict_codons(self):
        """Run the HYPHY script to predict the codons."""
        #   Get the base directory of the LRT package
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the pasta script
        hyphy_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'Prediction.sh')
        hyphy_path = check_modules.check_executable('HYPHYMP')
        #   Build the command for predictig
        cmd = [
            'bash',
            hyphy_script,
            self.hyphy_input,
            self.hyphy_output
            ]
        self.mainlog.debug(' '.join(cmd))
        #   Then run the command
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()

