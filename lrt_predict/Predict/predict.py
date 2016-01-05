#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import subprocess
import os
import tempfile

#   Import Biopython library
from Bio import AlignIO
from Bio import SeqIO

#   Import our helper scripts here
from ..General import parse_input
from ..General import set_verbosity
from ..General import check_modules


class LRTPredict(object):
    def __init__(
            self,
            hyphy_path,
            nuc_aln,
            treefile,
            query,
            substitutions,
            verbose):
        self.mainlog = set_verbosity.verbosity('LRT_Prediction', verbose)
        self.hyphy_path = check_modules.check_executable(hyphy_path)
        self.nmsa = AlignIO.read(open(nuc_aln, 'r'), 'fasta')
        self.phylogenetic = treefile
        self.query = SeqIO.read(query, 'fasta')
        self.substitutions = parse_input.parse_subs(substitutions, self.mainlog)
        self.query_pos = 0
        self.aligned_pos = None
        self.hyphy_input = None
        self.hyphy_output = None
        return

    def get_query_position(self):
        """Get the index of the query sequence in the Pasta alignment."""
        qname = self.query.id.replace(':', '_')
        #   And step through it, saving the position of the query
        for index, sequence in enumerate(self.nmsa):
            if sequence.id == qname:
                self.mainlog.debug(qname + ' is at position ' + str(index))
                self.query_pos = index
                break
        return

    def get_aligned_positions(self):
        """Get the position of the query codons in the alignment. Uses the
        query sequence, and counts gap characters to calculate the aligned
        positions to be computed."""
        #   Parse the alignment object and get the list of aligned positions
        #   that need to be predicted.
        self.mainlog.debug(
            'Searching for positions: ' +
            ', '.join([str(i) for i in self.substitutions]))
        self.aligned_pos = []
        real_position = 0
        self.mainlog.debug(self.nmsa[self.query_pos].seq)
        for index, column in enumerate(self.nmsa[self.query_pos].seq):
            if column == '-':
                continue
            else:
                real_position += 1
            if real_position % 3 == 0:
                if real_position / 3 in self.substitutions:
                    self.aligned_pos.append(index+1)
        self.mainlog.debug(
            'Aligned Pos: ' + ', '.join([str(i) for i in self.aligned_pos]))
        return

    def write_aligned_subs(self):
        """Write the aligned positions into a temporary file."""
        subsfile = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_HYPHY_Subs_',
            suffix='.txt',
            delete=False
            )
        subsfile.write('\n'.join([str(i) for i in self.aligned_pos]))
        return subsfile

    def prepare_hyphy_inputs(self):
        """Prepare the input files for the HYPHY prediction script. Writes the
        paths of the MSA, tree, and substitutions file into a plain text file.
        Also builds the name of the temporary output file to hold the
        predictions."""
        #   Write the substitutions file
        alignedsubs = self.write_aligned_subs()
        infile = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_HYHPY_In_',
            suffix='.txt'
            )
        #   We write the paths of the MSA, the tree, the positions, and the
        #   query name into the input file. But we need to put the full paths
        #   into the file.
        infile.write(os.path.abspath(self.nmsa) + '\n')
        infile.write(os.path.abspath(self.phylogenetic) + '\n')
        infile.write(alignedsubs.name + '\n')
        infile.write(self.query.id.replace(':', '_'))
        infile.flush()
        #   And then we create a name for the HYPHY output file
        outfile = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_HYPHY_Out_',
            suffix='.txt'
            )
        self.hyphy_input = infile
        self.hyphy_output = outfile
        return

    def predict_codons(self):
        """Run the HYPHY script to predict the codons."""
        #   Get the base directory of the LRT package
        lrt_path = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
        #   Then build the path to the hyphy script
        hyphy_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'Prediction.sh')
        #   And the actual HyPhy code that JCF wrote
        prediction_script = os.path.join(
            lrt_path,
            'Shell_Scripts',
            'LRT.hyphy')
        #   Build the command for predictig
        cmd = [
            'bash',
            hyphy_script,
            self.hyphy_path,
            prediction_script,
            self.hyphy_input.name,
            self.hyphy_output.name
            ]
        self.mainlog.debug(' '.join(cmd))
        #   Then run the command
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        self.mainlog.debug('stdout:\n' + out)
        self.mainlog.debug('stderr:\n' + err)
        #   Return the output file
        return self.hyphy_output
