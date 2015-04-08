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

class LRTPredict:
    def __init__(self, nuc_aln, pep_aln, treefile, query, substitutions, verbose):
        self.mainlog = set_verbosity.verbosity('LRT_Prediction', verbose)
        self.nmsa = nuc_aln
        self.pmsa = pep_aln
        self.phylogenetic = treefile
        self.query = query
        self.substitutions = parse_input.parse_subs(substitutions, self.mainlog)
        self.query_pos = 0
        return
    #   A function to get the index of the query sequence in the prank alignment
    def get_query_position(self):
        #   Get the name from the query sequence
        qseq = SeqIO.read(self.query, 'fasta')
        #   Read the alignment
        a = AlignIO.read(open(self.nmsa, 'r'), 'fasta')
        #   And step through it, saving the position of the query
        for index, sequence in enumerate(a):
            if sequence.id == qseq.id:
                self.mainlog.debug(qseq.id + ' is at position ' + str(index))
                self.query_pos = index
                break
        #   Return the alignment object, so we don't have to read it up again
        return a

    #   A function to extract the amino acid states from the prank alignment
    def get_codon_columns(self):
        pass

    #   A function to get the LRT predictions
    def predict_codons(self):
        pass
