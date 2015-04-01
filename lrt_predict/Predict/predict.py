#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import tempfile
import subprocess

#   Import our helper scripts here
from ..General import parse_input
from ..General import set_verbosity

class LRTPredict:
    def __init__(self, unaligned_sequences, query_sequence, substitutions, verbose):
        self.mainlog = set_verbosity('LRT Prediction', verbose)
        self.input_seq = unaligned_sequences
        self.query = query_sequence
        self.substitutions = parse_input.parse_subs(substitutions, self.mainlog)
        return

    #   A function to prepare the prank input file
    def add_query_to_seqlist(self):
        pass

    #   A function to call the prank alignment
    def prank_align(self):
        pass

    #   A function to extract the amino acid states from the prank alignment
    def get_aa_states(self):
        pass

    #   A function to get the LRT predictions
    def predict_codons(self):
        pass
