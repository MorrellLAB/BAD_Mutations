#!/usr/bin/env python

#   A script that performs the alignment and LRT prediction

#   Import standard library modules here
import subprocess
import os

#   Import our helper scripts here
from ..General import parse_input
from ..General import set_verbosity
from ..General import check_modules

class LRTPredict:
    def __init__(self, alignment, substitutions, verbose):
        self.mainlog = set_verbosity.verbosity('LRT_Prediction', verbose)
        self.msa = alignment
        self.substitutions = parse_input.parse_subs(substitutions, self.mainlog)
        return
    #   A function to extract the amino acid states from the prank alignment
    def get_aa_states(self):
        pass

    #   A function to get the LRT predictions
    def predict_codons(self):
        pass
