#!/usr/bin/env python

#   A script that contains functions to parse user input into the prediction
#   pipeline.

#   Import standard library modules here
import re

#   We need to handle sequence records
from Bio import SeqIO

#   Our helper scripts
import set_verbosity
from ..Fetch import file_funcs

#   A function to check the FASTA file input
def valid_fasta(f, log):
    #   Does the file exist?
    if not file_funcs.file_exists(f, log):
        log.error('File ' + f + ' does not exist.')
        return False
    else:
        #   Start checking it
        try:
            s = SeqIO.read(f, 'fatsta')
        except ValueError:
            log.error('Input file ' + f + ' has more than one record. This script only accepts single-record FASTA files.')
            return False
        return True
