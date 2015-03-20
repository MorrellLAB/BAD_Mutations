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
            s = SeqIO.read(f, 'fasta')
        except ValueError:
            log.error('Input file ' + f + ' has more than one record. This script only accepts single-record FASTA files.')
            return False
        return True

#   A function to parse the input substitutions file
def parse_subs(f, log):
    #   Does the file exist?
    if not file_funcs.file_exists(f, log):
        log.error('File ' + f + ' does not exist.')
        return False
    else:
        #   Begin parsing it
        subs_data = []
        with open(f, 'r') as subfile:
            for index, line in enumerate(subfile):
                tmp = line.strip().split('\t')
                #   Check the fields. The first one should be numeric
                try:
                    pos = int(tmp[0])
                except ValueError:
                    log.error('Line ' + str(index + 1) + ' of input file ' + f + ': First field is not an integer.')
                    exit(1)
                #   If we can sucessfully cast it to integer, then we continue
                #   If there is only one item in the list, then the SNP ID was not supplied
                #   we drop in the empty string
                if len(tmp) == 1:
                    snpid = ''
                    log.warning('Variant on line ' + str(index + 1) + ' of input file ' + f + ' does not have an ID. Using the empty string (\'\') as an ID.')
                else:
                    snpid = tmp[1]
                #   Return these as a tuple
                subs_data.append((pos, snpid))
    log.info('Input file ' + f + ' contains ' + str(index+1) + ' positions to predict.')
    return subs_data
