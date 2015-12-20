#!/usr/bin/env python

#   A script that contains functions to parse user input into the prediction
#   pipeline.

#   Import standard library modules here
import re

#   We need to handle sequence records, alignments, and phylogenies
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.NewickIO import NewickError

#   Our helper scripts
import set_verbosity
import file_funcs


def valid_tree(f, log):
    """Check that the phylogenetic tree is valid. This only checks the tree
    structure and doesn't check any of the branch lengths or names."""
    if not file_funcs.file_exists(f, log):
        log.error('File ' + f + ' does not exist')
        return False
    else:
        #   Phylo.read() raises a NewickError when the tree is not valid
        try:
            p = Phylo.read(f, 'newick')
        except NewickError:
            log.error(
                'Input file ' + \
                f + \
                ' is not a valid Newick tree file!')
            return False
        return True


def valid_msa(f, log):
    """Check if the MSA is a valid sequence alignment or not. All sequences
    should be the same length, and should be in FASTA format."""
    if not file_funcs.file_exists(f, log):
        log.error('File ' + f + ' does not exist.')
        return False
    else:
        #   AlignIO.read() raises a ValueError if the alignment is not in the
        #   right format, or if not all the sequences are the same length
        try:
            a = AlignIO.read(f, 'fasta')
        except ValueError:
            log.error(
                'Input file ' + \
                f + \
                ' is not a valid FASTA alignment!' + \
                ' Check the length of each sequence.')
            return False
        return True


def valid_fasta(f, log):
    """Check if the FASTA supplied is valid."""
    #   Does the file exist?
    if not file_funcs.file_exists(f, log):
        log.error('File ' + f + ' does not exist.')
        return False
    else:
        #   Start checking it
        try:
            s = SeqIO.read(f, 'fasta')
        except ValueError:
            log.error(
                'Input file ' + \
                f + \
                ' has more than one record. '+ \
                'This script only accepts single-record FASTA files.')
            return False
        return True


def parse_subs(f, log):
    """Parse the input substitutions file. Returns a list of integers."""
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
                    log.error(
                        'Line ' + \
                        str(index + 1) + \
                        ' of input file ' + f \
                        + ': First field is not an integer.')
                    exit(1)
                #   If we can sucessfully cast it to integer, then we continue
                #   If there is only one item in the list, then the SNP ID
                #   is abset. We drop in the empty string
                if len(tmp) == 1:
                    snpid = ''
                    log.warning(
                        'Variant on line ' + \
                        str(index + 1) + \
                        ' of input file ' + \
                        f + \
                        ' does not have an ID. ' + \
                        'Using the empty string (\'\') as an ID.')
                else:
                    snpid = tmp[1]
                #   Return these as a tuple
                subs_data.append(pos)
    log.info(
        'Input file ' + \
        f + \
        ' contains ' + \
        str(index+1) + \
        ' positions to predict.')
    return subs_data
