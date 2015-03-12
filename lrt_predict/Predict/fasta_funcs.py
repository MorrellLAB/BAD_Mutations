#!/usr/bin/env python

#   A script to parse the input fasta file

try:
    from Bio import SeqIO
except ImportError:
    print 'You need the Biopython libraries installed!'
    exit(1)

#   A function that parses out the relevant information from the FASTA file
#   This function expects a SeqRecord object, from Biopython
def get_fasta_info(sequence_record):
    #   Get the name of the sequence
    name = sequence_record.name
    #   And the sequence
    seq = str(sequence_record.seq)
    #   Get the qualifiers
    #   The description attribute contains anything after the > in the
    #   record identifier. We take the last elements of it
    qualifiers = sequence_record.description.split(' ')[1:]
    #   Then, traverse the list of qualifiers and ask if there are codons
    #   or SNP IDs, and populate as necessary
    codon_list = []
    id_list = []
    for q in qualifiers:
        #   Check for codons
        if q.startswith('codons'):
            #   Separate the keyword from the list of data
            keyword, c = q.split('=')
            #   Then convert them to integers
            codon_list = [int(x) for x in c.split(',')]
        #   Check for  SNP IDs
        if q.startswith('IDs'):
            #   separate it again
            keyword, i = q.split('=')
            id_list = i.split(',')
        #   Any unrecognized keywords are ignored.
    #   Build the dict to return the data with
    fasta_info = {
        'name': name,
        'seq': seq,
        'codons': codon_list,
        'ids': id_list}
    return fasta_info


#   A function that validates the fasta information. This is just checking that
#   the codons list isn't empty and that if the SNP ID list is present that it is
#   the same length as the codon list
def validate_fasta_info(fasta_info):
    #   Build a standard prefix for the error messages
    err_prefix = 'Error in defline of ' + fasta_info['name']
    #   first, we check if the codon list is a nonzero length
    #   Empty lists evaluate to Boolean False
    if not fasta_info['codons']:
        #   Build a nice error message
        msg = err_prefix + ': codon list is empty!'
        valid = False
    #   Check the lengths of SNP ID list and codon list
    elif len(fasta_info['codons']) != len(fasta_info['ids']):
        #   Again, build a nice error message
        msg = err_prefix + ': codon list and ID list are different lengths!'
        valid = False
    elif len(fasta_info['ids']) != len(set(fasta_info['ids'])):
        #   Make sure there aren't any duplicates!
        msg = err_prefix + ': duplicate IDs!'
        valid = False
    elif len(fasta_info['codons']) != len(set(fasta_info['codons'])):
        #   Same for codon list
        msg = err_prefix + ': duplicate codon positions!'
        valid = False
    else:
        #   It all checks out
        msg = None
        valid = True
    return (valid, msg)
