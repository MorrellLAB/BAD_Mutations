#!/usr/bin/env python
#   Helper script for Phytozome.py
#   Contains argument parsing code

try:
    import argparse
except ImportError:
    print 'Error! You need to have the argparse module installed.'
    exit(1)

import os


#   Import the helper script to validate arguments
import check_args

#   A function to actually parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description = 'LRT for deleterious SNP prediction in plants.',
        add_help=True)
    #   Define a sub-parser to handle the different actions we can perform
    #   we can either 'fetch' or 'predict'
    subparser = parser.add_subparsers(
        dest='action',
        title='Available actions',
        help='Sub-command help')
    #   Create a parser for 'fetch'
    fetch_args = subparser.add_parser(
        'fetch',
        help='Fetch CDS files from Phytozome and Ensembl')
    #   And give it some arguments
    fetch_args.add_argument(
        '--user',
        required=False,
        default=None,
        help='Username for jgi.doe.gov (For fetching from Phytozome)')
    fetch_args.add_argument(
        '--password',
        required=False,
        help='Password for jgi.doe.gov. If you are not comfortable supplying\
        this on the command-line in text, you can enter it on the prompt.',
        default=None)
    fetch_args.add_argument(
        '--base',
        '-b',
        required=False,
        help='Base directory for species databses. Defaults to .',
        default=os.getcwd())
    #   Create a new mutually exclusive group for deciding if we want to only
    #   fetch, or if we want to convert
    actions = fetch_args.add_mutually_exclusive_group(required=False)
    actions.add_argument(
        '--fetch-only',
        required=False,
        action='store_true',
        default=False,
        help='Do not convert CDS files to BLAST databases, just fetch.'
        )
    actions.add_argument(
        '--convert-only',
        required=False,
        action='store_true',
        default=False,
        help='Do not fetch new CDS from databases, just convert to BLAST db.'
        )
    #   Create a parser for 'predict'
    predict_args = subparser.add_parser(
        'predict',
        help='Run the LRT prediction pipeline')
    args = parser.parse_args()
    return args


#   Here we validate the arguments
def validate_args(args):
    #   We will convert it into a dictionary to access the data
    dict_args = vars(args)
    #   Check the base argument. If it starts with something other than a /
    #   then it is a relative path, and we should fix it
    if not dict_args['base'].startswith('/'):
        #   Add the cwd onto it, since the script fails otherwise
        dict_args['base'] = os.path.join(os.getcwd(), dict_args['base'])
    #   Then check the action
    #   If we are fetching, we have to check the username and base
    #   argparse should have checked for missing arguments by now
    #   If the arguments do not check out, return a message
    if dict_args['action'] == 'fetch':
        if dict_args['user']:
            if not check_args.valid_email(dict_args['user']):
                return (False, 'Username is not a valid e-mail address.')
        elif not check_args.valid_dir(dict_args['base']):
            return (False, 'Base directory is not readable/writable, or does not exist.')
        else:
            return (dict_args, None)
    #   The other subcommand isn't implemented yet
    else:
        return (dict_args, None)


#   This is just a simple function that shows when the user does not supply
#   any arguments. This is an issue with Python 2.*, and has been "fixed" in
#   Python 3+.
def usage():
    print '''Usage: LRT_Predict.py <subcommand> <options>

where <subcommand> is one of 'fetch' or 'predict.' This script will download
the necessary data to perform the likelihood ratio test (LRT) of deleterious
SNP prediction in plants, as described in Chun and Fay (2009) in Genome
Research. 

'fetch' will download gzipped CDS FASTA files from Phytozome, unzip them, and
convert them into BLAST databases. It requires a username and password for the
JGI Genomes Portal, which is free. Check with Phytozome for their data release
and usage policies. Use 'LRT_Predict.py fetch -h' for more infomration.

'predict' will run the LRT with a given query sequence and a list of affected
codons. This feature has not yet been implemented.

Dependencies:
    tblastx (NCBI BLAST executables)
    requests (Python HTTP requests module)
    prank-msa (Phylogeny-aware sequence alignment)'''
    return
