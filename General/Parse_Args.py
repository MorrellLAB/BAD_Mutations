#!/usr/bin/env python
#   Helper script for Phytozome.py
#   Contains argument parsing code

try:
    import argparse
except ImportError:
    print 'Error! You need to have the argparse module installed.'
    exit(1)

def parse_args():
    parser = argparse.ArgumentParser(
        description = 'LRT for deleterious SNP prediction in plants.')
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
        required=True,
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
        default='.')
    #   Create a parser for 'predict'
    predict_args = subparser.add_parser(
        'predict',
        help='Run the LRT prediction pipeline')
    args = parser.parse_args()
    return args
