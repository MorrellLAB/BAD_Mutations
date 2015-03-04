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
        description = 'Script to download CDS FASTA files from Phytozome.')
    parser.add_argument(
        '--user', 
        required=True, 
        help='User name for jgi.doe.gov')
    parser.add_argument(
        '--password', 
        required=True, 
        help='Password for jgi.doe.gov')
    parser.add_argument(
        '--base',
        '-b',
        help='Base directory for species databases. Defaults to .',
        default='.')
    args = parser.parse_args()
    return args
