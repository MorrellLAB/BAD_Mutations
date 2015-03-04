#!/usr/bin/env python
#   Helper script for Phytozome.py
#   Contains argument parsing code

import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description = "Script to download CDS FASTA files from Phytozome.")
    parser.add_argument(
        '--user', 
        required=True, 
        help='User name for jgi.doe.gov')
    parser.add_argument(
        '--password', 
        required=True, 
        help='Password for jgi.doe.gov')
    args = parser.parse_args()
    return args
