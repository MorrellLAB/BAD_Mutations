#!/usr/bin/env python

#   Script to perform LRT from Chun and Fay (2009) to predict deleterious SNPs
#   in plants.
#   Requires a user name and password for the JGI phyotzome portal (Free)
#   Thomas Kono
#       March 4, 2015
#       Saint Paul, MN

#   Dependencies:
#       1) requests http://www.python-requests.org
#       2) argparse https://code.google.com/p/argparse/
#           (Not required for Python 3 and Python 2 >= 2.7)
#       3) Biopython http://biopython.org/
#       4) BLAST+ executables from NCBI


#   To handle passwords
import getpass
#   To cd
import os
#   to check arguments
import sys
#   for verbosity messages
import logging

#   Import the dependency checking script
from lrt_predict.General import check_modules
#   Import the main fetching script
import lrt_predict.Fetch.phytozome as phytozome
#   Import the verbosity script
from lrt_predict.General import set_verbosity
#   Import our argument parsing script
from lrt_predict.General import parse_args
#   Import the BLAST search script
from lrt_predict.Blast import blast_search

#   A function to do the fetching
def fetch(arg, log):
    #   Create a new Phytozome instance that will handle our work with
    #   the JGI Genomes Portal.
    log.info('Creating a new instance to fetch data')
    #   We give it username, password, base directory, whether or not we have to log in and 
    p = phytozome.Phytozome(arg['user'], arg['password'], arg['base'], arg['convert_only'], arg['verbose'])
    if arg['convert_only']:
        log.info('Only converting files.')
        p.convert()
    elif arg['fetch_only']:
        log.info('Only downloading files.')
        p.get_xml_urls()
        p.fetch_cds()
    else:
        log.info('Downloading and converting files.')
        p.get_xml_urls()
        p.fetch_cds()
        p.convert()
    return


#   A function to do the BLAST searching
def blast(arg, log):
    log.info('Creating a new instance to BLAST.')
    b = blast_search.BlastSearch(arg['base'], arg['fasta'], arg['evalue'], arg['verbose'])
    b.blast_all()
    return b.homologues


#   A function to do the predicting
#   Main function
def main():
    #   Parse the arguments
    #   First, a check to see if any arguments were sent at all
    #   If not, then print the usage and exit
    if not sys.argv[1:]:
        parse_args.usage()
        exit(1)
    arguments = parse_args.parse_args()
    #   Pull out the verbosity switch right away
    verbose = set_verbosity.verbosity('LRT_Predict', arguments.verbose)
    arguments_valid, msg = parse_args.validate_args(arguments)
    #   If we got a return value that isn't False, then our arguments are good
    if arguments_valid:
        verbose.debug(arguments_valid['action'] + ' subcommand was invoked')
        #   Which command was invoked?
        if arguments_valid['action'] == 'fetch':
            #   Send it to the fetch command
            fetch(arguments_valid, verbose)
        elif arguments_valid['action'] == 'predict':
            homologues = blast(arguments_valid, verbose)
            print homologues
    else:
        verbose.error(msg)
    return


#   Do the work here
main()
