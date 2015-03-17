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

#   Import the main fetching script
import lrt_predict.Fetch.phytozome as phytozome
#   Import the verbosity script
from lrt_predict.General import set_verbosity
#   Import our argument parsing script
from lrt_predict.General import parse_args

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
            #   Create a new Phytozome instance that will handle our work with
            #   the JGI Genomes Portal.
            verbose.debug('Creating a new instance to fetch data')
            p = phytozome.Phytozome(
                arguments_valid['user'],
                arguments_valid['password'],
                arguments_valid['base'],
                arguments_valid['convert_only'],
                arguments_valid['verbose'])
            if arguments_valid['convert_only']:
                verbose.info('Only converting files.')
                p.convert()
            else:
                verbose.info('Downloading and converting files.')
                p.get_xml_urls()
                p.fetch_cds()
                p.convert()
        elif arguments_valid['action'] == 'predict':
            pass
    else:
        verbose.error(msg)
    return


#   Do the work here
main()
