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
#   Import the verbosity script
from lrt_predict.General import set_verbosity
#   Import our argument parsing script
from lrt_predict.General import parse_args


#   A function to do the predicting
#   Main function
def main():
    #   The very first thing we do is do a base check to make sure that we can
    #   parse arguments
    dep = check_modules.check_modules()
    if dep:
        check_modules.missing_mods(dep)
        exit(1)
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
            #   Next, we check the modules that are required by each subcommand
            fetchdeps = check_modules.check_modules(fetch=True)
            if fetchdeps:
                check_modules.missing_mods(fetchdeps)
                exit(1)
            #   Import the main fetching script
            #   We do it here, since we only want to import this if we are fetching
            import lrt_predict.Fetch.phytozome as phytozome
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
            #   Send it to the fetch command
            fetch(arguments_valid, verbose)
        elif arguments_valid['action'] == 'predict':
            #   Next, we check the modules that are required by each subcommand
            predictdeps = check_modules.check_modules(predict=True)
            if predictdeps:
                check_modules.missing_mods(predictdeps)
                exit(1)
            #   Import the BLAST search script
            #   Again, import here because we only want call these functions if we are running BLAST
            from lrt_predict.Blast import blast_search
            #   A function to do the BLAST searching
            def blast(arg, log):
                log.info('Creating a new instance to BLAST.')
                b = blast_search.BlastSearch(arg['base'], arg['fasta'], arg['evalue'], arg['verbose'])
                b.blast_all()
                hom = b.get_hit_seqs()
                hom.seek(0)
                print hom.read()
                hom.close()
                return
            blast(arguments_valid, verbose)
    else:
        verbose.error(msg)
    return


#   Do the work here
main()
