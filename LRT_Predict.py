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


#   Define a function for fetching
def fetch(arg, log):
    fetchdeps = check_modules.check_modules(fetch=True)
    if fetchdeps:
        check_modules.missing_mods(fetchdeps)
        exit(1)
    #   Next we check for the presence of the utilities we need to 
    #   do the fetching
    missing_reqs = check_modules.missing_executables(['bash', 'makeblastdb', 'gzip', 'sum'])
    if missing_reqs:
        log.error('Some required executables were not found on your system: ' + '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   Import the main fetching script
    #   We do it here, since we only want to import this if we are fetching
    import lrt_predict.Fetch.phytozome as phytozome
    import lrt_predict.Fetch.ensembl as ensembl
    #   Create a new Phytozome instance that will handle our work with
    #   the JGI Genomes Portal.
    log.info('Creating a new Phytozome instance to fetch data.')
    #   We give it username, password, base directory, whether or not we have to log in and 
    p = phytozome.Phytozome(arg['user'], arg['password'], arg['base'], arg['convert_only'], arg['loglevel'])
    log.info('Creating a new Ensembl instance to fetch data.')
    ens = ensembl.EnsemblPlants(arg['base'], arg['convert_only'], arg['loglevel'])
    if arg['convert_only']:
        log.info('Only converting files.')
        ens.convert()
        p.convert()
    elif arg['fetch_only']:
        log.info('Only downloading files.')
        log.info('Fetching from Ensembl Plants...')
        ens.get_ftp_urls()
        ens.download_files()
        log.info('Fetching from Phytozome...')
        p.get_xml_urls()
        p.fetch_cds()
    else:
        log.info('Downloading and converting Ensembl Plants files...')
        ens.get_ftp_urls()
        ens.download_files()
        log.info('Downloading and converting files.')
        p.get_xml_urls()
        p.fetch_cds()
        ens.convert()
        p.convert()
    return


#   Define a function for BLASTing
def blast(arg, log):
    blastdeps = check_modules.check_modules(predict=True)
    #   Check the module dependencies for the BLAST function
    if blastdeps:
        check_modules.missing_mods(blastdeps)
        exit(1)
    missing_reqs = check_modules.missing_executables(['bash', 'tblastx'])
    #   And then check the executable dependencies
    if missing_reqs:
        log.error('Some required executables were not found on your system: ' + '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   If all that checks out, import the BLAST class script
    from lrt_predict.Blast import blast_search
    log.info('Creating a new instance to BLAST.')
    b = blast_search.BlastSearch(arg['base'], arg['fasta'], arg['evalue'], arg['loglevel'])
    b.blast_all()
    #   hom contains the file object that has the unaligned sequence in it.
    hom = b.get_hit_seqs()
    return hom


#   Define a function for aligning
def align(arg, unaligned, log):
    aligndeps = check_modules.check_modules(predict=True)
    if aligndeps:
        check_modules.mossing_mods(aligndeps)
        exit(1)
    #   Check for the required executables
    missing_reqs = check_modules.missing_executables(['bash', 'prank'])
    if missing_reqs:
        log.error('Some required executables were not found on your system: ' + '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   Then we import the necessary modules
    from lrt_predict.Predict import align
    log.info('Creating a new instance of PrankAlign.')
    a = align.PrankAlign(unaligned, arg['fasta'], arg['loglevel'])
    #   Then align them
    a.add_query_to_seqlist()
    stdout, stderr, outfile = a.prank_align()
    log.debug('stdout: \n' + stdout)
    log.debug('stderr: \n' + stderr)
    return outfile

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
    loglevel = set_verbosity.verbosity('LRT_Predict', arguments.loglevel)
    arguments_valid, msg = parse_args.validate_args(arguments, loglevel)
    #   If we got a return value that isn't False, then our arguments are good
    if arguments_valid:
        loglevel.debug(arguments_valid['action'] + ' subcommand was invoked')
        #   Which command was invoked?
        if arguments_valid['action'] == 'fetch':
            fetch(arguments_valid, loglevel)
        elif arguments_valid['action'] == 'predict':
            #   We will return the filename that contains the unaligned
            #   sequences, as we will use these as inputs for prank
            unaligned_seqs = blast(arguments_valid, loglevel)
            #   Then add the query sequence and align them
            alignment = align(arguments_valid, unaligned_seqs, loglevel)
            #   prank creates files with a certain prefix
            alignment_nuc_file = alignment.name + '.best.nuc.fas'
            alignment_pep_file = alignment.name + '.best.pep.fas'
            alignment_tree_file = alignment.name + '.best.dnd'
            loglevel.info('Nucleotide alignment in ' + alignment_nuc_file)
            loglevel.info('Protein alignment in ' + alignment_pep_file)
            loglevel.info('Tree in ' + alignment_tree_file)
    else:
        loglevel.error(msg)
    return


#   Do the work here
main()
