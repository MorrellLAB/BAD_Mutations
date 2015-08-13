#!/usr/bin/env python
"""
Script to perform LRT from Chun and Fay (2009) to predict deleterious SNPs in
plants. Requires a user name and password for the JGI phyotzome portal (Free)
Author: Thomas Kono
        March 4, 2015
        Saint Paul, MN
"""
#   Dependencies:
#       1) requests http://www.python-requests.org
#       2) argparse https://code.google.com/p/argparse/
#           (Not required for Python 3 and Python 2 >= 2.7)
#       3) Biopython http://biopython.org/
#       4) BLAST+ executables from NCBI

#   to check arguments
import sys
#   To handle file copy
import shutil
import os

#   Import the dependency checking script
import lrt_predict.General.check_modules as check_modules
#   Import the verbosity script
import lrt_predict.General.set_verbosity as set_verbosity
#   Import our argument parsing script
import lrt_predict.General.parse_args as parse_args


def setup(arg):
    """A function to write the configuration file."""
    #   Check for the required modules for setup
    setupdeps = check_modules.check_modules(setup=True)
    if setupdeps:
        check_modules.missing_mods(setupdeps)
        exit(1)
    #   Import the setup script
    import lrt_predict.Setup.setup_env as setup_env
    #   Start a new instance of the configuration class
    s_env = setup_env.SetupEnv(
        arg['base'],
        arg['deps_dir'],
        arg['target'],
        arg['evalue'],
        arg['codon'],
        arg['missing_threshold'],
        arg['config'],
        arg['loglevel'])
    s_env.write_config()
    return


def fetch(arg, log):
    """A function to download the appropriate files from Ensembl Plants and
    Phytozome. Will also convert them to BLAST databases."""
    fetchdeps = check_modules.check_modules(fetch=True)
    if fetchdeps:
        check_modules.missing_mods(fetchdeps)
        exit(1)
    #   Next we check for the presence of the utilities we need to
    #   do the fetching
    missing_reqs = check_modules.missing_executables(
        ['bash',
         'makeblastdb',
         'gzip',
         'sum'])
    if missing_reqs:
        log.error(
            'Some required executables were not found on your system: ' +
            '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   Import the main fetching script
    #   We do it here, since we only want to import this if we are fetching
    import lrt_predict.Fetch.phytozome as phytozome
    import lrt_predict.Fetch.ensembl as ensembl
    #   Create a new Phytozome instance that will handle our work with
    #   the JGI Genomes Portal.
    log.debug('Creating a new Phytozome instance to fetch data.')
    phy = phytozome.Phytozome(
        arg['user'],
        arg['password'],
        arg['base'],
        arg['convert_only'],
        arg['loglevel'])
    log.debug('Creating a new Ensembl instance to fetch data.')
    ens = ensembl.EnsemblPlants(
        arg['base'],
        arg['convert_only'],
        arg['loglevel'])

    if arg['convert_only']:
        log.debug('Only converting files.')
        ens.convert()
        phy.convert()
    elif arg['fetch_only']:
        log.debug('Only downloading files.')
        log.info('Fetching from Ensembl Plants...')
        ens.get_ftp_urls()
        ens.download_files()
        log.info('Fetching from Phytozome...')
        phy.get_xml_urls()
        phy.fetch_cds()
    else:
        log.debug('Downloading and converting Ensembl Plants files...')
        ens.get_ftp_urls()
        ens.download_files()
        log.debug('Downloading and converting files.')
        phy.get_xml_urls()
        phy.fetch_cds()
        ens.convert()
        phy.convert()
    return


def blast(arg, log):
    """A function to search the databses with BLAST and collect the
    homologous sequences from them."""
    blastdeps = check_modules.check_modules(predict=True)
    #   Check the module dependencies for the BLAST function
    if blastdeps:
        check_modules.missing_mods(blastdeps)
        exit(1)
    missing_reqs = check_modules.missing_executables(['bash', 'tblastx'])
    #   And then check the executable dependencies
    if missing_reqs:
        log.error(
            'Some required executables were not found on your system: ' +
            '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   If all that checks out, import the BLAST class script
    import lrt_predict.Blast.blast_search as blast_search
    log.info('Creating a new instance to BLAST.')
    b_search = blast_search.BlastSearch(
        arg['base'],
        arg['fasta'],
        arg['evalue'],
        arg['loglevel'])
    b_search.blast_all()
    #   hom contains the file object that has the unaligned sequence in it.
    hom = b_search.get_hit_seqs()
    return hom


def align(arg, unaligned, log):
    """A function to align the homologous sequences with pasta, and return
    the aligned sequences and the phylogenetic tree."""
    aligndeps = check_modules.check_modules(predict=True)
    if aligndeps:
        check_modules.missing_mods(aligndeps)
        exit(1)
    #   Check for the required executables
    missing_reqs = check_modules.missing_executables(['bash', 'run_pasta.py'])
    if missing_reqs:
        log.error(
            'Some required executables were not found on your system: ' +
            '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   Then we import the necessary modules
    import lrt_predict.Predict.align as aligner
    log.info('Creating a new instance of PastaAlign.')
    aln = aligner.PastaAlign(
        unaligned,
        arg['fasta'],
        arg['loglevel'])
    #   Then align them
    stdout, stderr, alignment, tree = aln.pasta_align()
    log.debug('stdout: \n' + stdout)
    log.debug('stderr: \n' + stderr)
    return (alignment, tree)


def predict(arg, nuc, tree, log):
    """A function to run the HYPHY codon prediction model on each column of
    the alignment and return a score for each one."""
    predictdeps = check_modules.check_modules(predict=True)
    if predictdeps:
        check_modules.missing_mods(predictdeps)
        exit(1)
    #   Check for the required executables
    missing_reqs = check_modules.missing_executables(['bash', 'HYPHYSP'])
    if missing_reqs:
        log.error(
            'Some required executables were not found on your system: ' +
            '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   import the predict script
    import lrt_predict.Predict.predict as predictor
    #   Create a new instance of class LRTPredict
    lrt = predictor.LRTPredict(
        nuc,
        tree,
        arg['fasta'],
        arg['substitutions'],
        arg['loglevel'])
    position = lrt.get_query_position()
    return position


def main():
    """The main function."""
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
    #   Import the config handling script
    import lrt_predict.Setup.parse_config as config
    #   vars() will convert from the weird Namespace type to a dictionary
    arguments = vars(parse_args.parse_args())
    #   Pull out the verbosity switch right away
    loglevel = set_verbosity.verbosity('LRT_Predict', arguments['loglevel'])
    #   If the config variable was passed
    if arguments['config']:
        #   We ask then if the setup command was not passed
        #   If the user wants to setup, then don't bother trying to validate
        #   the config.
        if arguments['action'] != 'setup':
            cfg = config.ConfigHandler(
                arguments['config'],
                arguments,
                arguments['loglevel'])
            if cfg.is_valid():
                cfg.read_vars()
                config_opts = cfg.merge_options()
            else:
                loglevel.error('Config file is not valid!')
                exit(1)
        #   Else, just set it to the options that were passed
        else:
            config_opts = arguments
    else:
        config_opts = arguments
    arguments_valid, msg = parse_args.validate_args(config_opts, loglevel)
    #   If we got a return value that isn't False, then our arguments are good
    if arguments_valid:
        loglevel.debug(arguments_valid['action'] + ' subcommand was invoked')
        #   Which command was invoked?
        if arguments_valid['action'] == 'setup':
            setup(arguments_valid)
        elif arguments_valid['action'] == 'fetch':
            fetch(arguments_valid, loglevel)
        elif arguments_valid['action'] == 'predict':
            #   We will return the filename that contains the unaligned
            #   sequences, as we will use these as inputs for prank
            unaligned_seqs = blast(arguments_valid, loglevel)
            #   Then add the query sequence and align them
            alignment, tree_file = align(arguments_valid, unaligned_seqs, loglevel)
            loglevel.info('Nucleotide alignment in ' + alignment)
            loglevel.info('Tree in ' + tree_file)
            new_nuc = '/Users/tomkono/DataDisk/tmp/JCF_Barley_CSV_Testing/Pasta_Tests/' + file_funcs.local_name(alignment)
            new_tree = '/Users/tomkono/DataDisk/tmp/JCF_Barley_CSV_Testing/Pasta_Tests/' + file_funcs.local_name(tree_file)
            open(new_nuc, 'w').close()
            open(new_tree, 'w').close()
            shutil.copy2(alignment, new_nuc)
            shutil.copy2(tree_file, new_tree)
            #predict(arguments_valid, nuc_file, tree_file, loglevel)
    else:
        loglevel.error(msg)
    return

#   Do the work here
main()
