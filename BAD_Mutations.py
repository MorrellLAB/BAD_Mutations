#!/usr/bin/env python
"""Perform a likelihood ratio test from Chun and Fay (2009) for deleterious
SNP prediction in plants.

Authors: Thomas JY Kono and Paul J Hoffman
Paper: https://doi.org/10.1534/g3.118.200563
"""

#   to check arguments
import sys
#   To handle file copy
import shutil
import os
import pprint

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
        [arg['bash_path'],
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
    missing_reqs = check_modules.missing_executables(
        [
            arg['bash_path'],
            arg['tblastx_path']
        ])
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
        arg['target'],
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
    missing_reqs = check_modules.missing_executables(
        [
            arg['bash_path'],
            arg['pasta_path'],
            arg['clustalo_path'],
            arg['fasttree_path']
        ])
    if missing_reqs:
        log.error(
            'Some required executables were not found on your system: ' +
            '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   Then we import the necessary modules
    import lrt_predict.Predict.align as aligner
    log.info('Creating a new instance of PastaAlign.')
    aln = aligner.PastaAlign(
        arg['pasta_path'],
        arg['clustalo_path'],
        arg['fasttree_path'],
        unaligned,
        arg['fasta'],
        arg['loglevel'])
    #   Prepare the sequences for alignment:
    #       Check length is multiple of 3
    #       Translate to protein
    #       Remove STOP codons
    nseqs = aln.prepare_sequences()
    #   Then align them
    if nseqs > 2:
        stdout, stderr = aln.pasta_align()
        log.debug('stdout: \n' + stdout.decode('utf-8'))
        log.debug('stderr: \n' + stderr.decode('utf-8'))
    else:
        log.warning('Only two sequences; using clustal-omega for alignment.')
        stdout, stderr = aln.clustalo_align()
        log.debug('stdout: \n' + stdout.decode('utf-8'))
        log.debug('stderr: \n' + stderr.decode('utf-8'))
    #   Backtranslate the alignment
    aln.back_translate()
    #   Then sanitize the alignment and tree
    aln.sanitize_outputs()
    #   Then copy them over
    log.info('Nucleotide alignment in ' + aln.final_aln)
    log.info('Tree in ' + aln.tree_out)
    new_nuc = os.path.join(
        arg['output'],
        os.path.basename(
            arg['fasta'].replace(
                '.fasta',
                '_MSA.fasta')
            )
        )
    new_tree = os.path.join(
        arg['output'],
        os.path.basename(
            arg['fasta'].replace('.fasta', '.tree')
            )
        )
    open(new_nuc, 'w').close()
    open(new_tree, 'w').close()
    shutil.copy2(aln.final_aln, new_nuc)
    shutil.copy2(aln.tree_out, new_tree)
    log.info('MSA copied to ' + new_nuc)
    log.info('Tree copied to ' + new_tree)
    #   Cleanup the temporary file
    os.remove(aln.final_aln)
    return


def predict(arg, log):
    """A function to run the HYPHY codon prediction model on each column of
    the alignment and return a score for each one."""
    predictdeps = check_modules.check_modules(predict=True)
    if predictdeps:
        check_modules.missing_mods(predictdeps)
        exit(1)
    #   Check for the required executables
    missing_reqs = check_modules.missing_executables(
        [
            arg['bash_path'],
            arg['hyphy_path']
        ])
    if missing_reqs:
        log.error(
            'Some required executables were not found on your system: ' +
            '\n'.join(missing_reqs) + '\nPlease install them to continue.')
        exit(1)
    #   import the predict script
    import lrt_predict.Predict.predict as predictor
    #   Create a new instance of class LRTPredict
    lrt = predictor.LRTPredict(
        arg['hyphy_path'],
        arg['alignment'],
        arg['tree'],
        arg['fasta'],
        arg['substitutions'],
        arg['loglevel'])
    lrt.get_query_position()
    lrt.get_aligned_positions()
    lrt.write_aligned_subs()
    lrt.prepare_hyphy_inputs()
    outputfile = lrt.predict_codons()
    return outputfile


def compile_preds(arg, log):
    """A function to compile a directory full of HyPhy reports, and write a
    single report with all the necessary prediction information."""
    #   Import the HyPhyParser class source file
    import lrt_predict.Predict.hyphy_parser as hyphyparser
    #   Start a new hyphy parser object. We only need the directory containing
    #   the predictions.
    comp = hyphyparser.HyPhyParser(arg['pred_dir'], arg['loglevel'])
    #   Get all the HyPhy reports
    reports = comp.get_prediction_files()
    log.info('Found a total of ' + str(len(reports)) + ' reports.')
    #   Parse them into prediction data
    parsed_preds = [comp.parse_prediction(rep) for rep in reports]
    #   Read the long format substitutions file, keying on the same values as
    #   the prediction file (gene ID and postion). This file will also have
    #   the SNP ID in it.
    subs = arg['long_subs']
    alts = {}
    with open(subs, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                # Unpack the values.
                txid = tmp[0]
                aapos = tmp[1]
                alt_aa = tmp[2]
                snpid = tmp[3]
                # The HyPhy predictions are keyed on (txid, pos), but we have
                # to replace dot (.) with underscore because dots crash the
                # HyPhy sequence parser.
                key = (txid.replace('.', '_'), aapos)
                alts[key] = alt_aa
    #   Add P-values to the predictions
    logp_preds = []
    for genepred in parsed_preds:
        for snppred in genepred:
            a = alts.get((snppred[0], snppred[1]), 'NA')
            logp_preds.append(comp.add_regression(a, snppred))
        #   Then write them into the destination file
    comp.compile_predictions(logp_preds)
    return


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
    if 'config' in arguments:
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
        elif arguments_valid['action'] == 'align':
            #   We will return the filename that contains the unaligned
            #   sequences, as we will use these as inputs for pasta
            unaligned_seqs = blast(arguments_valid, loglevel)
            #   Then add the query sequence and align them
            align(
                arguments_valid,
                unaligned_seqs,
                loglevel)
        elif arguments_valid['action'] == 'predict':
            out = predict(arguments_valid, loglevel)
            #   copy the output file into the destination directory
            #   To build the output filename, we join the output directory
            #   with a new name based on the input filename
            out_fname = os.path.join(
                arguments_valid['output'],
                os.path.basename(
                    arguments_valid['fasta'].replace(
                        '.fasta',
                        '_Predictions.txt')
                    )
                )
            open(out_fname, 'w').close()
            shutil.copy2(out.name, out_fname)
            loglevel.info('Prediction in ' + out_fname)
        elif arguments_valid['action'] == 'compile':
            compile_preds(arguments_valid, loglevel)
            return
    else:
        loglevel.error(msg)
    return


main()
