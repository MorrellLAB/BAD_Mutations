#!/usr/bin/env python
"""Script to parse and validate arguments passed to BAD_Mutations.py"""

import argparse
import os
import getpass

#   Import the helper script to validate arguments
import lrt_predict.General.check_args as check_args
#   And the script to check the input files
import lrt_predict.General.parse_input as parse_input
#   And the script to operate on files
import lrt_predict.General.file_funcs as file_funcs
#   And the species lists
import lrt_predict.Fetch.ensembl_species as ensembl_species
import lrt_predict.Fetch.phytozome_species as phytozome_species

#   Create the list of allowable species by the combination of both
#   Phytozome and Ensembl. Currently this is restricted to angiosperms
SPECIES_LIST = ensembl_species.ensembl_fetch + phytozome_species.phyto_fetch


#   A function to actually parse the arguments
def parse_args():
    """Parse the arguments. We set up three subcommands here:
        setup:   Generate config file and download dependencies
        fetch:   Download CDS files from Phytozome and Ensembl
        predict: Align homologous sequences and predict deleterious subs

    See the manual or help message for detailed descriptions of all the
    arguments presented here."""
    parser = argparse.ArgumentParser(
        description='LRT for deleterious SNP prediction in plants.',
        add_help=True)
    #   Define a sub-parser to handle the different actions we can perform
    #   we can either 'fetch' or 'predict'
    subparser = parser.add_subparsers(
        dest='action',
        title='Available actions',
        help='Sub-command help')

    #   Create a parser for 'setup'
    setup_args = subparser.add_parser(
        'setup',
        help='Set up the runtime environment of BAD_Mutations')
    setup_args.add_argument(
        '--list-species',
        required=False,
        action='store_true',
        default=False,
        help='List the accepted names of query species (case sensitive).')
    setup_args.add_argument(
        '--config',
        '-c',
        required=False,
        help='Where to store the configuration file.',
        default=os.path.join(os.getcwd(), 'LRTPredict_Config.txt'))
    setup_args.add_argument(
        '--base',
        '-b',
        required=False,
        help='Base directory for species databases. Defaults to .',
        default=os.getcwd())
    setup_args.add_argument(
        '--deps-dir',
        '-d',
        required=False,
        help='Directory to house downloaded software. Defaults to .',
        default=os.getcwd())
    setup_args.add_argument(
        '--target',
        '-t',
        required=False,
        help=(
            'Which species are you predicting in (case sensitive)? Pass '
            '--list-species to see a full list of allowable species names.'
            ),
        default=None)
    setup_args.add_argument(
        '--evalue',
        '-e',
        required=False,
        default=0.05,
        type=float,
        help='E-value threshold for accepting sequences into the alignment.')
    setup_args.add_argument(
        '-m',
        '--min-seqs',
        required=False,
        type=int,
        default=10,
        help=(
            'Skip predictions for SNPs with fewer than this number of '
            'species represented in the multiple sequence alignment.'
            )
        )
    setup_args.add_argument(
        '--codon',
        required=False,
        action='store_const',
        const='codon',
        default='translate',
        help=('Use the codon alignment model for prank-msa, may give more '
              'accurate branch lengths but is much slower.'))

    #   Create a parser for 'fetch'
    fetch_args = subparser.add_parser(
        'fetch',
        help='Fetch CDS files from Phytozome and Ensembl')
    #   And give it some arguments
    fetch_args.add_argument(
        '--config',
        '-c',
        required=False,
        help='Use this configuration file.')
    fetch_args.add_argument(
        '--base',
        '-b',
        required=False,
        help='Base directory for species databses.')
    fetch_args.add_argument(
        '--user',
        '-u',
        required=False,
        default=None,
        help='Username for jgi.doe.gov (For fetching from Phytozome)')
    fetch_args.add_argument(
        '--password',
        '-p',
        required=False,
        help=(
            'Password for jgi.doe.gov. If you are not comfortable supplying '
            'this on the command-line in text, you can enter it on the prompt.'
            ),
        default=None)
    #   Create a new mutually exclusive group for deciding if we want to only
    #   fetch, or if we want to convert
    actions = fetch_args.add_mutually_exclusive_group(required=False)
    actions.add_argument(
        '--fetch-only',
        required=False,
        action='store_true',
        default=False,
        help='Do not convert CDS files to BLAST databases, just fetch.')
    actions.add_argument(
        '--convert-only',
        required=False,
        action='store_true',
        default=False,
        help='Do not fetch new CDS from databases, just convert to BLAST db.')

    #   Create a parser for 'align'
    align_args = subparser.add_parser(
        'align',
        help='Produce a multiple sequence alignment and phylogenetic tree.')
    align_args.add_argument(
        '--base',
        '-b',
        required=False,
        help='Base directory for species databses.')
    align_args.add_argument(
        '--config',
        '-c',
        required=False,
        help='Use this configuration file.')
    align_args.add_argument(
        '--evalue',
        '-e',
        required=False,
        default=0.05,
        type=float,
        help='E-value threshold for accepting sequences into the alignment.')
    align_args.add_argument(
        '--fasta',
        '-f',
        required=True,
        default=None,
        help='Path to the input FASTA file.')
    align_args.add_argument(
        '--output',
        '-o',
        required=False,
        default=os.getcwd(),
        help='Output directory.')

    #   Create a parser for 'predict'
    predict_args = subparser.add_parser(
        'predict',
        help='Run the LRT, generate HyPhy report.')
    predict_args.add_argument(
        '--config',
        '-c',
        required=False,
        help='Use this configuration file.')
    #   Give 'predict' some arguments
    predict_args.add_argument(
        '--fasta',
        '-f',
        required=True,
        default=None,
        help='Path to the input FASTA file.')
    predict_args.add_argument(
        '--alignment',
        '-a',
        required=True,
        default=None,
        help='Path to the input multiple sequence alignment.')
    predict_args.add_argument(
        '--tree',
        '-r',
        required=True,
        default=None,
        help='Path to the phylogenetic tree.')
    predict_args.add_argument(
        '--substitutions',
        '-s',
        required=True,
        default=None,
        help='Path to the input substitutions file.')
    predict_args.add_argument(
        '--output',
        '-o',
        required=False,
        default=os.getcwd(),
        help='Output directory.')

    #   Create a parser for 'compile'
    compile_args = subparser.add_parser(
        'compile',
        help=(
            'Compile a directory of HyPhy outputs into a single file with '
            'an evaluation of the impact for each SNP.'
            )
        )
    compile_args.add_argument(
        '--pred-dir',
        '-P',
        required=True,
        help='Directory where HyPhy outputs are stored.'
        )
    compile_args.add_argument(
        '--long-subs',
        '-S',
        required=True,
        default=None,
        help=(
            'Path to a long substitutions file, listing every SNP that is to '
            'be predicted. Same format as -s.'
            )
        )

    #   Add a switch for verbosity
    parser.add_argument(
        '--verbosity',
        '-v',
        required=False,
        dest='loglevel',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='Minimum verbosity level of messages printed.')
    args = parser.parse_args()
    return args


#   Here we validate the arguments
def validate_args(args, log):
    """A function that validates the arguments. For arguments that are
    filenames, it checks that they are readable. For directories, it checks
    that they are read/write. For usernames, it checks that they are valid
    email addresses. Prompts for username/password if they are not supplied
    on the command line. Validate input files for prediction."""
    #   Check the base argument. If it starts with something other than a /
    #   then it is a relative path, and we should fix it
    if not args['base'].startswith('/'):
        #   Add the cwd onto it, since the script fails otherwise
        args['base'] = os.path.join(os.getcwd(), args['base'])
    #   Then check the action
    #   If we are fetching, we have to check the username and base
    #   argparse should have checked for missing arguments by now
    #   If the arguments do not check out, return a message
    if args['action'] == 'setup':
        if args['list_species']:
            return (
                False,
                'The list of allowable species names is \n' +
                '\n'.join(SPECIES_LIST))
        if args['target'] not in SPECIES_LIST:
            return (
                False,
                ('The species name you provided is not in the list of '
                 'allowable species.'))
        #   Check the filename for the config file. It can be a relative path
        #   or start with a tilde.
        if not args['config'].startswith('/'):
            args['config'] = os.path.join(os.getcwd(), args['config'])
        elif args['config'].startswith('~'):
            #   os.path.expanduser() will transform ~/... into /home/user/...
            args['config'] = os.path.expanduser(args['config'])
        if not check_args.valid_dir(os.path.dirname(args['config'])):
            return (
                False,
                'You cannot create a configuration file in that directory.')
        if not check_args.valid_dir(args['base']):
            return (
                False,
                'Base directory is not readable/writable, or does not exist.')
    elif args['action'] == 'fetch':
        #   If config is suppled:
        if args['config']:
            if not file_funcs.file_exists(args['config'], log):
                return (
                    False,
                    'The specified configuration file does not exist!')
        #   If username is supplied:
        if args['user']:
            #   Check if it's valid
            if not check_args.valid_email(args['user']):
                return (
                    False,
                    'Username is not a valid e-mail address.')
        #   Username not supplied, and we need to access JGI
        elif not args['convert_only']:
            args['user'] = raw_input('Username for JGI Genomes Portal: ')
        #   Else, we only want to convert
        else:
            pass
        #   Same with password
        if args['password']:
            pass
        elif not args['convert_only']:
            args['password'] = getpass.getpass(
                'Password for JGI Genomes Portal: ')
        else:
            pass
        if not check_args.valid_dir(args['base']):
            return (
                False,
                'Base directory is not readable/writable, or does not exist.')
        else:
            pass
    #   Check the arguments passed to align
    elif args['action'] == 'align':
        #   If config is suppled:
        if args['config']:
            if not file_funcs.file_exists(args['config'], log):
                return (
                    False,
                    'The specified configuration file does not exist!')
        if not check_args.valid_dir(args['output']):
            return (
                False,
                'Output directory is not readable/writable, or does not exist.')
    #   Check arguments to predict
    elif args['action'] == 'predict':
        #   If config is suppled:
        if args['config']:
            if not file_funcs.file_exists(args['config'], log):
                return (
                    False,
                    'The specified configuration file does not exist!')
        if not check_args.valid_dir(args['output']):
            return (
                False,
                'Output directory is not readable/writable, or does not exist.')
        if not parse_input.valid_tree(args['tree'], log):
            return (
                False,
                'The input Newick tree is not valid.')
        if not parse_input.valid_msa(args['alignment'], log):
            return (
                False,
                'The input MSA file provided is not valid.')
        if not parse_input.parse_subs(args['substitutions'], log):
            return (
                False,
                'The input substitutions file provided is not valid.')
    return (args, None)


#   This is just a simple function that shows when the user does not supply
#   any arguments. This is an issue with Python 2.*, and has been "fixed" in
#   Python 3+.
def usage():
    """Print a usage message."""
    print '''Usage: BAD_Mutations.py <subcommand> <arguments>

where <subcommand> is one of 'setup', 'fetch', or 'predict.' This script will
download the necessary data to perform the likelihood ratio test (LRT) for
deleterious SNP prediction as described in Chun and Fay (2009) in Genome
Research. Because of the data sources used, this implementation is specific to
SNP annotation in plants.

The 'setup' subcommand will create a configuration file that contains paths to
requried executables and parameters for alignment. This is optional, but
recommended, as it makes batches of analysis much easier to standardize.

The 'fetch' subcommand will download gzipped CDS FASTA files from Phytozome,
unzip them, and convert them into BLAST databases. It requires a (free)
username and password for the JGI Genomes Portal. Check with Phytozome for
their data release and usage policies. Use 'BAD_Mutations.py fetch -h' for more
information.

The 'align' subcommand will run TBLASTX against all of the available BLAST
databases, collect orthologues that pass the E-value threshold, align them with
PASTA, and produce a phylogenetic tree.

The 'predict' subcommand will run the LRT with a given query sequence and a
list of affected codons.

Dependencies:
    Biopython
    tblastx (NCBI BLAST executables)
    requests (Python HTTP requests module)
    pasta (Phylogeny-aware sequence alignment)'''
    return
