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

###############################################################################
##### Future? #################################################################
###############################################################################
#   - Integrate creation of BLAST databases


#   To handle passwords
import getpass
#   To cd
import os
#   to check arguments
import sys
#   Import the main fetching script
import lrt_predict.Fetching.mainfetch as mainfetch
#   Import our argument parsing script
from lrt_predict.General import parse_args

#   A function to do the prediction
def predict():
    pass


#   Main function
def main():
    #   Parse the arguments
    #   First, a check to see if any arguments were sent at all
    #   If not, then print the usage and exit
    if not sys.argv[1:]:
        parse_args.usage()
        exit(1)
    arguments = parse_args.parse_args()
    arguments_valid, msg = parse_args.validate_args(arguments)
    #   If we got a return value that isn't False, then our arguments are good
    if arguments_valid:
        #   Which command was invoked?
        if arguments_valid.action == 'fetch':
            if arguments_valid.convert_only:
                mainfetch.convert(arguments_valid.base, None)
            else:
                updated_cds = mainfetch.fetch(
                    arguments_valid.base,
                    arguments_valid.user,
                    arguments_valid.password,
                    arguments_valid.fetch_only)
                mainfetch.convert(updated_cds)
        elif arguments_valid.action == 'predict':
            predict()
    else:
        print 'Error!', msg
    return


#   Do the work here
main()
