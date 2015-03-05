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
#   - Get the MD5 of the remote files and compare to MD5 of local files to
#     see if we need to actually download them again.
#   - Integrate creation of BLAST databases


#   To handle passwords
import getpass
#   To cd
import os
#   Import our argument parsing script
from General import parse_prgs
#   Import the directory handling script
from General import dir_funcs
#   Import our file handling script
from General import file_checks
#   Import the Phytozome signon and URL parsing script
from Fetching import get_urls
#   Import the downloading script
from Fetching import fetch_cds

#   A function to perform the fetching
def fetch(base, user, password):
    #   create the base for the LRT
    dir_funcs.makebase(base)
    #   cd into the base. This should work since we already checked the
    #   permissions with the function call above
    os.chdir(base)
    #   Check for an entered password, else take it
    if password:
        pass
    else:
        password = getpass.getpass('Password: ')
    #   Start a new login session to Phyotozome
    session = get_urls.signon(user, password)
    #   If session is NoneType, then there was a login problem
    if not session:
        print 'There was a problem logging in to JGI Genomes Portal. Check your\
username and password and try again.'
        exit(1)
    #   Get the list of all URLs and md5s from the downloads tree
    urls, md5s = get_urls.extract_all_urls(session)
    #   Get the CDS only
    cds = get_urls.extract_cds_urls(urls)
    #   What are the filenames?
    cds_names = []
    for c in cds:
        #   What is the local name?
        local_name = fetch_cds.localname(c)
        #   And the species name?
        species_name = fetch_cds.speciesname(local_name)
        #   Do these files exist already?
        spdir = dir_funcs.make_species_dir(base, species_name)
        #   Download the files
        fetch_cds.dl(session, c)
        print 'Downloaded ' + fname
        cds_names.append(fname)
    #   Now, move all the files into their correct directories
    #   Save a MD5 in each directory
    for c, m in zip(cds_names, md5s):
        #   Build the directory name from the filename
        #   the species is the first field, separated by _
        spname = c.split('_')[0]
        #   Create the species directory
        spdir = dir_funcs.make_species_dir(base, spname)
        #   And move the file into it
        dir_funcs.move_file(c, spdir)

    print "Done!"
    return


#   A function to do the prediction
def predict():
    pass


#   Main function
def main():
    #   Parse the arguments
    arguments = parse_args.parse_args()
    #   Which command was invoked?
    if arguments.action == 'fetch':
        fetch(arguments.base, arguments.user, arguments.password)
    elif arguments.action == 'predict':
        predict()
    return


#   Do the work here
main()
