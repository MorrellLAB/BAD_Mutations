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
#   Import our argument parsing script
from General import parse_args
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
        print 'There was a problem logging in to JGI Genomes Portal. Check your username and password and try again.'
        exit(1)
    #   Get the list of all URLs and md5s from the downloads tree
    urls, md5s = get_urls.extract_all_urls(session)
    #   Get the CDS only
    cds, md5s = get_urls.extract_cds_urls(urls, md5s)
    #   What are the filenames?
    cds_names = []
    for c, m in zip(cds, md5s):
        #   What is the local name?
        local_name = fetch_cds.local_name(c)
        #   And the species name?
        species_name = fetch_cds.species_name(local_name)
        print 'Fetching ' + species_name + ' ...'
        #   Create the species directory
        spdir = dir_funcs.make_species_dir(base, species_name)
        #   cd into that directory
        os.chdir(spdir)
        #   Do these files exist already?
        if file_checks.file_exists(local_name):
            #   If so, check for the md5 file
            if file_checks.file_exists(local_name + '.md5'):
                #   If it's there, compare the two MD5s
                md5s_same = file_checks.md5_is_same(local_name + '.md5', m)
            else:
                #   If it's not there, calculate it
                calc_md5 = file_checks.calculate_md5(local_name)
                #   And then compare
                md5s_same = file_checks.md5_is_same(calc_md5, m)
            #   Then we check if the md5s are the same. If they are the same,
            #   then we skip the file, and move on. Else, we download a new one
            if md5s_same:
                print local_name + ' already exists and MD5s are identical, skipping ...'
                continue
            else:
                print local_name + ' has different MD5 than server. Downloading ...'
                #   A variable to tell if we've sucessfully downloaded the file
                same = False
                while not same:
                    fetch_cds.dl(session, c, local_name)
                    r_md5 = m
                    l_md5 = file_checks.calculate_md5(local_name)
                    print r_md5, l_md5
                    same = file_checks.md5_is_same(l_md5, r_md5)
                #   Save the md5 for later
                handle = open(local_name + '.md5', 'w')
                handle.write(l_md5 + '\n')
                handle.close()
        #   If the files don't exist already, then we just download them
        else:
            #   a variable to tell if we've checked out
            same = False
            #   While making sure we downloaded reliably
            while not same:
                fetch_cds.dl(session, c, local_name)
                r_md5 = m
                l_md5 = file_checks.calculate_md5(local_name)
                same = file_checks.md5_is_same(l_md5, r_md5)
            #   Save the md5 for later
            handle = open(local_name + '.md5', 'w')
            handle.write(l_md5 + '\n')
            handle.close()
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
