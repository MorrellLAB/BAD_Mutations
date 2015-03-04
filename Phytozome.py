#!/usr/bin/env python

#   Script to download CDS files from phytozome.net
#   requires a user name and password for the JGI phyotzome portal (Free)
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


#   To take arguments
import sys
#   Import our argument parsing script
import Parse_Args as PA
#   Import the Phytozome signon and URL parsing script
import Get_URLs as GU
#   Import the downloading script
import Fetch_CDS as FC
#   Import the directory handling script
import Dir_Funcs as DF

#   Main function
def main():
    #   Parse the arguments
    arguments = PA.parse_args()
    #   create the base for the LRT
    DF.makebase(arguments.base)
    #   Start a new login session to Phyotozome
    session = GU.signon(arguments.user, arguments.password)
    #   Get the list of all URLs from the downloads tree
    urls = GU.extract_all_urls(session)
    #   Get the CDS only
    cds = GU.extract_cds_urls(urls)
    #   What are the filenames?
    cds_names = []
    for c in cds:
        #   Download the files
        fname = FC.dl(session, c)
        print 'Downloaded ' + fname
        cds_names.append(fname)
    #   Now, move all the files into their correct directories
    for c in cds_names:
        #   Build the directory name from the filename
        #   the species is the first field, separated by _
        spname = c.split('_')[0]
        #   Create the species directory
        spdir = DF.make_species_dir(arguments.base, spname)
        #   And move the file into it
        DF.move_file(c, spdir)
    print "Done!"
    return


#   Do the work here
main()
