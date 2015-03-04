#!/usr/bin/env python

#   Script to download CDS files from phytozome.net
#   requires a user name and password for the JGI phyotzome portal (Free)

#   To take arguments
import sys
#   Import our argument parsing script
import Parse_Args as PA
#   Import the Phytozome signon and URL parsing script
import Get_URLs as GU
#   Import the downloading script
import Fetch_CDS as FC

#   Main function
def main():
    #   Parse the arguments
    arguments = PA.parse_args()
    #   Start a new login session to Phyotozome
    session = GU.signon(arguments.user, arguments.password)
    #   Get the list of all URLs from the downloads tree
    urls = GU.extract_all_urls(session)
    #   Get the CDS only
    cds = GU.extract_cds_urls(urls)
    for c in cds:
        fname = FC.dl(session, c)
        print 'Downloaded ' + fname

#   Do the work here
main()
