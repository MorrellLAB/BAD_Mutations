#!/usr/bin/env python
#   This script will take the URLs of the CDS files and download them from
#   Phytozome.

#   To handle operating system commands like mkdir
import os

#   Define a base URL for downloading
DL_BASE = 'http://genome.jgi.doe.gov'

#   A function that generates a local filename from a remote one
def local_name(url):
    #   We assume that the remote url is separated by / and that the file name
    #   is the final part of the url
    return url.split('/')[-1]


#   A function that finds the species name given the filename
def species_name(fname):
    #   Again, assume this time that the species name is the first field in
    #   the _-delimited filename
    return fname.split('_')[0]

#   A function that downloads the CDS file
def dl(session, url, localname):
    #   With stream=True, it downloads the response right away
    r = session.get(DL_BASE + url, stream=True)
    #   Save the file
    with open(localname, 'wb') as f:
        #   Take the file in pieces
        for chunk in r.iter_content(chunk_size=1024):
            #   Empty chunks are for keepalive purposes, we don't save those
            if chunk:
                #   Write the file to disk
                f.write(chunk)
                #   and flush the buffer
                f.flush()
    return
