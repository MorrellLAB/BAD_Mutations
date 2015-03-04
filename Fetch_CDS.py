#!/usr/bin/env python
#   This script will take the URLs of the CDS files and download them from
#   Phytozome.

#   To handle operating system commands like mkdir
import os

#   Define a base URL for downloading
dl_base = 'http://genome.jgi.doe.gov'
#   The base directory
lrt_base = '.'

#   A function that downloads the CDS file
def dl(session, url):
    #   The local filename should be the basename of the remote file
    localname = url.split('/')[-1]
    #   With stream=True, it downloads the response right away
    r = session.get(dl_base + url, stream=True)
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
    return localname
