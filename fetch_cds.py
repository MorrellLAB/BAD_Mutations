#!/usr/bin/env python

import requests
from xml.etree import ElementTree

test_fil_url ='http://genome.jgi.doe.gov/PhytozomeV10/download/_JAMO/53112a6949607a1be0055904/Sitalica_164_v2.1.cds_primaryTranscriptOnly.fa.gz'

def download_file(url):
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter
    r = requests.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024): 
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    return local_filename

download_file(test_fil_url)
