#!/usr/bin/env python

import requests
import argparse
from xml.etree import ElementTree
from requests.auth import HTTPDigestAuth

def parse_args():
    print "parse args"
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--foo', help='foo help')
    args = parser.parse_args()
    parser.print_help()

def sign_on():
    login_data = {'commit':'1','login':'ribb0013','password':''}
    s = requests.Session()
    s.auth = ('user', 'pass')
    s.headers.update({'x-test': 'true'})

    signon_page ="https://signon.jgi.doe.gov/signon"
    requests.post(signon_page,login_data,verify=True)
    return s


def download_file(url,s):
    local_filename = url.split('/')[-1]
    r = s.get(url, stream=True)
    print r.headers
    print r.request.headers
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024): 
            if chunk:
                f.write(chunk)
                f.flush()
    return local_filename

def main():
    parse_args()
    test_fil_url ='http://genome.jgi.doe.gov/PhytozomeV10/download/_JAMO/53112a6949607a1be0055904/Sitalica_164_v2.1.cds_primaryTranscriptOnly.fa.gz'

    s = sign_on()
    download_file(test_fil_url,s)

parse_args()
#main()
