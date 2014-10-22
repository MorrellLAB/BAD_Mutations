#!/usr/bin/env python
import requests

def sign_on(user,password):
    s = requests.session()
    payload = {'login':user,'password':password}
    signon_page ="https://signon.jgi.doe.gov/signon/create"
    r = s.post(signon_page, data=payload)
    return s

def download_file(url,s):
    local_filename = url.split('/')[-1]
    r = s.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024): 
            if chunk:
                f.write(chunk)
                f.flush()
    return local_filename

def fetch_xml(s):
    url = "http://genome.jgi.doe.gov/ext-api/downloads/get-directory"
    payload = {'organism':'PhytozomeV10'}
    xml = s.get(url,params=payload)
    return xml

#def main():
#    args = parse_args()
#    s = sign_on(args)
#    fetch_xml(s)

#main()
