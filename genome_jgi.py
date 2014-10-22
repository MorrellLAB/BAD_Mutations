#!/usr/bin/env python
import requests
from xml.etree import ElementTree

class JGIUtils(object):

    def __init__(self,user,password):
        self.s = self.sign_on(user,password)

    def sign_on(self,user,password):
        s = requests.session()
        payload = {'login':user,'password':password}
        signon_page ="https://signon.jgi.doe.gov/signon/create"
        r = s.post(signon_page, data=payload)
        return s

    def download_file(self,url,s):
        local_filename = url.split('/')[-1]
        r = s.get(url, stream=True)
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024): 
                if chunk:
                    f.write(chunk)
                    f.flush()
        return local_filename

    def fetch_xml(self):
        url = "http://genome.jgi.doe.gov/ext-api/downloads/get-directory"
        payload = {'organism':'PhytozomeV10'}
        xml = self.s.get(url,params=payload)
        return xml.content

    def fetch_url_list(self):
        url_list = []
        xml = self.fetch_xml()
        tree = ElementTree.fromstring(xml)
        for elem in tree.findall('.//file'):
            url = elem.attrib.get('url')
            url_list.append(url)

        return url_list

    def fetch_cds_list(self):
        cds_url_list = []
        urls = self.fetch_url_list()
        for url in urls:
            part = url.split(".")
            for p in part:
                if p == "cds":
                    cds_url_list.append(url)
        return cds_url_list

#def main():
#    args = parse_args()
#    s = sign_on(args)
#    fetch_xml(s)

#main()
