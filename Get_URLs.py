#!/usr/bin/env python
#   A script to get URLs for the CDS files out of Phytozome XML

#   To handle HTTP requests
import requests
#   To parse XML files
from xml.etree import ElementTree

#   Some important URLs defined here
sign_on_page = 'https://signon.jgi.doe.gov/signon/create'
xml_dir = 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory'
#   These are the HTTP GET variables that are needed to get the full XML tree
xml_data = {'organism':'PhytozomeV10'}

#   A function to sign onto Phyotozome.net, and maintain an active session
def signon(username, password):
    #   We create a new session object, which maintains our login credentials
    s = requests.Session()
    #   Create a dictionary that contains the data to be sent to the login form
    d = {'user': username, 'password': password}
    #   We sign on to Phytozome with our username and password with HTTP POST
    #   This creates a Response object
    r = s.post(sign_on_page, data=d)
    #   Return the Response object so we can requests data from it later
    return(r)


#   A function to get all URLs out of the XML file
def extract_all_urls(response):
    #   Get the XML file itself
    xml = response.get(xml_dir, params=xml_data)
    #   Next, we create a tree of elements out of it
    #   Response.text contains the Unicode text of the response
    xmltree = ElementTree.fromstring(xml.text)
    #   We step through the tree and save all URLs from the XML
    urls = []
    for e in xmltree.findall('.//file'):
        urls.append(e.attrib.get('url'))
    #   And return the list of URLs
    return(urls)


#   A function to extract the URLs to the CDS files from the XML
def extract_cds_urls(url_list):
    #   The CDS files end in .cds.fa.gz
    #   A little messy, but it works...
    suffix = '.cds.fa.gz'
    cds = []
    for u in url_list:
        if u.endswith(suffix):
            cds.append(u)
    #   And return the list of cds URLs
    return(cds)
