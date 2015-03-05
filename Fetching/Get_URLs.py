#!/usr/bin/env python
#   A script to get URLs for the CDS files out of Phytozome XML

#   To handle HTTP requests
try:
    import requests
except ImportError:
    print 'Error! You need to have the requests module installed.'
    exit(1)
#   To parse XML files
from xml.etree import ElementTree

#   Some important URLs defined here
sign_on_page = 'https://signon.jgi.doe.gov/signon/create'
xml_dir = 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory'
#   These are the HTTP GET variables that are needed to get the full XML tree
xml_data = {'organism':'PhytozomeV10'}
#   Text to check for if a login was successful or not
#   This is VERY delicate, if JGI changes the wording on their form,
#   then this will have to be updated.
failed_login_text = 'Login and password do not match'

#   A function to sign onto Phyotozome.net, and maintain an active session
def signon(username, password):
    #   We create a new session object, which maintains our login credentials
    s = requests.Session()
    #   Create a dictionary that contains the data to be sent to the login form
    d = {'login': username, 'password': password}
    #   We sign on to Phytozome with our username and password with HTTP POST
    #   This creates a Response object
    r = s.post(sign_on_page, data=d)
    #   Check to see if there was a successful login
    #   This is VERY delicate, if JGI changes the wording on their form,
    #   then this will have to be updated.
    if failed_login_text in r.text:
        return None
    else:
        #   Return the Session object so we can use it to fetch data later
        return s


#   A function to get all URLs out of the XML file
#   Borrowed from YR's code
def extract_all_urls(session):
    #   Get the XML file itself, in a Response object
    xml = session.get(xml_dir, params=xml_data)
    #   Next, we create a tree of elements out of it
    #   Response.text contains the Unicode text of the response
    xmltree = ElementTree.fromstring(xml.text)
    #   We step through the tree and save all URLs from the XML
    urls = []
    #   And the md5 sums
    md5s = []
    for e in xmltree.findall('.//file'):
        urls.append((e.attrib.get('url'))
        md5s.append(e.attrib.get('md5'))
    #   And return the list of URLs
    return (urls, md5s)


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
    return cds
