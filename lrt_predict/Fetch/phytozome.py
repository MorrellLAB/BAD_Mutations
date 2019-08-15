#!/usr/bin/env python
"""Defines a class that is used to fetch data from Phytozome 10."""

#   Import python modules here
from xml.etree import ElementTree
import os
import subprocess
import tempfile

#   Import our helper scripts here
import lrt_predict.Fetch.phytozome_species as phytozome_species
import lrt_predict.Fetch.fetch as fetch
import lrt_predict.General.file_funcs as file_funcs


#   A new class to handle the requests to and responses from JGI
class Phytozome(fetch.Fetcher):
    """A class to handle fetching data from Phytozome 10, hosted by the Joint
    Genome Institute of the US Department of Energy.
    URL: http://phytozome.jgi.doe.gov

    Contains the following class attributes:
        JGI_LOGIN (str)       Sign in URL for JGI Genome Portal
        DL_BASE (str)         URL to the root of the downloads directory
        XML_URL (str)         URL to the XML file describing the downloads tree
        XML_DATA (dict)       HTTP GET variables used to request every species
                              for download
        FAILED_LOGIN (str)    Text to search for for identifying if a login
                              attempt has failed
        EXPIRED_ACCOUNT (str) Text to search for for identifying if an account
                              has an expired username/password

    Contains the following instance attributes:
        username (str)        User name for logging into JGI Genomes Portal
        password (str)        Password for the same
        session (Session)     HTTP Requests session for logging into JGI
        urls (list)           URLs to files to download
        md5s (list)           MD5 hashes of files for integrity checks

    Inherits the following attributes and methods from fetch.Fetcher:
        mainlog (logger)      Logging messages formatter and handler
        base (str)            Base directory for the data
        make_species_dir()    Create a directory for a species' data
        convert()             Convert the FASTA files to BLAST databases

    Contains the following methods:
        __init__(self, u, p, base, convertonly, verbose):
            Initialize the class with the username, password, base directory,
            whether or not to only convert to databases, and verbosity level.

        sign_on(self):
            Sign on to JGI Genomes Portal with the provided credentials.
            Searches for failed login attempts and expired username and
            password combinations. Returns a requests.Session object.

        get_xml_urls(self):
            Gets the URLs and th MD5s of the CDS files from the XML file from
            Phytozome. Stores these data in `urls' and `md5s' respectively.

        download_file(self, file):
            Fetches a remote file. Uses the requests.Session object to
            authenticate.

        fetch_cds(self):
            Iterates through the urls and md5s instance attributes and
            downloads the appropriate files. Checks the local MD5 against the
            remote MD5 and downloads the remote file if they differ. Appends
            the filenames of each updated file to the `to_convert' attribute.
    """
    #   These are common variables for every new Phytozome class
    #   They should not change from instance to instance
    JGI_LOGIN = 'https://signon-old.jgi.doe.gov/signon/create'
    DL_BASE = 'https://genome.jgi.doe.gov'
    XML_URL = 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=Phytozome'
    FAILED_LOGIN = 'Login and password do not match'
    EXPIRED_ACCOUNT = 'Sorry, your password has expired'
    TO_FETCH = phytozome_species.phyto_fetch

    def __init__(self, user, passwd, base,
                 convertonly, verbose):
        """Initialize the class with the username, password, base directory,
           whether or not to only convert to databases, and verbosity level."""
        fetch.Fetcher.__init__(self, base, verbose)
        self.username = user
        self.password = passwd
        self.mainlog.debug('Creating new instance of Phytozome')
        #   If we are only converting, then we don't have to sign on
        if convertonly:
            self.cookie = None
        else:
            self.cookie = self.sign_on()
        self.urls = []
        self.md5s = []
        return

    def sign_on(self):
        """Sign on to JGI Genomes Portal with the provided credentials.
           Searches for failed login attempts and expired username and
           password combinations. Returns a requests.Session object."""
        self.mainlog.debug('Logging in to JGI Genomes Portal with username ' +
                           self.username)
        #   Create a named temporary file to hold the cookie information
        cookie_file = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_JGI_Cookie_',
            suffix='.txt',
            delete=False)
        self.mainlog.debug('Cookies file: ' + cookie_file.name)
        #   Use cURL to log in to JGI, and create a cookie file that we will
        #   use to fetch all other files.
        cmd = [
            'curl',
            self.JGI_LOGIN,
            '--data-urlencode',
            'login=' + self.username,
            '--data-urlencode',
            'password=' + self.password,
            '-c',
            cookie_file.name
            ]
        #   Then we execute the command
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE
            )
        out, err = p.communicate()
        return cookie_file

    def get_xml_urls(self):
        """Gets the URLs and th MD5s of the CDS files from the XML file from
           Phytozome. Stores these data in `urls' and `md5s' respectively."""
        self.mainlog.debug('Fetching XML')
        #   Create another temporary named file for the XML output
        xml_out = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_JGI_XML_',
            suffix='.xml',
            delete=False)
        self.mainlog.debug('XML will be stored in ' + xml_out.name)
        #   Use cURL to download the XML, passing the cookies we generated
        #   earlier to authenticate.
        cmd = [
            'curl',
            self.XML_URL,
            '-b',
            self.cookie.name,
            '-o',
            xml_out.name
            ]
        #   Execute the command
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()
        self.mainlog.debug('cURL stdout: ' + out.decode('utf-8'))
        self.mainlog.debug('cURL stderr: ' + err.decode('utf-8'))
        #   Then, read the XML back from the file
        xml = xml_out.read()
        #   This suffix is what we want the filenames ending with
        #   this can change, depending on the target of the LRT
        suffix = '.cds.fa.gz'
        #   Use HTTP GET to fetch the XML from Phytozome's server
        #   This is also a response obkect
        self.mainlog.debug('The XML I got was \n\n' + xml)
        #   Create an element tree out of it, so we can easily step
        #   through the data
        xml_tree = ElementTree.fromstring(xml)
        #   Step through it and extract all CDS URLs
        for elem in xml_tree.findall('.//file'):
            #   if the URL ends in a certain suffix, then save it
            if elem.attrib.get('url').endswith(suffix):
                url = elem.attrib.get('url')
                md5 = elem.attrib.get('md5')
                #   Check to see that the file is in the list of
                #   species to download
                local_filename = file_funcs.local_name(url)
                species_name = file_funcs.species_name(local_filename)
                if species_name in self.TO_FETCH:
                    self.urls.append(url)
                    self.md5s.append(md5)
        self.mainlog.debug('Found ' + str(len(self.urls)) + ' files to fetch')
        return

    def download_file(self, url):
        """Fetches a remote file. Uses the cookies file from cURL to
           authenticate."""
        self.mainlog.debug('Fetching ' + url)
        #   Build the full URL to fetch
        full_url = self.DL_BASE + url
        #   And build the command to download it
        cmd = [
            'curl',
            full_url,
            '-b',
            self.cookie.name,
            '-o',
            file_funcs.local_name(url)
            ]
        #   Then download it
        p = subprocess.Popen(
            cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
            )
        out, err = p.communicate()
        self.mainlog.debug('Done fetching ' + url)
        return

    #    A function to fetch the CDS files
    def fetch_cds(self):
        """Iterates through the urls and md5s instance attributes and
           downloads the appropriate files. Checks the local MD5 against the
           remote MD5 and downloads the remote file if they differ. Appends
           the filenames of each updated file to the `to_convert' attribute."""
        self.mainlog.debug('Downloading files from ' +
                           str(len(self.urls)) +
                           ' species')
        for u, m in zip(self.urls, self.md5s):
            #   Get a local name of the CDS
            lname = file_funcs.local_name(u)
            target_dir = self.make_species_dir(u)
            os.chdir(target_dir)
            #   check to see if the file already exists
            if file_funcs.file_exists(lname, self.mainlog):
                #   Get the md5
                lmd5 = file_funcs.calculate_md5(lname, self.mainlog)
                #   Compare the MD5s
                md5s_same = file_funcs.checksum_is_same(lmd5, m, self.mainlog)
                #   If they are the same, skip it, and move on
                if md5s_same:
                    self.mainlog.info(lname + ' is current. Skipping.')
                    continue
                else:
                    self.mainlog.info(lname + ' is out of date. Downloading.')
                    #   Try to download it until the MD5s check out
                    same = False
                    while not same:
                        self.download_file(u)
                        new_lmd5 = file_funcs.calculate_md5(
                            lname,
                            self.mainlog)
                        same = file_funcs.checksum_is_same(
                            new_lmd5,
                            m,
                            self.mainlog)
                    #   Tack it onto the list of files to convert
                    self.to_convert.append(
                        os.path.join(
                            self.base,
                            target_dir,
                            lname)
                        )
            else:
                self.mainlog.info(lname + ' does not yet exist. Downloading.')
                #   And the same procedure as if the file were updated
                same = False
                while not same:
                    self.download_file(u)
                    new_lmd5 = file_funcs.calculate_md5(
                        lname,
                        self.mainlog)
                    same = file_funcs.checksum_is_same(
                        new_lmd5,
                        m,
                        self.mainlog)
                self.to_convert.append(
                    os.path.join(
                        self.base,
                        target_dir,
                        lname))
        self.mainlog.info('Done downloading CDS files from Phytozome.')
        return
