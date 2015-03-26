#!/usr/bin/env python

#   Import python modules here
import requests
from xml.etree import ElementTree
import logging
import os

#   Import our helper scripts here
import dir_funcs
import file_funcs
import format_blast
import phytozome_species
from ..General import set_verbosity
from ..General import check_modules

#   A new class to handle the requests to and responses from JGI
class Phytozome:
    #   These are common variables for every new Phytozome class
    #   They should not change from instance to instance
    JGI_LOGIN = 'https://signon.jgi.doe.gov/signon/create'
    DL_BASE = 'http://genome.jgi.doe.gov'
    XML_URL = 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory'
    XML_DATA = {'organism':'PhytozomeV10'}
    FAILED_LOGIN = 'Login and password do not match'
    TO_FETCH = phyozome_species.phyto_fetch
    #   When creating a new Phytozome class, we have these following pieces of
    #   data created and attached to the class
    def __init__(self, u, p, base, convertonly, verbose):
        self.username = u
        self.password = p
        self.mainlog = set_verbosity.verbosity(__name__, verbose)
        self.mainlog.debug('Creating new instance of Phytozome')
        #   If we are only converting, then we don't have to sign on
        if convertonly:
            self.session = None
        else:
            self.session = self.sign_on()
        self.urls = []
        self.md5s = []
        self.base = base
        self.to_convert = []
        #   Do stuff for the base directory
        self.dirlog = set_verbosity.verbosity('Dir_Funcs', verbose)
        dir_funcs.makebase(base, self.dirlog)
        #   And stuff for the file operations
        self.filelog = set_verbosity.verbosity('File_Funcs', verbose)

    #   A function to sign on
    def sign_on(self):
        self.mainlog.debug('Logging in to JGI Genomes Portal with username ' + self.username)
        #   Start a new session, this will store our login information
        s = requests.Session()
        #   This is the data we will send to phytozome
        d = {'login': self.username, 'password': self.password}
        #   Get the response data back - the results of our login attempt
        r = s.post(self.JGI_LOGIN, data=d)
        if self.FAILED_LOGIN in r.text:
            self.mainlog.critical('Could not log into JGI Genomes Portal. Check username and password')
            exit(1)
        self.mainlog.debug('Successfully logged in')
        return s

    #   A function to get all the CDS URLs out of the Phytozome XML
    def get_xml_urls(self):
        self.mainlog.debug('Fetching XML')
        #   This suffix is what we want the filenames ending with
        #   this can change, depending on the target of the LRT
        suffix = '.cds.fa.gz'
        #   Use HTTP GET to fetch the XML from Phytozome's server
        #   This is also a response obkect
        xml = self.session.get(self.XML_URL, params=self.XML_DATA)
        #   Create an element tree out of it, so we can easily step
        #   through the data
        xml_tree = ElementTree.fromstring(xml.text)
        #   Step through it and extract all CDS URLs
        for e in xml_tree.findall('.//file'):
            #   if the URL ends in a certain suffix, then save it
            if e.attrib.get('url').endswith(suffix):
                url = e.attrib.get('url')
                md5 = e.attrib.get('md5')
                #   Check to see that the file is in the list of species to download
                local_filename = file_funcs.local_name(url)
                species_name = file_funcs.species_name(local_filename)
                if species_name in self.TO_FETCH:
                    self.urls.append(url)
                    self.md5s.append(md5)
        self.mainlog.debug('Found ' + str(len(self.urls)) + ' files to fetch')
        return

    #   A function to download a file
    def download_file(self, url):
        self.mainlog.debug('Fetching ' + url)
        #   With stream=True, it downloads the response right away
        r = self.session.get(self.DL_BASE + url, stream=True)
        #   Save the file
        with open(file_funcs.local_name(url), 'wb') as f:
            #   Take the file in pieces
            for chunk in r.iter_content(chunk_size=1024):
                #   Empty chunks are for keepalive purposes, we don't save those
                if chunk:
                    #   Write the file to disk
                    f.write(chunk)
                    #   and flush the buffer
                    f.flush()
            self.mainlog.debug('Done fetching ' + url)
        return

    #    A function to fetch the CDS files
    def fetch_cds(self):
        self.mainlog.debug('Downloading files from ' + str(len(self.urls)) + ' species')
        for u, m in zip(self.urls, self.md5s):
            #   Get a local name of the CDS
            lname = file_funcs.local_name(u)
            #   Get the species name from the filename
            species_name = file_funcs.species_name(lname)
            #   Create the species directory
            target_dir = dir_funcs.make_species_dir(self.base, species_name, self.dirlog)
            #   cd into it
            os.chdir(target_dir)
            #   check to see if the file already exists
            if file_funcs.file_exists(lname, self.filelog):
                #   Get the md5
                lmd5 = file_funcs.calculate_md5(lname, self.filelog)
                #   Compare the MD5s
                md5s_same = file_funcs.checksum_is_same(lmd5, m, self.filelog)
                #   If they are the same, skip it, and move on
                if md5s_same:
                    self.mainlog.info(lname + ' already exists, and is current. Skipping.')
                    continue
                else:
                    self.mainlog.info(lname + ' exists, but it out of date. Downloading.')
                    #   Try to download it until the MD5s check out
                    same = False
                    while not same:
                        self.download_file(u)
                        new_lmd5 = file_funcs.calculate_md5(lname, self.filelog)
                        same = file_funcs.checksum_is_same(new_lmd5, m, self.filelog)
                    #   Tack it onto the list of files to convert
                    self.to_convert.append(os.path.join(self.base, target_dir, lname))
            else:
                self.mainlog.info(lname + ' does not yet exist. Downloading.')
                #   And the same procedure as if the file were updated
                same = False
                while not same:
                    self.download_file(u)
                    new_lmd5 = file_funcs.calculate_md5(lname, self.filelog)
                    same = file_funcs.checksum_is_same(new_lmd5, m, self.filelog)
                self.to_convert.append(os.path.join(self.base, target_dir, lname))
        self.mainlog.info('Done downloading CDS files from Phytozome.')
        return

    #   A function to convert downloaded files to BLAST databases
    def convert(self):
        #   What is the path to the makeblastdb executable?
        makeblastdb_path = check_modules.check_executable('makeblastdb')
        #   Check if the list of updated CDS files is empty or not
        if not self.to_convert:
            #   If it is empty, then populate it with all of them
            fname_list = file_funcs.get_file_by_ext(self.base, '.cds.fa.gz', self.filelog)
        else:
            fname_list = self.to_convert
        #   for each one
        for f in fname_list:
            out, error = format_blast.format_blast(makeblastdb_path, f)
            self.mainlog.info('stdout: \n' + out)
            self.mainlog.info('stderr: \n' + error)
        return
