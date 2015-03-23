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
from ..General import set_verbosity
from ..General import check_modules

#   A new class to handle the requests to and responses from Ensembl
class EnsemblPlants:
    #   When creating a new EnsemblPlants class, we have these following pieces of
    #   data created and attached to the class
    def __init__(self, base, convertonly, verbose):
        self.mainlog = set_verbosity.verbosity(__name__, verbose)
        self.mainlog.debug('Creating new instance of EnsemblPlants')
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

    #   A function to convert downloaded files to BLAST databases
    def convert(self):
        #   What is the path to the makeblastdb executable?
        makeblastdb_path = check_modules.check_executable('makeblastdb')
        #   Check if the list of updated CDS files is empty or not
        if not self.to_convert:
            #   If it is empty, then populate it with all of them
            fname_list = file_funcs.get_file_by_ext(self.base, '.fa.gz', self.filelog)
        else:
            fname_list = self.to_convert
        #   for each one
        for f in fname_list:
            out, error = format_blast.format_blast(makeblastdb_path, f)
            self.mainlog.info('stdout: \n' + out)
            self.mainlog.info('stderr: \n' + error)
        return
