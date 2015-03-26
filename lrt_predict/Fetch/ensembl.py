#!/usr/bin/env python

#   Import python modules here
import ftplib
import logging
import os
from StringIO import StringIO

#   Import our helper scripts here
from ..General import dir_funcs
from ..General import file_funcs
import format_blast
import ensembl_species
from ..General import set_verbosity
from ..General import check_modules

#   A new class to handle the requests to and responses from Ensembl
class EnsemblPlants:
    #   When creating a new EnsemblPlants class, we have these following pieces of
    #   data created and attached to the class
    ENSEMBL_URL = 'ftp.ensemblgenomes.org'
    ENSEMBL_PLANT_BASE = '/pub/plants/current/fasta/'
    ENSEMBL_TO_FETCH = ensembl_species.ensembl_fetch
    def __init__(self, base, convertonly, verbose):
        self.mainlog = set_verbosity.verbosity(__name__, verbose)
        self.mainlog.debug('Creating new instance of EnsemblPlants')
        #   If we are only converting, then we don't have to sign on
        if convertonly:
            self.session = None
        else:
            self.session = self.sign_on()
        self.urls = []
        self.cksums = []
        self.base = base
        self.to_convert = []
        #   Do stuff for the base directory
        self.dirlog = set_verbosity.verbosity('Dir_Funcs', verbose)
        dir_funcs.makebase(base, self.dirlog)
        #   And stuff for the file operations
        self.filelog = set_verbosity.verbosity('File_Funcs', verbose)

    #   A function to lot onto the FTP server
    def sign_on(self):
        ftp_session = ftplib.FTP(self.ENSEMBL_URL)
        #   login() without any args, use anonymous FTP by default
        ftp_session.login()
        ftp_session.cwd(self.ENSEMBL_PLANT_BASE)
        return ftp_session

    #   A function to just get a file.
    def get_file(self, fname):
        handle = open(file_funcs.local_name(fname), 'wb')
        self.session.retrbinary('RETR ' + fname, handle.write)
        handle.close()
        return

    #   A function to build the list of URLs to fetch
    def get_ftp_urls(self):
        #   nlst() returns a listing of the directory contents
        #   We are going to assume that each of these are directories, since
        #   that is the way Ensembl has their server set up...
        for d in self.session.nlst():
            #   If it's in our list of species to download...
            if d in self.ENSEMBL_TO_FETCH:
                self.mainlog.debug('Attempting to cd into ' + self.ENSEMBL_PLANT_BASE + d + '/cds/')
                #   cd into the CDS directory there
                self.session.cwd(self.ENSEMBL_PLANT_BASE + d + '/cds/')
                #   Get the contents
                listing = self.session.nlst()
                #   Then find the one that ends in .fa.gz
                for l in listing:
                    if l.endswith('.fa.gz'):
                        #   Tack the full path onto the list of URLS
                        self.mainlog.debug('Appending ' + self.session.pwd() + '/' + l + ' onto list of files to fetch.')
                        self.urls.append(self.session.pwd() + '/' + l)
                    #   get the checksum
                    elif l == 'CHECKSUMS':
                        #   We create a new StringIO instance, which can read/write
                        #   string data like file data
                        c = StringIO()
                        self.session.retrbinary('RETR ' + l, c.write)
                        #   Then we get the contents of the CHECKSUMS file
                        #   Split on newlines, since there is one line per file entry
                        c_str = c.getvalue()
                        c.close()
                        lines = c_str.strip().split('\n')
                        #   Check for the file that ends with .fa.gz 
                        for line in lines:
                            if line.endswith('.fa.gz'):
                                #   The first field is the 16-bit checksum
                                #   Cast to integer
                                self.cksums.append(int(line.split()[0]))
        return

    #   A function to download the files into the proper directories
    def download_files(self):
        #   For each URL we have:
        for u, c in zip(self.urls, self.cksums):
            #   We have to get the species directories as in the Phytozome class
            spname = file_funcs.ensembl_species_name(u)
            #   Then we have to make that directory
            target_dir = dir_funcs.make_species_dir(self.base, spname, self.dirlog)
            #   cd into it
            os.chdir(target_dir)
            #   What is the local file name?
            lname = file_funcs.local_name(u)
            #   If it exists, we check if the checksums are the same
            if file_funcs.file_exists(lname, self.filelog):
                local_cksum = file_funcs.calculate_crc32(lname, self.filelog)
                crc32_same = file_funcs.checksum_is_same(local_cksum, c, self.filelog)
                if crc32_same:
                    self.mainlog.info(lname + ' already exists and is current, skipping.')
                    continue
                else:
                    self.mainlog.info(lname + ' exists, but is out of date. Updating.')
                    same = False
                    while not same:
                        self.get_file(u)
                        new_local_cksum = file_funcs.calculate_crc32(lname, self.filelog)
                        same = file_funcs.checksum_is_same(new_local_cksum, c, self.filelog)
                    #   And save a record for those that need to be converted
                    self.to_convert.append(os.path.join(self.base, target_dir, lname))
            #   If the file doesn't exist, then it's the same as if the checksum were different
            else:
                self.mainlog.info(lname + ' does not exist. Downloading.')
                same = False
                while not same:
                    self.get_file(u)
                    new_local_cksum = file_funcs.calculate_crc32(lname, self.filelog)
                    same = file_funcs.checksum_is_same(new_local_cksum, c, self.filelog)
                self.to_convert.append(os.path.join(self.base, target_dir, lname))
        self.mainlog.info('Done downloading CDS files from Ensembl.')
        #   We are done with the FTP connection, log out
        self.session.quit()
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
