#!/usr/bin/env python
"""A class to handle downloads from Ensembl Plants."""

#   Import python modules here
import ftplib
import os
from StringIO import StringIO

#   Import our helper scripts here
import lrt_predict.Fetch.fetch as fetch
import lrt_predict.General.file_funcs as file_funcs
import lrt_predict.Fetch.ensembl_species as ensembl_species


#   A new class to handle the requests to and responses from Ensembl
class EnsemblPlants(fetch.Fetcher):
    """A class to handle fetching data from Ensembl Plants.
    URL: http://plants.ensembl.org/

    Contains the following class attributes:
        ENSEMBL_URL (str)        Base FTP URL for Ensembl
        ENSEMBL_PLANT_BASE (str) Directory on server for plant genome data
        ENSEMBL_TO_FETCH (list)  Species to download

    Contains the following instance attributes:
        session (Session)     HTTP Requests session for logging into JGI
        urls (list)           URLs to files to download
        cksums (list)         CRC sum of files for integrity checks
        to_covert (list)      Files to convert from FASTA to BLAST databases

    Inherits the following attributes and methods from fetch.Fetcher:
        mainlog (logger)      Logging messages formatter and handler
        base (str)            Base directory for the data
        make_species_dir()    Create a directory for a species' data
        convert()             Convert the FASTA files to BLAST databases

    Contains the following methods:
        __init__(self, base, convertonly, verbose):
            Initialize the class with the base directory

        sign_on(self):
            Initialize FTP connection to Ensembl Plants with ENSEMBL_URL
            and navigate to ENSEMBL_PLANT_BASE.

        get_file(self, fname):
            Download the file specified by `fname'

        get_ftp_urls(self):
            Build a list of FTP URLs to fetch. Currently searches for all
            files listed in a CDS directory. Also reads the CRC sums out of
            the CHECKSUMS file in each FTP directory. Save the URL in the
            `urls' attribute, and the checksums in the `cksums' attribute.

        dowload_files(self):
            Iterate through the list of URLs and download the appropriate
            files. Computes the CRC sum of existing files and compares them to
            the remote checksum to decide whether or not to to download.
    """

    ENSEMBL_URL = 'ftp.ensemblgenomes.org'
    ENSEMBL_PLANT_BASE = '/pub/plants/current/fasta/'
    ENSEMBL_TO_FETCH = ensembl_species.ensembl_fetch

    def __init__(self, base, convertonly, verbose):
        fetch.Fetcher.__init__(base, verbose)
        self.mainlog.debug('Creating new instance of EnsemblPlants')
        #   If we are only converting, then we don't have to sign on
        if convertonly:
            self.session = None
        else:
            self.session = self.sign_on()
        self.urls = []
        self.cksums = []
        self.to_convert = []
        return

    def sign_on(self):
        """Initialize FTP connection to Ensembl Plants with ENSEMBL_URL
        and navigate to ENSEMBL_PLANT_BASE."""
        ftp_session = ftplib.FTP(self.ENSEMBL_URL)
        #   login() without any args, use anonymous FTP by default
        ftp_session.login()
        ftp_session.cwd(self.ENSEMBL_PLANT_BASE)
        return ftp_session

    def get_file(self, fname):
        """Download the file specified by `fname'"""
        handle = open(file_funcs.local_name(fname), 'wb')
        self.session.retrbinary('RETR ' + fname, handle.write)
        handle.close()
        return

    def get_ftp_urls(self):
        """Build a list of FTP URLs to fetch. Currently searches for all
        files listed in a CDS directory. Also reads the CRC sums out of
        the CHECKSUMS file in each FTP directory. Save the URL in the
        `urls' attribute, and the checksums in the `cksums' attribute."""
        #   nlst() returns a listing of the directory contents
        #   We are going to assume that each of these are directories, since
        #   that is the way Ensembl has their server set up...
        for d in self.session.nlst():
            #   If it's in our list of species to download...
            if d in self.ENSEMBL_TO_FETCH:
                self.mainlog.debug(
                    'Attempting to cd into ' +
                    self.ENSEMBL_PLANT_BASE +
                    d +
                    '/cds/')
                #   cd into the CDS directory there
                self.session.cwd(self.ENSEMBL_PLANT_BASE + d + '/cds/')
                #   Get the contents
                listing = self.session.nlst()
                #   Then find the one that ends in .fa.gz
                for l in listing:
                    if l.endswith('.fa.gz'):
                        #   Tack the full path onto the list of URLS
                        self.mainlog.debug(
                            'Appending ' +
                            self.session.pwd() +
                            '/' +
                            l +
                            ' onto list of files to fetch.')
                        self.urls.append(self.session.pwd() + '/' + l)
                    #   get the checksum
                    elif l == 'CHECKSUMS':
                        #   We create a new StringIO instance, which can
                        #   read/write string data like file data
                        c = StringIO()
                        self.session.retrbinary('RETR ' + l, c.write)
                        #   Then we get the contents of the CHECKSUMS file
                        #   Split on newlines, since there is one line per
                        #   file entry
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

    def download_files(self):
        """Iterate through the list of URLs and download the appropriate
        files. Computes the CRC sum of existing files and compares them to
        the remote checksum to decide whether or not to to download."""
        #   For each URL we have:
        for u, c in zip(self.urls, self.cksums):
            target_dir = self.make_species_dir(u)
            #   cd into it
            os.chdir(target_dir)
            #   What is the local file name?
            lname = file_funcs.local_name(u)
            #   If it exists, we check if the checksums are the same
            if file_funcs.file_exists(lname, self.mainlog):
                local_cksum = file_funcs.calculate_crc32(lname, self.mainlog)
                crc32_same = file_funcs.checksum_is_same(
                    local_cksum,
                    c,
                    self.mainlog)
                if crc32_same:
                    self.mainlog.info(
                        lname +
                        ' already exists and is current, skipping.')
                    continue
                else:
                    self.mainlog.info(
                        lname +
                        ' exists, but is out of date. Updating.')
                    same = False
                    while not same:
                        self.get_file(u)
                        new_local_cksum = file_funcs.calculate_crc32(
                            lname,
                            self.mainlog)
                        same = file_funcs.checksum_is_same(
                            new_local_cksum,
                            c,
                            self.mainlog)
                    #   And save a record for those that need to be converted
                    self.to_convert.append(
                        os.path.join(
                            self.base,
                            target_dir, lname))
            #   If the file doesn't exist, then it's the same
            #   as if the checksum were different
            else:
                self.mainlog.info(lname + ' does not exist. Downloading.')
                same = False
                while not same:
                    self.get_file(u)
                    new_local_cksum = file_funcs.calculate_crc32(
                        lname,
                        self.mainlog)
                    same = file_funcs.checksum_is_same(
                        new_local_cksum,
                        c,
                        self.mainlog)
                self.to_convert.append(
                    os.path.join(
                        self.base,
                        target_dir,
                        lname))
        self.mainlog.info('Done downloading CDS files from Ensembl.')
        #   We are done with the FTP connection, log out
        self.session.quit()
        return
