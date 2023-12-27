#!/usr/bin/env python
"""A class that defines a general fetching class. Phytozome and Ensembl
classes will inherit from this one."""

#   Import our helper scripts
from lrt_predict.General import check_modules
from lrt_predict.General import file_funcs
from lrt_predict.General import set_verbosity
from lrt_predict.Fetch import format_blast
from lrt_predict.General import dir_funcs


class Fetcher(object):
    """A general class that contains methods and variables that are common to
    all of the fetching that we have to do."""

    def __init__(self, base, verbose):
        self.base = base
        self.to_convert = []
        self.mainlog = set_verbosity.verbosity(__name__, verbose)
        self.mainlog.debug('Creating new instance of Fetcher')
        dir_funcs.makebase(base, self.mainlog)
        return

    def make_species_dir(self, fname):
        """Creates a species directory based on a file name, and return the
        name of the created directory."""
        spname = file_funcs.ensembl_species_name(fname)
        #   Then we have to make that directory
        target_dir = dir_funcs.make_species_dir(
            self.base,
            spname,
            self.mainlog)
        return target_dir

    def convert(self):
        """Iterates through the `to_convert' attribute and converts each file
            from a gzipped FASTA file to a BLAST database."""
        #   What is the path to the makeblastdb executable?
        makeblastdb_path = check_modules.check_executable('makeblastdb')
        #   Check if the list of updated CDS files is empty or not
        if not self.to_convert:
            #   If it is empty, then populate it with all of them
            fname_list = file_funcs.get_file_by_ext(
                self.base,
                '.fa.gz',
                self.mainlog)
        else:
            fname_list = self.to_convert
        #   for each one
        for fname in fname_list:
            out, error = format_blast.format_blast(makeblastdb_path, fname)
            self.mainlog.info('stdout: \n' + out.decode('utf-8'))
            self.mainlog.info('stderr: \n' + error.decode('utf-8'))
        return
