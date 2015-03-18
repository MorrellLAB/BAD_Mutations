#!/usr/bin/env python

#   A script to perform BLAST searching


#   Import standard library modules here
import tempfile

#   Import the Biopython library
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast import NCBIXML

#   Import the script to give verbose messages
from ..General import set_verbosity
#   Check inputs and extract information from them
from ..General import parse_input
#   Operate on files
from ..Fetch import file_funcs

#   A class to handle our BLAST searches
#   Borrowing heavily from script written by Paul Hoffman
class BlastSearch:
    #   Initialize the class
    def __init__(self, base, query, evalue, verbose):
        self.query = query
        self.evalue = evalue
        self.homologues = []
        self.mainlog = set_verbosity.verbosity('BLAST_Search', verbose)
        self.basedir = base

    #   Define a function to create and manage output paths
    def gen_output(self):
        self.mainlog.debug('Creating named temporary file for BLAST output.')
        #   We use a python wrapper to mktemp - it automatically removes the
        #   temporary file when the handle is closed. We will just use it to
        #   store the BLAST output
        #   We use mode=w+t since we want read/write in text mode.
        #   We also need the filename, so used the NamedTemporaryFile method
        temp_output = tempfile.NamedTemporaryFile(mode='w+t')
        self.mainlog.debug('Temp file created with name ' + temp_output.name)
        #   Return the file-like object
        return temp_output

    #   Define a function to run the BLAST command
    def run_blast(self, db):
        #   Check to see if the input file is okay
        if not parse_input.valid_fasta(self.query, self.mainlog):
            self.mainlog.error('The input FASTA file provided is not valid.')
            exit(1)
        #   Create a temp file
        blastout = self.gen_output()
        #   Start building a command line
        cline = NcbitblastxCommandline(
            query=self.query,
            out=blastout.name,
            db=db,
            evalue=self.evalue,
            outfmt=5,
            max_target_seqs=5)
        self.mainlog.debug(str(cline))
        #   And then execute it
        cline()
        return blastout

    #   Define a function to get the best hit out of a BLAST record
    def get_seq_ID(self, out):
        #   Seek to the beginning of the temporary file to read its contents
        out.seek(0)
        #   Then read the xml out
        self.mainlog.debug('Beginning to parse the BLAST output for ' + self.query + '.')
        blast_records = NCBIXML.parse(out)
        #   Convert it to a list, since it should be relatively small...
        blast_records = list(blast_records)
        self.mainlog.info('We found ' + str(len(blast_records)) + ' hits for ' + self.query + '.')
        #   If the list is empty, then we have no hits
        if len(blast_records) > 0:
            #   For each record
            for r in blast_records:
                #   For each alignment in the record
                for a in r.alignments:
                    #   And then for each high-scoring pair, check the e-value
                    for hsp in a.hsps:
                        #   Check the e-value
                        if hsp.expect <= self.evalue:
                            #   If it's less than or equal to our threshold, keep it
                            homologue = a.title
                            self.mainlog.info('Saving ' + a.title + ' as best hit.')
                            break
        #   Close the temporary file to clean up
        #   it's automatically deleted
        out.close()
        return homologue

    #   Define a function to BLAST against every database
    def blast_all(self):
        databases = file_funcs.get_file_by_ext(self.basedir, '.cds.fa', self.mainlog)
        self.mainlog.info('Running BLAST on ' + str(len(databases)) + ' species databases.')
        self.mainlog.debug('BLAST databases:\n' + '\n'.join(databases))
        for d in databases:
            #   Get the BLAST output
            blast_output = self.run_blast(d)
            #   And parse it
            homologous_locus = self.get_seq_ID(blast_output)
            #   And then tack it onto the list of homologues
            self.homologues.append(homologous_locus)
        return
