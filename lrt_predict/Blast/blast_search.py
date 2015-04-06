#!/usr/bin/env python

#   A script to perform BLAST searching


#   Import standard library modules here
import tempfile

#   Import the Biopython library
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

#   Import the script to give verbose messages
from ..General import set_verbosity
#   Check inputs and extract information from them
from ..General import parse_input
#   Operate on files
from ..General import file_funcs
#   For checking for the presence of executables
from ..General import check_modules
#   For fetching sequences from the BLAST databases
import sequence_fetch

#   A class to handle our BLAST searches
#   Borrowing heavily from script written by Paul Hoffman
class BlastSearch:
    #   Initialize the class
    def __init__(self, base, query, evalue, verbose):
        self.query = query
        self.evalue = evalue
        self.homologues = {}
        self.mainlog = set_verbosity.verbosity('BLAST_Search', verbose)
        self.basedir = base

    #   Define a special function to get the best hit out of a BlastRecord
    #   We need to do this because of the ugly nested nature of the data structure
    #   'return' immediately stops iteration, whereas 'break' only works on one loop
    def best_hit(self, br):
        for a in br.alignments:
            for hsp in a.hsps:
                if hsp.expect <= self.evalue:
                    return a.title
        else:
            return None

    #   Define a function to create and manage output paths
    def gen_output(self):
        self.mainlog.debug('Creating named temporary file for BLAST output.')
        #   We use a python wrapper to mktemp - it automatically removes the
        #   temporary file when the handle is closed. We will just use it to
        #   store the BLAST output
        #   We use mode=w+t since we want read/write in text mode.
        #   We also need the filename, so used the NamedTemporaryFile method
        temp_output = tempfile.NamedTemporaryFile(mode='w+t', prefix='LRTPredict_BlastSearch_', suffix='_BLASTout.xml')
        self.mainlog.debug('Temp file created with name ' + temp_output.name)
        #   Return the file-like object
        return temp_output

    #   Define a function to run the BLAST command
    def run_blast(self, db):
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
                best = self.best_hit(r)
                if best:
                    self.mainlog.info('Saving ' + best + ' as best hit.')
                    break
        #   Close the temporary file to clean up
        #   it's automatically deleted
        out.close()
        return best

    #   Define a function to BLAST against every database
    def blast_all(self):
        databases = file_funcs.get_file_by_ext(self.basedir, '.cds.fa', self.mainlog)
        self.mainlog.info('Running BLAST on ' + str(len(databases)) + ' species databases.')
        self.mainlog.debug('BLAST databases:\n' + '\n'.join(databases))
        #   databases will always be at least length 1
        #   if there is nothing, then the list will just have the empty string
        if databases[0] == '':
            self.mainlog.error('The base directory ' + self.basedir + ' does not contain any BLAST databases!')
            exit(1)
        for d in databases:
            #   Get the BLAST output
            blast_output = self.run_blast(d)
            #   And parse it
            homologous_locus = self.get_seq_ID(blast_output)
            #   We do this check in case there is no match in a species
            #   Only save those that have a match
            if homologous_locus:
                #   We want the first and second parts, separated by a space
                fasta_info = homologous_locus.split(' ')
                seq_id = fasta_info[0]
                #   We need this part if we want to search by regex
                gb_id = fasta_info[1]
                #   And then tack it onto the list of homologues
                self.homologues[d] = (seq_id, gb_id)
        return

    #   Define a function to get the hit sequences out of the databses
    def get_hit_seqs(self):
        #   Create a temporary file for holding sequence information while we collect it
        self.mainlog.debug('Creating named temporary file for homologous sequences.')
        temp_output = tempfile.NamedTemporaryFile(mode='w+t', prefix='LRTPredict_BlastSearch_', suffix='_homologues.fasta')
        self.mainlog.debug('Created temporary file ' + temp_output.name)
        qseq = SeqIO.read(self.query, 'fasta')
        #    Start a new string to write the data into the file
        towrite = '>' + qseq.name + '\n' + str(qseq.seq) + '\n'
        #   Check to see if the blastdbcmd command is avilable
        blastdbcmd_path = check_modules.check_executable('blastdbcmd')
        if blastdbcmd_path:
            self.mainlog.debug('Using ' + blastdbcmd_path)
            for database, seqID in self.homologues.iteritems():
                fasta, error = sequence_fetch.blastdbcmd(blastdbcmd_path, database, seqID[0])
                towrite += fasta
        else:
            self.mainlog.debug('Using regex')
            for database, seqID in self.homologues.iteritems():
                fasta = sequence_fetch.get_seq_by_regex(database, seqID[1])
                towrite += fasta
        self.mainlog.debug('Writing sequences into ' + temp_output.name)
        temp_output.write(towrite)
        #   We flush() it so that there is no data left unwritten
        temp_output.flush()
        return temp_output
