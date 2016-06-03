#!/usr/bin/env python
"""A script to perform BLAST searching and sequence fetching."""

#   Import standard library modules here
import tempfile
import os
import re

#   Import the Biopython library
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

#   Import the script to give verbose messages
from lrt_predict.General import set_verbosity
#   Operate on files
from lrt_predict.General import file_funcs
#   For checking for the presence of executables
from lrt_predict.General import check_modules
#   For fetching sequences from the BLAST databases
from lrt_predict.Blast import sequence_fetch


#   A class to handle our BLAST searches
#   Borrowing heavily from script written by Paul Hoffman
class BlastSearch(object):
    """A class to handle BLAST searching and fetching of homologous sequences
    from plant CDS databases, fetched from Phytozome 10 and Ensembl Plants.

    The only functions that are designed to be called from outside the class
    are blast_all() and get_hit_seqs(). The others support the main two.

    best_hit():
        Iterates through a BlastRecord object, and returns the best (lowest
        E-value) hit.

    gen_output():
        Creates a temporary file for storing BLAST output. Returns a filename

    run_blast():
        Runs tblastx with a query sequence on a databse. Uses the filename
        from gen_output() to save the XML report. Returns the filename with
        the NCBI XML BLAST output.

    get_seq_id():
        Read through the BLAST XML and find the BlastRecord that has the
        lowest E-value. Feed that BlastRecord to best_hit(), and return the
        sequence ID of the hit.

    blast_all():
        BLAST search the provided query sequence against each species databse.
        Return the list of best BLAST hits.

    get_hit_seqs():
        Using the output from blast_all(), get the FASTA sequence of each of
        the homologous sequences and write them into a temporary file.
    """

    def __init__(self, base, target, query, evalue, verbose):
        """Initialize the class with base directory, query sequence, e-value
        threshold and verbosity level."""
        self.query = query
        self.evalue = evalue
        self.orthologues = {}
        self.mainlog = set_verbosity.verbosity('BLAST_Search', verbose)
        self.basedir = base
        self.target = target
        return

    def best_hit(self, brecord):
        """Define a special function to get the best hit out of a BlastRecord
        We need to do this because of the ugly nested nature of the data
        structure 'return' immediately stops iteration, whereas 'break'
        only works on one loop."""
        for aln in brecord.alignments:
            for hsp in aln.hsps:
                frames = [str(f) for f in hsp.frame]
                debug_msg = aln.title + ' Stats:\n' +\
                    'Bit Score: ' + str(hsp.bits) + '\n'\
                    'E-value: ' + str(hsp.expect) + '\n'\
                    'Identities ' + str(hsp.identities) + '\n'\
                    'Aln. Length: ' + str(hsp.align_length) + '\n'\
                    'Hit Start: ' + str(hsp.sbjct_start) + '\n'\
                    'Hit End: ' + str(hsp.sbjct_end) + '\n'\
                    'Frames: ' + ', '.join(frames)
                self.mainlog.debug(debug_msg)
                if hsp.expect <= self.evalue:
                    return aln.title
        else:
            return None

    def gen_output(self):
        """Define a function to create and manage output paths."""
        self.mainlog.debug('Creating named temporary file for BLAST output.')
        #   We use a python wrapper to mktemp - it automatically removes the
        #   temporary file when the handle is closed. We will just use it to
        #   store the BLAST output
        #   We use mode=w+t since we want read/write in text mode.
        #   We also need the filename, so used the NamedTemporaryFile method
        temp_output = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_BlastSearch_',
            suffix='_BLASTout.xml')
        self.mainlog.debug('Temp file created with name ' + temp_output.name)
        #   Return the file-like object
        return temp_output

    def run_blast(self, database):
        """Define a function to run the BLAST command."""
        #   Create a temp file
        blastout = self.gen_output()
        #   Start building a command line
        cline = NcbitblastxCommandline(
            query=self.query,
            out=blastout.name,
            db=database,
            evalue=self.evalue,
            outfmt=5,
            max_target_seqs=5)
        self.mainlog.debug(str(cline))
        #   And then execute it
        cline()
        return blastout

    def get_seq_id(self, out):
        """Define a function to get the best hit out of a BLAST record."""
        #   Seek to the beginning of the temporary file to read its contents
        out.seek(0)
        #   Then read the xml out
        self.mainlog.debug('Beginning to parse the BLAST output for ' +
                           self.query + '.')
        blast_records = NCBIXML.parse(out)
        #   Convert it to a list, since it should be relatively small...
        blast_records = list(blast_records)
        self.mainlog.info('We found ' + str(len(blast_records)) +
                          ' hits for ' + self.query + '.')
        #   If the list is empty, then we have no hits
        if len(blast_records) > 0:
            #   For each record
            for rec in blast_records:
                best = self.best_hit(rec)
                if best:
                    self.mainlog.info('Saving ' + best + ' as best hit.')
                    break
        #   Close the temporary file to clean up
        #   it's automatically deleted
        out.close()
        return best

    def blast_all(self):
        """Define a function to BLAST against every database."""
        databases = file_funcs.get_file_by_ext(self.basedir,
                                               '.fa',
                                               self.mainlog)
        self.mainlog.info('Running BLAST on ' +
                          str(len(databases)) +
                          ' species databases.')
        self.mainlog.debug('BLAST databases:\n' + '\n'.join(databases))
        #   databases will always be at least length 1
        #   if there is nothing, then the list will just have the empty string
        if databases[0] == '':
            self.mainlog.error('The base directory ' +
                               self.basedir +
                               ' does not contain any BLAST databases!')
            exit(1)
        for blast_db in databases:
            #   If the target species is in the filename of the FASTA sequence,
            #   we will skip it.
            if self.target.upper() in blast_db.upper():
                continue
            #   Get the BLAST output
            blast_output = self.run_blast(blast_db)
            #   And parse it
            homologous_locus = self.get_seq_id(blast_output)
            #   We do this check in case there is no match in a species
            #   Only save those that have a match
            if homologous_locus:
                #   We want the first and second parts, separated by a space
                fasta_info = homologous_locus.split(' ')
                seq_id = fasta_info[0]
                #   We need this part if we want to search by regex
                gb_id = fasta_info[1]
                #   And then tack it onto the list of orthologues
                self.orthologues[blast_db] = (seq_id, gb_id)
        return

    def get_hit_seqs(self):
        """Define a function to get the hit sequences out of the databses."""
        #   Create a temporary file for holding sequence information while we
        #   collect it
        if len(self.orthologues) == 0:
            self.mainlog.critical(
                'Could not find any BLAST hits! Try raising the E-value '
                'threshold for homology.')
            exit(2)
        self.mainlog.debug('Creating named tempfile for homologous sequences.')
        temp_output = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='BAD_Mutations_BlastSearch_',
            suffix='_orthologues.fasta')
        self.mainlog.debug('Created temporary file ' + temp_output.name)
        qseq = SeqIO.read(self.query, 'fasta')
        #    Start a new string to write the data into the file
        towrite = '>' + qseq.name + '\n' + str(qseq.seq) + '\n'
        #   Check to see if the blastdbcmd command is avilable
        blastdbcmd_path = check_modules.check_executable('blastdbcmd')
        if blastdbcmd_path:
            self.mainlog.debug('Using ' + blastdbcmd_path)
            for database, seqid in self.orthologues.iteritems():
                fasta, error = sequence_fetch.blastdbcmd(
                    blastdbcmd_path,
                    database,
                    seqid[0])
                self.mainlog.debug('Stdout:\n' + fasta)
                self.mainlog.debug('Stderr:\n' + error)
                #   If the sequence has ambiguous nucleotides (WRKYSMVBDHN),
                #   then we exclude it.
                db_seq = ''.join(fasta.split('\n')[1:])
                #   Then, search for the ambiguous nucleotides, skip if they are
                #   found.
                if re.search('S|W|R|K|Y|M|V|B|D|H|N', db_seq, re.I):
                    self.mainlog.warning(
                        'Removing sequence from ' +
                        os.path.basename(database) +
                        ' due to ambiguous nucleotides.')
                    continue
                #   We will use the name of the assembly as the species name
                spname = os.path.basename(database)
                #   Then split on . and take the first part
                spname = '>' + spname.split('.')[0]
                #   Then replace the weird ID with the species name
                fasta = re.sub('>.+', spname, fasta)
                towrite += fasta
        else:
            self.mainlog.debug('Using regex')
            for database, seqid in self.orthologues.iteritems():
                fasta = sequence_fetch.get_seq_by_regex(database, seqid[1])
                #   Perform the same check here for ambiguous nucleotides.
                db_seq = ''.join(fasta.split('\n')[1:])
                #   Then, search for the ambiguous nucleotides, skip if they are
                #   found.
                if re.search('S|W|R|K|Y|M|V|B|D|H|N', db_seq, re.I):
                    self.mainlog.warning(
                        'Removing sequence from ' +
                        os.path.basename(database) +
                        ' due to ambiguous nucleotides.')
                #   We will use the name of the assembly as the species name
                spname = os.path.basename(database)
                #   Then split on . and take the first part
                spname = '>' + spname.split('.')[0]
                #   Then replace the weird ID with the species name
                fasta = re.sub('>.+', spname, fasta)
                towrite += fasta
        self.mainlog.debug('Writing sequences into ' + temp_output.name)
        temp_output.write(towrite)
        #   We flush() it so that there is no data left unwritten
        temp_output.flush()
        #   Just as a test, try spitting the data back out
        temp_output.seek(0)
        self.mainlog.debug(temp_output.read())
        return temp_output
