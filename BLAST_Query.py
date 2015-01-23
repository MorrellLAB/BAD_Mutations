#!/usr/bin/env python
#   Takes a query sequence and runs BLAST against each of the species 
#   databases and returns a list of hits

#   Import the species databases file
import Species_Databases
import sys, os
try:
    #   Import our necessary Biopython libraries
    #       For running TBLASTX
    from Bio.Blast.Applications import NcbitblastxCommandline
    #       For parsing XML BLAST output
    from Bio.Blast import NCBIXML
    #       For handling/reading sequences
    from Bio import SeqIO
except ImportError:
    print "This script requires the Biopython libraries!"
    exit(1)

#   Define a function to take XML with BLAST results and return the best hit
def get_best_hit(results):
    pass


#   Define a function to run BLAST against a database
def run_blast(query, db):
    #   Takes a path to a BLAST database and a filename that contains just
    #   one sequence for a query
    #   Should add a way to check the number of sequences, as passing a
    #   multi-FASTA file as input probably breaks the script
    
    #   Check if the database exists and is readable
    if not os.access(db, os.R_OK):
        print "Warning! Database " + db + " is not readable or does not exist!"
        return(None)
