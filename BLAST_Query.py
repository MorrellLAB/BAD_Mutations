#!/usr/bin/env python
#   Takes a query sequence and runs BLAST against each of the species BLAST 
#   databases and returns a list of hits

import sys
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
