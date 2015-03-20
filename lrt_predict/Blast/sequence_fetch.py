#!/usr/bin/env python

#   A script to retrieve FASTA sequences from a BLAST database
#   Basically just a fancy wrapper around blastdbcmd

import subprocess
import re
import os

#   create a variable to hold the path to our installation directory
#   It is two levels above this one
#   os.path.realpath(__file__) is the full path to this script, sequence_fetch.py
#   rsplit() splits from the right end of the file, and the number argumetn
#   tells how many times to split. We take the first element of the resulting list.
LRT_PATH = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
#   The directory that contains shell scripts
SEQ_FETCH_SCRIPT = os.path.join(LRT_PATH, 'Shell_Scripts', 'Seq_From_BLASTdb.sh')

#   Easier way: just use the packaged blastdbcmd from NCBI
def blastdbcmd(path, db, seqID):
    #   The only thing is we have to be careful with pipes in the sequence ID
    #   escape them.
    esc_seqID = seqID.replace('|', '\|')
    #   The script is written in shell, so this function just calls it and
    #   checks the output
    #   Build the shell command
    cmd = ['sh', SEQ_FETCH_SCRIPT, path, db, seqID]
    #   Execute the script
    #   shell=False to ensure that we aren't executing commands from untrusted
    #   sources. We set out and err to subprocess.PIPE so we can save the actual
    #   output for later
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    #   return the output, which is just the fasta file.
    return (out, err)

#   Failing that, we can use a regular expression to get it
#   Regular expression written by Paul Hoffman.
def get_seq_by_regex(db, seqID):
    pass
