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
    #   The script is written in shell, so this function just calls it and
    #   checks the output
    #   Build the shell command. Apparently, we do no have to escape the pipe characters,
    #   the python warpper must handle that
    cmd = ['bash', SEQ_FETCH_SCRIPT, path, db, seqID]
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
#   Note that this solution is somewhat expensive, since it reads up
#   the entire sequence into memory
def get_seq_by_regex(db, seqID):
    #   Open a handle to the fasta file
    handle = open(db, 'r')
    #   And read its contents into memory.
    fasta = handle.read()
    #   Then we build the regex out of the sequence id
    #   We have to prepend a >, since that starts off the FASTA record
    s = '>' + seqID
    #   Then, we have to escape any dots in the name, since . is a special
    #   regex character
    s = s.replace('.', r'\.')
    #   Then build the full regex
    #   Match the start of the FASTA record, followed by any number of
    #   characters, followed by a newline, then ATCGN, across multiple lines
    #   Thanks to Paul Hoffman for the regex.
    pattern = re.compile('(^' + s + '.*\n[ATCGN\n]*)', re.M)
    #   Search the fasta for the right sequence and return it
    target_seq = pattern.search(fasta)
    #   Close this file handle
    handle.close()
    return fasta
