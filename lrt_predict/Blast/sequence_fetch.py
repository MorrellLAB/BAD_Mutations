#!/usr/bin/env python

#   A script to retrieve FASTA sequences from a BLAST database
#   Basically just a fancy wrapper around blastdbcmd

import subprocess
import re
import os

#   create a variable to hold the path to our installation directory
#   It is two levels above this one
LRT_PATH = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]
#   The directory that contains shell scripts
SEQ_FETCH_SCRIPT = os.path.join(
    LRT_PATH,
    'Shell_Scripts',
    'Seq_From_BLASTdb.sh')


def blastdbcmd(path, db, seqID):
    """Wrapper function for the blastdbcmd sequence fetch command."""
    cmd = ['bash', SEQ_FETCH_SCRIPT, path, db, seqID]
    #   Execute the script
    #   shell=False to ensure that we aren't executing commands from untrusted
    #   sources. We set out and err to subprocess.PIPE so we can save the
    #   output for later
    p = subprocess.Popen(
        cmd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = p.communicate()
    #   return the output, which is just the fasta file.
    return (out, err)


def get_seq_by_regex(db, seqID):
    """Use a regular expression to pull out FASTA sequences."""
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
    return target_seq
