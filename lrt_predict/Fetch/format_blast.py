#!/usr/bin/env python

#    A script to manage the creation of BLAST databases

#   To manage subprocesses
import subprocess
#   To handle paths
import os

#   create a variable to hold the path to our installation directory
#   It is two levels above this one
#   os.path.realpath(__file__) is the full path to this script, blast_databases.py
#   rsplit() splits from the right end of the file, and the number argumetn
#   tells how many times to split. We take the first element of the resulting list.
LRT_PATH = os.path.realpath(__file__).rsplit(os.path.sep, 3)[0]

#   The directory that contains shell scripts
DBFORMAT_SCRIPT = os.path.join(LRT_PATH, 'Shell_Scripts', 'Unzip_CDS.sh')

#   Function that just calls the shell script that handles the BLAST database
#   formatting
def format_blast(makeblastdb_path, fname):
    #   The script is written in shell, so this function just calls it and
    #   checks the output
    #   Build the shell command
    cmd = ['bash', DBFORMAT_SCRIPT, makeblastdb_path, fname]
    #   Execute the script
    #   shell=False to ensure that we aren't executing commands from untrusted
    #   sources
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)
