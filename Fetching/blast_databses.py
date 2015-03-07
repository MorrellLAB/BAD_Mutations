#!/usr/bin/env python

#    A script to manage the creation of BLAST databases

#   To manage subprocesses
import subprocess
#   To handle paths
import os

#   The directory that contains shell scripts
DBFORMAT_SCRIPT = os.path.join('.', 'Shell_Scripts', 'Unzip_CDS.sh')

#   Function that just calls the shell script that handles the BLAST database
#   formatting
def format_blast(basedir):
    #   The script is written in shell, so this function just calls it and
    #   checks the output
    #   Build the shell command
    cmd = [DBFORMAT_SCRIPT, basedir]
    #   Execute the script
    #   shell=False to ensure that we aren't executing commands from untrusted
    #   sources
    retval = subprocess.call(cmd, shell=False)
    return retval
