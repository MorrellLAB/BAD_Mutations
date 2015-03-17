#!/usr/bin/env python

#   A script that contains functions to operate on files/paths

import os
import hashlib
import subprocess
import logging


#   A simple wrapper around os.path.isfile(), since this name is easier to read
def file_exists(fname, l):
    l.debug('Checking if ' + fname + ' exists.')
    if os.path.isfile(fname):
        return True
    else:
        return False


#   A function that generates a local filename from a remote one
def local_name(url):
    #   We assume that the remote url is separated by / and that the file name
    #   is the final part of the url
    return url.split('/')[-1]


#   A function that finds the species name given the filename
def species_name(fname):
    #   Again, assume this time that the species name is the first field in
    #   the _-delimited filename
    return fname.split('_')[0]


#   A function to calculate md5
def calculate_md5(fname, l, blocksize=8192):
    l.debug('Calculating MD5 on ' + fname + ' with a blocksize of ' + str(blocksize))
    #   If the supplied block size isn't a multiple of 128, then we will exit
    #   with an error, since md5 has a digest block size of 128 bytes
    if (blocksize % 128) != 0:
        l.error('Block size should be an integer multiple of 128!')
        exit(1)
    #   Start the md5 digest object
    md5 = hashlib.md5()
    #   Next we read in the file, blocksize bytes at a time, and feed it to md5
    #   we have to use the binary read option here
    with open(fname, 'rb') as f:
        while True:
            chunk = f.read(blocksize)
            #   If this is empty, then we have read the end of the file
            if not chunk:
                break
            #   Tack the chunk onto the md5 object
            md5.update(chunk)
    #   Then, return the hexadecimal hash (this is the one that is transmitted
    #   through plain text)
    return md5.hexdigest()


#   A function to check if the md5 of our local file is the same as that of the
#   remote file
def md5_is_same(local_md5, remote_md5, l):
    l.debug('Comparing local MD5 of ' + local_md5 + ' and remote MD5 of ' + remote_md5)
    if local_md5 == remote_md5:
        return True
    else:
        return False


#   A function to find all files with a given suffix in a certain directory
#   This is essentially a wrapper around the unix find command. This probably
#   will not work on windows, but this shouldn't be running on Windows 
#   anyway......
def get_cds_files(basedir, l):
    l.debug('Finding all *.cds.fa.gz files in ' + basedir)
    #   cd into the base directory
    os.chdir(basedir)
    #   Then execute the find command to get all gzipped cds files
    #   Normally, shell=True is a security hazard, but since we aren't actually
    #   running user-fed commands, we should be okay here
    raw_files = subprocess.check_output('find . -name "*.cds.fa.gz"', shell=True)
    file_list = raw_files.strip().split('\n')
    l.debug('Found ' + str(len(file_list)) + ' .cds.fa.gz files in ' + basedir)
    return file_list
