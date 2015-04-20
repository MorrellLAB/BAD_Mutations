#!/usr/bin/env python

#   A script that contains functions to operate on files/paths

import os
import hashlib
import subprocess
import logging
import re

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
    species = fname.split('_')[0]
    #   Dots and underscores mess it up
    species = species.replace('.', '_')
    return species


#   One for ensembl, too
def ensembl_species_name(fname):
    #   The argument to this one is a whole path on the FTP server, so we have
    #   to split on / first
    lfile = local_name(fname)
    #   Ensembl filenames are like
    #   Genus_species.blahlblah.blahblah ... .fa.gz
    #   We want to match the genus_species part
    binomial = re.match('^[a-zA-Z]+_[a-zA-Z]+', lfile)
    return binomial.group()


#   A function to calculate md5
def calculate_md5(fname, l, blocksize=8192):
    l.info('Calculating MD5 on ' + fname + ' with a blocksize of ' + str(blocksize))
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


#   A function to calculate the CRC32 sum as implemented in the UNIX 'sum' cmd
def calculate_crc32(fname, l):
    l.info('Calculating CRC32 on ' + fname + ' with  `sum\' command.')
    cmd = ['sum', fname]
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    l.debug('stdout:\n' + out)
    l.debug('stderr:\n' + err)
    return int(out.strip().split()[0])


#   A function to check if the checksum of our local file is the same as that of the
#   remote file
def checksum_is_same(local_sum, remote_sum, l):
    l.info('Comparing local checksum of ' + str(local_sum) + ' and remote checksum of ' + str(remote_sum))
    if local_sum == remote_sum:
        return True
    else:
        return False


#   A function to find all files with a given suffix in a certain directory
#   This is essentially a wrapper around the unix find command. This probably
#   will not work on windows, but this shouldn't be running on Windows 
#   anyway......
def get_file_by_ext(basedir, suffix, l):
    l.info('Finding all files with suffix ' + suffix + ' in ' + basedir)
    #   Then execute the find command to get all gzipped cds files
    #   Normally, shell=True is a security hazard, but since we aren't actually
    #   running user-fed commands, we should be okay here
    #   Example, with basedir=/tmp and suffix=.fasta we have
    #       find /tmp -name "*.fasta"
    command = ' '.join(['find', basedir, '-name', '"*'+suffix+'"'])
    raw_files = subprocess.check_output(command, shell=True)
    file_list = raw_files.strip().split('\n')
    l.debug('Found ' + str(len(file_list)) +  ' ' + suffix + ' files in ' + basedir)
    return file_list
