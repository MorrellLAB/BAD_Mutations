#!/usr/bin/env python

#   A script that contains functions to operate on files/paths

import os
import hashlib
import subprocess
import logging
import re


def file_exists(fname, l):
    """A simple wrapper around os.path.isfile()."""
    l.debug('Checking if ' + fname + ' exists.')
    if os.path.isfile(fname):
        return True
    else:
        return False


def local_name(url):
    """Generate a local filename from a remote URL."""
    #   We assume that the remote url is separated by / and that the file name
    #   is the final part of the url
    return url.split('/')[-1]


def species_name(fname):
    """Find a species name given a filename."""
    #   Again, assume this time that the species name is the first field in
    #   the _-delimited filename
    species = fname.split('_')[0]
    #   Dots and underscores mess it up
    species = species.replace('.', '_')
    return species


def ensembl_species_name(fname):
    """Find a species name given an Ensembl filename."""
    #   The argument to this one is a whole path on the FTP server, so we have
    #   to split on / first
    lfile = local_name(fname)
    #   Ensembl filenames are like
    #   Genus_species.blahlblah.blahblah ... .fa.gz
    #   We want to match the genus_species part
    binomial = re.match('^[a-zA-Z]+(_[a-zA-Z]+)?', lfile)
    return binomial.group()


def calculate_md5(fname, l, blocksize=8192):
    """Calculate the MD5 sum."""
    l.info(
        'Calculating MD5 on ' + fname + ' with blocksize of ' + str(blocksize)
        )
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


def calculate_crc32(fname, l):
    """Calculate the CRC32 sum, as implemented in UNIX sum command."""
    l.info('Calculating CRC32 on ' + fname + ' with  `sum\' command.')
    cmd = ['sum', fname]
    p = subprocess.Popen(
        cmd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = p.communicate()
    l.debug('stdout:\n' + out.decode('utf-8'))
    l.debug('stderr:\n' + err.decode('utf-8'))
    return int(out.strip().split()[0])


def checksum_is_same(local_sum, remote_sum, l):
    """Compare the local and remote checksums."""
    l.info(
        'Comparing local checksum of ' +
        str(local_sum) +
        ' and remote checksum of ' +
        str(remote_sum))
    if local_sum == remote_sum:
        return True
    else:
        return False


def get_file_by_ext(basedir, suffix, l):
    """Find all files with a given suffix."""
    l.info('Finding all files with suffix ' + suffix + ' in ' + basedir)
    #   Then execute the find command to get all gzipped cds files
    #   Normally, shell=True is a security hazard, but since we aren't actually
    #   running user-fed commands, we should be okay here
    #   Example, with basedir=/tmp and suffix=.fasta we have
    #       find /tmp -name "*.fasta"
    command = ' '.join(['find', basedir, '-name', '"*'+suffix+'"'])
    raw_files = subprocess.check_output(command, shell=True).decode('utf-8')
    file_list = raw_files.strip().split('\n')
    l.debug(
        'Found ' +
        str(len(file_list)) +
        ' ' + suffix +
        ' files in ' +
        basedir)
    return file_list
