#!/usr/bin/env python

#   A script to perform some simple file checks

import os
import hashlib

#   A simple wrapper around os.path.isfile(), since this name is easier to read
def file_exists(fname):
    if os.path.isfile(fname):
        return True
    else:
        return False


#   A function to calculate md5
def calculate_md5(fname, blocksize=8192):
    #   If the supplied block size isn't a multiple of 128, then we will exit
    #   with an error, since md5 has a digest block size of 128 bytes
    if (blocksize % 128) != 0:
        print 'Block size should be an integer multiple of 128!'
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
def md5_is_same(local_md5, remote_md5):
    #   First, check if the local_md5 is a reference to a file or not
    #   If it is, read the contents
    if file_exists(local_md5 + '.md5'):
        handle = open(local_md5 + '.md5', 'r')
        #   This may break our program if someone has tampered with the
        #   md5 file. This shouldn't be the case, but it is possible.
        local_md5 = handle.read().strip()
        handle.close()
    else:
        pass
    if local_md5 == remote_md5:
        return True
    else:
        return False
