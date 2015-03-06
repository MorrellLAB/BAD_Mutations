#!/usr/bin/env python

#   A script to do some file handling, like making directories and moving files

import os

#   A function to create a base directory for the species
def makebase(basedir):
    #   Get the dirname of the LRT base
    base_parent = os.path.dirname(basedir)
    #   Make sure the directory is readable/writeable
    if not (os.access(base_parent, os.W_OK) and os.access(base_parent, os.R_OK)):
        print 'Error! Your LRT base directory ' + basedir + ' is not writable'
        exit(1)
    #    If the directory exists, then just skip it
    if os.path.exists(basedir):
        return
    else:
        print 'Creating base ' + basedir + ' ...'
        try:
            os.mkdir(basedir)
        except (OSError, IOError) as msg:
            print 'Could not create base', msg
    return


#   A general function to make a directory
def make_species_dir(basedir, spdir):
    #   Again make sure the directory is readable/writeable
    if not(os.access(basedir, os.W_OK) and os.access(basedir, os.R_OK)):
        print 'Error! Your LRT base directory ' + basedir + ' is not writable'
        exit(1)
    #   create the name of the target directory
    target_dir =os.path.join(basedir, spdir)
    #   Now do the checks
    if os.path.exists(target_dir):
        return target_dir
    else:
        print 'Creating ' + target_dir + ' ...'
        try:
            os.mkdir(target_dir)
        except (OSError, IOError) as msg:
            print 'Cound not create directory', msg
    return target_dir


#   A function for moving files
def move_file(filename, dirname):
    #   Check the target is read/write
    if not(os.access(dirname, os.W_OK) and os.access(dirname, os.R_OK)):
        print 'Error! ' + dirname + ' is not read/write'
        exit(1)
    #   Check that the file is read/write
    if not(os.access(filename, os.W_OK) and os.access(filename, os.R_OK)):
        print 'Error! ' + filename + ' is not read/write'
        exit(1)
    try:
        os.rename(filename, dirname + '/' + filename)
    except(OSError, IOError) as msg:
        print 'Could not move', filename, 'into', dirname, msg
    return
