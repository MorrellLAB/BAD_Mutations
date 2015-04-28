#!/usr/bin/env python

#   A script that contains functions to operate on directories

import os
import logging

#   A function to set the verbosity level
def verbose(v):
    #   Initialize the logger, with the current module name
    l = logging.getLogger(__name__)
    s = logging.StreamHandler()
    formatter = logging.Formatter('===%(asctime)s - %(name)s===\n%(levelname)s\t%(message)s')
    if v:
        #   If verbose, set the logger to DEBUG to emit all sorts of msgs
        l.setLevel(logging.DEBUG)
        s.setLevel(logging.DEBUG)
    else:
        #   Else, we set it to warning
        l.setLevel(logging.WARNING)
        s.setLevel(logging.WARNING)
    #   Get the logger ready to do stuff
    s.setFormatter(formatter)
    l.addHandler(s)
    #   and return it so we can throw messages
    return l


#   A function to create a base directory for the species
def makebase(basedir, l):
    l.debug('Base directory is ' + basedir)
    #   Get the dirname of the LRT base
    base_parent = os.path.dirname(basedir)
    #    If the directory exists, then just skip it
    if os.path.exists(basedir):
        return
    else:
        l.info('Creating base ' + basedir + ' ...')
        try:
            os.mkdir(basedir)
        except (OSError, IOError) as msg:
            l.error('Could not create base directory' + msg)
    return


#   A general function to make a directory
def make_species_dir(basedir, spdir, l):
    l.debug('Species directory is ' + os.path.join(basedir, spdir))
    #   create the name of the target directory
    target_dir = os.path.join(basedir, spdir)
    #   Now do the checks
    if os.path.exists(target_dir):
        return target_dir
    else:
        l.info('Creating ' + target_dir + ' ...')
        try:
            os.mkdir(target_dir)
        except (OSError, IOError) as msg:
            l.error('Could not create directory' + msg)
    return target_dir


#   A function for moving files
def move_file(filename, dirname, l):
    #   Check that the file is read/write
    if not(os.access(filename, os.W_OK) and os.access(filename, os.R_OK)):
        l.error('Error! ' + filename + ' is not read/write')
        exit(1)
    try:
        l.debug('Moving ' + filename + ' to ' + os.path.join(dirname, filename))
        os.rename(filename, os.path.join(dirname, filename))
    except(OSError, IOError) as msg:
        l.error('Could not move ' + filename + ' into ' + dirname + '. ' + msg)
    return
