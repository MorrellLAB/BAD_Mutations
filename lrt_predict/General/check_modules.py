#!/usr/bin/env python

#   A script to check for the presence of all the required modules

#   Gain access to the spawn.find_executable(), which works a lot like `which`
from distutils import spawn


def check_modules(setup=False, fetch=False, predict=False):
    """Function to try to import required modules, and return a list of modules
    that are not installed."""
    missing_modules = []
    #   Argparse
    try:
        import argparse
    except ImportError:
        missing_modules.append('argparse')
    #   If predict
    if predict:
        #   Biopython
        try:
            import Bio
        except ImportError:
            missing_modules.append('Biopython')
    return missing_modules


def missing_mods(modules):
    """Function to print a nice message about modules that are required."""
    msg = '''Some of the required modules were not found on your system. Please
install them and try again. The following modules were not found:'''
    print msg
    print '\n'.join(modules)
    return


def check_executable(exe):
    """Checks for a path to an executable."""
    path = spawn.find_executable(exe)
    if path:
        return path
    else:
        return False


def missing_executables(exelist):
    """Checks for the presence and execute permissions for all program names
    passed to it."""
    missing_programs = []
    for e in exelist:
        if check_executable(e):
            continue
        else:
            missing_programs.append(e)
    return missing_programs
