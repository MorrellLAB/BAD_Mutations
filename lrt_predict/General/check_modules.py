#!/usr/bin/env python

#   A script to check for the presence of all the required modules

#   Try to import them, return names of those that are not installed
def check_modules(fetch, predict):
    missing_modules = []
    #   Argparse
    try:
        import argparse
    except ImportError:
        missing_modules.append('argparse')
    #   If fetch:
    if fetch:
        #   Requests
        try:
            import requests
        except ImportError:
            missing_modules.append('requests')
    #   If predict
    if predict:
        #   Biopython
        try:
            import Bio
        except ImportError:
            missing_modules.append('Biopython')
    return missing_modules

#   A function to print a nice little message about missing modules
def missing_mods(modules):
    msg = '''Some of the required modules were not found on your system. Please
install them and try again. The following modules were not found:'''
    print msg
    print '\n'.join(modules)
    return

