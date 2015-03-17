#!/usr/bin/env python

import logging

#   A script to set the verbosity of the program
def verbosity(name, v):
    l = logging.getLogger(name)
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
    return l
