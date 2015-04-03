#!/usr/bin/env python

import logging

#   A script to set the verbosity of the program
def verbosity(name, v):
    l = logging.getLogger(name)
    s = logging.StreamHandler()
    formatter = logging.Formatter('===%(asctime)s - %(name)s===\n%(levelname)s\t%(message)s')
    l.setLevel(v)
    s.setLevel(v)
    #   Get the logger ready to do stuff
    s.setFormatter(formatter)
    l.addHandler(s)
    return l
