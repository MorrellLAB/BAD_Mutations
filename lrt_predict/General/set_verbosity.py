#!/usr/bin/env python

import logging


def verbosity(name, v):
    """Set the verbosity."""
    l = logging.getLogger(name)
    s = logging.StreamHandler()
    formatter = logging.Formatter(
        '===%(asctime)s - %(name)s===\n%(levelname)s\t%(message)s')
    l.setLevel(v)
    s.setLevel(v)
    #   Get the logger ready to do stuff
    s.setFormatter(formatter)
    l.addHandler(s)
    return l
