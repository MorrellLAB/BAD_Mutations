#!/usr/bin/env python

#   A script to check the validity of arguments
import os
import re

#   Does the directory exist, and is it readable and writable?
def valid_dir(dname):
    parent = os.path.dirname(dname)
    if os.path.isdir(dname):
        if (os.access(dname, os.R_OK) and os.access(dname, os.W_OK)):
            return True
        else:
            return False
    #   Check the parent directory
    elif os.path.isdir(parent):
        if (os.access(parent, os.R_OK) and os.access(dname, os.W_OK)):
            return True
        else:
            return False
    else:
        return False


#   Is the username an email (sorta?)
#       This is a weak regex to identify an email, there are better ones
#       but this will do.
def valid_email(email):
    #   This matches one or more non-whitespace, followed by @,
    #   followed by .-separated non-whitespace, and ending in non-whitespace
    email_regex = r'\S+@(\S+\.)+\S'
    if re.match(email_regex, email):
        return True
    else:
        return False
