#!/usr/bin/env python
#   A script to setup the LRT prediction environment
#   Will create and parse a configuration file

#   Import the file functions script
from ..General import file_funcs
#   Import the script to give verbose messages
from ..General import set_verbosity

class SetupEnv:
    #   These are all the variables that we will write into the config file
    def __init__(self, verbose):
        self.mainlog = set_verbosity.verbosity('Setup_Env', verbose)
        return

    #   A function to set the variables and write the config file
    def set_vars(self, base, user, password, evalue, 
                 model, missingness, cfg):
        self.base = base
        self.username = user
        self.password = password
        self.eval_thresh = evalue
        self.alignment_model = model
        self.miss_thresh = missingness
        self.config_file = cfg
        self.mainlog.debug('Setting variables: \n' +\
            'BASE=' + self.base +\
            'USERNAME=' + self.username +\
            'PASSWORD=' + self.password +\
            'EVAL_THRESHOLD=' + self.evalue +\
            'MISS_THRESHOLD=' + self.miss_thresh +\
            'PRANK_MODEL=' + self.alignment_model)
        return
