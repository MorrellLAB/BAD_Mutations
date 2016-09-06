#!/usr/bin/env python
#   A script to parse the configuration file

#   Import the verbosity script
from ..General import set_verbosity


class ConfigHandler(object):
    """A class to handle the reading and storing of configuration variables.
    Checks the validity of the configuration file as well as writes new
    configuration files."""
    #   Here are the keywords that we are accepting
    #   We define it as a dictionary with the config keywords as keys
    #   and the internal variable names as values.
    #   This is so that we can easily assign the config names to real
    #   variables that the program uses.
    KEYWORDS = {'BASE': 'base',
                'EVAL_THRESHOLD': 'evalue',
                'TARGET_SPECIES': 'target',
                'BASH': 'bash_path',
                'GZIP': 'gzip_path',
                'SUM': 'sum_path',
                'TBLASTX': 'tblastx_path',
                'PASTA': 'pasta_path',
                'HYPHY': 'hyphy_path'
                }
    #   Here is the string that prefixes a variable delcaration
    DECLR = '#define'

    def __init__(self, cfg, arg, verbose):
        self.config_file = cfg
        self.user_args = arg
        self.mainlog = set_verbosity.verbosity('Configuration_Handler', verbose)
        return

    def is_valid(self):
        """Checks if the configuration file is valid."""
        with open(self.config_file, 'r') as f:
            for index, line in enumerate(f):
                #   Check the variable declarations
                if line.startswith(self.DECLR):
                    #   If None is passed to str.split(), it will just split
                    #   on any whitespace. Pass 2, so that we split the string
                    #   into three elements:
                    #       #define KEYWORD VALUE
                    #   This allows VALUE to have spaces in it
                    #   (though not recommended)
                    tmp = line.strip().split(None, 2)
                    #   If we don't have three elements, the line is malformed
                    if len(tmp) != 3:
                        self.mainlog.error(
                            'Line ' +
                            str(index+1) +
                            ': Expected three fields, got ' +
                            str(len(tmp)))
                        return False
                    #   Then we get the keyword
                    k = tmp[1]
                    #   Check if the keyword is in the accepted keywords
                    if k not in self.KEYWORDS:
                        self.mainlog.warning(
                            'Line ' +
                            str(index+1) +
                            ': Unknown variable ' +
                            k)
                else:
                    continue
        return True

    def read_vars(self):
        """Step through the configuration file, and set variables as they are
        encountered."""
        conf_dict = {}
        with open(self.config_file, 'r') as f:
            #   Split up the lines as in the above function
            for line in f:
                if line.startswith(self.DECLR):
                    tmp = line.strip().split(None, 2)
                    k = tmp[1]
                    value = tmp[2]
                    #   We use a dictionary in lieu of a case/switch statement
                    if k in self.KEYWORDS:
                        self.mainlog.debug(
                            'Setting variable ' +
                            k +
                            ' to ' +
                            value)
                        conf_dict[self.KEYWORDS[k]] = value
                    else:
                        self.mainlog.warning('Unknown variable ' + k)
        self.config_vars = conf_dict
        return

    def merge_options(self):
        """Merge the two input dictionaries: the one recieved on the command
        line, and the one read out of the configuration file. Options given on
        the command line will take precedence over the ones in the
        configuration file."""
        #   List of keys that we do not want to clobber if they are not
        #   passed on the command line
        do_not_clobber = ['base']
        #   This takes up a bit of extra memory, but it will preserve the two
        #   dictionaries as separate variables.
        configs = self.config_vars.copy()
        #   We will iterate through the user-supplied arguments dictionary,
        #   and update those that need to be updated
        for option, val in self.user_args.iteritems():
            if option not in configs:
                configs[option] = val
            #   Check that the value was specified on the command line
            elif not val:
                configs.update({option: self.config_vars[option]})
            elif option not in do_not_clobber:
                self.mainlog.debug('OPTION: ' + option + '; VALUE: ' + str(val))
                configs.update({option: val})
        self.mainlog.debug(
            'Command line and config options merged. Values: ' + str(configs))
        #   Return them
        return configs
