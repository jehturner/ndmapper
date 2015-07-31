# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A simple package to automate pushing data through a series of processing steps.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

config = {'filename_regex' : '[NS][0-9]{8}S[0-9]{3,4}',
          'logfile' : None,
          'data_name' : 'SCI',
          'uncertainty_name' : 'VAR',
          'flags_name' : 'DQ'
         }

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:

    # from example_mod import *
    import re
    config['filename_regex'] = re.compile(config['filename_regex'])

    # Set PyRAF options in the environment before it gets imported.
    # This is to solve the plotting hang on Apple reported by Bryan M.
    import os
    os.environ['PYRAF_BETA_STATUS'] = '1'

