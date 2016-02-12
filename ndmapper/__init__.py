# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

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
          'labels' : {'data' : 'SCI', 'uncertainty' : 'VAR', 'flags' : 'DQ'},
          'reprocess' : None,
          'use_uncert' : True,
          'use_flags' : True,
          'interact' : False
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

    # Define a decorator that converts some package-API-standard arguments to
    # use package run-time default values, unless specified explicitly:
    from functools import wraps
    import inspect
    def ndprocess_defaults(proc_fn):

        @wraps(proc_fn)  # transfer target function name & docstring to wrapper

        # Wrapper function to be substituted for the original:
        def wrap_defaults(*args, **kwargs):

            # Evaluate what arguments the processing function we're wrapping
            # would receive if we called it directly:
            callargs = inspect.getcallargs(proc_fn, *args, **kwargs)

            # Convert any unspecified flags to the run-time defaults specified
            # in the package configuration dictionary:
            for param in ['interact', 'reprocess', 'use_uncert', 'use_flags']:
                if param in callargs and callargs[param] is None:
                    callargs[param] = config[param]

            # Call the processing function on our modified arguments:
            return proc_fn(**callargs)

        return wrap_defaults

