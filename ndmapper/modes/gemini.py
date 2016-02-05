# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from ndmapper import config

# Default cal dependencies are defined by instrument mode, with just an empty
# placeholder dict at this level:
CAL_DEPS = {'target' : []}


def gemini_iraf_helper():
    """
    Define a few common default values for Gemini IRAF tasks in one place
    (just to avoid more repetition of static values than necessary) and
    convert some Python defaults to be more appropriate for the IRAF tasks.
    """

    # Enable both VAR & DQ if either is specified, since most Gemini tasks
    # can't control them individually:
    if config['use_uncert'] or config['use_flags']:
        vardq = True
    else:
        vardq = False

    vals = {'gmosdata' : 'gmos$data/',
            'key_airmass' : 'AIRMASS',
            'key_biassec' : 'BIASSEC',
            'key_datasec' : 'DATASEC',
            'key_exptime' : 'EXPTIME',
            'key_gain' : 'GAIN',
            'key_ron' : 'RDNOISE',
            'observatory' : 'default',
            'vardq' : vardq,
            'verbose' : True
           }

    return vals

