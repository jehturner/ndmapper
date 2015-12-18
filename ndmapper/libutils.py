# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os

def splitext(path):
    """
    A version of splitext that splits at the first separator rather than the
    last one (so 'file.fits.gz' gives 'file' & 'fits.gz'). It also returns
    None for the extension value where there isn't one (instead of ''), just
    to avoid incorrect reconstruction of 'file.' as 'file' or vice versa.
    """

    components = path.split(os.extsep, 1) # always produces a root name

    return components[0], None if len(components) == 1 else components[1]


def addext(path, ext):
    """
    Reconstruct a filename from a (root, extension) tuple of the type
    produced by splitext().
    """
    return path + ('' if ext is None else os.extsep + ext)

