# Copyright(c) 2015-2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

"""
Some lower-level utility functions (with minimal dependencies) that are used
internally and are also made available as part of the public API.
"""

import os, os.path
import tempfile


def splitext(path):
    """
    A version of splitext that splits at the first separator rather than the
    last one (so 'file.fits.gz' gives 'file' & 'fits.gz'). It also returns
    None for the extension value where there isn't one (instead of ''), just
    to avoid incorrect reconstruction of 'file.' as 'file' or vice versa.
    """

    components = os.path.basename(path).split(os.extsep, 1)  # len always 1->2
    ext = None if len(components) == 1 else components[1]
    root = path if ext is None else path[:-len(os.extsep+ext)]

    return root, ext


def addext(path, ext):
    """
    Reconstruct a filename from a (root, extension) tuple of the type
    produced by splitext().
    """
    return path + ('' if ext is None else os.extsep + ext)


def new_filename(purpose='tmp', base='', ext='', full_path=False):
    """
    Generate a new filename string that is not already used in the current
    directory (beginning with 'tmp' by default, for use as a temporary file).
    Unlike Python's tempfile module, this function does not actually open the
    file, making the result suitable for passing to external programs, but as
    a result, a race condition may occur if the file is not created
    immediately, which is the user's responsibility.

    Parameters
    ----------

    purpose : str, optional
        Starting string, used to indicate the file's purpose (default 'tmp').

    base : convertible to str, optional
        A base name to add between "tmp_" and the last few, randomized
        characters, to help distinguish temporary filenames, eg. for
        troubleshooting purposes.

    ext : convertible to str, optional
        An file extension name to use (eg 'fits'). The leading dot is optional
        and will be added if needed.

    full_path : bool
        Return the full path to the file, rather than a (relative) filename in
        the current working directory (default False)?


    Returns
    -------

    str
        A filename that doesn't already exist in the current working directory.

    """
    base = str(base)
    ext = str(ext)

    # Add the leading dot to any specified file extension, if necessary
    # (checking type to produce a less obscure error below if not a string):
    if ext and not ext.startswith(os.extsep):
        ext = os.extsep + ext

    # Python doesn't provide a (non-deprecated) way to produce a temporary
    # filename without actually creating and opening the file (to avoid
    # possible race conditions & exploits). One can, however, let Python close
    # the file again and then recycle its name, saving the corresponding
    # DataFile immediately to avoid possible collisions.
    with tempfile.NamedTemporaryFile(
        prefix='{0}_{1}{2}'.format(purpose, base, '_' if base else ''),
        suffix=ext, dir='') as tmpfile:

        tmpname = tmpfile.name

    return tmpname if full_path else os.path.basename(tmpname)


def map_API_enum(name, value, map_dict):
    """
    Convert an enumerable parameter value from the Python API to its equivalent
    IRAF value in the supplied dictionary (where such a mapping exists),
    raising an appropriate exception if it's not recognized.
    """
    try:
        return map_dict[value]
    except KeyError:
        raise ValueError('unrecognized value for \'{0}\''.format(name))

