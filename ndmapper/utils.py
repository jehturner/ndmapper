# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

"""
Some high-level utilities for direct use in scripts & processing functions
(or which depend on the high-level data representations used in NDMapper).
"""

from ndmapper.data import FileName, DataFile, DataFileList

from .calibrations import K_CALIBRATIONS


def convert_region(region, convention):
    """
    Convert a NumPy- or FITS-style region string into a tuple of integers
    and Python slice objects that is suitable for subscripting arrays.

    Some corner cases (eg. involving '*', ':' or steps) currently behave
    differently from Python/IRAF subscripting but should be fairly harmless.

    Parameters
    ----------

    region : str
        Region string, eg. '100:110,50:60', '100,*' or ':,99'.

    convention : str
        Indexing convention used: 'NumPy' (default) or 'FITS'; case
        insensitive.

    """

    # Check arguments:
    if not isinstance(region, basestring):
        raise TypeError('region must be a string')
    
    convention = convention.lower()
    if convention not in ['numpy', 'fits']:
        raise ValueError('convention must be NumPy or FITS')

    # Apply appropriate syntax & indexing adjustments for convention used:
    if convention == 'numpy':
        nregion = region
        order = slice(None, None, None)
        orig = 0
    elif convention == 'fits':
        nregion = region.replace('*', ':')
        order = slice(None, None, -1)
        orig = 1
    else:
        raise ValueError('convention must be \'NumPy\' or \'FITS\'')

    # Split region into a range for each axis:
    axes = nregion.split(',')

    # Parse sub-string for each axis into a range & convert to slice object:
    slices = []
    for axis in axes[order]:
        err = False if axis else True    # disallow empty string
        vals = axis.split(':')
        nvals = len(vals)
        if nvals > 3:
            err = True                   # disallow more than start:stop:step
        elif nvals == 1:
            try:
                sliceobj = int(vals[0])-orig  # single row/column number
            except ValueError:
                err = True
        else:
            try:
                # Any adjustment for 1-based indexing is applied only to the
                # start of the range, since 1 has to be added to the ending
                # index to account for FITS/IRAF-style ranges being inclusive.
                sliceobj = slice(*(int(val)-adj if val else None \
                                   for val, adj in zip(vals, (orig, 0, 0))))
            except ValueError:
                err = True               # disallow non-numeric values etc.
        if err:
            raise ValueError('failed to parse region: [%s]' % region)

        slices.append(sliceobj)

    return tuple(slices)


def to_filename_strings(objects, strip_names=True, strip_dirs=True,
                        use_cal_dict=False):
    """
    Extract a list of filename strings from one or more `str`, `FileName` or
    `DataFile` objects (or a `DataFileList`), by default removing any path and
    processing suffix/prefixes. A calibration dictionary may also be given,
    in the format produced by calibrations.init_cal_dict(), if the relevant
    option is enabled. It is the caller's responsibility to ensure that the
    string values are actually valid filenames.

    This is typically used to reproduce base filenames for use in either
    downloading or looking up external information about the files.

    """
    # Convert any recognized single objects to a list (partly to ensure we
    # don't inadvertently iterate over the NDData instances of a DataFile):
    if isinstance(objects, (DataFile, FileName, basestring)):
        objects = [objects]

    # If the objects argument looks like a calibration dict, extract a list
    # of unique constituent filenames from all the calibrations:
    elif use_cal_dict and hasattr(objects, 'keys') and \
         K_CALIBRATIONS in objects:
        objects = list(set([fn for flist in \
                            objects[K_CALIBRATIONS].itervalues() \
                   for fn in (flist if hasattr(flist, '__iter__') else [])]))

    # This must be list-like after the above:
    if not hasattr(objects, '__iter__') or hasattr(objects, 'keys'):
        raise ValueError('objects parameter has an unexpected type')

    return [str(FileName(str(obj), strip=strip_names, \
                         dirname='' if strip_dirs else None)) \
            for obj in objects]


def to_datafilelist(arg, mode=None):
    """
    Derive a `DataFileList` object from one or a sequence of objects than can
    be converted to filename strings (eg. `str`, `FileName` or `DataFile`).
    The ``mode`` defaults to ``'read'`` when given one or more filenames and
    to the existing mode for `DataFile` and `DataFileList`. Where all of the
    inputs are existing `DataFile` instances, those are re-used by reference
    (``mode`` permitting), instead of opening new copies.

    Beware of feeding this inappropriate argument types, as almost anything can
    be converted to `str`...
    """

    # Ensure the input is in a list, for iteration:
    if isinstance(arg, DataFile) or not hasattr(arg, '__iter__'):
        arg = [arg]

    # If passed existing DataFile(s), use them as the data argument for
    # instantiation. Currently, DataFiles can only be re-used by reference if
    # all elements have that type. This case also covers DataFileList, which
    # needs re-instantating in case a new mode is specified.
    if arg and all([isinstance(df, DataFile) for df in arg]):
        outlist = DataFileList(data=arg, mode=mode)

    # Otherwise, convert inputs to filename strings and instantiate with those.
    # Should the acceptable types be enumerated to avoid strange results? That
    # might impact user creativity in both good and a bad ways...
    else:
        mode = 'read' if mode is None else mode
        outlist = DataFileList(filenames=[str(fn) for fn in arg], mode=mode)

    return outlist

