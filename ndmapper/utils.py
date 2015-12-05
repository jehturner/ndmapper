# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from ndmapper.io import FileName
from ndmapper.data import DataFile


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


def to_filename_strings(objects, strip=True):
    """
    Convert a list of str, FileName or DataFile objects (or a DataFileList)
    to their string representations, by default removing any path and
    processing suffix/prefixes. A single object of these types, or calibration
    dictionary in the format produced by calibrations.init_cal_dict(), may
    also be given. It is the caller's responsibility to ensure that the string
    values are actually valid filenames.

    """
    # Convert any single objects to a list (partly to ensure we don't
    # inadvertently iterate over the NDData instances of a DataFile):
    if isinstance(objects, (DataFile, FileName, basestring)):
        objects = [objects]

    # If the objects argument looks like a calibration dict, extract a list
    # of unique constituent filenames from all the calibrations:
    elif hasattr(objects, 'keys') and 'calibrations' in objects:
        objects = list(set([fn for flist in \
                            objects['calibrations'].itervalues() \
                   for fn in (flist if hasattr(flist, '__iter__') else [])]))

    if not hasattr(objects, '__iter__'):  # must be a list if converted above
        raise ValueError('objects parameter has an unexpected type')

    return [str(FileName(str(obj), strip=strip, \
                         dirname='' if strip else None)) \
            for obj in objects]

