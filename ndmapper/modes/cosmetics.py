# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import numpy as np

from ndmapper import config, ndprocess_defaults
from ndmapper.data import DataFile, NDLater
from ndmapper.utils import convert_region as _convert_region

__all__ = ['init_bpm']


@ndprocess_defaults
def init_bpm(reference, regions, convention='numpy', value=None,
             filename=None):
    """
    Initialize a bad pixel mask array from a string listing the corresponding
    regions in NumPy or IRAF/FITS syntax.


    Parameters
    ----------

    reference : DataFile
        DataFile whose NDLater instances the bad pixel mask should match
        (in number and array dimensions).

    regions : list of list of str
        A list containing one sub-list of bad regions per NDLater instance of
        the reference DataFile. Each bad region is specified as a string in
        NumPy/IRAF syntax, eg. '100:110,50:60', '100,*' or ':,99'. The indices
        can be specified in either NumPy or FITS ordering (see below).

    convention : str
        Indexing convention for the region strings: 'NumPy' (default) or
        'FITS' (case insensitive). The former uses 0-based indexing by
        (..., plane, row, column), similar to C or matrix notation, exclusive
        of the ending row/column, while the latter uses 1-based indexing by
        (column, row, plane, ...), similar to ForTran or IRAF and inclusive
        of the ending row/column (with the same origin in either case).

    value : int, optional
        Mask value for bad regions (default 1).

    filename : str or FileName, optional
        Filename for output DataFile.

    """

    # Use default bad pixel value of 1:
    if value is None:
        value = 1

    # Re-map the flags/data-quality array labels to the main science array for
    # the special purpose of keeping a master bad pixel mask (NDData doesn't
    # allow for .flags without associated .data). This method could do with
    # tidying up a bit, as simply giving labels=config['labels']['flags'] to
    # DataFile without also resetting labels['flags'] causes duplication, due
    # to the condition for including an array in the lists in DataFile.save():
    labels = config['labels'].copy()
    labels['data'] = labels['flags']
    labels['flags'] = None

    # New object for the output:
    output = DataFile(filename, mode='new', labels=labels)

    # Create an empty integer array of the same shape as each reference
    # NDData instance:
    for ndd in reference:
        # Do we want meta=ndd.meta.copy()? Probably not.
        output.append(NDLater(np.zeros(ndd.data.shape, dtype=np.uint16)))

    # Convert each region string from NumPy or FITS region string format to
    # a NumPy slice object and use it to set the specified flag value for the
    # corresponding pixels:
    for ndd, reglist in zip(output, regions):
        for region in reglist:
            slices = _convert_region(region, convention)
            ndd.data[slices] = value

    return output

