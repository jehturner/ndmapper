# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os.path
from copy import copy

import numpy as np

from ndmapper import config, ndprocess_defaults
from ndmapper.data import FileName, DataFile, DataFileList, NDLater
from ndmapper.utils import convert_region, to_datafilelist

__all__ = ['init_bpm', 'add_bpm']


@ndprocess_defaults
def init_bpm(reference, regions, convention='numpy', value=None,
             filename=None, reprocess=None):
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

    # Create (or read from disk, if not reprocessing) the output DataFile.
    # This is provisional, until a more definitive method is implemented,
    # possibly as part of the decorator.
    mode = 'new' if reprocess is None else \
           'update' if reprocess is False and os.path.exists(str(filename)) \
           else 'overwrite'
    output = DataFile(filename, mode=mode, labels=labels)

    if mode == 'update':
        return output

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
            slices = convert_region(region, convention)
            ndd.data[slices] = value

    return output


@ndprocess_defaults
def add_bpm(inputs, bpm, out_names=None, reprocess=None):
    """
    Incorporate an existing bad pixel mask into the corresponding flags
    attributes of the input DataFile instances.

    Parameters
    ----------

    inputs : `DataFileList` or `DataFile`
        Input images to which the bad pixel mask is to be added (via logical
        OR with any existing ``flags`` information).

    bpm : `DataFileList` or `DataFile`
        Bad pixel mask, as produced by `init_bpm`, with integer data quality
        flags stored in the main data array of each constituent NDLater
        instance. The length (number of NDLater instances) and shape of each
        constituent array must match the ``inputs``.

    out_names : convertible to `str`, optional
        Output filenames. If None (default), the names of the DataFile
        instances returned will be constructed from those of the input files,
        prefixed with 'b'.


    Returns
    -------

    DataFileList
        Copies of the ``inputs``, with the bad pixel mask values included in
        the flags attribute of each constituent `NDLater` instance.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """

    # Some of this bookeeping should be replaced by a standard wrapper.
    # Also, logging needs to be implemented, as for IRAF.
    inputs = to_datafilelist(inputs)
    if out_names is None:
        out_names = [FileName(df, prefix='b', dirname='') for df in inputs]

    mode = 'new' if reprocess is None else 'overwrite'
    outlist = DataFileList(mode=mode)

    # Here the wrapper would ensure that input/output lengths are the same.

    len_bpm = len(bpm)

    # Loop over the files:
    for in_df, out_name in zip(inputs, out_names):

        if len(in_df) != len_bpm:
            raise ValueError('input {0} and BPM have different lengths'\
                             .format(str(in_df)))

        # Read and use any existing file if not reprocessing:
        if reprocess is False and os.path.exists(str(out_name)):
            out_df = DataFile(out_name, mode='update')

        # Otherwise, OR each NDLater instance of the BPM with the corresponding
        # flags extension of the input file:
        else:
            # Make a copy the inputs, re-using the data by reference since
            # replacing the flags attribute won't affect the inputs.
            out_df = DataFile(out_name, data=in_df, mode=mode)

            # Set flags to the logical OR of the BPM & any existing flags:
            for out_ndd, bpm_ndd in zip(out_df, bpm):

                out_ndd.flags = copy(bpm_ndd.data) if out_ndd.flags is None \
                                else bpm_ndd.data | out_ndd.flags

        outlist.append(out_df)

    return outlist

