# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

"""
Public I/O routines used by NDMapper. These are wrappers for their
file-format-specific back ends (which call astropy.io.fits etc.). Users should
normally work with DataFile objects rather than calling these directly.
"""

# The functions defined in this file determine the API & doc-string for their
# format-specific counterparts, to which the work of loading/saving/etc. is
# delegated. Back-end look-up is done automatically by the _get_loader
# decorator, such that these can be minimal definitions, or can be cached by
# calling get_backend_fn() directly. Note that no common processing can be done
# here if the same behaviour is expected via either look-up method.

from .mapio import *            # want these in the public namespace
from .formats import formats
from ._util import _get_loader  # since private attributes get excluded from *
from ._util import *


@_get_loader
def load_common_meta(loader, filename):
    """
    Open an existing file and return any meta-data common to all data groups
    (eg. a FITS primary header) as a dict-like (eg. PyFITS header) object.

    Parameters
    ----------

    filename : `str`
        Name of an existing file from which to get the meta-data common to all
        groups of pixel data within that file.

    Returns
    -------

    `dict`-like
        A meta-data dictionary or compatible object.

    """
    return loader(filename)


@_get_loader
def load_array_meta(loader, filename, index):
    """
    Load the meta-data data associated with an array from the specified index
    within a file.
    """
    return loader(filename, index)


@_get_loader
def load_array(loader, filename, index):
    """
    Load a data array from the specified index within a file.
    """
    return loader(filename, index)


@_get_loader
def load_table_meta(loader, filename, index):
    """
    Load the meta-data data associated with a table from the specified index
    within a file.
    """
    return loader(filename, index)


@_get_loader
def load_table(loader, filename, index):
    """
    Load a table from the specified index within a file as an array-like
    object.
    """
    return loader(filename, index)


# This exists for completeness & experimentation but in practice we want to
# write the whole file at once from DataFile.
@_get_loader
def save_array(loader, filename, index, data, meta=None):
    """
    Save a data array to the specified index (array number) within an existing
    file. Associated meta-data may be saved at the same time, being needed to
    describe the array fully and to save writing to the file twice (unlike
    load_array_meta & load_array, which are used separately anyway when lazy
    loading and when using a single set of meta-data for data/uncertainty/flags
    arrays). If meta is not supplied, the minimal metadata needed to specify
    the array size, dimensionality & numeric type will be written. The index
    within the file must either exist already (for overwriting) or be the next
    one available (for appending).
    """
    return loader(filename, index, data, meta)


@_get_loader
def save_list(loader, filename, data, array_meta=None, identifiers=None,
              types=None, common_meta=None):
    """
    Save a list of data arrays and associated meta-data to a file (overwriting
    any existing data).

    Where data and array_meta are both None at a given list index, any data
    already present at the corresponding location within the file will be
    preserved, avoiding redundant writes where possible (it is the caller's
    responsibility to determine what needs re-writing and what can be kept).

    Parameters
    ----------

    filename : str
        Name of the (new or existing) file to which the data should be saved.

    data : list of ndarray or None
        Ndarray instances to be saved to each location within the file.

    array_meta : list of dict-like or None, optional
        Header/meta-data dictionaries associated with each data array.

    identifiers : list of (str or None, int or str or None), optional
        Name and group identifier for each array within the file (eg. FITS
        EXTNAME & EXTVER), where applicable overriding any such information
        in array_meta.

    types : list of str, optional
        Type of each data item to be saved: 'image' or 'table' (where the
        applicable file format makes such a distinction). Defaults to 'image'
        for each item.

    common_meta : dict-like, optional
        Header/meta-data dictionary common to all the arrays within the file.

    """
    return loader(filename, data, array_meta, identifiers, types, common_meta)


@_get_loader
def map_file(loader, filename, labels=None):
    """
    Open an existing file and return a list of NDMapIO instances corresponding
    to data/uncertainty/flags image array groups and a second list of TabMapIO
    instances corresponding to any additional data tables.

    Returns
    -------

    (list of NDMapIO, list of TabMapIO)

    """
    return loader(filename, labels)

