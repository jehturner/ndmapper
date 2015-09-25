# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# The functions defined in this file determine the API & doc-string for their
# format-specific counterparts, to which the work of loading/saving/etc. is
# delegated. Back-end look-up is done automatically by the _get_loader
# decorator, such that these can be minimal definitions.

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

    filename : FileName or str
        Name of an existing file from which to get the meta-data common to all
        groups of pixel data within that file.

    Returns
    -------

    dict-like
        A meta-data dictionary or compatible object.

    """
    return loader(filename)


@_get_loader
def load_array_meta(loader, filename, index):
    """
    Load the meta-data data associated with an array from the specified index
    within a file.
    """
    return loader(filename)


@_get_loader
def load_array(loader, filename, index):
    """
    Load a data array from the specified index within a file.
    """
    return loader(filename)


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
def map_file(loader, filename):
    """
    Open an existing file and return a list of corresponding NDMapIO instances.
    """
    return loader(filename)

