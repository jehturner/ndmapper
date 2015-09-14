# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# The functions defined in this file determine the API & doc-string for their
# format-specific counterparts, to which the work of loading/saving/etc. is
# delegated. Back-end look-up is done automatically by the _get_loader
# decorator, such that these can be minimal definitions.

from .formats import formats
from ._util import _get_loader  # since private attributes get excluded from *
from ._util import *            # these go into the public io namespace


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

