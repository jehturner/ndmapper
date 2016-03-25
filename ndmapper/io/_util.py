# Copyright(c) 2015-2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# These are helper functions for the I/O routines in ndmapper.io.

from functools import wraps

from ..libutils import splitext
from .formats import formats


__all__ = ['get_backend_fn']


def get_backend_fn(funcname, filename):
    """
    Given a filename string and the name of a loader function defined in
    ndmapper.io, return the implementation of the latter function from the
    sub-module appropriate for the file format.

    Currently we use only the file extension names, rather than more foolproof
    magic or try-loader-except to determine file types, avoiding unnecessary
    I/O overheads & complication.

    This function may be used either directly by applications wanting to cache
    look-ups when doing repeated I/O operations or, internally (when defining
    new generic functions in ndmapper.io), via the _get_loader decorator.

    """
    fext = (splitext(filename)[1] or '').lower()
    backend_fn = None
    for fmt, vals in formats.iteritems():
        if fext in vals:
            # Import back-end module if not done already; just assume it
            # exists if defined in formats dict, otherwise we have a bug.
            exec('from . import _{0} as {0}'.format(fmt))
            # Try to get the back-end function from the module:
            try:
                backend_fn = eval(fmt+'.'+funcname)
            except AttributeError:
                raise IOError('back end \'%s\' has no function \'%s\'' \
                              % (fmt, funcname))
            break
    if not backend_fn:  # no back-end for file extension
        raise IOError('unsupported file format \'%s\'' % fext)
    return backend_fn


def _get_loader(fn):
    """
    A decorator that calls get_backend_fn() to determine automatically the
    appropriate back-end function corresponding to the generic one called
    directly and provide it to the latter as an additional argument (similar
    to 'self' in classes). Intended for internal use within ndmapper.io.
    """
    @wraps(fn)  # use func's own name & docstring instead of the wrapper's
    def loader_wrapper(*args, **kwargs):
        loader = get_backend_fn(fn.__name__, args[0])
        return fn(loader, *args, **kwargs)
    return loader_wrapper

