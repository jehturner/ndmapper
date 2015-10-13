# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os.path
import string
import re
import hashlib

import numpy as np

from astropy.nddata import StdDevUncertainty
import astropy.io.fits as pyfits

from .. import config
from ._util import get_backend_fn


__all__ = ['FileName', 'NDMapIO']


class FileName(object):
    """
    A class for parsing filenames into components, reconstructing them &
    storing/sharing that information amongst interested parties.

    Parameters
    ----------

    path : str, FileName
        Single filename to parse into a FileName object representation.

    sep : str, None
        Separator for suffix components (defaults to "_").

    strip : bool
        Remove any existing prefix and suffixes from the supplied path
        (prior to adding any specified prefix & suffix)?

    prefix : str, None
        Prefix string to add before the base filename.

    suffix : str, None
        Suffix string to add after the base filename (including any initial
        separator).

    dirname : str, None
        Directory name to add to the path (replacing any existing directory).

    regex : str, re, None
        Regular expression matching root filename (without a file extension).
        By default this is None, causing the value of the package configuration
        variable "ndmapper.config['filename_regex']" to be used, which
        defaults to Gemini's "S20150101S0001"-style convention (thus allowing
        use of other conventions without having to override the regex every
        time a DataFile is instantiated, as well as optionally allowing the
        default to be pre-compiled).


    Attributes
    ----------

    dir : str
        Directory that the file resides in.

    prefix : str
        Sequence of characters preceding the base name.

    base : str
        Base filename in a standard format that can be recognized via the
        regex parameter, eg. S20150307S0001 for Gemini data.

    suffix : list
        List of one or more suffixes following the base name, including any
        separator character (eg. _forStack). 

    ext : str
        File extension(s), eg. "fits" or "fits.gz".

    sep : str, None
        One or more characters specified as a suffix separator.

    """

    def __init__(self, path=None, sep='_', strip=False, prefix=None, \
        suffix=None, dirname=None, regex=None):

        # Get default regex from package configuration so non-Gemini users
        # don't have to specify an alternative convention every time this
        # class is instantiated:
        if regex is None:
            regex = config['filename_regex']

        # Compile regular expression if supplied as a string:
        if isinstance(regex, basestring):
            self._re = re.compile(regex)
        else:
            self._re = regex

        # Record what separator we're using:
        self.sep = sep

        # If passed an existing instance, reconstruct and re-parse it, since
        # the regex or separator can differ (and it's simpler to do).
        if isinstance(path, FileName):
            path = str(path)
        elif path is not None and not isinstance(path, str):
            raise ValueError('path must be a str or %s instance' % \
                             str(self.__class__.__name__))

        # Actually parse the path or use placeholder attributes if it's None:
        if path is None:

            # In this case we want a completely empty string to represent the
            # file, so that (eg.) DataFile objects instantiated with None
            # won't get some anomalous default filename.
            self.dir = ''
            self.prefix = ''
            self.base = ''
            self.suffix = []
            self.ext = None

        else:
            # Separate directory, filename root & file extension:
            self.dir = os.path.dirname(path)
            # This splits at the first dot, unlike os.path.splitext:
            fncmp = os.path.basename(path).split(os.extsep, 1)
            froot = fncmp[0]  # split always produces a root/base name
            if len(fncmp) == 1:
                self.ext = None   # distinguish no ext at all from a dot
            else:
                self.ext = fncmp[1]

            # Separate any prefix and/or suffixes from the base name:
            match = self._re.search(froot)
            if match:
                self.standard = True
                self.base = match.group()
                if strip:
                    self.prefix = ''
                    self.suffix = []
                else:
                    self.prefix = froot[:match.start()]
                    self.suffix = self._split(froot[match.end():])
            else:
                self.standard = False
                self.base = froot
                self.prefix = ''
                self.suffix = []
            
        # Add on any specified prefix or suffix:
        if prefix:
            self.prefix = prefix + self.prefix
        if suffix:
            self.suffix.extend(self._split(suffix))

        # Add or replace any initial directory name if specified:
        if dirname is not None:
            self.dir = dirname

    # Split a string (suffix) into a list that includes the separator
    # character at the start of each element that originally had one (which
    # the first element is not bound to):
    def _split(self, suff):
        ls = suff.split(self.sep)  # always produces at least one element
        # If string starts with a separator, omit empty initial string
        # produced by the split:
        if not ls[0]:
            result = []
        # Otherwise, include the beginning of the string before the first
        # separator as-is:
        else:
            result = [ls[0]]
        # Restore the separator char before any subsequent elements:
        result += [self.sep + s for s in ls[1:]]
        return result

    @property
    def re(self):
        return self._re

    # Reconstruct the filename when printing the instance value (after any
    # user modifications to the individual components):
    def __str__(self):
        if self.ext is None:
            dotext = ''
        else:
            dotext = os.extsep+self.ext
        return (os.path.join(self.dir, (self.prefix+self.base+ \
                string.join(self.suffix, '')+dotext)))

    def __deepcopy__(self, memo):
        # One can't copy a regex, except by re.compile(_re.pattern), but they
        # look to be immutable anyway so this shouldn't be a problem; likewise
        # for strings:
        return FileName(path=str(self), sep=self.sep, regex=self._re)

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self == other


class NDMapIO(object):
    """
    Propagate additional information needed for NDData instances (or user code)
    to support lazy loading, allow saving only arrays/header attributes that
    have changed & report which FITS extensions they came from for IRAF etc.

    For lazy-loading or saving operations to succeed, the corresponding file
    must already exist. This class is intended to encapsulate bookkeeping
    within NDLater (managed by a DataFile instance) with reasonable overheads,
    rather than to provide a robust API: for the user-level interface, see
    DataFile instead.

    Attributes
    ----------

    filename : FileName
        The path to the file from which the data are to be mapped.

    group_id : int or str or None
        Group identifier appropriate for the file type (int EXTVER for FITS),
        which labels this particular NDData instance within a DataFile.

    data_idx : int
    uncertainty_idx : int or None
    flags_idx : int or None
        The original index of each constituent data/uncertainty/flags array
        within the host file (extension number for FITS).

    """

    _data_hash = None
    _uncertainty_hash = None
    _flags_hash = None

    def __init__(self, filename, group_id=None, data_idx=None, \
        uncertainty_idx=None, flags_idx=None):

        # Cast the filename to a new FileName instance whether or not it is
        # one already; we need an independent copy here, otherwise lazy
        # loading of data not yet in memory will fail when changing the
        # filename of a DataFile instance and trying to save the data.
        try:
            filename = FileName(filename)
        except ValueError:
            raise ValueError('must define filename as a str or FileName '
                             'object')

        self.filename = filename
        self.group_id = group_id
        self.data_idx = data_idx
        self.uncertainty_idx = uncertainty_idx
        self.flags_idx = flags_idx

        self._dloader = get_backend_fn('load_array', self.filename)
        self._mloader = get_backend_fn('load_array_meta', self.filename)
        self._saver = get_backend_fn('save_array', self.filename)

        # Consider automatically determining group_id from the input here
        # (ie. setting it to hdu.ver == EXTVER) if None.

    def load_data(self):
        data = self._dloader(self.filename, self.data_idx)
        # A NumPy array is directly hashable -- but doing so pulls memory-
        # mapped data entirely into memory, where they stay until unloaded
        # with "del ndd.data". A workaround of reading the file twice would
        # negate the intended benefit of being able to save it intelligently,
        # so just disable hashing in the first instance and ask Erik B. about
        # it later. It might be better to determine whether the copy is dirty
        # using object ids (weakref) and memory mapping instead, like PyFITS,
        # but that might mean re-reading the file after saving, to establish
        # memory mapping before we can mark the buffer clean.
        # self._data_hash = hashlib.sha1(data).hexdigest()
        return data

    # def save_data(self, data, header, force=False):
    #     # Should hash meta-data as well here, or else we'll lose changes that
    #     # aren't associated with changes to data.
    #     newhash = hashlib.sha1(data).hexdigest()
    #     if force or newhash != self._data_hash:
    #         self._data_hash = newhash
    #         self._saver(self.filename, self.data_idx, data, header)

    def load_uncertainty(self):
        if self.uncertainty_idx:
            uncert = self._dloader(self.filename, self.uncertainty_idx)
            # Presumably this kills any memory mapping? Worry about it later.
            # The sqrt is just a temporary hack until I write a Var subclass.
            # StdDevUncertainty isn't directly hashable so cast to str first
            # (also see load_data above for another reason).
            uncert = StdDevUncertainty(np.sqrt(uncert))
            # self._uncert_hash = hashlib.sha1(uncert).hexdigest()
            return uncert

    def load_flags(self):
        if self.flags_idx:
            flags = self._dloader(self.filename, self.flags_idx)
            # self._flags_hash = hashlib.sha1(flags).hexdigest()
            return flags

    def load_meta(self):
        meta = self._mloader(self.filename, self.data_idx)
        # This cast to str is a little bit slow, so let's see whether the hash
        # here turns out to be premature optimization before reinstating it:
        # self._meta_hash = hashlib.sha1(str(meta)).hexdigest()
        return meta

