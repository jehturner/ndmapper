# Copyright(c) 2015-2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import hashlib

import numpy as np

from astropy.nddata import StdDevUncertainty
from astropy.table import Table
import astropy.io.fits as pyfits

from ._util import get_backend_fn


__all__ = ['NDMapIO', 'TabMapIO']


class NDMapIO(object):
    """
    Propagate additional information needed for `NDLater` instances to support
    lazy loading, allow saving only arrays/header attributes that have changed
    & report which FITS extensions they came from for IRAF etc.

    For lazy-loading or saving operations to succeed, the corresponding file
    must already exist. This class is intended to encapsulate bookkeeping
    within `NDLater` (managed by a `DataFile` instance) with reasonable
    overheads, rather than to provide a robust API: for the user-level
    interface, see `DataFile` instead.

    Attributes
    ----------

    filename : `str`
        The path to the file from which the data are to be mapped.

    ident : `int` or `str` or `None`
        Group identifier appropriate for the file type (int EXTVER for FITS),
        which labels this particular `NDData` instance within a `DataFile`.

    data_idx : `int`
    uncertainty_idx : `int` or `None`
    flags_idx : `int` or `None`
        The original index of each constituent data/uncertainty/flags array
        within the host file (extension number for FITS).

    """

    _data_hash = None
    _uncertainty_hash = None
    _flags_hash = None

    def __init__(self, filename, ident=None, data_idx=None, \
        uncertainty_idx=None, flags_idx=None):

        # This must maintain a separate copy of the host object's filename,
        # otherwise lazy loading of data not yet in memory will fail when
        # changing the filename of a DataFile instance and trying to save it.

        # This should perhaps be changed to cache a reference to its data so
        # that one NDLater instance instantiated from another will share the
        # same data arrays independently of whether lazy loading is triggered
        # before or after instantiation. Once one of them is saved, it will
        # still get re-mapped independently.

        if not isinstance(filename, basestring):
            raise ValueError('filename must be supplied as a string')

        self.filename = filename
        self.ident = ident
        self.data_idx = data_idx
        self.uncertainty_idx = uncertainty_idx
        self.flags_idx = flags_idx

        self._dloader = get_backend_fn('load_array', self.filename)
        self._mloader = get_backend_fn('load_array_meta', self.filename)
        self._saver = get_backend_fn('save_array', self.filename)

        # Consider automatically determining ident from the input here
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


class TabMapIO(object):
    """
    A proxy object for lazily loading/saving AstroPy Table instances. This is
    similar to `NDMapIO`, but instead of being used by an `NDData` sub-class to
    load its own attributes lazily, `TabMapIO` is used to initialize a normal
    `Table` instance on demand, since the latter doesn't have several data
    arrays to load separately and sub-classing `Table` would likely prove more
    complicated with less benefit.

    At the user level, instances are managed by, and the corresponding table
    data accessed via, `DataFile` objects. For lazy-loading or saving
    operations to succeed, the corresponding file must already exist.

    Attributes
    ----------

    filename : `str`
        The path to the file from which the data are to be mapped.

    label : `str`
        Application-specific label/name identifying the type of `Table`
        (EXTNAME for FITS). Multiple tables of the same type can be
        distinguished via the ident parameter.

    ident : `int` or `str` or `None`
        Identifier appropriate for the file type (int EXTVER for FITS), which
        distinguishes this particular instance of a given type of Table within
        the applicable DataFile.

    idx : `int`
        The original array index/number within the host file (extension number
        for FITS).

    """

    _table = None

    def __init__(self, filename, idx, label=None, ident=None):

        # This must maintain a separate copy of the host object's filename,
        # otherwise lazy loading of data not yet in memory will fail when
        # changing the filename of a DataFile instance and trying to save it.

        if not isinstance(filename, basestring):
            raise ValueError('filename must be supplied as a string')

        self.filename = filename
        self.idx = idx
        self.label = label
        self.ident = ident

        self._dloader = get_backend_fn('load_table', self.filename)
        self._mloader = get_backend_fn('load_table_meta', self.filename)
        # self._saver = get_backend_fn('save_table', self.filename)

    def load_data(self):
        data = self._dloader(self.filename, self.idx)
        return data

    def load_meta(self):
        meta = self._mloader(self.filename, self.idx)
        return meta

    def load_table(self):
        meta = self.load_meta()
        data = self.load_data()
        self._table = Table(data=data, meta=meta, copy=False)

    @property
    def table(self):
        if not self._table:
            self.load_table()
        return self._table

    @table.setter
    def table(self, value):

        # Should this preserve the existing label & ident? Should it update
        # them in the new Table's meta (which would mean making a copy)?
        # EXTNAME & EXTVER should probably be removed while in memory instead.

        # Avoid converting existing Table instances to Table because that
        # converts .meta from an io.fits header to an OrderedDict, which it
        # turns out can choke on some odd values such as HISTORY.
        if not isinstance(value, Table):
            try:
                value = Table(value, copy=False)
            except ValueError:
                raise TypeError('value of .table must be convertible to Table')

        self._table = value

    def copy(self):
        """
        Generate a new instance that shares any already-loaded data but can
        be re-mapped independently.
        """
        newinst = TabMapIO(self.filename, self.idx, self.label, self.ident)
        newinst._table = self._table
        return newinst

