# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# Code & documentation based on NDData: Copyright (c) 2011-2015, Astropy
# Developers. Licensed under a 3-clause BSD style license - see ASTROPY.rst.

# Need to understand the differences between AstroPy and Gemini Python
# docstring conventions better, which are supposedly the same but one has
# extra fields.

import os.path
from copy import deepcopy

import numpy as np

#from astropy.units import Unit, Quantity
from astropy.nddata import NDDataBase, NDData, NDDataArray
from astropy.utils.compat.odict import OrderedDict

from . import config
from . import io as ndmio
from .io import FileName, NDMapIO


class DataFile(object):
    """
    A class representing a data-file-like object, including a filename and/or
    a list of associated NDData instances, along with any ancillary data such
    as a primary header, binary tables that describe the file as a whole and
    aren't part of the NDData structure itself or processing history.
    
    This class can simply store a collection of associated NDData instances in
    memory or it can simply store a filename, eg. for use with IRAF tasks that
    do disk-based I/O, or it can do both, handling any loading & saving. In
    principle this abstraction allows mixing Python code operating on NDData
    fairly seamlessly with steps that do disk I/O (eg. IRAF tasks, as long as
    expectations regarding metadata are compatible), as well as providing a
    convenient file-like way of organizing an NDData collection.

    Parameters
    ----------

    data : NDData or list of NDData or DataFile or None
        NDData instance(s) for the dataset(s) to be represented, or an
        existing DataFile instance (in which case the result will be a new
        DataFile referring to the same NDData instance(s) as the original,
        until such time as it is saved [not implemented]).

        If "filename" refers to an existing file and "data" is None, the data
        will be mapped lazily from that file (data=[] can be used to avoid
        this if the intention is to replace any existing data). Otherwise,
        the user-specified data will be used instead.

    meta : dict-like
        A meta-data dictionary / header describing the file as a whole
        (rather than any individual nddata object), replacing any existing
        meta-data if the file already exists on disk.

    filename : str or FileName, optional
        The filename on disk of the dataset(s) to be represented.

    mode : str
        'read' (default), 'new' or 'overwrite'
        Specifies whether the file is expected to exist on disk already and be
        read into this DataFile (only if a filename is provided). Since
        ndmapper performs reads and writes on request (rather than keeping the
        corresponding files open), DataFile instances can be read or written
        later irrespective of the mode, but "mode" controls whether to create
        this instance from an existing file and whether existence or
        non-existence are allowed (largely to avoid unexpected results, eg.
        when working in the wrong directory). With 'read', the specified file
        must already exist, with 'new', it must not exist and with 'overwrite',
        any existing file is ignored and will be replaced when writing to disk.
        The 'data' and 'filename' parameters always override whatever would
        otherwise be read from disk.

    strip : bool
        Remove any existing prefix and suffixes from the supplied filename
        (prior to adding any specified prefix & suffix)?

    prefix : str, None
        Prefix string to add before the base filename.

    suffix : str, None
        Suffix string to add after the base filename (including any initial
        separator).

    dirname : str, None
        Directory name to add to the filename (replacing any existing dir).

    Attributes
    ----------

    filename : FileName or None
        A filename-parsing object representing the path on disk.

    meta : dict-like
        The header/meta-data associated with the file as a whole (eg. the
        primary FITS header) rather than an individual nddata instance.

    data : list of NDData or None
        NDData instance(s) associated with the DataFile.

    """

    _filename = None
    meta = None
    log = ''

    _dirty = True    # has the file been modified since saved/loaded?


    def __init__(self, data=None, meta=None, filename=None, mode='read',
        strip=False, prefix=None, suffix=None, dirname=None):

        if isinstance(data, DataFile):  # create new copy of existing obj
            self.data = data.data
            self.meta = deepcopy(data.meta)
            self._filename = deepcopy(data.filename)
        elif isinstance(data, NDDataBase):
            self.data = [data]
        elif hasattr(data, '__iter__') and \
            all([isinstance(d, (NDDataBase)) for d in data]):  # True for []
            # Sequence of NDData (or empty to ignore existing file data):
            self.data = data
        elif data is None:
            self.data = None
        else:
            raise TypeError('DataFile: data parameter has an unexpected type')

        # Use any filename & meta-data inherited from input DataFile unless
        # specified:
        if not filename:
            filename = self._filename
        if meta is not None:
            self.meta = meta

        # Parse any filename into a FileName object:
        self._filename = FileName(filename, strip=strip, prefix=prefix,
            suffix=suffix, dirname=dirname)

        # Load the file if requested and it contents haven't been overridden
        # and check that it exists or doesn't, to match expectations:
        exists = os.path.exists(str(self.filename))
        read_file = False
        if str(self.filename):  # Ignore mode=='read' if there's no filename
            if mode == 'read':
                read_file = True
                if not exists:
                    raise IOError('%s not found' % str(self.filename))
            elif mode == 'new':
                if exists:
                    raise IOError('%s already exists' % str(self.filename))
            elif mode != 'overwrite':
                raise ValueError('unrecognized file mode, \'%s\'' % mode)

        if read_file and self.meta is None:
            self._load_meta()
        if read_file and data is None:
            self._load_data()

        # Ensure (meta-)data attributes have the right types & track how many
        # NDData objects this instance contains. An empty DataFile has length
        # zero.
        if self.data is None:
            self.data = []
        if self.meta is None:
            self.meta = {}
        self._len = len(self.data)

    # Make the filename attribute (not its component values) read-only, as the
    # object is shared between DataFile & NDLater to keep things in sync:
    @property
    def filename(self):
        return self._filename

    # Specify that class manages its own iteration with next() method:
    def __iter__(self):
        self._n = 0
        return self

    # Iteration over the DataFile instance returns ~NDData:
    def next(self):
        self._n += 1
        if self._n > self._len:
            raise StopIteration
        else:
            return self.data[self._n-1]

    # Allow subscripting this instance to get ~NDData:
    def __getitem__(self, key):
        return self.data[key]

    # Deletion of an item from the DataFile:
    # This is a simplistic implementation for now; it will probably need to
    # handle some cleaning up etc.
    def __delitem__(self, key):
        del self.data[key]

    # When printing the instance, show the filename:
    def __str__(self):
        return str(self.filename)

    # Need to fix this to do something sensible when there's no filename??
    # This produces an unquoted string in lists!
    def __repr__(self):
        return 'DataFile \'%s\' (len %d)' % (self.__str__(), self._len)

    # Facilitate appending new NDData instances to the DataFile, If the
    # input is another DataFile instance, any existing filename is dropped.
    def append(self, elements):
        if not isinstance(elements, DataFile):
            elements = DataFile(elements)
        if elements.data:
            if self.data is None:
                self.data = []
            self.data += elements.data
        if self.data is not None:
            self._len = len(self.data)

    def __len__(self):
        return self._len

    # Lazily-load the data from file. This is not closely tied to FITS but only
    # flat lists of NDData objects are supported, rather than any arbitrary
    # hierarchy supported by, say, HDF5.
    def _load_data(self):
        self.data = [NDLater(iomap=iomap) for iomap \
                     in ndmio.map_file(self.filename)]

    def _load_meta(self):
        self.meta = ndmio.load_common_meta(self.filename)

    def reload(self):
        """
        Re-load NDData instances & shared meta-data from the associated file
        on disk (to synchronize the DataFile instance with any changes made by
        external programs such as IRAF).

        When instantiating a new DataFile object, it is unnecessary to run
        reload() afterwards if the associated file already exists; it will be
        read automatically.

        Note that data arrays are not actually copied into memory here; they
        are re-mapped and still lazily-loaded once referenced (if applicable).
        """
        if not str(self.filename):
            raise IOError('Attempt to re-load DataFile object with no ' \
                          'associated file')

        # If the file doesn't exist etc., just pass through the IOError from
        # PyFITS, without masking the origin of any more obscure errors:
        self._load_meta()
        self._load_data()
        self._len = len(self.data)

    def save(self):
        """
        Save NDData instances & common meta-data to the associated file,
        creating it if it doesn't already exist. There should be no issue with
        clobbering existing files (unless they happen to be created by another
        program during execution) if the DataFile has been instantiated with
        the right mode, which involves a check for pre-existence.

        """
        # Consider saving the mode in self so we can re-check here that the
        # file (non-)existence is still as expected.

        # This code should be made to append None values for data groups that
        # have not changed since last loaded or saved, which the back-end will
        # then skip re-writing if it can. However, there are various ways for
        # this to backfire if not done very carefully (aka. premature
        # optimization) so this is left as an exercise for later.

        data_name = config['data_name']
        uncertainty_name = config['uncertainty_name']
        flags_name = config['flags_name']

        # Construct flat lists for the array, meta & (name, group_id) tuple,
        # to pass to the save_list function:

        data_list, meta_list, identifiers = [], [], []

        for ndd in self.data:

            group_id = ndd._io.group_id

            arr_group = (ndd.data,
                         ndd.uncertainty.array**2 if ndd.uncertainty else None,
                         ndd.flags)

            meta_group = (ndd.meta, None, None)

            id_group = ((data_name, group_id), (uncertainty_name, group_id),
                        (flags_name, group_id))

            # Include list entries for the main data array and only non-empty
            # uncertainty/flags (passing None values to save_list for those
            # would cause any existing information to be preserved at the
            # applicable location in the file, which isn't what we want).
            for arr, meta, ident in zip(arr_group, meta_group, id_group):
                if arr is not None or ident[0] == data_name:
                    data_list.append(arr)
                    meta_list.append(meta)
                    identifiers.append(ident)

        ndmio.save_list(self.filename, data_list, meta_list, identifiers,
                        self.meta)

        # Here we should re-index NDLater iomaps to reflect the actual order
        # in which the arrays got saved (as determined above) and update the
        # filename to switch lazy-loading to the recently-saved file.


class DataFileList(list):
    """
    A class that holds a list of DataFile objects, tracking filenames
    and/or NDData collections with ancillary information. This implementation 
    is pretty much a normal Python list but provides a more convenient
    interface for instantiating multiple DataFile objects from a list of
    filenames or nddata instances/lists. 

    Parameters
    ----------

    data : NDData or list of NDData or list of lists of NDData or DataFile
           or DataFileList, optional
        NDData/DataFile instance(s) for the dataset(s) to be represented (or
        an existing DataFileList instance). Any member DataFile instances will
        become new copies if any of the filename-modifying parameters are set,
        otherwise the DataFileList will simply hold references to the original
        instances (allowing manipulation of existing DataFiles via new lists).

    meta : dict-like or list of dict-like or None, optional
        The header/meta-data associated with each file as a whole (eg. the
        primary FITS header) rather than with individual nddata instances.
        There should be one instance per file or None (which preserves any
        information from an existing file).

    filenames : str or list of str, optional
        The filename(s) on disk of the dataset(s) to be represented. There
        should either be as many filenames as data or None (but lists of
        NDData can be nested to associate subsets with fewer filenames).

    mode : str
        'read', 'new' or 'overwrite'
        Specifies whether the files are expected to exist on disk already and
        be read into the corresponding DataFile instances (only if filenames
        are provided; also see DataFile). With 'read', the specified files must
        already exist, with 'new', they must not exist and with 'overwrite',
        any existing files are ignored and will be replaced when writing to
        disk. The 'data' and 'filename' parameters always override whatever
        would otherwise be read from disk.

    strip : bool
        Remove any existing prefix and suffixes from the supplied filename
        (prior to adding any specified prefix & suffix)?

    prefix : str, None
        Prefix string to add before the base filename.

    suffix : str, None
        Suffix string to add after the base filename (including any initial
        separator).

    dirname : str, None
        Directory name to add to the filename (replacing any existing dir).

    """

    filenames = []

    def __init__(self, data=None, meta=None, filenames=[], mode='read',
        strip=False, prefix=None, suffix=None, dirname=None):

        # This is a little fiddly (not too much considering) but handling
        # these cases here should make things conceptually simpler and/or
        # more readable from a user perspective...

        # First convert arguments to lists if they aren't already. There isn't
        # really a simple way of doing this with a function because for some
        # arguments we want to distinguish specific list-like objects from a
        # list of those objects, for others, we want to distinguish duck-typed
        # sequence objects from a list etc. -- ie. the criteria vary. Also, it
        # seems safer to require specific input types in some cases and warn
        # the user if something unexpected is received.

        # Deal with special cases of single objects that we accept for
        # convenience, although we are instantiating a list of things. These
        # cases are enumerated to distinguish them from container lists.
        if _compatible_data_obj(data):
            data = [data]
        elif data is not None and not hasattr(data, '__iter__'):
            # If the constituent elements of the list don't have the right
            # type, let DataFile complain rather than checking here.
            raise TypeError('DataFileList: data parameter has an unexpected ' \
                            'type')
        len_data = seqlen(data)

        if isinstance(filenames, basestring):
            filenames = [filenames]
        elif filenames is not None and not hasattr(filenames, '__iter__'):
            raise TypeError('DataFileList: filenames parameter has an ' \
                            ' unexpected type')
        len_fn = seqlen(filenames)

        if hasattr(meta, 'keys'):  # dict-like (incl. PyFITS headers).
            meta = [meta]
        len_meta = seqlen(meta)

        # Determine whether the filename is being modified, to help decide
        # below whether a new copy of the input is needed:
        fn_modified = filenames or strip or prefix or suffix \
            or dirname is not None

        # If instantiating with a list or DataFileList of unmodified
        # DataFiles (no filename prefix changes etc.), just reference
        # the existing instances by default, instead of creating copies,
        # to allow making new lists of existing DataFile instances:
        if isinstance(data, list) \
           and all([isinstance(item, DataFile) for item in data]) \
           and meta is None and not fn_modified:
            initlist = data

        # Otherwise create a new copy of each DataFile to hold mods. Here
        # we should have an existing DataFileList instance, a list of DataFile
        # objects, list of NDData or list of lists of NDData.
        else:
            # Expand out both data & filenames to lists of the same length,
            # to avoid repetition below.
            if len_data < 2:
                # This can produce either [None] or [[]], the latter of which
                # overrides any existing data when instantiating DataFile.
                if not data:
                    data = [data]
                # Only if both data & filenames are 0/None, expand them to
                # match meta:
                if meta and not filenames:
                    filenames = [None for mdict in meta]
                data = [data[0] for fn in filenames]
            # Filenames need to be unique so there must be None or as many
            # as there are items in data (but data can include singly-
            # nested lists of NDData, with each sub-list associated with
            # one filename and one resulting DataFile).
            elif not filenames:
                filenames = [None for obj in data]
            # Should we also be expanding out a single data object to the
            # number of filenames?
            elif len_fn != len_data:
                raise ValueError('DataFileList: data & filenames differ ' \
                    'in length')
            # If meta is None or a single dict, expand it to match the
            # length of data & filenames, which are by now the same. Meta
            # only determines the number of data files produced (above)
            # if both filenames and data are None.
            if len_meta < 2:
                if not meta:
                    meta = [None]
                meta = [meta[0] for fn in filenames]
            if len(meta) != len(filenames):
                raise ValueError('DataFileList: meta does not match ' \
                    'data/filenames in length')
            initlist = [DataFile(data=obj, meta=mdict, filename=fn,
                mode=mode, strip=strip, prefix=prefix, suffix=suffix,
                dirname=dirname) for obj, mdict, fn in \
                zip(data, meta, filenames)]

        # Do whatever initialization a list object normally does:
        list.__init__(self, initlist)

    # Wrap the normal list append in the same way as __init__:
    def append(self, data=None, meta=None, filename=None, strip=False, \
        prefix=None, suffix=None, dirname=None):

        if isinstance(data, DataFile) and not (meta or filename or strip or \
            prefix or suffix or dirname is not None):
            list.append(self, data)
        else:
            list.append(self, DataFile(data=data, meta=meta, filename=filename,
                mode=mode, strip=strip, prefix=prefix, suffix=suffix,
                dirname=dirname))

    # Wrap the normal list extend in the same way as __init__:
    def extend(self, data=None, meta=None, filenames=[], strip=False, \
        prefix=None, suffix=None, dirname=None):

        list.extend(self, DataFileList(data=data, meta=meta, \
            filenames=filenames, mode=mode, strip=strip, prefix=prefix,
            suffix=suffix, dirname=dirname))


def seqlen(arg, convert_empty=False):
    """
    Return the length of the argument if a sequence, otherwise 1 or None.
    """
    if hasattr(arg, '__iter__'):   # no strings etc.
        try:
            slen = len(arg)
        except TypeError:
            raise TypeError('seqlen: unexpectedly got an iterable object ' \
                            'with no length!')

        if convert_empty is True and slen == 0:
            slen = None

    elif arg is None:
        slen = None

    else:
        slen = 1

    return slen
    

def _compatible_data_obj(arg):
    # For now DataFile doesn't support ndarray but that's OK, as it reports
    # the same error DataFileList would for unsupported types and this will
    # work if and when it's added there:
    if isinstance(arg, DataFile) or isinstance(arg, NDDataBase) or \
       isinstance(arg, np.ndarray):
        return True
    else:
        return False


class NDLater(NDDataArray):
    """
    A compatible variant of NDDataArray that facilitates lazy loading of
    pixel data, allowing code to work freely with NDData-like objects
    (including just doing bookeeping with the headers) without using more
    memory than necessary.

    The main API difference from NDDataArray is that NDLater is initialized
    with a file name & indices (eg. FITS extension numbers) for getting/
    saving the data on demand. This very simple interface does not include
    any safety checks and it's the caller's responsibility to take care of
    managing the file structure, existence etc., which can be taken care of
    by the higher-level DataFile class. A subset of NDData.__init__()
    parameters are also accepted, allowing creation from an existing NDData
    instance or ndarrays (eg. for writing to a new file); this will override
    any existing data in the file at the specified indices.

    Any mask, wcs & unit values will be derived (once implemented) directly
    from the input meta-data, rather than specified here, and can then be
    overridden after instantiation if necessary.

    Parameters
    ----------

    data : `~numpy.ndarray` or `NDData`, optional.
        The main data array contained in this object, overriding any existing
        data in the mapped file, if applicable. If the intention is to use a
        new copy of the input object rather than a reference to it, the user
        should make that copy beforehand.

    uncertainty : `~astropy.nddata.NDUncertainty`, optional
        Uncertainties on the data.

    flags : `~numpy.ndarray`-like or `~astropy.nddata.FlagCollection`, optional
        Flags giving information about each pixel. These can be specified
        either as a Numpy array of any type (or an object which can be converted
        to a Numpy array) with a shape matching that of the data, or as a
        `~astropy.nddata.FlagCollection` instance which has a shape matching
        that of the data.

    meta : `dict`-like object, optional
        Metadata for this object.  "Metadata" here means all information that
        is included with this object but not part of any other attribute
        of this particular object. e.g., creation date, unique identifier,
        simulation parameters, exposure time, telescope name, etc.

    iomap : NDMapIO, optional
        An object that maps the data, uncertainty & flags attributes to a file
        name and indices within that file, enabling those attribute values to
        be lazily loaded or saved.

    (See NDDataArray doc string for methods & attributes.
     This is a Work in progress, to support DataFile.)

    """

    # After implementing lazy loading with NDLater instead of loading
    # everything in _load_nddata_from_FITS, 11 tests now take ~2.3s total to
    # run instead of ~1.2s (I think with 2 tests that read 2 SCI exts each).

    # This is based on the NDData & NDDataArray __init__ but avoids referencing
    # array attributes here, instead storing an obj that knows how to get them.
    def __init__(self, data=None, uncertainty=None, flags=None, meta=None,
                 iomap=None):

        if iomap and not isinstance(iomap, NDMapIO):
            raise TypeError('iomap must be an NDMapIO instance')

        if data is None and iomap is None:
            raise ValueError('must provide either data or iomap')

        # Remember our "parent class", for later use in getters/setters, where
        # to be on the safe side, we invoke the NDDataArray getter/setter logic
        # after actually loading the data array(s).
        self._parent = super(NDLater, self)

        # Add support for deriving these properly later, when needed. The
        # _unit must be set before the parent __init__ is called below (which
        # looks like a bug in NDDataArray), as the latter calls the uncertainty
        # setter, which uses the unit getter.
        self._mask = None
        self._wcs = None
        self._unit = None

        # Initializing the data (& uncertainty/flags) to None indicates that
        # the data haven't been loaded yet:
        self._data = None

        # Attach the object to which lazy loading is delegated (and which
        # tracks the mapping of attributes to extensions). This must be done
        # before setting some of the other attributes, whose setters need to
        # access the data (which may or may not be what we want for lazy
        # loading but is what currently happens in nddata).
        self._io = iomap

        if data is None:
            # When starting from scratch, the only initialization we can
            # inherit from our ancestors is the most basic stuff that happens
            # in the NDDataBase class (which doesn't do much at present but
            # just in case it does later...):
            super(NDData, self).__init__()

            # Normally self._data, self._uncertainty are set by NDData and
            # self._flags by NDDataArray. They are only used by the relevant
            # attribute getters & setters upstream.
            self._uncertainty = uncertainty
            self._flags = flags

            # If this remains undefined, we'll load it from file below:
            self._meta = meta

            # Initialize attributes via the upstream setter logic of the public
            # API where possible (wcs doesn't have one yet). Some of these
            # setters expect the private attribute version to be defined first:
            self.mask = self._mask
            self.unit = self._unit

        else:
            # When passed data as well as a filename, let our parent class
            # populate the class attributes and then we'll lazily load anything
            # that wasn't provided. It's the caller's responsibility to avoid
            # overriding inconsistent subsets of what's already in the file
            # (eg. data without the corresponding uncertainty) but the DataFile
            # class can help take care of that.
            self._parent.__init__(data, uncertainty=uncertainty, mask=None,
                flags=flags, wcs=None, meta=meta, unit=None)

            # NDDataArray doesn't copy flags from data when it should. We can't
            # duck type this because numpy also has a different "flags".
            if flags is None and isinstance(data, NDDataArray):
                self._flags = data.flags

            # Don't copy any existing NDLater _io attribute, require it to be
            # specified explicitly. If we're making a copy we probably want to
            # write somewhere different from the original and DataFile can
            # take care of this.

        # Don't bother loading the header lazily, but still get it via NDMapIO,
        # to avoid adding I/O logic in more places than necessary. We may also
        # need the meta-data here, eg. to determine things like units & WCS.
        if self._meta is None:
            if self._io:
                self._meta = self._io.load_meta()
            else:
                self._meta = OrderedDict()

        # Because NDDataArray's mask & flags setters refer to data.shape,
        # they will unintentionally trigger lazy loading if we're provided
        # explicitly with a flags argument. Consider overriding the shape
        # getter (etc.) to avoid this, using the metadata directly (with the
        # help of the io library?) or checking for _data==None and doing
        # del _data afterwards if necessary as a temporary workaround.

        # These setters only work after creating the _io attribute above.
        self.uncertainty = self._uncertainty
        self.flags = self._flags

    @property
    def data(self):
        if self._data is None and self._io:
            self._data = self._io.load_data()
        return self._data

    @data.setter
    def data(self, value):
        self._data = value

    @data.deleter
    def data(self):
        # This doesn't always free memory in practice, probably due to PyFITS's
        # default memory mapping, but either way it provides the indended means
        # of dropping our reference to the data.
        self._data = None

    @property
    def uncertainty(self):
        # This prevents resetting the value to None; use del instead
        if self._parent.uncertainty is None and self._io:
            self.uncertainty = self._io.load_uncertainty()
        return self._parent.uncertainty

    # Parent class setters seem not to get called automatically once a getter
    # is defined but unfortunately super() simply doesn't work as a proxy for
    # setters (Python issue 14965) so we have to do it the following way:

    @uncertainty.setter
    def uncertainty(self, value):
        NDDataArray.uncertainty.fset(self, value)

    @uncertainty.deleter
    def uncertainty(self):
        self.uncertainty = None

    @property
    def flags(self):
        # This prevents resetting the value to None; use del instead
        if self._parent.flags is None and self._io:
            self.flags = self._io.load_flags()
        return self._parent.flags

    @flags.setter
    def flags(self, value):
        NDDataArray.flags.fset(self, value)

    @flags.deleter
    def flags(self):
        self.flags = None


def load_datafilelist(filename, dirname=None, mode='read'):
    """
    Open a text file containing a list of filenames as a DataFileList object.

    Parameters
    ----------

    filename : str
        Name of a plain-text file, containing one filename per line.

    dirname : str, optional
        Directory where the files to be read/written are located, to be
        prefixed to each listed filename. This option cannot be used if any
        filenames already include a (full or relative) path.

    mode : str, optional
        Access mode parameter to pass to each instantiated DataFile object.


    Returns
    -------

    DataFileList
        A DataFileList object, mapped to the listed files.

    """

    f = open(filename, 'r')

    dfl = DataFileList()

    # Parse each filename, one per line, check it doesn't have a duplicate
    # path and append it to the DataFileList:
    for line in f:
        line = line.strip()         # remove new lines & trailing/leading space
        if line and line[0] != '#': # ignore empty lines & comments
            if dirname and os.path.dirname(line):
                raise IOError('Specified \'dirname\' when file already has '
                    'one:\n  %s' % line)
            dfl.append(DataFile(filename=line, mode=mode, dirname=dirname))

    f.close()

    return dfl

# To do:
# - Is the DataFile init logic needlessly re-reading any DataFile passed
#   as an argument? More specifically, I think this will be triggered when
#   adding a DataFile to a DataFileList with mode='read'.
# - Implement deepcopy methods?

