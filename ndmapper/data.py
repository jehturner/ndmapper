# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# Code & documentation based on NDData: Copyright (c) 2011-2015, Astropy
# Developers. Licensed under a 3-clause BSD style license - see ASTROPY.rst.

# Need to understand the differences between AstroPy and Gemini Python
# docstring conventions better, which are supposedly the same but one has
# extra fields.

import os.path
from copy import copy, deepcopy
import tempfile

import numpy as np

#from astropy.units import Unit, Quantity
from astropy.nddata import NDDataBase, NDData, NDDataArray
from astropy.table import Table
from astropy.utils.compat.odict import OrderedDict

from . import config
from . import io as ndmio
from .io import FileName, NDMapIO, TabMapIO


__all__ = ['DataFile', 'DataFileList', 'NDLater', 'load_file_list',
           'temp_saved_datafile']


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

    filename : str or FileName, optional
        The filename on disk of the dataset(s) to be represented.

    data : NDData or list of NDData or DataFile or None, optional
        NDData instance(s) for the dataset(s) to be represented, or an
        existing DataFile instance (in which case the result will be a new
        DataFile referring to the same NDData instance(s) as the original,
        until such time as it is saved [not implemented]).

        If "filename" refers to an existing file and "data" is None, the data
        will be mapped lazily from that file (data=[] can be used to avoid
        this if the intention is to replace any existing data). Otherwise,
        the user-specified data will be used instead.

    meta : dict-like, optional
        A meta-data dictionary / header describing the file as a whole
        (rather than any individual nddata object), replacing any existing
        meta-data if the file already exists on disk.

    mode : str, optional
        'read' (default), 'new', 'update' or 'overwrite'
        Specifies whether the file should exist on disk already and be used to
        initialize this DataFile (if a filename is provided) and whether it can
        later be written to disk. Although ndmapper does not hold file handles
        open with a fixed access mode, these options enforce the user's
        declared intention, to avoid mishaps (such as inadvertently overwriting
        input files or failures when working in the wrong directory). With
        'read' and 'update', the specified file must already exist, with 'new',
        it must not exist and with 'overwrite', any existing file is ignored
        and will be replaced when writing to disk. The 'data' and 'filename'
        parameters always override whatever would otherwise be read from disk.
        When no `filename` is given, the default is 'new' instead of 'read'.

    strip : bool, optional
        Remove any existing prefix and suffixes from the supplied filename
        (prior to adding any specified prefix & suffix)?

    prefix : str or None, optional
        Prefix string to add before the base filename.

    suffix : str or None, optional
        Suffix string to add after the base filename (including any initial
        separator).

    dirname : str or None, optional
        Directory name to add to the filename (replacing any existing dir).

    labels : str or dict of str : str, optional
        Naming of each NDData component array in the corresponding file on
        disk, overriding the package default values in config['labels']. Where
        a string is given, it overrides only the label of the main data array
        (ie. config['labels']['data']).

    Attributes
    ----------

    filename : FileName or None
        A filename-parsing object representing the path on disk.

    mode : str
        File access mode (see parameters). This initially reflects the
        corresponding parameter value / default and can change implicitly in
        certain circumstances when loading, saving or changing the filename.

    meta : dict-like
        The header/meta-data associated with the file as a whole (eg. the
        primary FITS header) rather than an individual nddata instance.

    cals : dict of str : DataFile
        A dictionary of associated processed calibrations, mapping calibration
        types (eg. 'bias') to DataFile instances that can be (or have been)
        used by the corresponding calibration steps. This attribute is
        currently not persistent when saving to disk, so any calibrations must
        be associated explicitly on reloading, if needed.

    The NDData instance(s) associated with the DataFile are accessed by
    iterating over or subscripting it like a list.

    """

    _filename = None
    _meta = None
    _tables = None    # change later to _extras?
    _labels = None
    _cals = None

    log = ''


    def __init__(self, filename=None, data=None, meta=None, mode=None,
        strip=False, prefix=None, suffix=None, dirname=None, labels=None):

        if isinstance(data, DataFile):  # create new copy of existing obj
            self._data = data._data
            self._tables = data._tables
            self._meta = deepcopy(data.meta)
            self._filename = deepcopy(data.filename)
            self._labels = copy(data._labels)
            self._cals = data._cals
        elif isinstance(data, NDDataBase):
            self._data = [NDLater(data=data)]
        elif hasattr(data, '__iter__') and \
            all([isinstance(d, (NDDataBase)) for d in data]):  # True for []
            # Sequence of NDData (or empty to ignore existing file data):
            self._data = [NDLater(data=d) for d in data]
        elif data is None:
            self._data = None
        else:
            raise TypeError('data parameter has an unexpected type')

        # Use any filename & meta-data inherited from input DataFile unless
        # specified:
        if not filename:
            filename = self._filename
        if meta is not None:
            self._meta = meta

        # Do likewise for component array labels, defaulting to the package
        # config values if unspecified and not copying another DataFile. Any
        # defaults not overridden still apply, to allow specifying a subset
        # without leaving things unlabelled.
        if not self._labels:
            self._labels = copy(config['labels'])
        if labels:
            if isinstance(labels, basestring):
                self._labels['data'] = labels
            elif hasattr(labels, 'keys'):
                self._labels.update(labels)
            else:
                raise TypeError('labels parameter has an unexpected type')

        # Parse any filename into a FileName object:
        self._filename = FileName(filename, strip=strip, prefix=prefix,
                                  suffix=suffix, dirname=dirname)

        # Load the file if requested and it contents haven't been overridden
        # and check that it exists or doesn't, to match expectations:
        exists = os.path.exists(str(self.filename))
        read_file = False
        if mode is None:
            mode = 'read' if str(self.filename) else 'new'
        if mode in ['read', 'update']:
            read_file = True
            if not exists:
                raise IOError('%s not found' % str(self.filename))
        elif mode == 'new':
            if exists:
                raise IOError('%s already exists' % str(self.filename))
        elif mode != 'overwrite':
            raise ValueError('unrecognized file mode, \'%s\'' % mode)
        self._mode = mode

        if read_file and self.meta is None:
            self._load_meta()
        if read_file and data is None:
            self._load_data()

        # Ensure (meta-)data attributes have the right types & track how many
        # NDData objects this instance contains. An empty DataFile has length
        # zero.
        if self._data is None:
            self._data = []
        if self._tables is None:
            self._tables = []
        if self.meta is None:
            self._meta = OrderedDict()
        if self._cals is None:
            self._cals = {}
        self._len = len(self._data)

        # Has the user overridden or accessed the file contents since
        # instantiation (if not, we can avoid saving without checking hashes
        # when supplying a DataFileList to an external program such as an
        # IRAF task via run_task)?
        self._unloaded = meta is None and data is None

    @property
    def filename(self):
        return self._filename

    # Re-parse any change of filename after instantiation. To override prefixes
    # etc., the user can supply a FileName instance as the argument. Also
    # update the mode to 'new' so that previously-read-only files may now be
    # saved if appropriate, without inadvertently clobbering any existing copy.
    @filename.setter
    def filename(self, value):
        if value != self._filename:
            self._filename = FileName(value)
            self._mode = 'new'

    @property
    def mode(self):
        # TO DO: make this update mode first if the filename attributes have
        # changed, eg. by storing a string of whatever it was last, otherwise
        # the mode doesn't change to 'new' as when replacing the filename.
        return self._mode

    @property
    def meta(self):
        self._unloaded = False
        return self._meta

    @property
    def cals(self):
        return self._cals

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
            return self._data[self._n-1]

    # Allow subscripting this instance to get ~NDData:
    def __getitem__(self, key):
        return self._data[key]

    # Deletion of an item from the DataFile:
    # This is a simplistic implementation for now; it will probably need to
    # handle some cleaning up etc.
    def __delitem__(self, key):
        del self._data[key]  # this can also be passed a slice
        self._len = len(self._data)

    # When printing the instance, show the filename:
    def __str__(self):
        return str(self.filename)

    # Need to fix this to do something sensible when there's no filename??
    # This produces an unquoted string in lists!
    def __repr__(self):
        return 'DataFile \'%s\' (len %d)' % (self.__str__(), self._len)

    # Append a new NDData instance to the DataFile:
    def append(self, item):
        # For the time being, restrict the items to append to NDData sub-
        # classes, otherwise quite strange types of arrays can result:
        if not isinstance(item, NDDataBase):
            raise ValueError('append argument should be NDData compatible')
        item = NDLater(data=item)
        # We pretty much have to recalculate the largest numeric identifier
        # here, since it would be fiddly to track when the user changes one.
        # Only auto-number new NDData if they don't already have identifiers.
        if item.ident is None:
            item.ident = self.next_ident
        self._data.append(item)
        self._len += 1
        self._unloaded = False

    # Extend the DataFile with the contents of another DataFile or list of
    # NDData-derived objects (or single object, like append). If the input is
    # another DataFile instance, any existing filename is dropped.
    def extend(self, items):
        if not isinstance(items, DataFile):
            items = DataFile(data=items)
        # Only auto-number new NDData if they don't already have identifiers.
        if all([ndd.ident is None for ndd in items]):
            for ident, item in enumerate(items, start=self.next_ident):
                item.ident = ident
        self._data += items._data
        # TO DO: re-number table idents if there's duplication between the dfs?
        # Maybe need "discard=True" option to eliminate duplicates?
        self._tables += items._tables
        self._len = len(self._data)
        self._unloaded = False

    def __len__(self):
        return self._len

    # Lazily-load the data from file. This is not closely tied to FITS but only
    # flat lists of NDData objects are supported, rather than any arbitrary
    # hierarchy supported by, say, HDF5.
    def _load_data(self):
        data_maps, table_maps = ndmio.map_file(self.filename,
                                               labels=self._labels)
        self._data = [NDLater(iomap=iomap) for iomap in data_maps]
        # Table proxy objects are kept directly in DataFile, rather than used
        # by a sub-class as for NDData, since there is little need to load
        # separate parts of a Table on demand (& it could be complicated):
        self._tables = table_maps

    def _load_meta(self):
        self._meta = ndmio.load_common_meta(self.filename)

    def reload(self):
        """
        Re-load NDData instances & shared meta-data from the associated file
        on disk (eg. to synchronize the DataFile instance with any changes made
        by external programs such as IRAF).

        When instantiating a new DataFile object, it is unnecessary to run
        reload() afterwards if the associated file already exists; it will be
        read automatically.

        Note that data arrays are not actually copied into memory here; they
        are re-mapped and still lazily-loaded once referenced (if applicable).
        """
        if not str(self.filename):
            raise IOError('Attempt to re-load DataFile object with no ' \
                          'associated file')

        # Currently if the user overwrites data in memory by reloading the
        # file, that's tough luck; safeguarding volatile memory is considered
        # less critical than data on disk.

        # If the file doesn't exist etc., just pass through the IOError from
        # PyFITS, without masking the origin of any more obscure errors:
        self._load_meta()
        self._load_data()
        self._len = len(self._data)

        # If an unsaved file with 'new' or 'overwrite' mode is reloaded after
        # some external process (IRAF) creates it, change the mode to 'update',
        # to reflect our copy now being based on the "hard" copy, just as at
        # instantiation:
        if self.mode in ['new', 'overwrite']:
            self._mode = 'update'

        # Since everything has been re-mapped from file here, we can reset
        # the unloaded flag (unlike when saving, where the user may still hold
        # a reference to the thing that was saved):
        self._unloaded = True

    def save(self):
        """
        Save NDData instances & common meta-data to the associated file,
        creating it if it doesn't already exist.

        """
        # Consider saving the mode in self so we can re-check here that the
        # file (non-)existence is still as expected.

        # Disallow overwriting stuff that the user didn't originally expect to.
        # This can be circumvented by changing the filename or casting to a
        # new DataFile copy with a different mode.
        if self.mode == 'read':
            raise IOError('attempted to save {0} with mode \'read\''\
                          .format(str(self.filename)))
        elif self.mode == 'new' and os.path.exists(str(self.filename)):
            raise IOError('file {0} with mode \'new\' would now overwrite an '\
                          'existing copy'.format(str(self.filename)))

        # This code should be made to append None values for data groups that
        # have not changed since last loaded or saved, which the back-end will
        # then skip re-writing if it can. However, there are various ways for
        # this to backfire if not done very carefully (aka. premature
        # optimization) so this is left as an exercise for later.

        data_label = self._labels['data']
        uncertainty_label = self._labels['uncertainty']
        flags_label = self._labels['flags']

        # Construct flat lists for the array, meta & (name, ident) tuple,
        # to pass to the save_list function. Also record the file location
        # index for each saved attribute, to allow remapping to the new file.

        # Currently it's assumed that image arrays & tables go in the same flat
        # list of arrays, which might be invalid for formats other than FITS.

        data_list, meta_list, identifiers, type_list = [], [], [], []
        imapidx, tmapidx = [], []

        idx = 0

        # Add any tables at the beginning of the file, to keep the NDData
        # arrays together if more get appended later:
        for tproxy in self._tables:
            idx += 1
            data_list.append(np.array(tproxy.table))
            meta_list.append(tproxy.table.meta)
            identifiers.append((tproxy.label, tproxy.ident))
            type_list.append('table')
            tmapidx.append(idx)

        # Add the NDData component arrays:
        for ndd in self._data:

            ident = ndd.ident

            arr_group = (ndd.data,
                         ndd.uncertainty.array**2 if ndd.uncertainty else None,
                         ndd.flags)

            meta_group = (ndd.meta, None, None)

            id_group = ((data_label, ident), (uncertainty_label, ident),
                        (flags_label, ident))

            type_group = ('image',) * 3

            # Include list entries for the main data array and only non-empty
            # uncertainty/flags (passing None values to save_list for those
            # would cause any existing information to be preserved at the
            # applicable location in the file, which isn't what we want).
            for arr, meta, arr_id, arr_type in \
                zip(arr_group, meta_group, id_group, type_group):

                if arr is not None or arr_id[0] == data_label:
                    idx += 1
                    data_list.append(arr)
                    meta_list.append(meta)
                    identifiers.append(arr_id)
                    type_list.append(arr_type)
                    imapidx.append(idx)
                else:
                    imapidx.append(None)

        ndmio.save_list(self.filename, data_list, meta_list, identifiers,
                        type_list, self.meta)

        # If the save succeeded without raising an exception, remap each
        # Table proxy & each NDLater's _io attribute to the newly-saved file.

        for (n, tproxy), idx in zip(enumerate(self._tables), tmapidx):
            self._tables[n] = TabMapIO(self.filename, idx=idx,
                                       label=tproxy.label, ident=tproxy.ident)

        for ndd, data_idx, uncertainty_idx, flags_idx in \
            zip(self._data, *[iter(imapidx)]*3):

            # Initialize a new _io instance in case it doesn't exist already:
            ndd._io = NDMapIO(FileName(self.filename),
                              ident=ndd.ident, data_idx=data_idx,
                              uncertainty_idx=uncertainty_idx,
                              flags_idx=flags_idx)

        # If the file mode was 'new', it needs changing to 'update' now it
        # has been saved, to allow saving further changes; likewise for
        # 'overwrite', to reflect any subsequent saves being based on the
        # existing file:
        if self.mode in ['new', 'overwrite']:
            self._mode = 'update'

    @property
    def unloaded(self):
        # To qualify as unloaded, the DataFile itself (ie. meta) has not to
        # have been touched, nor each constituent NDLater instance nor the
        # filename (WRT what NDLater is lazy-loading):
        return self._unloaded and all([ndd.unloaded and \
            ndd._io.filename == self.filename for ndd in self])

    @property
    def next_ident(self):
        """
        Return the next integer identifier greater than any existing integer
        identifier(s) or None if one or more existing identifier is undefined.
        """
        idents = [as_int_or_none(ndd.ident) for ndd in self]
        if not idents:
            ident = 1
        elif any([val is None for val in idents]):
            ident = None  # don't try to set ids if they're incomplete anyway
        else:
            ident = max(idents)
            ident = ident + 1 if ident else 1  # don't rely on False==0 FWIW
        return ident

    def renumber(self, idents=None):
        """
        Re-number/name the (ident attributes of) member NDData instances,
        either with user-supplied values or sequentially.

        Parameters
        ----------

        idents : list of int, optional
            The identifiers to use, if not numbering sequentially from 1.
            The intention is also to allow a list of str later (see NDLater).

        """

        if idents:
            if not isinstance(idents, list) or len(idents) != len(self):
                raise ValueError('idents must be a list matching the '\
                                 'DataFile length')
            pairs = zip(idents, self)
        else:
            pairs = enumerate(self, start=1)

        for ident, ndd in pairs:
            ndd.ident = ident

        self._unloaded = False


def as_int_or_none(val):
    """
    Convert int or str(int) to an integer, preserving None values and returning
    False for other types.
    """
    if val is None or isinstance(val, (int, long)):
        result = val
    elif isinstance(val, basestring):
        try:
            result = int(val)
        except ValueError, TypeError:
            result = False
    else:
        result = False
    return result


class DataFileList(list):
    """
    A class that holds a list of DataFile objects, tracking filenames
    and/or NDData collections with ancillary information. This implementation 
    is pretty much a normal Python list but provides a more convenient
    interface for instantiating multiple DataFile objects from a list of
    filenames or nddata instances/lists. 

    Parameters
    ----------

    filenames : str or list of str, optional
        The filename(s) on disk of the dataset(s) to be represented. There
        should either be as many filenames as data or None (but lists of
        NDData can be nested to associate subsets with fewer filenames).

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

    mode : str, optional
        'read' (default), 'new', 'update' or 'overwrite'
        Specifies whether the file should exist on disk already and be used to
        initialize this DataFile (if a filename is provided) and whether it can
        later be written to disk (also see DataFile). With 'read' and 'update',
        the specified file must already exist, with 'new', it must not exist
        and with 'overwrite', any existing file is ignored and will be replaced
        when writing to disk. The 'data' and 'filename' parameters always
        override whatever would otherwise be read from disk. If mode is None,
        it will be taken from any DataFile instances given as data, as long as
        they are all the same (raising an exception if not), defaulting to
        'read' if no DataFile instances are provided.

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

    def __init__(self, filenames=None, data=None, meta=None, mode=None,
        strip=False, prefix=None, suffix=None, dirname=None):

        # Check our args & if needed expand them out to match:
        filenames, data, meta, mode = self._expand_args(filenames=filenames,
            data=data, meta=meta, mode=mode, strip=strip, prefix=prefix,
            suffix=suffix, dirname=dirname)

        # If no mods are specified to the input DataFiles, use them directly:
        if filenames is None:  # means _expand_args just passed through data
            initlist = data
        # Otherwise, cast each data input to new DataFile with requested mods:
        else:
            initlist = [DataFile(filename=fn, data=obj, meta=mdict,
                mode=mode, strip=strip, prefix=prefix, suffix=suffix,
                dirname=dirname) for obj, mdict, fn in \
                zip(data, meta, filenames)]

        # Do whatever initialization a list object normally does:
        list.__init__(self, initlist)

        # Record the specified mode, to apply to any DataFiles added later:
        self._mode = mode

    def _expand_args(self, filenames, data, meta, mode, strip, prefix,
        suffix, dirname):

        # First convert filenames/data/meta to lists (if they aren't already)
        # or None. Although instantiating a list of things, single objects are
        # accepted for convenience and are expanded to match the lengths of
        # the other arguments if needed (except filenames). There isn't really
        # a simple way of doing this with a function because for some arguments
        # we want to distinguish specific list-like objects from a list of
        # those objects, for others, we want to distinguish duck-typed sequence
        # objects from a list etc. -- ie. the criteria vary. Also, it seems
        # safer to require specific input types in some cases and warn the user
        # if something unexpected is received. Handling all this is a little
        # fiddly but should make things conceptually simpler and/or more
        # readable from a user perspective.

        # Data can be a singly-nested list, with one sub-list per DataFile.
        # These cases are enumerated to distinguish them from container lists.
        if _compatible_data_obj(data) or data == []:
            data = [data]
        elif data is not None and not hasattr(data, '__iter__'):
            # If the constituent elements of the list don't have the right
            # type, let DataFile complain rather than checking here.
            raise TypeError('data parameter has an unexpected type')
        # How many data items do we have (None or 1 or length of list >=0)?
        len_data = seqlen(data)

        if isinstance(filenames, basestring):
            filenames = [filenames]
        elif filenames is not None and not hasattr(filenames, '__iter__'):
            raise TypeError('filenames parameter has an unexpected type')
        len_fn = seqlen(filenames)

        if hasattr(meta, 'keys'):  # dict-like (inc. PyFITS headers)
            meta = [meta]
        elif meta is not None and not hasattr(meta, '__iter__'):
            raise TypeError('meta parameter has an unexpected type')
        len_meta = seqlen(meta)

        # If the mode is undefined, use any existing DataFile modes if they're
        # all the same, otherwise default to 'read' as for DataFile.
        if mode is None:
            if isinstance(data, list):
                dfmodes = list(set([item.mode for item in data \
                                    if isinstance(item, DataFile)]))
                nmodes = len(dfmodes)
                if nmodes == 0:
                    mode = 'read'
                elif nmodes == 1:
                    mode = dfmodes[0]
                else:
                    raise ValueError('must specify DataFileList mode when ' \
                        'provided mixed DataFile modes')
            else:
                mode = 'read'

        # Determine whether the filename is being modified, to help decide
        # below whether a new copy of the input is needed:
        fn_modified = filenames or strip or prefix or suffix \
            or dirname is not None

        # If there are any modifiers to the input filename, meta or mode or
        # the data aren't already DataFiles, create a new copy of each DataFile
        # to hold those modifications without affecting existing objects
        # (otherwise, use what we were given directly, to allow making new
        # lists of existing DataFile instances):
        if not (isinstance(data, list) and meta is None and not fn_modified \
           and all([isinstance(item, DataFile) and item.mode == mode \
                    for item in data])):

            # Here we should have an existing DataFileList instance, a list of
            # DataFile objects, list of NDData or list of lists of NDData:

            # Expand out data, filenames & meta to lists of the same length,
            # to avoid repetition below.
            listlen = max(len_fn, len_data, len_meta)

            if listlen is None:
                filenames, data, meta = [], [], []
            else:
                listrange = range(listlen)
                # In the special case of filenames, we don't expand out a
                # single value to match the other lists, as any filenames are
                # expected to be unique within any given list:
                filenames = filenames if filenames else \
                            [None for item in listrange]
                data = data if len_data > 1 \
                       else [data[0] for item in listrange] if data \
                       else [data for item in listrange]
                       # This last line can produce [None] or [[]], the latter
                       # of which overrides any existing data in the DataFile.
                meta = meta if len_meta > 1 \
                       else [meta[0] for item in listrange] if meta \
                       else [meta for item in listrange]

            if not len(filenames) == len(data) == len(meta):
                raise ValueError('filenames, data & meta args are unmatched ' \
                    'in length')

        return filenames, data, meta, mode

    # Ensure the data to be added have a compatible mode. These rules might
    # still need a bit of tweaking but should be reasonable.
    def _check_mode(self, filename, data, strip, prefix, suffix, dirname):

        is_datafile = isinstance(data, DataFile)
        orig_mode = data.mode if is_datafile else self._mode

        if orig_mode == self._mode:
            return True

        # Disallow appending as-yet-unwritten files when the files in the list
        # are also expected to be on disk (eg. by IRAF):
        if self._mode in ['read', 'update'] and \
           orig_mode in ['new', 'overwrite']:
            raise ValueError('can\'t add unsaved file to DataFileList '\
                             'whose mode=\'{0}\''.format(self._mode))

        # When listing 'new' files, figure out the filename that would be added
        # and ensure it doesn't already exist on disk (simply instantiating
        # DataFile would produce a less clear error for this scenario).
        elif self._mode == 'new':
            fn = FileName(filename if filename else data.filename \
                          if is_datafile else None, strip=strip, prefix=prefix,
                          suffix=suffix, dirname=dirname)
            if os.path.exists(str(fn)):
                raise ValueError('can\'t add already-existing file(s) to '\
                                 'DataFileList whose mode=\'new\'')

        # Disallow overwriting files previously declared read-only by
        # accidentally appending them to a writeable list:
        elif orig_mode == 'read':  # list mode 'update'/'overwrite'
            raise ValueError('can\'t add read-only data to DataFileList '\
                             'whose mode=\'{0}\''.format(self._mode))

        return False

    # Wrap the normal list append:
    def append(self, data=None, meta=None, filename=None, strip=False, \
        prefix=None, suffix=None, dirname=None):

        # If the filename of the DataSet being appended is unmodified and its
        # existing mode is the same as the list's, we append that instance to
        # allow keeping the same DataFiles in multiple lists, otherwise we cast
        # to a new instance with the appropriate attributes modified. All
        # DataFiles inherit the file mode of the DataFileList. To avoid
        # surprises, the modes of existing DataFiles can only be overridden
        # arbitrarily at instantiation, when done explicitly, whereas here
        # certain combinations are disallowed.

        same_mode = self._check_mode(filename=filename, data=data, strip=strip,
            prefix=prefix, suffix=suffix, dirname=dirname)

        if isinstance(data, DataFile) and same_mode and not (meta or filename \
            or strip or prefix or suffix or dirname is not None):
            newdf = data
        else:
            newdf = DataFile(filename=filename, data=data, meta=meta,
                mode=self._mode, strip=strip, prefix=prefix, suffix=suffix,
                dirname=dirname)

        list.append(self, newdf)

    def extend(self, data=None, meta=None, filenames=None, strip=False, \
        prefix=None, suffix=None, dirname=None):

        # Check & if needed expand the argument lists, as in init.
        filenames, data, meta = self._expand_args(filenames=filenames,
            data=data, meta=meta, mode=self._mode, strip=strip, prefix=prefix,
            suffix=suffix, dirname=dirname)

        # If no mods are specified to the input DataFiles, use them directly:
        if filenames is None:
            initlist = data

        # Otherwise, cast each data input to new DataFile with requested mods:
        else:
            # This is a bit of a hack to run _check_mode as part of a list
            # comprehension; we only really care whether it raises an exception
            # but comparing its return value with None (which never happens)
            # conforms with the applicable syntax:
            initlist = [DataFile(filename=fn, data=obj, meta=mdict,
                mode=self._mode, strip=strip, prefix=prefix, suffix=suffix,
                dirname=dirname) for obj, mdict, fn in \
                zip(data, meta, filenames) \
                if self._check_mode(filename=fn, data=obj, strip=strip,
                    prefix=prefix, suffix=suffix, dirname=dirname) is not None]

        # Call the usual list.extend() method to finish the job:
        list.extend(self, initlist)


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

    ident : `int` (`str` to be supported later), optional
        File-independent identifier for this NDLater instance (eg. MOS slit
        number, CCD/amplifier number or object name). This is used to
        determine the correspondence of NDLater instances across multiple
        DataFile objects. [The use of string identifiers would currently
        cause incompatibility with IRAF run_task, until more bookkeeping is
        added to map them to numeric EXTVERs for the FITS kernel; also, the
        back-end loader/saver would also need modifications to handle them.]

    iomap : `NDMapIO`, optional
        An object that maps the data, uncertainty & flags attributes to a file
        name and indices within that file, enabling those attribute values to
        be lazily loaded. Internally, this iomap remains pointing to the
        existing file location with which the NDLater instance was last
        synchronized (loaded/saved), independently of whether the host DataFile
        (if any) changes, ensuring that the expected data are loaded on demand
        (as long as the original file still exists). The iomap.ident attribute
        is therefore independent of NDLater.ident, to allow renaming/numbering.

    Attributes
    ----------

    ident : `int` (`str` to be supported later)
        As described above. When the NDLater instance is read from a FITS file
        via DataFile, this value defaults to EXTVER (unless the file was saved
        from a previous instance with the value overridden). This attribute is
        an API short cut to meta['NDM_ID'], where the value persists.

    (See NDDataArray doc string for other methods & attributes.
     This is a Work in progress, to support DataFile.)

    """

    _id_key = 'NDM_ID'

    # After implementing lazy loading with NDLater instead of loading
    # everything in _load_nddata_from_FITS, 11 tests now take ~2.3s total to
    # run instead of ~1.2s (I think with 2 tests that read 2 SCI exts each).

    # This is based on the NDData & NDDataArray __init__ but avoids referencing
    # array attributes here, instead storing an obj that knows how to get them.
    def __init__(self, data=None, uncertainty=None, flags=None, meta=None,
                 ident=None, iomap=None):

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
        # to avoid adding I/O logic in more places than necessary. We also need
        # the meta-data here, to determine things like ident, units & WCS.
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

        # Where the file is mapped from disk but data & meta aren't loaded/
        # accessed yet, this flag remains false and we can avoid having to save
        # again prior to running some external program (ie. IRAF/run_task) on a
        # DataFileList, without first checking hashes. This is a one-time
        # optimization, as we don't know later (even after saving) whether the
        # user/app could modify a prior ref to the data, but it supports the
        # common case where DataFileList is only used to list input filenames.
        self._unloaded = data is None and uncertainty is None and \
                         flags is None and meta is None

        # Set or override identifier in order of precedence: 1. user parameter,
        # 2. value persisted in meta-data, 3. format-native iomap value.
        # Changing the value read from disk resets the above unloaded flag:
        if self.ident is None and self._io:
            self._meta[self._id_key] = self._io.ident
        if ident is not None:
            self.ident = ident

    @property
    def data(self):
        if self._data is None and self._io:
            self._data = self._io.load_data()
            self._unloaded = False
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
            self._unloaded = False
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
            self._unloaded = False
        return self._parent.flags

    @flags.setter
    def flags(self, value):
        NDDataArray.flags.fset(self, value)

    @flags.deleter
    def flags(self):
        self.flags = None

    @property
    def meta(self):
        self._unloaded = False
        return self._meta

    @property
    def ident(self):
        val = self._meta.get(self._id_key, None)
        return val or None  # PyFITS encodes None as ''; no need to distinguish

    @ident.setter
    def ident(self, value):
        if value != self.ident:
            self._unloaded = False
        self._meta[self._id_key] = value

    @property
    def unloaded(self):
        return self._unloaded


def load_file_list(filename):
    """
    Load a text file containing a list of filenames (or other strings) as a
    Python list.

    To obtain a list of DataFile objects, the result can easily be converted
    as in the example:

    >>> raw_biases = DataFileList(load_file_list('list_of_biases.txt'))

    If the listed files need downloading first, the usage would be similar to:

    >>> bias_list = load_file_list('list_of_biases.txt')
    >>> download_files(bias_list, server='gemini', dirname='raw')
    >>> raw_biases = DataFileList(bias_list, dirname='raw')

    (or it may be preferable to produce the initial list by other means, such
    as command-line arguments or a list definition in the user script).

    The DataFileList object can subsequently be used in place of the initial
    plain Python list.


    Parameters
    ----------

    filename : str
        Name of a plain-text file, containing one entry per line. Although
        the intention is mainly to work with filenames, any non-comment strings
        are valid. Lines whose first non-whitespace character is '#' are
        treated as comments.

    Returns
    -------

    list of str
        A list of filenames (or other strings), one per input line with any
        leading or trailing whitespace removed.

    """

    f = open(filename, 'r')

    flist = []

    # Append one filename per line to the list:
    for line in f:
        line = line.strip()         # remove new lines & trailing/leading space
        if line and line[0] != '#': # ignore empty lines & comments
            flist.append(line)

    f.close()

    return flist


def temp_saved_datafile(datafile):
    """
    Save a copy of a DataFile instance using a temporary filename, eg. for
    use by an external program, and return the copy. It is the caller's
    responsibility to delete the file once it is no longer needed.

    Although mapped to different files, both objects share the same pixel data
    until reloaded.
    """
    # Python doesn't provide a (non-deprecated) way to produce a temporary
    # filename without actually creating and opening the file (to avoid
    # possible race conditions & exploits). One can, however, let Python close
    # the file again and then recycle its name, saving the corresponding
    # DataFile immediately to avoid possible collisions.
    with tempfile.NamedTemporaryFile(
        prefix='tmp_{0}_'.format(datafile.filename.base),
        suffix='.'+datafile.filename.ext, dir='') as tmpfile:

        # Construct the new DataFile in memory from the old one:
        tdf = DataFile(data=datafile, filename=tmpfile.name, dirname='', \
                       mode='overwrite')

    # As soon as Python has closed its file handle, re-create the file by
    # saving the new DataFile object:
    tdf.save()

    # Pass the new temporary DataFile back to the caller (should it
    # automatically delete its own file when it goes out of scope?):
    return tdf


# To do:
# - Is the DataFile init logic needlessly re-reading any DataFile passed
#   as an argument? More specifically, I think this will be triggered when
#   adding a DataFile to a DataFileList with mode='read'.
# - Implement deepcopy methods?

