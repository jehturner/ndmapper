# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.

# Describe as a list or collection or files/datasets/exposures/nddata/...?
# Result/node/IO etc.??

# Need to figure out the overlap with Mark's nddata container object at some
# point, but it seems to have no concept of file association and for now that
# looks more work than getting started with something simple.

# Need to understand the differences between AstroPy and Gemini Python
# docstring conventions better, which are supposedly the same but one has
# extra fields.

import os.path
import string
import re
from copy import deepcopy
import numpy as np
from astropy.nddata import NDDataBase, NDDataArray
from astropy.nddata import NDUncertainty, StdDevUncertainty
import astropy.io.fits as pyfits
from . import config


class FileName(object):
    """
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
        variable "quadpype.config['filename_regex']" to be used, which
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

    ext : list
        File extension, eg. "fits".

    sep : str, None
        One or more characters specified as a separator.

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
            froot, fext = os.path.splitext(os.path.basename(path))
            if fext:
                self.ext = fext.lstrip('.')
            else:
                self.ext = None   # distinguish no ext at all from a dot

            # # Separate base filename into a prefix+base & a final suffix:
            # elements = froot.rsplit('_', 1)
            # froot = elements[0]
            # if len(elements) > 1:
            #     self.suffix=elements[1]
            # else:
            #     self.suffix=''

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
            dotext = '.'+self.ext
        return (os.path.join(self.dir, (self.prefix+self.base+ \
                string.join(self.suffix, '')+dotext)))

    def __deepcopy__(self, memo):
        # One can't copy a regex, except by re.compile(_re.pattern), but they
        # look to be immutable anyway so this shouldn't be a problem; likewise
        # for strings:
        return FileName(path=str(self), sep=self.sep, regex=self._re)


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

    filename : str or FileName, optional
        The filename on disk of the dataset(s) to be represented.

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

    data : list of NDData or None
        NDData instance(s) associated with the DataFile.

    """

    filename = None
    log = ''

    _dirty = True    # has the file been modified since saved/loaded?


    def __init__(self, data=None, filename=None, strip=False, prefix=None,
        suffix=None, dirname=None):

        if isinstance(data, DataFile):  # create new copy of existing obj
            self.data = data.data
            self.filename = deepcopy(data.filename)
        elif isinstance(data, NDDataBase):
            self.data = [data]
        elif hasattr(data, '__iter__'):  # sequence, presumably of NDData (or
            self.data = data             # empty to ignore existing file data)
        elif data is None:
            self.data = None
        else:
            raise TypeError('DataFile: data parameter has an unexpected type')

        # Use any filename inherited from input DataFile unless specified:
        if not filename:
            filename = self.filename

        # Parse any filename into a FileName object:
        self.filename = FileName(filename, strip=strip, prefix=prefix,
            suffix=suffix, dirname=dirname)

        # Does the file already exist on disk?
        if data is None and os.path.exists(str(self.filename)):
            self._load()

        # TO DO: should also read the primary header here as an attribute
        # of the data file.

        # Track how many NDData objects this instance contains. An empty
        # DataFile has length zero.
        if self.data is None:
            self._len = 0
        else:
            self._len = len(self.data)            

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

    # Lazily-load? the (meta-)data from file. For the time being we assume
    # the data live in a FITS file, since the NDData I/O functionality
    # (NDIOMixin & astropy.io.registry) is still being fleshed out and it's
    # unclear how it would handle mapping multiple NDData objects from a
    # file -- but things are otherwise not structured to assume FITS, except
    # in that we currently support only flat lists of NDData objects rather
    # than any arbitrary hierarchy supported by, say, HDF5.
    def _load(self):
        self.data = list(_load_nddata_from_FITS(str(self.filename)))


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

    filenames : str or list of str, optional
        The filename(s) on disk of the dataset(s) to be represented. There
        should either be as many filenames as data or None (but lists of
        NDData can be nested to associate subsets with fewer filenames).

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

    def __init__(self, data=None, filenames=[], strip=False, prefix=None,
        suffix=None, dirname=None):

        # This is a little fiddly (not too much considering) but handling
        # these cases here should make things conceptually simpler and/or
        # more readable from a user perspective...

        data_len = seqlen(data)
        fn_len = seqlen(filenames)

        fn_modified = filenames or strip or prefix or suffix \
            or dirname is not None
        
        # If we got a single NDData or DataFile instance, create an
        # DataFile directly from it, overriding specified prefixes etc.:
        if isinstance(data, DataFile) or isinstance(data, NDDataBase):
            if fn_len > 1:
                raise ValueError('DataFileList: got multiple filenames for ' \
                    + 'a single data object')
            elif fn_len == 1:
                filename = filenames[0]
            else:
                filename = filenames

            # If instantiating with an unmodified DataFile (no filename prefix
            # changes etc.), just reference the existing instance by default,
            # instead of making a copy, to allow making new lists of existing
            # DataFile instances:
            if isinstance(data, DataFile) and not fn_modified:
                initlist = [data]
            else:
                initlist = [DataFile(data=data, filename=filename, strip=strip,
                    prefix=prefix, suffix=suffix, dirname=dirname)]

        # If we got an existing DataFileList instance, re-create it from its
        # member DataFile objects to apply any override params. Likewise, if
        # we got a list of DataFile objects, list of NDData or list of lists
        # of NDData, construct a DataFileList from that.
        else:
            # If instantiating with a list or DataFileList of unmodified
            # DataFiles (no filename prefix changes etc.), just reference
            # the existing instances by default, instead of creating copies,
            # to allow making new lists of existing DataFile instances:
            if isinstance(data, list) \
                and all([isinstance(item, DataFile) for item in data]) \
                and not fn_modified:
                initlist = data

            # Otherwise create a new copy of each DataFile to hold mods:
            else:
                # Should have got something iterable or None for data. Expand
                # out both data & filenames to avoid repetition below.
                if not data:
                    if fn_len is None:
                        filenames = [filenames]
                    data = [None for fn in filenames]
                # Filenames need to be unique so there must be None or as many
                # as there are items in data (but data can include singly-
                # nested lists of NDData, with each sub-list associated with
                # one filename and one resulting DataFile).
                elif not filenames:
                    filenames = [None for obj in data]
                elif fn_len != data_len:
                    raise ValueError('DataFileList: data & filenames differ ' \
                        'in length')
                initlist = [DataFile(data=obj, filename=fn, strip=strip,
                    prefix=prefix, suffix=suffix, dirname=dirname) \
                    for obj, fn in zip(data, filenames)]

        # Do whatever initialization a list object normally does:
        list.__init__(self, initlist)

    # Wrap the normal list append in the same way as __init__:
    def append(self, data=None, filename=None, strip=False, prefix=None,
        suffix=None, dirname=None):

        if isinstance(data, DataFile) and not (filename or strip or prefix or \
            suffix or dirname is not None):
            list.append(self, data)
        else:
            list.append(self, DataFile(data=data, filename=filename,
                strip=strip, prefix=prefix, suffix=suffix, dirname=dirname))

    # Wrap the normal list extend in the same way as __init__:
    def extend(self, data=None, filenames=[], strip=False, prefix=None,
        suffix=None, dirname=None):

        list.extend(self, DataFileList(data=data, filenames=filenames, 
            strip=strip, prefix=prefix, suffix=suffix, dirname=dirname))


def seqlen(arg):
    """
    Return the argument's length only if it's a sequence, otherwise None.
    """
    if hasattr(arg, '__iter__'):    
        try:
            slen = len(arg)
        except TypeError:
            slen = None
    else:
        slen = None
    return slen
    

def _load_nddata_from_FITS(filename):
    """
    Open an existing FITS file and return the corresponding NDDataArray
    instances.

    TO DO: populate the mask attribute from flags.
    """

    hdulist = pyfits.open(filename, mode='readonly')

    # A dict of empty lists to sort recognized extensions into:
    ext_dict = {'data' : [], 'uncertainty' : [], 'flags' : [], 'undef' : []}

    # Sort any FITS image extensions by EXTNAME into SCI/VAR/DQ lists:
    have_names = False
    for hdu in hdulist:

        # Ignore any non-image extensions:
        if isinstance(hdu, pyfits.ImageHDU):

            # The name/ver attributes are semi-documented but seem to be
            # part of the public API now.
            if hdu.name:  # seems to default to ''; also works for None
                have_names = True

            if hdu.name == config['data_name']:
                ext_dict['data'].append(hdu)

            elif hdu.name == config['uncertainty_name']:
                ext_dict['uncertainty'].append(hdu)

            elif hdu.name == config['flags_name']:
                ext_dict['flags'].append(hdu)

            elif not hdu.name:
                ext_dict['undef'].append(hdu)

            # else:
            #     ignore any extensions with unrecognized names

            # print ext_dict

    # If there are no named image extensions, treat them as main data arrays:
    if not have_names:
        ext_dict['data'] = ext_dict['undef']

    # List existing data (SCI) array extension numbers for reference:
    extvers = [hdu.ver for hdu in ext_dict['data']]

    # Create the NDDataArray instances. Since the main data array is mandatory,
    # we just ignore any uncertainty (VAR) or flags (DQ) extensions without a
    # corresponding data (SCI) array and loop over the latter:
    lastver = 0
    ndlist = []
    for data_hdu in ext_dict['data']:

        # Give any unnumbered SCI extensions the next available EXTVER after
        # the last one used:
        if data_hdu.ver < 1:  # seems to default to -1; also works for None
            thisver = lastver + 1
            while thisver in extvers:
                thisver += 1
            data_hdu.ver = thisver
            uncert_hdu = None
            flags_hdu = None

        # Otherwise, if the EXTVER was defined to begin with, look for
        # associated uncertainty & flags (VAR/DQ) extensions:
        else:
            # Find uncertainty & flags HDUs matching this EXTVER:
            uncert_hdu = None
            for hdu in ext_dict['uncertainty']:
                if hdu.ver == data_hdu.ver:
                    uncert_hdu = hdu
                    break

            flags_hdu = None
            for hdu in ext_dict['flags']:
                if hdu.ver == data_hdu.ver:
                    flags_hdu = hdu
                    break

        lastver = data_hdu.ver

        if uncert_hdu:
            uncert_data = StdDevUncertainty(np.sqrt(uncert_hdu.data))
        else:
            uncert_data = None

        if flags_hdu:
            flags_data = flags_hdu.data
        else:
            flags_data = None

        # Instantiate the NDData instance:
        # TO DO: Fix this to use a VarUncertainty class instead of StdDev.
        ndlist.append(NDDataArray(data=data_hdu.data, uncertainty=uncert_data,
            mask=None, flags=flags_data, wcs=None, meta=data_hdu.header,
            unit=None))

    # We don't keep the file open continually, since it may get updated later
    # by IRAF or whatever (this means some trickery later on to keep io.fits
    # happy, since we haven't read in the lazy-loaded data arrays yet).
    hdulist.close()

    return ndlist


# To do:
# - In progress
#   - Add the primary header as a DataFile attribute.
#   - Add NDLater & NDFITSLoader/NDLoader classes for lazy-loaded NDData.
#     - Reasonably important for IRAF to avoid duplicating everything.
#   - Write a separate len method so data can be lazy-loaded??
#     - How to do lazy loading for NDData when main array is mandatory?
#       - Need a sub-class?? Replace an empty array?
#     - Ask Erik about it?
#   - Go back to looping over MEF files in iraf_task.
# - Moved onto iraf_task for now.
# - Implement deepcopy methods?
# - Start unit tests as soon as feasible.

