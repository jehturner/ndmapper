# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.

# Describe as a list or collection or files/datasets/exposures/nddata/...?
# Result/node/IO etc.??

# Need to figure out the overlap with Mark's nddata container object at some
# point, but it seems to have no concept of file association and for now that
# looks more work than getting started with something simple.

import os.path
import string
import re
from copy import deepcopy
from astropy.nddata import NDDataBase

class FileName(object):
    """
    Parameters
    ----------

    path : str, FileName
        Single filename to parse into a FileName object representation.

    regex : str, re
        Regular expression matching root filename (without a file extension;
        defaults to Gemini's "S20150101S0001"-style convention).

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

    re : re
        Regular expression used to match the base name.

    sep : str, None
        One or more characters specified as a separator.

    """

    def __init__(self, path=None, regex="[NS][0-9]{8}S[0-9]{3,4}", sep="_",
        strip=False, prefix=None, suffix=None, dirname=None):

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
            # file, so that (eg.) NDDataFile objects instantiated with None
            # won't get some anomalous default filename.
            self.dir = ''
            self.prefix = ''
            self.base = ''
            self.suffix = ''
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
                self.base = match.group()
                if strip:
                    self.prefix = ''
                    self.suffix = ['']
                else:
                    self.prefix = froot[:match.start()]
                    self.suffix = self._split(froot[match.end():])
            else:
                raise ValueError('failed to parse base name with regexp "%s"\n'
                    '            from %s' % (regex, path))

        # Add on any specified prefix or suffix:
        if prefix:
            self.prefix = prefix + self.prefix
        if suffix:
            self.suffix.extend(self._split(suffix))

        # Add or replace any initial directory name if specified:
        if dirname:
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
        return FileName(path=str(self), regex=self._re, sep=self.sep)

class NDDataFile(object):
    """
    A class representing a data-file-like object, including a filename and/or
    a list of associated NDData instances, along with any ancillary data such
    as a primary header or binary tables that describe the file as a whole and
    aren't part of the NDData structure itself.
    
    This class can simply store a collection of associated NDData instances in
    memory or it can simply store a filename, eg. for use with IRAF tasks that
    do disk-based I/O, or it can do both, handling any loading & saving. In
    principle this abstraction allows mixing Python code operating on NDData
    fairly seamlessly with steps that do disk I/O (eg. IRAF tasks), as well
    as providing a convenient file-like way of organizing an NDData collection.

    Parameters
    ----------

    filename : str or FileName, optional
        The filename on disk of the dataset(s) to be represented.

    data : NDData or list of NDData or NDDataFile, optional
        NDData instance(s) for the dataset(s) to be represented (or an
        existing NDDataFile instance to copy deeply).

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
        NDData instance(s) associated with the NDDataFile.

    """

    filename = None

    def __init__(self, data=None, filename=None, strip=False, prefix=None,
        suffix=None, dirname=''):

        if isinstance(data, NDDataFile):  # create new copy of existing obj
            self.data = deepcopy(data.data)
            self.filename = deepcopy(data.filename)
        elif isinstance(data, NDDataBase):
            self.data = [data]
        elif data and hasattr(data, '__iter__'):  # non-empty sequence,
            self.data = data                      # presumably of NDData
        elif data is None:
            self.data = None
        else:
            raise TypeError('NDDataFile: data parameter has an unexpected type')

        # Use any filename inherited from input NDDataFile unless specified:
        if not filename:
            filename = self.filename

        # Parse any filename into a FileName object:
        self.filename = FileName(filename, strip=strip, prefix=prefix,
            suffix=suffix, dirname=dirname)

        # Track how many NDData objects this instance contains. An empty
        # NDDataFile has length zero.
        if self.data is None:
            self._len = 0
        else:
            self._len = len(self.data)            

    # Specify that class manages its own iteration with next() method:
    def __iter__(self):
        self._n = 0
        return self

    # Iteration over the NDDataFile instance returns ~NDData:
    def next(self):
        self._n += 1
        if self._n > self._len:
            raise StopIteration
        else:
            return self.data[self._n-1]

    # Allow subscripting this instance to get ~NDData:
    def __getitem__(self, key):
        return self.data[key]

    # Deletion of an item from the NDDataFile:
    # This is a simplistic implementation for now; it will probably need to
    # handle some cleaning up etc.
    def __delitem__(self, key):
        del self.data[key]

    # When printing the instance, show the filename:
    def __str__(self):
        return str(self.filename)

    # Facilitate appending new NDData instances to the NDDataFile, If the
    # input is another NDDataFile instance, any existing filename is dropped.
    def append(self, elements):
        if not isinstance(elements, NDDataFile):
            elements = NDDataFile(elements)
        if elements.data:
            if self.data is None:
                self.data = []
            self.data += elements.data
        if self.data is not None:
            self._len = len(self.data)

    def __len__(self):
        return self._len

class DataList(object):
    """
    A class that represents a list of files and/or corresponding NDData
    objects in memory. It can also track the associated processing history. 

    For the time being, all this does is track lists of filenames and their
    prefixes, for running IRAF tasks, but the idea is eventually to be able to
    convert transparently between files on disk and NDData instances through
    this representation, allowing IRAF & Python steps to be mixed freely
    (*if* they have compatible expectations regarding metadata).

    # Need to understand the differences between AstroPy and Gemini Python
    # docstring conventions better, which are supposedly the same.

    Parameters
    ----------

    filenames : str or list of str, optional
        The filename(s) on disk of the dataset(s) to be represented.

    data : NDData or list of NDData or DataList, optional
        NDData instance(s) for the dataset(s) to be represented (or an
        existing DataList instance to copy deeply).

    """

    filenames = []

    # Support initializing with DataList instance(s) (as data)?

    def __init__(self, data=None, filenames=[], strip=False, prefix=None,
        suffix=None, dirname=''):

        # Can initialize with an existing DataList object as data, creating
        # a new copy (can also use copy function but not deepcopy):
        if isinstance(data, DataList):
            self.data = deepcopy(data.data)
            self.filenames = deepcopy(data.filenames)
        elif data is not None:
            raise NotImplementedError('DataList: parameter data only accepts ' \
                'a DataList for now')
        else:
            self.data = None

        # Get the default regexp from FileName class and cache it for
        # subsequent instances (without having to define multiple defaults):
        self._re = FileName().re

        # If the argument was a single filename, convert to a list:
        if isinstance(filenames, basestring):
            filenames = [filenames]

        # Use any filenames inherited from an input DataList if not specified:
        if not filenames:
            filenames = self.filenames

        # Parse the list of file names into a list of FileName objects:
        self.filenames = [FileName(fn, regex=self._re, strip=strip, \
            prefix=prefix, suffix=suffix, dirname=dirname) for fn in filenames]

        # When we init with nddata instances here later, we should ensure that
        # filenames is populated with the corresponding number of entries, even
        # if they are empty strings or None. Consider using some default names
        # until specified.

        self._len = len(self.filenames)

    # Specify that class manages its own iteration with next() method:
    def __iter__(self):
        self._n = 0
        return self

    # Iteration over the DataList instance:
    # What should this return, in order to retain the file/nddata
    # association??? A single-element DataList or a Data Set
    def next(self):
        self._n += 1
        if self._n > self._len:
            raise StopIteration
        else:
            return DataList(filenames=str(self.filenames[self._n-1]))

    # Subscripting the DataList instance:
    def __getitem__(self, key):
        return DataList(filenames=str(self.filenames[key]))

    # Deletion of an item from the DataList:
    # This is a simplistic implementation for now; it will probably need to
    # handle some cleaning up etc.
    def __delitem__(self, key):
        del self.filenames[key]

    # Return the first (and only, when iterating) filename from the list:
    @property
    def name(self):
        if self._len > 0:
            return self.filenames[0]
        else:
            return None

    # When printing the DataList instance, show the filename(s), which is
    # what a human being commonly uses to distinguish the items:
    # Print basenames without any path prefix?
    def __str__(self):
        return str([str(fn) for fn in self.filenames])

    def append(self, elements):
        if isinstance(elements, DataList):
            #dl = DataList(elements)
            if elements.data:
                self.data += deepcopy(elements.data)
            if elements.filenames:
                self.filenames += deepcopy(elements.filenames)
            self._len = len(self.filenames)
        # Could extend this to accept the same parameters as __init__ and
        # instantiate a new DataList before appending it.
        else:
            raise NotImplementedError('DataList: append() expects a DataList')

    # def __deepcopy__(self, memo):
    #     return DataList()  # fix this

    def __len__(self):
        return self._len

# eg flats = thingmy, a list

# To do:
# - Created NDDataFile class from DataList.
#   - Modify & simplify DataList to use it.
#     - Do we still need a special DataList object?
#   - Make DataList return NDDataFile when subscripting etc.
# - Implement deepcopy methods?
#   - Append a filename or list of filenames or a DataList
#     - Instantiate a DataList and append that?
# - Make base class for FileName with re.??
#   - But how to call that from NDDataFile etc.?
#   - Want some way to avoid re-compiling the re and specifying it every
#     time something is instantiated for non-Gemini data, eg. a global
#     that the user can modify.
# - Start unit tests as soon as feasible.

