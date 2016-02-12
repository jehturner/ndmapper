# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

"""
A module for executing Web queries on specific data archive services, for file
downloads, calibration file matching etc.
"""

import os, os.path
import sys
from contextlib import closing
from StringIO import StringIO
import urllib2
import xml.dom.minidom as xmd
import tarfile
import hashlib

from .utils import to_filename_strings
from .calibrations import init_cal_dict, save_cal_dict, recurse_file_cals

__all__ = ['look_up_cals', 'look_up_single_cal_gemini', 'download_files',
           'download_files_gemini', 'download_query_gemini']


def _get_back_end(server, fn_base):
    """
    Private function to get the appropriate back-end for performing a look-up
    on a specified server, given the back-end function's base name. If the
    first argument is a user-supplied function rather than a server name, it
    is used directly and the second argument is ignored.
    """

    if callable(server):
        back_end_fn = server
    else:
        try:
            back_end_fn = eval('{0}_{1}'.format(fn_base, server))
        except NameError:
            raise ValueError('unknown server \'{0}\' for {1}'\
                             .format(server, fn_base))

    return back_end_fn


def look_up_cals(filenames, dependencies, server, cache=None,
                 obs_type='target'):
    """
    Look up specified filenames on a specified server and return a dictionary
    of best-matching calibration files and their calibrations in turn. 

    WARNING: According to the Python docs, the minidom library used in
    back-ends for this should not be used to parse input from untrusted URLs.

    Parameters
    ----------

    filenames : list of (convertible to str) or convertible to str
        The filename(s) for which calibration files are to be looked up.

    dependencies : dict of str : list of str
        Dictionary defining which types of observations (potentially) require
        which other types as calibrations. Both the keys and list values are
        strings such as 'bias' or 'flat', as recognized by the calibration
        server. The first key has the special name 'target'. This dictionary
        is normally imported with the name CAL_DEPS from the instrument mode
        sub-package used, eg. modes.gmos.spec.ifu.

    server : str, function
        Name of the service to be used for the look-up (which determines the
        back-end function selected from this library). The currently-available
        options are: 'gemini' (default). Alternatively, the user may supply a
        look-up function/object, which must accept a filename and calibration
        type name (such as 'flat') and return a list of (filename, checksum)
        tuples, where the checksum may be None.

    cache : str, optional
        Filename of a JSON cache file to load and save, if any. This is
        essentially an editable Python dictionary. Any existing entries in
        this file take precedence over new look-up results, ensuring
        reproducibility and allowing the user to override matches iteratively.
        If the file does not already exist, it will be created from a new
        look-up, while the default of None produces a new look-up every time.

    obs_type : str
        The type of observation, corresponding to `filenames`, for which
        calibrations are to be looked up (default 'target'). This parameter
        exists to allow explicitly looking up a subset of calibrations for
        calibrations if needed, but normally it should be left at its default
        value and all those dependencies will be found recursively anyway,
        starting from the 'target'.


    Returns
    -------

    dict
        A dictionary of calibration files & associations, in the format
        produced by calibrations.init_cal_dict().

    """
    # Convert each filename object to a string and remove any path &
    # processing prefix/suffixes to get the original base names:
    filenames = to_filename_strings(filenames, strip_names=True,
                                    strip_dirs=True, use_cal_dict=False)

    # Check that the dependency dict looks as expected:
    if not hasattr(dependencies, 'keys') or 'target' not in dependencies or \
       not all([isinstance(val, list) for val in dependencies.itervalues()]):
        raise ValueError('malformed \'dependencies\' dict; see doc string')

    if obs_type not in dependencies:
        raise ValueError('obs_type \'{0}\' not found in dependencies dict'\
                         .format(obs_type))
    
    # Determine which server-specific look-up function to use:
    back_end_fn = _get_back_end(server, 'look_up_single_cal')

    # Load the cached/edited calibration dictionary if one already exists,
    # otherwise initialize a new one:
    priorcache = cache if (cache and os.path.exists(str(cache))) else None
    cal_dict = init_cal_dict(filename=priorcache)

    # Look up the calibrations for each input filename:
    for filename in filenames:

        # Look up each of the calibration types defined in the dependencies
        # dict recursively:
        recurse_file_cals(filename, obs_type, dependencies, back_end_fn,
                          cal_dict)

    # Save the updated cal dict to the cache file:
    if cache:
        save_cal_dict(cal_dict, cache)

    return cal_dict


def look_up_single_cal_gemini(filename, cal_type):
    """
    Query the Gemini science archive for the named input file and return a
    list of (filename, checksum) tuples corresponding to the requested type
    of calibration files.
    """

    service = 'http://archive.gemini.edu/calmgr'

    query = '/'.join((service, cal_type, filename))

    # Perform the Web query and parse the page contents, printing any failed
    # query details in case of HTTP errors. Apparently urllib doesn't support
    # "with" directly, so create the context manager with contextlib.closing.
    try:
        with closing(urllib2.urlopen(query)) as fileobj:
            xml_dom = xmd.parse(fileobj)
    except urllib2.HTTPError:
        sys.stderr.write('Failed query: {0}'.format(query))  # to do: log this
        raise

    # Also looked at astropy.utils.xml.iterparser here but it doesn't seem to
    # parse things in a convenient hierarchical way for this purpose.

    matches = []
   
    # Iterate over the list of matching calibration files and generate a
    # list of filename, checksum tuples:
    for cal_item in xml_dom.getElementsByTagName('calibration'):

        xml_caltype = parse_xml_value(cal_item, 'caltype')

        # Only include cals of the requested type. Others can be encountered
        # eg. when querying a non-existent type, which returns all of them.
        if xml_caltype == cal_type:
            filename = parse_xml_value(cal_item, 'filename')
            checksum = parse_xml_value(cal_item, 'md5') \
                       if filename else None
            
            # To do: log a warning message that we failed to parse the
            # calibration entry if filename is None here.

            # Append the filename & md5sum values or None values if they could
            # not be parsed:
            matches.append((filename, checksum))

    return matches


def parse_xml_value(parent, tag):

    tag_entries = parent.getElementsByTagName(tag)

    # Return a value only if a single entry is found for the specified tag:
    if len(tag_entries) == 1:

        contents = tag_entries[0].childNodes

        # Only if there's a single value for the tag, return it. This can be
        # False if there's no text or multiple paragraphs between the opening
        # & closing tags etc.
        if len(contents) == 1:
            return contents[0].data


def download_files(filenames, server, dirname=''):
    """
    Download a list of files from the specified server, skipping any that
    already exist in the target directory.

    User authentication is not supported yet (ie. can't get proprietary data).

    Parameters
    ----------

    filenames : list of (convertible to str) or str or calibration dict
        The names of the files to request from the server. If a calibration
        dictionary (in the format produced by calibrations.init_cal_dict())
        is supplied, all the raw calibration files listed in it will be
        downloaded.

    server : str, function
        Name of the service to be used for the look-up (which determines the
        back-end function selected from this library). The currently-available
        options are: 'gemini' (default). Alternatively, the user may supply a
        look-up function/object, which must accept a list of filename strings
        and a directory name as arguments (and should implement skipping
        existing files if needed).

    dirname : str, optional
        Writeable directory in which to place the files. This path may be
        absolute or relative to the current directory. If it doesn't exist, it
        will be created.

    Returns
    -------

    list of str
        A list of the file(s) downloaded, for syntactical convenience, eg.
        allowing "DataFileList(download_files(some_list, server='gemini'))".
        This is just a copy of the input name string(s) that are derived from
        `filenames`, prefixed with `dirname`.

    """
    # Convert each filename to a string or calibration dictionary to a list
    # of files and remove any path & processing prefix/suffixes to get the
    # original base names.
    filenames = to_filename_strings(filenames, strip_names=False,
                                    strip_dirs=True, use_cal_dict=True)

    # Determine which server-specific look-up function to use:
    back_end_fn = _get_back_end(server, 'download_files')

    # Make sure the download location exists and is writeable:
    if dirname:
        if not os.path.isdir(dirname):
            os.mkdir(dirname)  # fails if a file exists with the same name
        if not os.access(dirname, os.W_OK):
            raise IOError('directory not writeable: {0}'.format(dirname))

    back_end_fn(filenames, dirname)

    return [os.path.join(dirname, fn) for fn in filenames]


def download_files_gemini(filenames, dirname=''):
    """
    Query the Gemini science archive for the named input files and download
    them to the specified directory, skipping any that already exist there.

    User authentication is not supported yet (ie. can't get proprietary data).

    Parameters
    ----------

    filenames : list of str
        The names of the files to request from the server. Any that are
        already present in `dirname` will be ignored.

    dirname : str, optional
        The (absolute or relative) directory path in which to place the files.

    """
    # (The alternative to this query is 'file', which omits the checksum)
    service = 'download'

    # The Gemini archive currently doesn't support looking up a list of known
    # filenames, so request them one at a time:
    for filename in filenames:

        # Get only those files that don't already exist:
        if not os.path.exists(os.path.join(dirname, filename)):

            query = '/'.join((service, filename))
            try:
                download_query_gemini(query, dirname)
            except:
                # Pass through the exception after reporting which was the
                # failed query:
                sys.stderr.write('Failed query: {0}\n'.format(query))
                raise


def download_query_gemini(query, dirname=''):
    """
    Perform a user-specified Gemini science archive query and save the files
    returned to a specified directory.

    Any existing files with the same names are overwritten, since there's no
    way to avoid downloading them all as part of the query anyway. However,
    only files that pass their checksum comparison (or have no checksum, which
    shouldn't happen) get written to disk. In the event that one or more files
    are corrupt, the caller can therefore re-try the query, omitting any files
    that were fetched successfully the first time if the query supports it.

    The downloaded tar file is stored in memory while processing its contents,
    which should be optimal as long as the archive isn't unreasonably large
    (to do: consider adding an option to write it to a temporary file).

    Parameters
    ----------

    query : str
        The query URL (or just the path component) to request from the server.

    dirname : str, optional
        The (absolute or relative) directory path in which to place the files.

    """
    url = 'http://archive.gemini.edu'
    checksum_fn = 'md5sums.txt'
    aux_fn = [checksum_fn, 'README.txt']

    # Accept either a full query URL or just the slash-separated criteria:
    query = query if query.startswith(url) else '/'.join((url, query))

    # Perform Web query and download the tar file to a StringIO file object
    # in memory, passing through any HTTP errors.
    with closing(urllib2.urlopen(query)) as fileobj:
        fobj_buff = StringIO(fileobj.read())

    # Open the in-memory tar file & extract its contents.
    with tarfile.open(fileobj=fobj_buff) as tar_obj:

        # Read the filenames & checksums for each file:
        chk_dict = {}
        with closing(tar_obj.extractfile(checksum_fn)) as chk_fobj:
            for line in chk_fobj:
                try:
                    checksum, fn = line.split()
                except ValueError:
                    raise ValueError('failed to parse {0}'\
                                     .format(checksum_fn))
                chk_dict[fn] = checksum

        # Get a list of the remaining files in the archive, without making
        # assumptions as to how the name(s) match the input list (eg. they
        # have a compression suffix).
        tarlist = [tinfo for tinfo in tar_obj.getmembers() \
                   if tinfo.name not in aux_fn]

        nosum, corrupt = [], []

        # Verify each data file and extract it to disk:
        for member in tarlist:

            fn = member.name

            # Extract the (compressed) file from the archive to memory.
            # This probably duplicates the memory usage momentarily for
            # each file in turn but that's unlikely to be critical.
            with closing(tar_obj.extractfile(fn)) as fobj:
                data_file = fobj.read()

            # Compare the checksums and record any problems. If one fails,
            # the whole query will need re-trying, so we'll just raise an
            # exception at the end and let the caller try again (ideally
            # excluding those files we've already got successfully here).
            if fn in chk_dict:
                checksum = hashlib.md5(data_file).hexdigest()
                if checksum == chk_dict[fn]:
                    good = True
                else:
                    corrupt.append(fn)
                    good = False
            else:
                nosum.append(fn)
                good = False

            # Write the file only if it didn't have a bad checksum:
            if good:
                decompress_to_disk(data_file, fn, dirname)

        if nosum:
            sys.stderr.write(
                'Files downloaded with missing checksums: {0}\n\n'\
                .format(' '.join(nosum)))

        if corrupt:
            raise IOError('Corrupt files skipped: {0}'\
                          .format(' '.join(corrupt)))


def decompress_to_disk(data, filename, dirname=''):
    """
    Write the contents of a file in memory to the named file on disk,
    decompressing it first if the file extension is recognized as a
    compressed data format. The file will be overwritten if it already exists.
    """

    outname, ext = os.path.splitext(filename)
    cmp_format = ext.lstrip(os.extsep)

    # If the file format is recognized, decompress the data before writing to
    # disk, otherwise just write the file as it is already.
    if cmp_format == 'bz2':
        import bz2
        data = bz2.decompress(data)

    elif cmp_format == 'gz':
        import gzip
        with gzip.GzipFile(fileobj=StringIO(data)) as gzip_obj:
            data = gzip_obj.read()

    else:
        # Restore the original file extension since we're not decompressing:
        outname += ext

    path = os.path.join(dirname, outname)

    # Write the file, clobbering any existing copy:
    with open(path, 'w') as fobj:
        fobj.write(data)

