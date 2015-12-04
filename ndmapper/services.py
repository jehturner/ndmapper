# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os.path
from contextlib import closing
import urllib2
import xml.dom.minidom as xmd

from ndmapper.io import FileName
from ndmapper.data import DataFile
from .calibrations import init_cal_dict, recurse_file_cals

__all__ = ['look_up_cals', 'look_up_single_cal_gemini']


def look_up_cals(filenames, dependencies, server, cache=None,
                 obs_type='target'):
    """
    Look up specified filenames on a specified server and return a dictionary
    of best-matching calibration files and their calibrations in turn. 

    WARNING: According to the Python docs, the minidom library used in
    back-ends for this should not be used to parse input from untrusted URLs.

    Parameters
    ----------

    filenames : list of (str or DataFile or FileName)
                or str or DataFile or FileName
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
        [CACHE CURRENTLY UNIMPLEMENTED]

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

    # If given a single filename, make it into a list:
    if isinstance(filenames, (DataFile, FileName, basestring)):
        filenames = [filenames]
    elif not isinstance(filenames, list):
        raise ValueError('filenames parameter has an unexpected type')

    # Cast filenames to strings if they are DataFile instances and remove any
    # path & processing prefix/suffixes to get the original base name. This is
    # a bit convoluted to deal concisely with str & DataFile in the same way.
    # filenames = [FileName(str(fn)).base for fn in filenames]  # no: keep ext
    filenames = [str(FileName(str(fn), strip=True, dirname='')) \
                 for fn in filenames]

    # Check that the dependency dict looks as expected:
    if not hasattr(dependencies, 'keys') or 'target' not in dependencies or \
       not all([isinstance(val, list) for val in dependencies.itervalues()]):
        raise ValueError('malformed \'dependencies\' dict; see doc string')

    if obs_type not in dependencies:
        raise ValueError('obs_type \'{0}\' not found in dependencies dict'\
                         .format(obs_type))
    
    # Determine which server-specific look-up function to use:
    if callable(server):
        back_end_fn = server
    else:
        try:
            back_end_fn = eval('look_up_single_cal_'+server)
        except NameError:
            raise ValueError('unknown server \'{0}\''.format(server))

    # Load the cached/edited calibration dictionary if one already exists,
    # otherwise initialize a new one:
    priorcache = cache if os.path.exists(str(cache)) else None
    cal_dict = init_cal_dict(filename=priorcache)

    # Look up the calibrations for each input filename:
    for filename in filenames:

        # Look up each of the calibration types defined in the dependencies
        # dict recursively:
        recurse_file_cals(filename, obs_type, dependencies, back_end_fn,
                          cal_dict)

    # Save the dict to the cache here.

    return cal_dict


def look_up_single_cal_gemini(filename, cal_type):
    """
    Query the Gemini science archive for the named input file and return a
    list of (filename, checksum) tuples corresponding to the requested type
    of calibration files.
    """

    service = 'http://archive.gemini.edu/calmgr'

    query = '/'.join((service, cal_type, filename))

    # Perform the Web query and parse the page contents, passing back any HTTP
    # errors directly to the caller. Apparently urllib doesn't support "with"
    # directly so create the context manager with contextlib.closing.
    with closing(urllib2.urlopen(query)) as fileobj:
        xml_dom = xmd.parse(fileobj)

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

