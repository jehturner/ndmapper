# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os.path
from contextlib import closing
import urllib2
import xml.dom.minidom as xmd

from .utils import to_filename_strings
from .calibrations import init_cal_dict, save_cal_dict, recurse_file_cals

__all__ = ['look_up_cals', 'look_up_single_cal_gemini']


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
    filenames = to_filename_strings(filenames, strip=True, use_cal_dict=False)

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
        print 'Failed query: {0}'.format(query)  # to do: log this
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

