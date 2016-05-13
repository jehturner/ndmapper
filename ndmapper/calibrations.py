# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

"""
A module for management of calibration exposure dependencies & associations via
a Python dictionary (in conjuction with the "services" module).
"""

import os.path
import json

from ndmapper.libutils import splitext, addext
from ndmapper.data import FileName, DataFile, DataFileList

__all__ = ['init_cal_dict', 'save_cal_dict', 'add_cal_entry', 'cal_entries',
           'associate_cals']


K_ASSOCIATIONS = 'associations'
K_CALIBRATIONS = 'calibrations'
K_CHECKSUMS    = 'checksums'


def init_cal_dict(filename=None):
    """
    Load a calibration dictionary from a JSON file if one is specified,
    otherwise initialize a new dictionary in the same format.

    Parameters
    ----------

    filename : `str`, optional
        Name of the JSON cache file to load, if any.

    Returns
    -------

    `dict`
        An (empty or populated) dictionary of calibration files & associations,
        in the required format. This consists of 3 sub-dictionaries,
        'associations', mapping each input filename to a dictionary of
        calibration name labels, keyed by calibration type; 'calibrations',
        mapping each calibration name label to a list of constituent raw files
        and 'checksums', mapping each raw calibration file to an optional
        checksum. [To do: add an example.]

    """

    # To do: this should log a warning if the file doesn't exist, in case the
    # user has mistyped the name (there's no way to distinguish that from the
    # cache just not being created yet, though unintentional differences won't
    # ordinarily occur as both cases will involve the same line of user code).

    # If a cache file was specified and exists, read it:
    if filename and os.path.exists(filename):
        with open(filename) as fobj:
            jstr = fobj.read()

        # Convert JSON to Python dict, passing through any parsing errors:
        cal_dict = json.loads(jstr)

        # Catch any basic errors in the calibration dict format itself. This
        # could be made more bulletproof & informative but at least ensures
        # that the structures are nested as expected so look-ups will work:
        if isinstance(cal_dict, dict) and \
           sorted(cal_dict.iterkeys()) == \
               sorted([K_ASSOCIATIONS, K_CALIBRATIONS, K_CHECKSUMS]) and \
           all([isinstance(val, dict) \
                for val in cal_dict.itervalues()]) and \
           all([isinstance(val, dict) \
                for val in cal_dict[K_ASSOCIATIONS].itervalues()]) and \
           all([isinstance(val, list) \
                for val in cal_dict[K_CALIBRATIONS].itervalues()]):
            valid = True
        else:
            valid = False

        if not valid:
            raise ValueError('Bad calibration dict format in {0}'\
                             .format(filename))

    # Otherwise, initialize a new calibration dictionary:
    else:
        cal_dict = {
            K_ASSOCIATIONS : {},
            K_CALIBRATIONS : {},
            K_CHECKSUMS : {}
        }

    return cal_dict


def save_cal_dict(cal_dict, filename):
    """
    Save the provided calibration dictionary as a user-editable JSON file.

    Parameters
    ----------

    cal_dict : `dict`
        A dictionary of calibration files & associations, in the format
        produced by `init_cal_dict()`.

    filename : `str`
        Name of the JSON cache file to save the dictionary to.

    """
    # Create a pretty-printed JSON string that's easy to read & edit:
    jstr = json.dumps(cal_dict, sort_keys=True, indent=4)

    # Save it to the cache file, overwriting any existing copy:
    with open(filename, 'w') as fobj:
        fobj.write(jstr)


def recurse_file_cals(filename, obs_type, dependencies, lookup_fn, cal_dict):
    """
    Look up the calibrations on which a given file depends (including their
    calibrations in turn). This is a support function for the public user
    interface `services.look_up_cals()`.

    """

    # Get a reference to each sub-dictionary to avoid verbosity & repetition:
    associations = cal_dict[K_ASSOCIATIONS]
    calibrations = cal_dict[K_CALIBRATIONS]    

    # Loop over the cal types on which this observation type directly depends:
    for cal_type in dependencies[obs_type]:

        # Look up the applicable type of calibration for this file only if
        # there isn't already an association in the dictionary:
        if filename not in associations or \
           cal_type not in associations[filename]:

            # Look up the specified type of calibration for this file:
            matches = lookup_fn(filename, cal_type)

            # Populate sub-dicts with the retrieved filenames & checksums. If
            # no match was found, a placeholder None entry gets created:
            add_cal_entry(filename, cal_type, matches, cal_dict)

        # Now loop over the raw files for each matching calibration of the
        # applicable type and recursively look up any calibrations that they
        # need in turn. Typically, matches with their own dependencies will
        # either be single files / pairs (eg. a spectroscopic flat) or will
        # share the same dependencies (eg. a set of imaging flats with the
        # same bias), so the tree doesn't really grow exponentially.

        # Get the key for the group of files corresponding to this calibration:
        dep_label = associations[filename][cal_type]

        # Recursively call self on each of the files in the group (unless
        # the look-up failed to find any, in which case the label is None):
        if dep_label:
            for dep_file in calibrations[dep_label]:
                recurse_file_cals(dep_file, cal_type, dependencies, lookup_fn,
                                  cal_dict)

    return cal_dict


def add_cal_entry(filename, cal_type, matches, cal_dict):
    """
    Populate each component of the calibration dictionary with the specified
    list of calibration files matching a given filename.

    At this level, any existing association between the specified filename
    and the same type of calibration will be overwritten, but any existing
    lists of available calibration files and the corresponding checksums will
    be left unchanged. Normally, this function only gets called (eg. via
    `services.look_up_cals()`) for calibration entries that don't already exist
    in the dictionary, to preserve existing associations.

    Parameters
    ----------

    filename : `str`
        Name of the file for which the applicable calibration files are to
        be recorded.

    cal_type : `str`
        Type of calibration looked up (matching a name in the dependencies
        dictionary applicable to the instrument mode).

    matches : `list` of (`str`, `str` or `None`)
        List of (filename, checksum) pairs from a prior calibration look-up,
        to be included in the calibration dictionary.

    cal_dict : `dict`
        The calibration dictionary to populate, already initialized in the
        required format by `init_cal_dict()`.

    """

    # Complain if "matches" doesn't actually look like a list of matches.
    # Maybe this is more processing than one would habitually call from a
    # recursive loop but it should be fast compared with the look-up itself.
    if not hasattr(matches, '__iter__') or \
       not all([hasattr(match, '__getitem__') and len(match)==2 \
                for match in matches]):
        raise ValueError('matches should be a list of (filename, checksum)')

    # For now we assume the calibration dictionary is in the right format.
    # Checking that would involve writing a verification function.

    # Get a reference to each sub-dictionary to avoid verbosity & repetition:
    associations = cal_dict[K_ASSOCIATIONS]
    calibrations = cal_dict[K_CALIBRATIONS]
    checksums    = cal_dict[K_CHECKSUMS]

    # Add an entry for this filename in the calibration dict if not present:
    if filename not in associations:
        associations[filename] = {}

    # If there are no matches to add, just return a placeholder entry for
    # this calibration type, for the user to fill in manually:
    if not matches:
        associations[filename][cal_type] = None
        return

    # Sort the matches. This avoids duplicate lists if the same set of files
    # happens to be presented in a different order for multiple look-ups and
    # also reproduces Gemini's usual naming convention. Don't sort in place,
    # to avoid changing our input parameters.
    matches = sorted(matches)
    
    # Group calibration files under a label that's based on the first filename
    # in the sorted matches list; this label will later become the processed
    # calibration filename.
    ref_base, ref_ext = splitext(matches[0][0])
    label_base = '{0}_{1}'.format(ref_base, cal_type)
    label = addext(label_base, ref_ext)   # feasible to drop ext??

    # Extract the filenames from the matches:
    match_names = [match[0] for match in matches]

    # If there's already an entry for this calibration label in the dictionary
    # but the corresponding list of files is different, append a running number
    # to the label until it's either unique or does match an existing entry:
    n = 0
    while label in calibrations and calibrations[label] != match_names:
        n += 1
        label = addext('{0}_{1}'.format(label_base, n), ref_ext)

    # Add the list of calibration files to the dictionary if not already there:
    if label not in calibrations:
        calibrations[label] = match_names

    # Record each file's checksum, only if there is one. The condition here
    # ensures entries won't change, eg. if user has corrected a bad checksum:
    for key, val in matches:
        if val and key not in checksums:
            checksums[key] = val

    # Record the calibration association corresponding to this look-up:
    associations[filename][cal_type] = label


def cal_entries(cal_dict, cal_type, reference=None):
    """
    Extract entries of a specified type from the 'calibrations' section of
    a calibration dictionary (typically for iterating over and processing
    calibrations of a given type).

    This could also be done by the user in 1-2 lines of code, but involves
    some nested loops or (in)comprehensions that are less user-readable than
    a function call (to avoid duplication of type information in calibration
    dictionary, which would then have to be kept synchronized).

    Parameters
    ----------

    cal_dict : `dict`
        A dictionary of calibration files & associations, in the format
        produced by `init_cal_dict()` or `services.look_up_cals()`.

    cal_type : `str`
        Type of calibration to be looked up (matching a type name in the
        'associations' sub-dict).

    reference : `str` or `tuple` of `str`, optional
        One or more filenames with which matching calibrations must be
        associated, to limit the selection. By default, all available
        calibrations of type ``cal_type`` are selected.

    Returns
    -------

    `tuple` of (`str`, `list` of `str`)
        A tuple of (key, value) tuple pairs from the dictionary, where each
        first element is a calibration name label and each second element a
        list of constituent filenames. These are the ``items()`` of the input
        'calibrations' sub-dict that match ``cal_type`` (and, if specified,
        ``reference``). The result can be iterated over directly, for
        processing each entry in turn, or converted with ``dict()`` to
        reproduce the appropriate subset of ``cal_dict['calibrations']``.

    """

    # Convert reference to a tuple if a single name is provided:
    if isinstance(reference, basestring):
        reference = (reference,)

    # Extract the unique keys (calibration name labels) matching the specified
    # calibration type and reference files from the associations sub-dict:
    keys = {name for refkey, calmap in cal_dict[K_ASSOCIATIONS].iteritems() \
              if reference is None or refkey in reference \
            for typekey, name in calmap.iteritems() if typekey == cal_type}

    # Now extract the entries from the calibrations sub-dict corresponding to
    # the keys determined above. This could also be nested with the above set
    # comprehension if we want to enter the obfuscated 1-liner competition.
    # It would be preferable to leave this as a generator from a programmatic
    # viewpoint but a tuple is just more user-friendly.
    return tuple((key, cal_dict[K_CALIBRATIONS][key]) for key in keys)


def associate_cals(cals, inputs, cal_type, from_type=None, cal_dict=None):
    """
    Associate a given type of processed calibrations with the data they will
    subsequently be used to calibrate (both as `DataFile` instances).

    cals : `DataFileList` or `DataFile` or `dict`
        A list/dict of processed calibrations, available for association with
        zero or more matching ``inputs``. These are normally `DataFile`
        instances but, in the case of a dictionary, the values can also be
        objects of some other type (such as a list of background regions). If
        a `dict` is provided, with filename strings as the keys, calibrations
        are retrieved from it directly, using the names determined by
        calibration dictionary (or self) matching, otherwise such a dictionary
        is created internally from the supplied DataFile instances, keyed by
        their filenames. An entry must exist in ``cals`` for every match to
        ``inputs`` found in ``cal_dict``.

    inputs : `DataFileList` or `DataFile`
        DataFile instances with which processed calibrations are to be
        associated. If no calibration match can be determined from ``cal_dict``
        for a given input, the association is silently omitted (but it is an
        error for a calibration that is matched to be missing from ``cals``).

    cal_type : `str`
        Associated calibration type; this determines the name of the key
        (eg. 'flat' or 'trace') under which the calibrations get associated in
        the ``cals`` dictionaries of the ``inputs``. It must match a type name
        expected by later processing steps, rather than a type name defined in
        the calibration dictionary (though these are typically the same).

    from_type : `str`, optional
        Calibration dictionary type to associate, if different from the default
        of ``cal_type``. This can be used for multiple purposes, including:
        1. re-use of, for example, an associated 'flat' in the calibration
        dictionary (where available) as a 'trace' reference for spectral
        extraction; 2. when the type naming convention from a calibration
        look-up differs from that expected by subsequent processing steps
        (eg. due to different authors) and 3. associating calibration files
        with themselves (a case that may or may not exist explicitly in the
        calibration dictionary). The special parameter value of 'self' causes
        the cal_dict look-up to be skipped altogether and, for each input file,
        looks for a file in ``cals`` with the same original name. In case 1,
        where no substitute association exists, an entry ('trace') must be
        added manually to the calibration dictionary, pointing to either an
        existing calibration name label or a new one.

    cal_dict : `dict`, optional
        A populated calibration dictionary, in the specific format produced by
        `init_cal_dict()` or `services.look_up_cals()`, used to match ``cals``
        with the ``inputs``. Must be provided unless ``from_type== 'self'``.

    """

    # Ensure the inputs are DataFileLists if not already:
    inputs = DataFileList(data=inputs)

    # Are we associating a calibration with itself rather than some calibration
    # look-up result?
    from_self = from_type == 'self'

    # Convert cals to a dict keyed by name if not already in that format (to
    # avoid searching the list for matches repeatedly). Distinguish
    # DataFile(List) explicitly in case they have keys attributes added later.
    # When associating with 'self', convert keys to original filenames, as
    # otherwise found when looking up matches in the cal dict, to ignore
    # different processing stages. This is a bit fiddly to make more concise.
    if isinstance(cals, (DataFileList, DataFile)) or not hasattr(cals, 'keys'):
        if from_self:
            cals = {str(df.filename.orig) : df \
                    for df in DataFileList(data=cals)}
        else:
            cals = {str(df.filename) : df for df in DataFileList(data=cals)}

    # If passed a dict and associating with 'self', convert keys to original
    # filenames, as discussed above:
    elif from_self:
        cals = {FileName(key).orig : df for key, df in cals.iteritems()}

    # Default to using same cal_dict type convention as for associated type.
    # Could `from_type` be removed in favour of a brute-force look-up of files
    # associated with each input amongst the provided cals? Not if we want to
    # allow supplying a heterogenous mass of available cals, which might be
    # convenient from a user-bookeeping perspective.
    if from_type is None:
        from_type = cal_type

    # Calibration dict is only optional if self-associating:
    if cal_dict is None and not from_self:
        raise ValueError('cal_dict must be specified unless '\
                         'from_type is \'self\'')

    # To do: Add cal_dict validation, to produce a sensible error if needed.

    # Do the look-ups and perform the associations:
    for df in inputs:
        if from_self:
            match_name = df.filename.orig
        else:
            try:
                match_name = cal_dict[K_ASSOCIATIONS][df.filename.orig]\
                                     [from_type]
            except KeyError:
                match_name = None  # if no match is found, do nothing & move on

        if match_name is not None:
            try:
                df.cals[cal_type] = cals[match_name]
            except KeyError:
                raise KeyError('{0} missing from input cals'.format(match_name))

    # Return input files with cals now attached, for syntactic convenience:
    return inputs

