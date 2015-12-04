# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from ndmapper.libutils import splitext, addext

__all__ = ['init_cal_dict', 'add_cal_entry']


K_ASSOCIATIONS = 'associations'
K_CALIBRATIONS = 'calibrations'
K_CHECKSUMS    = 'checksums'


def init_cal_dict(filename=None):
    """
    Load a calibration dictionary from a JSON file if one is specified,
    otherwise initialize a new dictionary in the same format.

    Parameters
    ----------

    filename : str, optional
        Name of the JSON cache file to load, if any.

    Returns
    -------

    dict
        An (empty or populated) dictionary of calibration files & associations,
        in the required format. This consists of 3 sub-dictionaries,
        'associations', mapping each input filename to a dictionary of
        calibration name labels, keyed by calibration type, 'calibrations',
        mapping each calibration name label to a list of constituent raw files
        and 'checksums', mapping each raw calibration file to an optional
        checksum. [To do: add an example.]

    """

    # To do: this should log a warning if the file doesn't exist, in case the
    # user has mistyped the name (there's no way to distinguish that from the
    # cache just not being created yet, though unintentional differences won't
    # ordinarily occur as both cases will involve the same line of user code).

    # To do: cache results as JSON.

    # Read any existing file:
    # To do: remove the "and False" once persistence is implemented.
    if filename and False:
        cal_dict = None

    # Initialize a new calibration dictionary:
    else:
        cal_dict = {
            K_ASSOCIATIONS : {},
            K_CALIBRATIONS : {},
            K_CHECKSUMS : {}
        }

    return cal_dict


def recurse_file_cals(filename, obs_type, dependencies, lookup_fn, cal_dict):
    """
    Look up the calibrations on which a given file depends (including their
    calibrations in turn). This is a support function for the public user
    interface services.look_up_cals().

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
    services.look_up_cals()) for calibration entries that don't already exist
    in the dictionary, to preserve existing associations.

    Parameters
    ----------

    filename : str
        Name of the file for which the applicable calibration files are to
        be recorded.

    cal_type : str
        Type of calibration looked up (matching a name in the dependencies
        dictionary applicable to the instrument mode).

    matches : list of (str, str or None)
        List of (filename, checksum) pairs from a prior calibration look-up,
        to be included in the calibration dictionary.

    cal_dict : dict
        The calibration dictionary to populate, already initialized in the
        required format by init_cal_dict().

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

