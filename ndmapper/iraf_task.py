# Copyright(c) 2015-2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# Draft module to execute IRAF tasks with PyRAF.

# This works quite well but could do with refactoring for modularity after
# growing a bit ... erm ... organically as the the concept has been sketched
# out and features added. Also, an interface needs to be added for capturing
# the standard output and the logging needs tidying up significantly.

"""
A module to help execute IRAF tasks conveniently as part of an NDMapper data
reduction sequence, with DataFileList objects as inputs and outputs.
"""

import os
import os.path
import tempfile
import datetime
import traceback

from pyraf import iraf

from . import config
from .data import FileName, DataFile, DataFileList, temp_saved_datafile
from .utils import to_datafilelist


__all__ = ['run_task', 'get_extname_labels']


def run_task(taskname, inputs, outputs=None, prefix=None, suffix=None,
             comb_in=False, MEF_ext=True, path_param=None, reprocess=None,
             logfile=None, **params):
    """
    Wrapper to run an IRAF task on one or more `DataFile` objects and collect
    the results.

    Parameters
    ----------

    taskname : `str`
        Name of the IRAF task to run, optionally (if the package is already
        loaded) including the package hierarchy (eg. ``gemini.gmos.gfcube``). 

    inputs : `dict`
        Dictionary mapping task parameter names to one or more input `DataFile`
        instances to pass one at a time to the task (``comb_in is False``) or
        all together (``comb_in is True``). All the named files must already
        exist.

    outputs : `dict` or `None`, optional
        Specify output parameter name(s) and their filename value(s), if any,
        as dictionary keys & values. The files named must not already exist. 
        The equivalent dictionary is returned as output, after applying any
        automatic modifications. The `dict` values may have any type that can
        be converted to a string (eg. `FileName`), either individually or in
        a sequence.

        If the ``prefix`` and/or ``suffix`` parameter is set, the value(s) may
        name a parameter from the inputs dictionary, prefixed with '@'
        (eg. ``@infiles``), to create the output names based on the
        corresponding input names, or prefixed with '!' to create a single
        output name based on the first input name.

    prefix : `str` or `None`
        A default prefix to add to existing input filename(s) to form the
        output filename(s), if the output parameter value(s) specify this
        behaviour.

    suffix : `str` or `None`
        A suffix to add to existing input filename(s) to form the output
        filename(s), if the output parameter value(s) specify this behaviour.

    comb_in : `bool`
        Pass all the inputs to the task at once, in a single call (eg. for
        stacking), instead of the default behaviour of calling the task on
        each file in turn (per input parameter)? This parameter is named
        obscurely to avoid conflicts with IRAF tasks
        (eg. ``imcombine.combine``).

    MEF_ext : `bool`
        Specify and iterate over FITS image extensions, for tasks expecting
        simple FITS as input (eg. core IRAF tasks; default `True`)? This should
        be set to `False` for tasks that already handle multi-extension FITS
        files (eg. from Gemini IRAF) or when the input files are already simple
        FITS. The extension names to be iterated over are defined when the
        input DataFile instances are created, defaulting to values kept in the
        package ``config`` dictionary.

        The number of extensions named ``config['labels']['data']`` (eg. 'SCI')
        must be the same for every input file or one (in which case that
        single extension will be re-used for each iteration over the extensions
        of the other input files).

    path_param : `str` or `None`
        Name of a task parameter (eg. ``rawpath``) used to specify the location
        of the input files, instead of including the full path in each input
        filename as usual. The DataFile paths are automatically stripped from
        the inputs and supplied to the task via this parameter instead. Output
        files are still assumed to reside in the current working directory
        unless otherwise specified. To use this option, all inputs containing
        a directory path (other than '') must reside in the same directory --
        if this is not the case and the IRAF task does not understand paths in
        filenames then the user will need to copy the input files to the
        current directory before running it. The user must not supply filenames
        in another directory to input parameters for which the IRAF task does
        not apply the path_param prefix (usually input calibrations), or the
        task will fail to find some or all of the inputs.

    reprocess : `bool` or `None`, optional
        Overwrite or re-use existing output files? The default of `None`
        redirects to the value of the package configuration variable
        ``ndmapper.config['reprocess']``, which in turn defaults to `None`,
        causing an error to be raised where output files already exist before
        processing. If the value is set to `True`, any existing files will be
        deleted before the IRAF task is run, while `False` is a no-op for
        existing files, causing them to be re-used as output without repeating
        the processing. If the IRAF task produces any intermediate files that
        are not included in ``outputs`` (ie. that are unknown to run_task), it
        is the caller's responsibility to delete them before repeating any
        processing. The task is always (re-)run where there are no `outputs`.

    logfile : `str` or `dict` or `None`, optional
        Optional filename for logging output, which includes any IRAF log
        contents (delimited by run_task status lines) or Python exceptions.

        The default of `None` causes the value of the package configuration
        variable ``ndmapper.config['logfile']`` to be used, which itself
        defaults to `None` (in which case no log is written).

        Where only a filename string is provided and the IRAF task has a
        parameter named ``logfile``, the corresponding IRAF log will be captured
        automatically, otherwise only status information and Python exceptions
        will get recorded. Where a single-item dictionary is provided, the key
        specifies an alternative IRAF "log file" parameter name to use and the
        value again specifies the output filename [not implemented]. A special
        key string of ``STDOUT`` will cause the standard output to be captured
        in the log (instead of any IRAF log file contents) [unimplemented].

        The IRAF log contents relevant to each file are also appended to the
        corresponding output DataFile's log attribute (whether or not a log
        file is specified here and written to disk).

    params : `dict`
        Named IRAF task parameters. These may include ancillary input or
        output filenames that don't need to be tracked by the main inputs &
        outputs dictionaries.

    Returns
    -------

    outputs : `dict` of `str` : `DataFileList`
        The DataFile objects named by the parameter ``outputs``, containing
        the results from IRAF.

    Notes
    -----

    There is no support for mixing MEF- and simple FITS files in a single call.

    In principle, ``prefix`` & ``suffix`` could conflict with any like-named
    IRAF parameters that have a different meaning from the Gemini convention
    (ie. a string that is added to the start/end of the input filename to
    provide an output name), but there appears to be only one such case in
    Ureka (sqiid.getcoo); likewise for ``MEF_ext``, which has no known uses
    elsewhere, ``path_param`` and ``reprocess``. It is assumed that the
    widely-used ``logfile`` will only ever have the usual meaning.

    """

    # Print initial run_task() delimiter:
    dt = datetime.datetime.now()
    logstart = '-----\nSTART run_task(%s)  %s' % \
        (taskname, dt.strftime('%Y-%m-%d %H:%M:%S'))
    print(logstart)

    # Start the log file:
    if logfile is None:
        logfile = config['logfile']
    if logfile is None:
        userlog = None
    elif isinstance(logfile, basestring):
        userlog = open(logfile, 'a')
        userlog.write('%s\n' % logstart)
    else:
        # Dict, to be implemented:
        raise NotImplementedError('logfile must currently be str or None')

    # Determine from the config dict whether to reprocess data, if unspecified:
    if reprocess is None:
        reprocess = config['reprocess']

    # Keep a list of any temporary files that need cleaning up when finished:
    tmplist = []

    # This giant try-except block just exists to log any tracebacks before
    # re-raising the exception:
    try:

        # Ensure host package(s) is/are loaded:
        pkglist = taskname.split('.')[:-1]
        for pkg in pkglist:
            eval('iraf.'+pkg+'(_doprint=0, Stdout=1)') # capture+discard stdout

        # Get reference to the task as a function:
        task = eval('iraf.'+taskname)

        # Restore any old task parameter settings that aren't overridden in
        # "params" to their defaults, to ensure the results are reproducible:
        task.unlearn()

        # Check for a recognizeable "log file" task parameter that can be used
        # to capture IRAF logging output:
        logpar = None
        for parname in ['logfile']:
            if parname in task.getParDict().keys():
                logpar = parname
                break

        # Ensure the main inputs & outputs are dictionaries, otherwise we don't
        # know what IRAF parameters they go with. I think we want this to be
        # fairly brittle rather than duck-typed, to avoid unexpected behaviour.
        if outputs is None: outputs = {}
        if not isinstance(inputs, dict) or not isinstance(outputs, dict):
            raise TypeError('values of inputs/outputs must be parameter=' \
                            ' value dictionaries')

        # Make sure each input parameter is expressed as a filename list and
        # determine how many sets of input files there are to iterate over
        # (should be as many as the length of the longest list).
        inplen = conv_io_pars(inputs, mode=None)  # defaults to mode='read'
        nfiles = max(inplen)

        # Input files are no longer required already to exist on disk here, as
        # the unloaded flag will be False otherwise, which now causes temporary
        # copies to get saved below, at temp_saved_datafile().

        # Set the task's path_param if specified (but not overridden by the
        # user), after ensuring the path to the files is unique:
        if path_param and path_param not in params:
            paths=set()
            for dfl in inputs.itervalues():
                for df in dfl:
                    if df.filename.dir:  # don't count CWD ('')
                        paths.add(df.filename.dir)
            ndirs = len(paths)
            if ndirs == 0:
                path = ''
            elif ndirs == 1:
                (path,) = paths
            else:
                raise ValueError('inputs must all have the same path when ' \
                    '\'path_param\' is set')
            path = os.path.join(path, '')  # incl. trailing slash unless blank
            params[path_param] = path

        # Apply any specified prefix to the filenames of the reference input
        # parameter to form the corresponding output filenames (this usage is
        # a bit dodgy after re-doing DataFile modes but still works):
        if outputs is not None:
            for key, val in outputs.iteritems():
                if isinstance(val, basestring) and val and val[0] in '!@':
                    if prefix is None and suffix is None:
                        raise ValueError('output \'%s\' requires missing '
                            'suffix/prefix value' % key)
                    refpar = val[1:]
                    if val[0] == '!':
                        namerange = slice(0, 1)        # use first filename
                    else:
                        namerange = slice(None, None)  # use all filenames
                    preflist = DataFileList(mode='overwrite')
                    if refpar in inputs:
                        for datafile in inputs[refpar][namerange]:
                            newfile = DataFile(filename=datafile.filename,
                                               mode='overwrite')
                            newfile.filename.dir=''  # output goes in CWD
                            if prefix is not None:
                                newfile.filename.prefix = \
                                    prefix + newfile.filename.prefix
                            if suffix is not None:
                                newfile.filename.suffix.append(suffix)
                            preflist.append(newfile)
                    else:
                        raise ValueError('parameter name %s for prefix not '\
                                         'in inputs dictionary' % refpar)
                    outputs[key] = preflist

        # Make sure output parameters are DataFileLists, as for the input,
        # selecting the mode according to the reprocess parameter. Use
        # overwrite mode for reprocess==False until the files get reloaded, so
        # DataFile won't complain if they don't already exist.
        mode = 'new' if reprocess is None else 'overwrite'
        outplen = conv_io_pars(outputs, mode=mode)

        # Save temporary copies (to the current directory) of any input files
        # that could have changed in memory, having done what's needed with
        # the original input filenames above. Creating copies can slow things
        # down by a factor of ~2-3 (eg. from 2.3s to 5.7s when adding two 270M
        # FITS files with 3 SCI extensions each on an SSD) but avoids
        # unexpected results due to unsaved changes. Disk space usage could be
        # improved by copying only those files needed at each iteration but
        # that would complicate expansion of lists to the same length below etc.

        # When path_param is set and *any* of the inputs needs saving, copies
        # must be made of all the files, since they must reside in the same
        # directory (and the original location may not be writeable).
        if path_param and not all([df.unloaded for dfl in inputs.itervalues() \
                                   for df in dfl]):
            copyall = True
            params[path_param] = ''
        else:
            copyall = False

        # Substitute original DataFiles for temporary copies only where needed
        # (if we're not sure there's an up-to-date copy on disk already). The
        # method for deciding this is currently a bit primitive (optimizing it
        # is a bit of a minefield) but avoids routine copies in the common case
        # where DataFileList is simply used as a list of files for IRAF.
        for dfl in inputs.itervalues():
            for n, df in enumerate(dfl):
                if copyall or not df.unloaded:
                    tdf = temp_saved_datafile(df)
                    tmplist.append(tdf)
                    dfl[n] = tdf

        # Consider adding a section here that enables comb_in automatically
        # if the number of output files (for at least one output parameter?)
        # is 1 (in which case the following check won't run). The parameter
        # can't be eliminated entirely because we can't distinguish looping
        # over N inputs with N outputs from running a task once that processes
        # N inputs together and produces N outputs (eg. WCS updates onto a
        # common system, scaling to a common sky level etc.). Could add
        # comb_in="auto"/None option or make it into "force_comb" or
        # "separate" etc.
        # - At this point, the prefix parameter has already been applied to
        #   generate one or more output filenames, if applicable.
        # - Document decision in the log.

        # Now if we're iterating over the files and feeding them to the task
        # one at a time for each parameter, expand out any single filenames
        # implicitly to the length of any input list(s) so we can iterate over
        # everything together, complaining if given lists of different lengths
        # (ie. whose correspondence cannot be determined unambiguously).
        if not comb_in and nfiles > 1:
            for parset, parlens in [(inputs, inplen), (outputs, outplen)]:
                for param, n in zip(parset, parlens):
                    if n == 1:
                        parset[param] *= nfiles
                    elif n != nfiles:
                        raise ValueError('input/output file lists have ' \
                            'unmatched lengths and comb_in=False')

        # Create a list of inputs & outputs for each set of files on which the
        # task is run (just one if comb_in=True). At this point, the input
        # & output lists should all have the same lengths if comb_in=False.
        if comb_in:
            inlist = [inputs]
            outlist = [outputs]
        else:
            inlist = [{key : DataFileList(data=inputs[key][n]) for key in \
                       inputs.keys()} for n in range(nfiles)]
            outlist = [{key : DataFileList(data=outputs[key][n], mode=mode) \
                        for key in outputs.keys()} for n in range(nfiles)]

        # Define IRAF string format for any extension FITS extension numbers:
        in_extfmt = '[%s]'
        out_extfmt = '[%s,%s,append]'

        # To avoid obscure failures, do a pre-iteration over output parameter
        # set(s) and ensure there are no duplicate outputs between iterations
        # or, when reprocess is False, incomplete subsets of existing outputs
        # from any given iteration. Duplicate outputs are allowed within an
        # iteration (ie. task call) and if not legitimate should eventually be
        # caught when the task itself complains. Other errors like file
        # permissions are dealt with later, by checking that the files actually
        # get created by IRAF. While we're at it, add the run_task delimiter to
        # the DataFile log attributes (done separately from the IRAF log so
        # it's still present when not using the latter).
        prevnames = []
        for outpset in outlist:
            iternames = []
            existing = non_existing = False
            for dfl in outpset.itervalues():
                for df in dfl:
                    df.log += '\n%s\n' % logstart
                    # (This comparison should also deal automatically with any
                    # explicit .fits extensions once DataFile handles them:)
                    name = os.path.abspath(str(df))
                    if name in prevnames:
                        raise IOError('duplicate output file: %s' % str(df))
                    if reprocess is False:
                        if os.path.exists(name):
                            existing = True
                        else:
                            non_existing = True
                    iternames.append(name)
            prevnames.extend(iternames)
            if reprocess is False and existing and non_existing:
                raise IOError('reprocess is False & a subset of these outputs '\
                              'already exist:\n  {0}'
                              .format('\n  '.join([str(df) for dfl \
                                                   in outpset.itervalues() \
                                                   for df in dfl])))

        # Iterate over the parameter set(s) and run the task on each one:
        for inpset, outpset in zip(inlist, outlist):

            call_task = True

            # If re-processing, delete any existing files from this task call
            # and if not, check whether the call can be skipped. This check is
            # done here, separately from the above section, to avoid premature
            # removal of results in case of failure:
            if reprocess is not None:
                # This variable also gets re-used in the call_task else clause:
                names = [str(df) for dfl in outpset.itervalues() for df in dfl]
                for name in names:
                    if os.path.exists(name):
                        if reprocess:
                            os.remove(name)
                        else:
                            call_task = False
                            break  # either all or none exist after above sec.

            # Execute the task unless its outputs already exist and reprocess
            # is False, in which case we just re-use the files instead (via the
            # same reload() call at the end):
            if call_task:

                # When MEF_ext=True, we require all the inputs to have the same
                # or unit length at each iteration (or overall if comb_in=True)
                # and the same EXTVERs in order to match data extensions
                # unambiguously between the inputs. While one could envisage
                # more intelligent expansion behaviour than this for special
                # cases of comb_in=True (such as requiring the lengths to match
                # only between files at the same list positions) it would be
                # difficult to generalize without unnecessarily imposing fixed
                # relationships between sets of input files (eg. there's
                # nothing to stop an IRAF task from taking a different number
                # of files for each input or combining them in some way other
                # than one operation over all the inputs per list position).
                # The most likely case of iterating implicitly over MEF
                # extensions for multiple groups of files that only match
                # within each group can be handled using comb_in=False.
                if MEF_ext:
                    # List EXTVERs for each input file (matching inpset dict):
                    extdict = {param : [{ndd.ident : ndd._io.data_idx \
                                         for ndd in df} \
                                        for df in inpset[param]] \
                               for param in inpset}
                    # Also derive a flat list of all sorted EXTVER lists, to
                    # check easily that they match & get the nominal EXTVERs:
                    allvers = [sorted(extmap.iterkeys()) for extmaps in \
                               extdict.itervalues() for extmap in extmaps]
                    # Find longest extension list, which we'll iterate over if
                    # all goes well (others should be the same or unit length):
                    extvers = max(allvers, key=len)
                    # Fail if other non-unit-length EXTVER lists don't match:
                    if not all([dfvers == extvers for dfvers in allvers
                                if len(dfvers) > 1]):
                        raise ValueError('non-matching input MEF EXTVERs')

                # Not iterating explicitly over MEF extensions:
                else:
                    # Dummy dict to iterate over below, avoiding duplication:
                    extdict = {param : [None for df in inpset[param]] \
                               for param in inpset}
                    extvers = ['']

                # Run the task once on per MEF extension, if applicable,
                # otherwise just once in total:
                for ver in extvers:

                    # Complete the IRAF parameter set with input/output file
                    # lists for this iteration over the input files:
                    for param in inpset:

                        # IRAF filenames for this parameter:
                        fnlist = []

                        # Iterate over files for this param & their ext maps:
                        for df, dfextmap in zip(inpset[param], extdict[param]):

                            # OS filename, without any MEF ext (& without any
                            # path if task expects a separate path parameter):
                            if path_param:
                                fn = str(FileName(df.filename, dirname=''))
                            else:
                                fn = str(df)

                            # If iterating over FITS extensions, find the data
                            # extension number corresponding to this extver, or
                            # if there's only one data ext re-use that number:
                            if dfextmap:
                                if len(df) == 1:
                                    fn += in_extfmt % df[0]._io.data_idx
                                else:
                                    fn += in_extfmt % dfextmap[ver]

                            fnlist.append(fn)

                        # Convert filename list to IRAF comma-separated string
                        # and add the relevant task parameter/value:
                        params[param] = ','.join(fnlist)

                    # Similar IRAF comma-separated list for output files. Here
                    # we just give IRAF the extname/ver instead of the ext.
                    for param in outpset:
                        params[param] = ','.join( \
                            [str(df)+(out_extfmt % (df._labels['data'], ver)) \
                             if ver else str(df) for df in outpset[param]])

                    # Specify log file for IRAF. Even if the user specifies a
                    # name, use a temporary file and capture its contents
                    # before appending to the user-specified file.
                    if logpar is not None:
                        templog = tempfile.NamedTemporaryFile()
                        params[logpar] = templog.name
                    else:
                        templog = None

                    # print('pars', params)

                    # Execute with Python inputs converted to IRAF-style pars:
                    try:
                        task(**params)

                    # Note that PyRAF doesn't trap failures in IRAF tasks that
                    # accept input file lists and only issue a warning and
                    # carry on when an error specific to one of the files
                    # occurs, so we have to check separately that the expected
                    # output gets created to be confident it worked.

                    except (iraf.IrafError, KeyError), errstr:
                        # Currently just a placeholder for any clean-up.
                        raise

                    # Save any temporary IRAF log output whether or not the
                    # task succeeded:
                    finally:
                        if templog is not None:
                            logtext = templog.read()
                            # Copy temporary IRAF log into user-specified log:
                            if userlog:
                                userlog.write(logtext)
                            # Attach log text to all output DataFile objects
                            # since, where there's more than one, we don't know
                            # which if any is the main one and it may apply to
                            # them all:
                            # To do: currently get overwritten by reload below?
                            for dfl in outpset.itervalues():
                                for df in dfl:
                                    df.log += logtext
                            templog.close()

                    # Check that any non-blank output filenames got created:
                    for key, val in outpset.items():
                        for df in val:
                            namestr = str(df)
                            if namestr and not os.path.isfile(namestr):
                                raise RuntimeError(
                                    'No file %s after running %s' % \
                                    (namestr, taskname)
                                )

                # Here we would clean up any temp copies of input files from
                # this iteration over a given set of files, if and when the
                # copies are made per iteration instead of all at the start.

            # If processing was skipped, note that in the log:
            else:
                msg = 'Skip processing & re-use existing outputs for:'
                for name in names:   # from "if reprocess" at start of loop
                    msg += '\n  {0}'.format(name)
                for dfl in outpset.itervalues():
                    for df in dfl:
                        df.log += msg
                if userlog:
                    userlog.write('{0}\n'.format(msg))
                print(msg)

        # Print final run_task() delimiter:
        dt = datetime.datetime.now()
        logend='END   run_task(%s)  %s\n-----\n' % \
            (taskname, dt.strftime('%Y-%m-%d %H:%M:%S'))
        if userlog: userlog.write('%s' % logend)
        print(logend)

        # Add delimiter to individual DataFile log attributes as well:
        for dfl in outputs.itervalues():
            for df in dfl:
                df.log += '\n%s\n' % logend
                # print('dfl', df.log)

    except:
        if userlog:
            userlog.write(traceback.format_exc())  # save traceback
        raise

    finally:
        if userlog:
            userlog.close()
        for df in tmplist:
            os.remove(str(df))

    # Map data from files listed in the outputs after their creation by IRAF:
    for param in outputs:
        for df in outputs[param]:
            df.reload()

    # Return the outputs dictionary provided by the user, after expanding
    # any input parameter refs expanded to DataFileLists etc.:
    return outputs

    # TO DO:
    # - Finish logging behaviour as per the docstring(?).
    #   - Also send record of inputs/outputs to the log only in case the
    #     IRAF task doesn't do it?
    # - Want to ensure the name ends in .fits etc?
    #   - Would do this in DataFile, otherwise it would judge incorrectly
    #     whether the file already exists (and in any case DataFile needs to
    #     know how to load the file etc.).
    #     - Use a (prioritized?) list of all recognized file types rather than
    #       making DataFile FITS-specific.
    #       - Only when the file mode corresponds to input files(?). Otherwise
    #         there's no way to know which is the right implicit extension.
    #         Unless there's a package-configurable default?
    # - Capture any stdout or stderr as well as the log?
    # - Consider allowing params other than input & output to be DataFileLists
    #   and converting them to strings as needed, for convenience.
    # - Add IRAF tests with another 2 FITS files: with VAR/DQ, unnumbered.
    # - Improve error trapping further??
    #   - Check gemini status parameters?


def conv_io_pars(pardict, mode):
    """
    Convert `dict` of input or output Python file lists/names to `dict` of
    type `DataFileList` and return a `list` of the list lengths (private).

    """
    for param in pardict:
        pardict[param] = to_datafilelist(pardict[param], mode=mode)

    parlen = [len(val) for val in pardict.itervalues()]
    return parlen


def get_extname_labels(datafiles):
    """
    Ensure that all `DataFile` instances in datafiles use the same convention
    for labelling `NDData` component arrays (eg. 'SCI', 'VAR', 'DQ') and return
    the corresponding labels dictionary. The dictionary values can then be
    used to specify the EXTNAME conventions for IRAF tasks that operate on
    multi-extension FITS files (which is unnecessary in Python, where `NDData`
    defines which array is which). An empty dictionary is returned if the input
    list is empty.

    Raises `ValueError` if the constituent label dictionaries differ.
    """
    if isinstance(datafiles, DataFile):
        datafiles = [datafiles]
    elif not hasattr(datafiles, '__iter__') or \
         not all([isinstance(df, DataFile) for df in datafiles]):
        raise TypeError('datafiles must be a DataFile list or a DataFile')

    unique_items = set([tuple(df._labels.iteritems()) for df in datafiles])
    if len(unique_items) > 1:
        raise ValueError('datafiles must all have the same "labels" convention')

    labels = dict(unique_items.pop()) if unique_items else {}

    return labels

