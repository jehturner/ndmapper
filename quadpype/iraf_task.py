# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# Draft module to execute IRAF tasks with PyRAF.

import os.path
import tempfile
import datetime
import traceback
from pyraf import iraf
from . import config
from .data import DataFile, DataFileList


def run_task(taskname, inputs, outputs=None, prefix=None, combine=False, \
             MEF_ext=True, logfile=None, **params):
    """
    Wrapper to run an IRAF task on one or more DataFile objects and collect
    the results.

    Parameters
    ----------

    taskname : str
        Name of the IRAF task to run, optionally (if the package is already
        loaded) including the package hierarchy (eg. 'gemini.gmos.gfcube'). 

    inputs : dict of str : (DataFile or DataFileList)
        Dictionary mapping task parameter names to one or more input DataFile
        instances to pass one at a time to the task (combine=False) or all
        together (combine=True).

    outputs : (dict of str : (DataFile or DataFileList or str)) or None
        Specify output parameter name(s) and their filename value(s), if any.
        The same dictionary is returned as output after applying any automatic
        modifications.

        If the "prefix" parameter is set, the value(s) may name a parameter
        from the inputs dictionary, prefixed with '@' (eg. "@infiles"), to
        create the output names based on the corresponding input names, or
        prefixed with '!' to create a single output name based on the first
        input name.

    prefix : str or None
        A default prefix to add to the existing input name(s) if no output
        names are specified.

    combine : bool
        Pass all the inputs to the task in one call (eg. for stacking),
        instead of calling it on them one at a time by default (since not
        all tasks support input lists)?

    MEF_ext : bool
        Specify and iterate over FITS image extensions, for tasks expecting
        simple FITS as input (eg. core IRAF tasks; default = True)? This
        should be set to False for tasks that already handle multi-extension
        FITS files (eg. from Gemini IRAF) or when the input files are already
        simple FITS. The extension names iterated over are defined by the
        package "config" dictionary.

        The number of data_name extensions must be the same for every input file
        or one (in which case that single extension will be re-used for each
        iteration over the extensions of the other input files).

    logfile : str or {str : str} or None
        Optional filename for logging output, which includes any IRAF log
        contents (delimited by run_task status lines) or Python exceptions.

        The default of None causes the value of the package configuration
        variable "quadpype.config['logfile']" to be used, which itself
        defaults to None (in which case no log is written).

        Where only a filename string is provided and the IRAF task has a
        parameter named "logfile", the corresponding IRAF log will be captured
        automatically, otherwise only status information and Python exceptions
        will get recorded. Where single-item dictionary is provided, the key
        specifies an alternative IRAF "log file" parameter name to use and the
        value again specifies the output filename [not implemented]. A special
        key string of 'STDOUT' causes the standard output to be captured in
        the log (instead of any IRAF log file contents).

        The IRAF log contents relevant to each file are also appended to the
        corresponding output DataFile's log attribute (whether or not a log
        file is specified here and written to disk).

    params : dict
        Named IRAF task parameters. These may include ancillary input or
        output filenames that don't need to be tracked by the main inputs &
        outputs dictionaries.

    There is no support for mixing MEF- and simple FITS files in a single call.

    There is some possibility of the combine or extname parameters (possibly 
    also prefix if it has a different meaning from the Gemini one) conflicting
    with IRAF parameters, so these should perhaps be renamed more obscurely.


    Returns
    -------

    outputs : dict of str : (DataFile or DataFileList or str)
        The DataFile outputs specified by the input parameter "outputs".

    """

    # Print initial run_task() delimiter:
    dt = datetime.datetime.now()
    logstart = '-----\nSTART run_task(%s)  %s' % \
        (taskname, dt.strftime('%Y-%m-%d %H:%M:%S'))
    print logstart

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
        inplen = conv_io_pars(inputs)
        nfiles = max(inplen)

        # Apply any specified prefix to the filenames of the reference input
        # parameter to form the corresponding output filenames:
        if prefix is not None and outputs is not None:
            for key, val in outputs.iteritems():
                if isinstance(val, basestring) and val and val[0] in '!@':
                    refpar = val[1:]
                    if val[0] == '!':
                        namerange = slice(0, 1)        # use first filename
                    else:
                        namerange = slice(None, None)  # use all filenames
                    preflist = DataFileList()
                    if refpar in inputs:
                        for datafile in inputs[refpar][namerange]:
                            newfile = DataFile(filename=datafile.filename)
                            newfile.filename.dir=''  # output goes in CWD
                            newfile.filename.prefix = \
                                prefix + newfile.filename.prefix
                            preflist.append(newfile)
                    else:
                        raise ValueError('parameter name %s for prefix not '\
                                         'in inputs dictionary' % refpar)
                    outputs[key] = preflist

        # Make sure output parameters are filename lists, as for the input:
        outplen = conv_io_pars(outputs)

        # Now if we're iterating over the files and feeding them to the task
        # one at a time for each parameter, expand out any single filenames
        # implicitly to the length of any input list(s) so we can iterate over
        # everything together, complaining if given lists of different lengths
        # (ie. whose correspondence cannot be determined unambiguously).
        if not combine and nfiles > 1:
            for parset, parlens in [(inputs, inplen), (outputs, outplen)]:
                for param, n in zip(parset, parlens):
                    if n == 1:
                        parset[param] *= nfiles
                    elif n != nfiles:
                        raise ValueError('input/output file lists have ' \
                            'unmatched lengths and combine=False')

        # Create a list of inputs and outputs for each set of files on which
        # the task is run (just one if combine=True). At this point, the input
        # & output lists should all have the same lengths if combine=False.
        if combine:
            inlist = [inputs]
            outlist = [outputs]
        else:
            inlist = [{key : DataFileList(inputs[key][n]) for key in \
                       inputs.keys()} for n in range(nfiles)]
            outlist = [{key : DataFileList(outputs[key][n]) for key in \
                       outputs.keys()} for n in range(nfiles)]

        # Add run_task delimiter to individual DataFile log attributes (done
        # separately from IRAF log so it's still present if the latter isn't):
        for dfl in outputs.values():
            for df in dfl:
                df.log += '\n%s\n' % logstart

        # Define IRAF string format for any extension FITS extension numbers:
        in_extfmt = '[%s]'
        out_extfmt = '[%s,%s,append]' % (config['data_name'], '%s')

        # Iterate over the parameter set(s) and run the task on each one:
        for inpset, outpset in zip(inlist, outlist):

            # When MEF_ext=True, we require all the inputs to have the same or
            # unit length at each iteration (or overall if combine=True) and
            # the same EXTVERs in order to match data extensions unambiguously
            # between the inputs. While one could envisage more intelligent
            # expansion behaviour than this for special cases of combine=True
            # (such as requiring the lengths to match only between files at
            # the same list positions) it would be difficult to generalize
            # without unnecessarily imposing fixed relationships between sets
            # of input files (eg. there's nothing to stop an IRAF task from
            # taking a different number of files for each input or combining
            # them in some way other than one operation over all the inputs
            # per list position). The most likely case of iterating implicitly
            # over MEF extensions for multiple groups of files that only match
            # within each group can be handled using combine=False.
            if MEF_ext:
                # List EXTVERs for each input file (matching inpset dict):
                extdict = {param : [{ndd._io.group_id : ndd._io.data_idx \
                                     for ndd in df} for df in inpset[param]] \
                           for param in inpset}
                # Also derive a flat list of all sorted EXTVER lists, so we
                # can easily check that they match & get the nominal EXTVERs:
                allvers = [sorted(extmap.iterkeys()) for extmaps in \
                           extdict.itervalues() for extmap in extmaps]
                # Find the longest extension list, which we'll iterate over if
                # all goes well (the others should be the same or unit length):
                extvers = max(allvers, key=len)
                # Fail if the other non-unit-length EXTVER lists don't match:
                if not all([dfvers == extvers for dfvers in allvers
                            if len(dfvers) > 1]):
                    raise ValueError('non-matching input MEF EXTVERs')

            # Not iterating explicitly over MEF extensions:
            else:
                # Dummy dict to iterate over below, avoiding duplication:
                extdict = {param : [None for df in inpset[param]] \
                           for param in inpset}
                extvers = ['']

            # Run the task once on per MEF extension, if applicable, otherwise
            # just once in total:
            for ver in extvers:

                # Complete the IRAF parameter set with input/output file lists
                # for this iteration over the input files:
                for param in inpset:

                    # IRAF filenames for this parameter:
                    fnlist = []

                    # Iterate over files for this parameter & their ext maps:
                    for df, dfextmap in zip(inpset[param], extdict[param]):

                        # OS filename, without any ext:
                        fn = str(df)

                        # If iterating over FITS extensions, find the data
                        # extension number corresponding to this extver, or if
                        # there's only one data ext re-use that number:
                        if dfextmap:
                            if len(df) == 1:
                                fn += in_extfmt % df[0]._io.data_idx
                            else:
                                fn += in_extfmt % dfextmap[ver]

                        fnlist.append(fn)

                    # Convert filename list to IRAF comma-separated string
                    # and add the relevant task parameter/value:
                    params[param] = ','.join(fnlist)

                # Similar IRAF comma-separated list for the output files. Here
                # we just give IRAF the extname/ver to use instead of the ext.
                for param in outpset:
                    params[param] = ','.join( \
                        [str(df)+(out_extfmt % ver) if ver else str(df) \
                         for df in outpset[param]])

                # Specify log file for IRAF. Even if the user specifies a name,
                # use a temporary file and capture its contents before
                # appending to the user-specified file.
                if logpar is not None:
                    templog = tempfile.NamedTemporaryFile()
                    params[logpar] = templog.name
                else:
                    templog = None

                # print 'pars', params

                # Execute with our Python inputs converted to IRAF-style pars:
                try:
                    task(**params)

                # Note that PyRAF doesn't trap failures in IRAF tasks that
                # accept input file lists and only issue a warning and carry on
                # when an error specific to one of the files occurs, so we have
                # to check separately that the expected output gets created to
                # be confident it worked.

                except (iraf.IrafError, KeyError), errstr:
                    # Currently this is just a placeholder for any clean-up.
                    raise

                # Save any temporary IRAF log output whether or not the task
                # succeeded:
                finally:
                    if templog is not None:
                        logtext = templog.read()
                        # Copy the temporary IRAF log into user-specified log:
                        if userlog:
                            userlog.write(logtext)
                        # Attach log text to all output DataFile objects since,
                        # where there's more than one, we don't know which if
                        # any is the main one and it may apply to them all:
                        for dfl in outpset.values():
                            for df in dfl:
                                df.log += logtext
                        templog.close()

                # Check that any non-blank output filenames were created:
                for key, val in outpset.items():
                    for df in val:
                        namestr = str(df)
                        if namestr and not os.path.isfile(namestr):
                            raise RuntimeError('No file %s after running %s' % \
                                (namestr, taskname))

        # Print final run_task() delimiter:
        dt = datetime.datetime.now()
        logend='END   run_task(%s)  %s\n-----\n' % \
            (taskname, dt.strftime('%Y-%m-%d %H:%M:%S'))
        if userlog: userlog.write('%s' % logend)
        print logend

        # Add delimiter to individual DataFile log attributes as well:
        for dfl in outputs.values():
            for df in dfl:
                df.log += '\n%s\n' % logend
                # print 'dfl', df.log

    except:
        if userlog:
            userlog.write(traceback.format_exc())  # save traceback
        raise

    finally:
        if userlog:
            userlog.close()

    # Need to re-load any files listed in the outputs here, to figure out
    # their new lengths etc.

    # Return the outputs dictionary provided by the user, after expanding
    # any input parameter refs expanded to DataFileLists etc.:
    return outputs

    # TO DO:
    # - Check that output file doesn't already exist.
    # - Improve error trapping somehow.
    # - Finish logging behaviour as per the docstring(?).
    # - Scan for duplicates in output filenames? May be legitimate??
    # - Want to ensure the name ends in .fits etc?
    # - Copy files with path prefixes into the CWD under temporary filenames
    #   (to avoid conflicts) and log their correspondence to temporary names
    #   so the filenames aren't obfuscated too much in the log.
    #   - Don't cd to the data because it might break user & IRAF expectations
    #     about login.cl, database directories etc.
    # - Ensure check for existence of output files omits the [SCI].
    # - Capture any stdout or stderr as well as the log?
    # - Consider allowing params other than input & output to be DataFileLists
    #   and converting them to strings as needed, for convenience.
        
def conv_io_pars(pardict):
    """
    Convert dict of input or output Python file lists/names to dict of
    type DataFileList and return a list of the list lengths (private).
    """
    for param in pardict:
        # First cast any strings to a DataFile and then any DataFiles
        # to a DataFileList:
        parval = pardict[param]
        if isinstance(parval, basestring):
            parval = DataFile(filename=parval)
        try:
            pardict[param] = DataFileList(parval)
        except TypeError:
            raise TypeError('could not convert %s to DataFileList' \
                            % type(parval))
    parlen = [len(val) for val in pardict.values()]
    return parlen

