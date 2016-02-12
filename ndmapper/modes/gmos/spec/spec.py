# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os, os.path

from pyraf import iraf
from ndmapper import config, ndprocess_defaults
from ndmapper.data import DataFileList
from ndmapper.utils import to_datafilelist
from ndmapper.iraf_task import run_task, get_extname_labels
from ndmapper.modes.gemini import gemini_iraf_helper

from ..gmos import *
from ..gmos import __all__

__all__ = __all__ + ['CAL_DEPS', 'biases', 'traces', 'arcs', 'flats',
                     'standards', 'calibrate_flux', 'apply_flux_cal']


# Default calibration dependence for GMOS spectroscopy. The 'spectwilight'
# type is left out for now, as it's currently not used here.
CAL_DEPS = {'target' : ['specphot', 'flat', 'arc', 'bias'],
            'specphot' : ['flat', 'arc', 'bias'],
            'flat' : ['arc', 'bias'],
            'arc' : ['bias'],
            'bias' : []
           }

# For convenience & avoidance of visual noise, initialize these dictionaries
# that can be used for managing calibrations within user scripts:
biases, traces, arcs, flats, standards = {}, {}, {}, {}, {}


@ndprocess_defaults
def calibrate_flux(inputs, outputs=None, reference=None, lookup_dir=None,
                   reprocess=None, interact=None):
    """
    Generate an instrumental sensitivity spectrum from the 1D integrated
    spectrum of a spectrophotometric standard star.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing the 1D spectrum of a spectrophotometric
        standard star. Usually there will just be one input and output, but
        multiple files are accepted.

    outputs : DataFileList or DataFile, optional
        Output sensitivity spectra in magnitudes. If None (default), a new
        DataFileList will be returned, whose names are constructed from those
        of the input files, with '_sens' appended.

    reference : str or None, optional
        The name of a text file containing tabulated fluxes. This is usually
        the name of the standard star, in lower case. If None, a table matching
        the OBJECT header keyword will be sought. The precise bandpasses to be
        used can be adjusted by editing a copy of this file.

    lookup_dir : str or None, optional
        Directory name in which to find the `reference` file with tabulated
        fluxes, if not the current working directory. IRAF syntax may be used.

    See "help gsstandard" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The sensitivity spectra produced by gsstandard.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    interact : bool
        Enable interactive identification of bandpasses and fitting of the
        sensitivity curve (default False)? This may be overridden by the
        step's own "interact" parameter.

    """
    inputs = to_datafilelist(inputs)

    # Use default prefix if output filename unspecified & convert to a list:
    suffix = '_sens'
    if not outputs:
        outputs = ['@input'] * len(inputs)
    elif isinstance(outputs, basestring) or isinstance(outputs, DataFile):
        outputs = [outputs]
    if len(outputs) != len(inputs):
        raise ValueError('inputs & outputs have unmatched lengths')

    # Default to finding the look-up table in the CWD:
    if not lookup_dir:
        lookup_dir = './'

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # Loop over the files explicitly, rather than letting run_task do it,
    # since we need to modify or delete the "sfile" output each time:
    outlist = DataFileList(mode='update')
    for indf, outname in zip(inputs, outputs):

        # Default to obtaining the look-up table filename from the target
        # name in the header:
        status = indf.unloaded
        if not reference and 'OBJECT' in indf.meta:
            ref = str(indf.meta['OBJECT']).lower()
        else:
            ref = reference
        # Avoid run_task making an unnecessary temp copy pending a more
        # intelligent mechanism for discerning when data are modified:
        indf._unloaded = status

        # Generate a name for the 'sfile' output at each iteration. This is
        # not used directly but is kept for the user's information. Currently
        # there is no way to specify a different name.
        sfile = indf.filename.root + '_std'
        # This is insufficient to avoid the task failing if reprocess is False
        # and the main 'sfunction' output does not exist. It's probably best to
        # remove the file in that case too, given that this particular output
        # is only informational, but wait until we have more logic for
        # determining 'outputs' in pure Python steps.
        if reprocess and os.path.exists(sfile):
            os.remove(sfile)

        # Use a simple summation, at least for now, as getting the rejection
        # etc. right with the level of contrast involved can be tricky.
        result = run_task(
            'gemini.gmos.gsstandard',
            inputs={'input' : indf}, outputs={'sfunction' : outname},
            prefix=None, suffix=suffix, comb_in=False, MEF_ext=False,
            path_param=None, reprocess=reprocess, sfile=sfile,
            sci_ext=labels['data'], var_ext=labels['uncertainty'],
            dq_ext=labels['flags'], key_airmass=gemvars['key_airmass'],
            key_exptime=gemvars['key_exptime'], fl_inter=interact,
            starname=ref, samestar=True, apertures='', beamswitch=False,
            bandwidth=iraf.INDEF, bandsep=iraf.INDEF, fnuzero=3.68e-20,
            caldir=lookup_dir, observatory=gemvars['observatory'], mag='',
            magband='', teff='', ignoreaps=True, extinction='',
            out_extinction='extinct.dat', function='spline3', order=3,
            graphs='sr', marks='plus cross box', colors='2 1 3 4',
            verbose=gemvars['verbose']
        )

        outlist.append(result['sfunction'][0])

    return outlist


@ndprocess_defaults
def apply_flux_cal(inputs, outputs=None, reprocess=None, interact=None):
    """
    Apply a previously-measured sensitivity spectrum to each input to
    convert the data values to spectrophotometric flux units.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing spectra to be converted to physical flux
        units. Each input DataFile must already have an associated entry named
        'specphot' in its `cals` dictionary, corresponding to the measured
        instrumental sensitivity spectrum.

    outputs : DataFileList or DataFile, optional
        Output spectra with units of flux. If None (default), a new
        DataFileList will be returned, whose names are constructed from those
        of the input files, prefixed with 'c' as in the Gemini IRAF package.


    Returns
    -------

    output : DataFileList
        The flux-calibrated spectra produced by gscalibrate.


    Package 'config' options
    ------------------------

    use_uncert : bool
        Enable NDData 'uncertainty' (variance) propagation (default True)?

    use_flags : bool
        Enable NDData 'flags' (data quality) propagation (default True)?

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """

    # Use default prefix if output filename unspecified:
    prefix = 'c'
    if not outputs:
        outputs = '@input'

    # Get a list of sensitivity spectra to use from the input file "cals"
    # dictionaries.
    try:
        sensspec = DataFileList(data=[df.cals['specphot'] for df in inputs])
    except KeyError:
        raise KeyError('one or more inputs is missing an associated '\
                       'sensitivity spectrum')

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # To do: deal with extinction correction?
    result = run_task(
        'gemini.gmos.gscalibrate',
        inputs={'input' : inputs, 'sfunction' : sensspec},
        outputs={'output' : outputs}, prefix=prefix, suffix=None,
        comb_in=False, MEF_ext=False, path_param=None, reprocess=reprocess,
        sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], key_airmass=gemvars['key_airmass'],
        key_exptime=gemvars['key_exptime'], fl_vardq=gemvars['vardq'],
        fl_ext=False, fl_flux=True, fl_scale=True, fluxscale=1.0e15,
        ignoreaps=True, fl_fnu=False, extinction='',
        observatory=gemvars['observatory'], verbose=gemvars['verbose']
    )

    return result['output']

