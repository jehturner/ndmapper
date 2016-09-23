# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os, os.path

from pyraf import iraf
from ndmapper import config, ndprocess_defaults
from ndmapper.data import DataFileList
from ndmapper.utils import to_datafilelist
from ndmapper.iraf_task import run_task, get_extname_labels
from ndmapper.iraf_db import add_db_entry
from ndmapper.lib.gemini import gemini_iraf_helper

from ..gmos import *
from ..gmos import __all__

__all__ = __all__ + ['CAL_DEPS', 'biases', 'traces', 'arcs', 'flats', 'bg_reg',
                     'standards', 'calibrate_flux', 'apply_flux_cal',
                     'normalize_QE', 'shift_spectra']


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
biases, traces, arcs, flats, bg_reg, standards = {}, {}, {}, {}, {}, {}


@ndprocess_defaults
def calibrate_flux(inputs, out_names=None, reference=None, lookup_dir=None,
                   reprocess=None, interact=None):
    """
    Generate an instrumental sensitivity spectrum in magnitudes from the 1D
    integrated spectrum of a spectrophotometric standard star.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing the 1D spectrum of a spectrophotometric
        standard star. Usually there will just be one input and output, but
        multiple files are accepted.

    out_names : `str`-like or list of `str`-like, optional
        Names of output sensitivity spectra. If None (default), the names of
        the DataFile instances returned will be constructed from those of the
        input files, with '_sens' appended.

    reference : str or None, optional
        The name of a text file containing tabulated fluxes. This is usually
        the name of the standard star, in lower case. If None, a table matching
        the OBJECT header keyword will be sought (with any 'LTT' prefix
        abbreviated to 'L', to match the "noao$onedstds" convention). The
        precise bandpasses to be  used can be adjusted by editing a copy of
        this file.

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
    if not out_names:
        out_names = ['@input'] * len(inputs)
    elif isinstance(out_names, basestring) or isinstance(out_names, DataFile):
        out_names = [out_names]
    if len(out_names) != len(inputs):
        raise ValueError('inputs & out_names have unmatched lengths')

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
    for indf, outname in zip(inputs, out_names):

        # Default to obtaining the look-up table filename from the target
        # name in the header:
        status = indf.unloaded
        if not reference and 'OBJECT' in indf.meta:
            ref = str(indf.meta['OBJECT']).lower()
            # Map normal/Gemini naming to NOAO onedstds convention:
            if ref.startswith('ltt'):
                ref = 'l' + ref[3:]
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
        # determining 'out_names' in pure Python steps.
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
def apply_flux_cal(inputs, out_names=None, reprocess=None, interact=None):
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

    out_names : `str`-like or list of `str`-like, optional
        Names of flux-calibrated output spectra. If None (default), the names
        of the DataFile instances returned will be constructed from those of
        the input files, prefixed with 'c' as in the Gemini IRAF package.


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
    if not out_names:
        out_names = '@input'

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
        outputs={'output' : out_names}, prefix=prefix, suffix=None,
        comb_in=False, MEF_ext=False, path_param=None, reprocess=reprocess,
        sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], key_airmass=gemvars['key_airmass'],
        key_exptime=gemvars['key_exptime'], fl_vardq=gemvars['vardq'],
        fl_ext=False, fl_flux=True, fl_scale=True, fluxscale=1.0e15,
        ignoreaps=True, fl_fnu=False, extinction='',
        observatory=gemvars['observatory'], verbose=gemvars['verbose']
    )

    return result['output']


@ndprocess_defaults
def normalize_QE(inputs, out_names=None, reprocess=None, interact=None):
    """
    Correct for relative differences in quantum efficiency as a function of
    wavelength between the constituent GMOS CCDs (to avoid discontinuities
    between the detectors etc.).

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Bias-subtracted, unmosaicked input images, containing spectra whose
        fluxes are to be  normalized to a smooth function of QE(lambda) across
        the CCDs (equivalent to the response of the central CCD). Each input
        DataFile must already have an associated entry named 'arc' in its
        `cals` dictionary, identifying a corresponding arc exposure that has
        had its individual fibre spectra extracted and calibrated in wavelength
        (this information can be used in conjunction with meta-data to evaluate
        a pre-determined QE(wavelength) function on the original pixel grid).

    out_names : `str`-like or list of `str`-like, optional
        Names of corrected output images. If None (default), the names of the
        DataFile instances returned will be constructed from those of the input
        files, prefixed with 'q', as in the Gemini IRAF package.


    Returns
    -------

    output : DataFileList
        The corrected images produced by gqecorr.


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
    prefix = 'q'
    if not out_names:
        out_names = '@inimages'

    # Get a list of arcs from each input's "cals" dictionary.
    try:
        arcs = DataFileList(data=[df.cals['arc'] for df in inputs])
    except KeyError:
        raise KeyError('one or more inputs is missing an associated arc')

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # We don't need to define corrimages unless fl_keep=True, in which case it
    # would need making into a run_task output so reprocess handles it properly.
    # The corrections should be derivable by dividing the input & output anyway.
    result = run_task(
        'gemini.gmos.gqecorr',
        inputs={'inimages' : inputs, 'refimages' : arcs},
        outputs={'outimages' : out_names}, prefix=prefix,
        suffix=None, comb_in=False, MEF_ext=False, path_param=None,
        reprocess=reprocess, fl_correct=True, fl_keep=False, corrimages="",
        ifu_preced=1, fl_vardq=gemvars['vardq'], sci_ext=labels['data'],
        var_ext=labels['uncertainty'], dq_ext=labels['flags'],
        mdf_ext=gemvars['mask_table_name'], key_detsec=gemvars['key_detsec'],
        key_ccdsec=gemvars['key_ccdsec'], key_datasec=gemvars['key_datasec'],
        key_biassec=gemvars['key_biassec'], key_ccdsum=gemvars['key_ccdsum'],
        qecorr_data=gemvars['gmosdata']+'gmosQEfactors.dat',
        database='database/', verbose=gemvars['verbose']
    )

    return result['outimages']


@ndprocess_defaults
def shift_spectra(inputs, shift=None):
    """
    Apply a zero-point shift to the wavelength calibrations of the input
    spectra, eg. to correct the solution from a day-time arc spectrum for
    flexure with respect to science exposures, as determined from sky lines.

    Currently, this is a bit of a placeholder step that just corrects zero
    points in an IRAF database and must be run prior to ``rectify_wavelength``.
    This is useful for making a first-order correction to the arc, eg. so that
    QE correction will be performed accurately. It may additionally be
    necessary to correct for (a smaller amount of) mutual flexure between
    science exposures; at some point, this function will likely be generalized
    to allow adjusting the WCS as well, after rectification to linear
    wavelength.

    The database file(s) to which the shifts are added will be those
    corresponding to the 'REFSPEC1' header keyword, if present, otherwise each
    database file is expected to have the same name as the input file, prefixed
    with 'id' and with its filename extension replaced by a 3-digit integer
    suffix (eg. '_001'), corresponding to the slit number. Where a reference
    spectrum is defined, corrections can be applied to it via any spectra that
    reference it, but doing so will overwrite any shift applied previously for
    other spectra that reference the same arc.

    [To do: consider having this clone the solution and use the copy?]


    Parameters
    ----------

    inputs : DataFileList or DataFile
        Spectra on which ``calibrate_wavelength`` (or IRAF gswavelength) has
        been run, to generate a wavelength database file.

    shift : float or None, optional
        Shift to apply to each aperture, in the units of the database file 
        (normally Angstroms). If None, the 'shift' value from each input's
        ``cals`` dictionary is used, where present, otherwise the value
        defaults to 0.


    Returns
    -------

    DataFileList
        The inputs, with adjusted wavelength zero points (currently just in the
        IRAF database).

    """

    inputs = to_datafilelist(inputs)

    for df in inputs:

        if shift is None:
            shift = float(df.cals['shift']) if 'shift' in df.cals else 0.0

        for n, ndd in enumerate(df, start=1):

            unloaded = ndd.unloaded

            # If a reference spectrum is defined, edit its database, otherwise
            # look for one named after the input file:
            if 'REFSPEC1' in ndd.meta:
                slitdbname = ndd.meta['REFSPEC1']
            else:
                slitdbname = '{0}_{1:03}'.format(df.filename.root, n)

            dbfile = os.path.join('database', 'id' + slitdbname)
            add_db_entry(dbfile, 'shift', shift)

            # Until we have more intelligent checksumming, avoid leaving the
            # input files flagged as dirty just due to the above meta look-up:
            ndd._unloaded = unloaded

    return inputs

