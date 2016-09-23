# Copyright(c) 2015-2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os
from pyraf import iraf

from ndmapper import config, ndprocess_defaults
from ndmapper.data import FileName, DataFile, DataFileList
from ndmapper.iraf_task import run_task, get_extname_labels
from ndmapper.libutils import new_filename
from ndmapper.utils import to_datafilelist

from ...gemini import gemini_iraf_helper
from .spec import *
from .spec import __all__

__all__ = __all__ + ['prepare', 'subtract_bias', 'extract_spectra',
                     'calibrate_wavelength', 'rectify_wavelength', 'make_flat',
                     'subtract_sky', 'resample_to_cube', 'sum_spectra',
                     'background_regions', 'subtract_bg', 'align_wcs',
                     'mosaic', 'log_rebin']


@ndprocess_defaults
def prepare(inputs, out_names=None, mdf=None, reprocess=None):
    """
    Update the meta-data from raw GMOS images in preparation for subsequent
    processing.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images to be "prepared" by attaching a "mask definition file"
        (MDF) extension, describing the slit mapping to sky, and performing
        meta-data updates.

    out_names : `str`-like or list of `str`-like, optional
        Names of the output "prepared" files. If None (default), the names
        of the DataFile instances returned will be constructed from those of
        the input files, prefixed with 'g' as in the Gemini IRAF package.

    mdf : str, optional
        Name of a FITS mask definition file (in the current working directory)
        to be attached as a FITS extension of each output file. If None and all
        the inputs have mdf associations in their `cals` dictionaries, those
        mdfs are used, otherwise, if there are no such associations, the GMOS
        package default is used.

    See "help gfreduce" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The "prepared" images produced by gfreduce.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """

    # Use default prefix if output filename unspecified:
    prefix = 'g'
    if not out_names:
        out_names = '@inimages'

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    task_inputs = {'inimages' : inputs}
    optargs = {}

    # Get MDF from calibration associations if present, otherwise look in
    # gmos$data by default or the CWD if a file is specified via the parameter:
    if mdf is None:
        mdfs = [df.cals['mdf'] if 'mdf' in df.cals else None for df in inputs]
        missing = [entry is None for entry in mdfs]
        if all(missing):
            optargs['mdffile'] = 'default'
            optargs['mdfdir'] = gemvars['gmosdata']
        elif any(missing):
            raise KeyError('all or no inputs must have associated mdf tables')
        else:
            task_inputs['mdffile'] = DataFileList(data=mdfs)
    else:
        optargs['mdffile'] = mdf
        optargs['mdfdir'] = ''

    # Here some "unnecessary" parameters are defined just to be certain they
    # don't change anything and so we can easily copy this when wrapping
    # other operations with gfreduce later. Using gfreduce for this ensures
    # that the MDF is determined in the same way as in my CL example, though
    # that logic seems to be duplicated in gprepare (argh).
    result = run_task('gemini.gmos.gfreduce', inputs=task_inputs,
        outputs={'outimages' : out_names}, prefix=prefix, comb_in=False,
        MEF_ext=False, path_param='rawpath', reprocess=reprocess,
        outpref='default', slits='header', exslits='*', fl_nodshuffl=False,
        fl_inter=False, fl_vardq=gemvars['vardq'], fl_addmdf=True,
        fl_over=False, fl_trim=False, fl_bias=False, fl_gscrrej=False,
        fl_gnsskysub=False, fl_extract=False, fl_gsappwave=False,
        fl_wavtran=False, fl_skysub=False, fl_fluxcal=False, fl_fulldq=True,
        dqthresh=0.1, key_mdf='', key_biassec=gemvars['key_biassec'],
        key_datasec=gemvars['key_datasec'], 
        bpmfile=gemvars['gmosdata']+'chipgaps.dat', grow=1.5, reference='',
        response='', wavtraname='', sfunction='', extinction='',
        fl_fixnc=False, fl_fixgaps=True, fl_novlap=True, perovlap=10.0,
        nbiascontam='default', biasrows='default', order='default',
        low_reject=3.0, high_reject=3.0, niterate=2, line=iraf.INDEF, nsum=10,
        trace=False, recenter=False, thresh=200., function='chebyshev',
        t_order=21, t_nsum=10, weights='none',
        gratingdb=gemvars['gmosdata']+'GMOSgratings.dat',
        filterdb=gemvars['gmosdata']+'GMOSfilters.dat', xoffset=iraf.INDEF,
        expr='default', sepslits=False, w1=iraf.INDEF, w2=iraf.INDEF,
        dw=iraf.INDEF, nw=iraf.INDEF, observatory=gemvars['observatory'],
        sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], verbose=gemvars['verbose'], **optargs)

    return result['outimages']


@ndprocess_defaults
def subtract_bias(inputs, out_names=None, ovs_function='spline3', ovs_order=1,
                  ovs_lsigma=2.0, ovs_hsigma=2.0, ovs_niter=5, reprocess=None,
                  interact=None):

    """
    Subtract overscan level & pixel-to-pixel variations in zero point. This
    also (as an artifact of using gfreduce) converts the units to electrons.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images from which to subtract the bias & overscan levels and
        trim off the overscan region. Currently these must reside in the
        current working directory (which will normally be the case after
        running prepare). Each input image should already have an entry named
        'bias' in its dictionary of associated calibration files (`cals`
        attribute); these bias images should already have their overscan levels
        subtracted and overscan columns removed.

    out_names : `str`-like or list of `str`-like, optional
        Names of output bias-subtracted files. If None (default), the names
        of the DataFile instances returned will be constructed from those of
        the input files, prefixed with 'r' as in the Gemini IRAF package.

    ovs_function : str
        Function to use for fitting the overscan region in IRAF (default
        'spline3'; may also be 'chebyshev, 'legendre' or 'spline1').

    ovs_order : int
        Order of the overscan fitting function (default 1).

    ovs_lsigma : float
        Negative sigma rejection threshold for overscan fitting (default 2.0).

    ovs_hsigma : float
        Positive sigma rejection threshold for overscan fitting (default 2.0).

    ovs_niter : int
        Number of rejection iterations for overscan fitting (default 5).

    interact : bool, None
        Fit the overscan region interactively in IRAF? If None (default),
        interactivity is instead controlled by the package configuration
        dictionary (see below).

    See "help gfreduce" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The bias-subtracted images produced by gfreduce.


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

    interact : bool
        Enable interactive plotting (default False)? This may be overridden
        by the task's own "interact" parameter.

    """

    # Use default prefix if output filename unspecified:
    prefix = 'r'
    if not out_names:
        out_names = '@inimages'

    # Get a list of biases to use from the input file "cals" dictionaries.
    # If the entry exists, assume the value is already a valid DataFile.
    try:
        bias = DataFileList(data=[df.cals['bias'] for df in inputs])
    except KeyError:
        raise KeyError('one or more inputs is missing an associated bias')

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # Need to add QE correction parameters when copying this, once available.

    # Here some "unnecessary" parameters are defined just to be certain they
    # don't change anything and so we can easily copy this when wrapping
    # other operations with gfreduce later:
    result = run_task('gemini.gmos.gfreduce', inputs={'inimages' : inputs,
        'bias' : bias}, outputs={'outimages' : out_names}, prefix=prefix,
        comb_in=False, MEF_ext=False, path_param='rawpath',
        reprocess=reprocess, outpref='default', slits='header', exslits='*',
        fl_nodshuffl=False, fl_inter=interact, fl_vardq=gemvars['vardq'],
        fl_addmdf=False, fl_over=True, fl_trim=True, fl_bias=True,
        fl_gscrrej=False, fl_gnsskysub=False, fl_extract=False,
        fl_gsappwave=False, fl_wavtran=False, fl_skysub=False,
        fl_fluxcal=False, fl_fulldq=True, dqthresh=0.1, key_mdf='',
        mdffile='default', mdfdir=gemvars['gmosdata'],
        key_biassec=gemvars['key_biassec'], key_datasec=gemvars['key_datasec'],
        bpmfile=gemvars['gmosdata']+'chipgaps.dat', grow=1.5, reference='',
        response='', wavtraname='', sfunction='', extinction='',
        fl_fixnc=False, fl_fixgaps=True, fl_novlap=True, perovlap=10.0,
        nbiascontam='default', biasrows='3:64', order=ovs_order,
        low_reject=ovs_lsigma, high_reject=ovs_hsigma, niterate=ovs_niter,
        line=iraf.INDEF, nsum=10, trace=False, recenter=False, thresh=200.,
        function=ovs_function, t_order=21, t_nsum=10, weights='none',
        gratingdb=gemvars['gmosdata']+'GMOSgratings.dat',
        filterdb=gemvars['gmosdata']+'GMOSfilters.dat', xoffset=iraf.INDEF,
        expr='default', sepslits=False, w1=iraf.INDEF, w2=iraf.INDEF,
        dw=iraf.INDEF, nw=iraf.INDEF, observatory=gemvars['observatory'],
        sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], verbose=gemvars['verbose'])

    return result['outimages']


@ndprocess_defaults
def extract_spectra(inputs, out_names=None, startpos=None, threshold=1000.,
                    order=21, overlap_buffer=0.1, reprocess=None,
                    interact=None):

    """
    Extract one spectrum per IFU fibre from each 2D spectrogram to produce
    row-stacked 1D spectra, mosaicking the separate GMOS CCDs beforehand and
    adding an approximate wavelength solution to the headers afterwards. The
    results are also divided by any associated flat spectra, where available.

    The combination of multiple operations here is a side effect of wrapping
    gfextract; at least the mosaicking & flat fielding are likely to be
    separated into their own steps in a later version.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images from which to extract spectra (in which the individual
        detectors are not yet mosaicked). If flat fielding is to be performed,
        every input DataFile must already have an entry named 'flat' in its
        `cals` dictionary, corresponding to a processed detector+IFU flat
        (otherwise, flat fielding is omitted if there are no such
        associations). Likewise, if the fibre traces are to be taken from one
        or more reference files(s) (usually flats), rather than re-measured,
        every input must have an entry named 'trace' in its `cals` dictionary,
        pointing to a previous output file from this step.

    out_names : `str`-like or list of `str`-like, optional
        Names of output files containing extracted & flat-fielded spectra. If
        None (default), the names of the DataFile instances returned will be
        constructed from those of the input files, prefixed with 'e' as in the
        Gemini IRAF package.

    startpos : int, optional
        Starting column in which the peaks corresponding to each fibre are
        identified (when not using a reference image), before tracing the
        corresponding spectra in either direction along the wavelength axis.
        This is defined relative to the equally-sized sub-region(s) cut out by
        gfextract for the IFU slits (see the log file), rather than the full
        input image. The default is the middle column of each region (ie. of
        the useful wavelength range). This parameter is typically useful for
        avoiding localized artifacts, such as cosmic rays, that can be mistaken
        by apall for fibres when bright enough. The definition of this
        parameter is likely to change in a future version, such that the value
        will be an offset relative to the default column (currently awkward to
        implement, as the regions cut out by gfextract are unknown beforehand).

    threshold : float, optional
        Feature detection threshold for centring around initial peak positions.
        This is the minimum difference required between the lowest & highest
        pixel values in the surrounding region in order for the peak to be
        considered a feature (see IRAF center1d).

    order : int, optional
        Order of the fit used to trace each fibre centre as a function of
        position along the dispersion axis.

    overlap_buffer : float, optional
        Amount by which to extend the removal of any region where the slits
        overlap in 2-slit mode, as a fraction of the overlap (gfextract
        perovlap parameter). The default of 0.1 will remove 110% of the nominal
        overlap, while negative values will include some of the overlap region.

    interact : bool, None
        Identify the fibres interactively in IRAF? If None (default),
        interactivity is instead controlled by the package configuration
        dictionary (see below).

    See "help gfextract" in IRAF for more detailed information. Note that the
    behaviour of this step can depend on whether an IRAF database file already
    exists from a previous run, causing the earlier results to be re-used.


    Returns
    -------

    outimage : DataFileList
        The extracted spectra produced by gfextract.


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

    interact : bool
        Enable interactive plotting (default False)? This may be overridden
        by the step's own "interact" parameter.

    """

    # Initialize the task inputs dict from the main list of input files; any
    # associated "reference" or "response" file lists then get added below,
    # only if applicable:
    task_inputs = {'inimage' : inputs}

    # Use default prefix if output filename unspecified:
    prefix = 'e'
    if not out_names:
        out_names = '@inimage'

    # Use apall starting column default if not specified:
    if startpos is None:
        startpos = iraf.INDEF

    # Convert overlap parameter to IRAF convention, deliberately allowing
    # negative values etc.
    perovlap = 10. * overlap_buffer

    # Make a list of flats from each input's "cals" dictionary. If NO files
    # have flats associated, the response parameter is omitted and no flat
    # fielding is done. Otherwise, at present, they must all have one.
    flats = [df.cals['flat'] if 'flat' in df.cals else None for df in inputs]
    if None in flats:
        if not all([entry is None for entry in flats]):
            raise KeyError('one or more inputs is missing an associated flat')
    else:
        task_inputs['response'] = DataFileList(data=flats)

    # Do likewise for the reference traces.
    traces = [df.cals['trace'] if 'trace' in df.cals else None for df in inputs]
    if None in traces:
        trace = True
        if not all([entry is None for entry in traces]):
            raise KeyError('one or more inputs is missing an associated '\
                           'reference trace')
    else:  # have existing traces
        trace = False
        task_inputs['reference'] = DataFileList(data=traces)

    # To do:
    # - Replace the above replicated blocks with a DataFileList cals attribute
    #   that returns another DataFileList? Allow for None values in
    #   DataFileList so files can remain associated positionally when there are
    #   missing associations (probably non-trivial)?
    # - Consider adapting run_task to accept None values for inputs and set
    #   the corresponding filename parameter to "" in that case?
    #   - An alternative might be to add a cals parameter to run_task, which
    #     would tell it to derive each named parameter from a specified cal
    #     type in the inputs cal dicts -- but this could be messy, as it means
    #     mapping each cal parameter to both an inputs parameter and a cal type
    #     rather than a single DataFileList.
    # - Some additional parameters probably want exposing, such as perovlap.

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # The apextract trace parameter seems to control both whether a trace in
    # wavelength is used at all (as opposed to a fixed range of rows around the
    # measured starting centre, with no tilt or curvature) and whether it is
    # repeated for the input spectra when using a reference image (leading to
    # artifacts in science data), so the appropriate setting depends on whether
    # a reference is used. The recenter parameter controls whether a fixed
    # adjustment is applied to the fibre centre zero points when using a
    # reference image -- but that would break the correspondence with the pixel
    # flat, so leave it disabled for now (probably better to lose a bit of
    # signal for small shifts than to introduce unquantified systematics).

    # Change the default function to spline3 (etc.), like in the CL example?

    result = run_task(
        'gemini.gmos.gfextract', inputs=task_inputs,
        outputs={'outimage' : out_names}, prefix=prefix, comb_in=False,
        MEF_ext=False, path_param=None, reprocess=reprocess, outpref='e',
        title='', exslits='*', line=startpos, nsum=10, trace=trace,
        recenter=False, thresh=threshold, function='spline3', order=order,
        t_nsum=10, weights='none', bpmfile=gemvars['gmosdata']+'chipgaps.dat',
        grow=1.5, gaindb='default',
        gratingdb=gemvars['gmosdata']+'GMOSgratings.dat',
        filterdb=gemvars['gmosdata']+'GMOSfilters.dat', xoffset=iraf.INDEF,
        perovlap=perovlap, sci_ext=labels['data'],
        var_ext=labels['uncertainty'], dq_ext=labels['flags'],
        fl_inter=interact, fl_vardq=gemvars['vardq'], fl_novlap=True,
        fl_gnsskysub=False, fl_fixnc=False, fl_fixgaps=True, fl_gsappwave=True,
        fl_fulldq=True, dqthresh=0.1, verbose=gemvars['verbose']
    )

    return result['outimage']


@ndprocess_defaults
def calibrate_wavelength(inputs, order=4, line_list=None, interact=None):
    """
    Determine an absolute wavelength solution for each row (fibre) of the
    input spectra. Currently, this information is saved only in an IRAF
    database, whence it is picked up directly by the rectification step.

    In future, this will be replaced by two Python steps that use AstroPy
    gWCS to propagate a two-component wavelength solution, one mapping the
    differences between fibres and the other a reference spectrum from pixels
    to wavelength units.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images with extracted, row-stacked spectra, from which to
        determine a wavelength fit by comparision with the `line_list`.

    order : int, optional
        Order of the one-dimensional Chebyshev fitting function, used to
        describe wavelength as a function of pixel index for each row
        (default 4).

    line_list : str, optional
        Name of a text file from which reference wavelengths are to be read
        for the spectral lines present in the source spectrum (see the help
        for IRAF task "identify"). The default is "gmos$data/CuAr_GMOS.dat".

    interact : bool, None
        Inspect & modify the line identifications and fitting functions
        interactively in IRAF? If None (default), interactivity is instead
        controlled by the package configuration  dictionary (see below).

    See "help gswavelength" in IRAF for more detailed information.


    Returns
    -------

    DataFileList
        The (unmodified) `inputs`, with newly-associated wavelength solution
        databases in IRAF.


    Package 'config' options
    ------------------------

    interact : bool
        Enable interactive plotting (default False)? This may be overridden
        by the step's own "interact" parameter.

    """

    # For the time being, the output from this only exists in the IRAF
    # database on disk and is thereby communicated between gswavelength &
    # gftransform. When the step is eventually replaced with Python code,
    # the solution will persist in the gWCS attribute(s) of the output
    # NDLater instances.

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    if line_list is None:
        line_list = gemvars['gmosdata']+'CuAr_GMOS.dat'

    fl_inter = 'yes' if interact else 'NO'  # enumerated value in this task

    result = run_task(
        'gemini.gmos.gswavelength', inputs={'inimages' : inputs}, outputs=None,
        prefix=None, suffix=None, comb_in=False, MEF_ext=False, path_param=None,
        crval='CRVAL1', cdelt='CD1_1', crpix='CRPIX1', key_dispaxis='DISPAXIS',
        dispaxis=1, database='database', coordlist=line_list,
        gratingdb=gemvars['gmosdata']+'GMOSgratings.dat',
        filterdb=gemvars['gmosdata']+'GMOSfilters.dat', fl_inter=fl_inter,
        section='default', nsum=1, ftype='emission', fwidth=10., gsigma=0.,
        cradius=10., threshold=0., minsep=2.5, match=-6., function='chebyshev',
        order=order, sample='*', niterate=10., low_reject=2.5, high_reject=2.5,
        grow=0., refit=True, step=1, trace=True, nlost=10, maxfeatures=150,
        ntarget=30, npattern=5, fl_addfeat=True, aiddebug='s', fl_dbwrite='YES',
        fl_overwrite=True, fl_gsappwave=False, fitcfunc='chebyshev',
        fitcxord=4, fitcyord=4, verbose=gemvars['verbose']
    )

    # Propagate the task inputs unmodified since the solution currently isn't
    # being attached to them but will be in future.
    return inputs


@ndprocess_defaults
def rectify_wavelength(inputs, out_names=None, start_wl=None, end_wl=None,
                       delta_wl=None, npix=None, reprocess=None):
    """
    Resample ("transform") fibre spectra onto a pixel array with a linear
    wavelength increment.

    This step currently depends on a version of gftransform that has been
    patched with respect to the public version 1.13.1, to accept ".fits" in
    the filenames passed to its `wavtraname` parameter. This change is
    available in the "ifudrgmos" version of the GMOS package (r145+) and is
    likely to be included in the public Gemini IRAF package in 2017.


    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing extracted fibre spectra, to be interpolated
        and resampled to linear World co-ordinates along the wavelength axis.
        Each input must already have an associated entry named 'arc' in its
        `cals` dictionary (whose name is used by gftransform, in the current
         implementation, to find the corresponding IRAF database files).

    out_names : `str`-like or list of `str`-like, optional
        Names of output images with a constant increment in wavelength per
        pixel. If None (default), the names of the DataFile instances returned
        will be constructed from those of the input files, prefixed with 't' as
        in the Gemini IRAF package.

    start_wl, end_wl, delta_wl : float, optional
        Starting & ending wavelengths (at the centre of the first/last column)
        and wavelength increment per pixel to use for the output grid, in the
        same units as the database (? Usually Angstroms). By default, the
        output spans the same average wavelength range as the input, with the
        same number of pixels. Any combination of these parameters and `npix`
        can be overridden (precedence if all four specified TBC), leaving the
        remainder to be determined based on the input.

    npix : int, optional
        Length in pixels of the wavelength axis of the output image(s), by
        default the same as in the input (used in conjuction with `start_wl`,
        `end_wl` & `delta_wl` to constrain the wavelength increment & range).

    See "help gftransform" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The resampled spectra produced by gftransform, with a linear
        wavelength increment.


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
    prefix = 't'
    if not out_names:
        out_names = '@inimages'

    # Make a list of arcs from each input's "cals" dictionary. Every file must
    # have one associated. The wavtraname parameter for gftransform is an
    # oddball case because the files it names don't actually need to exist
    # (just the corresponding database entries) and the original task fails if
    # .fits extensions are included -- had to fix it in order to use a
    # DataFileList, otherwise run_task won't correlate multiple inputs properly.
    try:
        arcs = DataFileList(data=[df.cals['arc'] for df in inputs])
    except KeyError:
        raise KeyError('one or more inputs is missing an associated arc')

    # The reference arcs must all have been saved to disk and not be flagged
    # as dirty, otherwise run_task will save them with a temp filename that
    # doesn't match the database. We could guard against this by saving the
    # files and resetting the "unloaded" (dirty) flag here, though it would
    # want restoring again afterwards, since external references that can be
    # used to modify the data may already exist. This problem will go away as
    # and when the "dirty flag logic" can be made more sophisticated.

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # Convert unspecified values to IRAF syntax:
    if start_wl is None: start_wl = iraf.INDEF
    if end_wl is None: end_wl = iraf.INDEF
    if delta_wl is None: delta_wl = iraf.INDEF
    if npix is None: npix = iraf.INDEF

    result = run_task(
        'gemini.gmos.gftransform',
        inputs={'inimages' : inputs, 'wavtraname' : arcs},
        outputs={'outimages' : out_names}, prefix=prefix, suffix=None,
        comb_in=False, MEF_ext=False, path_param=None, reprocess=reprocess,
        database='database', w1=start_wl, w2=end_wl, dw=delta_wl, nw=npix,
        dqthresh=0.1, fl_vardq=gemvars['vardq'], fl_flux='yes',
        sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], verbose=gemvars['verbose']
    )

    return result['outimages']


@ndprocess_defaults
def make_flat(inputs, flat_names=None, order=45, sample='*', reprocess=None,
              interact=None):
    """
    Generate a normalized flat field calibration spectrum that includes both
    fibre-to-fibre and pixel-to-pixel variations, by fitting the average
    continuum and dividing each fibre spectrum by the result (scaled to match
    its wavelength solution).

    This requires the modified version of the IRAF task gfresponse found in
    the "ifudrgmos" version of the GMOS package (r145+), published at
    http://drforum.gemini.edu/topic/gmos-ifu-data-reduction-scripts/.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, each containing extracted fibre spectra of a spatially
        flat, spectrally smooth source (normally the GCal Quartz-Halogen lamp),
        with known wavelength solutions. Each input must already have an
        associated entry named 'arc' in its `cals` dictionary (whose name is
        currently used to find the corresponding IRAF database files).

    flat_names : `str`-like or list of `str`-like, optional
        Names of output normalized flat-field images, one per input. If None
        (default), the names of the DataFile instances returned will be
        constructed from those of the input files by appending '_flat'.
        [To do: should this combine the inputs & produce only one output?]

    order : int, optional
        Order of the Chebyshev continuum fit (default 45).

    sample : str, optional
        The range of sample pixels in the wavelength dimension to use for
        continuum fitting (currently IRAF convention, default '*').

    See "help gfresponse" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The normalized flat field produced by gfresponse.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    interact : bool
        Enable interactive continuum fitting (default False)? This may be
        overridden by the step's own "interact" parameter.

    """

    # Default to appending "_flat" if an output filename is not specified:
    if not flat_names:
        flat_names = '@inimage'

    # Make a list of arcs from each input's "cals" dictionary. Every file must
    # have one associated.
    try:
        arcs = DataFileList(data=[df.cals['arc'] for df in inputs])
    except KeyError:
        raise KeyError('one or more inputs is missing an associated arc')

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    result = run_task(
        'gemini.gmos.gfresponse',
        inputs={'inimage' : inputs, 'wavtraname' : arcs},
        outputs={'outimage' : flat_names}, prefix=None, suffix='_flat',
        comb_in=False, MEF_ext=False, path_param=None, reprocess=reprocess,
        title='', skyimage='', database='database', fl_inter=interact,
        fl_fit=False, function='chebyshev', order=order, sample=sample,
        sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], verbose=gemvars['verbose']
    )

    return result['outimage']


@ndprocess_defaults
def subtract_sky(inputs, out_names=None, reprocess=None):
    """
    Create an average sky spectrum over rows corresponding to the background
    IFU field and subtract it from each fibre spectrum (image row). This is
    done separately for each half of the fields (IFU "slit") in 2-slit mode,
    to match the CCD gaps, wavelength range, detector regions & slit
    characteristics optimally.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing extracted, row-stacked fibre spectra with
        linearized wavelength co-ordinates.

    out_names : `str`-like or list of `str`-like, optional
        Names of output images, containing the sky-subtracted spectra and the
        1D sky spectrum used for each slit [currently only on disk until
        DataFile propagates "extras" properly]. If None (default), the names
        of the DataFile instances returned will be constructed from those of
        the input files, prefixed with 's' as in the Gemini IRAF package.

    See "help gfskysub" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The sky-subtracted spectra produced by gfskysub.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """
    # Use default prefix if output filename unspecified:
    prefix = 's'
    if not out_names:
        out_names = '@inimages'

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # Disable interactivity as it seems slow and not that useful.
    # Note that running gfskysub via gfreduce has been observed to corrupt
    # the MDF in the past, requiring it to be re-copied from the input to the
    # output, but so far so good here.
    result = run_task(
        'gemini.gmos.gfskysub',
        inputs={'inimages' : inputs}, outputs={'outimages' : out_names},
        prefix=prefix, suffix=None, comb_in=False, MEF_ext=False,
        path_param=None, reprocess=reprocess, apertures="", expr='default',
        combine='average', reject='avsigclip', scale='none', zero='none',
        weight='none', sepslits=True, lthreshold=iraf.INDEF,
        hthreshold=iraf.INDEF, nlow=1, nhigh=1, nkeep=0, mclip=True, lsigma=3.,
        hsigma=3., key_ron=gemvars['key_ron'], key_gain=gemvars['key_gain'],
        snoise="0.0", sigscale=0.1, pclip=-0.5, grow=0.0, blank=0.0,
        sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], fl_inter=False, verbose=gemvars['verbose']
    )

    return result['outimages']


@ndprocess_defaults
def resample_to_cube(inputs, out_names=None, bitmask=8, use_uncert=None,
                     reprocess=None):
    """
    Resample fibre spectra onto a 3D "data cube", interpolating spatially over
    the hexagonal/triangular grid of the IFU. The alignment between wavelength
    planes is also corrected for estimated atmospheric dispersion as part of
    the process, using a model from SLALIB.


    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing extracted, row-stacked fibre spectra with
        linearized wavelength co-ordinates.

    out_names : `str`-like or list of `str`-like, optional
        Names of output images. If None (default), the names of the DataFile
        instances returned will be constructed from those of the input files,
        prefixed with 'd' as in the Gemini IRAF package.

    bitmask : int, optional
        Data quality bits used to exclude bad pixels from the interpolation,
        where available in the input file (default 8, for cosmic rays). Where
        multiple data quality bits are summed to produce a bitmask that is not
        a power of two, any one or more of those bits in the data quality plane
        will cause the corresponding pixel to be excluded. A value of 0 causes
        all the input pixels to be used. Since this step interpolates only in
        the spatial directions, features of approximately constant wavelength,
        such as chip gaps (16) and bad columns (1) are better dealt with by
        spectral interpolation at an earlier reduction step.

    See "help gfcube" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The 3D images produced by gfcube, each containing a 3D "data cube" with
        0.1" pixels spatially and the flux units of the input converted to
        values per square arcsecond.


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
    prefix = 'd'
    if not out_names:
        out_names = '@inimage'

    # Always generate output DQ since it tracks the interpolation bounds.
    result = run_task(
        'gemini.gmos.gfcube',
        inputs={'inimage' : inputs}, outputs={'outimage' : out_names},
        prefix=prefix, suffix=None, comb_in=False, MEF_ext=False,
        path_param=None, reprocess=reprocess, ssample=0.1, bitmask=bitmask,
        fl_atmdisp=True, fl_flux=True, fl_var=use_uncert, fl_dq=True
    )

    return result['outimage']


@ndprocess_defaults
def sum_spectra(inputs, out_names=None, reprocess=None):
    """
    Sum over fibre spectra in the main "object" IFU field to produce a 1D
    output spectrum.


    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing extracted, row-stacked fibre spectra with
        linearized wavelength co-ordinates.

    out_names : `str`-like or list of `str`-like, optional
        Names of output images, each containing a 1D spectrum. If None
        (default), the names of the DataFile instances returned will be
        constructed from those of the corresponding input files, prefixed with
        'a' as in the Gemini IRAF package.

    See "help gfapsum" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The 1D spectra produced by gfapsum.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """
    # Use default prefix if output filename unspecified:
    prefix = 'a'
    if not out_names:
        out_names = '@inimages'

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # Use a simple summation, at least for now, as getting the rejection etc.
    # right with the level of contrast involved can be tricky.
    result = run_task(
        'gemini.gmos.gfapsum',
        inputs={'inimages' : inputs}, outputs={'outimages' : out_names},
        prefix=prefix, suffix=None, comb_in=False, MEF_ext=False,
        path_param=None, reprocess=reprocess, apertures='', expr='default',
        combine='sum', reject='none', scale='none', zero='none', weight='none',
        lthreshold=iraf.INDEF, hthreshold=iraf.INDEF, nlow=1, nhigh=1, nkeep=0,
        mclip=True, lsigma=3., hsigma=3., key_ron=gemvars['key_ron'],
        key_gain=gemvars['key_gain'], snoise='0.0', sigscale=0.1, pclip=-0.5,
        grow=0.0, blank=0.0, sci_ext=labels['data'],
        var_ext=labels['uncertainty'], dq_ext=labels['flags'], fl_inter=False,
        verbose=gemvars['verbose']
    )

    return result['outimages']


@ndprocess_defaults
def background_regions(input_ref):
    """
    Determine unilluminated regions of the detector, which can subsequently
    be used (with optional user modifications) to estimate an instrumental
    background or scattered light level as a function of position.

    Parameters
    ----------

    input_ref : DataFileList or DataFile
        Input image with the dimensions of the raw data. This must have an
        entry named 'trace' in its `cals` dictionary, from which the
        illuminated and unilluminated detector regions can be determined.

    See "help gffindblocks" in IRAF for more detailed information.


    Returns
    -------

    mask : list of list of int
        List of rectangular regions in the format [[x1, x2, y1, y2], ...]
        (currently FITS convention).

    """

    # Get the reference traces from the input.
    try:
        trace = input_ref.cals['trace']
    except KeyError:
        raise KeyError('input_ref is missing an associated reference trace')

    # Generate a temp filename in which to store the output mask:
    outfn = new_filename(base=trace.filename.base+'_gaps', ext='')

    run_task(
        'gemini.gmos.gffindblocks',
        inputs={'image' : input_ref, 'extspec' : trace}, outputs=None,
        prefix=None, suffix=None, comb_in=False, MEF_ext=False,
        path_param=None, reprocess=None, mask=outfn
    )

    # Read & parse the regions from the output text file:
    regions = []
    with open(outfn, 'r') as outfile:
        for line in outfile:
            line = line.strip()
            if not line or line[0] == '#':
                continue
            # fields = [int(field)-1 for field in line.split()]  # NumPy conv?
            fields = [int(field) for field in line.split()]
            if len(fields) != 4:
                raise ValueError('unexpected format for gffindblocks output '\
                                 '{0}'.format(outfn))

            # regions.append(fields[2:4] + fields[0:2])  # use NumPy conv?
            regions.append(fields)

    # Don't need the output file after reading it into memory:
    os.remove(outfn)

    return regions


@ndprocess_defaults
def subtract_bg(inputs, out_names=None, x_order=None, y_order=None,
                reprocess=None):
    """
    Model the instrumental background or scattered light level as a function
    of position in the input files (based on the counts within specified
    nominally-unilluminated regions) and subtract the result from each input
    to remove the estimated contamination.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input, bias-subtracted images, in the raw data format, each of which
        must have an entry named 'bg_reg' in its `cals` dictionary, specifying
        the unilluminated detector regions to use for background estimation;
        see ``background_regions``.

    out_names : `str`-like or list of `str`-like, optional
        Names of output images containing the background-subtracted spectra. If
         None (default), the names of the DataFile instances returned will be
        constructed from those of the corresponding input files, prefixed with
        'b' as in the Gemini IRAF package.

    x_order, y_order : int or list of int, optional
        Order of the Legendre surface fit along rows and columns, respectively,
        for each CCD (or all CCDs if a single integer). With the default of
        None, orders of [5,9,5] or [5,5,9,5,5,5] are used for x and [5,7,5] or
        [5,5,7,5,5,5] for columns, as appropriate. The index of the higher
        number may need adjusting by the user to match the CCD where the IFU
        slits overlap (if applicable). This logic will probably be made a bit
        more intelligent in a future version.

    See "help gfscatsub" in IRAF for more detailed information.


    Returns
    -------

    outimage : DataFileList
        The background-subtracted images produced by gfscatsub.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """

    # Here we have to expand out the output filenames instead of letting
    # run_task do it because it currently doesn't recognize text files as main
    # inputs. This should be replaced by run-task-like decorator functionality
    # in the longer run.

    # Convert inputs to a DataFileList if needed:
    inputs = to_datafilelist(inputs)

    # Use default prefix if the output filenames are unspecified:
    prefix = 'b'
    if not out_names:
        out_names = [FileName(indf, prefix=prefix) for indf in inputs]
    elif len(out_names) != len(inputs):
        raise ValueError('inputs & out_names have unmatched lengths')

    # Get lists of bg regions to use from the input file "cals" dictionaries:
    try:
        bg_reg_list = [df.cals['bg_reg'] for df in inputs]
    except KeyError:
        raise KeyError('one or more inputs is missing associated list of '\
                       'background regions')

    # Avoid raising obscure errors if the wrong thing gets attached as bg_reg.
    # To do: consider writing a more generic type-checking function.
    if not all(bg_reg and hasattr(bg_reg, '__iter__') and \
               all(reg and hasattr(reg, '__iter__') and \
                   all(isinstance(n, (int, long, basestring)) for n in reg) \
                 for reg in bg_reg \
               ) for bg_reg in bg_reg_list
           ):
        raise ValueError('cals[\'bg_reg\'] should be a list of limit lists')

    # Loop over the inputs explicitly, since run_task currently can't recognize
    # lists of text files as inputs:
    mode = 'update' if not reprocess else 'overwrite'
    outputs = DataFileList(mode=mode)
    for indf, bg_reg, outname in zip(inputs, bg_reg_list, out_names):

        # Save the background regions for each instance as a temporary file
        # for IRAF:
        gapfn = new_filename(base=indf.filename.base+'_gaps', ext='')
        with open(gapfn, 'w') as gapfile:
            for reg in bg_reg:
                gapfile.write('{0}\n'.format(' '.join(str(n) for n in reg)))

        # Generate default orders appropriate for the number of detectors in
        # each DataFile, if unspecified:
        len_df = len(indf)
        if x_order is None:
            xorder = [5] * len_df
            xorder[(len_df-1)//2] = 9
        else:
            xorder = x_order if hasattr(x_order, '__iter__') else (x_order,)
        if y_order is None:
            yorder = [5] * len_df
            yorder[(len_df-1)//2] = 7
        else:
            yorder = y_order if hasattr(y_order, '__iter__') else (y_order,)

        # Convert list of orders to comma-separated IRAF syntax:
        xorder = ','.join(str(n) for n in xorder)
        yorder = ','.join(str(n) for n in yorder)

        result = run_task(
            'gemini.gmos.gfscatsub',
            inputs={'image' : indf}, outputs={'outimage' : outname},
            prefix=None, suffix=None, comb_in=False, MEF_ext=False,
            path_param=None, reprocess=reprocess, mask=gapfn,
            xorder=xorder, yorder=yorder, cross=True
        )

        # Finished with the temporary file:
        os.remove(gapfn)

        # Accumulate the output DataFileList, copying the dictionary of cals
        # from each input to the output until persistence is implemented, since
        # the same ones are usually needed at the next step:
        outdf = result['outimage'][0]
        outdf.cals.update(indf.cals)
        outputs.append(outdf)

    return outputs


@ndprocess_defaults
def align_wcs(inputs, method='correlate'):
    """
    Measure spatial offsets between IFU data cubes, by comparing their image
    features (summed over wavelength), and update their WCS zero points
    accordingly, to match the first cube in ``inputs``. This provides the
    alignment information needed for subsequent co-addition with ``mosaic``.

    Users should compare the registration of output WCSs carefully (for which
    the ``mosaic.separate`` option may be useful) and adjust the CRVAL keywords
    manually in case of errors.

    This function depends on the ``pyfu`` PyRAF/Python package, which must be
    installed separately.


    Parameters
    ----------

    inputs : DataFileList or DataFile
        Reduced data cubes (normally from ``resample_to_cube`` or Gemini IRAF's
        gfcube), with any cosmic ray flux removed. For meaningful results,
        their fields of view must overlap sufficiently to identify one or more
        spatial peaks in common.

    method : {'correlate', 'centroid'}
        Spatial registration algorithm to use (after collapsing each cube in
        wavelength to form an image). The default of 'correlate' should always
        be used where image structure is nebulous or multi-peaked. The simple
        'centroid' algorithm can be used for sources that have a single 
        well-defined peak (but is unreliable in the presence of cosmic rays).
        In limited testing to date, 'correlate' has been found to work
        comparably well to the original 'centroid' for single-peaked sources,
        so has been made the default in this wrapper.

    Processing is currently performed using the PyFU function "pyfalign".


    Returns
    -------

    DataFileList
        The input images with their WCS zero-points adjusted. Pyfalign modifies
        its inputs, rather than creating a new copy (since no information is
        lost in the process, with the first file retaining the original WCS).
        The files are consequently always re-processed. It is possible that
        this will change in future.

    """

    import pyfu

    inputs = to_datafilelist(inputs)
    names = [str(df) for df in inputs]

    pyfu.pyfalign(names, method=method)

    inputs.reload()  # sync with output saved by pyfalign

    return inputs


@ndprocess_defaults
def mosaic(inputs, out_name=None, separate=False, use_uncert=None,
           reprocess=None):
    """
    Resample and co-add IFU data cubes onto a single, mosaicked output cube
    (with registration determined by their WCS zero-point differences).

    This function depends on the ``pyfu`` PyRAF/Python package, which must be
    installed separately.


    Parameters
    ----------

    inputs : DataFileList or DataFile
        Reduced data cubes (normally from ``resample_to_cube`` or Gemini IRAF's
        gfcube) whose WCS zero points have been adjusted onto a common system
        (normally by ``align_wcs``, PyFU's "pyfalign" or manual determination).

    out_name : `str`-like or list of `str`-like, optional
        Names of output file, containing the datacube mosaic. If None
        (default), the name of the DataFile instance returned will be
        constructed from that of the first input file, with '_add' appended.

    separate : `bool`
        Write one output cube per input to separate NDData arrays in the
        output (separate FITS extensions), instead of co-adding the cubes onto
        a single array (default False)? This option is used for inspecting the
        registration of the resampled component cubes (or to allow co-addition
        with another program).

    Processing is currently performed using the PyFU function "pyfmosaic".


    Returns
    -------

    DataFileList
        The output datacube mosaic.


    Package 'config' options
    ------------------------

    use_uncert : bool
        Enable NDData 'uncertainty' (variance) propagation (default True)?

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """

    import pyfu

    suffix='_add'

    inputs = to_datafilelist(inputs)
    inlist = [str(df) for df in inputs]
    if len(inputs) < 1:
        raise ValueError('inputs cannot be empty')

    if out_name is None:
        out_name = FileName(inlist[0], suffix=suffix, dirname='')
    elif isinstance(out_name, list):
        if len(out_name) == 1:
            out_name = out_name[0]
        else:
            raise ValueError('out_name should have one value, if specified')
    out_name = str(out_name)

    mode = 'new' if reprocess is None else \
           'update' if reprocess is False and os.path.exists(out_name) \
           else 'overwrite'

    outdf = DataFile(out_name, mode=mode)

    if mode != 'update':

        # NB. The posangle option seems to have been broken for a while
        # (rotation about the wrong centre).
        pyfu.pyfmosaic(inlist, outimage=out_name, posangle=None,
                       separate=separate, propvar=use_uncert)

        outdf.reload()

    return outdf


@ndprocess_defaults
def log_rebin(inputs, out_names=None, conserve_flux=False, use_uncert=None,
              reprocess=None):
    """
    Resample the wavelength axis of each x-y-lambda data cube onto constant
    increments in log(wavelength), ie. constant redshift per pixel. The output
    WCS is defined according to the FITS paper III convention (equation 5).

    This function depends on the ``pyfu`` PyRAF/Python package, which must be
    installed separately.


    Parameters
    ----------

    inputs : DataFileList or DataFile
        Data cubes with linear WCSs in wavelength (each one mapped to a FITS
        extension named 'SCI' for compatibility with PyFU), as produced by
        ``resample_to_cube`` or  ``mosaic``.

    out_names : `str`-like or list of `str`-like, optional
        Names of output files containing the re-binned datacubes. If None
        (default), the names of the DataFile instances returned will be
        constructed from those of the corresponding input files, prefixed
        with 'l'.

    conserve_flux : `bool`
        Conserve flux per pixel? The default of False conserves flux density,
        which is appropriate for data that have been calibrated in physical
        flux units (per Angstrom, nm or micron), as opposed to counts or
        electrons per pixel.


    Processing is currently performed using the PyFU function "pyflogbin".


    Returns
    -------

    DataFileList
        The output data cubes, resampled to logarithmic increments in
        wavelength.


    Package 'config' options
    ------------------------

    use_uncert : bool
        Enable NDData 'uncertainty' (variance) propagation (default True)?

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """

    import pyfu

    inputs = to_datafilelist(inputs)  # To do: check if these are unsaved
    if len(inputs) < 1:
        raise ValueError('inputs cannot be empty')

    # Use default prefix if the output filenames are unspecified:
    prefix = 'l'
    if out_names is None:
        out_names = [FileName(indf, prefix=prefix) for indf in inputs]
    elif len(out_names) != len(inputs):
        raise ValueError('inputs & out_names have unmatched lengths')

    # Create a list to hold the outputs:
    outputs = DataFileList(mode='overwrite')

    # Loop over the input/output pairs:
    for indf, outname in zip(inputs, out_names):

        outname = str(outname)

        mode = 'new' if reprocess is None else \
               'update' if reprocess is False and os.path.exists(outname) \
                else 'overwrite'

        outdf = DataFile(outname, mode=mode)

        if mode != 'update':

            # Use PyFU to resample the saved cube:
            pyfu.pyflogbin(str(indf), outname, fluxcons=conserve_flux,
                           propvar=use_uncert)

            outdf.reload()

        outputs.append(outdf)

    return outputs

