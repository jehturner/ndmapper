# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from pyraf import iraf
from ndmapper import config, ndprocess_defaults
from ndmapper.data import DataFileList
from ndmapper.iraf_task import run_task, get_extname_labels
from ndmapper.modes.gemini import gemini_iraf_helper

from .spec import *

@ndprocess_defaults
def prepare(inputs, outputs=None, mdf=None):
    """
    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images to be "prepared" by attaching a "mask definition file"
        (MDF) extension, describing the slit mapping to sky, and performing
        meta-data updates.

    outputs: DataFileList or DataFile, optional
        Output "prepared" files. If None (default), a new DataFileList
        will be returned, whose names are constructed from those of the input
        files, prefixed with 'g' as in the Gemini IRAF package.

    mdf : str, optional
        Name of a FITS mask definition file to be attached as a FITS extension
        of each output file.

    See "help gfreduce" in IRAF for more detailed information.


    Returns
    -------

    outimages : DataFileList
        The "prepared" images produced by gfreduce.

    """

    # Use default prefix if output filename unspecified:
    prefix = 'g'
    if not outputs:
        outputs = '@inimages'

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # Get MDF from gmos$data by default or, if specified, from the CWD:
    if mdf is None:
        mdf = 'default'
        mdfdir = gemvars['gmosdata']
    else:
        mdfdir = ''

    # Here some "unnecessary" parameters are defined just to be certain they
    # don't change anything and so we can easily copy this when wrapping
    # other operations with gfreduce later. Using gfreduce for this ensures
    # that the MDF is determined in the same way as in my CL example, though
    # that logic seems to be duplicated in gprepare (argh).
    result = run_task('gemini.gmos.gfreduce', inputs={'inimages' : inputs},
        outputs={'outimages' : outputs}, prefix=prefix, comb_in=False,
        MEF_ext=False, path_param='rawpath', outpref='default', slits='header',
        exslits='*', fl_nodshuffl=False, fl_inter=False,
        fl_vardq=gemvars['vardq'], fl_addmdf=True, fl_over=False, fl_trim=False,
        fl_bias=False, fl_gscrrej=False, fl_gnsskysub=False, fl_extract=False,
        fl_gsappwave=False, fl_wavtran=False, fl_skysub=False,
        fl_fluxcal=False, fl_fulldq=True, dqthresh=0.1, key_mdf='',
        mdffile=mdf, mdfdir=mdfdir, key_biassec=gemvars['key_biassec'],
        key_datasec=gemvars['key_datasec'],
        bpmfile=gemvars['gmosdata']+'chipgaps.dat', grow=1.5,
        reference='', response='', wavtraname='', sfunction='', extinction='',
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
        dq_ext=labels['flags'], verbose=gemvars['verbose'])

    return result['outimages']


@ndprocess_defaults
def subtract_bias(inputs, outputs=None, ovs_function='spline3', ovs_order=1,
                  ovs_lsigma=2.0, ovs_hsigma=2.0, ovs_niter=5, interact=None):

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

    outputs: DataFileList or DataFile, optional.
        Output bias subtracted files. If None (default), a new DataFileList
        will be returned, whose names are constructed from those of the input
        files, prefixed with 'r' as in the Gemini IRAF package.

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

    data_name : str
        Name identifying main NDData data (science) arrays within a file
        (default 'SCI' = Gemini convention).

    uncertainty_name : str
        Name identifying NDData uncertainty (variance) arrays within a file
        (default 'VAR' = Gemini convention).

    flags_name : str
        Name identifying NDData flags (data quality) arrays within a file
        (default 'DQ' = Gemini convention).

    interact : bool
        Enable interactive plotting (default False)? This may be overridden
        by the task's own "interact" parameter.

    """

    # Use default prefix if output filename unspecified:
    prefix = 'r'
    if not outputs:
        outputs = '@inimages'

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
        'bias' : bias}, outputs={'outimages' : outputs}, prefix=prefix,
        comb_in=False, MEF_ext=False, path_param='rawpath', outpref='default',
        slits='header', exslits='*', fl_nodshuffl=False, fl_inter=interact,
        fl_vardq=gemvars['vardq'], fl_addmdf=False, fl_over=True, fl_trim=True,
        fl_bias=True, fl_gscrrej=False, fl_gnsskysub=False, fl_extract=False,
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
def extract_spectra(inputs, outputs=None, startpos=None, interact=None):

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

    outputs: DataFileList or DataFile, optional
        Output files containing extracted & flat-fielded spectra. If None
        (default), a new DataFileList will be returned, whose names are
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
    if not outputs:
        outputs = '@inimage'

    # Use apall starting column default if not specified:
    if startpos is None:
        startpos = iraf.INDEF

    # Make a list of flats from each input's "cals" dictionary. If NO files
    # have flats associated, the response parameter is omitted and no flat
    # fielding is done. Otherwise, at present, they must all have one.
    flats = [df.cals['flat'] if 'flat' in df.cals else None for df in inputs]
    if None in flats:
        print [entry for entry in flats]
        if not all([entry is None for entry in flats]):
            raise KeyError('one or more inputs is missing an associated flat')
    else:
        task_inputs['response'] = DataFileList(data=flats)

    # Do likewise for the reference traces.
    traces = [df.cals['trace'] if 'trace' in df.cals else None for df in inputs]
    if None in traces:
        if not all([entry is None for entry in traces]):
            raise KeyError('one or more inputs is missing an associated '\
                           'reference trace')
    else:
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

    # The apextract trace parameter only seems to control whether a trace in
    # wavelength is used at all (as opposed to a fixed range of rows around the
    # measured starting centre, with no tilt or curvature); thus, there is no
    # reason to set it to False. In either case, where a reference is used, the
    # existing trace terms (but not the fibre centre / zero point) get copied
    # from it verbatim. The recenter parameter controls whether a fixed
    # adjustment is applied to the fibre centre zero points when using a
    # reference image -- but that would break the correspondence with the pixel
    # flat, so leave it disabled for now (probably better to lose a bit of
    # signal for small shifts than to introduce unquantified systematics).

    result = run_task(
        'gemini.gmos.gfextract', inputs=task_inputs,
        outputs={'outimage' : outputs}, prefix=prefix, comb_in=False,
        MEF_ext=False, path_param=None, outpref='e', title='', exslits='*',
        line=startpos, nsum=10, trace=True, recenter=False, thresh=200.,
        function='chebyshev', order=5, t_nsum=10, weights='none',
        bpmfile=gemvars['gmosdata']+'chipgaps.dat', grow=1.0, gaindb='default',
        gratingdb=gemvars['gmosdata']+'GMOSgratings.dat',
        filterdb=gemvars['gmosdata']+'GMOSfilters.dat', xoffset=iraf.INDEF,
        perovlap=10., sci_ext=labels['data'], var_ext=labels['uncertainty'],
        dq_ext=labels['flags'], fl_inter=interact, fl_vardq=gemvars['vardq'],
        fl_novlap=True, fl_gnsskysub=False, fl_fixnc=False, fl_fixgaps=True,
        fl_gsappwave=True, fl_fulldq=True, dqthresh=0.1,
        verbose=gemvars['verbose']
    )

