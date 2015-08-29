# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from pyraf import iraf
from ndmapper import config, ndprocess_defaults
from ndmapper.iraf_task import run_task
from ndmapper.modes.gemini import gemini_iraf_helper

@ndprocess_defaults
def prepare(inputs, outputs=None, mdf=None):
    """
    Parameters
    ----------

    images : DataFileList or DataFile
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


    Package 'config' options
    ------------------------

    data_name : str
        Name identifying main NDData data (science) arrays within a file
        (default 'SCI' = Gemini convention).

    """

    # Use default prefix if output filename unspecified:
    prefix = 'g'
    if not outputs:
        outputs = '@inimages'

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Get MDF from gmos$data by default or, if specified, from the CWD:
    if mdf is None:
        mdf = 'default'
        mdfdir = gemvars['gmosdata']
    else:
        mdfdir = ''

    # Here some "unnecessary" parameters are defined just to be certain they
    # don't change anything and so we can easily copy this when wrapping
    # other operations with gfreduce later:
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
        sci_ext=config['data_name'], var_ext=config['uncertainty_name'],
        dq_ext=config['flags_name'], verbose=gemvars['verbose'])

    return result['outimages']


@ndprocess_defaults
def subtract_bias(inputs, bias, outputs=None, ovs_function='spline3',
    ovs_order=1, ovs_lsigma=2.0, ovs_hsigma=2.0, ovs_niter=5, interact=None):

    """
    Parameters
    ----------

    images : DataFileList or DataFile
        Input images from which to subtract the bias & overscan levels and
        trim off the overscan region. Currently these must reside in the
        current working directory (which will normally be the case after
        running prepare).

    bias : DataFile or DataFileList
        Processed bias image to be used for subtracting pixel-to-pixel
        variations in zero point (which should likewise have its overscan
        level subtracted and the overscan region removed).

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

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

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
        sci_ext=config['data_name'], var_ext=config['uncertainty_name'],
        dq_ext=config['flags_name'], verbose=gemvars['verbose'])

    return result['outimages']

