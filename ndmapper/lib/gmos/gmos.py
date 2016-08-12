# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from pyraf import iraf
from ndmapper import config, ndprocess_defaults
from ndmapper.data import FileName, DataFile, DataFileList
from ndmapper.iraf_task import run_task, get_extname_labels
from ndmapper.utils import to_datafilelist

from ..gemini import *

# Redefine explicitly what to expose from our parent module here, rather than
# appending to its __all__ (as at lower levels of the tree), since it contains
# observatory-specific helper functions that are not processing steps.
__all__ = ['CAL_DEPS', 'make_bias', 'clean_pixels']

# These functions are intended to represent logical processing steps, rather
# than strict one-to-one wrappers for existing IRAF tasks; the aim is not to
# replicate the IRAF user interface exactly. Each "wrapper" function must
# therefore deal with any options exposed to the user explicitly, rather than
# allowing IRAF task parameters to be specified arbitrarily via **params.

# At least to begin with, the aim is to simplify things aggressively and not
# expose every option for atypical usage until a use case arises. When that
# happens, we'll weigh up passing additional options directly vs indirectly vs
# adding more wrappers to support different common usage scenarios.

# NB. PyRAF accepts True/False, iraf.yes/iraf.no & 'yes'/'no' interchangeably.

@ndprocess_defaults
def make_bias(inputs, bias=None, bpm=None, ovs_function='spline3',
    ovs_order=1, ovs_lsigma=2.0, ovs_hsigma=2.0, ovs_niter=5,
    comb_lsigma=2.0, comb_hsigma=2.0, reprocess=None, interact=None):

    """
    Combine individual bias exposures to produce a reference bias frame for
    calibrating other data.


    Parameters
    ----------

    inputs : DataFileList
        Input raw bias images.

    bias : str-like, optional
        Output combined bias image name. If None (default), a new DataFile
        will be returned whose name is constructed from that of the first input
        file, by appending '_bias'.

    bpm : DataFile or DataFileList, optional
        A bad pixel mask, used if 'use_uncert' and/or 'use_flags' is enabled.
        This can be created with the Gemini IRAF task GBPM.

    ovs_function : str
        Function to use for fitting the overscan region in IRAF (default
        'chebyshev'; may also be 'legendre', 'spline1' or 'spline3').

    ovs_order : int
        Order of the overscan fitting function (default 1).

    ovs_lsigma : float
        Negative sigma rejection threshold for overscan fitting (default 2.0).

    ovs_hsigma : float
        Positive sigma rejection threshold for overscan fitting (default 2.0).

    ovs_niter : int
        Number of rejection iterations for overscan fitting (default 5).

    comb_lsigma : float
        Negative sigma rejection threshold for averaging biases (default 2.0).

    comb_hsigma : float
        Positive sigma rejection threshold for averaging biases (default 2.0).

    interact : bool, None
        Fit the overscan region interactively in IRAF? If None (default),
        interactivity is instead controlled by the package configuration
        dictionary (see below).

    See "help gbias" in IRAF for more detailed information.


    Returns
    -------

    outbias : DataFile
        The combined bias image produced by gbias.


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

    # Some candidate parameters to open up to the UI:
    verbose = True

    # Default to appending "_bias" if an output filename is not specified:
    if not bias:
        bias = '!inimages'

    # Determine input DataFile EXTNAME convention, to pass to the task:
    labels = get_extname_labels(inputs)

    # Insert a BPM in the task inputs if supplied by the user
    # (NB. Use of this BPM parameter is untested at the time of writing; it
    # would need a multi-extension FITS BPM in place of the pixel list files
    # distributed with the package):
    inputs = {'inimages' : inputs}
    if bpm:
        inputs['bpm'] = bpm

    # Most of the IRAF package tasks don't have the granularity to control
    # VAR & DQ propagation separately, so just turn them both on if either
    # is specified. This isn't handled by ndprocess_defaults since the
    # Python & IRAF APIs are different (two parameters vs one):
    if config['use_uncert'] or config['use_flags']:
        vardq = True
    else:
        vardq = False

    # Wrap gbias, defining the parameters reproducibly (for a given version)
    # but omitting inapplicable parameters such as minmax options. Certain
    # parameters, such as logfile & rawpath, are set directly by run_task.
    result = run_task('gemini.gmos.gbias', inputs=inputs,
        outputs={'outbias' : bias}, suffix='_bias', comb_in=True,
        MEF_ext=False, path_param='rawpath', reprocess=reprocess, fl_over=True,
        fl_trim=True, key_biassec='BIASSEC', key_datasec='DATASEC',
        key_ron='RDNOISE', key_gain='GAIN', ron=3.5, gain=2.2,
        gaindb='default', sci_ext=labels['data'],
        var_ext=labels['uncertainty'], dq_ext=labels['flags'], sat='default',
        nbiascontam='default', biasrows='default', fl_inter=interact,
        median=False, function=ovs_function, order=ovs_order,
        low_reject=ovs_lsigma, high_reject=ovs_hsigma, niterate=ovs_niter,
        combine='average', reject='avsigclip', lthreshold=iraf.INDEF,
        hthreshold=iraf.INDEF, masktype='goodvalue', maskvalue=0.0,
        scale='none', zero='none', weight='none', statsec='[*,*]',
        key_exptime='EXPTIME', nkeep=1, mclip=True, lsigma=comb_lsigma,
        hsigma=comb_hsigma, sigscale=0.1, grow=0.0, fl_vardq=vardq,
        verbose=verbose)

    # Return the only DataFile instance from the output DataFileList
    # corresponding to the task's "outbias" parameter:
    return result['outbias'][0]


