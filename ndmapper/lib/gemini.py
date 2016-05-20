# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os.path

from ndmapper import config, ndprocess_defaults
from ndmapper.data import FileName
from ndmapper.libutils import map_API_enum
from ndmapper.iraf_task import run_task

__all__ = ['CAL_DEPS', 'gemini_iraf_helper', 'clean_pixels']


# Default cal dependencies are defined by instrument mode, with just an empty
# placeholder dict at this level:
CAL_DEPS = {'target' : []}


def gemini_iraf_helper():
    """
    Define a few common default values for Gemini IRAF tasks in one place
    (just to avoid more repetition of static values than necessary) and
    convert some Python defaults to be more appropriate for the IRAF tasks.

    Returns
    -------

    dict
        Common parameter values suitable for passing to Gemini IRAF tasks.

    """

    # Enable both VAR & DQ if either is specified, since most Gemini tasks
    # can't control them individually:
    if config['use_uncert'] or config['use_flags']:
        vardq = True
    else:
        vardq = False

    vals = {'gmosdata' : 'gmos$data/',
            'key_airmass' : 'AIRMASS',
            'key_biassec' : 'BIASSEC',
            'key_ccdsec'  : 'CCDSEC',
            'key_ccdsum'  : 'CCDSUM',
            'key_datasec' : 'DATASEC',
            'key_detsec'  : 'DETSEC',
            'key_exptime' : 'EXPTIME',
            'key_gain' : 'GAIN',
            'key_ron' : 'RDNOISE',
            'mask_table_name' : 'MDF',
            'observatory' : 'default',
            'vardq' : vardq,
            'verbose' : True
           }

    return vals


@ndprocess_defaults
def clean_pixels(inputs, out_names=None, method='global', grow=1.5,
                 bitmask=65535, axis=1, order=None, low_reject=3.,
                 high_reject=2.3, iterations=5, reprocess=None, interact=None):
    """
    Replace each pixel whose corresponding `flags` value matches one or more
    non-zero bits of the `bitmask` with a (locally- or globally-) interpolated
    estimate of its true, uncontaminated value.

    This currently works on 2D images.

    Parameters
    ----------

    inputs : DataFileList or DataFile
        Input images, containing extracted, row-stacked fibre spectra with
        linearized wavelength co-ordinates.

    out_names : `str`-like or list of `str`-like, optional
        Names of output images, each containing a 1D spectrum. If None
        (default), the names of the DataFile instances returned will be
        constructed from those of the corresponding input files, prefixed
        with 'p'.

    method : {'local', 'global'}
        Method to use for interpolating good data to generate replacement
        values; either 'local' interpolation along the narrowest dimension of
        each contiguous bad region using IRAF's ``proto.fixpix`` (with
        ``linterp=INDEF`` and ``cinterp=INDEF``) or a 'global' fit to each row
        or column using ``fit1d``. The default is 'global'.

    grow : float, optional
        The radius in pixels (default 1.5) by which to expand the rejection of
        regions matching the ``bitmask`` in the ``flags`` array when
        generating replacement values. This does not cause additional pixels
        to be replaced, it merely avoids basing replacement values on the
        immediately-surrounding pixels, where those are also contaminated at
        a lower level.

    bitmask : int, optional
        The bit-wise OR of ``flags`` bits used to trigger pixel replacement.
        The default of 65535 causes all pixels with DQ > 0 in the input to be
        replaced in the output, while bitmask=0 would copy the input unchanged,
        with intermediate values used to reject some defects but not others
        (eg. 9 flags cosmic rays (8) and detector defects (1) in Gemini IRAF).

    axis : {0, 1}
        The image axis (Python convention) along which to fit 1D Chebyshev
        models when using method 'global' (default 1, ie. rows).

    order : int or None, optional
        The order of the 1D Chebyshev fits when using method 'global'. With the
        default of None, a value is selected automatically by 'gemfix'.

    low_reject, high_reject : int, optional
        Lower and upper thresholds, in standard deviations, for rejection of
        nominally-good pixels (in addition to  those excluded by the ``flags``
        array) when performing 'global' fits, to help avoid any residual
        contamination when generating replacement values. These default to
        3.0 and 2.3, respectively (avoiding the statistically high values
        associated with common defects such as cosmic rays more aggressively
        than low ones).

    iterations : int, optional
        The number of pixel rejection and re-fitting iterations to perform to
        obtain a final fit when using the 'global' method (default 5).


    See "help gemfix" in IRAF for more detailed information.


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

    interact : bool
        Enable interactive fitting for the 'global' method (default False)?
        This may be overridden by the step's own "interact" parameter.

    """
    # Use default prefix if output filename unspecified. We have to do some
    # extra work here, normally done by run_task, to allow for the IRAF bug
    # fix below to behave according to the reprocess flag.
    prefix = 'p'
    if not out_names:
        out_names = [FileName(indf, prefix=prefix) for indf in inputs]
    elif len(out_names) != len(inputs):
        raise ValueError('inputs & out_names have unmatched lengths')

    # Map some Python API values to equivalent IRAF ones:
    method = map_API_enum('method', method, \
                          {'local' : 'fixpix', 'global' : 'fit1d'})
    axis = map_API_enum('axis', axis, {0 : 2, 1 : 1})
    order = order if order else 0

    # Get a few common Gemini IRAF defaults.:
    gemvars = gemini_iraf_helper()

    # Keep copies of any input MDF tables, because an IRAF bug mangles them
    # when running gemfix, causing subsequent steps to fail (eg. with only one
    # of the two slits getting processed).
    mdfs = {}
    for indf, outname in zip(inputs, out_names):
        filename = indf.filename.orig
        mdfs[filename] = None
        if reprocess or not os.path.exists(str(outname)):
            for tp in indf._tables:  # to do: add/use a public API
                if str(tp.label).upper() == 'MDF':
                    mdfs[filename] = tp.table
                    break

    # Use a simple summation, at least for now, as getting the rejection etc.
    # right with the level of contrast involved can be tricky.
    result = run_task(
        'gemini.gemtools.gemfix',
        inputs={'inimages' : inputs}, outputs={'outimages' : out_names},
        prefix=prefix, suffix=None, comb_in=False, MEF_ext=False,
        path_param=None, reprocess=reprocess, method=method, grow=grow,
        bitmask=bitmask, axis=axis, order=order, low_reject=low_reject,
        high_reject=high_reject, niterate=iterations, fl_inter=interact,
        verbose=gemvars['verbose']
    )

    outlist = result['outimages']

    # Restore the original MDF extensions to the output files, if applicable:
    for outdf in outlist:
        filename = outdf.filename.orig
        mdf = mdfs[filename]
        if mdf:
            for tp in outdf._tables:
                if str(tp.label).upper() == 'MDF':
                    tp.table = mdfs[filename]
                    outdf.save()
                    break

    return outlist

