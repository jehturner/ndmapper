# Copyright(c) 2015-2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os.path
from copy import copy, deepcopy

import numpy as np
from astropy.modeling import models, fitting

from ndmapper import config, ndprocess_defaults
from ndmapper.data import FileName, DataFile, DataFileList, NDLater
from ndmapper.utils import convert_region, to_datafilelist
from .fitting import fit_1D

__all__ = ['init_bpm', 'add_bpm', 'lacosmic_spec', 'clean_cosmic_rays']


@ndprocess_defaults
def init_bpm(reference, regions, convention='numpy', value=None,
             filename=None, reprocess=None):
    """
    Initialize a bad pixel mask array from a string listing the corresponding
    regions in NumPy or IRAF/FITS syntax.


    Parameters
    ----------

    reference : DataFile
        DataFile whose NDLater instances the bad pixel mask should match
        (in number and array dimensions).

    regions : list of list of str
        A list containing one sub-list of bad regions per NDLater instance of
        the reference DataFile. Each bad region is specified as a string in
        NumPy/IRAF syntax, eg. '100:110,50:60', '100,*' or ':,99'. The indices
        can be specified in either NumPy or FITS ordering (see below).

    convention : str
        Indexing convention for the region strings: 'NumPy' (default) or
        'FITS' (case insensitive). The former uses 0-based indexing by
        (..., plane, row, column), similar to C or matrix notation, exclusive
        of the ending row/column, while the latter uses 1-based indexing by
        (column, row, plane, ...), similar to ForTran or IRAF and inclusive
        of the ending row/column (with the same origin in either case).

    value : int, optional
        Mask value for bad regions (default 1).

    filename : str or FileName, optional
        Filename for output DataFile.

    """

    # Use default bad pixel value of 1:
    if value is None:
        value = 1

    # Re-map the flags/data-quality array labels to the main science array for
    # the special purpose of keeping a master bad pixel mask (NDData doesn't
    # allow for .flags without associated .data). This method could do with
    # tidying up a bit, as simply giving labels=config['labels']['flags'] to
    # DataFile without also resetting labels['flags'] causes duplication, due
    # to the condition for including an array in the lists in DataFile.save():
    labels = config['labels'].copy()
    labels['data'] = labels['flags']
    labels['flags'] = None

    # Create (or read from disk, if not reprocessing) the output DataFile.
    # This is provisional, until a more definitive method is implemented,
    # possibly as part of the decorator.
    mode = 'new' if reprocess is None else \
           'update' if reprocess is False and os.path.exists(str(filename)) \
           else 'overwrite'
    output = DataFile(filename, mode=mode, labels=labels)

    if mode == 'update':
        return output

    # Create an empty integer array of the same shape as each reference
    # NDData instance:
    for ndd in reference:
        # Do we want meta=ndd.meta.copy()? Probably not.
        output.append(NDLater(np.zeros(ndd.data.shape, dtype=np.uint16)))

    # Convert each region string from NumPy or FITS region string format to
    # a NumPy slice object and use it to set the specified flag value for the
    # corresponding pixels:
    for ndd, reglist in zip(output, regions):
        for region in reglist:
            slices = convert_region(region, convention)
            ndd.data[slices] = value

    return output


@ndprocess_defaults
def add_bpm(inputs, bpm, out_names=None, reprocess=None):
    """
    Incorporate an existing bad pixel mask into the corresponding flags
    attributes of the input DataFile instances.

    Parameters
    ----------

    inputs : `DataFileList` or `DataFile`
        Input images to which the bad pixel mask is to be added (via logical
        OR with any existing ``flags`` information).

    bpm : `DataFileList` or `DataFile`
        Bad pixel mask, as produced by `init_bpm`, with integer data quality
        flags stored in the main data array of each constituent NDLater
        instance. The length (number of NDLater instances) and shape of each
        constituent array must match the ``inputs``.

    out_names : convertible to `str`, optional
        Output filenames. If None (default), the names of the DataFile
        instances returned will be constructed from those of the input files,
        prefixed with 'k'.


    Returns
    -------

    DataFileList
        Copies of the ``inputs``, with the bad pixel mask values included in
        the flags attribute of each constituent `NDLater` instance.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """

    # Some of this bookeeping should be replaced by a standard wrapper.
    # Also, logging needs to be implemented, as for IRAF.
    inputs = to_datafilelist(inputs)
    if out_names is None:
        out_names = [FileName(df, prefix='k', dirname='') for df in inputs]

    mode = 'new' if reprocess is None else 'overwrite'
    outlist = DataFileList(mode=mode)

    # Here the wrapper would ensure that input/output lengths are the same.

    len_bpm = len(bpm)

    # Loop over the files:
    for in_df, out_name in zip(inputs, out_names):

        if len(in_df) != len_bpm:
            raise ValueError('input {0} and BPM have different lengths'\
                             .format(str(in_df)))

        # Read and use any existing file if not reprocessing:
        if reprocess is False and os.path.exists(str(out_name)):
            out_df = DataFile(out_name, mode='update')

        # Otherwise, OR each NDLater instance of the BPM with the corresponding
        # flags extension of the input file:
        else:
            # Make a copy the inputs, re-using the data by reference since
            # replacing the flags attribute won't affect the inputs.
            out_df = DataFile(out_name, data=in_df, mode=mode)

            # Set flags to the logical OR of the BPM & any existing flags:
            for out_ndd, bpm_ndd in zip(out_df, bpm):

                out_ndd.flags = copy(bpm_ndd.data) if out_ndd.flags is None \
                                else bpm_ndd.data | out_ndd.flags

        outlist.append(out_df)

    return outlist


def lacosmic_spec(input_ndd, x_order=None, y_order=None, sigclip=4.5,
                  sigfrac=0.32, objlim=1.0, niter=5, sepmed=True,
                  cleantype='meanmask'):
    """
    Detect and clean cosmic rays in a 2D wavelength-dispersed image, using the
    well-known LA Cosmic algorithm of van Dokkum (2001)*, as implemented in
    McCully's optimized version for Python, "lacosmicx"+.

    * LA Cosmic: http://www.astro.yale.edu/dokkum/lacosmic
    + lacosmicx: https://github.com/cmccully/lacosmicx

    lacosmicx is an optional dependency, whose absence will cause this
    function to fail with an ImportError. For the time being, the slightly
    modified fork at https://github.com/jehturner/lacosmicx must be used.

    Currently, input and output bad pixel masks are expected to be found in
    the NDDataArray `flags` attribute, but this will likely change to `mask`
    in future (throughout NDMapper).

    Parameters
    ----------

    input_ndd : NDDataArray-like
        Input image in which cosmic rays are to be detected.

    x_order, y_order : int or None, optional
        Order for fitting and subtracting object continuum and sky line models,
        prior to running the main cosmic ray detection algorithm. When None,
        defaults are used, according to the image size (as in the IRAF task
        gemcrspec). When 0, no fit is done.

    sigclip : float, optional
        Laplacian-to-noise limit for cosmic ray detection. Lower values will
        flag more pixels as cosmic rays. Default: 4.5.

    sigfrac : float, optional
        Fractional detection limit for neighboring pixels. For cosmic ray
        neighbor pixels, a lapacian-to-noise detection limit of
        sigfrac * sigclip will be used. Default: 0.32.

    objlim : float, optional
        Minimum contrast between Laplacian image and the fine structure image.
        Increase this value if cores of bright stars are flagged as cosmic
        rays. Default: 1.0.

    niter : int, optional
        Number of iterations of the LA Cosmic algorithm to perform. Default: 5.

    sepmed : boolean, optional
        Use the separable median filter instead of the full median filter.
        The separable median is not identical to the full median filter, but
        they are approximately the same and the separable median filter is
        significantly faster and still detects cosmic rays well. Default: True

    cleantype : {'median', 'medmask', 'meanmask', 'idw'}, optional
        Set which clean algorithm is used:
        'median': An umasked 5x5 median filter
        'medmask': A masked 5x5 median filter
        'meanmask': A masked 5x5 mean filter
        'idw': A masked 5x5 inverse distance weighted interpolation
        Default: "meanmask".

    Returns
    -------

    NDLater
        A cleaned copy of the input, with cosmic ray detections added to
        its `flags` array (with a value of 8).

    """

    # This can fail if the optional dependency is missing:
    from lacosmicx import lacosmicx

    # For this function to be general-purpose, it should use a meta-data
    # abstraction for the gain, read noise & saturation, based on a
    # configuration database of recognized instruments, but until that's
    # implemented they are required to be called GAIN and RNOISE in meta-data
    # and the assumed saturation level is fixed at 65k.
    gain = input_ndd.meta['GAIN']
    read_noise = input_ndd.meta['RDNOISE']
    saturation = 65535.

    # Convert DQ to a boolean array of pixels for lacosmicx to treat as bad
    # from the beginning:
    bitmask = 65535  # do we want to include everything?
    if input_ndd.flags is None:
        inmask = None
    else:
        inmask = (input_ndd.flags & bitmask) > 0

    # Use default orders from gemcrspec (from Bryan):
    ny, nx = input_ndd.shape
    if x_order is None:
        x_order = 9
    if y_order is None:
        y_order = 2 if ny < 50 else 3 if ny < 80 else 5

    # To do: use dispaxis below, rather than -1.

    # Fit and subtract the object spectrum. For some reason, subtracting the
    # model directly from the NDData instance here resets .flags to None.
    if x_order > 0:
        objfit = fit_1D(input_ndd.data, function='legendre', order=x_order,
                        axis=1, lsigma=4.0, hsigma=4.0, iterations=3)
        input_ndd.data -= objfit
    else:
        objfit = np.zeros_like(input_ndd.data)

    # Fit sky lines:
    if y_order > 0:
        skyfit = fit_1D(input_ndd.data, function='legendre', order=y_order,
                        axis=0, lsigma=4.0, hsigma=4.0, iterations=3)
        input_ndd.data -= skyfit
        objfit += skyfit  # keep combined fits for later restoration
        del skyfit

    # Delegate all the actual identification and cleaning to lacosmicx (a
    # version to which I've added a bgsub parameter that allows for a
    # previously-subtracted spectroscopic object+sky model to be included in
    # the noise estimates, as in the original):
    cr_mask, clean_data = lacosmicx(
        input_ndd.data, inmask=inmask, bgsub=objfit, sigclip=sigclip,
        sigfrac=sigfrac, objlim=objlim, gain=gain, readnoise=read_noise,
        satlevel=saturation, pssl=0.0, niter=niter, sepmed=sepmed,
        cleantype=cleantype, fsmode='median', verbose=True
    )

    # Add object & sky signal back in, after cleaning what structure is left:
    clean_data += objfit

    # Obey config options for whether to propagate uncertainty & flags:
    if input_ndd.uncertainty is not None and config['use_uncert']:
        uncertainty = deepcopy(input_ndd.uncertainty)
    else:
        uncertainty = None
    if input_ndd.flags is not None and config['use_flags']:
        # To do: abstract the DQ convention somewhere instead of using "8":
        flags = input_ndd.flags | (np.array(cr_mask, dtype=np.uint16) * 8)
    else:
        flags = None

    # Construct an NDData-like object from the lacosmicx output, plus the
    # original variance and bad pixel mask, and return it,
    return NDLater(data=clean_data, uncertainty=uncertainty, flags=flags,
                   meta=input_ndd.meta)


@ndprocess_defaults
def clean_cosmic_rays(inputs, out_names=None, x_order=None, y_order=None,
                      sigma=4.5, sigfrac=0.32, objlim=1.0, iterations=5,
                      reprocess=None):
    """
    Identify pixels contaminated by cosmic ray flux, using McCully's version
    of the LACosmic algorithm (see `lacosmic_spec`), record the detections in
    a copy of the input `flags` array and replace bad values with an estimate
    of the clean value, based on a masked mean of surrounding pixels. This
    works on 2D images (with fitting disabled) and spectra.


    Parameters
    ----------

    inputs : DataFileList or DataFile
        2D bias-subtracted input images to be cleaned.

    out_names : `str`-like or list of `str`-like, optional
        Names of cleaned output images. If None (default), the names of the
        DataFile instances returned will be constructed from those of the
        corresponding input files, prefixed with 'x', as in Gemini IRAF.

    x_order, y_order : int or None, optional
        Order for fitting and subtracting object continuum and sky line models,
        prior to running the main cosmic ray detection algorithm. When None,
        defaults are used, according to the image size (as in the IRAF task
        gemcrspec). When 0, no fit is done.

    sigma : float, optional
        Laplacian-to-noise limit for cosmic ray detection. Lower values will
        flag more pixels as cosmic rays. Default: 4.5.

    sigfrac : float, optional
        Fractional detection limit for neighboring pixels. For cosmic ray
        neighbour pixels, a lapacian-to-noise detection limit of
        sigfrac * sigclip will be used. Default: 0.32.

    objlim : float, optional
        Minimum contrast between Laplacian image and the fine structure image.
        Increase this value if cores of bright stars are flagged as cosmic
        rays. Default: 1.0.

    iterations : int, optional
        Number of iterations of the detection algorithm to perform. Default: 5.


    Returns
    -------

    DataFileList
        A copy of the input data with cosmic ray detections included in its
        `flags` arrays and the corresponding `data` values replaced by local
        masked mean values.


    Package 'config' options
    ------------------------

    reprocess : bool or None
        Re-generate and overwrite any existing output files on disk or skip
        processing and re-use existing results, where available? The default
        of None instead raises an exception where outputs already exist
        (requiring the user to delete them explicitly). The processing is
        always performed for outputs that aren't already available.

    """
    # Default prefix if output filenames unspecified:
    prefix = 'x'

    inputs = to_datafilelist(inputs)
    if not out_names:
        out_names = [FileName(indf, prefix=prefix) for indf in inputs]
    elif len(out_names) != len(inputs):
        raise ValueError('inputs & out_names have unmatched lengths')

    mode = 'new' if reprocess is None else 'overwrite'
    outlist = DataFileList(mode=mode)

    # For now, GAIN & RDNOISE are hard-wired into lacosmic_spec, but should
    # be abstracted at some point. We don't use gemvars in Python and the
    # config options to use uncertainty & flags are picked up directly by
    # lacosmic_spec().

    # Loop over the files:
    for in_df, out_name in zip(inputs, out_names):

        # Read and use any existing file if not reprocessing:
        # To do: move this sort of logic into shared infrastructure.
        if reprocess is False and os.path.exists(str(out_name)):
            out_df = DataFile(out_name, mode='update')

        # Otherwise do/repeat the cosmic ray detection:
        else:
            # Make a new copy of the input file (copied by reference, as the
            # actual NDData instances get overwritten below and higher-level
            # attributes get deep-copied anyway -- though not any MDF tables):
            out_df = DataFile(out_name, data=in_df, mode=mode)

            # Run lacosmic_spec() on each NDData instance of the current file:
            for n, in_ndd in enumerate(in_df):

                out_df[n] = lacosmic_spec(in_ndd, x_order=x_order,
                                          y_order=y_order, sigclip=sigma,
                                          sigfrac=sigfrac, objlim=objlim,
                                          niter=iterations, sepmed=True,
                                          cleantype='meanmask')

        outlist.append(out_df)

    return outlist

