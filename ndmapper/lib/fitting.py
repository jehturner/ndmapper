# Copyright(c) 2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import numpy as np
from astropy.modeling import models, fitting

function_map = {
    'chebyshev' : models.Chebyshev1D,
    'legendre' : models.Legendre1D,
    'polynomial' : models.Polynomial1D,
}

def fit_1D(arr, function='legendre', order=1, axis=-1, lsigma=3.0, hsigma=3.0,
           iterations=0):
    """
    An routine for evaluating the result of fitting 1D polynomials to each
    vector along some axis of an N-dimensional image array, with (not yet)
    iterative pixel rejection and re-fitting, similar to IRAF's fit1d.

    Only a subset of fit1d functionality is currently supported (not
    including interactivity).

    Parameters
    ----------

    arr : `ndarray`
        N-dimensional input array containing the values to be fitted.

    function : {'legendre'}, optional
        Fitting function/model type to be used (current default 'legendre').

    order : `int`, optional
        Order (number of terms or degree+1) of the fitting function
        (default 1).

    axis : `int`, optional
        Array axis along which to perform fitting (Python convention;
        default -1, ie. rows).

    lsigma, hsigma : `float`, optional
        Rejection threshold in standard deviations below and above the mean,
        respectively (default 3.0).

    iterations : `int`, optional
        Number of rejection and re-fitting iterations (default 0, ie. a single
        fit with no iteration).

    Returns
    -------

    `ndarray`
        An array of the same shape as the input, whose values are evaluated
        from the polynomial fits to each 1D vector.

    """

    # Determine how many pixels we're fitting each vector over:
    try:
        npix = arr.shape[axis]
    except IndexError:
        raise ValueError('axis={0} out of range for input shape {1}'\
                         .format(axis, arr.shape))

    # To support fitting any axis of an N-dimensional array, first convert the
    # input array to a 2D stack of 1D rows (a no-op for 2D images with axis=-1)
    # -- make a list of dimension numbers, move the specified dim to the end,
    # transpose the array accordingly, and then flatten all dimensions prior to
    # the last one into a single axis:
    ndim = len(arr.shape)
    dims = range(ndim)
    dims.append(dims.pop(axis))
    arr = arr.transpose(dims)
    newshape = arr.shape
    arr = arr.reshape(-1, npix)

    # Perform the fits "simultaneously" with AstroPy's modelling module:
    nfits = arr.shape[0]
    points = np.arange(npix)
    models = function_map[function](degree=order-1, n_models=nfits)
    fitter = fitting.LinearLSQFitter()
    models = fitter(models, points, arr)

    # Evaluate the fits at each pixel:
    points_2D = np.tile(points, (nfits,1))
    arr = models(points_2D)

    # Restore the original ordering & shape of the (substitute) array:
    arr = arr.reshape(newshape)
    dims = range(ndim)
    dims.insert(axis, dims[-1])  # insert before pop to handle [-1] properly
    del dims[-1]
    arr = arr.transpose(dims)

    return arr

