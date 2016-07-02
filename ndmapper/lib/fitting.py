# Copyright(c) 2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import math
import numpy as np
from astropy.modeling import models, fitting

__all__ = ['fit_1D']

function_map = {
    'chebyshev' : models.Chebyshev1D,
    'legendre' : models.Legendre1D,
    'polynomial' : models.Polynomial1D,
}

def fit_1D(image, function='legendre', order=1, axis=-1, lsigma=3.0, hsigma=3.0,
           iterations=0):
    """
    An routine for evaluating the result of fitting 1D polynomials to each
    vector along some axis of an N-dimensional image array, with iterative
    pixel rejection and re-fitting, similar to IRAF's fit1d.

    Only a subset of fit1d functionality is currently supported (not
    including interactivity).

    Parameters
    ----------

    image : `ndarray`
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
        npix = image.shape[axis]
    except IndexError:
        raise ValueError('axis={0} out of range for input shape {1}'\
                         .format(axis, image.shape))

    # Record input dtype so we can cast the evaluated fits back to it, since
    # modelling always seems to return float64:
    intype = image.dtype

    # To support fitting any axis of an N-dimensional array, first convert the
    # input array to a 2D stack of 1D rows (a no-op for 2D images with axis=-1)
    # -- make a list of dimension numbers, move the specified dim to the end,
    # transpose the array accordingly, and then flatten all dimensions prior to
    # the last one into a single axis:
    ndim = len(image.shape)
    dims = range(ndim)
    dims.append(dims.pop(axis))
    image = image.transpose(dims)
    newshape = image.shape
    image = image.reshape(-1, npix)

    # If applicable, make a temp copy of the image for cleaned pixels:
    clean = image if iterations == 0 else image.copy()

    # Prepare to perform the fits "simultaneously" with AstroPy modelling:
    nfits = image.shape[0]
    points = np.arange(npix, dtype=np.int16)
    points_2D = np.tile(points, (nfits,1))
    models = function_map[function](degree=order-1, n_models=nfits)
    fitter = fitting.LinearLSQFitter()

    # Create a mask for tracking rejected pixel values:
    mask_2D = np.zeros_like(image, dtype=np.bool)

    # (Re-)fit the data with rejection of outlying points:
    for n in range(iterations+1):

        # Fit the pixel data:
        models = fitter(models, points, clean)

        # Evaluate the fits at each pixel:
        fitvals = models(points_2D).astype(intype)

        # Replace deviant pixels in each row with the fitted values:
        for n, (row, fit, mask) in enumerate(zip(clean, fitvals, mask_2D)):
            diff = row - fit
            stddev = math.sqrt(np.mean(diff*diff))
            rejpix = np.where((diff > hsigma*stddev) | (diff < -lsigma*stddev))
            mask[rejpix] = True
            row[mask] = fit[mask]  # replace cumulative rejections with new fit
            lastdev = stddev

    # # TEST: Plot the fit:
    # import pylab
    # nrow = 3087
    # mask = mask_2D[nrow]
    # pylab.plot(points, image[nrow], 'k.')
    # pylab.plot(points, fitvals[nrow])
    # pylab.plot(points[mask], image[nrow][mask], 'rx')
    # pylab.show()

    # Restore the original ordering & shape of the (substitute) array:
    fitvals = fitvals.reshape(newshape)
    dims = range(ndim)
    dims.insert(axis, dims[-1])  # insert before pop to handle [-1] properly
    del dims[-1]
    fitvals = fitvals.transpose(dims)

    return fitvals

