import numpy as np
import scipy.interpolate


def resample_response_curves(
        wavelength, rcurves, resolution, kind='slinear'):
    """Resample the given response curve to specified spectral resolution

    Parameters
    ----------
    wavelength : ndarray shape(nvalues)
        wavelength
    rcurve : ndarray shape(nbands, nvalues)
        sensor response curve
    resolution : float
        resolution to interpolate to
    kind : str
        interpolation algorithm for
        scipy.interpolate.interp1d

    Returns
    -------
    wavelength : ndarray
        new wavelength
    rcurves : ndarray
        resampled rcurves
    """
    f = scipy.interpolate.interp1d(
            wavelength, rcurves, kind=kind, axis=1,
            bounds_error=False, fill_value=0)
    start_wv = wavelength[0]
    end_wv = wavelength[-1]
    nsteps = round((end_wv - start_wv) / resolution) + 1
    xnew = np.linspace(start_wv, end_wv, nsteps)
    return xnew, f(xnew)
