from __future__ import division
import numpy as np
import scipy.interpolate


def imresize(a, newshape):
    """Interpolate an a to a new shape"""

    # If there is only one value then assign it to each cell in the new array
    if a.size == 1:
        out = np.zeros(newshape)
        out[:] = a[0]
        return out

    # Otherwise interpolate
    # At least two points in each dimension are needed for linear interpolation
    if a.shape[0] == 1:
        a = np.vstack((a, a))
    if a.shape[1] == 1:
        a = np.hstack((a, a))

    # Assume that the values are regularly spaced within the new array
    dx = newshape[1] / a.shape[1]
    xx = (np.arange(a.shape[1]) + 0.5) * dx
    dy = newshape[0] / a.shape[0]
    yy = (np.arange(a.shape[0]) + 0.5) * dy

    xxnew = np.arange(newshape[0])
    yynew = np.arange(newshape[1])

    f = scipy.interpolate.interp2d(xx, yy, a)
    return f(yynew, xxnew)
