from __future__ import division

import affine
import pyproj
import numpy as np

from atmcorr import resampling


def get_tiled_transform_shape(src_transform, src_shape, dst_res):
    """Get transform and shape of tile grid with resolution dst_res

    Paramters
    ---------
    src_transform : affine.Affine
        source transform
    src_shape : int, int
        source shape
    dst_res : float or tuple (float, float)
        destination resolution

    Returns
    -------
    affine.Affine
        target transform
    tuple (int, int)
        target shape
    """
    src_res = np.array((src_transform.a, src_transform.e))
    scale = np.abs(dst_res / src_res)
    dst_transform = src_transform * affine.Affine(scale[0], 0, 0, 0, scale[1], 0)
    dst_shape = tuple(np.ceil(np.array(src_shape) / scale).astype(int))
    return dst_transform, dst_shape


def _get_corner_coordinates(transform, height, width):
    """Get coordinates of all four pixel corners of an image of given transform and shape

    Parameters
    ----------
    transform : affine.Affine
        image transform
    height, width : int
        image shape

    Returns
    -------
    ndarray of shape (2, 4, height, width)
        x, y corner coordinates
        ul, ur, lr, ll
    """
    # j index top-first to get bottom-up image with negative transform.e
    i = np.arange(width + 1)
    j = np.arange(height + 1)[::-1]
    jj, ii = np.meshgrid(j, i, indexing='ij')
    xx, yy = transform * (ii, jj)
    ul = np.stack((xx[:-1, :-1], yy[:-1, :-1]), axis=0)
    ur = np.stack((xx[:-1, 1:], yy[:-1, 1:]), axis=0)
    lr = np.stack((xx[1:, 1:], yy[1:, 1:]), axis=0)
    ll = np.stack((xx[1:, :-1], yy[1:, :-1]), axis=0)
    corners = np.zeros((2, 4, height, width))
    corners[:, 0, ...] = ul
    corners[:, 1, ...] = ur
    corners[:, 2, ...] = lr
    corners[:, 3, ...] = ll
    return corners


def _transform_corners(corners, src_crs, dst_crs):
    """Transform corners from array indices to dst_crs coordinates

    Parameters
    ----------
    corners : ndarray shape(2, N, ...) dtype(int)
        x,y pairs for N corners
    src_crs : dict or rasterio.crs.CRS
        source coordinate reference system
    dst_crs : dict or rasterio.crs.CRS
        destination coordinate reference system

    Returns
    -------
    ndarray, ndarray
        projected coordinates
    """
    src_proj = pyproj.Proj(src_crs)
    dst_proj = pyproj.Proj(dst_crs)
    xs, ys = corners
    xout, yout = pyproj.transform(src_proj, dst_proj, xs, ys)
    return xout, yout


def _corners_to_extents(xs, ys):
    """Convert arrays of corner coordinates to an extent record array

    Parameters
    ----------
    xs, ys : ndarray shape(N, ...)
        x and y coordinates of N corners

    Returns
    -------
    np.recarray shape(...)
        xmin, xmax, ymin, ymax
    """
    extent_rec = np.core.records.fromarrays(
            [np.min(xs, axis=0), np.max(xs, axis=0), np.min(ys, axis=0), np.max(ys, axis=0)],
            names=['xmin', 'xmax', 'ymin', 'ymax'])
    return extent_rec


def get_projected_extents(transform, height, width, src_crs, dst_crs={'init': 'epsg:4326'}):
    """Get extents of pixels in WGS84 or other projection

    Parameters
    ----------
    transform : affine.Affine
        image transform
    height, width : int
        image shape
    src_crs : dict or rasterio.crs.CRS
        source coordinate reference system
    dst_crs : dict or rasterio.crs.CRS
        destination coordinate reference system
        default: WGS84 (lon, lat)

    Returns
    -------
    np.recarray shape(...)
        xmin, xmax, ymin, ymax
    """
    corners = _get_corner_coordinates(transform, height, width)
    xproj, yproj = _transform_corners(corners, src_crs, dst_crs=dst_crs)
    return _corners_to_extents(xproj, yproj)


def bounds_to_projected_extents(left, bottom, right, top, src_crs, dst_crs={'init': 'epsg:4326'}):
    """Get extents record array from bounds

    Parameters
    ----------
    left, bottom, right, top : float
        extents
    src_crs, dst_crs : dict
        source and destination coordinate reference systems

    Returns
    -------
    np.recarray shape (1, 1)
        with names xmin, xmax, ymin, ymax
    """
    src_proj = pyproj.Proj(src_crs)
    dst_proj = pyproj.Proj(dst_crs)
    xs = np.array([left, left, right, right])
    ys = np.array([bottom, top, top, bottom])
    xproj, yproj = pyproj.transform(src_proj, dst_proj, xs, ys)
    return _corners_to_extents(xproj, yproj)[np.newaxis, np.newaxis]


def get_projected_image_extent(transform, height, width, src_crs, dst_crs={'init': 'epsg:4326'}):
    """Get extents of a whole image in WGS84 or other projection

    Parameters
    ----------
    transform : affine.Affine
        image transform
    height, width : int
        image shape
    src_crs : dict or rasterio.crs.CRS
        source coordinate reference system
    dst_crs : dict or rasterio.crs.CRS
        destination coordinate reference system
        default: WGS84 (lon, lat)

    Returns
    -------
    np.recarray shape (1, 1)
        with names xmin, xmax, ymin, ymax
    """
    left, top = transform * (0, 0)
    right, bottom = transform * (height, width)
    return bounds_to_projected_extents(
        left, bottom, right, top, src_crs, dst_crs=dst_crs)


def recarr_take_dict(a, *idx):
    return dict(zip(a.dtype.names, a[idx]))


def resample_correction_params(corrparams, src_transform, src_crs, dst_transform, dst_shape):
    """Resampled 6S correction parameters from one per tile to one per image pixel

    Parameters
    ----------
    corrparams : np.recarray shape(nbands, njtiles, nitiles)
        correction parameters
    src_transform : affine.Affine
        tiling grid transform
    src_crs : dict
        source (and destination) spatial reference system
    dst_transform : affine.Affine
        destination image transform
    dst_shape : tuple (int, int)
        destination shape

    Returns
    -------
    np.recarray
        resampled correction parameters
    """
    nj, ni = corrparams.shape[1:]
    ntiles = nj * ni

    resampled = np.empty(dst_shape, dtype=corrparams.dtype)
    for field in resampled.dtype.names:
        if ntiles == 1:
            # take single value
            resampled[field][:] = corrparams[field][0, 0]
        else:
            # interpolate
            resampled[field] = resampling.resample(
                source=corrparams[field],
                src_transform=src_transform,
                src_crs=src_crs,
                dst_shape=dst_shape,
                dst_transform=dst_transform,
                dst_crs=src_crs,
                resampling=None)
    return resampled
