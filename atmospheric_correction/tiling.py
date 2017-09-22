from __future__ import division

import numpy as np
import pyproj


def get_tile_corners(height, width, tilesize):
    """Tile image with width and height into windows of size tilesize

    Paramters
    ---------
    height, width : int
        image dimensions
    tilesize : int
        symmetric tile size

    Returns
    -------
    ndarray shape(ny, nx) of rasterio.windows.Window
        windows
    """
    nx = int(width / tilesize + 0.5)
    ny = int(height / tilesize + 0.5)

    jj = np.arange(ny)
    ii = np.arange(nx)
    jmesh, imesh = np.meshgrid(jj, ii, indexing='ij')
    assert jmesh.shape == (ny, nx)

    # tile grid dims, 4 corners, (x, y)
    corners = np.zeros((2, 4, ny, nx))

    bottom = jmesh * tilesize
    left = imesh * tilesize
    top = np.minimum((bottom + tilesize - 1), ny)
    right = np.minimum((left + tilesize - 1), nx)
    # now create four corner points
    corners[:, 0, ...] = np.stack((left, bottom), axis=0)
    corners[:, 1, ...] = np.stack((left, top), axis=0)
    corners[:, 2, ...] = np.stack((right, top), axis=0)
    corners[:, 3, ...] = np.stack((right, bottom), axis=0)
    return corners


def transform_corners(corners, src_transform, src_crs, dst_crs={'init': 'epsg:4326'}):
    """Transform corners from array indices to dst_crs coordinates

    Parameters
    ----------
    corners : ndarray shape(2, ...) dtype(int)
        x,y pairs for N corners
    src_transform : affine.Affine
        source image transform
    src_crs : dict or rasterio.crs.CRS
        source coordinate reference system
    dst_crs : dict or rasterio.crs.CRS
        destination coordinate reference system
        default: WGS84 (lon, lat)

    Returns
    -------
    ndarray, ndarray
        projected coordinates
    """
    src_proj = pyproj.Proj(src_crs)
    dst_proj = pyproj.Proj(dst_crs)
    xs, ys = src_transform * corners
    xout, yout = pyproj.transform(src_proj, dst_proj, xs, ys)
    return xout, yout


def corners_to_extents(xs, ys):
    """Convert arrays of corner coordinates to an extent record array

    Parameters
    ----------
    xs, ys : ndarray shape(N, ...)
        x and y coordinates

    Returns
    -------
    np.recarray shape(...)
        xmin, xmax, ymin, ymax
    """
    extent_rec = np.core.records.fromarrays(
            [xs.min(axis=0), xs.max(axis=0), ys.min(axis=0), ys.max(axis=0)],
            names=['xmin', 'xmax', 'ymin', 'ymax'])
    return extent_rec


def get_tile_extents(height, width, tilesize, src_transform, src_crs):
    """For an image of height, width get extents for tiles of size tilesize

    Parameters
    ----------
    height, width : int
        image shape
    tilesize : int
        symmetric tile size
    src_transform : affine.Affine
        image transform
    src_crs : dict or rasterio.crs.CRS
        source coordinate reference system

    Returns
    -------
    np.recarray shape(...)
        xmin, xmax, ymin, ymax
    """
    corners = get_tile_corners(height, width, tilesize)
    lon, lat = transform_corners(corners, src_transform, src_crs)
    extents_rec = corners_to_extents(lon, lat)
    return extents_rec
