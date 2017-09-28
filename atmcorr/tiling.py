from __future__ import division

import numpy as np
import pyproj


def get_tile_corners_ij(height, width, xtilesize, ytilesize):
    """Tile image with width and height into windows of size tilesize

    Paramters
    ---------
    height, width : int
        image dimensions
    xtilesize, ytilesize : int
        tile sizes

    Returns
    -------
    ndarray shape(ny, nx) of rasterio.windows.Window
        windows
    """
    nx = int(width / xtilesize + 0.5)
    ny = int(height / ytilesize + 0.5)

    jj = np.arange(ny)
    ii = np.arange(nx)
    jmesh, imesh = np.meshgrid(jj, ii, indexing='ij')
    assert jmesh.shape == (ny, nx)

    # tile grid dims, 4 corners, (x, y)
    corners = np.zeros((2, 4, ny, nx))

    bottom = jmesh * ytilesize
    left = imesh * xtilesize
    top = np.minimum((bottom + ytilesize - 1), (height - 1))
    right = np.minimum((left + xtilesize - 1), (width - 1))
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
            [np.min(xs, axis=0), np.max(xs, axis=0), np.min(ys, axis=0), np.max(ys, axis=0)],
            names=['xmin', 'xmax', 'ymin', 'ymax'])
    return extent_rec


def get_tile_extents(height, width, src_transform, src_crs, xtilesize, ytilesize):
    """For an image of height, width get extents for tiles of size xtilesize, ytilesize

    Parameters
    ----------
    height, width : int
        image shape
    src_transform : affine.Affine
        image transform
    src_crs : dict or rasterio.crs.CRS
        source coordinate reference system
    xtilesize, ytilesize : int
        tile sizes

    Returns
    -------
    np.recarray shape(...)
        xmin, xmax, ymin, ymax
    """
    corners = get_tile_corners_ij(height, width, xtilesize, ytilesize)
    lon, lat = transform_corners(corners, src_transform, src_crs)
    extents_rec = corners_to_extents(lon, lat)
    return extents_rec


def extents_from_bounds(left, bottom, right, top, src_crs, dst_crs={'init': 'epsg:4326'}):
    """Get extents record array from bounds

    Parameters
    ----------
    left, bottom, right, top : float
        extents
    src_crs, dst_crs : dict
        source and destination coordinate reference systems
    """
    src_proj = pyproj.Proj(src_crs)
    dst_proj = pyproj.Proj(dst_crs)
    xs = [left, left, right, right]
    ys = [bottom, top, top, bottom]
    xout, yout = pyproj.transform(src_proj, dst_proj, xs, ys)
    return corners_to_extents(xs, ys)
