from __future__ import division

import numpy as np

import rasterio.windows
import rasterio.warp


def get_tiled_windows(height, width, tilesize):
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
    windows = np.empty((ny, nx), dtype=object)
    for j in range(ny):
        for i in range(nx):
            col_offset = i * tilesize
            row_offset = j * tilesize
            w = rasterio.windows.Window(col_offset, tilesize, row_offset, tilesize)
            if j == (ny == 1) or i == (nx - 1):
                w = rasterio.windows.crop(w, height, width)
            windows[j, i] = w
    return windows


def _transform_bounds(w, src_transform, src_crs, dst_crs):
    src_bounds = rasterio.windows.bounds(w, src_transform)
    return rasterio.warp.transform_bounds(src_crs, dst_crs, *src_bounds)


def get_transformed_bounds(windows, src_transform, src_crs, dst_crs):
    """Get bounds for a sequence of windows in dst_crs

    Parameters
    ----------
    windows : ndarray of rasterio.windows.Window
        windows
    src_transform : Affine
        source image transform
    src_crs, dst_crs : rasterio.crs.CRS or similar
        source and destination coordinate
        reference systems
    """
    newshape = (4, ) + windows.shape
    bounds = np.zeros(newshape, int)
    winiter = windows.flat
    for _ in range(windows.size):
        w = next(winiter)
        j, i = winiter.coords
        bounds[:, j, i] = _transform_bounds(w, src_transform, src_crs, dst_crs)
    return bounds


def getTileExtents(width, height, transform, tileSize):
    """ Split the image into square tiles of tileSize pixels

    Parameters
    ----------
    width, height : int
        nx, ny shape of image
    transform : affine.Affine
        image transform
    tileSize : int
        number of pixels in each square tile

    Returns
    -------
    extents : ([minX, maxY, maxX, minY])
        extents of each tile in the projected units.
    """
    gt = transform.to_gdal()

    nx = int(width // tileSize) + 1
    ny = int(height // tileSize) + 1
    rows = [int(i * tileSize) for i in range(ny)]
    cols = [int(i * tileSize) for i in range(nx)]
    rows.append(height)
    cols.append(width)

    nrows = len(rows)
    ncols = len(cols)
    tileExtents = np.zeros((nrows, ncols, 4))
    for j in range(1, nrows):
        for i in range(1, ncols):
            ext = getExtent(gt, cols[i]-1, rows[j]-1, cols[i-1], rows[j-1])
            tileExtents[j-1, i, :] = [ext[0][0], ext[0][1], ext[2][0], ext[2][1]]

    return tileExtents


def getExtent(gt, endCol, endRow, startCol=0, startRow=0):
    """Get list of corner coordinates from geotransform

    Parameters
    ----------
    gt : sequence of float
        GDAL geotransform
    endCol : int
        ending column of the subset relative to gt origin
    endRow : int
        ending row of the subset relative to gt origin
    startCol : int
        starting column of the subset relative to gt origin
    startRow : int
        starting row of the subset relative to gt origin

    Returns
    -------
    list of (float, float) : x,y coordinates of each corner
    """
    ext = []
    xarr = [startCol, endCol]
    yarr = [startRow, endRow]
    for px in xarr:
        for py in yarr:
            x = gt[0] + (px * gt[1]) + (py * gt[2])
            y = gt[3] + (px * gt[4]) + (py * gt[5])
            ext.append([x, y])
        yarr.reverse()
    return ext
