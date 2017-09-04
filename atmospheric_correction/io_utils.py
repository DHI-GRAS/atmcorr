from __future__ import division
import numpy as np


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
            ext.append([x,y])
        yarr.reverse()
    return ext
