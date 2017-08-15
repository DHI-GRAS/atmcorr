

def getTileExtents(img, tileSize):
    # Split the image into square tiles of tileSize pixels and returns the extents
    # ([minX, maxY, maxX, minY]) of each tile in the projected units.
    tileExtents = []
    gt = img.GetGeoTransform()

    rows = [int(i*tileSize) for i in range(int(img.RasterYSize/tileSize)+1)]
    rows.append(img.RasterYSize)
    cols = [int(i*tileSize) for i in range(int(img.RasterXSize/tileSize)+1)]
    cols.append(img.RasterXSize)

    for y in range(1, len(rows)):
        tileExtents.append([])
        for x in range(1, len(cols)):
            ext = getExtent(gt, cols[x]-1, rows[y]-1, cols[x-1], rows[y-1])
            tileExtents[y-1].append([ext[0][0], ext[0][1], ext[2][0], ext[2][1]])

    return tileExtents


def getExtent(gt, endCol, endRow, startCol=0, startRow=0):
    ''' Return list of corner coordinates from a geotransform

        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type endCol:   C{int}
        @param endCol: ending column of the subset relative to gt origin
        @type endRow:   C{int}
        @param endRow: ending row of the subset relative to gt origin
        @type startCol:   C{int}
        @param startCol: starting column of the subset relative to gt origin
        @type startRow:   C{int}
        @param startRow: starting row of the subset relative to gt origin
        @rtype:    C{[float,...,float]}
        @return:   coordinates of each corner
    '''
    ext = []
    xarr = [startCol,endCol]
    yarr = [startRow,endRow]
    for px in xarr:
        for py in yarr:
            x = gt[0]+(px*gt[1])+(py*gt[2])
            y = gt[3]+(px*gt[4])+(py*gt[5])
            ext.append([x,y])
        yarr.reverse()
    return ext
