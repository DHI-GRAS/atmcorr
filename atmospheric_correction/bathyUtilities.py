import os

from osgeo import gdal, ogr, osr
import numpy as np

driverOptionsGTiff = ['PREDICTOR=1', 'BIGTIFF=IF_SAFER']  # 'COMPRESS=DEFLATE']


def check_gdal_success(outfile, cmd):
    """Make sure GDAL command `cmd` succeeded in creating `outfile`"""
    if not os.path.isfile(outfile):
        raise RuntimeError('GDAL command \'{}\' did not produce the '
                'expected output file {}.'.format(cmd, outfile))


def mask(img, mask, maskValue=0):
    maskData = mask.GetRasterBand(1).ReadAsArray()
    maskedData = np.zeros((img.RasterYSize, img.RasterXSize, img.RasterCount))
    for band in range(1,img.RasterCount+1):
        maskedData[:,:,band-1] = img.GetRasterBand(band).ReadAsArray()
    maskedData[maskData == maskValue,:] = 0
    res = saveImg(maskedData, img.GetGeoTransform(), img.GetProjection(), "MEM")
    return res


def rasterize(in_vector, out_raster, pixel_size=25):
    # Input: Path to a shape file, filename of output raster and pixel size
    # Output: GDAL object with rasterised shape file
    # Define pixel_size and NoData value of new raster
    NoData_value = -9999

    # Filename of the raster Tiff that will be created
    raster_fn = out_raster

    # Open the data source and read in the extent
    source_ds = in_vector
    source_layer = source_ds.GetLayer()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()

    # Create the destination data source
    x_res = int((x_max - x_min) / pixel_size)
    y_res = int((y_max - y_min) / pixel_size)
    if out_raster != 'MEM':
        target_ds = gdal.GetDriverByName('GTiff').Create(raster_fn, x_res, y_res, 1, gdal.GDT_Byte)
    else:
        target_ds = gdal.GetDriverByName('MEM').Create(raster_fn, x_res, y_res, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(NoData_value)

    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])

    source_ds = None
    return target_ds


def calcMinCoveringExtent(imgA, imgB):
    # Input: two GDALDataset
    # Output: minimum covering extent of the two datasets
    # [minX, maxY, maxX, minY] (UL, LR)
    aGeoTrans = imgA.GetGeoTransform()
    bGeoTrans = imgB.GetGeoTransform()
    minX = max(aGeoTrans[0], bGeoTrans[0])
    maxY = min(aGeoTrans[3], bGeoTrans[3])
    maxX = min(aGeoTrans[0]+imgA.RasterXSize*aGeoTrans[1], bGeoTrans[0]+imgB.RasterXSize*bGeoTrans[1])
    minY = max(aGeoTrans[3]+imgA.RasterYSize*aGeoTrans[5], bGeoTrans[3]+imgB.RasterYSize*bGeoTrans[5])
    return [minX, maxY, maxX, minY]


def clipRasterWithShape(rasterImg, shapeImg):
    # Input: raster file and shape file
    # Output: GDAL object with the clipped raster
    # get the raster pixels size and extent
    #rasterImg = gdal.Open(raster,gdal.GA_ReadOnly)
    rasterGeoTrans = rasterImg.GetGeoTransform()

    # rasterize the shape and get its extent
    shapeRasterImg = rasterize(shapeImg, 'MEM', rasterGeoTrans[1])
    shapeRasterGeoTrans = shapeRasterImg.GetGeoTransform()

    # make sure that raster and shapeRaster pixels are co-aligned
    ulX = rasterGeoTrans[0] + round((shapeRasterGeoTrans[0]-rasterGeoTrans[0])/rasterGeoTrans[1])*rasterGeoTrans[1]
    ulY = rasterGeoTrans[3] + round((shapeRasterGeoTrans[3]-rasterGeoTrans[3])/rasterGeoTrans[5])*rasterGeoTrans[5]
    shapeRasterGeoTrans = (ulX, shapeRasterGeoTrans[1], shapeRasterGeoTrans[2], ulY, shapeRasterGeoTrans[4], shapeRasterGeoTrans[5])

    # get minimum covering extent for the raster and the shape and covert them
    # to pixels
    [minX, maxY, maxX, minY] = calcMinCoveringExtent(rasterImg, shapeRasterImg)
    rasterSubsetPixs         = world2Pixel(rasterGeoTrans, minX, maxY)      + world2Pixel(rasterGeoTrans, maxX, minY)
    shapeRasterSubsetPixs    = world2Pixel(shapeRasterGeoTrans, minX, maxY) + world2Pixel(shapeRasterGeoTrans, maxX, minY)


    # clip the shapeRaster to min covering extent
    shapeRasterData = shapeRasterImg.GetRasterBand(1).ReadAsArray()
    shapeRasterClipped = shapeRasterData[shapeRasterSubsetPixs[1]:shapeRasterSubsetPixs[3], shapeRasterSubsetPixs[0]:shapeRasterSubsetPixs[2]]

    # go through the raster bands, clip the to the minimum covering extent and mask out areas not covered by vector
    maskedData = np.zeros((np.shape(shapeRasterClipped)[0], np.shape(shapeRasterClipped)[1], rasterImg.RasterCount))
    bandNum = rasterImg.RasterCount
    for band in range(bandNum):
        rasterData = rasterImg.GetRasterBand(band+1).ReadAsArray()
        clippedData = rasterData[rasterSubsetPixs[1]:rasterSubsetPixs[3], rasterSubsetPixs[0]:rasterSubsetPixs[2]]
        maskedData[:, :, band] = np.where(shapeRasterClipped > 0, clippedData, np.NaN)

    # get the geotransform array for the masekd array
    maskedGeoTrans = (ulX, rasterGeoTrans[1], rasterGeoTrans[2], ulY, rasterGeoTrans[4], rasterGeoTrans[5])

    # save the masked img to memory and return it
    maskedImg = saveImg(maskedData, maskedGeoTrans, rasterImg.GetProjection(), "MEM", np.NaN)
    return maskedImg


def openAndClipRaster(inFilename, shapeRoiFilename):
    inImg = gdal.Open(inFilename, gdal.GA_ReadOnly)

    # If the ROI file is not specified or does not exist then return
    # unclipped image
    if not shapeRoiFilename or not os.path.exists(shapeRoiFilename):
        return inImg

    shapeRoi = ogr.Open(shapeRoiFilename)
    clippedImg = clipRasterWithShape(inImg, shapeRoi)
    inImg = None
    shapeRoi = None
    return clippedImg


def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    pixel = int(round((x - ulX) / xDist))
    line = int(round((ulY - y) / yDist))
    return [pixel, line]


# save the data to geotiff or memory
def saveImg(data, geotransform, proj, outPath, noDataValue=np.nan):
    # Start the gdal driver for GeoTIFF
    if outPath == "MEM":
        driver = gdal.GetDriverByName("MEM")
        driverOpt = []
    else:
        driver = gdal.GetDriverByName("GTiff")
        driverOpt = driverOptionsGTiff

    shape=data.shape
    if len(shape) > 2:
        ds = driver.Create(outPath, shape[1], shape[0], shape[2], gdal.GDT_Float32, driverOpt)
        ds.SetProjection(proj)
        ds.SetGeoTransform(geotransform)
        for i in range(shape[2]):
            ds.GetRasterBand(i+1).WriteArray(data[:,:,i])
            ds.GetRasterBand(i+1).SetNoDataValue(noDataValue)
    else:
        ds = driver.Create(outPath, shape[1], shape[0], 1, gdal.GDT_Float32)
        ds.SetProjection(proj)
        ds.SetGeoTransform(geotransform)
        ds.GetRasterBand(1).WriteArray(data)
        ds.GetRasterBand(1).SetNoDataValue(noDataValue)

    if outPath != "MEM":
        check_gdal_success(outPath, cmd='saveImg')

    print('Saved ' + outPath)
    return ds


def saveImgByCopy(outImg, outPath):
    driver = gdal.GetDriverByName( "GTiff" )
    driver.CreateCopy(outPath, outImg , 0 , driverOptionsGTiff)
    check_gdal_success(outPath, cmd='saveImgByCopy')
    print('Saved ' + outPath)


def getTileExtents(inImg, tileSize):
    # Split the image into square tiles of tileSize pixels and returns the extents
    # ([minX, maxY, maxX, minY]) of each tile in the projected units.
    tileExtents = []
    gt = inImg.GetGeoTransform()

    rows = [int(i*tileSize) for i in range(int(inImg.RasterYSize/tileSize)+1)]
    rows.append(inImg.RasterYSize)
    cols = [int(i*tileSize) for i in range(int(inImg.RasterXSize/tileSize)+1)]
    cols.append(inImg.RasterXSize)

    for y in range(1, len(rows)):
        tileExtents.append([])
        for x in range(1, len(cols)):
            ext = getExtent(gt, cols[x]-1, rows[y]-1, cols[x-1], rows[y-1])
            tileExtents[y-1].append([ext[0][0], ext[0][1], ext[2][0], ext[2][1]])

    return tileExtents


def getExtent(gt,endCol, endRow, startCol = 0, startRow = 0):
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
    ext=[]
    xarr=[startCol,endCol]
    yarr=[startRow,endRow]

    for px in xarr:
        for py in yarr:
            x=gt[0]+(px*gt[1])+(py*gt[2])
            y=gt[3]+(px*gt[4])+(py*gt[5])
            ext.append([x,y])
            #print x,y
        yarr.reverse()
    return ext


def reprojectCoords(coords,src_srs,tgt_srs):
    ''' Reproject a list of x,y coordinates.

        @type geom:     C{tuple/list}
        @param geom:    List of [[x,y],...[x,y]] coordinates
        @type src_srs:  C{osr.SpatialReference}
        @param src_srs: OSR SpatialReference object
        @type tgt_srs:  C{osr.SpatialReference}
        @param tgt_srs: OSR SpatialReference object
        @rtype:         C{tuple/list}
        @return:        List of transformed [[x,y],...[x,y]] coordinates
    '''
    trans_coords=[]
    transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords
