# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 09:30:18 2014

@author: rmgu
"""

# A quick fix for PyInstaller problem with matplotlib fonts (https://github.com/pyinstaller/pyinstaller/issues/885)
# A permament solution was pulled to GitHub on 18th November but should wait untill new release of PyInstaller is out
import sys
import os 
if not sys.getfilesystemencoding(): 
    sys.getfilesystemencoding = lambda: 'UTF-8'

import re
from math import *
# !!
# NOTE: Two line below must be uncommented if compiling with pyinstaller
#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
import numpy as np
import scipy.stats as st
import csv
from scipy.interpolate import interp1d
import scipy.interpolate
import scipy.ndimage


driverOptionsGTiff = ['COMPRESS=DEFLATE', 'PREDICTOR=1', 'BIGTIFF=IF_SAFER']

###############################################################################################
# Utility functions

# Band numbers are in numpy convention (starting from 0). GDAL starts the numbering from 1. 
# Therefore if a band is read using GDAL GetRasterBand the band number should be incremented by 1. 
def getSensorBandNumber(sensor):
    if sensor == "WV2" or sensor == "WV3":
        # mapping of band names to band number for WV2
        coastal = 0; blue = 1; green = 2; yellow = 3; red = 4; redEdge = 5; nir1 = 6; nir2 = 7;
    elif sensor == "GE1":
        # mapping of band names to band number for GE1
        coastal = 0; blue = 0; green = 1; yellow = 1; red = 2; redEdge = 2; nir1 = 3; nir2 = 3;
    elif sensor == "PHR1A" or sensor == "PHR1B" or sensor == "SPOT6":
        # mapping of band names to band number for Pleiades
        coastal = 2; blue = 2; green = 1; yellow = 1; red = 0; redEdge = 0; nir1 = 3; nir2 = 3;
    elif sensor == "L8":
        # mapping of band names to band number for Landsat 8
        coastal = 0; blue = 1; green = 2; yellow = 2; red = 3; redEdge = 3; nir1 = 4; nir2 = 4;
    elif sensor == "L7":
        # mapping of band names to band number for Landsat 7
        coastal = 0; blue = 0; green = 1; yellow = 1; red = 2; redEdge = 2; nir1 = 3; nir2 = 3;
    elif sensor == "S2A_10m":
        # mapping of band names to band number for 10 m Sentinel-2A
        coastal = 0; blue = 0; green = 1; yellow = 1; red = 2; redEdge = 2; nir1 = 7; nir2 = 7;
    elif sensor == "S2A_60m":
        # mapping of band names to band number for 10 m Sentinel-2A
        # assumes that bands with higher spatial resolution have been resampled to 60 m
        coastal = 0; blue = 1; green = 2; yellow = 2; red = 3; redEdge = 4; nir1 = 7; nir2 = 9; 
    else:
        raise Exception("Unknown sensor!")        
        
    return coastal, blue, green, yellow, red, redEdge, nir1, nir2


# read the sun zenith angle and off nadir viewing angle from the WV2 metadata file 
# and return them in radians
def viewingGeometryWV2(metadataFile):
    meanSunElRegex = "\s*meanSunEl\s*=\s*(.*);"
    meanOffNadriView = "\s*meanOffNadirViewAngle\s*=\s*(.*);"
    sunEl = 0.0
    offNadirView = 90.0
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(meanSunElRegex, line)
            if match:
                sunEl = float(match.group(1))
            match = re.match(meanOffNadriView, line)
            if match:
                offNadirView = float(match.group(1))
                
    sunZen = 90.0 - sunEl
    return radians(sunZen), radians(offNadirView)  

def viewingGeometryL8(metadataFile):
    
    # read viewing gemotery from Landsat metadata file
    sunElRegex = "\s*SUN_ELEVATION\s*=\s*(.*)\s*"
    sunEl = 0.0
    satZen = 0.0
        
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(sunElRegex , line)
            if match:
                sunEl = float(match.group(1))
        
    sunZen = 90 - sunEl
    
    return radians(sunZen), radians(satZen)

def viewingGeometryS2A(metadataFile):
    return None, None

def viewingGeometry(metadataFile, sensor):
    if sensor == "WV2" or sensor == "WV3" or sensor == "GE1":
        return viewingGeometryWV2(metadataFile)
    elif sensor == "L8" or sensor == "L7":
        return viewingGeometryL8(metadataFile)
    elif sensor == "S2A_10m" or sensor == "S2A_60m":
        return viewingGeometryS2A(metadataFile)

def mask(img, mask, maskValue = 0):
    maskData = mask.GetRasterBand(1).ReadAsArray()
    maskedData = np.zeros((img.RasterYSize, img.RasterXSize, img.RasterCount))    
    for band in range(1,img.RasterCount+1): 
        maskedData[:,:,band-1] = img.GetRasterBand(band).ReadAsArray()
    maskedData[maskData==maskValue,:] = 0
    res = saveImg (maskedData, img.GetGeoTransform(), img.GetProjection(), "MEM")
    return res    

# Input: Path to a shape file, filename of output raster and pixel size
# Output: GDAL object with rasterised shape file
def rasterize(in_vector, out_raster, pixel_size=25):
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

# Input: two GDALDataset 
# Output: minimum covering extent of the two datasets
# [minX, maxY, maxX, minY] (UL, LR)
def calcMinCoveringExtent(imgA, imgB):
    aGeoTrans = imgA.GetGeoTransform()
    bGeoTrans = imgB.GetGeoTransform()
    minX = max(aGeoTrans[0], bGeoTrans[0])
    maxY = min(aGeoTrans[3], bGeoTrans[3])
    maxX = min(aGeoTrans[0]+imgA.RasterXSize*aGeoTrans[1], bGeoTrans[0]+imgB.RasterXSize*bGeoTrans[1])
    minY = max(aGeoTrans[3]+imgA.RasterYSize*aGeoTrans[5], bGeoTrans[3]+imgB.RasterYSize*bGeoTrans[5])
    return [minX, maxY, maxX, minY]
        
# Input: raster file and shape file
# Output: GDAL object with the clipped raster
def clipRasterWithShape(rasterImg, shapeImg):
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
    rasterSubsetPixs = world2Pixel(rasterGeoTrans, minX, maxY) +  world2Pixel(rasterGeoTrans, maxX, minY)
    shapeRasterSubsetPixs = world2Pixel(shapeRasterGeoTrans, minX, maxY) + world2Pixel(shapeRasterGeoTrans, maxX, minY)     
    

    # clip the shapeRaster to min covering extent
    shapeRasterData = shapeRasterImg.GetRasterBand(1).ReadAsArray()
    shapeRasterClipped = shapeRasterData[shapeRasterSubsetPixs[1]:shapeRasterSubsetPixs[3], shapeRasterSubsetPixs[0]:shapeRasterSubsetPixs[2]]    
    
    # go through the raster bands, clip the to the minimum covering extent and mask out areas not covered by vector
    maskedData = np.zeros((np.shape(shapeRasterClipped)[0], np.shape(shapeRasterClipped)[1], rasterImg.RasterCount))
    bandNum = rasterImg.RasterCount    
    for band in range(1,bandNum+1):
        rasterData = rasterImg.GetRasterBand(band).ReadAsArray()
        clippedData = rasterData[rasterSubsetPixs[1]:rasterSubsetPixs[3], rasterSubsetPixs[0]:rasterSubsetPixs[2]]
        maskedData[:,:,band-1] = np.where(shapeRasterClipped > 0, clippedData, np.NaN)
        
    # get the geotransform array for the masekd array
    maskedGeoTrans = (ulX, rasterGeoTrans[1], rasterGeoTrans[2], ulY, rasterGeoTrans[4], rasterGeoTrans[5])
    
    # save the masked img to memory and return it
    maskedImg = saveImg(maskedData, maskedGeoTrans, rasterImg.GetProjection(), "MEM", np.NaN)    
    return maskedImg


def openAndClipRaster(inFilename, shapeRoiFilename):
    inImg = gdal.Open(inFilename , gdal.GA_ReadOnly)

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
  rtnX = geoMatrix[2]
  rtnY = geoMatrix[4]
  pixel = int(round((x - ulX) / xDist))
  line = int(round((ulY - y) / xDist))
  return [pixel, line]
   
# save the data to geotiff or memory    
def saveImg (data, geotransform, proj, outPath, noDataValue = np.nan):
    
    # Start the gdal driver for GeoTIFF
    if outPath == "MEM":
        driver = gdal.GetDriverByName("MEM")
        driverOpt = []
    else:
        driver = gdal.GetDriverByName("GTiff")
        driverOpt = driverOptionsGTiff
        #driverOpt = []    
    
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
            
    print('Saved ' +outPath )

    return ds

def saveImgByCopy(outImg, outPath):
        
    driver = gdal.GetDriverByName ( "GTiff" )
    savedImg = driver.CreateCopy(outPath, outImg , 0 , driverOptionsGTiff)
    savedImg = None
    outImg = None
    print('Saved ' +outPath )

def plotRasterPointsFromVector(inImg, inVector, fieldName, invertMeasuredPoints, limits = [-30, 0], tide = 0, title = "", qualityImg = None, qualityThresholds = None, figureFilename = ""):

    print("Plotting")
    gt=inImg.GetGeoTransform()
    rb=inImg.GetRasterBand(1).ReadAsArray()
    
    # Set quality class according to quality bands and thresholds
    if qualityImg and qualityThresholds and qualityImg.RasterCount == len(qualityThresholds) == 2:
        # First band contains euclaidian distance between modelled and measure spectra, second band contains the contribution of bottom signal to total modelled reflectance        
        qualityMap = np.zeros((qualityImg.RasterYSize, qualityImg.RasterXSize, 2))        
        qualityMap[:,:,0] = qualityImg.GetRasterBand(1).ReadAsArray()
        qualityMap[:,:,1] = qualityImg.GetRasterBand(2).ReadAsArray()
        qualityClassMap = np.zeros((qualityImg.RasterYSize, qualityImg.RasterXSize))
        # good closure and optically shallow 
        qualityClassMap[np.logical_and(qualityMap[:,:,0] <= qualityThresholds[0], qualityMap[:,:,1] >= qualityThresholds[1])] = 1
        # good closure and optically deep        
        qualityClassMap[np.logical_and(qualityMap[:,:,0] <= qualityThresholds[0], qualityMap[:,:,1] < qualityThresholds[1])] = 2
        # bad closure and optically shallow
        qualityClassMap[np.logical_and(qualityMap[:,:,0] > qualityThresholds[0], qualityMap[:,:,1] >= qualityThresholds[1])] = 3
        # bad closure and optically deep
        qualityClassMap[np.logical_and(qualityMap[:,:,0] > qualityThresholds[0], qualityMap[:,:,1] < qualityThresholds[1])] = 4
    else:
        qualityClassMap = None
        
    modelledValues = []
    measuredValues = []
    qualityValues = []    
    
    lyr=inVector.GetLayer()
    for feat in lyr:
        geom = feat.GetGeometryRef()
        mx,my=geom.GetX(), geom.GetY()  #coord in map units
    
        # Convert from map to pixel coordinates.
        [px, py] = world2Pixel(gt, mx, my)

        # pick the pixels corresponding to known depth points 
        if px >= 0 and px < inImg.RasterXSize and py >= 0 and py < inImg.RasterYSize :       
            pointValue = rb[py,px]+tide
            measuredValue = float(feat.GetField(fieldName))
            if invertMeasuredPoints:
                measuredValue = -measuredValue
            if not np.isnan(pointValue) and not np.isnan(measuredValue) and pointValue < limits[1] and pointValue >= limits[0] and measuredValue < limits[1] and measuredValue >= limits[0]:
                modelledValues.append(pointValue)
                measuredValues.append(measuredValue)
                if qualityClassMap is not None:
                    qualityValues.append(qualityClassMap[py,px])
    
    # calculate statistics
    a = np.array(measuredValues)
    b = np.array(modelledValues)    
    rmse = sqrt(np.mean((a-b)**2))
    bias = np.mean(b) - np.mean(a)
    r, _ = st.pearsonr(a,b)
    print("RMSE: "+str(rmse))
    print("Bias: "+str(bias))    
    print("r: "+str(r))
    print("Number of points: "+str(len(measuredValues))) 
    
    # plot
    plt.figure()
    if qualityValues:
        if len(measuredValues) > 300:
            marker = '.'
            s = 3
        else:
            marker = 'x'
            s = 30
        measuredValues = np.array(measuredValues)
        modelledValues = np.array(modelledValues)
        qualityValues = np.array(qualityValues)
        
        plt.scatter(measuredValues[qualityValues == 2], modelledValues[qualityValues == 2], s = s, marker=marker, c = 'b', edgecolors='none')
        plt.scatter(measuredValues[qualityValues == 3], modelledValues[qualityValues == 3], s = s, marker=marker, c = 'r', edgecolors='none')
        plt.scatter(measuredValues[qualityValues == 4], modelledValues[qualityValues == 4], s = s, marker=marker, c = 'k', edgecolors='none')
        plt.scatter(measuredValues[qualityValues == 1], modelledValues[qualityValues == 1], s = s, marker=marker, c = 'g', edgecolors='none')
    else:
        if len(measuredValues) > 300:    
            plt.hexbin(measuredValues, modelledValues, gridsize=(500,100), cmap=plt.cm.RdBu_r, mincnt=1)
        else:
            plt.scatter(measuredValues, modelledValues, s = 30, marker='x')    
    plt.axis([limits[0], limits[1], limits[0], limits[1]])    
    plt.grid(True)
    plt.xlabel('Measurements')
    plt.ylabel('Bathymetry')
    plt.title(title)
    plt.text(plt.xlim()[0]+0.5, plt.ylim()[1]+(plt.ylim()[0]/20.0)*1.0, "RMSE:   %.2f"% rmse)
    plt.text(plt.xlim()[0]+0.5, plt.ylim()[1]+(plt.ylim()[0]/20.0)*2.0, "Bias:     %.2f"% bias)
    plt.text(plt.xlim()[0]+0.5, plt.ylim()[1]+(plt.ylim()[0]/20.0)*3.0, "r:          %.2f"% r)
    plt.text(plt.xlim()[0]+0.5, plt.ylim()[1]+(plt.ylim()[0]/20.0)*4.0, "Points:          %d"% len(measuredValues)) 
    if figureFilename:
        try:
            plt.savefig(figureFilename)
            plt.close()
        except:
            plt.close()
    else:
        plt.show()
        
# resample the given band filter to specified spectral resolution
def resampleBandFilters(bandFilter, startWV, endWV, resolution):
    x = np.linspace(startWV, endWV, len(bandFilter))
    f = interp1d(x, bandFilter, kind='slinear')

    xnew = np.linspace(startWV, endWV, (endWV-startWV)/resolution + 1)

    return f(xnew)

# Read band filters from CSV file and assign them to a given sensor
def readBandFiltersFromCSV(csvFilename, sensor, isPan):
    # create empty lists for all the possible bands    
    pan = []; coastal = []; blue = []; green = []; yellow = []; red = []; rededge = []; nir1 = []; nir2 = []; wavelength = []    
    # S2 has some extra bands
    rededge2 = []; rededge3 = []#; nir3= []; swir1 = []; swir2 = []; swir3 = []  
    
    # read in the data from CSV file    
    with open(csvFilename, 'r') as csvFile:
        reader = csv.DictReader(csvFile)
        for line in reader:
            if sensor == "WV2" or sensor == "WV3":
                wavelength.append(float(line["Wavelength"]))
                pan.append(float(line["Panchromatic"]))
                coastal.append(float(line["Coastal"]))
                blue.append(float(line["Blue"]))
                green.append(float(line["Green"]))
                yellow.append(float(line["Yellow"]))
                red.append(float(line["Red"]))
                rededge.append(float(line["Red Edge"]))
                nir1.append(float(line["NIR1"]))
                nir2.append(float(line["NIR2"]))
            elif sensor == "PHR1A" or sensor == "PHR1B" or sensor == "SPOT6":
                wavelength.append(float(line["Wavelength"]))
                blue.append(float(line["B1Blue"]))
                green.append(float(line["B2Green"]))
                red.append(float(line["B3Red"]))
                nir1.append(float(line["B4NIR"]))
            elif sensor == "L8":
                wavelength.append(float(line["Wavelength"]))
                coastal.append(float(line["L8B1Coast"]))
                blue.append(float(line["L8B2Blue"]))
                green.append(float(line["L8B3Green"]))
                red.append(float(line["L8B4Red"]))
                nir1.append(float(line["L8B5NIR"]))
                pan.append(float(line["L8B8Pan"]))
            elif sensor == "L7":
                wavelength.append(float(line["Wavelength"]))
                blue.append(float(line["L7B1Blue"]))
                green.append(float(line["L7B2Green"]))
                red.append(float(line["L7B3Red"]))
                nir1.append(float(line["L7B4NIR"]))
            elif sensor == "S2A_10m" or sensor == "S2A_60m":
                wavelength.append(float(line["SR_WL"]))
                coastal.append(float(line["SR_AV_B1"]))
                blue.append(float(line["SR_AV_B2"]))
                green.append(float(line["SR_AV_B3"]))
                red.append(float(line["SR_AV_B4"]))
                rededge.append(float(line["SR_AV_B5"]))
                rededge2.append(float(line["SR_AV_B6"]))
                rededge3.append(float(line["SR_AV_B7"]))
                nir1.append(float(line["SR_AV_B8"]))
                nir2.append(float(line["SR_AV_B8A"]))                                
#                nir3.append(float(line["SR_AV_B9"]))
#                swir1.append(float(line["SR_AV_B10"]))
#                swir2.append(float(line["SR_AV_B11"]))
#                swir3.append(float(line["SR_AV_B12"]))
    
    # collect bands specific for each sensor and start and end wavelenghts            
    startWV = wavelength[0]
    endWV = wavelength[-1]
    if sensor == "WV2" or sensor == "WV3":
        if not isPan:        
            bandFilters = [coastal, blue, green, yellow, red, rededge, nir1, nir2]
        else:
            bandFilters = [pan]
    elif sensor == "PHR1A" or sensor == "PHR1B" or sensor == "SPOT6":
        # the order of Pleadis bands is like below (RGBN), not like indicated in the metadata file (BGRN)
        bandFilters = [red, green, blue, nir1]
    elif sensor == "L8":
        if not isPan:
            bandFilters = [coastal, blue, green, red, nir1]
        else:
            bandFilters = [pan]
    elif sensor == "L7":
        if not isPan:
            bandFilters = [blue, green, red, nir1]
        else:
            bandFilters = [pan]
    elif sensor == "S2A_10m" or sensor == "S2A_60m":
        bandFilters = [coastal, blue, green, red, rededge, rededge2, rededge3, nir1, nir2]#, nir3, swir1, swir2, swir3]
            
    return startWV, endWV, bandFilters

# Split the image into square tiles of tileSize pixels and returns the extents
# ([minX, maxY, maxX, minY]) of each tile in the projected units.
def getTileExtents(inImg, tileSize):
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
    
# Return list of corner coordinates from a geotransform
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
    
# Reproject a list of x,y coordinates.
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
    
