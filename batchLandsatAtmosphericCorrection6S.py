import os
import fnmatch
import re
import glob

from osgeo import ogr, gdal

from atmParametersMODIS import getAtmParametersMODIS
from atmCorr6S import atmCorr6S

from bathymetry.preprocessing import taoRadiance


def batchLandsatAtmosphericCorrection6S(landsatDataDir, modisAtmDir, roiDir, outDir):
    """
    landsat data folder must have following structure:
    landsatDataFolder
          landsatImageFolder
              ImageFiles

    modis atm folder must have MOD04, MOD05 and MOD07 files

    Roi folder must contain shape files for each landsat tile with the tile
    number present in the filename. The roi must be in same projection as
    the Landsat scene.
    """
    roiFiles = glob.glob(os.path.join(roiDir, "*.shp"))
    extent = None

    for root, dirnames, filenames in os.walk(landsatDataDir):

        # For each landsat scene find the roi file for the same tile and the roi's extent
        for filename in fnmatch.filter(filenames, '*_B1.tif'):
            roi = None
            filename = filename.replace('TIF','tif')
            match = re.search("L[A-Z]\d(\d{6})(\d{7}).*\.tif", filename)
            if match:
                print ("=========== "+filename+"===========" )
                tile = match.group(1)
                yearDoy = match.group(2)
                for roiFile in roiFiles:
                    if tile in roiFile:
                       roi = ogr.Open(roiFile)
                       extent = roi.GetLayer().GetExtent()
                       filename = os.path.join(root, filename)
                       break
                if not extent:
                    continue

                # Find the overpass time
                overpassTime = -1
                metadataFile = filename.replace("B1.tif", "MTL.txt")
                with open(metadataFile, 'r') as metadata:
                    for line in metadata:
                        match = re.search("SCENE_CENTER_TIME\s*=\s*(\d{2}):(\d{2}):.*Z", line)
                        if match:
                            overpassTime = float(match.group(1)) + float(match.group(2))/60.0
                            break
                if overpassTime < 0:
                    continue

                # With the extent and landsat filename call atmParametersMODIS to get
                # the atmospheric parameters for 6S
                aot, wv, ozone = getAtmParametersMODIS(filename, extent, modisAtmDir, yearDoy, overpassTime, roi)

                # Finally perform the atmospheric correction
                inImg = gdal.Open(filename)
                radianceImg = taoRadiance(inImg, metadataFile, sensor = "L8", doDOS = False)
                reflectanceImg = atmCorr6S(metadataFile, radianceImg, atm = {'AOT':aot, 'PWV':wv, 'ozone':ozone}, sensor = "L8", isPan = False, adjCorr = False)

                # Save the results and clean up
                outFile = os.path.join(outDir, os.path.basename(filename).replace("B1.tif", "6S.tif"))
                driver = gdal.GetDriverByName ( "GTiff" )
                savedImg = driver.CreateCopy(outFile, reflectanceImg , 0 )
                savedImg = None
                inImg = None
                radianceImg = None
                reflectanceImg = None
                #saveImg(outImg, inImg.GetGeoTransform(), inImg.GetProjection(), outFile)
                with open(outFile.replace("6S.tif", "6S_params.txt"), 'w') as paramFile:
                    paramFile.write(str({'AOT':aot, 'PWV':wv, 'ozone':ozone}))


if __name__ == "__main__":
    landsatDataDir = r'Q:\Landsat_Data\tiles'
    modisAtmDir = r'Q:\MODIS_ATM'
    roiDir = r'Q:\Shapefiles\6S_polygons_tiles'
    outDir = r'V:\Landsat_6S_temp'
    batchLandsatAtmosphericCorrection6S(landsatDataDir, modisAtmDir, roiDir, outDir)
