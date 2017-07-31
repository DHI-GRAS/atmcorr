import os
import re
import logging
from math import cos, radians, pi
from xml.etree import ElementTree as ET

from osgeo import gdal
import numpy as np
from tqdm import trange
from gdal_utils.gdal_utils import array_to_gtiff

from atmospheric_correction.sensors import sensor_is

logger = logging.getLogger(__name__)


def toaRadiance(inImg, metadataFile, sensor, doDOS, isPan=False):
    if sensor_is(sensor, 'WV'):
        return toaRadianceWV(inImg, metadataFile, doDOS, isPan, sensor)
    elif sensor_is(sensor, 'PHR'):
        return toaRadiancePHR1(inImg, metadataFile, doDOS)
    elif sensor_is(sensor, 'L7L8'):
        return toaRadianceL8(inImg, metadataFile, doDOS, isPan, sensor)
    elif sensor_is(sensor, 'S2'):
        return toaRadianceS2(inImg, metadataFile)


def toaRadianceWV(inImg, metadataFile, doDOS, isPan, sensor):
    absCalFactorRegex = "\s*absCalFactor\s*=\s*(.*);"
    effectiveBandwidthRegex = "\s*effectiveBandwidth\s*=\s*(.*);"

    # get the correction factors from the metadata file, assuming the number and
    # order of bands is the same in the image and the metadata file
    effectiveBandwidth, absCalFactor = [], []
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(absCalFactorRegex, line)
            if match:
                absCalFactor.append(float(match.group(1)))
            match = re.match(effectiveBandwidthRegex, line)
            if match:
                effectiveBandwidth.append(float(match.group(1)))

    bandNum = inImg.RasterCount
    if sensor == "WV3":
        # WV3 calibration factors
        if isPan:
            gain = [1.045]
            bias = [2.22]
        else:
            gain = [1.157, 1.07, 1.082, 1.048, 1.028, 0.979, 1.006,
                    0.975]  # as per calibration from DG released 3/6/2015
            bias = [7.07, 4.253, 2.633, 2.074, 1.807, 2.633, 3.406,
                    2.258]  # as per calibration from DG relased 3/6/2015
    else:
        # WV2 calibration factors
        if isPan:
            gain = [1.0264]
            bias = [6.5783]
        else:
            gain = [0.8632, 1.001, 1.0436, 1.0305, 1.0249, 0.9779, 0.981, 0.9217]
            bias = [6.6863, 2.399, 0.3973, 0.7744, -0.1495, 2.0383, 1.859, 2.0357]

    # perform dark object substraction
    dosDN = bandNum * [0]
    if doDOS:
        dosDN = darkObjectSubstraction(inImg)

    # apply the radiometric correction factors to input image
    radiometricData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))
    for band in trange(bandNum, desc='Radiometric correction WV', unit='band'):
        rawData = inImg.GetRasterBand(band + 1).ReadAsArray()
        radiometricData[:, :, band] = np.where(np.logical_and((rawData - dosDN[band]) > 0, rawData < 65536),
                                               rawData - dosDN[band], 0) * absCalFactor[band] / effectiveBandwidth[
                                          band] * (2 - gain[band]) - bias[band]

    validMask = np.sum(radiometricData, axis=2)
    # Mark the pixels which have all radiances of 0 as invalid
    invalidMask = np.where(validMask > 0, False, True)
    radiometricData[invalidMask, :] = 0

    return array_to_gtiff(radiometricData, "MEM", inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)


def toaRadiancePHR1(inImg, metadataFile, doDOS=False):
    """Apply radiometric correction to Pleadis image,
       with DOS atmospheric correction optional.
    """
    # get correction factors
    gain = [0, 0, 0, 0]
    bias = [0, 0, 0, 0]

    # In the XML file the band order is specified as BGRN.
    # However in reality it is RGBN. Therefore mapping is required
    # between the two
    bandMapping = {'B2': 0, 'B1': 1, 'B0': 2, 'B3': 3}

    # get down to the appropirate node
    tree = ET.parse(metadataFile)
    root = tree.getroot()
    Radiometric_Data = root.findall('Radiometric_Data')[0]
    Radiometric_Calibration = Radiometric_Data.findall('Radiometric_Calibration')[0]
    Instrument_Calibration = Radiometric_Calibration.findall('Instrument_Calibration')[0]
    Band_Measurement_List = Instrument_Calibration.findall('Band_Measurement_List')[0]
    for Band_Radiance in Band_Measurement_List.findall('Band_Radiance'):
        band = Band_Radiance.findall('BAND_ID')[0].text
        gain[bandMapping[band]] = float(Band_Radiance.findall('GAIN')[0].text)
        bias[bandMapping[band]] = float(Band_Radiance.findall('BIAS')[0].text)

    bandNum = inImg.RasterCount

    # perform dark object substraction
    if doDOS:
        dosDN = darkObjectSubstraction(inImg)
    else:
        dosDN = list(np.zeros(bandNum))

        # apply the radiometric correction factors to input image
    logger.info("Radiometric correctionPHR1")
    radiometricData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))
    validMask = np.zeros((inImg.RasterYSize, inImg.RasterXSize))
    for band in range(bandNum):
        logger.info(band + 1)
        rawData = np.int_(inImg.GetRasterBand(band + 1).ReadAsArray())
        radiometricData[:, :, band] = np.where(np.logical_and((rawData - dosDN[band]) > 0, rawData != 65536),
                                               (rawData - dosDN[band]) / gain[band] + bias[band], 0)
        validMask += radiometricData[:, :, band]

    # Mark the pixels which have all radiances of 0 as invalid
    invalidMask = np.where(validMask > 0, False, True)
    radiometricData[invalidMask, :] = np.nan

    return array_to_gtiff(radiometricData, "MEM", inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)


def toaRadianceL8(inImg, metadataFile, doDOS, isPan, sensor):
    multFactorRegex = "\s*RADIANCE_MULT_BAND_\d\s*=\s*(.*)\s*"
    addFactorRegex = "\s*RADIANCE_ADD_BAND_\d\s*=\s*(.*)\s*"

    if inImg.RasterCount == 1:
        if sensor == "L8":
            # The first 5 bands in L8 are VIS/NIR
            visNirBands = range(1,6)
        elif sensor == "L7":
            # The first 4 bands in L7 are VIS/NIR
            visNirBands = range(1,5)
        rawData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, len(visNirBands)))

        # Raw Landsat 8/7 data has each band in a separate image. Therefore first open images with all the required band data.
        imgDir = os.path.dirname(inImg.GetFileList()[0])
        for _, _, files in os.walk(imgDir):
            for name in sorted(files):
                match = re.search('(([A-Z]{2}\d).+)_B(\d+)\.TIF$', name)
                if match and int(match.group(3)) in visNirBands:
                    band = int(match.group(3))
                    logger.info(band)
                    rawImg = gdal.Open(os.path.join(imgDir, name), gdal.GA_ReadOnly)
                    rawData[:, :, band-1] = rawImg.GetRasterBand(1).ReadAsArray()
                    rawData = np.int_(rawData)
                    rawImg = None

    # Panchromatic should only be one band but this way the isPan option can also
    # be used to processed L8 images which are stacked in one file.
    else:
        rawData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))
        for i in range(inImg.RasterCount):
            rawData[:, :, i] = inImg.GetRasterBand(i+1).ReadAsArray()

    # get the correction factors from the metadata file, assuming the number and
    # order of bands is the same in the image and the metadata file
    multFactor = []
    addFactor = []
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(multFactorRegex, line)
            if match:
                multFactor.append(float(match.group(1)))
            match = re.match(addFactorRegex, line)
            if match:
                addFactor.append(float(match.group(1)))

    # perform dark object substraction
    if doDOS:
        rawImg = array_to_gtiff(rawData, "MEM", inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)
        dosDN  = darkObjectSubstraction(rawImg)
        rawImg = None
    else:
        dosDN = list(np.zeros(rawData.shape[2]))

    # apply the radiometric correction factors to input image
    logger.info("Radiometric correctionL8")
    radiometricData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, rawData.shape[2]))
    validMask       = np.zeros((inImg.RasterYSize, inImg.RasterXSize))
    for band in range(rawData.shape[2]):
        logger.info(band + 1)
        radiometricData[:, :, band] = np.where((rawData[:, :, band] - dosDN[band]) > 0,
                                               (rawData[:, :, band] - dosDN[band]) * multFactor[band] + addFactor[band],
                                               0)
        validMask += radiometricData[:, :, band]

    # Mark the pixels which have all radiances of 0 as invalid
    invalidMask = np.where(validMask > 0, False, True)
    radiometricData[invalidMask, :] = np.nan

    return array_to_gtiff(radiometricData, "MEM", inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)


def toaRadianceS2(inImg, metadataFile):
    """Method taken from the bottom of http://s2tbx.telespazio-vega.de/sen2three/html/r2rusage.html

    Note
    ----
    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    rc = float(metadataFile['reflection_conversion'])
    u = float(metadataFile['quantification_value'])
    e0 = [float(e) for e in metadataFile['irradiance_values']]
    z = float(metadataFile[metadataFile['current_granule']]['sun_zenit'])

    # Convert to radiance
    logger.info("Radiometric correction")
    radiometricData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, 9))
    for i in range(9):
        rToa = (inImg.GetRasterBand(i + 1).ReadAsArray().astype(float)) / rc
        radiometricData[:, :, i] = (rToa * e0[i] * cos(radians(z))) / (pi * u)
    return array_to_gtiff(
            radiometricData, "MEM",
            inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)


def darkObjectSubstraction(inImg):
    logger.info("DOS correction")
    dosDN = []
    tempData = inImg.GetRasterBand(1).ReadAsArray()
    numElements = np.size(tempData[tempData != 0])
    for band in range(1,inImg.RasterCount+1):
        # use histogram with 2048 bins since WV2 has 11 bit radiometric resolution
        hist, edges = np.histogram(
                inImg.GetRasterBand(band).ReadAsArray(),
                bins=2048, range=(1, 2048), density=False)
        for i in range(1, len(hist)):
            if hist[i] - hist[i-1] > (numElements-numElements*0.999999):
                dosDN.append(i-1)
                break
    return dosDN
