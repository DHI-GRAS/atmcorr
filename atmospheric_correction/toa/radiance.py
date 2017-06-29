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
from atmospheric_correction import dos

logger = logging.getLogger(__name__)


def toa_radiance(img, mtdfile, sensor, doDOS, isPan=False):
    if sensor_is(sensor, 'WV'):
        return toa_radiance_WV(img, mtdfile, doDOS, isPan, sensor)
    elif sensor_is(sensor, 'PHR'):
        return toa_radiance_PHR1(img, mtdfile, doDOS)
    elif sensor_is(sensor, 'L7L8'):
        return toa_radiance_L8(img, mtdfile, doDOS, isPan, sensor)
    elif sensor_is(sensor, 'S2'):
        return toa_radiance_S2(img, mtdfile)


def get_effectivebw_abscalfactor_WV(mtdfile):
    # get the correction factors from the metadata file, assuming the number and
    # order of bands is the same in the image and the metadata file
    abscalfactor_rgx = "\s*absCalFactor\s*=\s*(.*);"
    effectivebw_rgx = "\s*effectiveBandwidth\s*=\s*(.*);"

    effectivebw = []
    abscalfactor = []
    with open(mtdfile) as mf:
        for line in mf:
            match = re.match(abscalfactor_rgx, line)
            if match:
                abscalfactor.append(float(match.group(1)))
            match = re.match(effectivebw_rgx, line)
            if match:
                effectivebw.append(float(match.group(1)))
    if not effectivebw:
        raise ValueError('Unable to get effective bandwidth from mtdfile.')
    if not abscalfactor:
        raise ValueError('Unable to get abs cal factor from mtdfile.')
    return np.array(effectivebw), np.array(abscalfactor)


def get_gain_bias_WV(sensor, isPan):
    """Get WV calibration factors"""
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
    return np.array(gain), np.array(bias)


def toa_radiance_WV(img, mtdfile, doDOS, isPan, sensor):
    """Compute TOA radiance for WV"""

    gain, bias = get_gain_bias_WV(sensor, isPan)
    effectivebw, abscalfactor = get_effectivebw_abscalfactor_WV(mtdfile)
    scalefactor = abscalfactor / effectivebw * (2 - gain) - bias

    nbands = img.RasterCount

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        dosDN = dos.do_dos(img)
    else:
        dosDN = np.zeros(nbands)

    # apply the radiometric correction factors to input image
    radiance = np.zeros((img.RasterYSize, img.RasterXSize, img.RasterCount))
    for b in trange(nbands, desc='Radiometric correction WV', unit='b'):
        rawData = img.GetRasterBand(b + 1).ReadAsArray()
        mask = (rawData - dosDN[b]) > 0
        mask &= rawData < 65536

        radiance_band = rawData - dosDN[b]
        radiance[mask, b] = radiance_band[mask]
        radiance[:, :, b] *= scalefactor[b]

    # Mark the pixels which have all radiances of 0 as invalid
    mask = np.sum(radiance, axis=2) > 0
    radiance[~mask, :] = 0

    return array_to_gtiff(
            radiance, "MEM", img.GetProjection(), img.GetGeoTransform(), banddim=2)


def toa_radiance_PHR1(img, mtdfile, doDOS=False):
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
    tree = ET.parse(mtdfile)
    root = tree.getroot()
    Radiometric_Data = root.findall('Radiometric_Data')[0]
    Radiometric_Calibration = Radiometric_Data.findall('Radiometric_Calibration')[0]
    Instrument_Calibration = Radiometric_Calibration.findall('Instrument_Calibration')[0]
    Band_Measurement_List = Instrument_Calibration.findall('Band_Measurement_List')[0]
    for Band_Radiance in Band_Measurement_List.findall('Band_Radiance'):
        band = Band_Radiance.findall('BAND_ID')[0].text
        gain[bandMapping[band]] = float(Band_Radiance.findall('GAIN')[0].text)
        bias[bandMapping[band]] = float(Band_Radiance.findall('BIAS')[0].text)

    nbands = img.RasterCount

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        dosDN = dos.do_dos(img)
    else:
        dosDN = list(np.zeros(nbands))

        # apply the radiometric correction factors to input image
    logger.info("Radiometric correctionPHR1")
    radiance = np.zeros((img.RasterYSize, img.RasterXSize, img.RasterCount))
    validMask = np.zeros((img.RasterYSize, img.RasterXSize))
    for band in range(nbands):
        logger.info(band + 1)
        rawData = np.int_(img.GetRasterBand(band + 1).ReadAsArray())
        radiance[:, :, band] = np.where(np.logical_and((rawData - dosDN[band]) > 0, rawData != 65536),
                                               (rawData - dosDN[band]) / gain[band] + bias[band], 0)
        validMask += radiance[:, :, band]

    # Mark the pixels which have all radiances of 0 as invalid
    invalidMask = np.where(validMask > 0, False, True)
    radiance[invalidMask, :] = np.nan

    return array_to_gtiff(radiance, "MEM", img.GetProjection(), img.GetGeoTransform(), banddim=2)


def toa_radiance_L8(img, mtdfile, doDOS, isPan, sensor):
    multFactorRegex = "\s*RADIANCE_MULT_BAND_\d\s*=\s*(.*)\s*"
    addFactorRegex = "\s*RADIANCE_ADD_BAND_\d\s*=\s*(.*)\s*"

    if img.RasterCount == 1:
        if sensor == "L8":
            # The first 5 bands in L8 are VIS/NIR
            visNirBands = range(1,6)
        elif sensor == "L7":
            # The first 4 bands in L7 are VIS/NIR
            visNirBands = range(1,5)
        rawData = np.zeros((img.RasterYSize, img.RasterXSize, len(visNirBands)))

        # Raw Landsat 8/7 data has each band in a separate image. Therefore first open images with all the required band data.
        imgDir = os.path.dirname(img.GetFileList()[0])
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
        rawData = np.zeros((img.RasterYSize, img.RasterXSize, img.RasterCount))
        for i in range(img.RasterCount):
            rawData[:, :, i] = img.GetRasterBand(i+1).ReadAsArray()

    # get the correction factors from the metadata file, assuming the number and
    # order of bands is the same in the image and the metadata file
    multFactor = []
    addFactor = []
    with open(mtdfile, 'r') as metadata:
        for line in metadata:
            match = re.match(multFactorRegex, line)
            if match:
                multFactor.append(float(match.group(1)))
            match = re.match(addFactorRegex, line)
            if match:
                addFactor.append(float(match.group(1)))

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        rawImg = array_to_gtiff(rawData, "MEM", img.GetProjection(), img.GetGeoTransform(), banddim=2)
        dosDN = dos(rawImg)
        rawImg = None
    else:
        dosDN = list(np.zeros(rawData.shape[2]))

    # apply the radiometric correction factors to input image
    logger.info("Radiometric correctionL8")
    radiance = np.zeros((img.RasterYSize, img.RasterXSize, rawData.shape[2]))
    validMask = np.zeros((img.RasterYSize, img.RasterXSize))
    for band in range(rawData.shape[2]):
        logger.info(band + 1)
        mask = (rawData[:, :, band] - dosDN[band]) > 0
        radiance[~mask, band] = 0
        radiance[mask, band] = (
                (rawData[mask, band] - dosDN[band]) *
                multFactor[band] + addFactor[band])
        validMask += radiance[:, :, band]

    # Mark the pixels which have all radiances of 0 as invalid
    invalidMask = validMask <= 0
    radiance[invalidMask, :] = np.nan

    return array_to_gtiff(
            radiance, "MEM",
            img.GetProjection(), img.GetGeoTransform(), banddim=2)


def toa_radiance_S2(img, metadata):
    """Method taken from the bottom of http://s2tbx.telespazio-vega.de/sen2three/html/r2rusage.html

    Parameters
    ----------
    img : gdal image
        input data
    metadata : dict
        meta data

    Metadata contents
    -----------------
    metadata['band_ids'] : sequence, optional
        band IDs (0-based)
        default: 0-9

    Note
    ----
    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    rc = float(metadata['reflection_conversion'])
    u = float(metadata['quantification_value'])
    irradiance = [float(e) for e in metadata['irradiance_values']]
    sun_zenith = float(metadata[metadata['current_granule']]['sun_zenit'])
    band_ids = metadata.get('band_ids', None)
    if band_ids is None:
        logger.info('Assuming 9 bands')
        band_ids = list(range(9))

    # Convert to radiance
    logger.info("Radiometric correction")
    radiance = np.zeros((img.RasterYSize, img.RasterXSize, len(band_ids)))
    for i, band_id in enumerate(band_ids):
        rToa = img.GetRasterBand(i + 1).ReadAsArray().astype(float)
        rToa /= rc
        factor = irradiance[band_id] * cos(radians(sun_zenith)) / (pi * u)
        radiance[:, :, i] = rToa * factor

    return array_to_gtiff(
            radiance, "MEM",
            img.GetProjection(), img.GetGeoTransform(), banddim=2)
