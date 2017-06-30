import os
import re
import glob
import logging
from math import cos, radians, pi

from osgeo import gdal
import numpy as np
from tqdm import trange
from gdal_utils.gdal_utils import array_to_gtiff

from atmospheric_correction.sensors import sensor_is
from atmospheric_correction import metadata as metamod
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
        return toa_radiance_S2(img, mtd_dict=mtdfile)


def toa_radiance_WV(img, mtdfile, doDOS, isPan, sensor):
    """Compute TOA radiance for WV"""

    gain, bias = metamod.wv.get_gain_bias_WV(sensor, isPan)
    effectivebw, abscalfactor = metamod.wv.get_effectivebw_abscalfactor_WV(mtdfile)
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
        rawdata = img.GetRasterBand(b + 1).ReadAsArray()
        mask = (rawdata - dosDN[b]) > 0
        mask &= rawdata < 65536

        radiance_band = rawdata - dosDN[b]
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
    gain, bias = metamod.phr1.get_gain_bias_PHR1(mtdfile)

    nbands = img.RasterCount
    ny, nx = img.RasterYSize, img.RasterXSize

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        dosDN = dos.do_dos(img)
    else:
        dosDN = list(np.zeros(nbands))

        # apply the radiometric correction factors to input image
    logger.info("Radiometric correctionPHR1")
    radiance = np.zeros((ny, nx, nbands))
    for band in range(nbands):
        logger.info(band + 1)
        rawdata = np.int_(img.GetRasterBand(band + 1).ReadAsArray())
        mask = (rawdata - dosDN[band]) > 0
        mask &= rawdata != 65536
        radiance[mask, band] = (rawdata - dosDN[band]) / gain[band] + bias[band]

    # Mark the pixels which have all radiances of 0 as invalid
    allzero = np.all((radiance == 0), axis=-1)
    radiance[allzero, :] = np.nan

    return array_to_gtiff(
            radiance, "MEM", img.GetProjection(), img.GetGeoTransform(), banddim=2)


def read_landsat(img, sensor):

    if img.RasterCount > 1:
        # Panchromatic should only be one band but this way the isPan option can also
        # be used to processed L8 images which are stacked in one file.
        rawdata = np.zeros((img.RasterYSize, img.RasterXSize, img.RasterCount))
        for i in range(img.RasterCount):
            rawdata[:, :, i] = img.GetRasterBand(i+1).ReadAsArray()
    else:
        visnirbands = metamod.l78.get_visnirbands(sensor)
        rawdata = np.zeros((img.RasterYSize, img.RasterXSize, len(visnirbands)))

        # Raw Landsat 8/7 data has each band in a separate image.
        # Therefore first open images with all the required band data.
        imgdir = os.path.dirname(img.GetFileList()[0])
        pattern = os.path.join(imgdir, '*.TIF')
        bandfiles = glob.glob(pattern)
        for bf in sorted(bandfiles):
            try:
                bandstr = re.search('_B(\d+)\.TIF$', os.path.basename(bf)).group(1)
                band = int(bandstr)
            except AttributeError:
                continue
            if band in visnirbands:
                logger.info(band)
                bandimg = gdal.Open(os.path.join(imgdir, bf), gdal.GA_ReadOnly)
                rawdata[:, :, band-1] = bandimg.GetRasterBand(1).ReadAsArray()
                rawdata = np.int_(rawdata)
    return rawdata


def toa_radiance_L8(img, mtdfile, doDOS, isPan, sensor):
    mult_factor, add_factor = metamod.l78.get_correction_factors(mtdfile)

    rawdata = read_landsat(img, sensor)
    ny, nx, nbands = rawdata.shape

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        bandimg = array_to_gtiff(
                rawdata, "MEM", img.GetProjection(), img.GetGeoTransform(), banddim=2)
        dosDN = dos(bandimg)
        bandimg = None
    else:
        dosDN = np.zeros(nbands)

    # apply the radiometric correction factors to input image
    logger.info("Radiometric correctionL8")
    radiance = np.zeros(rawdata.shape)
    for band in range(nbands):
        logger.info(band + 1)
        mask = (rawdata[:, :, band] - dosDN[band]) > 0
        radiance[~mask, band] = 0
        radiance[mask, band] = (
                (rawdata[mask, band] - dosDN[band]) * mult_factor[band] +
                add_factor[band])

    # Mark the pixels which have all radiances of 0 as invalid
    allzero = np.all((radiance == 0), axis=-1)
    radiance[allzero, :] = np.nan

    return array_to_gtiff(
            radiance, "MEM",
            img.GetProjection(), img.GetGeoTransform(), banddim=2)


def toa_radiance_S2(img, mtd_dict):
    """Method taken from the bottom of http://s2tbx.telespazio-vega.de/sen2three/html/r2rusage.html

    Parameters
    ----------
    img : gdal image
        input data
    mtd_dict : dict
        meta data

    Metadata contents
    -----------------
    mtd_dict['band_ids'] : sequence, optional
        band IDs (0-based)
        default: 0-9

    Note
    ----
    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    rc = float(mtd_dict['reflection_conversion'])
    u = float(mtd_dict['quantification_value'])
    irradiance = [float(e) for e in mtd_dict['irradiance_values']]
    sun_zenith = float(mtd_dict[mtd_dict['current_granule']]['sun_zenit'])
    band_ids = mtd_dict.get('band_ids', None)
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
