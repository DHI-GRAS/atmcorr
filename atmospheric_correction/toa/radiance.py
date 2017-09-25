import logging

import numpy as np

from atmospheric_correction.sensors import sensor_is
from atmospheric_correction import metadata as metamod
from atmospheric_correction import dos

logger = logging.getLogger(__name__)


def toa_radiance(
        data, sensor, mtdFile, band_ids,
        doDOS=False, mtdFile_tile=None):
    """Compute TOA radiance

    Parameters
    ----------
    data : ndarray shape(nbands, ny, nx)
        input data
    sensor : str
        sensor name
    mtdFile : str
        path to metadata file
    band_ids : list of int
        band IDs of original product contained in array
        0-based
    doDOS : bool
        do a dark object subtraction
    mtdFile_tile : str
        tile metadata file
        required for Sentinel 2
    """
    commonkw = dict(
            data=data,
            mtdFile=mtdFile,
            doDOS=doDOS,
            band_ids=band_ids)
    if sensor_is(sensor, 'WV'):
        return toa_radiance_WV(sensor=sensor, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        return toa_radiance_PHR1(**commonkw)
    elif sensor_is(sensor, 'L7L8'):
        return toa_radiance_L8(sensor=sensor, **commonkw)
    elif sensor_is(sensor, 'S2'):
        commonkw.pop('doDOS')
        return toa_radiance_S2(mtdFile_tile=mtdFile_tile, **commonkw)


def toa_radiance_WV(data, mtdFile, sensor, band_ids, doDOS=False):
    """Compute TOA radiance for WV"""

    gain = metamod.wv.GAIN[sensor]
    bias = metamod.wv.BIAS[sensor]
    effectivebw, abscalfactor = metamod.wv.get_effectivebw_abscalfactor_WV(mtdFile)
    scalefactor = abscalfactor / effectivebw * (2 - gain)
    bias_bands = bias[band_ids]
    scalefactor_bands = scalefactor[band_ids]

    nbands = data.shape[0]

    # perform dark object substraction
    if doDOS:
        logger.info('DOS correction')
        dosDN = dos.do_dos(data)
        logger.info('Done.')
    else:
        dosDN = np.zeros(nbands)

    # apply the radiometric correction factors to input image
    logger.info('Radiometric correction')
    radiance = np.zeros(data.shape, dtype='f4')
    for i in range(nbands):
        rawdata = data[i]
        good = rawdata > dosDN[i]
        good &= rawdata < 65536
        gooddata = (rawdata[good] - dosDN[i]) * scalefactor_bands[i] - bias_bands[i]
        gooddata[gooddata < 0] = 0
        radiance[i, good] = gooddata

    logger.info('Done with radiometric correction.')
    return radiance


def toa_radiance_PHR1(data, mtdFile, band_ids, doDOS=False):
    """Apply radiometric correction to Pleadis image,
       with DOS atmospheric correction optional.
    """
    gain, bias = metamod.phr1.get_gain_bias_PHR1(mtdFile)

    gain = [gain[i] for i in band_ids]
    bias = [bias[i] for i in band_ids]

    nbands = data.shape[0]

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        dosDN = dos.do_dos(data)
    else:
        dosDN = np.zeros(nbands)

        # apply the radiometric correction factors to input image
    logger.info("Radiometric correction PHR1")
    radiance = np.zeros(data.shape)
    for i in range(nbands):
        logger.info(i + 1)
        rawdata = data[i]
        mask = (rawdata - dosDN[i]) > 0
        mask &= rawdata != 65536
        radiance[mask, i] = (rawdata - dosDN[i]) / gain[i] + bias[i]

    # Mark the pixels which have all radiances of 0 as invalid
    allzero = np.all((radiance == 0), axis=-1)
    radiance[allzero, :] = np.nan

    logger.info('Done with radiometric correction.')
    return radiance


def toa_radiance_L8(data, mtdFile, sensor, band_ids, doDOS=False):
    to_multiply, to_add = metamod.l78.get_correction_factors(mtdFile)

    # subset to band IDs
    to_multiply = [to_multiply[i] for i in band_ids]
    to_add = [to_add[i] for i in band_ids]

    radiance = data.copy()

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        dosDN = dos.do_dos(data)
        radiance -= dosDN[:, None, None]

    # apply the radiometric correction factors to input image
    logger.info("Radiometric correction L8")
    radiance[radiance < 0] = 0
    radiance *= to_multiply[:, None, None]
    radiance += to_add[:, None, None]

    # Mark the pixels which have all radiances of 0 as invalid
    allzero = np.all((radiance == 0), axis=0)
    radiance[allzero, :] = np.nan

    logger.info('Done with radiometric correction.')
    return radiance


def toa_radiance_S2(data, mtdFile, mtdFile_tile, band_ids):
    """Method taken from the bottom of http://s2tbx.telespazio-vega.de/sen2three/html/r2rusage.html

    Parameters
    ----------
    data : ndarray shape(nbands, ny, nx)
        input data
    mtdFile : str
        path to metadata file
    mtdFile_tile : str
        path to granule metadata file
    band_ids : sequence of int
        band IDs (0-based index of bands in img)

    Note
    ----
    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    if mtdFile_tile is None:
        raise ValueError('Tile metadata file required!')

    metadata = metamod.s2.parse_mtdfile(mtdFile, mtdFile_tile=mtdFile_tile)
    tile = list(metadata['granules'])[0]
    logger.debug('Tile is \'%s\'.', tile)
    rc = metadata['reflectance_conversion']
    qv = metadata['quantification_value']
    irradiance = np.array(metadata['irradiance_values'])
    sun_zenith = metadata['granules'][tile]['sun_zenith']

    irradiance = irradiance[band_ids]

    # Convert to radiance
    logger.info('Radiometric correction')
    factor = irradiance * np.cos(np.radians(sun_zenith)) / (np.pi * qv)
    radiance = data.astype('f4')
    for i in range(data.shape[0]):
        radiance[i] /= rc
        radiance[i] *= factor[i]
    logger.info('Done with radiometric correction.')
    return radiance
