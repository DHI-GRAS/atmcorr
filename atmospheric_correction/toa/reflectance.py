import logging

import numpy as np

from atmospheric_correction.metadata import wv as wvmeta
from atmospheric_correction.sensors import sensor_is

logger = logging.getLogger(__name__)


def toa_reflectance(data, mtdfile, sensor, band_ids):
    commonkw = dict(data=data, mtdfile=mtdfile)
    if sensor_is(sensor, 'WV'):
        res = toa_reflectance_WV(band_ids=band_ids, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        res = toa_reflectance_PHR1(**commonkw)
    elif sensor_is(sensor, 'L7L8'):
        res = toa_reflectance_L8(**commonkw)
    elif sensor_is(sensor, 'S2'):
        res = toa_reflectance_S2(**commonkw)
    return res


def toa_reflectance_WV(data, mtdfile, band_ids):
    """Estimate toa reflectance of radiometric WV data ignoric atmospheric, topographic and
       BRDF effects.

    Notes
    -------
    Based on http://www.digitalglobe.com/sites/default/files/
    Radiometric_Use_of_WorldView-2_Imagery%20%281%29.pdf

    Also works with GeoEye-1 and might work with other Digital Globe
    providers after a small modification
    """
    des, ssi, sza = wvmeta.get_earth_sun_distance(mtdfile)

    # apply the radiometric correction factors to input image
    logger.info("TOA reflectance")
    reflectance = np.zeros(data.shape)
    nbands = data.shape[0]
    for i in range(nbands):
        reflectance[i] = (data[i] * des**2 * np.pi) / (ssi[i] * np.cos(sza))

    return reflectance


def toa_reflectance_PHR1(data, mtdfile):
    # TODO implement
    raise NotImplementedError()


def toa_reflectance_L8(data, mtdfile):
    # TODO Implement
    raise NotImplementedError()


def toa_reflectance_S2(data, mtdfile):
    """Convert to TOA reflectance

    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    rc = float(mtdfile['reflection_conversion'])
    logger.info("TOA reflectance")
    reflectance = data.astype('f4')
    reflectance /= rc
    return reflectance
