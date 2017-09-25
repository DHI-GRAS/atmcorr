import logging

import numpy as np

logger = logging.getLogger(__name__)


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
