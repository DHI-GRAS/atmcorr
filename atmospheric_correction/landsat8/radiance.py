import logging

import numpy as np

from atmospheric_correction import dos
from atmospheric_correction.landsat8 import calibration

logger = logging.getLogger(__name__)


def toa_radiance(data, mtdFile, sensor, band_ids, doDOS=False):
    to_multiply, to_add = calibration.get_correction_factors(mtdFile)

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
