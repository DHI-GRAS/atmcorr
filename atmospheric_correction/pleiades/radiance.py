import logging

import numpy as np

from atmospheric_correction import dos
from atmospheric_correction.pleiades import calibration

logger = logging.getLogger(__name__)


def toa_radiance(data, mtdFile, band_ids, doDOS=False):
    """Apply radiometric correction to Pleadis image,
       with DOS atmospheric correction optional.
    """
    gain, bias = calibration.get_gain_bias_PHR1(mtdFile)

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
