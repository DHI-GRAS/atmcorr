import logging

import numpy as np

from atmcorr import dos
from atmcorr.pleiades import calibration

logger = logging.getLogger(__name__)


def dn_to_radiance(dndata, mtdFile, band_ids, doDOS=False):
    """Apply radiometric correction to Pleadis image,
       with DOS atmospheric correction optional.
    """
    gain, bias = calibration.get_gain_bias_PHR1(mtdFile)

    gain = [gain[i] for i in band_ids]
    bias = [bias[i] for i in band_ids]

    nbands = dndata.shape[0]

    # perform dark object substraction
    if doDOS:
        logger.info("DOS correction")
        dosDN = dos.do_dos(dndata)
    else:
        dosDN = np.zeros(nbands)

        # apply the radiometric correction factors to input image
    logger.info("Radiometric correction PHR1")
    radiance = np.zeros(dndata.shape)
    for i in range(nbands):
        logger.info(i + 1)
        rawdata = dndata[i]
        mask = (rawdata - dosDN[i]) > 0
        mask &= rawdata != 65536
        radiance[mask, i] = (rawdata - dosDN[i]) / gain[i] + bias[i]

    # Mark the pixels which have all radiances of 0 as invalid
    allzero = np.all((radiance == 0), axis=-1)
    radiance[allzero, :] = np.nan

    logger.info('Done with radiometric correction.')
    return radiance
