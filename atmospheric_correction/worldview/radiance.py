import logging

import numpy as np

from atmospheric_correction import dos

logger = logging.getLogger(__name__)


def toa_radiance_WV(data, mtdFile, sensor, band_ids, doDOS=False):
    """Compute TOA radiance for WV"""

    gain = metamod.wv.GAIN[sensor][:-1]
    bias = metamod.wv.BIAS[sensor][:-1]
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
