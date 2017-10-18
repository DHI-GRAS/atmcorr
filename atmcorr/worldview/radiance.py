import numpy as np

from atmcorr import dos
from atmcorr.worldview import calibration


def toa_radiance(data, mtdFile, sensor, band_ids, doDOS=False):
    """Compute TOA radiance for WV"""

    gain = calibration.GAIN[sensor][:-1]
    bias = calibration.BIAS[sensor][:-1]
    effectivebw, abscalfactor = calibration.get_effectivebw_abscalfactor_WV(mtdFile)
    scalefactor = abscalfactor / effectivebw * (2 - gain)
    bias_bands = bias[band_ids]
    scalefactor_bands = scalefactor[band_ids]

    nbands = data.shape[0]

    # perform dark object substraction
    if doDOS:
        dosDN = dos.do_dos(data)
    else:
        dosDN = np.zeros(nbands)

    # apply the radiometric correction factors to input image
    radiance = np.zeros(data.shape, dtype='f4')
    for i in range(nbands):
        rawdata = data[i]
        good = rawdata > dosDN[i]
        good &= rawdata < 65536
        gooddata = (rawdata[good] - dosDN[i]) * scalefactor_bands[i] - bias_bands[i]
        gooddata[gooddata < 0] = 0
        radiance[i, good] = gooddata

    return radiance
