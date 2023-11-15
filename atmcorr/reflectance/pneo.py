import numpy as np

import satmeta.pneo.parser as plparser


def toa_reflectance(dndata, mtdFile):
    """Convert to TOA reflectance"""
    metadata = plparser.parse_metadata(mtdFile)
    gain, bias = 10000, 0
    radata = np.zeros(dndata.shape, dtype='float32')
    with np.errstate(invalid='ignore'):
        for i in range(radata.shape[0]):
            radata[i, ...] = dndata[i] / gain[i] + bias[i]
        radata[radata < 0] = 0
    return radata