import numpy as np

import satmeta.pneo.parser as plparser


def dn_to_radiance(dndata, mtdFile, band_ids):
    """Apply radiometric correction to Pleiades Neo image"""
    metadata = plparser.parse_metadata(mtdFile)
    gain_bias = metadata['calibration_values']
    gain, bias = (np.array(gain_bias[key])[band_ids] for key in ['gain', 'bias'])
    radata = np.zeros(dndata.shape, dtype='float32')
    with np.errstate(invalid='ignore'):
        for i in range(radata.shape[0]):
            radata[i, ...] = (dndata[i] / 10000) / gain[i] + bias[i]
        radata[radata < 0] = 0
    return radata