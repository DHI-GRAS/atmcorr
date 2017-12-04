import numpy as np

from atmcorr.metadata import pleiades as metadata


def dn_to_radiance(dndata, mtdFile, band_ids):
    """Apply radiometric correction to Pleadis image"""
    gain, bias = metadata.get_gain_bias_PHR1(mtdFile)
    gain = [gain[i] for i in band_ids]
    bias = [bias[i] for i in band_ids]
    radata = np.zeros(dndata.shape, dtype='float32')
    radata[dndata == 65536] = np.nan
    for i in range(radata.shape[0]):
        radata[i, ...] = dndata[i] / gain[i] + bias[i]
    with np.errstate(invalid='ignore'):
        radata[radata < 0] = 0
    return radata
