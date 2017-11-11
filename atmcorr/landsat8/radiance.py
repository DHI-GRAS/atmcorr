import numpy as np

from satmeta.l8 import meta as l8meta


def dn_to_radiance(dndata, mtdfile, band_ids=slice(None, None)):
    """Compute radiance from digital numbers

    Parameters
    ----------
    dndata : ndarray shape (nbands, ny, nx)
        digital numbers data
    mtdfile : str
        path to metadata file
    band_ids : list of int, optional
        band IDs (0-based index of bands) in dndata

    Returns
    -------
    ndarray
        radiance data

    Source
    ------
    https://yceo.yale.edu/how-convert-landsat-dns-top-atmosphere-toa-reflectance
    """
    with open(mtdfile) as fin:
        metadata = l8meta.parse_metadata(fin)
    bias = np.asarray(metadata['RADIOMETRIC']['ADD'], 'float32')[band_ids]
    gain = np.asarray(metadata['RADIOMETRIC']['MULT'], 'float32')[band_ids]
    radiance = dndata.astype('float32') * gain[:, np.newaxis, np.newaxsi] + bias
    return radiance
