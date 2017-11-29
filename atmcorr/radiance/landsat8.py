import numpy as np

from satmeta.l8 import meta as l8meta


def dn_to_radiance(dndata, mtdFile, band_ids=None):
    """Compute radiance from digital numbers

    Parameters
    ----------
    dndata : ndarray shape (nbands, ny, nx)
        digital numbers data
    mtdFile : str
        path to metadata file
    band_ids : list of int, optional
        band IDs (0-based index of bands) in data

    Returns
    -------
    ndarray
        radiance data

    Source
    ------
    https://yceo.yale.edu/how-convert-landsat-dns-top-atmosphere-toa-reflectance
    """
    if band_ids is None:
        band_ids = slice(None, None)
    with open(mtdFile) as fin:
        metadata = l8meta.parse_metadata(fin)
    rescaling = metadata['rescaling']
    bias = np.asarray(rescaling['RADIANCE']['ADD'], 'float32')[band_ids]
    gain = np.asarray(rescaling['RADIANCE']['MULT'], 'float32')[band_ids]
    radiance = (
        dndata.astype('float32') *
        gain[:, np.newaxis, np.newaxis] +
        bias[:, np.newaxis, np.newaxis])
    return radiance
