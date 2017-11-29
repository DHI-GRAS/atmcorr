import logging

import numpy as np

import satmeta.s2.meta as s2meta

logger = logging.getLogger(__name__)


def toa_reflectance_to_radiance(dndata, mtdFile, mtdFile_tile, band_ids):
    """Method taken from the bottom of http://s2tbx.telespazio-vega.de/sen2three/html/r2rusage.html

    Parameters
    ----------
    dndata : ndarray shape(nbands, ny, nx)
        input data
    mtdFile : str
        path to metadata file
    mtdFile_tile : str
        path to granule metadata file
    band_ids : sequence of int
        band IDs (0-based index of bands in img)

    Note
    ----
    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    if mtdFile_tile is None:
        raise ValueError('Tile metadata file required!')

    meta = s2meta.parse_metadata(mtdFile)
    gmeta = s2meta.parse_granule_metadata(mtdFile_tile)
    rc = meta['reflectance_conversion']
    qv = meta['quantification_value']
    irradiance = np.array(meta['irradiance_values'])
    sun_zenith = gmeta['sun_zenith']

    irradiance = irradiance[band_ids]

    # Convert to radiance
    factor = irradiance * np.cos(np.radians(sun_zenith)) / (np.pi * qv)
    radiance = dndata.astype('f4')
    for i in range(dndata.shape[0]):
        radiance[i] /= rc
        radiance[i] *= factor[i]
    return radiance