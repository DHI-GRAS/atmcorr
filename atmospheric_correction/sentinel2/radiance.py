import logging

import numpy as np

logger = logging.getLogger(__name__)


def toa_radiance_S2(data, mtdFile, mtdFile_tile, band_ids):
    """Method taken from the bottom of http://s2tbx.telespazio-vega.de/sen2three/html/r2rusage.html

    Parameters
    ----------
    data : ndarray shape(nbands, ny, nx)
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

    metadata = metamod.s2.parse_mtdfile(mtdFile, mtdFile_tile=mtdFile_tile)
    tile = list(metadata['granules'])[0]
    logger.debug('Tile is \'%s\'.', tile)
    rc = metadata['reflectance_conversion']
    qv = metadata['quantification_value']
    irradiance = np.array(metadata['irradiance_values'])
    sun_zenith = metadata['granules'][tile]['sun_zenith']

    irradiance = irradiance[band_ids]

    # Convert to radiance
    logger.info('Radiometric correction')
    factor = irradiance * np.cos(np.radians(sun_zenith)) / (np.pi * qv)
    radiance = data.astype('f4')
    for i in range(data.shape[0]):
        radiance[i] /= rc
        radiance[i] *= factor[i]
    logger.info('Done with radiometric correction.')
    return radiance
