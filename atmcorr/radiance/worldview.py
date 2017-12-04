import numpy as np

from dg_calibration import radiance


def dn_to_radiance(dndata, mtdFile, band_ids):
    """Compute TOA radiance from DigitalGlobe digital numbers

    Parameters
    ----------
    dndata : ndarray shape (nbands, ny, nx)
        digital numbers data
    mtdFile : str
        path to IMD metadata file
    band_ids : sequence of int
        band IDs
    """
    radata = radiance.dn_to_radiance(dndata, mtdFile, band_ids=band_ids)
    radata[dndata == 65536] = np.nan
    with np.errstate(invalid='ignore'):
        radata[radata < 0] = 0
    return radata
