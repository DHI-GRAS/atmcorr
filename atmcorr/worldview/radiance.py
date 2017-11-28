import numpy as np

from atmcorr import dos

from dg_calibration import radiance


def toa_radiance(dndata, mtdFile, band_ids, doDOS=False):
    """Compute TOA radiance from DigitalGlobe digital numbers

    Parameters
    ----------
    dndata : ndarray shape (nbands, ny, nx)
        digital numbers data
    mtdFile : str
        path to IMD metadata file
    band_ids : sequence of int
        band IDs
    doDOS : bool
        perform dark-object-subtraction before
    """
    nbands = dndata.shape[0]
    if nbands != len(band_ids):
        raise ValueError(
            'First dimension in dndata ({}) must have same size as band_ids ({}).'
            .format(nbands, len(band_ids)))

    # roll band axis last for broadcasting
    dndata = np.rollaxis(dndata, 0, 3)

    # perform dark object substraction
    if doDOS:
        dosDN = dos.do_dos(dndata)
    else:
        dosDN = np.zeros(nbands)

    dndata = dndata.astype('float32')

    good = dndata > dosDN
    good &= dndata < 65536

    dndata[~good] = np.nan

    if doDOS:
        dndata -= dosDN

    radata = radiance.dn_to_radiance(dndata, mtdFile, band_ids=band_ids)
    radata[radata < 0] = np.nan
    return np.rollaxis(radata, 2, 0)
