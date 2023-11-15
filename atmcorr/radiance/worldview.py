import numpy as np

from ..metadata.worldview import parse_metadata
from ..wv_meta.gain_offset import get_gain_values, get_offset_values
from ..wv_meta.bands import get_values_sorted


def calculate_radiance(dn, gain, offset, absCalFactor, effectiveBandwidth):
    """Compute radiance from digital numbers

    Parameters
    ----------
    dn : ndarray (nbands, ny, nx)
        digital numbers data
    gain, offset : ndarray (nbands)
        gain and offset for all bands in dn data
    absCalFactor, effectiveBandwidth : ndarray (nbands)
        coefficients for all bands in dn data

    Returns
    -------
    ndarray
        radiance
    """
    radiance = np.zeros(dn.shape, dtype='float32')
    with np.errstate(invalid='ignore'):
        for b in range(radiance.shape[0]):
            radiance[b] = gain[b] * dn[b] * absCalFactor[b] / effectiveBandwidth[b] + offset[b]
    return radiance


def get_rad_params(mtd, band_ids=None):
    if mtd['bandId'] != 'Multi':
        raise NotImplementedError(
            'Currently only supporting \'Multi\' (multispectral) metadata. Got \'{}\'.'
            .format(mtd['bandId']))
    if band_ids is None:
        band_ids = slice(None, None)
    sat_id = mtd['satId']
    calvals = {k: get_values_sorted(d, sat_id=sat_id) for k, d in mtd['calibration'].items()}
    calkw = {k: np.array(calvals[k])[band_ids] for k in calvals}
    gain = np.array(get_gain_values(sat_id))[band_ids]
    offset = np.array(get_offset_values(sat_id))[band_ids]
    return dict(calkw, gain=gain, offset=offset)


def dn_to_radiance(dndata, mtdFile, band_ids=None):
    """Get radiance from digital numbers and meta data

    Parameters
    ----------
    dn : ndarray (nbands, ny, nx) or (ny, nx, nbands)
        digital numbers data
    imdfile_or_str : str
        path to IMD file or full contents of one
    band_ids : sequence of int, optional
        0-based index of bands in dn data

    Returns
    -------
    ndarray
        radiance
    """
    mtd = parse_metadata(mtdFile)
    kw = get_rad_params(mtd, band_ids)
    radata = calculate_radiance(dndata, **kw)
    with np.errstate(invalid='ignore'):
        radata[radata < 0] = 0
    return radata
