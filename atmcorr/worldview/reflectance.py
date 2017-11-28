from dg_calibration import reflectance


def toa_reflectance(radata, mtdFile, band_ids):
    """Estimate toa reflectance of radiometric data ignoring atmospheric, topographic and BRDF effects

    Parameters
    ----------
    radata : ndarray shape (nbands, ny, nx)
        radiance data
    mtdFile : str
        path to IMD metadata file
    band_ids : sequence of int
        band IDs

    Returns
    -------
    ndarray
        reflectance
    """
    return reflectance.radiance_to_reflectance(radata, mtdFile, band_ids=band_ids)
