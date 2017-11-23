import numpy as np
import scipy.ndimage


def adjacency_correction(refl, view_z, tau, T_dir, T_dif, pixel_size, radius=1.0):
    """Adjacency correction

    Parameters
    ----------
    refl : ndarray
        reflectance data
    u_v : ndarray
        np.cos(np.radians(sixs.geometry.view_z))
    tau : ndarray
        sixs.outputs.optical_depth_total.total
    T_dir : ndarray
        sixs.outputs.transmittance_global_gas.upward
    T_dif : ndarray
        sixs.outputs.transmittance_total_scattering.upward
    pixel_size : int
        pixel size
    radius : int
        radius

    Sources
    -------
    Following Ouaidrari & Vermote 1999: Operational Atmospheric Correction of Landsat TM Data
    """
    # TODO: NEEDS TO BE DOUBLE CHECKED
    rolled = False
    if np.ndim(view_z) == 1:
        rolled = True
        refl = np.rollaxis(refl, 0, 3)

    # definition below eq (4)
    u_v = np.cos(np.radians(view_z))
    T = 1 - ((1 - T_dif) + (1 - T_dir))

    # Fill in any NaN values, particularly at the edges of the image
    mask = np.isnan(refl)
    refl[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), refl[~mask])

    # Calculate average reflectance of adjacent pixels
    # The adjacency effect can come from pixels within 1km
    # of the central pixel (Verhoef et al., 2003) so
    # sigma should be half of that in gaussian filter
    sigma = radius / pixel_size
    adjRefl = scipy.ndimage.filters.gaussian_filter(refl, sigma)

    # eq (8)
    t_d = T_dif - np.exp(-tau / u_v)
    refl = (refl * T - adjRefl * t_d) / np.exp(-tau / u_v)

    # http://www.cesbio.ups-tlse.fr/multitemp/?p=2277
    # albedo = sixs.outputs.spherical_albedo.total
    # refl = ( refl*T*(1-refl*sixs)/(1-adjRefl*sixs) - adjRefl*t_d ) / exp(-tau/u_v)
    # T = 1 - ((1-T_dif) + (1-T_dir))
    # refl = (refl*T*(1-refl*albedo)/(1-adjRefl*albedo) - adjRefl*T_dif) / T_dir

    # Clean up
    refl[refl < 0.0] = 0.0
    refl[mask] = np.nan
    if rolled:
        return np.rollaxis(refl, 2, 0)
    else:
        return refl
