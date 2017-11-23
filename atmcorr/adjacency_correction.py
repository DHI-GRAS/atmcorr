import numpy as np
import scipy.ndimage


def adjacency_correction(refl, pixel_size, params, radius=1.0):
    """Adjacency correction

    Sources
    -------
    Following Ouaidrari & Vermote 1999: Operational Atmospheric Correction of Landsat TM Data
    """
    # TODO: NEEDS TO BE DOUBLE CHECKED

    # definition below eq (4)
    u_v = np.cos(np.radians(params[0]))
    tau = params[1]
    T_dir = params[2]
    T_dif = params[3]
    T = 1 - ((1 - T_dif) + (1 - T_dir))
    """
    u_v = np.cos(np.radians(params.geometry.view_z))
    tau = params.outputs.optical_depth_total.total
    T_dir = params.outputs.transmittance_global_gas.upward
    T_dif = params.outputs.transmittance_total_scattering.upward
    T= 1 - ((1-T_dif) + (1-T_dir))
    """

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
    # albedo = params.outputs.spherical_albedo.total
    # refl = ( refl*T*(1-refl*params)/(1-adjRefl*params) - adjRefl*t_d ) / exp(-tau/u_v)
    # T = 1 - ((1-T_dif) + (1-T_dir))
    # refl = (refl*T*(1-refl*albedo)/(1-adjRefl*albedo) - adjRefl*T_dif) / T_dir

    # Clean up
    refl[refl < 0.0] = 0.0
    refl[mask] = np.nan
    return refl
