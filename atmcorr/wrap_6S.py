from __future__ import division
import copy
import logging

import numpy as np
import Py6S

logger = logging.getLogger(__name__)

GEOMETRY_ATTRS = {
        'solar_z': 'sun_zenith',
        'solar_a': 'sun_azimuth',
        'view_z': 'sensor_zenith',
        'view_a': 'sensor_azimuth',
        'day': 'day',
        'month': 'month'}


AERO_PROFILES = [
        'NoAerosols',
        'Continental',
        'Maritime',
        'Urban',
        'Desert',
        'BiomassBurning',
        'Stratospheric']


def setup_sixs(
        atm, geometry, aeroProfile='Continental',
        reflectance_constant=10, is_ocean=False):
    """Set up 6S instance

    Parameters
    ----------
    atm : dict
        AOT, PWV, ozone : float
            atmospheric constituents
    aeroProfile : str
        name of profile
        must be attribute of Py6S.AeroProfile
    geometry : dict
        from viewing_geometry
    reflectance_constant : float
        magic number for AtmosCorr<name>FromReflectance
    is_ocean : bool
        use ocean-optimized settings

    Returns
    -------
    initialized SixS instance
    """
    mysixs = Py6S.SixS()

    # set ground reflectance
    if is_ocean:
        mysixs.ground_reflectance = Py6S.GroundReflectance.HomogeneousOcean(1.0, 0, -1, 0.5)

    # set atmospheric profile
    mysixs.atmos_profile = Py6S.AtmosProfile.UserWaterAndOzone(atm['PWV'], atm['ozone'])
    mysixs.aot550 = atm['AOT']

    # set aero profile
    try:
        aeroProfile_attr = getattr(Py6S.AeroProfile, aeroProfile)
    except (TypeError, AttributeError):
        raise ValueError(
            'aeroProfile \'{}\' not recognized. Try one of {}'
            .format(aeroProfile, AERO_PROFILES))
    else:
        mysixs.aero_profile = Py6S.AeroProfile.PredefinedType(aeroProfile_attr)

    # Set 6S altitude
    mysixs.altitudes.set_target_sea_level()
    mysixs.altitudes.set_sensor_satellite_level()

    # set 6S atmospheric correction method
    mysixs.atmos_corr = (
        Py6S.AtmosCorr.AtmosCorrBRDFFromReflectance(reflectance_constant) if is_ocean else
        Py6S.AtmosCorr.AtmosCorrLambertianFromReflectance(reflectance_constant))

    # Set 6S geometry
    mysixs.geometry = Py6S.Geometry.User()
    for attrname, key in GEOMETRY_ATTRS.items():
        value = geometry[key]
        if value is not None:
            setattr(mysixs.geometry, attrname, value)

    return mysixs


def generate_jobs(rcurves_dict, sixs_params):
    """Sets up Py6S instance and returns job for run_sixs_for_wavelength

    Parameters
    ----------
    rcurves_dict : dict
        sensor response curve parameters
    sixs_params : dict
        keyword arguments for setup_sixs

    Returns
    -------
    list of tuples (SixS, float, float, ndarray)
        SixS instance
        start and env wavelengths
        sensor response curve
    """
    # set up SixS instance once
    mysixs = setup_sixs(**sixs_params)

    # Run 6S for each spectral band
    jobs = [
        (copy.deepcopy(mysixs), rcurves_dict['start_wv'], rcurves_dict['end_wv'], rcurve)
        for rcurve in rcurves_dict['rcurves']]
    return jobs


def run_sixs_job(args):
    """Run sixs for a specific wavelength

    Parameters
    ----------
    args : sequence
        mysixs, start_wv, end_wv, rcurve

    Elements
    --------
    mysixs : SixS instance
        initialized 6S instance
    start_wv, end_wv : float
        wave length range
    rcurve : ndarray
        sensor response curve

    Returns
    -------
    dict
        correction parameters
    list of float
        adjacency correction parameters
    list
        arguments passed through
    """
    mysixs, start_wv, end_wv, rcurve = args[:4]
    moreargs = args[4:]
    mysixs.wavelength = Py6S.Wavelength(start_wv, end_wv, rcurve)
    mysixs.run()
    xdict = {
        'xa': mysixs.outputs.coef_xa,  # inverse of transmitance
        'xb': mysixs.outputs.coef_xb,  # scattering term of the atmosphere
        'xc': mysixs.outputs.coef_xc}  # reflectance of atmosphere for isotropic light (albedo)
    adjcorr_params = [
        mysixs.geometry.view_z,
        mysixs.outputs.optical_depth_total.total,
        mysixs.outputs.transmittance_global_gas.upward,
        mysixs.outputs.transmittance_total_scattering.upward]
    return (xdict, adjcorr_params, moreargs)


def perform_correction(radiance, corrparams):
    """Perform atmospheric correction

    Parameters
    ----------
    radiance : ndarray shape(nbands, ny, nx)
        input data
    corrparams : recarray or dict of float
        correction parameters
        fields/keys: xa, xb, xc
    """
    xa, xb, xc = (corrparams[field] for field in ['xa', 'xb', 'xc'])
    nbands = radiance.shape[0]
    reflectance = np.full(radiance.shape, np.nan, dtype='f4')

    for i in range(nbands):
        # Perform the atmospheric correction
        y = xa[i] * radiance[i] - xb[i]
        refl_band = y / (xc[i] * y + 1.0)
        with np.errstate(invalid='ignore'):
            # safe to ignore: NaN elements remain NaN
            refl_band[refl_band < 0] = 0
        reflectance[i] = refl_band
    return reflectance
