from __future__ import division
import logging

import numpy as np
from Py6S import SixS
from Py6S import AtmosProfile
from Py6S import AeroProfile
from Py6S import AtmosCorr
from Py6S import Wavelength
from Py6S import Geometry

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
        AOT, PWV, ozone,
        aeroProfile, geometry_dict):
    """Set up 6S instance

    Parameters
    ----------
    AOT, PWV, ozone : float
        atmospheric constituents
    aeroProfile : str
        name of profile
        must be attribute of Py6S.AeroProfile
    geometry_dict : dict
        from viewing_geometry
    """
    mysixs = SixS()

    # Set 6S BRDF model to 1 m/mysixs wind ocean with typical salinity and pigment concentration
    # mysixs.ground_reflectance = GroundReflectance.HomogeneousOcean(1.0, 0, -1, 0.5)
    mysixs.atmos_profile = AtmosProfile.UserWaterAndOzone(PWV, ozone)

    try:
        aeroProfile_attr = getattr(AeroProfile, aeroProfile)
    except (TypeError, AttributeError):
        raise ValueError(
                'aeroProfile \'{}\' not recognized. Try one of {}'
                ''.format(aeroProfile, AERO_PROFILES))

    mysixs.aero_profile = AeroProfile.PredefinedType(aeroProfile_attr)
    # get from MOD04 or MODATML2
    mysixs.aot550 = AOT

    # Set 6S altitude
    mysixs.altitudes.set_target_sea_level()
    mysixs.altitudes.set_sensor_satellite_level()

    # Set 6S atmospheric correction
    mysixs.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(10)

    # Set 6S geometry
    mysixs.geometry = Geometry.User()
    for attrname in GEOMETRY_ATTRS:
        key = GEOMETRY_ATTRS[attrname]
        value = geometry_dict[key]
        if value is not None:
            setattr(mysixs.geometry, attrname, value)

    return mysixs


def generate_jobs(
        atm,
        rcurves_dict,
        geometry_dict,
        aeroProfile="Continental"):
    """Sets up Py6S instance and returns job for run_sixs_for_wavelength

    Parameters
    ----------
    atm : dict
        atmospheric parameters
    rcurves_dict : dict
        sensor response curve parameters
    geometry_dict : dict
        viewing / Sun angles parameters
    aeroProfile : str
        aero profile for 6S
    """
    # set up SixS instance once
    mysixs = setup_sixs(
            atm['AOT'],
            atm['PWV'],
            atm['ozone'],
            aeroProfile,
            geometry_dict)

    # Run 6S for each spectral band
    jobs = [
            (mysixs, rcurves_dict['start_wv'], rcurves_dict['end_wv'], rcurve)
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
    mysixs.wavelength = Wavelength(start_wv, end_wv, rcurve)
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
        with np.errstate(invalid='ignore'):
            # safe to ignore: failing elements already NaN
            mask_nan = radiance[i] == 0
        y = xa[i] * radiance[i] - xb[i]
        y[mask_nan] = np.nan
        refl_band = y / (xc[i] * y + 1.0)
        with np.errstate(invalid='ignore'):
            # safe to ignore: NaN elements remain NaN
            refl_band[refl_band < 0] = 0
        reflectance[i] = refl_band
    return reflectance
