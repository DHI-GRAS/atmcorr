from __future__ import division
import logging

import numpy as np
from Py6S import SixS
from Py6S import AtmosProfile
from Py6S import AeroProfile
from Py6S import AtmosCorr
from Py6S import Wavelength
from Py6S import Geometry

from atmcorr import utils
from atmcorr.adjacency_correction import adjacency_correction

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
        sensor,
        AOT, PWV, ozone,
        aeroProfile, geometry_dict):
    """Set up 6S instance

    Parameters
    ----------
    sensor : str
        sensor name
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


def run_sixs_for_wavelength(args):
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
    """
    mysixs, start_wv, end_wv, rcurve = args
    mysixs.wavelength = Wavelength(start_wv, end_wv, rcurve)
    mysixs.run()
    xdict = {
            'xa': mysixs.outputs.coef_xa,  # inverse of transmitance
            'xb': mysixs.outputs.coef_xb,  # scattering term of the atmosphere
            'xc': mysixs.outputs.coef_xc}  # reflectance of atmosphere for isotropic light (albedo)
    retvals = [
            mysixs.geometry.view_z,
            mysixs.outputs.optical_depth_total.total,
            mysixs.outputs.transmittance_global_gas.upward,
            mysixs.outputs.transmittance_total_scattering.upward]
    return (xdict, retvals)


def get_correction_params(
        processor_pool,
        sensor,
        atm,
        rcurves_dict,
        geometry_dict,
        aeroProfile="Continental"):
    """Get correction parameters

    Parameters
    ----------
    processor_pool : multiprocessing.Pool
        processor pool
    sensor : str
        sensor name
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
            sensor,
            atm['AOT'],
            atm['PWV'],
            atm['ozone'],
            aeroProfile,
            geometry_dict)

    # Run 6S for each spectral band
    jobs = [
            (mysixs, rcurves_dict['start_wv'], rcurves_dict['end_wv'], rcurve)
            for rcurve in rcurves_dict['rcurves']]
    output = []
    mysixs = None
    for res in processor_pool.imap(run_sixs_for_wavelength, jobs):
        output.append(res[0])
        mysixs = res[1]
    return mysixs, output


def perform_correction(data, corrparams, pixel_size, radius=1, adjCorr=False, mysixs=None):
    """Perform atmospheric correction

    Parameters
    ----------
    data : ndarray shape(nbands, ny, nx)
        input data
    corrparams : record array shape=(nbands, nj, ni) fields=['xa', 'xb', 'xc']
        correction parameters (can be tiled)
        nj, ni are the tile indices
    pixel_size : int
        pixel size
    radius : int
        radius
    adjCorr : bool
        do adjacency correct
    mysixs : SixS instance, optional
        SixS model
        required for adjCorr
    """
    if adjCorr and mysixs is None:
        raise ValueError('adjCorr requires 6S instance')

    nbands, nj, ni = corrparams.shape
    ntiles = nj * ni

    if not nbands == data.shape[0]:
        raise ValueError(
                'First dimension of corrparams should correspond to '
                'first (band) dimension of data.')

    reflectance = np.full(data.shape, np.nan, dtype='f4')

    for i in range(nbands):
        corrparams_band = corrparams[i]
        # Read uncorrected radiometric data and correct
        radiance = data[i]

        # Interpolate the 6S correction parameters from one per image tile to
        # one per image pixel
        if ntiles == 1:
            # take single value
            xa = corrparams_band['xa'][0, 0]
            xb = corrparams_band['xb'][0, 0]
            xc = corrparams_band['xc'][0, 0]
        else:
            # interpolate
            xa = utils.imresize(corrparams_band['xa'], radiance.shape)
            xb = utils.imresize(corrparams_band['xb'], radiance.shape)
            xc = utils.imresize(corrparams_band['xc'], radiance.shape)

        # Perform the atmospheric correction
        with np.errstate(invalid='ignore'):
            # safe to ignore: failing elements already NaN
            mask_nan = radiance == 0
        y = xa * radiance - xb
        y[mask_nan] = np.nan
        refl_band = y / (xc * y + 1.0)
        with np.errstate(invalid='ignore'):
            # safe to ignore: NaN elements remain NaN
            refl_band[refl_band < 0] = 0
        reflectance[i] = refl_band

    # Perform adjecency correction if required
    if adjCorr:
        for i in range(nbands):
            reflectance[i] = adjacency_correction(
                    reflectance[i],
                    pixel_size,
                    mysixs,
                    radius)
    return reflectance
