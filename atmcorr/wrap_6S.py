from __future__ import division
import logging
import multiprocessing

import numpy as np
import scipy.ndimage
from Py6S import SixS
from Py6S import AtmosProfile
from Py6S import AeroProfile
from Py6S import AtmosCorr
from Py6S import Wavelength
from Py6S import Geometry

import sensor_response_curves as srcurves
import sensor_response_curves.resample as srcresample

from atmcorr import utils
from atmcorr import exception_handling
from atmcorr import viewing_geometry as vg

logger = logging.getLogger(__name__)

try:
    import numexpr as ne
    logger.debug('Using numexpr.')
except ImportError:
    ne = None
    logger.debug('Not using numexpr.')


_geometry_attrs_to_keys = {
        'solar_z': 'sun_zenith',
        'solar_a': 'sun_azimuth',
        'view_z': 'sensor_zenith',
        'view_a': 'sensor_azimuth',
        'day': 'day',
        'month': 'month'}


_aeroProfile_names = [
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
                ''.format(aeroProfile, _aeroProfile_names))

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
    for attrname in _geometry_attrs_to_keys:
        key = _geometry_attrs_to_keys[attrname]
        setattr(mysixs.geometry, attrname, geometry_dict[key])

    return mysixs


def run_sixs(args):
    """
    Parameters
    ----------
    args : sequence
        mysixs, start_wv, end_wv, rcurve

    where

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
        sensor,
        mtdFile,
        atm,
        band_ids,
        mtdFile_tile=None,
        aeroProfile="Continental",
        extent=None,
        nprocs=None):
    """Get correction parameters

    Parameters
    ----------
    sensor : str
        sensor name
    mtdFile : str
        path to metadata file
    atm : dict
        atmospheric parameters
    band_ids : list of int
        bands in input data
        0-based index wrt. original product
    mtdFile_tile : str
        path to tile metadata file
        required for Sentinel 2
    aeroProfile : str
        aero profile for 6S
    extent : list of float
        image extent
    nprocs : int
        number of processors to use
    """
    logger.info('Getting correction parameters')
    if nprocs is None:
        nprocs = multiprocessing.cpu_count()

    # Set 6S band filters
    wavelength, rcurves = srcurves.get_response_curves(
            sensor, band_ids=band_ids)
    wavelength = wavelength.astype('float') / 1e3
    # Also need to resample the band filters from 1nm to 2.5nm
    # as this is the highest spectral resolution supported by 6S
    wavelength, rcurves = srcresample.resample_response_curves(
                wavelength, rcurves, resolution=0.0025)

    start_wv = wavelength[0]
    end_wv = wavelength[-1]

    geometry_dict = vg.get_geometry(sensor, mtdFile, mtdFile_tile)

    nprocs = min((nprocs, len(rcurves)))

    # set up SixS instance once
    mysixs = setup_sixs(
            sensor,
            atm['AOT'],
            atm['PWV'],
            atm['ozone'],
            aeroProfile,
            geometry_dict)

    # Run 6S for each spectral band
    pool = multiprocessing.Pool(nprocs)
    jobs = [
            (mysixs, start_wv, end_wv, rcurve)
            for rcurve in rcurves]

    logger.info(
            'Running %d 6S jobs on %d processors',
            len(jobs), nprocs)
    output = []
    mysixs = None
    for res in pool.imap(run_sixs, jobs):
        output.append(res[0])
        mysixs = res[1]
    pool.close()
    pool.join()
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
        with exception_handling.ignore_nan_warnings():
            if ne is not None:
                y = ne.evaluate('xa * radiance - xb')
                mask = radiance == 0
                y[mask] = np.nan
                refl_band = ne.evaluate('y / (1.0 + xc * y)')
                mask = ne.evaluate('refl_band < 0')
                refl_band[mask] = 0
            else:
                y = xa * radiance
                y -= xb
                mask = radiance == 0
                y[mask] = np.nan
                denom = xc * y
                denom += 1.0
                refl_band = y / denom
                refl_band[refl_band < 0] = 0
        reflectance[i] = refl_band

    # Perform adjecency correction if required
    if adjCorr:
        logger.info('Performing adjacency correction')
        for i in range(nbands):
            reflectance[i] = adjacency_correction(
                    reflectance[i],
                    pixel_size,
                    mysixs,
                    radius)
    return reflectance


def adjacency_correction(refl, pixel_size, mysixs, radius=1.0):
    """Adjacency correction

    Sources
    -------
    Following Ouaidrari & Vermote 1999: Operational Atmospheric Correction of Landsat TM Data
    """
    # TODO: NEEDS TO BE DOUBLE CHECKED

    # definition below eq (4)
    u_v, tau, T_dir, T_dif = np.cos(np.radians(mysixs[0])), mysixs[1], mysixs[2], mysixs[3]
    T = 1 - ((1 - T_dif) + (1 - T_dir))
    """
    u_v = np.cos(np.radians(mysixs.geometry.view_z))
    tau = mysixs.outputs.optical_depth_total.total
    T_dir = mysixs.outputs.transmittance_global_gas.upward
    T_dif = mysixs.outputs.transmittance_total_scattering.upward
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
    refl = (refl*T - adjRefl*t_d) / np.exp(-tau/u_v)

    # http://www.cesbio.ups-tlse.fr/multitemp/?p=2277
    # albedo = mysixs.outputs.spherical_albedo.total
    # refl = ( refl*T*(1-refl*mysixs)/(1-adjRefl*mysixs) - adjRefl*t_d ) / exp(-tau/u_v)
    # T = 1 - ((1-T_dif) + (1-T_dir))
    # refl = (refl*T*(1-refl*albedo)/(1-adjRefl*albedo) - adjRefl*T_dif) / T_dir

    # Clean up
    refl[mask] = np.NaN
    refl[refl < 0.0] = 0.0
    return refl