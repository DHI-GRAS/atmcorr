import os
import logging
import multiprocessing

import numpy as np
import gdal_utils.gdal_utils as gu
from scipy.ndimage import filters
from scipy import interpolate
from tqdm import trange
from tqdm import tqdm
from Py6S import SixS
from Py6S import AtmosProfile
from Py6S import AeroProfile
from Py6S import AtmosCorr
from Py6S import Wavelength
from Py6S import Geometry
import srcurves

from atmospheric_correction import viewing_geometry as vg

logger = logging.getLogger(__name__)

PATH_6S = os.path.join(os.path.dirname(__file__), 'dependency', "sixsV1.1")

_geometry_attrs_to_keys = {
        'solar_z': 'sun_zenith',
        'solar_a': 'sun_azimuth',
        'view_z': 'sensor_zenith',
        'view_a': 'sensor_azimuth',
        'day': 'day',
        'month': 'month'}


def setup_SixS(
        sensor,
        AOT, PWV, ozone, bandFilter,
        aeroProfile, geometry_dict, start_wv, end_wv):

    mysixs = SixS(PATH_6S)

    # Set 6S BRDF model to 1 m/mysixs wind ocean with typical salinity and pigment concentration
    # mysixs.ground_reflectance = GroundReflectance.HomogeneousOcean(1.0, 0, -1, 0.5)
    mysixs.atmos_profile = AtmosProfile.UserWaterAndOzone(PWV, ozone)

    aeroProfileDict = {
            "No Aerosols": AeroProfile.NoAerosols,
            "Continental": AeroProfile.Continental,
            "Maritime": AeroProfile.Maritime,
            "Urban": AeroProfile.Urban,
            "Desert": AeroProfile.Desert,
            "BiomassBurning": AeroProfile.BiomassBurning,
            "Stratospheric": AeroProfile.Stratospheric}

    mysixs.aero_profile = AeroProfile.PredefinedType(aeroProfileDict[aeroProfile])
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

    mysixs.wavelength = Wavelength(start_wv, end_wv, bandFilter)
    return mysixs


def fun_SixS(setup_args):
    mysixs = setup_SixS(setup_args)
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


def getCorrectionParams6S(
        sensor,
        mtdfile,
        atm,
        isPan=False,
        aeroProfile="Continental",
        extent=None,
        nprocs=None):

    if nprocs is None:
        nprocs = multiprocessing.cpu_count()

    # Set 6S band filters
    start_wv, end_wv, rcurves = srcurves.get_response_curves(sensor, isPan)
    start_wv /= 1000.0
    end_wv /= 1000.0

    # Also need to resample the band filters from 1nm to 2.5nm
    # as this is the highest spectral resolution supported by 6S
    for i, band in enumerate(rcurves):
        rcurves[i] = srcurves.resample_response_curves(
                band, start_wv, end_wv, 0.0025)

    geometry_dict = vg.get_geometry(sensor, mtdfile)

    # Run 6S for each spectral band
    pool = multiprocessing.Pool(nprocs)
    jobs = [(
        sensor,
        atm['AOT'],
        atm['PWV'],
        atm['ozone'],
        bandFilter,
        aeroProfile,
        geometry_dict,
        start_wv,
        end_wv)
        for bandFilter in rcurves]
    logger.info(
            'Running %d atmospheric correction jobs on %d processors',
            len(jobs), nprocs)
    output = []
    mysixs = None
    for res in tqdm(
            pool.imap(fun_SixS, jobs),
            desc='Atmospheric Correction 6S',
            unit='job'):
        output.append(res[0])
        mysixs = res[1]
    pool.close()
    pool.join()
    return mysixs, output


def performAtmCorrection(img, correctionParams6S, radius=1, mysixs=None):
    shape = (img.RasterYSize, img.RasterXSize, img.RasterCount)
    refl = np.zeros(shape, dtype='float32')
    # assume same horizontal and vertical resolution
    pixelSize = img.GetGeoTransform()[1]

    nbands = len(correctionParams6S)

    for b in trange(nbands, desc='atmcorr', unit='band'):
        correctionParam = correctionParams6S[b]
        # Read uncorrected radiometric data and correct
        radiance = img.GetRasterBand(b + 1).ReadAsArray()

        # Interpolate the 6S correction parameters from one per image tile to
        # one per image pixel
        xa = explodeCorrectionParam(correctionParam['xa'], radiance.shape)
        xb = explodeCorrectionParam(correctionParam['xb'], radiance.shape)
        xc = explodeCorrectionParam(correctionParam['xc'], radiance.shape)

        # Perform the atmospheric correction
        y = xa * radiance - xb
        y[(radiance == 0)] = np.nan
        refl_band = y / (1.0 + xc * y)
        refl_band = np.maximum(refl_band, 0.0)
        refl_band[np.isnan(refl_band)] = 0.0
        refl[:, :, b] = refl_band

    # Perform adjecency correction if required
    if mysixs is not None:
        logger.info(
                'Performing adjacency correction')
        for b in trange(nbands, desc='adjcorr', unit='band'):
            refl[:, :, b] = adjacencyCorrection(
                    refl[:, :, b],
                    pixelSize,
                    mysixs,
                    radius)

    return gu.array_to_gtiff(
            refl, "MEM",
            img.GetProjection(), img.GetGeoTransform(),
            banddim=2)


def explodeCorrectionParam(param, newshape):
    # Interpolate the an array with parameters into the new shape
    data = np.data(param)
    ny, nx = data.shape

    # If there is only one value then assign it to each cell in the new data
    if data.size == 1:
        return np.zeros(newshape) + data[0]

    # Otherwise interpolate
    # At least two points in each dimension are needed for linerar interpolation
    if ny == 1:
        data = np.vstack((data, data))
    if nx == 1:
        data = np.hstack((data, data))

    # Assume that the parameters are regularly spaced within the new data
    dx = newshape[1] / nx
    xx = (np.arange(nx) + 0.5) * dx

    dy = newshape[0] / ny
    yy = (np.arange(ny) + 0.5) * dy

    f = interpolate.interp2d(xx, yy, data, fill_value=None)
    explodedParam = f(range(newshape[1]), range(newshape[0]))
    return explodedParam


def adjacencyCorrection(refl, pixelSize, mysixs, radius=1.0):
    # NEEDS TO BE DOUBLE CHECKED
    # Following Ouaidrari & Vermote 1999: Operational Atmospheric Correction of Landsat TM Data

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
    sigma = radius/pixelSize
    adjRefl = filters.gaussian_filter(refl, sigma)

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
