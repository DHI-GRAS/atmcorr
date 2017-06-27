import os
import sys
import logging
import multiprocessing
from math import exp, cos, radians

import numpy as np
from scipy.ndimage import filters
from scipy import interpolate
from tqdm import trange
from tqdm import tqdm

from gdal_utils.gdal_utils import array_to_gtiff

import bathyUtilities
from Py6S import SixS, AtmosProfile, AeroProfile, AtmosCorr, Wavelength

from . import band_filters

logger = logging.getLogger(__name__)

PATH_6S = os.path.join(os.path.dirname(__file__), 'dependency', "sixsV1.1")


def setup_SixS(args):
    AOT, PWV, ozone, bandFilter, aeroProfile, metadataFile, startWV, endWV = args

    # Have different paths to 6S and spectral response curves on Windows where,
    # I run the code mostly through Spyder and on Linux (CentOS/RedHat) where
    # I run mostly the complied program

    s = SixS(PATH_6S)

    #########################################################
    # Set 6S BRDF model to 1 m/s wind ocean with typical salinity and pigment concentration
    # s.ground_reflectance = GroundReflectance.HomogeneousOcean(1.0, 0, -1, 0.5)

    s.atmos_profile = AtmosProfile.UserWaterAndOzone(PWV, ozone)

    aeroProfileDict = {"No Aerosols": AeroProfile.NoAerosols,
                       "Continental": AeroProfile.Continental,
                       "Maritime": AeroProfile.Maritime,
                       "Urban": AeroProfile.Urban,
                       "Desert": AeroProfile.Desert,
                       "BiomassBurning": AeroProfile.BiomassBurning,
                       "Stratospheric": AeroProfile.Stratospheric}

    s.aero_profile = AeroProfile.PredefinedType(aeroProfileDict[aeroProfile])
    # get from MOD04 or MODATML2
    s.aot550 = AOT

    #########################################################
    # Set 6S altitude
    s.altitudes.set_target_sea_level()
    s.altitudes.set_sensor_satellite_level()

    #########################################################
    # Set 6S atmospheric correction
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(10)

    #########################################################
    # Set 6S geometry
    readGeometryWV2(metadataFile, s)

    s.wavelength = Wavelength(startWV, endWV, bandFilter)
    return s


def fun_SixS(args):
    s = setup_SixS(args)
    s.run()
    xa = s.outputs.coef_xa  # inverse of transmitance
    xb = s.outputs.coef_xb  # scattering term of the atmosphere
    xc = s.outputs.coef_xc  # reflectance of atmosphere for isotropic light (albedo)

    return {'xa': xa, 'xb': xb, 'xc': xc}, [s.geometry.view_z, s.outputs.optical_depth_total.total,
            s.outputs.transmittance_global_gas.upward, s.outputs.transmittance_total_scattering.upward]


def getCorrectionParams6S(metadataFile, atm={'AOT': -1, 'PWV': -1, 'ozone': -1}, sensor="WV2", isPan=False,
                          aeroProfile="Continental", extent=None, nprocs=False):

    if not nprocs:
        nprocs = multiprocessing.cpu_count()

    if getattr(sys, 'frozen', False):
        application_path = os.path.dirname(sys.executable)
    elif __file__:
        application_path = os.path.dirname(__file__)

    path = os.path.join(application_path, 'dependency', 'sensorResponseCurves', sensor + ".txt")

    # Set 6S band filters
    startWV, endWV, bandFilters = band_filters.readBandFiltersFromCSV(path, sensor, isPan)
    startWV /= 1000.0
    endWV /= 1000.0

    # Also need to resample the band filters from 1nm to 2.5nm as this is the highest spectral resolution supported by 6S
    for i, band in enumerate(bandFilters):
        bandFilters[i] = bathyUtilities.resampleBandFilters(band, startWV, endWV, 0.0025)

    # Run 6S for each spectral band
    pool = multiprocessing.Pool(nprocs)
    jobs = [(atm['AOT'], atm['PWV'], atm['ozone'], bandFilter, aeroProfile, metadataFile, startWV, endWV)
            for bandFilter in bandFilters]
    logger.info('Running {} atmospheric correction jobs on {} processors'.format(len(jobs), nprocs))
    output = []
    s = None
    for res in tqdm(pool.imap(fun_SixS, jobs), desc='Atmospheric Correction 6S', unit='job'):
        output.append(res[0])
        s = res[1]
    pool.close()
    pool.join()
    return s, output


def performAtmCorrection(inImg, correctionParams6S, radius=1, s=None):
    refl = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount), dtype='float32')
    pixelSize = inImg.GetGeoTransform()[1]  # assume same horizontal and vertical resolution
    for bandNum in trange(len(correctionParams6S), desc='performAtmCorrection', unit='band'):
        correctionParam = correctionParams6S[bandNum]
        # Read uncorrected radiometric data and correct
        radianceData = inImg.GetRasterBand(bandNum + 1).ReadAsArray()

        # Interpolate the 6S correction parameters from one per image tile to
        # one per image pixel
        xa = explodeCorrectionParam(correctionParam['xa'], radianceData.shape)
        xb = explodeCorrectionParam(correctionParam['xb'], radianceData.shape)
        xc = explodeCorrectionParam(correctionParam['xc'], radianceData.shape)

        # Perform the atmospheric correction
        y = np.where(radianceData == 0, np.nan, xa * radianceData - xb)
        refl[:, :, bandNum] = np.where(np.isnan(y), 0, np.maximum(y / (1.0 + xc * y), 0.0))
        # Perform adjecency correction if required
        if s is not None:
            logger.info('Performing adjacency correction for band {} ...'.format(bandNum+1))
            refl[:, :, bandNum] = adjacencyCorrection(refl[:, :, bandNum], pixelSize, s, radius)
    return array_to_gtiff(refl, "MEM", inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)


# Interpolate the an array with parameters into the new shape
def explodeCorrectionParam(param, newShape):
    array = np.array(param)

    # If there is only one value then assign it to each cell in the new array
    if array.size == 1:
        return np.zeros(newShape)+array[0]
    # Otherwise interpolate
    else:
        # At least two points in each dimension are needed for linerar interpolation
        if array.shape[0] == 1:
            array = np.vstack((array, array))
        if array.shape[1] == 1:
            array = np.hstack((array, array))

        # Assume that the parameters are regularly spaced within the new array
        xStep = newShape[1]/(array.shape[1])
        xLoc = np.array([(x+0.5)*xStep for x in range(array.shape[1])])
        yStep = newShape[0]/(array.shape[0])
        yLoc = np.array([(y+0.5)*yStep for y in range(array.shape[0])])

        f = interpolate.interp2d(xLoc, yLoc, array, fill_value = None)
        explodedParam = f(range(newShape[1]), range(newShape[0]))

        return explodedParam

# NEEDS TO BE DOUBLE CHECKED
# Following Ouaidrari & Vermote 1999: Operational Atmospheric Correction of Landsat TM Data
def adjacencyCorrection(refl, pixelSize, s, radius = 1.0):

    # definition below eq (4)
    u_v, tau,  T_dir, T_dif = cos(radians(s[0])), s[1], s[2], s[3]
    T = 1 - ((1 - T_dif) + (1 - T_dir))
    """
    u_v = cos(radians(s.geometry.view_z))
    tau = s.outputs.optical_depth_total.total
    T_dir = s.outputs.transmittance_global_gas.upward
    T_dif = s.outputs.transmittance_total_scattering.upward
    T= 1 - ((1-T_dif) + (1-T_dir))
    """

    # Fill in any NaN values, particularly at the edges of the image
    mask = np.isnan(refl)
    refl[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), refl[~mask])

    # Calculate average reflectance of adjacent pixels
    # The adjacency effect can come from pixels within 1km of the central pixel (Verhoef et al., 2003) so  sigma should be half of that in gaussian filter
    sigma = radius/pixelSize
    adjRefl = filters.gaussian_filter(refl, sigma)

    # eq (8)
    t_d = T_dif - exp(-tau/u_v)
    refl = (refl*T - adjRefl*t_d)/exp(-tau/u_v)

    # http://www.cesbio.ups-tlse.fr/multitemp/?p=2277
    #albedo = s.outputs.spherical_albedo.total
    #refl = ( refl*T*(1-refl*s)/(1-adjRefl*s) - adjRefl*t_d ) / exp(-tau/u_v)
    #T = 1 - ((1-T_dif) + (1-T_dir))
    #refl = (refl*T*(1-refl*albedo)/(1-adjRefl*albedo) - adjRefl*T_dif) / T_dir

    # Clean up
    refl[mask] = np.NaN
    refl[refl < 0.0] = 0.0
    return refl
