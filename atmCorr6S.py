import os
import re
import sys
import time
import logging
import multiprocessing
from math import exp, cos, radians
from xml.etree import ElementTree as ET

import numpy as np
from scipy.ndimage import filters
from scipy import interpolate

from graspy.gdal_utils import array_to_gtiff

import bathyUtilities
from Py6S import SixS, AtmosProfile, AeroProfile, AtmosCorr, Wavelength, Geometry


logger = logging.getLogger('atmProcessing.atmCorr6S')

def setup_SixS(args):
    AOT, PWV, ozone, bandFilter, aeroProfile, metadataFile, startWV, endWV = args

    # Have different paths to 6S and spectral response curves on Windows where,
    # I run the code mostly through Spyder and on Linux (CentOS/RedHat) where
    # I run mostly the complied program

    wd = os.path.dirname(__file__)
    PATH_6S = os.path.join(wd, 'dependency', "sixsV1.1")

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
    startWV, endWV, bandFilters = bathyUtilities.readBandFiltersFromCSV(path, sensor, isPan)
    startWV /= 1000.0
    endWV /= 1000.0

    # Also need to resample the band filters from 1nm to 2.5nm as this is the highest spectral resolution supported by 6S
    for i, band in enumerate(bandFilters):
        bandFilters[i] = bathyUtilities.resampleBandFilters(band, startWV, endWV, 0.0025)

    # Run 6S for each spectral band
    pool = multiprocessing.Pool(nprocs)
    jobArgs = [(atm['AOT'], atm['PWV'], atm['ozone'], bandFilter, aeroProfile, metadataFile, startWV, endWV)
               for bandFilter in bandFilters]
    output = []
    s = None

    start = time.time()
    sys.stdout.write('\r  {0:8.2f}% Atmospheric correction 6S. time: {1:8.2f}'.format(0.0, time.time() - start))
    for i, res in enumerate(pool.imap(fun_SixS, jobArgs)):
        sys.stdout.write('\r  {0:8.2f}% Atmospheric correction 6S. time: {1:8.2f}'.format(100 * i / float(len(jobArgs)-1),
                                                                           time.time() - start))
        output.append(res[0])
        s = res[1]

    pool.close()
    pool.join()
    print ""
    return s,  output


def performAtmCorrection(inImg, correctionParams6S, radius=1, s=None):
    refl  = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount), dtype='float32')
    pixelSize = inImg.GetGeoTransform()[1]  # assume same horizontal and vertical resolution
    start = time.time()
    sys.stdout.write('\r  {0:8.2f}% Adjacency   correction 6S. time: {1:8.2f}'.format(0.0, time.time() - start))
    for bandNum, correctionParam in enumerate(correctionParams6S):
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
            refl[:, :, bandNum] = adjacencyCorrection(refl[:, :, bandNum], pixelSize, s, radius)
        sys.stdout.write('\r  {0:8.2f}% Adjacency   correction 6S. time: {1:8.2f}'.format(
            100 * bandNum / float(len(correctionParams6S) - 1),
            time.time() - start))


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


def readGeometryPHR1(metadataFile, model6S):

    s = model6S

    tree = ET.parse(metadataFile)

    # get down to the appropirate node
    root = tree.getroot()
    Geometric_Data = root.findall('Geometric_Data')[0]
    Use_Area = Geometric_Data.findall('Use_Area')[0]
    for LGV in Use_Area.findall('Located_Geometric_Values'):
        # get angles for centre of the image
        if LGV.findall('LOCATION_TYPE')[0].text == "Center":
            Acquisition_Angles = LGV.findall('Acquisition_Angles')[0]
            satAz = float(Acquisition_Angles.findall('AZIMUTH_ANGLE')[0].text)
            satZen = float(Acquisition_Angles.findall('INCIDENCE_ANGLE')[0].text)
            Solar_Incidences = LGV.findall('Solar_Incidences')[0]
            sunAz = float(Solar_Incidences.findall('SUN_AZIMUTH')[0].text)
            sunEl = float(Solar_Incidences.findall('SUN_ELEVATION')[0].text)

            # get month and day
            timeStr = LGV.findall('TIME')[0].text
            dateRegex = '\d{4}-(\d{2})-(\d{2})T.*'
            match = re.match(dateRegex, timeStr)
            if match:
                month = int(match.group(1))
                day = int(match.group(2))

            break

    sunZen = 90.0 - sunEl

    s.geometry = Geometry.User()
    s.geometry.solar_z = sunZen
    s.geometry.solar_a = sunAz
    s.geometry.view_z = satZen
    s.geometry.view_a = satAz
    s.geometry.day = day
    s.geometry.month = month


def readGeometryWV2(metadataFile, s):
    # read viewing gemotery from WV2 metadata file
    meanSunElRegex = "\s*meanSunEl\s*=\s*(.*);"
    meanSunAzRegex = "\s*meanSunAz\s*=\s*(.*);"
    meanSatElRegex = "\s*meanSatEl\s*=\s*(.*);"
    meanSatAzRegex = "\s*meanSatAz\s*=\s*(.*);"
    # depending on the product type there can be either firstLineTime or earliestAcqTime in the metadata file
    firstLineTimeRegex = "\s*firstLineTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"
    earliestAcqTimeRegex = "\s*earliestAcqTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"

    month, day, sunEl, sunAz, satEl, satAz = 0, 0, 0, 0, 0, 0

    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(firstLineTimeRegex, line)
            if not match:
                 match = re.match(earliestAcqTimeRegex, line)

            if match:
                 month = int(match.group(2))
                 day   = int(match.group(3))

            match = re.match(meanSunElRegex, line)
            if match:
                sunEl = float(match.group(1))

            match = re.match(meanSunAzRegex, line)
            if match:
                sunAz = float(match.group(1))

            match = re.match(meanSatElRegex, line)
            if match:
                satEl = float(match.group(1))

            match = re.match(meanSatAzRegex, line)
            if match:
                satAz = float(match.group(1))

    sunZen = 90.0 - sunEl
    satZen = 90.0 - satEl

    s.geometry = Geometry.User()
    s.geometry.solar_z = sunZen
    s.geometry.solar_a = sunAz
    s.geometry.view_z = satZen
    s.geometry.view_a = satAz
    s.geometry.day = day
    s.geometry.month = month


def readGeometryL8(metadataFile, model6S, extent):

    s = model6S

    # read viewing gemotery from Landsat metadata file
    sunElRegex = "\s*SUN_ELEVATION\s*=\s*(.*)\s*"
    sunAzRegex = "\s*SUN_AZIMUTH\s*=\s*(.*)\s*"
    dateAcquiredRegex = "\s*DATE_ACQUIRED\s*=\s*\d{4}-(\d{2})-(\d{2})\s*"
    minXRegex = "\s*CORNER_LL_PROJECTION_X_PRODUCT\s*=\s*(\d+\.\d+)\s*"
    maxXRegex = "\s*CORNER_UR_PROJECTION_X_PRODUCT\s*=\s*(\d+\.\d+)\s*"

    month = 0; day = 0; sunEl = 0.0; sunAz = 0.0; satZen = 0.0; satAz = 0.0;
    minX = 0; maxX = 0

    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(dateAcquiredRegex, line)
            if match:
                month = int(match.group(1))
                day = int(match.group(2))
            match = re.match(sunElRegex , line)
            if match:
                sunEl = float(match.group(1))
            match = re.match(sunAzRegex , line)
            if match:
                sunAz = float(match.group(1))
            match = re.match(minXRegex, line)
            if match:
                minX = float(match.group(1))
            match = re.match(maxXRegex, line)
            if match:
                maxX = float(match.group(1))

    sunZen = 90 - sunEl

    s.geometry = Geometry.User()
    s.geometry.solar_z = sunZen
    s.geometry.solar_a = sunAz
    if extent is None:
        s.geometry.view_z = satZen
    else:
        # extent is [minX, maxY, maxX, minY]
        extentMiddleX = (extent[0]+extent[2])/2
        imageMiddleX = (minX + maxX)/2
        # L8 has 15 deg field of view so VZA ranges between 0 and 7.5 degrees
        s.geometry.view_z = abs(imageMiddleX - extentMiddleX)/(imageMiddleX - minX) * 7.5
    s.geometry.view_a = satAz
    s.geometry.day = day
    s.geometry.month = month


def readGeometryS2(metadataFile, model6S):
    # Get the granule metadata from metadata dictionary
    sunZen = float(metadataFile[metadataFile['current_granule']]['sun_zenit'])
    sunAz = float(metadataFile[metadataFile['current_granule']]['sun_azimuth'])
    sensorZen = float(metadataFile[metadataFile['current_granule']]['sensor_zenit'])
    sensorAz = float(metadataFile[metadataFile['current_granule']]['sensor_azimuth'])

    dateTime = metadataFile['product_start']
    dayMonthRegex = "\d{4}-(\d{2})-(\d{2})T"
    match = re.match(dayMonthRegex, dateTime)
    if match:
        month = int(match.group(1))
        day = int(match.group(2))

    s = model6S
    s.geometry = Geometry.User()
    s.geometry.solar_z = sunZen
    s.geometry.solar_a = sunAz
    s.geometry.view_z = sensorZen
    s.geometry.view_a = sensorAz
    s.geometry.day = day
    s.geometry.month = month
