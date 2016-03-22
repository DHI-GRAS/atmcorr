# -*- coding: utf-8 -*-
"""
Created on Mon Jul 07 12:47:31 2014

@author: rmgu
"""

import re
import platform
from osgeo import gdal, ogr
import numpy as np
import numpy.ma as ma
from scipy.ndimage import filters
from scipy import interpolate
import bathyUtilities
import os
from math import *
from xml.etree import ElementTree as ET
import sys

# import Py6S - should be implmented in a better way
wd = os.path.dirname(__file__)
sys.path.append(os.path.join(wd, "dependency"))
from Py6S import *


def getCorrectionParams6S(metadataFile, inImg, atm = {'AOT':-1, 'PWV':-1, 'ozone':-1}, sensor="WV2", isPan=False, aeroProfile="Continental", extent = None, adjCorrOutput=False):

    # Have different paths to 6S and spectral response curves on Windows where, 
    # I run the code mostly through Spyder and on Linux (CentOS/RedHat) where
    # I run mostly the complied program
    if platform.system() == "Windows":
        wd = os.path.dirname(__file__)
        PATH_6S = os.path.join(wd, 'dependency', "sixsV1.1")
        bandFiltersPath = os.path.join(wd, 'dependency', 'sensorResponseCurves')
    else:
        wd = os.path.dirname(__file__)
        PATH_6S = os.path.join(wd, 'dependency', "sixsV1.1")
        bandFiltersPath = os.path.join(wd, 'dependency', 'sensorResponseCurves')
    
    s = SixS(PATH_6S)
      
    #########################################################
    # Set 6S BRDF model to 1 m/s wind ocean with typical salinity and pigment concentration
    #s.ground_reflectance = GroundReflectance.HomogeneousOcean(1.0, 0, -1, 0.5)     
    
    #########################################################
    # Set 6S atmospheric and aerosol profile
    #s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
    # get from MOD05, MOD07 or MODATML2
    PWV = atm['PWV'] if atm['PWV'] >= 0 else 1.80
    ozone = atm['ozone'] if atm['ozone'] >= 0 else 0.30    
    s.atmos_profile =  AtmosProfile.UserWaterAndOzone(PWV, ozone)
    
    aeroProfileDict = {"No Aerosols": AeroProfile.NoAerosols,
                       "Continental": AeroProfile.Continental,
                       "Maritime": AeroProfile.Maritime,
                       "Urban": AeroProfile.Urban,
                       "Desert": AeroProfile.Desert,
                       "BiomassBurning": AeroProfile.BiomassBurning,
                       "Stratospheric": AeroProfile.Stratospheric}
                       
    s.aero_profile = AeroProfile.PredefinedType(aeroProfileDict[aeroProfile]) 
    # get from MOD04 or MODATML2
    s.aot550 = atm['AOT'] if atm['AOT'] >= 0 else 0.10     
    
    #########################################################
    # Set 6S altitude  
    s.altitudes.set_target_sea_level()    
    s.altitudes.set_sensor_satellite_level()
    
    #########################################################
    # Set 6S atmospheric correction
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(10)
    
    #########################################################
    # Set 6S geometry    
    if sensor == "WV2" or sensor == "WV3":
        readGeometryWV2(metadataFile, s)
    elif sensor == "PHR1A" or sensor == "PHR1B" or sensor == "SPOT6":
        readGeometryPHR1(metadataFile, s)
    elif sensor == "L8" or sensor == "L7":
        readGeometryL8(metadataFile, s, extent)
    elif sensor == "S2A_10m" or sensor == "S2A_60m":
        readGeometryS2(metadataFile, s)

    ##############################################################
    # Set 6S band filters   
    startWV, endWV, bandFilters = bathyUtilities.readBandFiltersFromCSV(os.path.join(bandFiltersPath, sensor+".txt"), sensor, isPan)
    # Convert from nanometers to micrometers since this is what 6S needs     
    startWV = startWV/1000.0
    endWV = endWV/1000.0
    # Also need to resample the band filters from 1nm to 2.5nm as this is the highest spectral resolution supported by 6S        
    for i, band in enumerate(bandFilters):
        bandFilters[i] = bathyUtilities.resampleBandFilters(band, startWV, endWV, 0.0025)
    #########################################################
    # Run 6S for each spectral band
    
    bandNum = 1
    outputParams = [] 
    print "Estimating 6S correction parameters..."
    for bandFilter in bandFilters:
        print(bandNum)       
        # run 6S and get correction factors        
        s.wavelength = Wavelength(startWV, endWV, bandFilter)
        ################
        # print inputs to 6S:
#        print "atm profile: ", s.atmos_profile
#        print "aero profile: ", s.aero_profile
#        print "altitude: ", s.altitudes
#        print "atoms corr: ", s.atmos_corr
#        print "geometry: ", s.geometry
#        print "geometry.solar_z: ", s.geometry.solar_z
#        print "geometry.solar_a: ", s.geometry.solar_a
#        print "geometry.view_z: ", s.geometry.view_z
#        print "geometry.view_a: ", s.geometry.view_a
#        print "geometry.day: ", s.geometry.day
#        print "geometry.month: ", s.geometry.month
#        print "wavelength: ", s.wavelength
        ################
        s.run()
        xa = s.outputs.coef_xa # inverse of transmitance
        xb = s.outputs.coef_xb # scattering term of the atmosphere
        xc = s.outputs.coef_xc # reflectance of atmosphere for isotropic light (albedo)
        
        outputParams.append({'xa': xa, 'xb': xb, 'xc': xc})
        bandNum += 1
    # If the function is being used in the usual fashion, it returns outputParams, if it is being run in preparation for adjencency correction, the s-class is returned
    if adjCorrOutput == False:
        s = None
        return outputParams
    else:
        return s
    
def performAtmCorrection(inImg, correctionParams6S, adjCorr=False, s=None):
    refl = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount)) 
    pixelSize = inImg.GetGeoTransform()[1] # assume same horizontal and vertical resolution
    
    for bandNum, correctionParam in enumerate(correctionParams6S):
        # Read uncorrected radiometric data and correct
        radianceData = inImg.GetRasterBand(bandNum+1).ReadAsArray()

        # Interpolate the 6S correction parameters from one per image tile to
        # one per image pixel
        xa = explodeCorrectionParam(correctionParam['xa'], radianceData.shape)         
        xb = explodeCorrectionParam(correctionParam['xb'], radianceData.shape)
        xc = explodeCorrectionParam(correctionParam['xc'], radianceData.shape)                   
        
        # Perform the atmospheric correction        
        y = np.where(np.isnan(radianceData), np.nan, xa*radianceData - xb)
        refl[:,:,bandNum] = np.where(np.isnan(y), 0, np.maximum(y/(1.0+xc*y),0.0))
        
        # Perform adjecency correction if required
        if adjCorr:
            try:
                radius = float(adjCorr)
            except:
                radius = 300    
            refl[:,:,bandNum-1] = adjacencyCorrection(refl[:,:,bandNum-1], pixelSize, s, radius)
     
            
    res = bathyUtilities.saveImg (refl, inImg.GetGeoTransform(), inImg.GetProjection(), "MEM")
    refl = None
    radianceData = None    
    
    return res

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
def adjacencyCorrection(refl, pixelSize, s, radius = 1000.0):
    print("Adjacency correction")    
    # definition below eq (4)
    u_v = cos(radians(s.geometry.view_z))

    tau = s.outputs.optical_depth_total.total
    T_dir = s.outputs.transmittance_global_gas.upward
    T_dif = s.outputs.transmittance_total_scattering.upward  
    T= 1 - ((1-T_dif) + (1-T_dir))           
    
    # Fill in any NaN values, particularly at the edges of the image
    mask = np.isnan(refl)
    refl[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), refl[~mask])
    
    # Calculate average reflectance of adjacent pixels
    # The adjacency effect can come from pixels within 1km of the central pixel (Verhoef et al., 2003) so  sigma should be half of that in gaussian filter
    sigma = int(radius/pixelSize)
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
    refl[refl<0.0] = 0.0
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
    
    
    
def readGeometryWV2(metadataFile, model6S):

    s = model6S    
    
    # read viewing gemotery from WV2 metadata file
    meanSunElRegex = "\s*meanSunEl\s*=\s*(.*);"
    meanSunAzRegex = "\s*meanSunAz\s*=\s*(.*);"
    meanSatElRegex = "\s*meanSatEl\s*=\s*(.*);"
    meanSatAzRegex = "\s*meanSatAz\s*=\s*(.*);"
    # depending on the product type there can be either firstLineTime or earliestAcqTime in the metadata file
    firstLineTimeRegex = "\s*firstLineTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"   
    earliestAcqTimeRegex = "\s*earliestAcqTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"
    
    
    month = 0; day = 0; sunEl = 0.0; sunAz = 0.0; satEl = 0.0; satAz = 0.0;
    
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(firstLineTimeRegex, line)
            if not match:
                 match = re.match(earliestAcqTimeRegex, line)
            if match:
                 month = int(match.group(2))
                 day = int(match.group(3))    
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
            
if __name__ == "__main__":
    atmCorr6S_WV2("\\\\DKCPH1-STOR.DHI.DK\\Projects\\18800137\\Metohi_data\\053660695010_01\\053660695010_01_P001_MUL\\13SEP19094612-M2AS-053660695010_01_P001.IMD")
