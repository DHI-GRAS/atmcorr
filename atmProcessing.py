# -*- coding: utf-8 -*-
"""
Created on Tue Sep 09 12:54:33 2014

@author: rmgu
"""
import numpy as np
import re
import os
from xml.etree import ElementTree as ET
from osgeo import gdal
from atmCorr6S import *
from bathyUtilities import *
from atmParametersMODIS import *
###############################################################################
def atmProcessingMain(options):
    # Set the band numbers to the appropriate sensor
    sensor = options['sensor']
    coastal, blue, green, yellow, red, redEdge, nir1, nir2 = getSensorBandNumber(sensor)

    # Commonly used filenames
    dnFile = options["dnFile"]
    metadataFile = options["metadataFile"] 
    reflectanceFile = options["reflectanceFile"]
    roiFile = options["roiFile"] if "roiFile" in options else ""
    
    # 
    atmCorrMethod = options["atmCorrMethod"]
    atm = options["atm"]
    isPan = options["isPan"]
    adjCorr = 0 # this parameter is related to diffuse light. According to Mikkel the parameter is not important and is therefore set to 0.
    aeroProfile = options["aeroProfile"]

    # DN -> Radiance -> Reflectance        
    if atmCorrMethod == "6S":
        doDOS = False
        inImg = openAndClipRaster(dnFile, roiFile)            
        radianceImg = toaRadiance(inImg, metadataFile, sensor, doDOS=doDOS, isPan=isPan)
        inImg = None
        # If atmospheric parameters needed by 6S are not specified then
        # donwload and use MODIS atmopsheric products
        if not (atm['AOT'] and atm['PWV'] and atm['ozone']):
            aot, pwv, ozone = downloadAtmParametersMODIS(dnFile, metadataFile, sensor)
            if not atm['AOT']:
                atm['AOT'] = aot[0]
            if not atm['PWV']:
                atm['PWV'] = pwv[0]
            if not atm['ozone']:
                atm['ozone'] = ozone[0]
        print("AOT: "+str(atm['AOT']))
        print("Water Vapour: "+str(atm['PWV']))
        print("Ozone: "+str(atm['ozone']))
        
        reflectanceImg = atmCorr6S(metadataFile, radianceImg, atm=atm, sensor=sensor, isPan=isPan, adjCorr=adjCorr, aeroProfile=aeroProfile)
        radianceImg = None
        
    elif atmCorrMethod in ["DOS", "TOA"]:
        if atmCorrMethod == "DOS":
            doDOS = True
        else:
            doDOS = False 
        inImg = openAndClipRaster(dnFile, roiFile)
        radianceImg = toaRadiance(inImg, metadataFile, sensor, doDOS=doDOS)
        inImg = None
        reflectanceImg = toaReflectance(radianceImg, metadataFile, sensor)    
        radianceImg = None
        
    elif atmCorrMethod == "RAD":
       doDOS = False 
       inImg = openAndClipRaster(dnFile, roiFile)
       radianceImg = toaRadiance(inImg, metadataFile, sensor, doDOS=doDOS)
       reflectanceImg = radianceImg
        
    outImg = reflectanceImg
    outFile = reflectanceFile
    
    reflectanceImg = None        
    print("Saving corrected image")    
    driver = gdal.GetDriverByName ( "GTiff" )
    savedImg = driver.CreateCopy(outFile, outImg , 0 , driverOptionsGTiff)
    savedImg = None
    outImg = None

################################################################################################

def toaRadiance(inImg, metadataFile, sensor, doDOS, isPan = False):
    if sensor == "WV2" or sensor == "WV3":
        res = toaRadianceWV(inImg, metadataFile, doDOS, isPan, sensor)
    elif sensor == "PHR1A" or sensor == "PHR1B" or sensor == "SPOT6":
        res = toaRadiancePHR1(inImg, metadataFile, doDOS) 
    elif sensor == "L8" or sensor == "L7":
        res = toaRadianceL8(inImg, metadataFile, doDOS, isPan, sensor)
    return res

def toaReflectance(inImg, metadataFile, sensor):
    if sensor == "WV2" or sensor == "WV3":
        res = toaReflectanceWV2(inImg, metadataFile)
    elif sensor == "PHR1A" or sensor == "PHR1B" or sensor == "SPOT6":
        res = toaReflectancePHR1(inImg, metadataFile) 
    elif sensor == "L8" or sensor == "L7":
        res = toaReflectanceL8(inImg, metadataFile)
    return res

# Apply radiometric correction to WV2 image, with DOS atmospheric correction optional.
# Based on http://www.digitalglobe.com/sites/default/files/Radiometric_Use_of_WorldView-2_Imagery%20%281%29.pdf.
# [doesn't anymore after addition of gain and bias] Also works with GeoEye-1 (and possibly other DigitalGlobe satellties) since they have the same metadata file structure.
def toaRadianceWV(inImg, metadataFile, doDOS, isPan, sensor):

    absCalFactorRegex = "\s*absCalFactor\s*=\s*(.*);"
    effectiveBandwidthRegex = "\s*effectiveBandwidth\s*=\s*(.*);" 
    
    # get the correction factors from the metadata file, assuming the number and 
    # order of bands is the same in the image and the metadata file
    absCalFactor = []
    effectiveBandwidth = []
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(absCalFactorRegex, line)
            if match:
                absCalFactor.append(float(match.group(1)))
            match = re.match(effectiveBandwidthRegex, line)
            if match:
                effectiveBandwidth.append(float(match.group(1))) 

    bandNum = inImg.RasterCount 
    if sensor == "WV3":
        # WV3 calibration factors
        if isPan:
            gain = [1.045]
            bias = [2.22]
        else:
            gain = [1.157, 1.07, 1.082, 1.048, 1.028, 0.979, 1.006, 0.975] #as per calibration from DG released 3/6/2015
            bias = [7.07, 4.253, 2.633, 2.074, 1.807, 2.633, 3.406, 2.258] #as per calibration from DG relased 3/6/2015
    else:
        # WV2 calibration factors 
        if isPan:
            gain = [1.0264]
            bias = [6.5783]
        else:
            gain = [0.8632, 1.001, 1.0436, 1.0305, 1.0249, 0.9779, 0.981, 0.9217]
            bias = [6.6863, 2.399, 0.3973, 0.7744, -0.1495, 2.0383, 1.859, 2.0357]
    
    # perform dark object substraction          
    if doDOS:
        dosDN = darkObjectSubstraction(inImg)
    else:
        dosDN = list(np.zeros(bandNum))            
    
    # apply the radiometric correction factors to input image
    print("Radiometric correction")
    radiometricData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))       
    validMask = np.zeros((inImg.RasterYSize, inImg.RasterXSize))    
    for band in range(1,bandNum+1):
        print(band)
        rawData = np.int_(inImg.GetRasterBand(band).ReadAsArray())
        #radiometricData[:,:,band-1] = np.where(np.logical_and((rawData-dosDN[band-1])>0, rawData != 65536),rawData-dosDN[band-1],0)*absCalFactor[band-1]/effectiveBandwidth[band-1]
        radiometricData[:,:,band-1] = np.where(np.logical_and((rawData-dosDN[band-1])>0, rawData != 65536),rawData-dosDN[band-1],0)*absCalFactor[band-1]/effectiveBandwidth[band-1]*(2-gain[band-1])-bias[band-1]
        validMask = validMask + radiometricData[:,:,band-1]
    
    # Mark the pixels which have all radiances of 0 as invalid    
    invalidMask = np.where(validMask > 0, False, True)    
    radiometricData[invalidMask,:] = np.nan    
    
    res = saveImg (radiometricData, inImg.GetGeoTransform(), inImg.GetProjection(), "MEM")
    return res
   
# Estimate toa reflectance of radiometric WV2 data ignoric atmospheric, topographic and 
# BRDF effects.
# Based on http://www.digitalglobe.com/sites/default/files/Radiometric_Use_of_WorldView-2_Imagery%20%281%29.pdf
# Also works with GeoEye-1 and might work with other Digital Globe providers after a small modification
def toaReflectanceWV2(inImg, metadataFile):
     
    # Band averaged solar spectral irradiances at 1 AU Earth-Sun distance. The first one is for panchromatic band.  
    # For WV2 coming from Table 4 from the document in units of (W/m^2/μm/str).
    # GE01 irradiance is from https://apollomapping.com/wp-content/user_uploads/2011/09/GeoEye1_Radiance_at_Aperture.pdf
    # and is in units of (mW/cm^2/μm/str)
    ssi = {"WV02":[1580.8140, 1758.2229, 1974.2416, 1856.4104, 1738.4791, 1559.4555, 1342.0695, 1069.7302, 861.2866],
           "WV03":[1574.41, 1757.89, 2004.61, 1830.18, 1712.07, 1535.33, 1348.08, 1055.94, 858.77], # Thuillier 2003
           # "WV03":[1578.28, 1743.9, 1974.53, 1858.1, 1748.87, 1550.58, 1303.4, 1063.92, 858.632], # ChKur
           # "WV03":[1583.58, 1743.81, 1971.48, 1856.26, 1749.4, 1555.11, 1343.95, 1071.98, 863.296], # WRC
           "GE01":[161.7, 196.0, 185.3, 150.5, 103.9]} 
        
    # depending on the product type there can be either firstLineTime or earliestAcqTime in the metadata file
    firstLineTimeRegex = "\s*firstLineTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"   
    earliestAcqTimeRegex = "\s*earliestAcqTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"      
    meanSunElRegex = "\s*meanSunEl\s*=\s*(.*);"
    satIdRegex = "\s*satId\s*=\s*\"(.*)\";"
     
    # get year, month, day and time and sun zenith angle from the metadata file
    with open(metadataFile, 'r') as metadata:
         for line in metadata:
             match = re.match(firstLineTimeRegex, line)
             if not match:
                 match = re.match(earliestAcqTimeRegex, line)
             if match:
                 year = int(match.group(1))
                 month = int(match.group(2))
                 day = int(match.group(3))
                 UT = float(match.group(4)) + float(match.group(5))/60 + float(match.group(6))/3600
             match = re.match(meanSunElRegex, line)
             if match:
                 sza = radians(90-float(match.group(1)))
             match = re.match(satIdRegex, line)
             if match:
                 ssi = ssi[match.group(1)]

    # get actual Earth-Sun distance follwoing equations from Radiometric Use Of WorldView-2 Imagery - Technical note
    if month < 3:
        year = year-1.0 
        month = month+12.0               
    A = int(year/100.0)
    B = 2.0-A+int(A/4.0)
    JD = int(365.25*(year+4716.0)) + int(30.6001*(month+1)) + day + UT/24.0 + B - 1524.5 
    D = JD - 2451545.0
    g = 357.529+0.98560025*D
    des = 1.00014 - 0.01671*cos(radians(g))-0.00014*cos(radians(2*g))
    
    # apply the radiometric correction factors to input image
    print("TOA reflectance")
    reflectanceData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))
    bandNum = inImg.RasterCount    
    for band in range(1,bandNum+1):
        print(band)
        radData = inImg.GetRasterBand(band).ReadAsArray()
        reflectanceData[:,:,band-1] = np.where(np.isnan(radData), np.nan, (radData*des**2*pi)/(ssi[band]*cos(sza)))

    res = saveImg (reflectanceData, inImg.GetGeoTransform(), inImg.GetProjection(), "MEM")
    return res    

# Apply radiometric correction to Pleadis image, with DOS atmospheric correction optional.
def toaRadiancePHR1(inImg, metadataFile, doDOS):
    
    # get correction factors
    gain = [0, 0, 0, 0]
    bias = [0, 0, 0, 0]
    
    # In the XML file the band order is specified as BGRN. However in reality it is RGBN. Therefore mapping is required 
    # between the two
    bandMapping = {'B2': 0, 'B1': 1, 'B0': 2, 'B3': 3}    
    
    # get down to the appropirate node
    tree = ET.parse(metadataFile)
    root = tree.getroot()
    Radiometric_Data = root.findall('Radiometric_Data')[0]
    Radiometric_Calibration = Radiometric_Data.findall('Radiometric_Calibration')[0]
    Instrument_Calibration = Radiometric_Calibration.findall('Instrument_Calibration')[0]
    Band_Measurement_List = Instrument_Calibration.findall('Band_Measurement_List')[0]
    for Band_Radiance in Band_Measurement_List.findall('Band_Radiance'):
        band = Band_Radiance.findall('BAND_ID')[0].text        
        gain[bandMapping[band]] = float(Band_Radiance.findall('GAIN')[0].text)
        bias[bandMapping[band]] = float(Band_Radiance.findall('BIAS')[0].text)
    
    bandNum = inImg.RasterCount 

    # perform dark object substraction          
    if doDOS:
        dosDN = darkObjectSubstraction(inImg)
    else:
        dosDN = list(np.zeros(bandNum))            
    
    # apply the radiometric correction factors to input image
    print("Radiometric correction")
    radiometricData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))
    validMask = np.zeros((inImg.RasterYSize, inImg.RasterXSize))
    for band in range(1,bandNum+1):
        print(band)
        rawData = np.int_(inImg.GetRasterBand(band).ReadAsArray())
        radiometricData[:,:,band-1] = np.where(np.logical_and((rawData-dosDN[band-1])>0, rawData != 65536),(rawData-dosDN[band-1])/gain[band-1] + bias[band-1],0)
        validMask = validMask + radiometricData[:,:,band-1]
    
    # Mark the pixels which have all radiances of 0 as invalid
    invalidMask = np.where(validMask > 0, False, True)
    radiometricData[invalidMask,:] = np.nan  
    res = saveImg (radiometricData, inImg.GetGeoTransform(), inImg.GetProjection(), "MEM")
    return res    

def toaReflectancePHR1(inImg, metadataFile):
    # for now just do nothing
    return inImg


def toaRadianceL8(inImg, metadataFile, doDOS, isPan, sensor):
    
    multFactorRegex = "\s*RADIANCE_MULT_BAND_\d\s*=\s*(.*)\s*"
    addFactorRegex = "\s*RADIANCE_ADD_BAND_\d\s*=\s*(.*)\s*"     
    
    if inImg.RasterCount==1:
        if sensor == "L8":
            # The first 5 bands in L8 are VIS/NIR    
            visNirBands = range(1,6)
        elif sensor == "L7":
            # The first 4 bands in L7 are VIS/NIR
            visNirBands = range(1,5)
        rawData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, len(visNirBands)))
        
        # Raw Landsat 8/7 data has each band in a separate image. Therefore first open images with all the required band data.
        imgDir = os.path.dirname(inImg.GetFileList()[0])
        for _, _, files in os.walk(imgDir):
            for name in sorted(files):
                match = re.search('(([A-Z]{2}\d).+)_B(\d+)\.TIF$', name)
                if match and int(match.group(3)) in visNirBands:
                    band = int(match.group(3))
                    print(band)
                    rawImg = gdal.Open(os.path.join(imgDir, name), gdal.GA_ReadOnly)
                    rawData[:, :, band-1] = rawImg.GetRasterBand(1).ReadAsArray()
                    rawData = np.int_(rawData)
                    rawImg = None
                    
    # Panchromatic should only be one band but this way the isPan option can also
    # be used to processed L8 images which are stacked in one file.
    else:
        rawData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))
        for i in range(inImg.RasterCount):
            rawData[:,:,i] = inImg.GetRasterBand(i+1).ReadAsArray()
                
    # get the correction factors from the metadata file, assuming the number and 
    # order of bands is the same in the image and the metadata file
    multFactor = []
    addFactor = []
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(multFactorRegex, line)
            if match:
                multFactor.append(float(match.group(1)))
            match = re.match(addFactorRegex, line)
            if match:
                addFactor.append(float(match.group(1))) 
 
    # perform dark object substraction          
    if doDOS:
        rawImg  = saveImg (rawData, inImg.GetGeoTransform(), inImg.GetProjection(), "MEM")
        dosDN = darkObjectSubstraction(rawImg)
        rawImg = None
    else:
        dosDN = list(np.zeros(rawData.shape[2])) 
 
    # apply the radiometric correction factors to input image
    print("Radiometric correction")
    radiometricData = np.zeros((inImg.RasterYSize, inImg.RasterXSize, rawData.shape[2]))
    validMask = np.zeros((inImg.RasterYSize, inImg.RasterXSize))       
    for band in range(1, rawData.shape[2]+1):
        print(band)
        radiometricData[:,:,band-1] = np.where((rawData[:,:,band-1]-dosDN[band-1])>0,(rawData[:,:,band-1]-dosDN[band-1])*multFactor[band-1] + addFactor[band-1], 0)  
        validMask = validMask + radiometricData[:,:,band-1]
    
    # Mark the pixels which have all radiances of 0 as invalid    
    invalidMask = np.where(validMask > 0, False, True)    
    radiometricData[invalidMask,:] = np.nan  
            
    res = saveImg (radiometricData, inImg.GetGeoTransform(), inImg.GetProjection(), "MEM")
    return res
    
    
def toaReflectanceL8(inImg, metadataFile):
    # for now just do nothing
    return inImg    

def darkObjectSubstraction(inImg):
    print("DOS correction")
    dosDN = []
    tempData = inImg.GetRasterBand(1).ReadAsArray()
    numElements = np.size(tempData[tempData!=0])
    tempData = None
    for band in range(1,inImg.RasterCount+1):
        # use histogram with 2048 bins since WV2 has 11 bit radiometric resolution
        hist, edges = np.histogram(inImg.GetRasterBand(band).ReadAsArray(), bins=2048, range=(1,2048), density = False)
        for i in range(1, len(hist)):
            if hist[i] - hist[i-1] > (numElements-numElements*0.999999):
                dosDN.append(i-1)
                break
    return dosDN

###############################################################################
