import re
from math import cos, radians, pi

import numpy as np
from gdal_utils.gdal_utils import array_to_gtiff

from ..sensors import sensor_is


def toaReflectance(inImg, metadataFile, sensor):
    if sensor_is(sensor, 'WV'):
        res = toaReflectanceWV2(inImg, metadataFile)
    elif sensor_is(sensor, 'PHR'):
        res = toaReflectancePHR1(inImg, metadataFile)
    elif sensor_is(sensor, 'L7L8'):
        res = toaReflectanceL8(inImg, metadataFile)
    elif sensor_is(sensor, 'S2'):
        res = toaReflectanceS2(inImg, metadataFile)
    return res


def toaReflectanceWV2(inImg, metadataFile):
    """Estimate toa reflectance of radiometric WV2 data ignoric atmospheric, topographic and
       BRDF effects.

    Notes
    -------
    Based on http://www.digitalglobe.com/sites/default/files/Radiometric_Use_of_WorldView-2_Imagery%20%281%29.pdf
    Also works with GeoEye-1 and might work with other Digital Globe providers after a small modification
    """

    # Band averaged solar spectral irradiances at 1 AU Earth-Sun distance.
    # The first one is for panchromatic band.
    # For WV2 coming from Table 4 from the document in units of (W/m^2/μm/str).
    # GE01 irradiance is from https://apollomapping.com/wp-content/user_uploads/2011/09/GeoEye1_Radiance_at_Aperture.pdf
    # and is in units of (mW/cm^2/mum/str)
    ssi = {"WV02":[1580.8140, 1758.2229, 1974.2416, 1856.4104, 1738.4791, 1559.4555, 1342.0695, 1069.7302, 861.2866],
           "WV03":[1574.41, 1757.89, 2004.61, 1830.18, 1712.07, 1535.33, 1348.08, 1055.94, 858.77], # Thuillier 2003
           # "WV03":[1578.28, 1743.9, 1974.53, 1858.1, 1748.87, 1550.58, 1303.4, 1063.92, 858.632], # ChKur
           # "WV03":[1583.58, 1743.81, 1971.48, 1856.26, 1749.4, 1555.11, 1343.95, 1071.98, 863.296], # WRC
           "GE01":[161.7, 196.0, 185.3, 150.5, 103.9]}

    # depending on the product type there can be either firstLineTime or earliestAcqTime in the metadata file
    firstLineTimeRegex   = "\s*firstLineTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"
    earliestAcqTimeRegex = "\s*earliestAcqTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"
    meanSunElRegex       = "\s*meanSunEl\s*=\s*(.*);"
    satIdRegex           = "\s*satId\s*=\s*\"(.*)\";"

    # get year, month, day and time and sun zenith angle from the metadata file
    with open(metadataFile, 'r') as metadata:
        for line in metadata:
            match = re.match(firstLineTimeRegex, line)
            if not match:
                match = re.match(earliestAcqTimeRegex, line)
            if match:
                year =  int(match.group(1))
                month = int(match.group(2))
                day =   int(match.group(3))
                UT =    float(match.group(4)) + float(match.group(5))/60 + float(match.group(6))/3600
            match = re.match(meanSunElRegex, line)
            if match:
                sza = radians(90-float(match.group(1)))
            match = re.match(satIdRegex, line)
            if match:
                ssi = ssi[match.group(1)]

    # get actual Earth-Sun distance follwoing equations from Radiometric Use Of WorldView-2 Imagery - Technical note
    if month < 3:
        year -= 1.0
        month += 12.0
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
    for band in range(bandNum):
        print(band + 1)
        radData = inImg.GetRasterBand(band + 1).ReadAsArray()
        reflectanceData[:, :, band] = np.where(np.isnan(radData), np.nan,
                                               (radData * des ** 2 * pi) / (ssi[band + 1] * cos(sza)))

    return array_to_gtiff(reflectanceData, "MEM", inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)


def toaReflectancePHR1(inImg, metadataFile):
    # for now just do nothing
    return inImg


def toaReflectanceL8(inImg, metadataFile):
    # todo Implement
    return inImg


def toaReflectanceS2(inImg, metadataFile):
    """Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    rc = float(metadataFile['reflection_conversion'])
    # Convert to TOA reflectance
    print("TOA reflectance")
    rToa = np.zeros((inImg.RasterYSize, inImg.RasterXSize, inImg.RasterCount))
    for i in range(inImg.RasterCount):
        print(i)
        rToa[:, :, i] = inImg.GetRasterBand(i + 1).ReadAsArray().astype(float) / rc

    return array_to_gtiff(rToa, "MEM", inImg.GetProjection(), inImg.GetGeoTransform(), banddim=2)
