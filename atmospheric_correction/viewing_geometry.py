import re
from xml.etree import ElementTree as ET

from Py6S import Geometry


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
                day = int(match.group(3))
            else:
                raise ValueError('This did not work so well.')

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
