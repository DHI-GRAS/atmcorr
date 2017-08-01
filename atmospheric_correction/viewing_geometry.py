import re
from xml.etree import ElementTree as ET

from atmospheric_correction.sensors import sensor_is

_required = {'sun_zenith', 'sun_azimuth', 'sensor_zenith', 'sensor_azimuth', 'month', 'day'}


def _check_gdict(d):
    dkeys = set(d)
    if not dkeys >= _required:
        raise ValueError(
                'Geometry dict is lacking the following keys: {}'.format(_required - dkeys))


def get_geometry(sensor, mtdfile):
    if sensor_is(sensor, 'WV'):
        gdict = get_geometry_WV2(mtdfile)
    elif sensor_is(sensor, 'PHR'):
        gdict = get_geometry_PHR1(mtdfile)
    elif sensor_is(sensor, 'L7L8'):
        gdict = get_geometry_L8(mtdfile)
    elif sensor_is(sensor, 'S2'):
        gdict = get_geometry_S2(mtd_dict=mtdfile)
    _check_gdict(gdict)
    return gdict


def get_geometry_PHR1(mtdfile):

    gdict = {}

    tree = ET.parse(mtdfile)

    # get down to the appropirate node
    root = tree.getroot()
    Geometric_Data = root.findall('Geometric_Data')[0]
    Use_Area = Geometric_Data.findall('Use_Area')[0]
    for LGV in Use_Area.findall('Located_Geometric_Values'):
        # get angles for centre of the image
        if LGV.findall('LOCATION_TYPE')[0].text == "Center":
            Acquisition_Angles = LGV.findall('Acquisition_Angles')[0]
            gdict['sensor_azimuth'] = float(Acquisition_Angles.findall('AZIMUTH_ANGLE')[0].text)
            gdict['sensor_zenith'] = float(Acquisition_Angles.findall('INCIDENCE_ANGLE')[0].text)
            Solar_Incidences = LGV.findall('Solar_Incidences')[0]
            gdict['sun_azimuth'] = float(Solar_Incidences.findall('SUN_AZIMUTH')[0].text)
            sun_elevation = float(Solar_Incidences.findall('SUN_ELEVATION')[0].text)
            gdict['sun_zenith'] = 90.0 - sun_elevation

            # get month and day
            timeStr = LGV.findall('TIME')[0].text
            date_regex = '\d{4}-(\d{2})-(\d{2})T.*'
            match = re.match(date_regex, timeStr)
            try:
                gdict['month'] = int(match.group(1))
                gdict['day'] = int(match.group(2))
            except AttributeError:
                raise ValueError('Unable to get month and day')
            break

    return gdict


def get_geometry_WV2(mtdfile):
    # read viewing gemotery from WV2 metadata file
    meanSunEl_regex = "\s*meanSunEl\s*=\s*(.*);"
    meanSunAz_regex = "\s*meanSunAz\s*=\s*(.*);"
    meanSatEl_regex = "\s*meanSatEl\s*=\s*(.*);"
    meanSatAz_regex = "\s*meanSatAz\s*=\s*(.*);"
    # depending on the product type there can be either
    # firstLineTime or earliestAcqTime in the metadata file
    firstLineTime_regex = (
            "\s*firstLineTime\s*=\s*(\d{4})[-_](\d{2})[-_]"
            "(\d{2})T(\d{2}):(\d{2}):(.*)Z;")
    earliestAcqTime_regex = (
            "\s*earliestAcqTime\s*=\s*(\d{4})[-_](\d{2})[-_]"
            "(\d{2})T(\d{2}):(\d{2}):(.*)Z;")

    gdict = {}
    with open(mtdfile, 'r') as metadata:
        for line in metadata:
            match = re.match(firstLineTime_regex, line)
            if not match:
                match = re.match(earliestAcqTime_regex, line)
            if match:
                gdict['month'] = int(match.group(2))
                gdict['day'] = int(match.group(3))
                continue
            else:
                raise ValueError('This did not work so well.')

            match = re.match(meanSunEl_regex, line)
            if match:
                sun_elevation = float(match.group(1))
                gdict['sun_zenith'] = 90.0 - sun_elevation
                continue

            match = re.match(meanSunAz_regex, line)
            if match:
                gdict['sun_azimuth'] = float(match.group(1))
                continue

            match = re.match(meanSatEl_regex, line)
            if match:
                sensor_elevation = float(match.group(1))
                gdict['sensor_zenith'] = 90.0 - sensor_elevation
                continue

            match = re.match(meanSatAz_regex, line)
            if match:
                gdict['sensor_azimuth'] = float(match.group(1))
                continue
    _check_gdict(gdict)
    return gdict


def get_geometry_L8(mtdfile, extent):
    # read viewing gemotery from Landsat metadata file

    sun_elevation_regex = "\s*SUN_ELEVATION\s*=\s*(.*)\s*"
    sun_azimuth_regex = "\s*SUN_AZIMUTH\s*=\s*(.*)\s*"
    dateAcquired_regex = "\s*DATE_ACQUIRED\s*=\s*\d{4}-(\d{2})-(\d{2})\s*"
    minX_regex = "\s*CORNER_LL_PROJECTION_X_PRODUCT\s*=\s*(\d+\.\d+)\s*"
    maxX_regex = "\s*CORNER_UR_PROJECTION_X_PRODUCT\s*=\s*(\d+\.\d+)\s*"

    gdict = {}
    with open(mtdfile, 'r') as metadata:
        for line in metadata:
            match = re.match(dateAcquired_regex, line)
            if match:
                gdict['month'] = int(match.group(1))
                gdict['day'] = int(match.group(2))
            match = re.match(sun_elevation_regex, line)
            if match:
                sun_elevation = float(match.group(1))
                gdict['sun_zenith'] = 90 - sun_elevation
            match = re.match(sun_azimuth_regex, line)
            if match:
                gdict['sun_azimuth'] = float(match.group(1))
            match = re.match(minX_regex, line)
            if match:
                minX = float(match.group(1))
            match = re.match(maxX_regex, line)
            if match:
                maxX = float(match.group(1))

    if extent is None:
        gdict['sensor_zenith'] = 0.0
    else:
        # extent is [minX, maxY, maxX, minY]
        extentMiddleX = (extent[0]+extent[2])/2
        imageMiddleX = (minX + maxX)/2
        # L8 has 15 deg field of view so VZA ranges between 0 and 7.5 degrees
        gdict['sensor_zenith'] = abs(imageMiddleX - extentMiddleX)/(imageMiddleX - minX) * 7.5
    _check_gdict(gdict)
    return gdict


def get_geometry_S2(mtd_dict):
    """Get geometry dictionaty for S2

    Parameters
    ----------
    mtd_dict : dict
        meta data dict from `satmeta`
    """
    granule = mtd_dict['current_granule']
    granule_meta = mtd_dict[granule]
    gdict = {}
    # Get the granule metadata from metadata dictionary
    copy_keys = ['sun_zenith', 'sun_azimuth', 'sensor_zenith', 'sensor_azimuth']
    for key in copy_keys:
        gdict[key] = granule_meta[key]
    d = mtd_dict['sensing_time']
    gdict['month'] = d.month
    gdict['day'] = d.day
    _check_gdict(gdict)
    return gdict
