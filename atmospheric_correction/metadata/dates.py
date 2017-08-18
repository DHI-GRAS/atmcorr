import re
import datetime
from xml.etree import ElementTree as ET

from . import s2 as s2meta
from atmospheric_correction.sensors import sensor_is


def get_date_l7l8(mtdFile):
    DATEACQUIREDRegex = "\s*DATE_ACQUIRED\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})"
    with open(mtdFile, 'r') as metadata:
        for line in metadata:
            match = re.match(DATEACQUIREDRegex, line)
            if match is None:
                continue
            year = int(match.group(1))
            month = int(match.group(2))
            day = int(match.group(3))
            return datetime.date(year, month, day)
    raise ValueError('Unable to get date from file \'{}\''.format(mtdFile))


def get_date_wv(mtdFile):
    firstLineTimeRegex = (
            "\s*firstLineTime\s*=\s*(\d{4})[-_](\d{2})"
            "[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;")
    earliestAcqTimeRegex = (
            "\s*earliestAcqTime\s*=\s*(\d{4})[-_](\d{2})"
            "[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;")
    with open(mtdFile, 'r') as metadata:
        for line in metadata:
            match = re.match(firstLineTimeRegex, line)
            if match is None:
                match = re.match(earliestAcqTimeRegex, line)
            if match is None:
                continue
            year = int(match.group(1))
            month = int(match.group(2))
            day = int(match.group(3))
            return datetime.date(year, month, day)
    raise ValueError('Unable to get date from file \'{}\''.format(mtdFile))


def get_date_phr(mtdFile):
    tree = ET.parse(mtdFile)
    # get down to the appropirate node
    root = tree.getroot()
    Geometric_Data = root.findall('Geometric_Data')[0]
    Use_Area = Geometric_Data.findall('Use_Area')[0]
    for LGV in Use_Area.findall('Located_Geometric_Values'):
        # get angles for centre of the image
        if LGV.findall('LOCATION_TYPE')[0].text == "Center":
            # get year month and day
            timeStr = LGV.findall('TIME')[0].text
            dateRegex = '(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})\.*'
            match = re.match(dateRegex, timeStr)
            year = int(match.group(1))
            month = int(match.group(2))
            day = int(match.group(3))
            return datetime.date(year, month, day)
    raise ValueError('Unable to get date from file \'{}\''.format(mtdFile))


def get_sensing_date(sensor, mtdFile):
    if sensor_is(sensor, 'WV'):
        return get_date_wv(mtdFile)
    elif sensor_is(sensor, 'L7L8'):
        return get_date_l7l8(mtdFile)
    elif sensor_is(sensor, "PHR"):
        return get_date_phr(mtdFile)
    elif sensor_is(sensor, 'S2'):
        metadata = s2meta.parse_mtdfile(mtdFile)
        return metadata['sensing_time']
    else:
        raise ValueError('Unknown sensor.')
