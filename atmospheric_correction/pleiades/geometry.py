import re
from xml.etree import ElementTree as ET


def get_geometry_PHR1(mtdFile):

    gdict = {}

    tree = ET.parse(mtdFile)

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
