from ..metadata.worldview import parse_metadata


def get_geometry(mtdFile):
    mtd = parse_metadata(mtdFile)
    angles = mtd['angles']
    gdict = {}
    gdict['day'] = mtd['sensing_time'].day
    gdict['month'] = mtd['sensing_time'].month
    gdict['sun_zenith'] = 90.0 - angles['meanSunEl']
    gdict['sun_azimuth'] = angles['meanSunAz']
    gdict['sensor_zenith'] = 90.0 - angles['meanSatEl']
    gdict['sensor_azimuth'] = angles['meanSatAz']
    return gdict
