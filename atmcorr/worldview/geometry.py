from dg_calibration import metadata


def get_geometry(mtdFile):
    mtd = metadata.parse_metadata(mtdFile)
    angles = mtd['angles']
    gdict = {}
    gdict['day'] = mtd['sensing_date'].day
    gdict['month'] = mtd['sensing_date'].month
    gdict['sun_zenith'] = 90.0 - angles['meanSunEl']
    gdict['sun_azimuth'] = angles['meanSunAz']
    gdict['sensor_zenith'] = 90.0 - angles['meanSatEl']
    gdict['senzor_azimuth'] = angles['meanSatAz']
    return gdict
