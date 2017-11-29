import re


def get_geometry(mtdFile, extent=None):
    """

    extent is [minX, maxY, maxX, minY]
    """
    # read viewing gemotery from Landsat metadata file

    sun_elevation_regex = "\s*SUN_ELEVATION\s*=\s*(.*)\s*"
    sun_azimuth_regex = "\s*SUN_AZIMUTH\s*=\s*(.*)\s*"
    dateAcquired_regex = "\s*DATE_ACQUIRED\s*=\s*\d{4}-(\d{2})-(\d{2})\s*"
    minX_regex = "\s*CORNER_LL_PROJECTION_X_PRODUCT\s*=\s*(\d+\.\d+)\s*"
    maxX_regex = "\s*CORNER_UR_PROJECTION_X_PRODUCT\s*=\s*(\d+\.\d+)\s*"

    gdict = {}
    with open(mtdFile, 'r') as metadata:
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

    gdict['sensor_azimuth'] = None

    if extent is None:
        gdict['sensor_zenith'] = 0.0
    else:
        extentMiddleX = (extent[0] + extent[2]) / 2
        imageMiddleX = (minX + maxX) / 2
        # L8 has 15 deg field of view so VZA ranges between 0 and 7.5 degrees
        gdict['sensor_zenith'] = abs(imageMiddleX - extentMiddleX)/(imageMiddleX - minX) * 7.5
    return gdict
