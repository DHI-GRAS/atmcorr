import re


def get_geometry(mtdFile):
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
    with open(mtdFile, 'r') as metadata:
        for line in metadata:
            match = re.match(firstLineTime_regex, line)
            if not match:
                match = re.match(earliestAcqTime_regex, line)
            if match:
                gdict['month'] = int(match.group(2))
                gdict['day'] = int(match.group(3))
                continue

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
    return gdict
