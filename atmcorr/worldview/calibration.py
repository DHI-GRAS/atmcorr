import re
import datetime

import numpy as np

GAIN = {
        'WV2': np.array([
            0.8632, 1.001, 1.0436, 1.0305, 1.0249, 0.9779, 0.981, 0.9217, 1.0264]),
        'WV3': np.array([
            1.157, 1.07, 1.082, 1.048, 1.028, 0.979, 1.006,
            0.975, 1.045])}  # as per calibration from DG released 3/6/2015

BIAS = {
        'WV2': np.array([
            6.6863, 2.399, 0.3973, 0.7744, -0.1495, 2.0383, 1.859, 2.0357, 6.5783]),
        'WV3': np.array([
            7.07, 4.253, 2.633, 2.074, 1.807, 2.633, 3.406,
            2.258, 2.22])}  # as per calibration from DG relased 3/6/2015

SSI = {
        'WV02': [
            1758.2229, 1974.2416,
            1856.4104, 1738.4791, 1559.4555,
            1342.0695, 1069.7302, 861.2866,
            1580.8140],
        'WV03': [
            1757.89, 2004.61,
            1830.18, 1712.07, 1535.33,
            1348.08, 1055.94, 858.77,
            1574.41],  # Thuillier 2003
        'WV03_ChKur': [
             1743.9, 1974.53,
             1858.1, 1748.87, 1550.58,
             1303.4, 1063.92, 858.632,
             1578.28],  # ChKur
        'WV03_WRC': [
            1743.81, 1971.48,
            1856.26, 1749.4, 1555.11,
            1343.95, 1071.98, 863.296,
            1583.58],  # WRC
        'GE01': [
            196.0, 185.3, 150.5, 103.9, 161.7]}


def get_effectivebw_abscalfactor_WV(mtdfile):
    # get the correction factors from the metadata file, assuming the number and
    # order of bands is the same in the image and the metadata file
    abscalfactor_rgx = "\s*absCalFactor\s*=\s*(.*);"
    effectivebw_rgx = "\s*effectiveBandwidth\s*=\s*(.*);"

    effectivebw = []
    abscalfactor = []
    with open(mtdfile) as mf:
        for line in mf:
            match = re.match(abscalfactor_rgx, line)
            if match:
                abscalfactor.append(float(match.group(1)))
            match = re.match(effectivebw_rgx, line)
            if match:
                effectivebw.append(float(match.group(1)))
    if not effectivebw:
        raise ValueError('Unable to get effective bandwidth from mtdfile.')
    if not abscalfactor:
        raise ValueError('Unable to get abs cal factor from mtdfile.')
    return np.array(effectivebw), np.array(abscalfactor)


def _get_meta(mtdfile):

    # depending on the product type there can be either
    # firstLineTime or earliestAcqTime in the metadata file
    firstLineTimeRegex = "\s*firstLineTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"
    earliestAcqTimeRegex = "\s*earliestAcqTime\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})T(\d{2}):(\d{2}):(.*)Z;"
    meanSunElRegex = "\s*meanSunEl\s*=\s*(.*);"
    satIdRegex = "\s*satId\s*=\s*\"(.*)\";"

    satID = None
    sza = None
    date = None
    hours = None
    # get year, month, day and time and sun zenith angle from the metadata file
    with open(mtdfile, 'r') as metadata:
        for line in metadata:
            match = re.match(firstLineTimeRegex, line)
            if not match:
                match = re.match(earliestAcqTimeRegex, line)
            if match:
                year = int(match.group(1))
                month = int(match.group(2))
                day = int(match.group(3))
                hours = (
                        float(match.group(4)) +
                        float(match.group(5))/60 +
                        float(match.group(6))/3600)
                date = datetime.datetime(year, month, day)
            match = re.match(meanSunElRegex, line)
            if match:
                sun_el = float(match.group(1))
                sza = np.radians(90 - sun_el)
            match = re.match(satIdRegex, line)
            if match:
                satID = match.group(1)

    if None in (satID, sza, date, hours):
        raise ValueError('This did not work, as expected.')

    return satID, sza, date, hours


def get_earth_sun_distance(mtdfile):

    # Band averaged solar spectral irradiances at 1 AU Earth-Sun distance.
    # The last one is for panchromatic band.
    # For WV2 coming from Table 4 from the document in units of (W/m^2/\mum/str).
    # GE01 irradiance is from
    # https://apollomapping.com/wp-content/user_uploads/2011/09/GeoEye1_Radiance_at_Aperture.pdf
    # and is in units of (mW/cm^2/mum/str)
    satID, sza, date, hours = _get_meta(mtdfile)

    ssi = SSI[satID]

    # get actual Earth-Sun distance follwoing equations from
    # Radiometric Use Of WorldView-2 Imagery - Technical note
    year = date.year
    month = date.month
    day = date.day
    if month < 3:
        year -= 1.0
        month += 12.0
    A = int(year/100.0)
    B = 2.0-A+int(A/4.0)
    JD = (
            int(365.25 * (year + 4716.0)) +
            int(30.6001 * (month + 1)) +
            day + hours / 24.0 + B - 1524.5)
    D = JD - 2451545.0
    g = 357.529 + 0.98560025 * D
    des = 1.00014 - 0.01671 * np.cos(np.radians(g)) - 0.00014 * np.cos(np.radians(2 * g))

    return des, ssi, sza
