import os
import csv

import numpy as np
from scipy.interpolate import interp1d

from atmospheric_correction.sensors import sensor_is

_here = os.path.abspath(os.path.dirname(__file__))
_csvdir = os.path.join(_here, 'data')


def read_band_filters(sensor, isPan):
    if sensor_is(sensor, 'S2'):
        fname = 'S2.txt'
    else:
        sensor + '.txt'
    csvfile = os.path.join(_csvdir, fname)

    return _read_csv(csvfile, sensor, isPan)


def _read_csv(csvFilename, sensor, isPan):
    """Read band filters from CSV file and assign them to a given sensor"""
    # create empty lists for all the possible bands
    pan = []; coastal = []; blue = []; green = []; yellow = []; red = []; rededge = []; nir1 = []; nir2 = []; wavelength = []
    # S2 has some extra bands
    rededge2 = []; rededge3 = []#; nir3= []; swir1 = []; swir2 = []; swir3 = []

    # read in the data from CSV file
    with open(csvFilename, 'r') as csvFile:
        reader = csv.DictReader(csvFile)
        for line in reader:
            if sensor == "WV2" or sensor == "WV3":
                wavelength.append(float(line["Wavelength"]))
                pan.append(float(line["Panchromatic"]))
                coastal.append(float(line["Coastal"]))
                blue.append(float(line["Blue"]))
                green.append(float(line["Green"]))
                yellow.append(float(line["Yellow"]))
                red.append(float(line["Red"]))
                rededge.append(float(line["Red Edge"]))
                nir1.append(float(line["NIR1"]))
                nir2.append(float(line["NIR2"]))
            elif sensor == "PHR1A" or sensor == "PHR1B" or sensor == "SPOT6":
                wavelength.append(float(line["Wavelength"]))
                blue.append(float(line["B1Blue"]))
                green.append(float(line["B2Green"]))
                red.append(float(line["B3Red"]))
                nir1.append(float(line["B4NIR"]))
            elif sensor == "L8":
                wavelength.append(float(line["Wavelength"]))
                coastal.append(float(line["L8B1Coast"]))
                blue.append(float(line["L8B2Blue"]))
                green.append(float(line["L8B3Green"]))
                red.append(float(line["L8B4Red"]))
                nir1.append(float(line["L8B5NIR"]))
                pan.append(float(line["L8B8Pan"]))
            elif sensor == "L7":
                wavelength.append(float(line["Wavelength"]))
                blue.append(float(line["L7B1Blue"]))
                green.append(float(line["L7B2Green"]))
                red.append(float(line["L7B3Red"]))
                nir1.append(float(line["L7B4NIR"]))
            elif sensor == "S2A_10m" or sensor == "S2A_60m":
                wavelength.append(float(line["SR_WL"]))
                coastal.append(float(line["SR_AV_B1"]))
                blue.append(float(line["SR_AV_B2"]))
                green.append(float(line["SR_AV_B3"]))
                red.append(float(line["SR_AV_B4"]))
                rededge.append(float(line["SR_AV_B5"]))
                rededge2.append(float(line["SR_AV_B6"]))
                rededge3.append(float(line["SR_AV_B7"]))
                nir1.append(float(line["SR_AV_B8"]))
                nir2.append(float(line["SR_AV_B8A"]))
#                nir3.append(float(line["SR_AV_B9"]))
#                swir1.append(float(line["SR_AV_B10"]))
#                swir2.append(float(line["SR_AV_B11"]))
#                swir3.append(float(line["SR_AV_B12"]))

    # collect bands specific for each sensor and start and end wavelenghts
    startWV = wavelength[0]
    endWV = wavelength[-1]

    if isPan:
        bandFilters = [pan]
        return startWV, endWV, bandFilters

    if sensor_is(sensor, 'WV'):
        bandFilters = [coastal, blue, green, yellow, red, rededge, nir1, nir2]
    elif sensor_is(sensor, 'PHR'):
        # the order of Pleadis bands is like below (RGBN),
        # not like indicated in the metadata file (BGRN)
        bandFilters = [red, green, blue, nir1]
    elif sensor_is(sensor, "L8"):
        bandFilters = [coastal, blue, green, red, nir1]
    elif sensor_is(sensor, "L7"):
        bandFilters = [blue, green, red, nir1]
    elif sensor_is(sensor, 'S2'):
        bandFilters = [
                coastal, blue, green, red,
                rededge, rededge2, rededge3,
                nir1, nir2]  # , nir3, swir1, swir2, swir3]

    return startWV, endWV, bandFilters


def resample_band_filters(bandFilter, startWV, endWV, resolution):
    """Resample the given band filter to specified spectral resolution"""
    x = np.linspace(startWV, endWV, len(bandFilter))
    f = interp1d(x, bandFilter, kind='slinear')
    xnew = np.linspace(startWV, endWV, (endWV-startWV)/resolution + 1)
    return f(xnew)
