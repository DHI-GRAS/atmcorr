from xml.etree import ElementTree as ET

import numpy as np
import dateutil.parser


def get_date(mtdFile):
    tree = ET.parse(mtdFile)
    # get down to the appropirate node
    root = tree.getroot()
    Geometric_Data = root.findall('Geometric_Data')[0]
    Use_Area = Geometric_Data.findall('Use_Area')[0]
    for LGV in Use_Area.findall('Located_Geometric_Values'):
        # get angles for centre of the image
        if LGV.findall('LOCATION_TYPE')[0].text == "Center":
            # get year month and day
            datestr = LGV.findall('TIME')[0].text
            return dateutil.parser.parse(datestr)
    raise ValueError('Unable to get date from file \'{}\''.format(mtdFile))


def get_gain_bias(mtdfile):

    # get correction factors
    gain = np.zeros(4)
    bias = np.zeros(4)

    # In the XML file the band order is specified as BGRN.
    # However in reality it is RGBN. Therefore mapping is required
    # between the two
    bandMapping = {'B2': 0, 'B1': 1, 'B0': 2, 'B3': 3}

    # get down to the appropirate node
    tree = ET.parse(mtdfile)
    root = tree.getroot()
    Radiometric_Data = root.findall('Radiometric_Data')[0]
    Radiometric_Calibration = Radiometric_Data.findall('Radiometric_Calibration')[0]
    Instrument_Calibration = Radiometric_Calibration.findall('Instrument_Calibration')[0]
    Band_Measurement_List = Instrument_Calibration.findall('Band_Measurement_List')[0]
    for Band_Radiance in Band_Measurement_List.findall('Band_Radiance'):
        band = Band_Radiance.findall('BAND_ID')[0].text
        gain[bandMapping[band]] = float(Band_Radiance.findall('GAIN')[0].text)
        bias[bandMapping[band]] = float(Band_Radiance.findall('BIAS')[0].text)
    return gain, bias
