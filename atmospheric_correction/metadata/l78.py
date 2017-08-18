import re

import numpy as np


def get_correction_factors(mtdfile):
    # get the correction factors from the metadata file, assuming the number and
    # order of bands is the same in the image and the metadata file
    mult_factorRegex = "\s*RADIANCE_MULT_BAND_\d\s*=\s*(.*)\s*"
    add_factorRegex = "\s*RADIANCE_ADD_BAND_\d\s*=\s*(.*)\s*"

    mult_factor = []
    add_factor = []
    with open(mtdfile, 'r') as metadata:
        for line in metadata:
            match = re.match(mult_factorRegex, line)
            if match:
                mult_factor.append(float(match.group(1)))
            match = re.match(add_factorRegex, line)
            if match:
                add_factor.append(float(match.group(1)))
    return np.array(mult_factor), np.array(add_factor)


def get_visnirbands(sensor):
    if sensor == "L8":
        # The first 5 bands in L8 are VIS/NIR
        return list(range(1, 6))
    elif sensor == "L7":
        # The first 4 bands in L7 are VIS/NIR
        return list(range(1, 5))
