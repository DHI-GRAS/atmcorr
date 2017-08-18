import re

import numpy as np


def get_gain_bias_WV(sensor, isPan):
    """Get WV calibration factors"""
    if sensor == "WV3":
        # WV3 calibration factors
        if isPan:
            gain = [1.045]
            bias = [2.22]
        else:
            gain = [1.157, 1.07, 1.082, 1.048, 1.028, 0.979, 1.006,
                    0.975]  # as per calibration from DG released 3/6/2015
            bias = [7.07, 4.253, 2.633, 2.074, 1.807, 2.633, 3.406,
                    2.258]  # as per calibration from DG relased 3/6/2015
    else:
        # WV2 calibration factors
        if isPan:
            gain = [1.0264]
            bias = [6.5783]
        else:
            gain = [0.8632, 1.001, 1.0436, 1.0305, 1.0249, 0.9779, 0.981, 0.9217]
            bias = [6.6863, 2.399, 0.3973, 0.7744, -0.1495, 2.0383, 1.859, 2.0357]
    return np.array(gain), np.array(bias)


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
