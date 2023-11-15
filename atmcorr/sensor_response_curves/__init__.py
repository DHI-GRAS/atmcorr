from __future__ import division
import os

import numpy as np


CSVDIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), 'data'
)

SENSOR_GROUPS = {
    'WV2': 'WV',
    'WV3': 'WV',
    'WV4': 'WV_4band',
    'GE1': 'WV_4band',
    'PHR1A': 'PHR',
    'PHR1B': 'PHR',
    'SPOT6': 'PHR',
    'PNEO3': 'PNEO',
    'PNEO4': 'PNEO',
    'L7': 'L7',
    'L8': 'L8',
    'S2A': 'S2',
    'S2B': 'S2'
}

SUPPORTED_SENSORS = sorted(list(SENSOR_GROUPS))

BANDS_TO_COLS = {
    'S2': {
        'wavelength': 'SR_WL',
        'coastal': 'SR_AV_B1',
        'blue': 'SR_AV_B2',
        'green': 'SR_AV_B3',
        'red': 'SR_AV_B4',
        'rededge': 'SR_AV_B5',
        'rededge2': 'SR_AV_B6',
        'rededge3': 'SR_AV_B7',
        'nir1': 'SR_AV_B8',
        'nir2': 'SR_AV_B8A',
        'nir3': 'SR_AV_B9',
        'swir1': 'SR_AV_B10',
        'swir2': 'SR_AV_B11',
        'swir3': 'SR_AV_B12'},
    'WV': {
        'wavelength': 'Wavelength',
        'pan': 'Panchromatic',
        'coastal': 'Coastal',
        'blue': 'Blue',
        'green': 'Green',
        'yellow': 'Yellow',
        'red': 'Red',
        'rededge': 'Red_Edge',
        'nir1': 'NIR1',
        'nir2': 'NIR2',
        'swir1': 'SWIR1',
        'swir2': 'SWIR2',
        'swir3': 'SWIR3',
        'swir4': 'SWIR4',
        'swir5': 'SWIR5',
        'swir6': 'SWIR6',
        'swir7': 'SWIR7',
        'swir8': 'SWIR8',
    },
    'WV_4band': {
        'wavelength': 'Wavelength',
        'pan': 'PAN',
        'blue': 'BLUE',
        'green': 'GREEN',
        'red': 'RED',
        'nir1': 'NIR'
    },
    'PHR': {
        'wavelength': 'Wavelength',
        'green': 'B2Green',
        'blue': 'B1Blue',
        'red': 'B3Red',
        'nir1': 'B4NIR'
    },
    'PNEO': {
        'wavelength': 'Wavelength',
        'green': 'Green',
        'blue': 'Blue',
        'red': 'Red',
        'coastal': 'BlueCoastal',
        'pan' : 'Pan',
        'nir1' : 'NIR',
        'rededge': 'RedEdge'
    },
    'L8': {
        'wavelength': 'Wavelength',
        'coastal': 'L8B1Coast',
        'blue': 'L8B2Blue',
        'green': 'L8B3Green',
        'red': 'L8B4Red',
        'nir1': 'L8B5NIR',
        'pan': 'L8B8Pan'
    },
    'L7': {
        'wavelength': 'Wavelength',
        'blue': 'L7B1Blue',
        'green': 'L7B2Green',
        'red': 'L7B3Red',
        'nir1': 'L7B4NIR',
        'nir2': 'L7B5NIR',
        'pan': 'L7B8Pan'
    }
}

COLS_TO_BANDS = {sk: {v: k for k, v in BANDS_TO_COLS[sk].items()} for sk in BANDS_TO_COLS}

BAND_SEQUENCE = {
    'S2': [
        'coastal', 'blue', 'green', 'red',
        'rededge', 'rededge2', 'rededge3',
        'nir1', 'nir2', 'nir3',
        'swir1', 'swir2', 'swir3',
    ],
    'WV': [
        'coastal', 'blue', 'green', 'yellow',
        'red', 'rededge', 'nir1', 'nir2',
        'swir1', 'swir2', 'swir3', 'swir4',
        'swir5', 'swir6', 'swir7', 'swir8',
    ],
    'WV_4band': [
        'pan', 'blue', 'green', 'red', 'nir1',
    ],
    'PHR': [
        'red', 'blue', 'green', 'nir1',
    ],
    'PNEO': [
        'red', 'blue', 'green', 'coastal', 'nir1', 'rededge'
    ],
    'L7': [
        'blue', 'green', 'red', 'nir1',
        # 'swir1', 'thermal', 'swir2', 'pan'],  # not defined in csv
    ],
    'L8': [
        'coastal', 'blue', 'green', 'red', 'nir1', 'pan',
        # 'swir1', 'swir2', 'cirrus', 'thermal1', 'thermal2']}  # not defined in csv
    ]
}


def _check_supported_sensor(sensor):
    if sensor not in SUPPORTED_SENSORS:
        raise ValueError(
            'Sensor \'{}\' is not supported. Choose from {}.'
            ''.format(sensor, SUPPORTED_SENSORS))


def _parse_csv(infile):
    return np.genfromtxt(infile, delimiter=',', names=True, dtype='float', comments='#', encoding="utf_8_sig")


def _rename_fields_inplace(a, name_map):
    a.dtype.names = [name_map[name] for name in a.dtype.names if name in a.dtype.names]


def _rename_sensor_fields(data, sensor):
    sensorgroup = SENSOR_GROUPS[sensor]
    colmap = COLS_TO_BANDS[sensorgroup]
    _rename_fields_inplace(data, colmap)


def _get_default_bands(sensor):
    sensorgroup = SENSOR_GROUPS[sensor]
    bandkeys = BAND_SEQUENCE[sensorgroup]
    return bandkeys


def _get_csv_file(sensor):
    """Get CSV file for given sensor

    Parameters
    ----------
    sensor : str in SUPPORTED_SENSORS
        sensor name
    """
    _check_supported_sensor(sensor)
    fname = sensor + '.txt'
    csvfile = os.path.join(CSVDIR, fname)
    if not os.path.isfile(csvfile):
        raise ValueError(
            'CSV file missing from package: \'{}\'.'.format(csvfile))
    return csvfile


def get_data_raw(sensor):
    """Get sensor response curve data for sensor

    Parameters
    ----------
    sensor : str in SUPPORTED_SENSORS
        sensor name

    Returns
    -------
    pandas.DataFrame
    """
    _check_supported_sensor(sensor)
    infile = _get_csv_file(sensor)
    return _parse_csv(infile)


def get_data_standard_names(sensor):
    """Get sensor response curve data for sensor with renamed fields

    Parameters
    ----------
    sensor : str in SUPPORTED_SENSORS
        sensor name

    Returns
    -------
    pandas.DataFrame : data with standard names
        column order corresponds to band order
    """
    data = get_data_raw(sensor)
    _rename_sensor_fields(data, sensor)
    if 'wavelength' not in data.dtype.names:
        raise RuntimeError(
                'Renaming seems to have failed.')
    return data


def get_response_curves(
        sensor, pan_only=False, bandkeys=None, band_ids=None):
    """Read response curves for given sensor

    Parameters
    ----------
    sensor : str in SUPPORTED_SENSORS
        sensor name
    pan_only : bool
        whether to get pan band only
    bandkeys : list of str, optional
        list of bands to get curves for
        default: BAND_SEQUENCE for the
                 given sensor group
    band_ids : list of int, optional
        list of band numbers to get curves for
        default: as bandkeys

    Returns
    -------
    wavelength : ndarray shape(nvals)
        wavelength column
    rcurves : ndarray shape(nvals, nbands)
        sensor response curves for all (selected) bands
    """
    data = get_data_standard_names(sensor)
    wavelength = data['wavelength']
    # collect bands specific for each sensor and start and end wavelenghts
    if bandkeys is not None:
        pass
    elif band_ids is not None:
        all_bands = _get_default_bands(sensor)
        bandkeys = [all_bands[i] for i in band_ids]
    elif pan_only:
        bandkeys = ['pan']
    else:
        bandkeys = _get_default_bands(sensor)
    data_stacked = np.vstack([data[name] for name in bandkeys if name in data.dtype.names])
    return wavelength, data_stacked
