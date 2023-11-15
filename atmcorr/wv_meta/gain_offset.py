from . import bands

SOURCE = (
    'https://dg-cms-uploads-production.s3.amazonaws.com/'
    'uploads/document/file/209/ABSRADCAL_FLEET_2016v0_Rel20170606.pdf')

GAIN = {
    'WV03': {  # 2016v0.Int
        'PAN': 0.950,
        'COASTAL': 0.905,
        'BLUE': 0.940,
        'GREEN': 0.938,
        'YELLOW': 0.962,
        'RED': 0.964,
        'REDEDGE': 1.000,
        'NIR1': 0.961,
        'NIR2': 0.978,
        'SIWR1': 1.200,
        'SIWR2': 1.227,
        'SIWR3': 1.199,
        'SIWR4': 1.196,
        'SIWR5': 1.262,
        'SIWR6': 1.314,
        'SIWR7': 1.346,
        'SIWR8': 1.376},
    'WV02': {  # 2016v0.Int
        'PAN': 0.942,
        'COASTAL': 1.151,
        'BLUE': 0.988,
        'GREEN': 0.936,
        'YELLOW': 0.949,
        'RED': 0.952,
        'REDEDGE': 0.974,
        'NIR1': 0.961,
        'NIR2': 1.002},
    'WV04': {  # Awaiting calibration results from DG RSS team
        'PAN': 1.0,
        'BLUE': 1.0,
        'GREEN': 1.0,
        'RED': 1.0,
        'NIR1': 1.0},
    'GE01': {  # 2016v3.Int
        'PAN': 0.970,
        'BLUE': 1.053,
        'GREEN': 0.994,
        'RED': 0.998,
        'NIR1': 0.994},
    'QUICKBIRD': {  # 2016v0.Int
        'PAN': 0.870,
        'BLUE': 1.105,
        'GREEN': 1.071,
        'RED': 1.060,
        'NIR1': 1.020},
    'WV01': {  # 2016v0.Int
        'PAN': 1.016},
    'IKONOS': {  # 2014v3
        'PAN': 0.907,
        'BLUE': 1.073,
        'GREEN': 0.990,
        'RED': 0.940,
        'NIR1': 1.043}}

OFFSET = {
    'WV03': {  # 2016v0.Int
        'PAN': -3.629,
        'COASTAL': -8.604,
        'BLUE': -5.809,
        'GREEN': -4.996,
        'YELLOW': -3.649,
        'RED': -3.021,
        'REDEDGE': -4.521,
        'NIR1': -5.522,
        'NIR2': -2.992,
        'SIWR1': -5.546,
        'SIWR2': -2.600,
        'SIWR3': -2.309,
        'SIWR4': -1.676,
        'SIWR5': -0.705,
        'SIWR6': -0.669,
        'SIWR7': -0.512,
        'SIWR8': -0.512},
    'WV02': {  # 2016v0.Int
        'PAN': -2.704,
        'COASTAL': -7.478,
        'BLUE': -5.736,
        'GREEN': -3.546,
        'YELLOW': -3.564,
        'RED': -2.512,
        'REDEDGE': -4.120,
        'NIR1': -3.300,
        'NIR2': -2.891},
    'WV04': {  # Awaiting calibration results from DG RSS team
        'PAN': 0.0,
        'BLUE': 0.0,
        'GREEN': 0.0,
        'RED': 0.0,
        'NIR1': 0.0},
    'GE01': {  # 2016v3.Int
        'PAN': -1.926,
        'BLUE': -4.537,
        'GREEN': -4.175,
        'RED': -3.754,
        'NIR1': -3.870},
    'QUICKBIRD': {  # 2016v0.Int
        'PAN': -1.491,
        'BLUE': -2.820,
        'GREEN': -3.338,
        'RED': -2.954,
        'NIR1': -4.722},
    'WV01': {  # 2016v0.Int
        'PAN': -1.824},
    'IKONOS': {  # 2014v3
        'PAN': -4.461,
        'BLUE': -9.699,
        'GREEN': -7.937,
        'RED': -4.767,
        'NIR1': -8.869}}


def get_gain_values(sat_id):
    """Returns values for multispectral bands (see bands.BANDS_MULT)"""
    try:
        return bands.get_values_sorted(GAIN[sat_id], sat_id=sat_id)
    except KeyError:
        raise ValueError(
            'Satellite ID \'{}\' not supported. Choose from {}.'
            .format(sat_id, list(GAIN)))


def get_offset_values(sat_id):
    """Returns values for multispectral bands (see bands.BANDS_MULT)"""
    try:
        return bands.get_values_sorted(OFFSET[sat_id], sat_id=sat_id)
    except KeyError:
        raise ValueError(
            'Satellite ID \'{}\' not supported. Choose from {}.'
            .format(sat_id, list(OFFSET)))
