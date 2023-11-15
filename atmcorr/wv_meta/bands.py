BANDS_MULT = ['COASTAL', 'BLUE', 'GREEN', 'YELLOW', 'RED', 'REDEDGE', 'NIR1', 'NIR2']

NBANDS_MULT = {
    'GE01': 4,
    'WV02': 8,
    'WV03': 8,
    'WV04': 4,
    'QUICKBIRD': 4,
    'WV01': 0,
    'IKONOS': 4}


def get_values_sorted(bandsdict, sat_id):
    """Get values of multispectral bands in bandsdict"""
    values = [bandsdict[band] for band in BANDS_MULT if band in bandsdict]
    n_expected = NBANDS_MULT[sat_id]
    if len(values) != n_expected:
        raise ValueError(
            'Expecting {} multispectral values for sat ID \'{}\'. Found {}.'
            .format(n_expected, sat_id, len(values)))
    return values
