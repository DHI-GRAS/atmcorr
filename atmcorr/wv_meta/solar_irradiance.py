from . import bands


CSV = {
    'WV02': """
PAN;1571.36;1575.38;1580.76
COASTAL;1773.81;1759.24;1757.77
BLUE;2007.27;1977.4;1974.29
GREEN;1829.62;1857.89;1856.03
YELLOW;1701.85;1738.11;1738.59
RED;1538.85;1554.95;1559.35
REDEDGE;1346.09;1302.19;1342.05
NIR1;1053.21;1061.4;1069.59
NIR2;856.599;856.816;861.201
    """,
    'WV03': """
PAN;1574.41;1578.28;1583.58
COASTAL;1757.89;1743.9;1743.81
BLUE;2004.61;1974.53;1971.48
GREEN;1830.18;1858.1;1856.26
YELLOW;1712.07;1748.87;1749.4
RED;1535.33;1550.58;1555.11
REDEDGE;1348.08;1303.4;1343.95
NIR1;1055.94;1063.92;1071.98
NIR2;858.77;858.632;863.296
SWIR1;479.019;478.873;494.595
SWIR2;263.797;257.55;261.494
SWIR3;225.283;221.448;230.518
SWIR4;197.552;191.583;196.766
SWIR 5;90.4178;86.5651;80.365
SWIR6;85.0642;82.0035;74.7211
SWIR7;76.9507;74.7411;69.043
SWIR8;68.0988;66.3906;59.8224
    """,
    'GE01': """
PAN;1610.73;1614.88;1619.49
BLUE;1993.18;1966.03;1963.53
GREEN;1828.83;1857.12;1855.25
RED;1491.49;1500.38;1506.29
NIR1;1022.58;1029.61;1037.7
    """,
    'QUICKBIRD': """
PAN;1370.92;1376.3;1381.72
BLUE;1949.59;1926.55;1924.62
GREEN;1823.64;1844.26;1842.81
RED;1553.78;1571.58;1574.65
NIR1;1102.85;1107.47;1113.72
    """,
    'WV01': """
PAN;1478.62;1481.48;1487.92
    """,
    'WV04': """
PAN;1608.01;nan;nan
BLUE;2009.45;nan;nan
GREEN;1831.88;nan;nan
RED;1492.12;nan;nan
NIR1;937.80;nan;nan
    """,
    'IKONOS': """
PAN;1353.25;1358.59;1364.06
BLUE;1921.26;1902.54;1901.19
GREEN;1803.28;1827.32;1826.04
RED;1517.76;1526.48;1532.48
NIR1;1145.8;1150.51;1155.37
    """}

SOURCES = ['Thuillier2003', 'ChKur', 'WRC']

DEFAULT_SOURCE = 'Thuillier2003'


def _parse_csv(sat_id):
    csv = CSV[sat_id]
    data = {}
    for line in csv.splitlines():
        if not line.strip():
            continue
        ss = line.split(';')
        if len(ss) != 4:
            raise ValueError('Data line corrupted: {} --> {}'.format(line, ss))
        band = ss[0]
        values = [float(v) for v in ss[1:]]
        data[band] = values
    return data


def get_solar_irradiance_dict(sat_id, source=DEFAULT_SOURCE):
    try:
        i = SOURCES.index(source)
    except ValueError:
        raise ValueError('Source must be one of {}. Got {}.'.format(SOURCES, source))
    data = _parse_csv(sat_id)
    return {band: data[band][i] for band in data}


def get_solar_irradiance_values(sat_id, source=DEFAULT_SOURCE):
    """Returns values for multispectral bands (see bands.BANDS_MULT)"""
    data = get_solar_irradiance_dict(sat_id=sat_id, source=source)
    return bands.get_values_sorted(data, sat_id=sat_id)
