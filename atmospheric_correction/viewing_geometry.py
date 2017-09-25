from atmospheric_correction.sensors import sensor_is

REQUIRED_KEYS = {'sun_zenith', 'sun_azimuth', 'sensor_zenith', 'sensor_azimuth', 'month', 'day'}


def _check_gdict(d):
    dkeys = set(d)
    if not dkeys >= REQUIRED_KEYS:
        raise ValueError(
                'Geometry dict is lacking the following keys: {}'.format(REQUIRED_KEYS - dkeys))


def get_geometry(sensor, mtdFile, mtdFile_tile=None):
    if sensor_is(sensor, 'WV'):
        from atmospheric_correction import worldview
        gdict = worldview.geometry.get_geometry(mtdFile)
    elif sensor_is(sensor, 'PHR'):
        from atmospheric_correction import pleiades
        gdict = pleiades.geometry.get_geometry(mtdFile)
    elif sensor_is(sensor, 'L7L8'):
        from atmospheric_correction import landsat8
        gdict = landsat8.geometry.get_geometry(mtdFile)
    elif sensor_is(sensor, 'S2'):
        from atmospheric_correction import sentinel2
        gdict = sentinel2.geometry.get_geometry(mtdFile, mtdFile_tile)
    _check_gdict(gdict)
    return gdict
