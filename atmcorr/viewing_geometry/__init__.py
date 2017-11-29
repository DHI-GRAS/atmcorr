from atmcorr.sensors import sensor_is
from atmcorr.sensors import sensor_is_any

REQUIRED_KEYS = {'sun_zenith', 'sun_azimuth', 'sensor_zenith', 'sensor_azimuth', 'month', 'day'}


def _check_gdict(d):
    dkeys = set(d)
    if not dkeys >= REQUIRED_KEYS:
        raise ValueError(
                'Geometry dict is lacking the following keys: {}'.format(REQUIRED_KEYS - dkeys))


def get_geometry(sensor, mtdFile, **kwargs):
    if sensor_is_any(sensor, 'WV', 'WV_4band'):
        from . import worldview
        gdict = worldview.get_geometry(mtdFile)
    elif sensor_is(sensor, 'PHR'):
        from . import pleiades
        gdict = pleiades.get_geometry(mtdFile)
    elif sensor_is(sensor, 'L7L8'):
        from . import landsat8
        gdict = landsat8.get_geometry(mtdFile)
    elif sensor_is(sensor, 'S2'):
        from . import sentinel2
        gdict = sentinel2.get_geometry(mtdFile, **kwargs)
    _check_gdict(gdict)
    return gdict
