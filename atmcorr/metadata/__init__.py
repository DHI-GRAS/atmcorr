from atmcorr.sensors import sensor_is
from atmcorr.sensors import sensor_is_any


def get_date(sensor, mtdFile):
    if sensor_is_any(sensor, 'WV', 'WV_4band'):
        from . import worldview
        return worldview.get_date(mtdFile)
    elif sensor_is(sensor, 'L7L8'):
        from . import landsat8
        return landsat8.get_date(mtdFile)
    elif sensor_is(sensor, "PHR"):
        from . import pleiades
        return pleiades.get_date(mtdFile)
    elif sensor_is(sensor, 'S2'):
        from . import sentinel2
        return sentinel2.get_date(mtdFile)
    else:
        raise ValueError('Unknown sensor.')
