from atmcorr.sensors import sensor_is
from atmcorr.sensors import sensor_is_any


def _toa_reflectance(data, mtdfile, sensor, band_ids):
    commonkw = dict(data=data, mtdfile=mtdfile)
    if sensor_is_any(sensor, 'WV', 'WV_4band'):
        from . import worldview
        return worldview.toa_reflectance(band_ids=band_ids, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        from . import pleiades
        return pleiades.toa_reflectance(**commonkw)
    elif sensor_is(sensor, 'PNEO'):
        from . import pneo
        return pneo.toa_reflectance(**commonkw)
    elif sensor_is(sensor, 'S2'):
        from . import sentinel2
        return sentinel2.calibrate_reflectance(**commonkw)
    else:
        raise NotImplementedError(sensor)
