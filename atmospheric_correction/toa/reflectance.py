from atmospheric_correction.sensors import sensor_is


def toa_reflectance(data, mtdfile, sensor, band_ids):
    commonkw = dict(data=data, mtdfile=mtdfile)
    if sensor_is(sensor, 'WV'):
        res = toa_reflectance_WV(band_ids=band_ids, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        res = toa_reflectance_PHR1(**commonkw)
    elif sensor_is(sensor, 'L7L8'):
        res = toa_reflectance_L8(**commonkw)
    elif sensor_is(sensor, 'S2'):
        res = toa_reflectance_S2(**commonkw)
    return res
