from atmospheric_correction.sensors import sensor_is


def get_sensing_date(sensor, mtdFile):
    if sensor_is(sensor, 'WV'):
        return get_date_wv(mtdFile)
    elif sensor_is(sensor, 'L7L8'):
        return get_date_l7l8(mtdFile)
    elif sensor_is(sensor, "PHR"):
        return get_date_phr(mtdFile)
    elif sensor_is(sensor, 'S2'):
        metadata = s2meta.parse_mtdfile(mtdFile)
        return metadata['sensing_time']
    else:
        raise ValueError('Unknown sensor.')
