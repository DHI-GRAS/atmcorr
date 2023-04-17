SUPPORTED_SENSORS = [
    'WV2', 'WV3', 'WV4', 'GE1',
    'PHR1A', 'PHR1B', 'SPOT6',
    'L7', 'L8', 'S2A', 'S2B',
    'PNEO3', 'PNEO4']


def check_sensor_supported(sensor):
    if sensor not in SUPPORTED_SENSORS:
        raise ValueError(
            'Sensor \'{}\' not among supported ({}).'
            .format(sensor, SUPPORTED_SENSORS))


def sensor_is(sensor, target):
    if sensor == target:
        return True
    elif sensor in ['WV2', 'WV3']:
        return target == 'WV'
    elif sensor in ['WV4', 'GE1']:
        return target == 'WV_4band'
    elif sensor in ['PHR1A', 'PHR1B', 'SPOT6']:
        return target == 'PHR'
    elif sensor in ['L8', 'L7']:
        return target in ['L7L8', 'L7', 'L8']
    elif sensor in ['S2A', 'S2B']:
        return target == 'S2'
    elif sensor.startswith("PNEO"):
        return target == 'PNEO'
    raise ValueError('Unknown sensor \'{}\'.'.format(sensor))


def sensor_is_any(sensor, *targets):
    res = False
    for target in targets:
        res = res or sensor_is(sensor, target)
    return res


def sensor_group_bands(sensor):
    """Get target of sensors that have similar sets of bands"""
    for target in ['WV', 'PHR', 'L7', 'L8', 'S2', 'WV_4band', 'PNEO']:
        if sensor_is(sensor, target):
            return target
    raise ValueError('Unable to get sensor group for \'%s\'.'.format(sensor))
