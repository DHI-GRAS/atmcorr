def sensor_is(sensor, key):
    if sensor == key:
        return True
    elif sensor in ['WV2', 'WV3']:
        return key == 'WV'
    elif sensor in ['PHR1A', 'PHR1B', 'SPOT6']:
        return key == 'PHR'
    elif sensor in ['L8', 'L7']:
        return key in ['L7L8', 'L7', 'L8']
    elif sensor.startswith('S2'):
        if sensor.endswith('10m'):
            return key == 'S2_10m'
        elif sensor.endswith('60m'):
            return key == 'S2_60m'
        else:
            return key == 'S2'
    raise ValueError('Unknown sensor \'{}\'.'.format(sensor))
