import satmeta.pneo.parser as plparser


def get_geometry(mtdFile):
    copy_keys = ['sun_azimuth', 'sensor_azimuth', 'sensor_zenith']
    metadata = plparser.parse_metadata(mtdFile)
    gdict = {key: metadata['angles'][key] for key in copy_keys}
    gdict['sun_zenith'] = 90.0 - metadata['angles']['sun_elevation']
    gdict['month'] = metadata['sensing_time'].month
    gdict['day'] = metadata['sensing_time'].day
    return gdict
