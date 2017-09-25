import numpy as np

from atmospheric_correction.sentinel2 import metadata as metamod


def get_geometry(mtdFile, mtdFile_tile):
    """Get geometry dictionaty for S2

    Parameters
    ----------
    mtfile_file
    """
    if mtdFile_tile is None:
        raise ValueError('Tile metadata file is required.')
    metadata = metamod.parse_mtdfile(mtdFile, mtdFile_tile=mtdFile_tile)
    tile = list(metadata['granules'])[0]
    metadata_granule = metadata['granules'][tile]
    gdict = {}
    # Get the granule metadata from metadata dictionary
    copy_keys = ['sun_zenith', 'sun_azimuth']
    for key in copy_keys:
        gdict[key] = metadata_granule[key]
    average_keys = ['sensor_zenith', 'sensor_azimuth']
    for key in average_keys:
        gdict[key] = np.mean(metadata_granule[key])
    d = metadata['sensing_time']
    gdict['month'] = d.month
    gdict['day'] = d.day
    return gdict
