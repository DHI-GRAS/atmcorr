import numpy as np
import satmeta.s2.meta as s2meta
from satmeta.s2 import angles_2d as s2angles


def _get_angles_means(gmeta):
    gdict = {}
    copy_keys = ['sun_zenith', 'sun_azimuth']
    for key in copy_keys:
        gdict[key] = gmeta[key]
    average_keys = ['sensor_zenith', 'sensor_azimuth']
    for key in average_keys:
        gdict[key] = np.mean(gmeta[key])
    return gdict


def _reformat_angles_dict(angles):
    gdict = {}
    for ksrc, kdst in zip(['Viewing_Incidence', 'Sun'], ['sensor', 'sun']):
        for kang in ['Zenith', 'Azimuth']:
            key = '{}_{}'.format(kdst, kang.lower())
            gdict[key] = angles[ksrc][kang]
    return gdict


def _get_angles_2d(mtdfile_tile, dst_transform, dst_shape):
    angles = s2angles.parse_resample_angles(
        mtdfile_tile,
        dst_transform=dst_transform,
        dst_shape=dst_shape,
        resample_method='rasterio')
    return _reformat_angles_dict(angles)


def get_geometry(mtdFile, mtdFile_tile, dst_transform=None, dst_shape=None):
    """Get geometry dictionaty for S2"""
    if mtdFile_tile is None:
        raise ValueError('Tile metadata file is required.')
    meta = s2meta.parse_metadata(mtdFile)
    if dst_transform is None and dst_shape is None:
        gmeta = s2meta.parse_granule_metadata(mtdFile_tile)
        angles = _get_angles_means(gmeta)
    else:
        angles = _get_angles_2d(
            mtdFile_tile, dst_transform=dst_transform, dst_shape=dst_shape)
    gdict = {}
    gdict.update(angles)
    # Get the granule metadata from metadata dictionary
    d = meta['sensing_time']
    gdict['month'] = d.month
    gdict['day'] = d.day
    return gdict
