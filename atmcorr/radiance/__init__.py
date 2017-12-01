from atmcorr.sensors import sensor_is
from atmcorr.sensors import sensor_is_any


def dn_to_radiance(
        dndata, sensor, mtdFile, band_ids, mtdFile_tile=None):
    """Compute TOA radiance

    Parameters
    ----------
    dndata : ndarray shape(nbands, ny, nx)
        input data: digital numbers
        except for Sentinel2
    sensor : str
        sensor name
    mtdFile : str
        path to metadata file
    band_ids : list of int
        band IDs of original product contained in array
        0-based
    mtdFile_tile : str
        tile metadata file
        required for Sentinel 2
    """
    commonkw = dict(
        dndata=dndata,
        mtdFile=mtdFile,
        band_ids=band_ids)
    if sensor_is_any(sensor, 'WV', 'WV_4band'):
        from . import worldview
        return worldview.dn_to_radiance(**commonkw)
    elif sensor_is(sensor, 'PHR'):
        from . import pleiades
        return pleiades.dn_to_radiance(**commonkw)
    elif sensor_is(sensor, 'S2'):
        from . import sentinel2
        return sentinel2.toa_reflectance_to_radiance(
            mtdFile_tile=mtdFile_tile, **commonkw)
    elif sensor_is(sensor, 'L8'):
        from . import landsat8
        return landsat8.dn_to_radiance(**commonkw)
    else:
        raise NotImplementedError(sensor)
