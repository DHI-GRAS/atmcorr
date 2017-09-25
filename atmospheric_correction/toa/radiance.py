from atmospheric_correction.sensors import sensor_is


def toa_radiance(
        data, sensor, mtdFile, band_ids,
        doDOS=False, mtdFile_tile=None):
    """Compute TOA radiance

    Parameters
    ----------
    data : ndarray shape(nbands, ny, nx)
        input data
    sensor : str
        sensor name
    mtdFile : str
        path to metadata file
    band_ids : list of int
        band IDs of original product contained in array
        0-based
    doDOS : bool
        do a dark object subtraction
    mtdFile_tile : str
        tile metadata file
        required for Sentinel 2
    """
    commonkw = dict(
            data=data,
            mtdFile=mtdFile,
            doDOS=doDOS,
            band_ids=band_ids)
    if sensor_is(sensor, 'WV'):
        return toa_radiance_WV(sensor=sensor, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        return toa_radiance_PHR1(**commonkw)
    elif sensor_is(sensor, 'L7L8'):
        return toa_radiance_L8(sensor=sensor, **commonkw)
    elif sensor_is(sensor, 'S2'):
        commonkw.pop('doDOS')
        return toa_radiance_S2(mtdFile_tile=mtdFile_tile, **commonkw)
