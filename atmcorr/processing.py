import os
import copy
import logging

import numpy as np
import dateutil
import rasterio

from atmcorr import wrap_6S
from atmcorr import sensors
from atmcorr import tiling
from atmcorr.sensors import sensor_is

try:
    import modis_atm.params
    HAS_MODIS = True
except ImportError:
    HAS_MODIS = False

logger = logging.getLogger(__name__)

MODIS_ATM_DIR = os.path.expanduser(os.path.join('~', 'MODIS_ATM'))


def main_optdict(options):
    main(**options)


def main(
        sensor, mtdFile, method,
        aeroProfile, atm={},
        dnFile=None, data=None, profile=None,
        tileSizePixels=0,
        adjCorr=True,
        aotMultiplier=1.0, nprocs=None,
        mtdFile_tile=None, band_ids=None, date=None,
        outfile=None, return_profile=False,
        use_modis=False, modis_atm_dir=MODIS_ATM_DIR,
        earthdata_credentials={}):
    """Main workflow function for atmospheric correction

    Parameters
    ----------
    sensor : str
        S2, WV2, WV3, PHR1A, PHR1B, SPOT6,
        L7, L8, S2
    mtdFile : str
        path to mtdFile file
    method : str
        6S, RAD, DOS, TOA (same as DOS)
    atm : dict, optional
        atmospheric parameters
    aeroProfile : str
        aero profile name
    dnFile : str, optional
        path to digital numbers input file
        instead, you can also provide the data
    data : ndarray, optional shape(nbands, ny, nx)
        digital numbers input data
        can be used instead of dnFile
    profile : dict, optional
        rasterio file profile
        required when using data
    tileSizePixels : int
        tile size in pixels
    adjCorr : bool
        perform adjacency correction
    aotMultiplier : float
        Atmospheric Optical Depth
        scale factor
    nprocs : int
        number of processors to use
        default: all available
    mtdFile_tile : str
        path to tile mtdFile
        required for S2
    band_ids : list of int or str
        band IDs (0-based) from complete
        sensor band set
        required if not full se is used
    date : datetime.datetime, optional
        image date
        will be retrieved from metadata
        if not specified
    outfile : str
        path to output file
    return_profile : bool
        when returning data,
        also return profile
    use_modis : bool
        download atm parameters from MODIS
    modis_atm_dir : str
        path where to store / cache modis
        files for re-use
    earthdata_credentials : dict
        username, password
        required for use_modis

    Returns
    -------
    None
        if outfile is specified
    data
        if not return_profile
        nodata are masked as NaN
    data, profile
        if return_profile
        where profile is the rasterio file profile
        and data as above
    """
    if data is not None and profile is None:
        raise ValueError('Data and profile must be provided together.')

    if data is not None and data.shape[0] != profile['count']:
        raise ValueError(
                'Data array must have bands as its first dimension.'
                'Number of bands was {} but data has shape {}.'
                .format(profile['count'], data.shape))

    if use_modis:
        if not HAS_MODIS:
            raise ValueError(
                    'To download MODIS parameters, you need to '
                    'install the `modis_atm` package.')
        if not earthdata_credentials:
            raise ValueError(
                    'To download MODIS parameters, you need to '
                    'provide your earthdata_credentials (https://earthdata.nasa.gov/).')
        if not os.path.isdir(modis_atm_dir):
            os.makedirs(modis_atm_dir)

    if date is None:
        date = _get_sensing_date(sensor, mtdFile)
    else:
        date = dateutil.parser.parse(date)

    if data is None:
        with rasterio.open(dnFile) as src:
            data = src.read()
            profile = copy.copy(src.profile)

    nbands = profile['count']
    if len(band_ids) != nbands:
        raise ValueError(
                'Number of band IDs ({}) does not correspond to number of bands ({}).'
                .format(len(band_ids), nbands))

    if band_ids is None:
        band_ids = np.arange(nbands)

    kwargs_toa_radiance = dict(
            mtdFile=mtdFile,
            mtdFile_tile=mtdFile_tile,
            band_ids=band_ids)

    nodata = profile['nodata']
    if nodata is not None:
        if nodata not in [0, 65536]:
            data[data == nodata] = 0
            profile['nodata'] = 0

    # keep unchanged copy of atm dict
    if atm is None:
        atm = {'AOT': None, 'PWV': None, 'ozone': None}

    # DN -> Radiance -> Reflectance
    if method == '6S':
        data = _main_6S(
                data, profile, band_ids, sensor, date,
                mtdFile, mtdFile_tile,
                atm, kwargs_toa_radiance, tileSizePixels,
                adjCorr, aotMultiplier, aeroProfile,
                use_modis, modis_atm_dir, earthdata_credentials,
                nprocs)

    elif method in ['DOS', 'TOA']:
        if method == 'DOS':
            doDOS = True
        else:
            doDOS = False

        if sensors.sensor_is(sensor, 'S2'):
            # S2 data is provided in L1C meaning in TOA reflectance
            data = _toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)
        else:
            data = _toa_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)
            data = _toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)

    elif method == 'RAD':
        doDOS = False
        data = _toa_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)

    else:
        raise ValueError('Unknown method \'{}\'.'.format(method))

    profile['dtype'] = 'float32'
    profile['nodata'] = np.nan
    logger.debug(profile)
    if outfile is not None:
        logger.info('Saving output to \'%s\'', outfile)
        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(data)
    else:
        if return_profile:
            return data, profile
        else:
            return data


def _main_6S(
        data, profile, band_ids, sensor, date,
        mtdFile, mtdFile_tile,
        atm, kwargs_toa_radiance, tileSizePixels,
        adjCorr, aotMultiplier, aeroProfile,
        use_modis, modis_atm_dir, earthdata_credentials,
        nprocs):

    res = profile['transform'].a

    logger.info('Computing TOA radiance ...')
    data = _toa_radiance(data, sensor, doDOS=False, **kwargs_toa_radiance)
    if not np.any(data):
        raise RuntimeError('Data is all zeros.')

    nbands, height, width = data.shape
    profile.update(width=width, height=height)

    if tileSizePixels > 0:
        tilegrid_kw = dict(
                xtilesize=tileSizePixels,
                ytilesize=tileSizePixels)
    else:
        tilegrid_kw = dict(
                xtilesize=width,
                ytilesize=height)
    tile_extents = tiling.get_tile_extents(
            height=height,
            width=width,
            src_transform=profile['transform'],
            src_crs=profile['crs'],
            **tilegrid_kw)
    if tile_extents.shape == ():
        tile_extents = tile_extents.reshape((1, 1))

    # Structure holding the 6S correction parameters has, for each band in
    # the image, arrays of values (one for each tile)
    # of the three correction parameters
    njtiles, nitiles = tile_extents.shape
    correctionParams = np.zeros(
            (nbands, njtiles, nitiles),
            dtype=dict(names=['xa', 'xb', 'xc'], formats=(['f4'] * 3)))

    nruns = njtiles * nitiles
    if nruns > 1 and adjCorr:
        raise ValueError('Adjacency correction only works for un-tiled image (tileSize=0)')

    # Get 6S correction parameters for an extent of each tile
    mysixs = None
    for j in range(njtiles):
        for i in range(nitiles):
            extent = dict(zip(tile_extents.dtype.names, tile_extents[j, i]))
            logger.debug('tile %d,%d extent is %s', j, i, extent)
            atm = copy.copy(atm)
            # If MODIS atmospheric data was downloaded then use it to set
            # different atmospheric parameters for each tile
            if use_modis:
                logger.info('Retrieving MODIS atmospheric parameters')
                atm_modis = modis_atm.params.retrieve_parameters(
                        date=date,
                        extent=extent,
                        credentials=earthdata_credentials,
                        download_dir=modis_atm_dir)
                logger.info(atm_modis)
                atm = atm_modis
                for k in atm:
                    if atm[k] is None:
                        raise ValueError(
                                'One or more atm parameters could not be retrieved: {}.'
                                .format(atm))

            atm['AOT'] *= aotMultiplier
            logger.debug('AOT: %s', atm['AOT'])
            logger.debug('Water Vapour: %s', atm['PWV'])
            logger.debug('Ozone: %s', atm['ozone'])

            mysixs, tilecp = wrap_6S.get_correction_params(
                    sensor=sensor,
                    mtdFile=mtdFile,
                    mtdFile_tile=mtdFile_tile,
                    atm=atm,
                    band_ids=band_ids,
                    aeroProfile=aeroProfile,
                    extent=extent,
                    nprocs=nprocs)

            nbands = len(tilecp)
            for b in range(nbands):
                correctionParams[b, j, i]['xa'] = tilecp[b]['xa']
                correctionParams[b, j, i]['xb'] = tilecp[b]['xb']
                correctionParams[b, j, i]['xc'] = tilecp[b]['xc']
    logger.debug('Correction parameters: %s', correctionParams)

    logger.info('Apply parameters to image')
    data = wrap_6S.perform_correction(
            data, correctionParams, pixel_size=res, adjCorr=adjCorr,
            mysixs=mysixs)

    return data


def _toa_radiance(
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
        from atmcorr import worldview
        return worldview.radiance.toa_radiance(sensor=sensor, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        from atmcorr import pleiades
        return pleiades.radiance.toa_radiance(**commonkw)
    elif sensor_is(sensor, 'L7L8'):
        from atmcorr import landsat8
        return landsat8.radiance.toa_radiance(sensor=sensor, **commonkw)
    elif sensor_is(sensor, 'S2'):
        commonkw.pop('doDOS')
        from atmcorr import sentinel2
        return sentinel2.radiance.toa_radiance(mtdFile_tile=mtdFile_tile, **commonkw)


def _toa_reflectance(data, mtdfile, sensor, band_ids):
    commonkw = dict(data=data, mtdfile=mtdfile)
    if sensor_is(sensor, 'WV'):
        from atmcorr import worldview
        res = worldview.reflectance.toa_reflectance(band_ids=band_ids, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        from atmcorr import pleiades
        res = pleiades.reflectance.toa_reflectance(**commonkw)
    elif sensor_is(sensor, 'L7L8'):
        from atmcorr import landsat8
        res = landsat8.reflectange.toa_reflectance(**commonkw)
    elif sensor_is(sensor, 'S2'):
        from atmcorr import sentinel2
        res = sentinel2.reflectance.toa_reflectance(**commonkw)
    return res


def _get_sensing_date(sensor, mtdFile):
    if sensor_is(sensor, 'WV'):
        from atmcorr import worldview
        return worldview.metadata.get_date(mtdFile)
    elif sensor_is(sensor, 'L7L8'):
        from atmcorr import landsat8
        return landsat8.metadata.get_date(mtdFile)
    elif sensor_is(sensor, "PHR"):
        from atmcorr import pleiades
        return pleiades.metadata.get_date(mtdFile)
    elif sensor_is(sensor, 'S2'):
        from atmcorr import sentinel2
        return sentinel2.metadata.get_date(mtdFile)
    else:
        raise ValueError('Unknown sensor.')
