import os
import copy
import logging

import numpy as np
import dateutil
import rasterio

from atmospheric_correction import wrap_6S
from atmospheric_correction import sensors
from atmospheric_correction import io_utils
from atmospheric_correction.toa.radiance import toa_radiance
from atmospheric_correction.toa.reflectance import toa_reflectance
import atmospheric_correction.metadata.bands as meta_bands
import atmospheric_correction.metadata.dates as meta_dates

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
        atm, aeroProfile,
        dnFile=None, data=None, profile=None,
        tileSizePixels=0,
        isPan=False, adjCorr=True,
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
    atm : dict
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
    isPan : bool
        use Pan band only
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
                ''.format(profile['count'], data.shape))

    if use_modis:
        if not HAS_MODIS:
            raise ValueError(
                    'To download MODIS parameters, you need to '
                    'install the `modis_atm` package.')
        if not earthdata_credentials:
            raise ValueError(
                    'To download MODIS parameters, you need to '
                    'provide your earthdata_credentials (https://earthdata.nasa.gov/).')

    sensor_group = sensors.sensor_group_bands(sensor)
    if band_ids is None:
        try:
            band_ids = meta_bands.default_band_ids[sensor_group]
            logger.info('Assuming original set of %d bands.', len(band_ids))
        except KeyError:
            pass

    kwargs_toa_radiance = dict(
            isPan=isPan,
            mtdFile=mtdFile,
            mtdFile_tile=mtdFile_tile,
            band_ids=band_ids)

    if date is None:
        date = meta_dates.get_sensing_date(sensor, mtdFile)
    else:
        date = dateutil.parser.parse(date)

    if data is None:
        with rasterio.open(dnFile) as src:
            data = src.read()
            profile = src.profile.copy()

    nbands = profile['count']
    if len(band_ids) != nbands:
        raise ValueError(
                'Number of band IDs {} does not correspond to number of bands {}.',
                ''.format(len(band_ids), nbands))

    nodata = profile['nodata']
    if nodata is not None:
        if nodata not in [0, 65536]:
            data[data == nodata] = 0
            profile['nodata'] = 0

    # keep unchanged copy of atm dict
    if atm is not None:
        atm_original = copy.deepcopy(atm)
        atm = None
    else:
        atm_original = {'AOT': None, 'PWV': None, 'ozone': None}

    # DN -> Radiance -> Reflectance
    if method == "6S":
        doDOS = False
        res = profile['transform'].a

        logger.info('Computing TOA radiance ...')
        data = toa_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)
        if not np.any(data):
            raise RuntimeError('Data is all zeros.')

        if tileSizePixels > 0:
            tileExtents = io_utils.getTileExtents(
                    height=profile['height'],
                    width=profile['width'],
                    transform=profile['transform'],
                    tileSize=tileSizePixels)
        else:
            tileExtents = np.array([[None]])

        # Structure holding the 6S correction parameters has, for each band in
        # the image, arrays of values (one for each tile)
        # of the three correction parameters
        nrows, ncols = tileExtents.shape[:2]
        correctionParams = np.zeros(
                (nbands, nrows, ncols),
                dtype=dict(names=['xa', 'xb', 'xc'], formats=(['f4'] * 3)))

        nruns = nrows * ncols
        if nruns > 1 and adjCorr:
            raise ValueError('Adjacency correction only works for un-tiled image (tileSize=0)')

        # Get 6S correction parameters for an extent of each tile
        mysixs = None
        for j in range(nrows):
            for i in range(ncols):
                extent = tileExtents[j, i]
                logger.debug('tile %d,%d extent is %s', j, i, extent)
                atm = atm_original.copy()
                # If MODIS atmospheric data was downloaded then use it to set
                # different atmospheric parameters for each tile
                if use_modis:
                    raise NotImplementedError()
                    atm_modis = modis_atm.params.retrieve_parameters(
                            date=date,
                            extent=extent,
                            credentials=earthdata_credentials,
                            download_dir=modis_atm_dir)
                    atm = atm_original.copy()
                    for key in atm_modis:
                        if key not in atm or atm[key] is None:
                            atm[key] = atm_modis[key]

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
                        isPan=isPan,
                        aeroProfile=aeroProfile,
                        extent=extent,
                        nprocs=nprocs)

                nbands = len(tilecp)
                for b in range(nbands):
                    correctionParams[b, j, i]['xa'] = tilecp[b]['xa']
                    correctionParams[b, j, i]['xb'] = tilecp[b]['xb']
                    correctionParams[b, j, i]['xc'] = tilecp[b]['xc']

        logger.debug('Correction parameters: %s', np.squeeze(correctionParams))

        data = wrap_6S.perform_correction(
                data, correctionParams, pixel_size=res, adjCorr=adjCorr,
                mysixs=mysixs)

    elif method in ['DOS', 'TOA']:
        if method == 'DOS':
            doDOS = True
        else:
            doDOS = False

        if sensors.sensor_is(sensor, 'S2'):
            # S2 data is provided in L1C meaning in TOA reflectance
            data = toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)
        else:
            data = toa_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)
            data = toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)

    elif method == 'RAD':
        doDOS = False
        data = toa_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)

    else:
        raise ValueError('Unknown method \'{}\'.'.format(method))

    profile['dtype'] = 'float32'
    profile['nodata'] = np.nan
    if outfile is not None:
        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(data)
    else:
        if return_profile:
            return data, profile
        else:
            return data
