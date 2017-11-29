import os
import copy
import logging
import datetime
import itertools
import concurrent.futures
import multiprocessing

import numpy as np
import dateutil
import tqdm

from atmcorr import wrap_6S
from atmcorr import tiling
from atmcorr import metadata
from atmcorr import radiance
from atmcorr.sensors import sensor_is
from atmcorr.sensors import sensor_is_any
from atmcorr.sensors import check_sensor_supported
from atmcorr import viewing_geometry
from atmcorr import response_curves
from atmcorr.adjacency_correction import adjacency_correction
from atmcorr import resampling

try:
    import modis_atm.params
    HAS_MODIS = True
except ImportError:
    HAS_MODIS = False

logger = logging.getLogger(__name__)

MODIS_ATM_DIR = os.path.expanduser(os.path.join('~', 'MODIS_ATM'))
NUM_PROCESSES = None


def _earthdata_credentials_from_env():
    auth = dict(
        username=os.environ.get('EARTHDATA_USERNAME'),
        password=os.environ.get('EARTHDATA_PASSWORD'))
    if None in list(auth.values):
        raise ValueError(
            'Environment variables ERTHDATA_USERNAME and EARTHDATA_PASSWORD not set.')
    return auth


def main_optdict(options):
    main(**options)


def main(
        data, profile,
        sensor, mtdFile, method,
        aeroProfile, atm={},
        tileSize=0,
        adjCorr=True,
        aotMultiplier=1.0,
        mtdFile_tile=None, band_ids=None, date=None,
        use_modis=False, modis_atm_dir=MODIS_ATM_DIR,
        earthdata_credentials=None):
    """Main workflow function for atmospheric correction

    Parameters
    ----------
    data : ndarray, shape(nbands, ny, nx)
        digital numbers input data
    profile : dict
        rasterio file profile
    sensor : str
        see atmcorr.sensors.SUPPORTED_SENSORS
    mtdFile : str
        path to mtdFile file
    method : str
        6S, RAD, DOS, TOA (same as DOS)
    atm : dict, optional
        atmospheric parameters
    aeroProfile : str
        aero profile name
    tileSize : int
        tile size in pixels
    adjCorr : bool
        perform adjacency correction
    aotMultiplier : float
        Atmospheric Optical Depth
        scale factor
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
    data, profile
        nodata are masked as NaN
        where profile is the modified rasterio file profile
    """
    check_sensor_supported(sensor)

    profile = _copy_check_profile(profile)

    if use_modis:
        if not HAS_MODIS:
            raise ValueError(
                'To download MODIS parameters, you need to '
                'install the `modis_atm` package.')
        if not earthdata_credentials:
            try:
                earthdata_credentials = _earthdata_credentials_from_env()
            except ValueError:
                raise ValueError(
                    'To download MODIS parameters, you need to '
                    'provide your earthdata_credentials (https://earthdata.nasa.gov/).')
        if not os.path.isdir(modis_atm_dir):
            os.makedirs(modis_atm_dir)

    if date is None:
        date = metadata.get_date(sensor, mtdFile)
    elif isinstance(date, datetime.datetime):
        pass
    else:
        date = dateutil.parser.parse(date)

    nbands = data.shape[0]
    if len(band_ids) != nbands:
        raise ValueError(
            'Number of band IDs ({}) does not correspond to number of bands ({}).'
            .format(len(band_ids), nbands))

    if band_ids is None:
        band_ids = np.arange(nbands)

    if atm is None:
        atm = {'AOT': None, 'PWV': None, 'ozone': None}

    kwargs_toa_radiance = dict(
            mtdFile=mtdFile,
            mtdFile_tile=mtdFile_tile,
            band_ids=band_ids)

    nodata = profile['nodata']
    if nodata is not None:
        if nodata not in [0, 65536]:
            data[data == nodata] = 0
            profile['nodata'] = 0

    # DN -> Radiance -> Reflectance
    if method == '6S':
        data = _main_6S(
                data, profile, band_ids, sensor, date,
                mtdFile, mtdFile_tile,
                atm, kwargs_toa_radiance, tileSize,
                adjCorr, aotMultiplier, aeroProfile,
                use_modis, modis_atm_dir, earthdata_credentials)
    elif method in ['DOS', 'TOA']:
        doDOS = (method == 'DOS')
        if sensor_is(sensor, 'S2'):
            # S2 data is provided in L1C meaning in TOA reflectance
            data = _toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)
        else:
            data = radiance.dn_to_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)
            data = _toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)
    elif method == 'RAD':
        doDOS = False
        data = radiance.dn_to_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)
    else:
        raise ValueError('Unknown method \'{}\'.'.format(method))

    profile['dtype'] = 'float32'
    profile['nodata'] = np.nan
    return data, profile


def _main_6S(
        data, profile, band_ids, sensor, date,
        mtdFile, mtdFile_tile,
        atm, kwargs_toa_radiance, tileSize,
        adjCorr, aotMultiplier, aeroProfile,
        use_modis, modis_atm_dir, earthdata_credentials):

    # convert to radiance
    data = _toa_radiance(data, sensor, doDOS=False, **kwargs_toa_radiance)
    if not np.any(data):
        raise RuntimeError('Data is all zeros.')

    nbands, height, width = data.shape

    if tileSize:
        tiling_transform, tiling_shape = tiling.get_tiled_transform_shape(
            src_transform=profile['transform'], src_shape=(height, width), dst_res=tileSize)
        tile_extents_wgs = tiling.get_projected_extents(
            transform=tiling_transform,
            height=tiling_shape[0], width=tiling_shape[1],
            src_crs=profile['crs'])
    else:
        tile_extents_wgs = tiling.get_projected_image_extent(
            transform=profile['transform'],
            height=height, width=width,
            src_crs=profile['crs'])
        tiling_transform = None
        tiling_shape = (1, 1)

    geometry_dict = viewing_geometry.get_geometry(
        sensor, mtdFile,
        mtdFile_tile=mtdFile_tile,
        dst_transform=tiling_transform,
        dst_shape=tiling_shape)

    rcurves_dict = response_curves.get_response_curves(sensor, band_ids)

    # Structure holding the 6S correction parameters has, for each band in
    # the image, arrays of values (one for each tile)
    # of the three correction parameters
    correction_params = np.zeros(
        ((nbands, ) + tiling_shape),
        dtype=dict(names=['xa', 'xb', 'xc'], formats=(['f4'] * 3)))
    if adjCorr:
        adjacency_params = np.zeros(((4, nbands) + tiling_shape), dtype='f4')

    def _get_tile_index_iter():
        return itertools.product(range(tiling_shape[0]), range(tiling_shape[1]))

    # retrieve MODIS parameters for all tiles
    all_atm = []
    if not use_modis:
        all_atm.append(atm)
    else:
        for j, i in _get_tile_index_iter():
            extent = tiling.recarr_take_dict(tile_extents_wgs, j, i)
            logger.debug('tile %d,%d extent is %s', j, i, extent)
            # If MODIS atmospheric data was downloaded then use it to set
            # different atmospheric parameters for each tile
            logger.info('Retrieving MODIS atmospheric parameters')
            atm_modis = modis_atm.params.retrieve_parameters(
                date=date,
                extent=extent,
                credentials=earthdata_credentials,
                download_dir=modis_atm_dir)
            logger.info(atm_modis)
            # check missing
            missing = []
            for key in atm_modis:
                if atm_modis[key] is None:
                    missing.append(key)
            if missing:
                raise ValueError(
                    'Some atm parameters could not be retrieved ({}): {}.'
                    .format(missing, atm_modis))
            all_atm.append(atm_modis)

    # apply aotMultiplier
    for atm_tile in all_atm:
        atm_tile['AOT'] *= aotMultiplier

    # extract geometry dict for each tile
    all_geom = []
    if tiling_shape == (1, 1):
        # get angles for whole image
        all_geom.append(geometry_dict)
    else:
        for j, i in _get_tile_index_iter():
            # extract angles for this tile
            geometry_dict_tile = {}
            for key, value in geometry_dict.items():
                geometry_dict_tile[key] = (
                    value[j, i] if isinstance(value, np.ndarray) else value)
            all_geom.append(geometry_dict_tile)

    # Get 6S correction parameters for each tile
    def _job_generator():
        for ktile, (j, i) in enumerate(_get_tile_index_iter()):
            atm_tile = all_atm[min(ktile, len(all_atm)-1)]
            geom_tile = all_geom[ktile]
            jobs = wrap_6S.generate_jobs(
                atm=atm_tile,
                geometry_dict=geom_tile,
                rcurves_dict=rcurves_dict,
                aeroProfile=aeroProfile)
            for b, job in enumerate(jobs):
                yield tuple(job) + (b, j, i)

    njobs = tiling_shape[0] * tiling_shape[1] * nbands
    pbarkw = dict(total=njobs, desc='Getting 6S params', unit='job', smoothing=0)

    # initialize processing pool
    nprocs = NUM_PROCESSES
    if nprocs is None:
        nprocs = multiprocessing.cpu_count()
    # execute 6S jobs
    with concurrent.futures.ProcessPoolExecutor(nprocs) as executor:
        for tilecp, adjcoef, idx in tqdm.tqdm(
                executor.map(wrap_6S.run_sixs_job, _job_generator()),
                **pbarkw):
            b, j, i = idx
            for field in correction_params.dtype.names:
                correction_params[field][b, j, i] = tilecp[field]
            if adjCorr:
                adjacency_params[:, b, j, i] = adjcoef

    # reproject parameters
    if tiling_shape == (1, 1):
        corrparams = {
            field: correction_params[field][:, 0, 0]
            for field in correction_params.dtype.names}
        if adjCorr:
            adjparams = adjacency_params[:, :, 0, 0]
    else:
        logger.debug('resampling correction parameters')
        # resample 6S correction parameters to image
        corrparams = np.empty(data.shape, dtype=correction_params.dtype)
        for field in corrparams.dtype.names:
            corrparams[field] = resampling.resample(
                source=correction_params[field],
                src_transform=tiling_transform,
                src_crs=profile['crs'],
                dst_transform=profile['transform'],
                dst_shape=data.shape)
        if adjCorr:
            # resample adjacency correction parameters to image
            # collapse first two dimensions
            collapsed_in = adjacency_params.reshape((-1, ) + adjacency_params.shape[2:])
            dst_shape = (collapsed_in.shape[0], ) + data.shape[1:]
            collapsed_out = resampling.resample(
                source=collapsed_in,
                src_transform=tiling_transform,
                src_crs=profile['crs'],
                dst_transform=profile['transform'],
                dst_shape=dst_shape)
            # unpack dimensions again
            final_shape = (adjacency_params.shape[0], ) + data.shape
            adjparams = collapsed_out.reshape(final_shape)

    # apply 6s correction parameters
    data = wrap_6S.perform_correction(data, corrparams)

    if adjCorr:
        # perform adjecency correction
            data = adjacency_correction(
                data,
                *adjparams,
                pixel_size=profile['transform'].a)

    return data


def _toa_reflectance(data, mtdfile, sensor, band_ids):
    commonkw = dict(data=data, mtdfile=mtdfile)
    if sensor_is_any(sensor, 'WV', 'WV_4band'):
        from atmcorr import worldview
        return worldview.reflectance.toa_reflectance(band_ids=band_ids, **commonkw)
    elif sensor_is(sensor, 'PHR'):
        from atmcorr import pleiades
        return pleiades.reflectance.toa_reflectance(**commonkw)
    elif sensor_is(sensor, 'S2'):
        from atmcorr import sentinel2
        return sentinel2.reflectance.toa_reflectance(**commonkw)
    else:
        raise NotImplementedError(sensor)


def _copy_check_profile(profile):
    profile = copy.deepcopy(profile)
    required = {'nodata', 'transform', 'crs'}
    missing = required - set(profile)
    if missing:
        raise ValueError(
            'Metadata dict `profile` is missing information: \'{}\''.format(missing))
    return profile
