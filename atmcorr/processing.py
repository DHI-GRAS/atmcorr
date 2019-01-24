import os
import copy
import logging
import datetime
import itertools
import functools
import concurrent.futures
import multiprocessing
import platform

import numpy as np
import dateutil
import tqdm

from atmcorr import tiling
from atmcorr import wrap_6S
from atmcorr import metadata
from atmcorr import radiance
from atmcorr import resampling
from atmcorr import response_curves
from atmcorr import viewing_geometry
from atmcorr.sensors import check_sensor_supported
from atmcorr.adjacency_correction import adjacency_correction

try:
    import modis_atm.params
    HAS_MODIS = True
except ImportError:
    HAS_MODIS = False

logger = logging.getLogger(__name__)

MODIS_ATM_DIR = os.path.expanduser(os.path.join('~', 'MODIS_ATM'))


def _earthdata_credentials_from_env():
    auth = dict(
        username=os.environ.get('EARTHDATA_USERNAME'),
        password=os.environ.get('EARTHDATA_PASSWORD'))
    if None in list(auth.values()):
        raise ValueError(
            'Environment variables ERTHDATA_USERNAME and EARTHDATA_PASSWORD not set.')
    return auth


def main_optdict(options):
    main(**options)


def main(
        data, profile,
        sensor, mtdFile,
        band_ids,
        sixs_params,
        tileSize=0,
        adjCorr=True,
        aotMultiplier=1.0,
        mtdFile_tile=None, date=None,
        use_modis=False, modis_atm_dir=MODIS_ATM_DIR,
        earthdata_credentials=None, nthreads=None):
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
    band_ids : list of int
        bands of data
        0-based indices of full band stack
    sixs_params : dict
        6S parameters passed to setup_sixs
        keys: aeroProfile, atm, is_ocean
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
    nthreads : int, optional
        number of threads to use
        default: cpu count * 2

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
    else:
        if not sixs_params.get('atm', None):
            raise ValueError(
                'If not using MODIS, you must provide atmospheric '
                'composition in sixs_params[\'atm\'].')

    if date is None:
        date = metadata.get_date(sensor, mtdFile)
    elif not isinstance(date, datetime.datetime):
        date = dateutil.parser.parse(date)

    nbands, height, width = datashape = data.shape

    if len(band_ids) != nbands:
        raise ValueError(
            'Number of band IDs ({}) does not correspond to number of bands ({}).'
            .format(len(band_ids), nbands))

    # convert input data to float
    nodata = profile['nodata']
    if nodata is not None:
        data_float = data.astype('float32')
        data_float[data == nodata] = np.nan
        data = data_float

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
        all_atm.append(sixs_params['atm'])
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
            _sixs_params = sixs_params.copy()
            _sixs_params['atm'] = all_atm[min(ktile, len(all_atm)-1)]
            _sixs_params['geometry'] = all_geom[ktile]
            jobs = wrap_6S.generate_jobs(
                rcurves_dict=rcurves_dict,
                sixs_params=_sixs_params)
            for b, job in enumerate(jobs):
                yield tuple(job) + (b, j, i)

    # prepare parallel execution
    jobgen = _job_generator()
    njobs = tiling_shape[0] * tiling_shape[1] * nbands
    pbar = functools.partial(
        tqdm.tqdm, total=njobs, desc='Getting 6S params', unit='job', smoothing=0)

    # initialize processing pool
    if nthreads is None:
        nthreads = multiprocessing.cpu_count() * 2
    nthreads = min((nthreads, njobs))
    if platform.system() == 'Windows':
        # Windows doesn't seem to like too many starting processes at a time
        nthreads = min(nthreads, 36)
    # execute 6S jobs
    with concurrent.futures.ThreadPoolExecutor(nthreads) as executor:
        for tilecp, adjcoef, idx in pbar(executor.map(wrap_6S.run_sixs_job, jobgen)):
            b, j, i = idx
            for field in correction_params.dtype.names:
                correction_params[field][b, j, i] = tilecp[field]
            if adjCorr:
                adjacency_params[:, b, j, i] = adjcoef

    # reproject parameters
    if tiling_shape == (1, 1):
        # extract params from recarray
        corrparams = {
            field: correction_params[field][:, 0, 0]
            for field in correction_params.dtype.names}
        if adjCorr:
            adjparams = adjacency_params[:, :, 0, 0]
    else:
        # resample 6S correction parameters to image
        logger.info('Resampling correction parameters')
        corrparams = np.empty(datashape, dtype=correction_params.dtype)
        for field in corrparams.dtype.names:
            corrparams[field] = resampling.resample(
                source=correction_params[field],
                src_transform=tiling_transform,
                src_crs=profile['crs'],
                dst_transform=profile['transform'],
                dst_shape=datashape)
        if adjCorr:
            # resample adjacency correction parameters to image
            logger.info('Resampling adjacency correction parameters')
            # collapse first two dimensions
            collapsed_in = adjacency_params.reshape((-1, ) + adjacency_params.shape[2:])
            dst_shape = (collapsed_in.shape[0], ) + datashape[1:]
            collapsed_out = resampling.resample(
                source=collapsed_in,
                src_transform=tiling_transform,
                src_crs=profile['crs'],
                dst_transform=profile['transform'],
                dst_shape=dst_shape)
            # unpack dimensions again
            final_shape = (adjacency_params.shape[0], ) + datashape
            adjparams = collapsed_out.reshape(final_shape)

    if not adjCorr:
        adjparams = None

    # DN -> Radiance
    data = radiance.dn_to_radiance(
        data, sensor,
        mtdFile=mtdFile,
        mtdFile_tile=mtdFile_tile,
        band_ids=band_ids,
        dst_transform=profile['transform']
    )

    # Radiance -> Reflectance
    # apply 6s correction parameters
    data = wrap_6S.perform_correction(data, corrparams)

    if adjCorr:
        # perform adjecency correction
        data = adjacency_correction(
            data, *adjparams,
            pixel_size=profile['transform'].a
        )

    profile['dtype'] = 'float32'
    profile['nodata'] = np.nan
    return data, profile


def _copy_check_profile(profile):
    profile = copy.deepcopy(profile)
    required = {'nodata', 'transform', 'crs'}
    missing = required - set(profile)
    if missing:
        raise ValueError(
            'Metadata dict `profile` is missing information: \'{}\''.format(missing))
    return profile
