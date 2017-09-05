import os
import copy
import logging

import numpy as np
import gdal_utils.gdal_binaries as gbin
import dateutil
import rasterio

from atmospheric_correction import wrap_6S
from atmospheric_correction import sensors
from atmospheric_correction import io_utils
from atmospheric_correction.toa.radiance import toa_radiance
from atmospheric_correction.toa.reflectance import toa_reflectance
import atmospheric_correction.metadata.bands as meta_bands
import atmospheric_correction.metadata.dates as meta_dates

# import modis_atm

logger = logging.getLogger(__name__)


def main_optdict(options):
    main(**options)


def main(
        sensor, mtdFile, method,
        atm, aeroProfile,
        dnFile=None, data=None, profile=None,
        tileSizePixels=0,
        isPan=False, adjCorr=True, use_modis=False,
        aotMultiplier=1.0, roiFile=None, nprocs=None,
        mtdFile_tile=None, band_ids=None, date=None,
        outfile=None):
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
    aeroProfile : dict
        TODO: what is this?
    dnFile : str, optional
        path to digital numbers input file
        instead, you can also provide the data
    data : ndarray, optional
        digital numbers input data
        can be used instead of dnFile
    profile : dict, optional
        rasterio file profile
        required when using data
    tileSizePixels : int
        tile size in pixels
    isPan : bool
        TODO: Who is Pan?
    adjCorr : bool
        perform adjacency correction
    use_modis : bool
        download atm parameters from MODIS
    aotMultiplier : float
        Atmospheric Optical Depth
        scale factor
    roiFile : str, optional
        path to ROI file to clip the image with
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
    """
    if data is not None and profile is None:
        raise ValueError('Data and profile must be provided together.')

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
        if roiFile is not None:
            logger.info('Clipping image to ROI ...')
            dnFile_clipped = '{}_clipped{}'.format(*os.path.splitext(dnFile))
            gbin.cutline(dnFile, inshp=roiFile, outfile=dnFile_clipped)
            dnFile = dnFile_clipped
            logger.info('Done clipping.')

        # load data
        with rasterio.open(dnFile) as src:
            data = src.read()
            profile = src.profile.copy()

    nbands = profile['count']
    if len(band_ids) != nbands:
        raise ValueError(
                'Number of band IDs {} does not correspond to number of bands {}.',
                ''.format(len(band_ids), nbands))

    # keep unchanged copy
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

        if tileSizePixels > 0:
            tileExtents = io_utils.getTileExtents(
                    height=profile['height'],
                    width=profile['width'],
                    transform=profile['transform'],
                    tileSize=tileSizePixels)
        else:
            tileExtents = np.array([[None]])

        # If atmospheric parameters needed by 6S are not specified then
        # donwload and use MODIS atmopsheric products
        if use_modis:
            logger.info('Retrieving MODIS atmospheric parameters ...')
            raise NotImplementedError()
            # modisAtmDir = modis_atm.get_atm(extent, date)

        # Structure holding the 6S correction parameters has for each band in
        # the image a dictionary with arrays of values (one for each tile)
        # of the three correction parameters
        nrows, ncols = tileExtents.shape[:2]
        _empty = np.zeros(
                (nrows, ncols),
                dtype=dict(names=['xa', 'xb', 'xc'], formats=(['f8'] * 3)))
        correctionParams = [_empty.copy() for _ in range(nbands)]

        # Get 6S correction parameters for an extent of each tile
        nruns = nrows * ncols
        for j in range(nrows):
            for i in range(ncols):
                extent = tileExtents[j, i]
                atm = atm_original.copy()
                # If MODIS atmospheric data was downloaded then use it to set
                # different atmospheric parameters for each tile
                if use_modis:
                    raise NotImplementedError()
                    atm_modis = None  # modis_atm.get_atm(date, extent)
                    atm = atm_original.copy()
                    for key in atm_modis:
                        if key not in atm or atm[key] is None:
                            atm[key] = atm_modis[key]

                atm['AOT'] *= aotMultiplier
                logger.debug("AOT: %s", atm['AOT'])
                logger.debug("Water Vapour: %s", atm['PWV'])
                logger.debug("Ozone: %s", atm['ozone'])

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
                    correctionParams[b][j, i]['xa'] = tilecp[b]['xa']
                    correctionParams[b][j, i]['xb'] = tilecp[b]['xb']
                    correctionParams[b][j, i]['xc'] = tilecp[b]['xc']

                if tileSizePixels == 0 and adjCorr:
                    if nruns != 1:
                        raise RuntimeError('Should be only one run with tileSizePixels=0')
                    data = wrap_6S.perform_correction(
                            data, correctionParams, adjCorr=True, mysixs=mysixs)
                    break

        if tileSizePixels > 0 or not adjCorr:
            logger.info('Perform atm correction')
            data = wrap_6S.perform_correction(
                    data, correctionParams, pixel_size=res, adjCorr=False)

    elif method in ["DOS", "TOA"]:
        if method == "DOS":
            doDOS = True
        else:
            doDOS = False

        if sensors.sensor_is(sensor, 'S2'):
            # S2 data is provided in L1C meaning in TOA reflectance
            data = toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)
        else:
            data = toa_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)
            data = toa_reflectance(data, mtdFile, sensor, band_ids=band_ids)

    elif method == "RAD":
        doDOS = False
        data = toa_radiance(data, sensor, doDOS=doDOS, **kwargs_toa_radiance)

    else:
        raise ValueError('Unknown method \'{}\'.'.format(method))

    if outfile is not None:
        profile['dtype'] = 'float32'
        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(data)
    else:
        return data
