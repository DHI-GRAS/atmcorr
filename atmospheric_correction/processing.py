import copy
import logging

import numpy as np
import gdal
import gdal_utils.gdal_utils as gu
import dateutil

from atmospheric_correction import wrap_6S
from atmospheric_correction import sensors
from atmospheric_correction import io_utils
from atmospheric_correction import toa
import atmospheric_correction.metadata.bands as meta_bands
import atmospheric_correction.metadata.dates as meta_dates

# import modis_atm

logger = logging.getLogger(__name__)


def main_optdict(options):
    main(**options)


def main(
        sensor, dnFile, mtdFile, method,
        atm, aeroProfile, tileSizePixels=0,
        isPan=False, adjCorr=True, use_modis=False,
        aotMultiplier=1.0, roiFile=None, nprocs=None,
        mtdFile_tile=None, band_ids=None, date=None):
    """Main workflow function for atmospheric correction

    Parameters
    ----------
    sensor : str
        S2, WV2, WV3, PHR1A, PHR1B, SPOT6,
        L7, L8, S2
    dnFile : str
        path to digital numbers input file
    mtdFile : str
        path to mtdFile file
    method : str
        6S, RAD, DOS, TOA (same as DOS)
    atm : dict
        atmospheric parameters
    aeroProfile : dict
        TODO: what is this?
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
    """
    # keep unchanged copy
    if atm is not None:
        atm_original = copy.deepcopy(atm)
        atm = None
    else:
        atm_original = {'AOT': None, 'PWV': None, 'ozone': None}

    kwargs_toa_radiance = dict(
            isPan=isPan,
            mtdFile=mtdFile,
            mtdFile_tile=mtdFile_tile,
            band_ids=band_ids)

    sensor_group = sensors.sensor_group_bands(sensor)
    if band_ids is None:
        try:
            band_ids = meta_bands.default_band_ids[sensor_group]
        except KeyError:
            pass

    if date is None:
        date = meta_dates.get_sensing_date(sensor, mtdFile)
    else:
        date = dateutil.parser.parse(date)

    # DN -> Radiance -> Reflectance
    reflectanceImg = None
    if method == "6S":
        doDOS = False

        if roiFile is not None:
            logger.info('Clipping image to ROI ...')
            img = gu.cutline_to_shape_name(dnFile, roiFile)
        else:
            img = gu.gdal_open(dnFile)

            if band_ids is None:
                band_ids = list(range(img.GetRasterCount()))

        logger.info('Computing TOA radiance ...')
        radianceImg = toa.toa_radiance(img, sensor, doDOS=doDOS, **kwargs_toa_radiance)
        img = None

        if tileSizePixels > 0:
            tileExtents = io_utils.getTileExtents(radianceImg, tileSizePixels)
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
        correctionParams = [_empty.copy() for _ in range(radianceImg.RasterCount)]

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
                    reflectanceImg = wrap_6S.perform_correction(
                            radianceImg, correctionParams, adjCorr=True, mysixs=mysixs)
                    break

        if tileSizePixels > 0 or not adjCorr:
            logger.info('Perform atm correction')
            reflectanceImg = wrap_6S.perform_correction(
                    radianceImg, correctionParams, adjCorr=False)

    elif method in ["DOS", "TOA"]:
        if method == "DOS":
            doDOS = True
        else:
            doDOS = False
        if roiFile is not None:
            img = gu.cutline_to_shape_name(dnFile, roiFile)
        else:
            img = gdal.open(dnFile)
        if img is None:
            raise RuntimeError('Unable to read dnFile.')

        if sensors.sensor_is(sensor, 'S2'):
            # S2 data is provided in L1C meaning in TOA reflectance
            reflectanceImg = toa.toa_reflectance(img, mtdFile, sensor)
        else:
            radianceImg = toa.toa_radiance(img, sensor, doDOS=doDOS, **kwargs_toa_radiance)
            reflectanceImg = toa.toa_reflectance(radianceImg, mtdFile, sensor)
            radianceImg = None
        img = None

    elif method == "RAD":
        doDOS = False
        img = gu.cutline_to_shape_name(dnFile, roiFile)
        radianceImg = toa.toa_radiance(img, sensor, doDOS=doDOS, **kwargs_toa_radiance)
        reflectanceImg = radianceImg

    else:
        raise ValueError('Unknown method \'{}\'.'.format(method))

    return reflectanceImg
