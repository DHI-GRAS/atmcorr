import copy
import logging

import gdal
import gdal_utils.gdal_utils as gu

from . import modis_params
from . import wrap_6S
from . import metadata
from .sensors import sensor_is
from .io_utils import getTileExtents
from .toa.radiance import toa_radiance
from .toa.reflectance import toa_reflectance

logger = logging.getLogger(__name__)


def main_optdict(options):
    main(**options)


def main(
        sensor, dnFile, mtdfile, method,
        atm, aeroProfile, tileSizePixels,
        isPan=False, adjCorr=True,
        aotMultiplier=1.0, roiFile=None, nprocs=None,
        mtdfile_tile=None, tile=None, band_ids=None):
    """Main workflow function for atmospheric correction

    Parameters
    ----------
    sensor : str
        S2, WV2, WV3, PHR1A, PHR1B, SPOT6,
        L7, L8, S2
    dnFile : str
        path to digital numbers input file
    mtdfile : str
        path to mtdfile file
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
    aotMultiplier : float
        Atmospheric Optical Depth
        scale factor
    roiFile : str, optional
        path to ROI file to clip the image with
    nprocs : int
        number of processors to use
        default: all available
    mtdfile_tile : str
        path to tile mtdfile
        required for S2
    tile : str
        tile name for S2
        e.g. 09ABC
    band_ids : list of int or str
        band IDs (0-based) from complete
        sensor band set
        required if not full se is used
    """
    # keep unchanged copy
    atm_original = copy.deepcopy(atm)

    mtd_dict = mtdfile
    if sensor_is(sensor, 'S2'):
        mtd_dict = metadata.readMetadataS2L1C(
                mtdfile=mtdfile,
                mtdfile_tile=mtdfile_tile)
        logger.debug('S2 mtd dict:\n%s', mtd_dict)
        mtd_dict.update({
            'current_granule': tile,
            'band_ids': band_ids})

    # DN -> Radiance -> Reflectance
    if method == "6S":
        doDOS = False

        if roiFile is not None:
            logger.info('Clipping image to ROI ...')
            img = gu.cutline_to_shape_name(dnFile, roiFile)
        else:
            img = gdal.open(dnFile)
        if img is None:
            raise RuntimeError('Unable to read dnFile.')

        logger.info('Computing TOA radiance ...')
        radianceImg = toa_radiance(img, mtd_dict, sensor, doDOS=doDOS, isPan=isPan)

        img = None
        reflectanceImg = None

        tileExtents = [[None]]
        if tileSizePixels > 0:
            tileExtents = getTileExtents(radianceImg, tileSizePixels)

        # If atmospheric parameters needed by 6S are not specified then
        # donwload and use MODIS atmopsheric products
        modisAtmDir = None
        if not (atm['AOT'] and atm['PWV'] and atm['ozone']):
            logger.info('Retrieving MODIS atmospheric parameters ...')
            modisAtmDir = modis_params.downloadAtmParametersMODIS(dnFile, mtd_dict, sensor)

        # Structure holding the 6S correction parameters has for each band in
        # the image a dictionary with arrays of values (one for each tile)
        # of the three correction parameter
        n_first = len(tileExtents[0])
        n_extents = len(tileExtents)
        correctionParams = (
                [{
                    'xa': [[0] * n_first] * n_extents,
                    'xb': [[0] * n_first] * n_extents,
                    'xc': [[0] * n_first] * n_extents}] *
                radianceImg.RasterCount)

        # Get 6S correction parameters for an extent of each tile
        for y, tileRow in enumerate(tileExtents):
            for x, extent in enumerate(tileRow):
                # If MODIS atmospheric data was downloaded then use it to set
                # different atmospheric parameters for each tile
                if modisAtmDir:
                    atm = atm_original
                    aot, pwv, ozone = modis_params.estimateAtmParametersMODIS(
                            dnFile, modisAtmDir, extent=extent, yearDoy="",
                            time=-1, roiShape=None)

                    if not atm['AOT']:
                        atm['AOT'] = aot
                    if not atm['PWV']:
                        atm['PWV'] = pwv
                    if not atm['ozone']:
                        atm['ozone'] = ozone

                atm['AOT'] *= aotMultiplier
                logger.debug("AOT: " + str(atm['AOT']))
                logger.debug("Water Vapour: " + str(atm['PWV']))
                logger.debug("Ozone: " + str(atm['ozone']))

                s, tileCorrectionParams = wrap_6S.getCorrectionParams6S(
                        sensor=sensor, mtd_dict=mtd_dict, atm=atm, isPan=isPan,
                        aeroProfile=aeroProfile, extent=extent, nprocs=nprocs)

                for band, bandCorrectionParams in enumerate(tileCorrectionParams):
                    correctionParams[band]['xa'][y][x] = bandCorrectionParams['xa']
                    correctionParams[band]['xb'][y][x] = bandCorrectionParams['xb']
                    correctionParams[band]['xc'][y][x] = bandCorrectionParams['xc']
                if tileSizePixels == 0 and adjCorr:
                    reflectanceImg = wrap_6S.performAtmCorrection(
                            radianceImg, correctionParams, adjCorr, s)

        if tileSizePixels > 0 or not adjCorr:
            logger.info('Perform atm correction')
            reflectanceImg = wrap_6S.performAtmCorrection(radianceImg, correctionParams, s=None)

        if modisAtmDir:
            logger.info('MODIS cleanup')
            modis_params.deleteDownloadedModisFiles(modisAtmDir)
        radianceImg = None

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

        if sensor_is(sensor, 'S2'):
            # S2 data is provided in L1C meaning in TOA reflectance
            reflectanceImg = toa_reflectance(img, mtd_dict, sensor)
        else:
            radianceImg = toa_radiance(img, mtd_dict, sensor, doDOS=doDOS)
            reflectanceImg = toa_reflectance(radianceImg, mtd_dict, sensor)
            radianceImg = None
        img = None

    elif method == "RAD":
        doDOS = False
        img = gu.cutline_to_shape_name(dnFile, roiFile)
        radianceImg = toa_radiance(img, mtd_dict, sensor, doDOS=doDOS)
        reflectanceImg = radianceImg

    return reflectanceImg
