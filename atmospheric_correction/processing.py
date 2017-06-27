import os
import re
import logging
import multiprocessing

from gdal_utils.gdal_utils import cutline_to_shape_name

from . import modis_params
from .sensors import sensor_is
from .satmeta import readMetadataS2L1C
from .wrap_6S import getCorrectionParams6S
from .wrap_6S import performAtmCorrection
from .io_utils import getTileExtents
from .toa.radiance import toaRadiance
from .toa.reflectance import toaReflectance

logger = logging.getLogger('atmProcessing')


def _tile_from_fname(fname):
    fname = os.path.basename(fname)
    tile = fname[len(fname) - 10:-4]
    if re.search('\d{2}[A-Z]{3}', tile) is not None:
        return tile
    else:
        raise ValueError(
                'Unable to get tile from fname \'{}\'.'
                ''.format(fname))


def atmProcessingMain(options):
    # Set the band numbers to the appropriate sensor
    sensor = options['sensor']

    # Commonly used filenames
    dnFile = options["dnFile"]
    metadata = options["metadataFile"]  # can be dictionary or file path
    roiFile = options.get("roiFile", "")

    # Correction options
    atmCorrMethod = options["atmCorrMethod"]

    # Make a copy of atm since it will be changing
    # in the code but the original is still required
    atm = options["atm"].copy()
    isPan = options["isPan"]
    adjCorr = options["adjCorr"]
    aeroProfile = options["aeroProfile"]
    tileSize = options["tileSizePixels"]
    aotMultiplier = options.get("aotMultiplier", 1.0)

    # special case for Sentinel-2 - read metadata in to dictionary
    if sensor_is(sensor, 'S2'):
        tile = options.get('tile', None)
        if tile is None:
            tile = _tile_from_fname(dnFile)
        metadata = readMetadataS2L1C(
                mtdfile=options["metadataFile"],
                mtdfile_tile=options.get('metadataFile_tile', None),
                band_ids=options.get('band_ids', None),
                tile=tile)
        logger.debug('S2 metadata:\n%s', metadata)
        # Add current granule (used to extract relevant metadata later...)
        metadata.update({'current_granule': tile})

    # DN -> Radiance -> Reflectance
    if atmCorrMethod == "6S":
        doDOS = False

        logger.info('Clipping image to ROI ...')
        inImg = cutline_to_shape_name(dnFile, roiFile)

        logger.info('Computing TOA radiance ...')
        radianceImg = toaRadiance(inImg, metadata, sensor, doDOS=doDOS, isPan=isPan)

        inImg = None
        reflectanceImg = None

        tileExtents = [[None]]
        if tileSize > 0:
            tileExtents = getTileExtents(radianceImg, tileSize)

        # If atmospheric parameters needed by 6S are not specified then
        # donwload and use MODIS atmopsheric products
        modisAtmDir = None
        if not (atm['AOT'] and atm['PWV'] and atm['ozone']):
            logger.info('Retrieving MODIS atmospheric parameters ...')
            modisAtmDir = modis_params.downloadAtmParametersMODIS(dnFile, metadata, sensor)

        # Structure holding the 6S correction parameters has for each band in
        # the image a dictionary with arrays of values (one for each tile)
        # of the three correction parameter
        correctionParams = [{'xa': [[0] * len(tileExtents[0])] * len(tileExtents),
                             'xb': [[0] * len(tileExtents[0])] * len(tileExtents),
                             'xc': [[0] * len(tileExtents[0])] * len(tileExtents)} for _ in
                            range(radianceImg.RasterCount)]

        # Get 6S correction parameters for an extent of each tile
        for y, tileRow in enumerate(tileExtents):
            for x, extent in enumerate(tileRow):
                # If MODIS atmospheric data was downloaded then use it to set
                # different atmospheric parameters for each tile
                if modisAtmDir:
                    atm = options["atm"].copy()
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
                s, tileCorrectionParams = getCorrectionParams6S(
                        metadata, atm=atm, sensor=sensor, isPan=isPan,
                        aeroProfile=aeroProfile, extent=extent,
                        nprocs=options.get("nprocs", multiprocessing.cpu_count()))

                for band, bandCorrectionParams in enumerate(tileCorrectionParams):
                    correctionParams[band]['xa'][y][x] = bandCorrectionParams['xa']
                    correctionParams[band]['xb'][y][x] = bandCorrectionParams['xb']
                    correctionParams[band]['xc'][y][x] = bandCorrectionParams['xc']
                if tileSize == 0 and adjCorr:
                    reflectanceImg = performAtmCorrection(
                            radianceImg, correctionParams, adjCorr, s)

        if tileSize > 0 or not adjCorr:
            logger.info('Perform atm correction')
            reflectanceImg = performAtmCorrection(radianceImg, correctionParams, s=None)

        if modisAtmDir:
            logger.info('MODIS cleanup')
            modis_params.deleteDownloadedModisFiles(modisAtmDir)
        radianceImg = None

    elif atmCorrMethod in ["DOS", "TOA"]:
        if atmCorrMethod == "DOS":
            doDOS = True
        else:
            doDOS = False
        inImg = cutline_to_shape_name(dnFile, roiFile)
        if not sensor_is(sensor, 'S2'):
            radianceImg = toaRadiance(inImg, metadata, sensor, doDOS=doDOS)
            inImg = None
            reflectanceImg = toaReflectance(radianceImg, metadata, sensor)
            radianceImg = None
        # S2 data is provided in L1C meaning in TOA reflectance
        else:
            reflectanceImg = toaReflectance(inImg, metadata, sensor)
            inImg = None

    elif atmCorrMethod == "RAD":
        doDOS = False
        inImg = cutline_to_shape_name(dnFile, roiFile)
        radianceImg = toaRadiance(inImg, metadata, sensor, doDOS=doDOS)
        reflectanceImg = radianceImg

    return reflectanceImg
