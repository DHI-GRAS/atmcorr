import copy
import logging

import numpy as np
import gdal
import gdal_utils.gdal_utils as gu

from atmospheric_correction import modis_params
from atmospheric_correction import wrap_6S
from atmospheric_correction import sensors
from atmospheric_correction import io_utils
from atmospheric_correction import toa
from atmospheric_correction import metadata as metamod

logger = logging.getLogger(__name__)


def main_optdict(options):
    main(**options)


def main(
        sensor, dnFile, mtdFile, method,
        atm, aeroProfile, tileSizePixels,
        isPan=False, adjCorr=True, use_modis=False,
        aotMultiplier=1.0, roiFile=None, nprocs=None,
        mtdFile_tile=None, band_ids=None):
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
    """
    # keep unchanged copy
    if atm is not None:
        atm_original = copy.deepcopy(atm)
        atm = None  # unassign atm

    kwargs_toa_radiance = dict(
            isPan=isPan,
            mtdFile=mtdFile,
            mtdFile_tile=mtdFile_tile,
            band_ids=band_ids)

    sensor_group = sensors.sensor_group_bands(sensor)
    if band_ids is None:
        band_ids = metamod.bands.default_band_ids[sensor_group]

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
        radianceImg = toa.toa_radiance(img, sensor, doDOS=doDOS, **kwargs_toa_radiance)

        img = None
        reflectanceImg = None

        tileExtents = [[None]]
        if tileSizePixels > 0:
            tileExtents = io_utils.getTileExtents(radianceImg, tileSizePixels)

        # If atmospheric parameters needed by 6S are not specified then
        # donwload and use MODIS atmopsheric products
        modisAtmDir = None
        if atm_original is None:
            if use_modis:
                logger.info('Retrieving MODIS atmospheric parameters ...')
                modisAtmDir = modis_params.downloadAtmParametersMODIS(dnFile, mtdFile, sensor)
            else:
                atm_original = {'AOT': -1, 'PWV': -1, 'ozone': -1}

        # Structure holding the 6S correction parameters has for each band in
        # the image a dictionary with arrays of values (one for each tile)
        # of the three correction parameter
        n_first = len(tileExtents[0])
        n_extents = len(tileExtents)
        empty = np.zeros(
                (n_extents, n_first),
                dtype=dict(names=['xa', 'xb', 'xc'], formats=(['f8'] * 3)))
        correctionParams = [empty.copy() for _ in range(radianceImg.RasterCount)]

        # Get 6S correction parameters for an extent of each tile
        for j, tileRow in enumerate(tileExtents):
            for i, extent in enumerate(tileRow):
                # If MODIS atmospheric data was downloaded then use it to set
                # different atmospheric parameters for each tile
                if modisAtmDir:
                    atm = atm_original.copy()
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
                logger.debug("AOT: %s", atm['AOT'])
                logger.debug("Water Vapour: %s", atm['PWV'])
                logger.debug("Ozone: %s", atm['ozone'])

                s, tileCorrectionParams = wrap_6S.getCorrectionParams6S(
                        sensor=sensor, mtdFile=mtdFile, mtdFile_tile=mtdFile_tile,
                        atm=atm, band_ids=band_ids, isPan=isPan,
                        aeroProfile=aeroProfile, extent=extent, nprocs=nprocs)

                for band, bandCorrectionParams in enumerate(tileCorrectionParams):
                    correctionParams[band]['xa'][j, i] = bandCorrectionParams['xa']
                    correctionParams[band]['xb'][j, i] = bandCorrectionParams['xb']
                    correctionParams[band]['xc'][j, i] = bandCorrectionParams['xc']
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

    return reflectanceImg
