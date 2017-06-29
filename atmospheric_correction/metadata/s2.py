import os
import re
import glob
from xml.etree import ElementTree as ET


def tile_from_fname(fname):
    """Get S2 tile from file name"""
    fname = os.path.basename(fname)
    try:
        re.search('\d{2}[A-Z]{3}', fname).group(0)
    except AttributeError:
        raise ValueError(
                'Unable to get tile from fname \'{}\'.'
                ''.format(fname))


def _find_granule_metadata_relpath(mtdfile, granule):
    granuleDir = os.path.join(
            os.path.dirname(mtdfile), "GRANULE", granule)
    pattern = os.path.join(granuleDir, '*.xml')
    try:
        return glob.glob(pattern)[0]
    except IndexError:
        raise RuntimeError(
                'Unable to find granule metadata files with pattern \'{}\'.'
                ''.format(pattern))


def parse_granule_mtdfile(mtdfile_tile, granulename=None, granulekey=None):
    # read metadata of tile
    tree = ET.parse(mtdfile_tile)
    root = tree.getroot()
    namespace = root.tag.split('}')[0]+'}'
    # Get sun geometry - use the mean
    baseNodePath = "./"+namespace+"Geometric_Info/Tile_Angles/"
    sunGeometryNodeName = baseNodePath+"Mean_Sun_Angle/"
    sunZen = root.find(sunGeometryNodeName+"ZENITH_ANGLE").text
    sunAz = root.find(sunGeometryNodeName+"AZIMUTH_ANGLE").text
    # Get sensor geometry - assume that all bands have the same angles
    # (they differ slightly)
    sensorGeometryNodeName = baseNodePath+"Mean_Viewing_Incidence_Angle_List/Mean_Viewing_Incidence_Angle/"
    sensorZen = root.find(sensorGeometryNodeName+"ZENITH_ANGLE").text
    sensorAz = root.find(sensorGeometryNodeName+"AZIMUTH_ANGLE").text
    EPSG = tree.find("./"+namespace+"Geometric_Info/Tile_Geocoding/HORIZONTAL_CS_CODE").text
    cldCoverPercent = tree.find("./"+namespace+"Quality_Indicators_Info/Image_Content_QI/CLOUDY_PIXEL_PERCENTAGE").text
    for elem in tree.iter(tag='Size'):
        if elem.attrib['resolution'] == '10':
            rows_10 = int(elem[0].text)
            cols_10 = int(elem[1].text)
        if elem.attrib['resolution'] == '20':
            rows_20 = int(elem[0].text)
            cols_20 = int(elem[1].text)
        if elem.attrib['resolution'] == '60':
            rows_60 = int(elem[0].text)
            cols_60 = int(elem[1].text)
    for elem in tree.iter(tag='Geoposition'):
        if elem.attrib['resolution'] == '10':
            ULX_10 = int(elem[0].text)
            ULY_10 = int(elem[1].text)
        if elem.attrib['resolution'] == '20':
            ULX_20 = int(elem[0].text)
            ULY_20 = int(elem[1].text)
        if elem.attrib['resolution'] == '60':
            ULX_60 = int(elem[0].text)
            ULY_60 = int(elem[1].text)

    # save to dictionary
    if granulekey is None:
        if granulename is None:
            raise ValueError('granulename must be specified if granulekey is None')
        granulekey = granulename[len(granulename)-13:-7]
    return {
        granulekey: {
            'sun_zenit': sunZen,
            'sun_azimuth': sunAz,
            'sensor_zenit': sensorZen,
            'sensor_azimuth': sensorAz,
            'projection': EPSG,
            'cloudCoverPercent': cldCoverPercent,
            'rows_10': rows_10,
            'cols_10': cols_10,
            'rows_20': rows_20,
            'cols_20': cols_20,
            'rows_60': rows_60,
            'cols_60': cols_60,
            'ULX_10': ULX_10,
            'ULY_10': ULY_10,
            'ULX_20': ULX_20,
            'ULY_20': ULY_20,
            'ULX_60': ULX_60,
            'ULY_60': ULY_60}}


def parse_mtdfile_S2L1C(mtdfile, mtdfile_tile=None, tile=None):

    # Get parameters from main metadata file
    ProductName = os.path.basename(os.path.dirname(mtdfile))
    tree = ET.parse(mtdfile)
    root = tree.getroot()
    namespace = root.tag.split('}')[0]+'}'

    baseNodePath = "./"+namespace+"General_Info/Product_Info/"
    dateTimeStr = root.find(baseNodePath+"PRODUCT_START_TIME").text
    procesLevel = root.find(baseNodePath+"PROCESSING_LEVEL").text
    spaceCraft = root.find(baseNodePath+"Datatake/SPACECRAFT_NAME").text
    orbitDirection = root.find(baseNodePath+"Datatake/SENSING_ORBIT_DIRECTION").text

    baseNodePath = "./"+namespace+"General_Info/Product_Image_Characteristics/"
    quantificationVal = root.find(baseNodePath+"QUANTIFICATION_VALUE").text
    reflectConversion = root.find(baseNodePath+"Reflectance_Conversion/U").text
    irradianceNodes = root.findall(baseNodePath+"Reflectance_Conversion/Solar_Irradiance_List/SOLAR_IRRADIANCE")
    e0 = []
    for node in irradianceNodes:
        e0.append(node.text)

    # save to dictionary
    metadict = {}
    metadict.update({'product_name':ProductName,
                     'product_start':dateTimeStr,
                     'processing_level':procesLevel,
                     'spacecraft':spaceCraft,
                     'orbit_direction':orbitDirection,
                     'quantification_value':quantificationVal,
                     'reflection_conversion':reflectConversion,
                     'irradiance_values': e0})
    # granules
    metadict['granules'] = []

    if mtdfile_tile is not None:
        gdict = parse_granule_mtdfile(mtdfile_tile, granulekey=tile)
        metadict.update(gdict)
        metadict['granules'].append(tile)
        return metadict

    glist = list(tree.iter(tag='Granules'))
    for elem in glist:
        granulename = elem.attrib['granuleIdentifier']
        mtdfile_tile = _find_granule_metadata_relpath(mtdfile, granulename=granulename)

        gdict = parse_granule_mtdfile(mtdfile_tile, granulename=granulename)
        metadict.update(gdict)
        granulekey = list(metadict.keys())[0]
        metadict['granules'].append(granulekey)
    return metadict


def parse_mtdfile_S2L2A(mtdfile, mtdfile_tile):
    # Get parameters from main metadata file
    ProductName = os.path.basename(os.path.dirname(mtdfile))
    tree = ET.parse(mtdfile)
    root = tree.getroot()
    namespace = root.tag.split('}')[0]+'}'

    baseNodePath = "./"+namespace+"General_Info/L2A_Product_Info/"
    dateTimeStr = root.find(baseNodePath+"PRODUCT_START_TIME").text
    procesLevel = root.find(baseNodePath+"PROCESSING_LEVEL").text
    spaceCraft = root.find(baseNodePath+"Datatake/SPACECRAFT_NAME").text
    orbitDirection = root.find(baseNodePath+"Datatake/SENSING_ORBIT_DIRECTION").text

    baseNodePath = "./"+namespace+"General_Info/L2A_Product_Image_Characteristics/"
    quantificationValBOA = root.find(baseNodePath+"L1C_L2A_Quantification_Values_List/L2A_BOA_QUANTIFICATION_VALUE").text
    quantificationValAOT = root.find(baseNodePath+"L1C_L2A_Quantification_Values_List/L2A_AOT_QUANTIFICATION_VALUE").text
    quantificationValWVP = root.find(baseNodePath+"L1C_L2A_Quantification_Values_List/L2A_WVP_QUANTIFICATION_VALUE").text
    reflectConversion = root.find(baseNodePath+"Reflectance_Conversion/U").text
    irradianceNodes = root.findall(baseNodePath+"Reflectance_Conversion/Solar_Irradiance_List/SOLAR_IRRADIANCE")
    e0 = []
    for node in irradianceNodes:
        e0.append(node.text)

    # save to dictionary
    metadict = {}
    metadict.update({'product_name':ProductName,
                     'product_start':dateTimeStr,
                     'processing_level':procesLevel,
                     'spacecraft':spaceCraft,
                     'orbit_direction':orbitDirection,
                     'quantification_value_boa':quantificationValBOA,
                     'quantification_value_aot':quantificationValAOT,
                     'quantification_value_wvp':quantificationValWVP,
                     'reflection_conversion':reflectConversion,
                     'irradiance_values': e0})
    # granules
    glist = list(tree.iter(tag='Granules'))
    metadict['granules'] = []
    for elem in glist:
        granule = elem.attrib['granuleIdentifier']

        if mtdfile_tile is not None:
            if len(glist) > 1:
                raise ValueError(
                        'mtdfile_tile should only be used for single-tile '
                        'products. Fount {}.'.format(len(glist)))
        else:
            mtdfile_tile = _find_granule_metadata_relpath(mtdfile, granule)

        # read metadata of tile
        tree = ET.parse(mtdfile_tile)
        root = tree.getroot()
        namespace = root.tag.split('}')[0]+'}'
        # Get sun geometry - use the mean
        baseNodePath = "./"+namespace+"Geometric_Info/Tile_Angles/"
        sunGeometryNodeName = baseNodePath+"Mean_Sun_Angle/"
        sunZen = root.find(sunGeometryNodeName+"ZENITH_ANGLE").text
        sunAz = root.find(sunGeometryNodeName+"AZIMUTH_ANGLE").text
        # Get sensor geometry - assume that all bands have the same angles
        # (they differ slightly)
        sensorGeometryNodeName = baseNodePath+"Mean_Viewing_Incidence_Angle_List/Mean_Viewing_Incidence_Angle/"
        sensorZen = root.find(sensorGeometryNodeName+"ZENITH_ANGLE").text
        sensorAz = root.find(sensorGeometryNodeName+"AZIMUTH_ANGLE").text
        EPSG = tree.find("./"+namespace+"Geometric_Info/Tile_Geocoding/HORIZONTAL_CS_CODE").text
        for elem in tree.iter(tag='Size'):
            if elem.attrib['resolution'] == '10':
                rows_10 = int(elem[0].text)
                cols_10 = int(elem[1].text)
            if elem.attrib['resolution'] == '20':
                rows_20 = int(elem[0].text)
                cols_20 = int(elem[1].text)
            if elem.attrib['resolution'] == '60':
                rows_60 = int(elem[0].text)
                cols_60 = int(elem[1].text)
        for elem in tree.iter(tag='Geoposition'):
            if elem.attrib['resolution'] == '10':
                ULX_10 = int(elem[0].text)
                ULY_10 = int(elem[1].text)
            if elem.attrib['resolution'] == '20':
                ULX_20 = int(elem[0].text)
                ULY_20 = int(elem[1].text)
            if elem.attrib['resolution'] == '60':
                ULX_60 = int(elem[0].text)
                ULY_60 = int(elem[1].text)

        # save to dictionary
        granulekey = granule[len(granule)-13:-7]
        metadict['granules'].append(granulekey)
        metadict.update({
            granulekey: {
                'sun_zenit': sunZen,
                'sun_azimuth': sunAz,
                'sensor_zenit': sensorZen,
                'sensor_azimuth': sensorAz,
                'projection': EPSG,
                'rows_10': rows_10,
                'cols_10': cols_10,
                'rows_20': rows_20,
                'cols_20': cols_20,
                'rows_60': rows_60,
                'cols_60': cols_60,
                'ULX_10': ULX_10,
                'ULY_10': ULY_10,
                'ULX_20': ULX_20,
                'ULY_20': ULY_20,
                'ULX_60': ULX_60,
                'ULY_60': ULY_60}})
    return metadict
