import os
import re
from xml.etree import ElementTree as ET

import dateutil.parser


def find_tile_name(fname):
    """Get S2 tile from file name"""
    fname = os.path.basename(fname)
    try:
        return re.search('\d{2}[A-Z]{3}', fname).group(0)
    except AttributeError:
        raise ValueError('Unable to find tile name in \'{}\'.'.format(fname))


def parse_granule_mtdfile(mtdFile_tile):
    gdict = {}
    # read metadata of tile
    tree = ET.parse(mtdFile_tile)
    root = tree.getroot()
    namespace = root.tag.split('}')[0]+'}'
    # find tile name
    baseNodePath = "./"+namespace+"General_Info/"
    tile_id = root.find(baseNodePath+"TILE_ID").text
    tile_name = find_tile_name(tile_id)
    # Get sun geometry - use the mean
    baseNodePath = "./"+namespace+"Geometric_Info/Tile_Angles/"
    sunGeometryNodeName = baseNodePath+"Mean_Sun_Angle/"
    sun_zenith = float(root.find(sunGeometryNodeName+"ZENITH_ANGLE").text)
    sun_azimuth = float(root.find(sunGeometryNodeName+"AZIMUTH_ANGLE").text)
    # Get sensor geometry - assume that all bands have the same angles
    # (they differ slightly)
    sensorGeometryNodeName = baseNodePath+"Mean_Viewing_Incidence_Angle_List/Mean_Viewing_Incidence_Angle/"
    sensor_zenith = float(root.find(sensorGeometryNodeName+"ZENITH_ANGLE").text)
    sensor_azimuth = float(root.find(sensorGeometryNodeName+"AZIMUTH_ANGLE").text)
    EPSG = tree.find("./"+namespace+"Geometric_Info/Tile_Geocoding/HORIZONTAL_CS_CODE").text
    cloud_cover_percentage = float(tree.find("./"+namespace+"Quality_Indicators_Info/Image_Content_QI/CLOUDY_PIXEL_PERCENTAGE").text)
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

    gdict.update({
        'sun_zenith': sun_zenith,
        'sun_azimuth': sun_azimuth,
        'sensor_zenith': sensor_zenith,
        'sensor_azimuth': sensor_azimuth,
        'projection': EPSG,
        'cloud_cover_percentage': cloud_cover_percentage,
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
        'ULY_60': ULY_60})

    # save to dictionary
    return {tile_name: gdict}


def parse_mtdfile(mtdFile, mtdFile_tile=None):

    # Get parameters from main metadata file
    ProductName = os.path.basename(os.path.dirname(mtdFile))
    tree = ET.parse(mtdFile)
    root = tree.getroot()
    namespace = root.tag.split('}')[0]+'}'

    baseNodePath = "./"+namespace+"General_Info/Product_Info/"
    sensing_time = dateutil.parser.parse(root.find(baseNodePath+"PRODUCT_START_TIME").text)
    procesLevel = root.find(baseNodePath+"PROCESSING_LEVEL").text
    spaceCraft = root.find(baseNodePath+"Datatake/SPACECRAFT_NAME").text
    orbitDirection = root.find(baseNodePath+"Datatake/SENSING_ORBIT_DIRECTION").text

    baseNodePath = "./"+namespace+"General_Info/Product_Image_Characteristics/"
    quantificationVal = float(root.find(baseNodePath+"QUANTIFICATION_VALUE").text)
    reflectConversion = float(root.find(baseNodePath+"Reflectance_Conversion/U").text)
    irradianceNodes = root.findall(baseNodePath+"Reflectance_Conversion/Solar_Irradiance_List/SOLAR_IRRADIANCE")
    e0 = []
    for node in irradianceNodes:
        e0.append(float(node.text))
    metadict = {
            'product_name': ProductName,
            'sensing_time': sensing_time,
            'processing_level': procesLevel,
            'spacecraft': spaceCraft,
            'orbit_direction': orbitDirection,
            'quantification_value': quantificationVal,
            'reflection_conversion': reflectConversion,
            'irradiance_values': e0}
    if mtdFile_tile is not None:
        metadict['granules'] = parse_granule_mtdfile(mtdFile_tile)
    return metadict
