import os
from xml.etree import ElementTree as ET
import glob

def readMetadataS2L1C(metadataFile):
    # Get parameters from main metadata file   
    ProductName = os.path.split(os.path.dirname(metadataFile))[1]
    tree = ET.parse(metadataFile)
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
    metaDict = {}
    metaDict.update({'product_name':ProductName,
                     'product_start':dateTimeStr,
                     'processing_level':procesLevel,
                     'spacecraft':spaceCraft,
                     'orbit_direction':orbitDirection,
                     'quantification_value':quantificationVal,
                     'reflection_conversion':reflectConversion,
                     'irradiance_values': e0})
    # granules
    XML_mask = 'S2?_*.xml'
    for elem in tree.iter(tag='Granules'):
        granule = elem.attrib['granuleIdentifier']
        granuleDir = os.path.join(os.path.dirname(metadataFile), "GRANULE", 
                                  granule)
        globlist = granuleDir + '/' + XML_mask
        metadataTile = glob.glob(globlist)[0]
        # read metadata of tile
        tree = ET.parse(metadataTile)
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
        metaDict.update({granule[len(granule)-13:-7]:{'sun_zenit':sunZen,
                                  'sun_azimuth':sunAz,
                                  'sensor_zenit':sensorZen,
                                  'sensor_azimuth':sensorAz,
                                  'projection':EPSG,
                                  'cloudCoverPercent':cldCoverPercent,
                                  'rows_10':rows_10,
                                  'cols_10':cols_10,
                                  'rows_20':rows_20,
                                  'cols_20':cols_20,
                                  'rows_60':rows_60,
                                  'cols_60':cols_60,
                                  'ULX_10':ULX_10,
                                  'ULY_10':ULY_10,
                                  'ULX_20':ULX_20,
                                  'ULY_20':ULY_20,
                                  'ULX_60':ULX_60,
                                  'ULY_60':ULY_60,}})
    return metaDict

def readMetadataS2L2A(metadataFile):
    # Get parameters from main metadata file   
    ProductName = os.path.split(os.path.dirname(metadataFile))[1]
    tree = ET.parse(metadataFile)
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
    metaDict = {}
    metaDict.update({'product_name':ProductName,
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
    XML_mask = 'S2?_*.xml'
    for elem in tree.iter(tag='Granules'):
        granule = elem.attrib['granuleIdentifier']
        granuleDir = os.path.join(os.path.dirname(metadataFile), "GRANULE", 
                                  granule)
        globlist = granuleDir + '/' + XML_mask
        metadataTile = glob.glob(globlist)[0]
        # read metadata of tile
        tree = ET.parse(metadataTile)
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
        metaDict.update({granule[len(granule)-13:-7]:{'sun_zenit':sunZen,
                                  'sun_azimuth':sunAz,
                                  'sensor_zenit':sensorZen,
                                  'sensor_azimuth':sensorAz,
                                  'projection':EPSG,
                                  'rows_10':rows_10,
                                  'cols_10':cols_10,
                                  'rows_20':rows_20,
                                  'cols_20':cols_20,
                                  'rows_60':rows_60,
                                  'cols_60':cols_60,
                                  'ULX_10':ULX_10,
                                  'ULY_10':ULY_10,
                                  'ULX_20':ULX_20,
                                  'ULY_20':ULY_20,
                                  'ULX_60':ULX_60,
                                  'ULY_60':ULY_60,}})
    return metaDict