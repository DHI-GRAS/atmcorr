# -*- coding: utf-8 -*-

"""
/***************************************************************************
 AtmosphericCorrection
                                 A QGIS plugin
 Use 6S module to perform atmospheric correction on satellite imagery
                              -------------------
        begin                : 2015-12-17
        copyright            : (C) 2015 by DHI-GRAS
        email                : rfn@dhi-gras.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'DHI-GRAS'
__date__ = '2015-12-17'
__copyright__ = '(C) 2015 by DHI-GRAS'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from processing.core.GeoAlgorithm import GeoAlgorithm
from processing.core.outputs import OutputRaster
from processing.core.parameters import ParameterSelection
from processing.core.parameters import ParameterBoolean
from processing.core.parameters import ParameterRaster
from processing.core.parameters import ParameterFile
from processing.core.parameters import ParameterNumber

from atmProcessing import atmProcessingMain

class AtmosphericCorrectionAlgorithm(GeoAlgorithm):

    SATELLITE = 'SATELLITE'
    SATELLITES = ['WorldView-2', 'WorldView-3', 'Landsat-8', 'Landsat-7',
                  'Pleiades-A', 'Pleiades-B', 'SPOT-6']
    DN_FILE = 'DN_FILE'
    METAFILE = 'METAFILE'
    PANCHROMATIC = 'PANCHROMATIC'
    METHOD = 'METHOD'
    METHODS = ['6S', 'DOS', 'TOA', 'RAD']
    ATMOSPHERIC_PROFILE = 'ATMOSPHERIC_PROFILE'
    ATMOSPHERIC_PROFILES = ['No aerosols', 'Continental', 'Maritime', 'Urban',
                            'Desert', 'Biomass burning', 'Stratospheric']
    DOWNLOAD_MODIS = 'DOWNLOAD_MODIS'
    AOT = 'AOT'
    PWV = 'PWV'
    OZONE = 'OZONE'
    OUTPUT_FILE = 'OUTPUT_FILE' 


    def defineCharacteristics(self):
        # The name that the user will see in the toolbox
        self.name = 'Atmospheric correction'
        # The branch of the toolbox under which the algorithm will appear
        self.group = 'Tools'
        
        self.addParameter(ParameterSelection(self.SATELLITE, 'Satellite', self.SATELLITES))
        self.addParameter(ParameterRaster(self.DN_FILE, 'DN file', showSublayersDialog=False))
        self.addParameter(ParameterFile(self.METAFILE, 'Metafile', optional=False))
        self.addParameter(ParameterBoolean(self.PANCHROMATIC, 'Panchromatic', default=False))
        self.addParameter(ParameterSelection(self.METHOD, 'Method', self.METHODS))
        self.addParameter(ParameterSelection(self.ATMOSPHERIC_PROFILE, 'Atmospheric profile', self.ATMOSPHERIC_PROFILES))
#        self.addParameter(ParameterBoolean(self.DOWNLOAD_MODIS, 'download modis'))
        self.addParameter(ParameterNumber(self.AOT, 'aot'))
        self.addParameter(ParameterNumber(self.PWV, 'pwv'))
        self.addParameter(ParameterNumber(self.OZONE, 'ozone'))
        self.addOutput(OutputRaster(self.OUTPUT_FILE, 'output file'))

    def processAlgorithm(self, progress):
        """Here is where the processing itself takes place."""
        # The first thing to do is retrieve the values of the parameters
        # entered by the user
        
        sensorList = ["WV2", "WV3", "L8", "L7", "PHR1A", "PHR1B", "SPOT6"]
        methodList = ["6S", "DOS", "TOA", "RAD"]
        atmProfList = ["No Aerosols", "Continental", "Maritime", "Urban", "Desert", "BiomassBurning", "Stratospheric"]
        
        options = {}
        # input/output parameters
        options["sensor"] = sensorList[self.getParameterValue(self.SATELLITE)]
        options["dnFile"] = self.getParameterValue(self.DN_FILE)
        options["metadataFile"] = self.getParameterValue(self.METAFILE)
        options["reflectanceFile"] = self.getOutputValue(self.OUTPUT_FILE)
        # Atmospheric correction parameters
        options["atmCorrMethod"] = methodList[self.getParameterValue(self.METHOD)]
        options["atm"] = self.atmParam()
        options["isPan"] = self.getParameterValue(self.PANCHROMATIC)
        options["aeroProfile"] = atmProfList[self.getParameterValue(self.ATMOSPHERIC_PROFILE)]
        
        atmProcessingMain(options)

    def atmParam(self):
        atm = {'AOT': self.getParameterValue(self.AOT), 'PWV': self.getParameterValue(self.PWV), 'ozone': self.getParameterValue(self.OZONE)}
        return atm
           
#    def getCustomParametersDialog(self):
#        return FieldsCalculatorDialog(self)
