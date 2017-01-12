# -*- coding: utf-8 -*-
"""
/***************************************************************************
 atmCorrectionDialog
                                 A QGIS plugin
 Use 6S module to perform atmospheric correction on satellite imagery
                             -------------------
        begin                : 2015-12-09
        git sha              : $Format:%H$
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

import os
from PyQt4 import QtGui, uic
from atmProcessing import atmProcessingMain
from bathyUtilities import saveImgByCopy

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'atmospheric_correction_dialog_base.ui'))

class atmCorrectionDialog(QtGui.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(atmCorrectionDialog, self).__init__(parent)
        self.setupUi(self)

        self.toolButton_DN.clicked.connect(self.selectDN)
        self.toolButton_meta.clicked.connect(self.selectMeta)
        self.checkBox_modis.stateChanged.connect(self.modisCheckBox)
        self.toolButton_output.clicked.connect(self.selectOutput)
        self.pushButton_cancel.clicked.connect(self.closeWindow)
        self.pushButton_save.clicked.connect(self.runAtmCorrection)
        self.comboBox_method.currentIndexChanged.connect(self.methodChange)

    def selectDN(self):
        self.lineEdit_DN.setText(QtGui.QFileDialog.getOpenFileName(self,
                    "Select DN file", ""))

    def selectMeta(self):
        self.lineEdit_meta.setText(QtGui.QFileDialog.getOpenFileName(self,
                    "Select metadata file", ""))

    def selectOutput(self):
        self.lineEdit_output.setText(QtGui.QFileDialog.getSaveFileName(self,
                    "Save corrected file", "", 'Image (*.tif)'))

    def closeWindow(self):
        self.close()

    def satellite(self):
        index = self.comboBox_satellite.currentIndex()
        sensorList = ["WV2", "WV3", "L8", "L7", "PHR1A", "PHR1B", "SPOT6", "S2A_10m", "S2A_60m"]
        return sensorList[index]

    def method(self):
        index = self.comboBox_method.currentIndex()
        methodList = ["DOS", "TOA", "RAD"]
        return methodList[index]

    def atmProfile(self):
        index = self.comboBox_atmProfile.currentIndex()
        atmProfList = ["No Aerosols", "Continental", "Maritime", "Urban", "Desert", "BiomassBurning", "Stratospheric"]
        return atmProfList[index]

    def atmParamSwitch(self, switch):
        if switch == "on":
            self.lineEdit_aot.setEnabled(True)
            self.lineEdit_pwv.setEnabled(True)
            self.lineEdit_ozone.setEnabled(True)
            self.label_Ozone.setEnabled(True)
            self.label_PWV.setEnabled(True)
            self.label_AOT.setEnabled(True)
        elif switch == "off":
            self.lineEdit_aot.setEnabled(False)
            self.lineEdit_pwv.setEnabled(False)
            self.lineEdit_ozone.setEnabled(False)
            self.label_Ozone.setEnabled(False)
            self.label_PWV.setEnabled(False)
            self.label_AOT.setEnabled(False)

    def modisCheckBox(self):
        if self.checkBox_modis.isChecked():
            self.atmParamSwitch("off")
        else:
            self.atmParamSwitch("on")

    def methodChange(self):
        self.atmParamSwitch("off")
        self.label_atmParameters.setEnabled(False)
        self.label_atmProfile.setEnabled(False)
        self.comboBox_atmProfile.setEnabled(False)
        self.checkBox_modis.setEnabled(False)


    def atmParam(self):
        if self.checkBox_modis.isChecked():
            atm = {'AOT': "", 'PWV': "", 'ozone': ""}
        else:
            atm = {'AOT': float(self.lineEdit_aot.text()), 'PWV': float(self.lineEdit_pwv.text()), 'ozone': float(self.lineEdit_ozone.text())}
        return atm

    def isPan(self):
        if self.checkBox_ispan.isChecked():
            isPanbool = True
        else:
            isPanbool = False
        return isPanbool

    def tileSizePixels(self):
        return self.spinBox_tileSize.value()

    def runAtmCorrection(self):
        options = {}
        # input/output parameters
        options["sensor"] = self.satellite()
        options["dnFile"] = self.lineEdit_DN.text()
        options["metadataFile"] = self.lineEdit_meta.text()
        options["reflectanceFile"] = self.lineEdit_output.text()
        # Atmospheric correction parameters
        options["atmCorrMethod"] = self.method()
        options["atm"] = self.atmParam()
        options["isPan"] = self.isPan()
        options["aeroProfile"] = self.atmProfile()
        options["adjCorr"] = 0
        options["tileSizePixels"]=self.tileSizePixels()

        reflectanceImg = atmProcessingMain(options)
        saveImgByCopy(reflectanceImg, options["reflectanceFile"])
        self.closeWindow()
