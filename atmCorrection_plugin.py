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

import os
import sys
import inspect

from processing.core.Processing import Processing
from atmCorrection_provider import AtmosphericCorrectionProvider
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from PyQt4.QtGui import QAction, QIcon
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from atmospheric_correction_dialog import atmCorrectionDialog
import os.path


cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]

if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)


class AtmosphericCorrectionPlugin:

    def __init__(self, iface):
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'atmCorrection_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = atmCorrectionDialog()

        # get reference to the QGIS message bar
        self.msg_bar = self.iface.messageBar()

        # Create the provider for Processing Toolbox
        self.provider = AtmosphericCorrectionProvider()


    def initGui(self):
        """
        Create action that will start plugin configuration
        """
        self.action = QAction(
            QIcon(":/plugins/atmCorrection/icon.png"),
            u"Atmospheric Correction", self.iface.mainWindow())
        # connect the action to the run method
        self.action.triggered.connect(self.run)

        # Add toolbar button and menu item
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToVectorMenu(u"&Atmospheric Correction", self.action)

        # Add algorithms to Processing Toolbox
        Processing.addProvider(self.provider)

    def unload(self):
        """
        Remove the plugin menu item and icon
        """
        self.iface.removePluginVectorMenu(u"&Atmospheric Correction", self.action)
        self.iface.removeToolBarIcon(self.action)
        Processing.removeProvider(self.provider)
        Processing.removeProvider(self.provider)

    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass