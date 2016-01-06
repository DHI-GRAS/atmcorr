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

from processing.core.AlgorithmProvider import AlgorithmProvider
from processing.core.ProcessingConfig import Setting, ProcessingConfig
from atmCorrection_algorithm import AtmosphericCorrectionAlgorithm
from PyQt4.QtGui import QIcon
import os

class AtmosphericCorrectionProvider(AlgorithmProvider):

#    MY_DUMMY_SETTING = 'MY_DUMMY_SETTING'

    def __init__(self):
        AlgorithmProvider.__init__(self)

        # Deactivate provider by default
#        self.activate = False

        # Load algorithms
        self.alglist = [AtmosphericCorrectionAlgorithm()]
        for alg in self.alglist:
            alg.provider = self

#    def initializeSettings(self):
#        """In this method we add settings needed to configure our
#        provider.
#
#        Do not forget to call the parent method, since it takes care
#        or automatically adding a setting for activating or
#        deactivating the algorithms in the provider.
#        """
#        AlgorithmProvider.initializeSettings(self)
#        #ProcessingConfig.addSetting(Setting('Atmospheric correction - 1',
#        #    AtmosphericCorrectionProvider.MY_DUMMY_SETTING,
#        #    'Example setting', 'Default value'))
#        ProcessingConfig.addSetting(Setting(self.getDescription(), "atmCorrection", "BEAM install directory", "default_value"))

    def unload(self):
        """Setting should be removed here, so they do not appear anymore
        when the plugin is unloaded.
        """
        AlgorithmProvider.unload(self)
#        ProcessingConfig.removeSetting(
#            AtmosphericCorrectionProvider.MY_DUMMY_SETTING)

    def getName(self):
        """This is the name that will appear on the toolbox group.

        It is also used to create the command line name of all the
        algorithms from this provider.
        """
        return 'DHIGRAS'

    def getDescription(self):
        """This is the provired full name.
        """
        return 'DHI-GRAS algorithms'

    def getIcon(self):
        """We return the default icon.
        """
        #return AlgorithmProvider.getIcon(self)
        return QIcon(os.path.dirname(__file__) + "/icon.png")

    def _loadAlgorithms(self):
        """Here we fill the list of algorithms in self.algs.

        This method is called whenever the list of algorithms should
        be updated. If the list of algorithms can change (for instance,
        if it contains algorithms from user-defined scripts and a new
        script might have been added), you should create the list again
        here.

        In this case, since the list is always the same, we assign from
        the pre-made list. This assignment has to be done in this method
        even if the list does not change, since the self.algs list is
        cleared before calling this method.
        """
        self.algs = self.alglist
