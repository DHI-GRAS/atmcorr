# -*- coding: utf-8 -*-
"""
/***************************************************************************
 atmCorrection
                                 A QGIS plugin
 Use 6S module to perform atmospheric correction on satellite imagery
                             -------------------
        begin                : 2015-12-09
        copyright            : (C) 2015 by DHI-GRAS
        email                : rfn@dhi-gras.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load atmCorrection class from file atmCorrection.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .atmCorrection_plugin import AtmosphericCorrectionPlugin
    return AtmosphericCorrectionPlugin(iface)
