from xml.etree import ElementTree as ET

import dateutil.parser


def get_date(mtdFile):
    tree = ET.parse(mtdFile)
    # get down to the appropirate node
    root = tree.getroot()
    Geometric_Data = root.findall('Geometric_Data')[0]
    Use_Area = Geometric_Data.findall('Use_Area')[0]
    for LGV in Use_Area.findall('Located_Geometric_Values'):
        # get angles for centre of the image
        if LGV.findall('LOCATION_TYPE')[0].text == "Center":
            # get year month and day
            datestr = LGV.findall('TIME')[0].text
            return dateutil.parser.parse(datestr)
    raise ValueError('Unable to get date from file \'{}\''.format(mtdFile))
