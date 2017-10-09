import re

import dateutil.parser


def get_date(mtdFile):
    """
    DATE_ACQUIRED = 2017-09-09
    SCENE_CENTER_TIME = "10:32:22.4278750Z"
    """
    date_acquired_pattern = re.compile('DATE_ACQUIRED\s=\s([\d-]*)')
    scene_center_time_pattern = re.compile('SCENE_CENTER_TIME\s=\s"([\d:\.]*Z)"')

    datestr = None
    timestr = None
    with open(mtdFile, 'r') as metadata:
        for line in metadata:
            match = date_acquired_pattern.search(line)
            if match is not None:
                datestr = match.group(1)
            match = scene_center_time_pattern.search(line)
            if match is not None:
                timestr = match.group(1)

    if datestr is None or timestr is None:
        raise ValueError(
                'Unable to get datestr ({}) and/or timestr ({}) from MTL file \'{}\'.'
                .format(datestr, timestr, mtdFile))

    fullstr = datestr + 'T' + timestr
    return dateutil.parser.parse(fullstr)
