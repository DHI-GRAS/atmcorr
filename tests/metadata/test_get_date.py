import datetime

import atmcorr.metadata

from ..data import MTDFILES


def test_get_date():
    for sensor, kw in MTDFILES.items():
        date = atmcorr.metadata.get_date(sensor, kw['mtdFile'])
        assert isinstance(date, datetime.datetime)
