import datetime

import atmcorr.metadata

from .data import METADATA


def test_get_date():
    for sensor, kw in METADATA.items():
        date = atmcorr.metadata.get_date(sensor, kw['mtdFile'])
        assert isinstance(date, datetime.datetime)
