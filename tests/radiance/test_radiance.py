import numpy as np

import atmcorr.radiance

from ..data import MTDFILES


def test_dn_to_radiance():
    dndata = np.random.randint(0, np.iinfo('int16').max, size=(3, 200, 300), dtype='int16')

    for sensor, kwargs in MTDFILES.items():
        rad = atmcorr.radiance.dn_to_radiance(
            dndata=dndata,
            sensor=sensor,
            band_ids=[0, 1, 2],
            **kwargs)
        assert rad.shape == dndata.shape
        assert rad.dtype == np.dtype('float32')
