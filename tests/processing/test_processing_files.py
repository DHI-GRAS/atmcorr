import pytest
import rasterio
import numpy as np

try:
    from atmcorr_testdata import wv2 as wv2data
except ImportError:
    pytestmark = pytest.mark.skip(reason='atmcorr_testdata not installed')

from atmcorr import processing
import ruamel.yaml


@pytest.mark.slow
def test_processing_wv2():

    config = ruamel.yaml.safe_load(open(wv2data.DATAFILES['config.yaml']))
    config.update(
        mtdFile=wv2data.DATAFILES['mtdfile.imd'])

    with rasterio.open(wv2data.DATAFILES['dn_allbands.tif']) as src:
        data = src.read()
        profile = src.profile
    new, _ = processing.main(data=data, profile=profile, **config)

    with rasterio.open(wv2data.DATAFILES['out_expected.tif']) as src:
        expected = src.read()

    assert new.shape == expected.shape

    assert not np.all(np.isnan(new))
    assert not np.all(np.isnan(expected))

    assert np.nanmax(np.abs(new - expected)) < 1e-6
