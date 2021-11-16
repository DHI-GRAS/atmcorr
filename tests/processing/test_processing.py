import affine
import numpy as np
import rasterio.crs
import pytest

from atmcorr import processing

from .data import MTDFILES


def _get_inputs():
    band_ids = [0, 1, 2]
    nbands = len(band_ids)
    ny, nx = 200, 100
    data = np.random.randint(
            0, np.iinfo('uint16').max,
            size=(nbands, ny, nx),
            dtype='uint16')
    mtdFile = MTDFILES['WV02.imd']
    profile = dict(
            width=nx, height=ny,
            crs=rasterio.crs.CRS('epsg:32640'),
            transform=affine.Affine(
                2.0, 0.0, 364466.0808031342,
                0.0, -2.0, 2836767.9090107735),
            nodata=0,
            count=nbands)
    inputs = dict(
        data=data,
        profile=profile,
        sensor='WV2',
        mtdFile=mtdFile,
        sixs_params=dict(
            aeroProfile='Maritime',
            atm={
                'AOT': 0.35,
                'PWV': 1.0,
                'ozone': 0.15}),
        band_ids=band_ids,
        adjCorr=True,
        aotMultiplier=1.0,
        mtdFile_tile=None,
        date=None,
        use_modis=False,
        modis_atm_dir=None,
        earthdata_credentials={})
    return inputs


@pytest.mark.slow
def test_workflow_wv2_adjCorr():
    inputs = _get_inputs()
    inputs.update(
        adjCorr=True)
    outdata, _ = processing.main(**inputs)
    assert outdata.shape == inputs['data'].shape


@pytest.mark.slow
def test_workflow_wv2_adjCorr_tiled():
    inputs = _get_inputs()
    inputs.update(
        adjCorr=True,
        tileSize=100)
    outdata, _ = processing.main(**inputs)
    assert outdata.shape == inputs['data'].shape
