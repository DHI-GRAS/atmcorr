import affine
import numpy as np
import rasterio.crs
import pytest

from atmcorr import processing

from . import testdata


@pytest.mark.slow
def test_workflow_wv2():
    band_ids = [0, 1, 2]
    nbands = len(band_ids)
    ny, nx = 200, 100
    data = np.random.randint(
            0, np.iinfo('uint16').max,
            size=(nbands, ny, nx),
            dtype='uint16')

    mtdFile = testdata.mtdfiles['WV2']

    profile = dict(
            width=nx, height=ny,
            crs=rasterio.crs.CRS({'init': 'epsg:32640'}),
            transform=affine.Affine(
                2.0, 0.0, 364466.0808031342,
                0.0, -2.0, 2836767.9090107735),
            nodata=0,
            count=nbands)

    outdata = processing.main(
            sensor='WV2',
            mtdFile=mtdFile,
            method='6S',
            atm={
                'AOT': 0.35,
                'PWV': 1.0,
                'ozone': 0.15},
            aeroProfile='Maritime',
            dnFile=None,
            data=data,
            profile=profile,
            tileSizePixels=0,
            band_ids=band_ids,
            adjCorr=True,
            aotMultiplier=1.0,
            nprocs=None,
            mtdFile_tile=None,
            date=None,
            outfile=None,
            return_profile=False,
            use_modis=False,
            modis_atm_dir=None,
            earthdata_credentials={})

    assert outdata.shape == data.shape
