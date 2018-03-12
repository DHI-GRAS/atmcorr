import affine
import numpy as np

from atmcorr import tiling


def test_get_tiled_transform_shape():
    src_res = 10
    src_transform = affine.Affine(src_res, 0, 36000, 0, -src_res, 18000)
    src_shape = (1000, 2000)
    dst_res = src_res * 10
    dst_transform, dst_shape = tiling.get_tiled_transform_shape(src_transform, src_shape, dst_res)
    expected_transform = src_transform * affine.Affine(10, 0, 0, 0, 10, 0)
    assert dst_transform == expected_transform
    assert dst_shape == (100, 200)


def test_get_projected_extents():
    extent_rec = tiling.get_projected_extents(
        transform=affine.Affine(10, 0, 36000, 0, -10, 18000),
        height=100,
        width=200,
        src_crs={'init': 'epsg:4032'})
    assert extent_rec.shape == (100, 200)
    assert 'xmin' in extent_rec.dtype.names


def test_bounds_to_projected_extents():
    bounds = (364466.0808031342, 2829605.9090107735, 365034.0808031342, 2836767.9090107735)
    src_crs = {'init': 'epsg:32640'}
    extent_rec = tiling.bounds_to_projected_extents(*bounds, src_crs=src_crs)
    assert extent_rec.shape == (1, 1)
    assert 'xmin' in extent_rec.dtype.names
    assert np.issubdtype(extent_rec['xmin'][0, 0], np.floating)


def test_recarr_take_dict():
    a = tiling.get_projected_extents(
        transform=affine.Affine(10, 0, 36000, 0, -10, 18000),
        height=100,
        width=200,
        src_crs={'init': 'epsg:4032'})
    d = tiling.recarr_take_dict(a, 0, 0)
    assert type(d) == dict
    assert np.issubdtype(d['xmin'], np.floating)
