import affine

from atmcorr import tiling


def test_get_tile_corners_ij():
    corners = tiling.get_tile_corners_ij(height=1000, width=2000, xtilesize=100, ytilesize=100)
    assert corners.shape == (2, 4, 10, 20)


def test_get_tile_corners_ij_round():
    corners = tiling.get_tile_corners_ij(height=999, width=1999, xtilesize=100, ytilesize=100)
    assert corners.shape == (2, 4, 10, 20)


def test_transform_corners():
    corners = tiling.get_tile_corners_ij(height=1000, width=2000, xtilesize=100, ytilesize=100)
    lon, lat = tiling.transform_corners(
            corners,
            src_transform=affine.Affine(10, 0, 36000, 0, 10, 18000),
            src_crs={'init': 'epsg:4032'})
    assert lon.shape == (4, 10, 20)


def test_get_tile_extents():
    extent_rec = tiling.get_tile_extents(
            height=1000,
            width=2000,
            xtilesize=100,
            ytilesize=100,
            src_transform=affine.Affine(10, 0, 36000, 0, 10, 18000),
            src_crs={'init': 'epsg:4032'})
    assert extent_rec.shape == (10, 20)
    assert 'xmin' in extent_rec.dtype.names


def test_get_tile_extents_single_tile():
    extent_rec = tiling.get_tile_extents(
            height=1000,
            width=2000,
            xtilesize=2000,
            ytilesize=1000,
            src_transform=affine.Affine(10, 0, 36000, 0, 10, 18000),
            src_crs={'init': 'epsg:4032'})
    assert extent_rec.shape == (1, 1)


def test_corners_from_bounds():
    bounds = (364466.0808031342, 2829605.9090107735, 365034.0808031342, 2836767.9090107735)
    src_crs = {'init': 'epsg:32640'}
    extent_rec = tiling.extents_from_bounds(*bounds, src_crs=src_crs)
    assert extent_rec.shape == ()
    assert 'xmin' in extent_rec.dtype.names
    assert extent_rec['xmin']
