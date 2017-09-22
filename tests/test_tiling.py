import affine

from atmospheric_correction import tiling


def test_get_tile_corners():
    corners = tiling.get_tile_corners(height=1000, width=2000, tilesize=100)
    assert corners.shape == (2, 4, 10, 20)


def test_get_tile_corners_round():
    corners = tiling.get_tile_corners(height=999, width=1999, tilesize=100)
    assert corners.shape == (2, 4, 10, 20)


def test_transform_corners():
    corners = tiling.get_tile_corners(height=1000, width=2000, tilesize=100)
    lon, lat = tiling.transform_corners(
            corners,
            src_transform=affine.Affine(10, 0, 36000, 0, 10, 18000),
            src_crs={'init': 'epsg:4032'})
    assert lon.shape == (4, 10, 20)


def test_get_tile_extents():
    extent_rec = tiling.get_tile_extents(
            height=1000,
            width=2000,
            tilesize=100,
            src_transform=affine.Affine(10, 0, 36000, 0, 10, 18000),
            src_crs={'init': 'epsg:4032'})
    assert extent_rec.shape == (10, 20)
    assert 'xmin' in extent_rec.dtype.names
