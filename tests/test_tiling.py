import affine

from atmospheric_correction import tiling


def test_get_tiled_windows():
    windows = tiling.get_tiled_windows(height=1000, width=2000, tilesize=10)
    assert windows.shape == (100, 200)


def test_get_tiled_windows_round():
    windows = tiling.get_tiled_windows(height=999, width=1999, tilesize=10)
    assert windows.shape == (100, 200)


def test_get_transformed_bounds():
    windows = tiling.get_tiled_windows(height=1000, width=2000, tilesize=10)
    bb = tiling.get_transformed_bounds(
            windows,
            src_transform=affine.Affine(10, 0, 36000, 0, 10, 18000),
            src_crs={'init': 'epsg:4032'},
            dst_crs={'init': 'epsg:4326'})

    assert bb.shape == (4, ) + windows.shape
