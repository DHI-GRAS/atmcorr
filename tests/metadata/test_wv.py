import atmospheric_correction.viewing_geometry as vg

from . import testdata_wv as testdata


def test_viewing_geometry():
    vg.get_geometry_WV2(testdata.mtdfile)
