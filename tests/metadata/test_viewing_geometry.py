import atmcorr.viewing_geometry as vg

from ..data import MTDFILES


def test_viewing_geometry():
    for sensor, kw in MTDFILES.items():
        gdict = vg.get_geometry(sensor=sensor, **kw)
        assert isinstance(gdict, dict)
        assert 'sun_zenith' in gdict
