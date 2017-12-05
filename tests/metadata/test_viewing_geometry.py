import atmcorr.viewing_geometry as vg

from .data import METADATA


def test_viewing_geometry():
    for sensor, kw in METADATA.items():
        gdict = vg.get_geometry(sensor=sensor, **kw)
        assert isinstance(gdict, dict)
        assert 'sun_zenith' in gdict
