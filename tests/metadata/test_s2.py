import atmospheric_correction.metadata.s2 as s2meta

from . import testdata_s2l1c as testdata


def test_parse_mtdfile():
    metadata = s2meta.parse_mtdfile(testdata.mtdfile)
    assert 'reflection_conversion' in metadata
    assert 'irradiance_values' in metadata
    assert 'quantification_value' in metadata


def test_parse_granule_mtdfile():
    gmeta = s2meta.parse_granule_mtdfile(testdata.mtdfile_tile)
    tile = list(gmeta)[0]
    assert tile == testdata.tile
    assert 'sun_zenith' in gmeta[tile]
    assert isinstance(gmeta[tile]['sensor_azimuth'], float)
