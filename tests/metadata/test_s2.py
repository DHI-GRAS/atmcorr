import atmospheric_correction.metadata.s2 as s2meta

from . import testdata_s2l1c as testdata


def test_parse_mtdfile():
    metadata = s2meta.parse_mtdfile(testdata.mtdfile)
    assert 'irradiance_values' in metadata
    assert 'quantification_value' in metadata
