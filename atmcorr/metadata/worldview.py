from ..wv_meta.parsing import parse_metadata_raw as _parse_metadata_raw
from ..wv_meta.postprocessing import postprocess_metadata


def get_date(mtdFile):
    mtd = parse_metadata(mtdFile)
    return mtd['sensing_time']


def _tastes_like_imd(s):
    return 'BEGIN_GROUP' in s


def parse_metadata_raw(imdfile_or_str):
    """Parse metadata from IMD

    Parameters
    ----------
    imdfile_or_str : str
        path to IMD file
        or string with complete IMD file contents

    Returns
    -------
    dict
        parsed metadata
    """
    if _tastes_like_imd(imdfile_or_str):
        lines = imdfile_or_str.splitlines()
    else:
        with open(imdfile_or_str) as fin:
            lines = fin.readlines()
    return _parse_metadata_raw(lines)


def parse_metadata(imdfile_or_str):
    """Parse metadata from IMD including derived attributes

    Parameters
    ----------
    imdfile_or_str : str
        path to IMD file
        or string with complete IMD file contents

    Returns
    -------
    dict
        parsed metadata
    """
    mtd = parse_metadata_raw(imdfile_or_str)
    mtd_postproc = postprocess_metadata(mtd)
    mtd.update(mtd_postproc)
    return mtd
