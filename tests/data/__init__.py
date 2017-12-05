import os

HERE = os.path.abspath(os.path.dirname(__file__))

MTDFILES = {
    'WV2': {'mtdFile': '13SEP20073259-M2AS-055783286010_01_P001.imd'},
    'PHR1A': {'mtdFile': 'DIM_PHR1A_PMS_201710070757436_PRJ_2540070101-001.XML'},
    'S2A': {
        'mtdFile': 'MTD_MSIL1C.xml',
        'mtdFile_tile': 'MTD_TL.xml'}}


def _join_paths(d):
    for key, path in d.items():
        d[key] = os.path.join(HERE, path)


[_join_paths(d) for d in MTDFILES.values()]
