import satmeta.s2.meta as s2meta


def parse_mtdfile(mtdFile, mtdFile_tile=None):
    meta = s2meta.parse_metadata(mtdFile)
    if mtdFile_tile is not None:
        gmeta = s2meta.parse_granule_metadata(mtdFile_tile)
        meta['granules'] = {gmeta['tile_name']: gmeta}
    return meta


def get_date(mtdFile):
    meta = s2meta.parse_metadata(mtdFile)
    return meta['sensing_time']
