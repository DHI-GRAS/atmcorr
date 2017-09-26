import satmeta.s2.meta as s2meta


def parse_mtdfile(mtdFile, mtdFile_tile=None):
    metadict = s2meta.parse_metadata(mtdFile)
    if mtdFile_tile is not None:
        gmeta = s2meta.parse_granule_metadata(mtdFile_tile)
        metadict['granules'] = {gmeta['tile_name']: gmeta}
    return metadict


def get_date(mtdFile):
    metadata = parse_mtdfile(mtdFile)
    return metadata['sensing_time']
