import satmeta.s2.meta as s2meta


def get_date(mtdFile):
    meta = s2meta.parse_metadata(mtdFile)
    return meta['sensing_time']
