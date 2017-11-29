import satmeta.l8.meta as l8meta


def get_date(mtdFile):
    with open(mtdFile) as fin:
        metadata = l8meta.parse_metadata(fin)
    return metadata['sensing_time']
