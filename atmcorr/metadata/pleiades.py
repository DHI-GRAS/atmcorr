from satmeta.pleiades import parser as plparser


def get_date(mtdFile):
    metadata = plparser.parse_metadata(mtdFile)
    return metadata['sensing_time']
