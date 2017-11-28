from dg_calibration import metadata


def get_date(mtdFile):
    mtd = metadata.parse_metadata(mtdFile)
    return mtd['sensing_time']
