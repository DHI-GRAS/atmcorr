import re
import datetime


def get_date(mtdFile):
    DATEACQUIREDRegex = "\s*DATE_ACQUIRED\s*=\s*(\d{4})[-_](\d{2})[-_](\d{2})"
    with open(mtdFile, 'r') as metadata:
        for line in metadata:
            match = re.match(DATEACQUIREDRegex, line)
            if match is None:
                continue
            year = int(match.group(1))
            month = int(match.group(2))
            day = int(match.group(3))
            return datetime.date(year, month, day)
    raise ValueError('Unable to get date from file \'{}\''.format(mtdFile))
