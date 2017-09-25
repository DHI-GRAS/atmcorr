import dateutil


def get_date_wv(mtdFile):
    timekeys = ['firstLinetime', 'earliestAcqTime']
    with open(mtdFile, 'r') as fin:
        for line in fin:
            linestrip = line.strip()
            for timekey in timekeys:
                if linestrip.startswith(timekey):
                    datestr = linestrip.split('=')[1]
                    datestr = datestr.strip().rstrip(';')
                    return dateutil.parser.parse(datestr)
    raise ValueError('Unable to get date from file \'{}\''.format(mtdFile))
