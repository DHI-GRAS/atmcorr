from atmcorr.landsat8 import unzip


def get_bandfiles(infiles, bands):
    if isinstance(infiles, (list, tuple)):
        if not len(infiles) == len(bands):
            raise ValueError(
                    'Number of band files {} does not match number of bands {}'
                    .format(infiles, bands))
        return infiles
    else:
        # assume str
        tarfile = infiles
        return unzip.get_bandfile_urls(tarfile, bands)
