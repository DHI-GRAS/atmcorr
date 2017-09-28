import tarfile
import fnmatch


def _get_band_fnpattern(band):
    return '*_B{}.TIF'.format(band)


def _find_band_file(names, band):
    fnpattern = _get_band_fnpattern(band)
    try:
        return fnmatch.filter(names, fnpattern)[0]
    except IndexError:
        raise RuntimeError('Unable to find file for band {} in list {}'.format(band, names))


def _get_names_in_file(infile):
    with tarfile.open(infile) as tar:
        return tar.getnames()


def _generate_member_url(infile, memberpath):
    return 'tar://' + infile + '!' + memberpath


def get_bandfile_urls(infile, bands):
    """Get URLs to band files in TAR archive

    Parameters
    ----------
    infile : str
        path to .tar(.gz) file
    bands : list of str
        bands to get
    """
    names = _get_names_in_file(infile)
    urls = []
    for band in bands:
        bfpath = _find_band_file(names, band=band)
        url = _generate_member_url(infile, bfpath)
        urls.append(url)
    return urls
