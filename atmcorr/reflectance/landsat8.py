from rio_toa import reflectance

from atmcorr.landsat8 import utils as l8utils


def toa_reflectance(infiles, outfile, mtdfile, bands):
    """Calculate TOA reflectance

    Parameters
    ----------
    infiles : list of str or str
        paths to Landsat 8 band files
        or URIs for members in TAR file
        or path to TAR file
    outfile : str
        path to save output to
    mtdfile : str
        path to metadata file
    bands : list of int
        bands to extract from TAR file
        or bands that the URIs correspond to
    """
    bandfiles = l8utils.get_bandfiles(infiles, bands)
    reflectance.calculate_landsat_reflectance(
            src_path=bandfiles,
            src_mtl=mtdfile,
            dst_path=outfile,
            rescale_factor=None,
            creation_options={},
            bands=bands,
            dst_dtype='float32',
            processes=1,
            pixel_sunangle=True)
