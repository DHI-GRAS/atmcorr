import logging

logger = logging.getLogger(__name__)


def toa_reflectance_S2(data, metadata):
    """Convert to TOA reflectance

    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    rc = float(metadata['reflection_conversion'])
    logger.info("TOA reflectance")
    reflectance = data.astype('f4')
    reflectance /= rc
    return reflectance
