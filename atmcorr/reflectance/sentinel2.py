
def toa_reflectance(data, metadata):
    """Convert to TOA reflectance

    Assumes a L1C product which contains TOA reflectance:
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-2-msi/product-types
    """
    rc = float(metadata['reflection_conversion'])
    reflectance = data.astype('f4')
    reflectance /= rc
    return reflectance
