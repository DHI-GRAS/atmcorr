

def resample(
        source, src_transform, src_crs, dst_shape,
        dst_transform=None, dst_crs=None,
        resampling=None, **reprojectkw):
    """Resample data

    Parameters
    ----------
    source : ndarray
        source data
    src_transform : affine.Affine
        source transformation
    src_crs : dict or rasterio.crs.CRS
        source coordinate reference system
    dst_shape : tuple
        output shape
    dst_transform : affine.Affine, optional
        destination transform
        default: same as src_transform
    dst_crs : dict or rasterio.crs.CRS, optional
        destination coordinate reference system
        default: same as src_crs
    resampling : int
        resampling method
        see rasterio.warp.Resampling
        default: Resampling.bilinear
    **reprojectkw : additional keyword arguments
        passed to rasterio.warp.reproject

    Returns
    -------
    ndarray : resampled data
    """
    import rasterio.warp
    import numpy as np

    if dst_transform is None:
        dst_transform = src_transform
    if dst_crs is None:
        dst_crs = src_crs

    if resampling is None:
        resampling = rasterio.warp.Resampling.bilinear

    destination = np.zeros(dst_shape, 'f8')
    rasterio.warp.reproject(
        source=source,
        destination=destination,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        resampling=resampling,
        **reprojectkw)
    return destination
