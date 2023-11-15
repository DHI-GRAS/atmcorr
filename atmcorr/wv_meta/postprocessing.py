import re

import shapely.geometry
import affine


def _get_angles(image_meta):
    angles = {}
    for key in image_meta:
        if re.match(r'(?:min|max|mean)(?:Sun|Sat|InTrack|CrossTrack|OffNadir)', key) is not None:
            angles[key] = image_meta[key]
    return angles


def _get_calibration_constants(band_meta):
    calibration_constants = {}
    for key in ['absCalFactor', 'effectiveBandwidth']:
        calibration_constants[key] = {band: band_meta[band][key] for band in band_meta}
    return calibration_constants


def _get_points_lonlat(cdict):
    points = []
    for sj, si in zip('UULL', 'LRRL'):
        sbase = sj + si
        xkey = sbase + 'Lon'
        ykey = sbase + 'Lat'
        points.append((cdict[xkey], cdict[ykey]))
    assert len(points) == 4
    points.append(points[0])
    return points


def _get_points_xy(cdict):
    points = []
    for sj, si in zip('UULL', 'LRRL'):
        sbase = sj + si
        xkey = sbase + 'X'
        ykey = sbase + 'Y'
        points.append((cdict[xkey], cdict[ykey]))
    assert len(points) == 4
    points.append(points[0])
    return points


def _points_to_polygon(points):
    return shapely.geometry.Polygon(points)


def _get_transform(projection_meta):
    coefs = (
        projection_meta[key]
        for key in [
                'colSpacing', 'orientationAngle', 'originX',
                'orientationAngle', 'rowSpacing', 'originY'])
    a, b, c, d, e, f = coefs
    return affine.Affine(a, b, c, d, -e, f)


def postprocess_metadata(mtd_raw):
    """Further derive types from raw metadata

    Parameters
    ----------
    mtd_raw : dict
        metadata parsed with parsing.parse_metadata_raw

    Returns
    -------
    dict
        postprocessed metadata
    """
    mtd_postproc = {}
    bm = next(iter(mtd_raw['band_meta'].values()))
    imgm = next(iter(mtd_raw['image_meta'].values()))
    copy_imgm = ['satId']
    for key in copy_imgm:
        mtd_postproc[key] = imgm[key]
    mtd_postproc['angles'] = _get_angles(imgm)
    mtd_postproc['footprint'] = _points_to_polygon(_get_points_lonlat(bm))
    mtd_postproc['sensing_time'] = (
        imgm['firstLineTime'] if 'firstLineTime' in imgm
        else mtd_raw['earliestAcqTime'])
    mtd_postproc['calibration'] = _get_calibration_constants(mtd_raw['band_meta'])
    mtd_postproc['transform'] = _get_transform(mtd_raw['projection_meta'])
    mtd_postproc['footprint_projected'] = _points_to_polygon(
        _get_points_xy(mtd_raw['projection_meta']))
    return mtd_postproc
