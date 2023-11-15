import re
import dateutil.parser

RENAME_BAND_GROUPS = {
    'P': 'PAN',
    'C': 'COASTAL',
    'B': 'BLUE',
    'G': 'GREEN',
    'Y': 'YELLOW',
    'R': 'RED',
    'RE': 'REDEDGE',
    'N': 'NIR1',
    'N2': 'NIR2',
    'S1': 'SWIR1',
    'S2': 'SWIR2',
    'S3': 'SWIR3',
    'S4': 'SWIR4',
    'S5': 'SWIR5',
    'S6': 'SWIR6',
    'S7': 'SWIR7',
    'S8': 'SWIR8',
    }


def parse_metadata_raw(lines):
    """Parse lines from an IMD file into Python types, renaming bands

    Parameters
    ----------
    lines : iterable of str
        lines from IMD file

    Returns
    -------
    dict
        metadata dict with int, float, str, datetime
        and groups 'band_meta', 'image_meta', 'projection_meta'
    """
    r_begin_group = re.compile(r'^BEGIN_GROUP\ =\ (.*)')
    r_end_group = re.compile(r'^END_GROUP\ =\ (.*)')
    r_floats = re.compile(r'\s*(.*?)\s*=\s*([\+\-]?[\d\.eE\-\+]+)')
    r_strings = re.compile(r'\s*(.*?)\s*=\s*"(.*)"\s*')
    r_dates = re.compile(r'\s*(.*?)\s*=\s*(\d{4}\-\d{2}\-\d{2}T[\d\:\.]+[zZ]?)')
    r_integers = re.compile(r'\s*(.*?)\s*=\s*([\+\-]?\d+);')

    root = {}
    band_meta = {}
    image_meta = {}
    projection_meta = {}
    in_group = False
    gname_begin = None
    for line in lines:
        begin_group = r_begin_group.search(line)
        if begin_group is not None:
            gname_begin = begin_group.group(1)
            in_group = True
            continue

        end_group = r_end_group.search(line)
        if end_group is not None:
            gname_end = end_group.group(1)
            if gname_end != gname_begin:
                raise ValueError('{} != {}'.format(gname_begin, gname_end))
            in_group = False
            gname_begin = None
            continue

        if in_group:
            gname = gname_begin
            if gname.startswith('BAND_'):
                key = gname.replace('BAND_', '')
                band = RENAME_BAND_GROUPS[key]
                if band not in band_meta:
                    band_meta[band] = {}
                g = band_meta[band]
            elif gname.startswith('IMAGE_'):
                if gname not in image_meta:
                    image_meta[gname] = {}
                g = image_meta[gname]
            elif gname == 'MAP_PROJECTED_PRODUCT':
                g = projection_meta
            else:
                # skip everything else
                continue
        else:
            g = root

        strings = r_strings.search(line)
        if strings is not None:
            key, value = strings.groups()
            g[key] = value
            continue

        integers = r_integers.search(line)
        if integers is not None:
            key, value = integers.groups()
            g[key] = int(value)
            continue

        dates = r_dates.search(line)
        if dates is not None:
            key, value = dates.groups()
            g[key] = dateutil.parser.parse(value)
            continue

        floats = r_floats.search(line)
        if floats is not None:
            key, value = floats.groups()
            g[key] = float(value)
            continue

    full = root.copy()
    full['band_meta'] = band_meta
    full['image_meta'] = image_meta
    full['projection_meta'] = projection_meta
    return full
