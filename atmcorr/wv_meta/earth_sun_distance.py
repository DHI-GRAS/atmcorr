from __future__ import division
import numpy as np


def get_earth_sun_distance(sensing_time):
    """Get actual Earth-Sun distance following equations from
       Radiometric Use Of WorldView-2 Imagery - Technical note

    Parameters
    ----------
    sensing_time : datetime.datetime
        sensing date-time

    Returns
    -------
    float
        Earth-Sun distance in m
    """
    year = sensing_time.year
    month = sensing_time.month
    day = sensing_time.day
    hours = sensing_time.hour
    if month < 3:
        year -= 1
        month += 12
    A = int(year / 100.0)
    B = 2.0 - A + int(A / 4.0)
    julian_day = (
            int(365.25 * (year + 4716.0)) +
            int(30.6001 * (month + 1)) +
            day + hours / 24.0 + B - 1524.5)
    D = julian_day - 2451545.0
    g = 357.529 + 0.98560025 * D
    des = 1.00014 - 0.01671 * np.cos(np.radians(g)) - 0.00014 * np.cos(np.radians(2 * g))
    return des
