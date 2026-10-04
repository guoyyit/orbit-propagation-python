"""Reference-frame transformations: J2000 inertial <-> ECEF, topocentric.

Inertial frame: J2000 (EME2000).  Earth-fixed positions use an ECEF frame
realised by a single z-axis rotation through Greenwich Mean Sidereal Time;
precession/nutation/polar motion are neglected (docs/conventions.md).
Ground stations use WGS-84 geodetic coordinates.
"""
from __future__ import annotations

import math
from datetime import datetime

from .elements import EARTH_FLATTENING, EARTH_RADIUS_M
from .timeutil import datetime_to_jd, greenwich_sidereal_time


def _gmst(when: datetime) -> float:
    return greenwich_sidereal_time(datetime_to_jd(when))


def j2000_to_ecef(
    position: tuple[float, float, float],
    when: datetime,
) -> tuple[float, float, float]:
    """Rotate a J2000 position into ECEF (z-axis rotation by GMST)."""
    gmst = _gmst(when)
    c = math.cos(gmst)
    s = math.sin(gmst)
    x, y, z = position
    return (c * x + s * y, -s * x + c * y, z)


def ecef_to_j2000(
    position: tuple[float, float, float],
    when: datetime,
) -> tuple[float, float, float]:
    """Rotate an ECEF position into J2000 (inverse GMST rotation)."""
    gmst = _gmst(when)
    c = math.cos(gmst)
    s = math.sin(gmst)
    x, y, z = position
    return (c * x - s * y, s * x + c * y, z)


def geodetic_to_ecef(
    lat_rad: float,
    lon_rad: float,
    alt_m: float,
) -> tuple[float, float, float]:
    """Geodetic (WGS-84) latitude/longitude/altitude -> ECEF metres."""
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    e2 = EARTH_FLATTENING * (2.0 - EARTH_FLATTENING)
    prime_vertical = EARTH_RADIUS_M / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    x = (prime_vertical + alt_m) * cos_lat * math.cos(lon_rad)
    y = (prime_vertical + alt_m) * cos_lat * math.sin(lon_rad)
    z = (prime_vertical * (1.0 - e2) + alt_m) * sin_lat
    return (x, y, z)


def ecef_to_enu_rotation(lat_rad: float, lon_rad: float) -> list[list[float]]:
    """3x3 rotation matrix mapping ECEF relative vectors to local ENU."""
    sl = math.sin(lat_rad)
    cl = math.cos(lat_rad)
    so = math.sin(lon_rad)
    co = math.cos(lon_rad)
    return [
        [-so, co, 0.0],
        [-sl * co, -sl * so, cl],
        [cl * co, cl * so, sl],
    ]


def ecef_to_enu(
    delta_ecef: tuple[float, float, float],
    lat_rad: float,
    lon_rad: float,
) -> tuple[float, float, float]:
    """Rotate an ECEF relative vector into local ENU components."""
    rot = ecef_to_enu_rotation(lat_rad, lon_rad)
    return tuple(
        sum(rot[row][col] * delta_ecef[col] for col in range(3))
        for row in range(3)
    )
