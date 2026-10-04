"""Physical constants and documented tolerances.

Units: SI (metres, seconds, radians) unless a name explicitly says otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BodyConstants:
    mu: float                # gravitational parameter, m^3/s^2
    radius_equator: float    # equatorial radius, m
    flattening: float        # geometric flattening
    j2: float                # second zonal harmonic
    omega: float             # sidereal rotation rate, rad/s


WGS84 = BodyConstants(
    mu=3.986004418e14,
    radius_equator=6378137.0,
    flattening=1.0 / 298.257223563,
    j2=1.08262668e-3,
    omega=7.2921151467e-5,
)

SUN_ANGULAR_RADIUS_DEG = 0.2666
AU_METRES = 1.495978707e11

DEFAULT_ROOT_TOL_SECONDS = 1e-4

# Documented agreement tolerances (see docs/accuracy.md).
TLE_COE_AGREEMENT_METRES = 1.0
TIME_INPUT_AGREEMENT_US = 1.0
J2_REFERENCE_METRES = 50.0
TWOBODY_ENERGY_DRIFT = 1e-8
TWOBODY_ANGMOM_DRIFT = 1e-10
WINDOW_STABLE_TOL_SECONDS = 1e-3
