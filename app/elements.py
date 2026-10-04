"""Classical orbital elements and Kepler/state conversions.

Angles are in radians, lengths in metres, times in seconds.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime

from .errors import ErrorCode, OrbitError


# WGS-84 / EGM96 constants used by the service.
EARTH_MU = 3.986004418e14  # m^3 / s^2
EARTH_RADIUS_M = 6_378_137.0  # equatorial radius, m
EARTH_FLATTENING = 1.0 / 298.257223563
EARTH_J2 = 1.08262668e-3
EARTH_ROTATION_RATE = 7.2921151467e-5  # rad/s (WGS-84 omega_e)
# WGS-72 constants, used internally for TLE conversions.
WGS72_MU = 3.986008e14
WGS72_RADIUS_M = 6_378_135.0

TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class ClassicalElements:
    """Classical (Keplerian) elements referenced to the epoch.

    Frame: J2000 (EME2000) geocentric inertial.
    """

    epoch: datetime
    semi_major_axis_m: float
    eccentricity: float
    inclination: float
    raan: float
    arg_perigee: float
    mean_anomaly: float
    mu: float = EARTH_MU


@dataclass(frozen=True)
class State:
    """Position/velocity state in the J2000 inertial frame (m, m/s)."""

    epoch: datetime
    position: tuple[float, float, float]
    velocity: tuple[float, float, float]


def validate_elements(el: ClassicalElements) -> None:
    """Validate ranges and physical feasibility; raise OrbitError otherwise.

    An orbit is feasible only if perigee is above the WGS-84 ellipsoid surface
    (the service deliberately fails non-physical subterranean perigees instead
    of silently propagating them).
    """

    def fail(message: str) -> None:
        raise OrbitError(ErrorCode.NON_PHYSICAL_ORBIT, message)

    values = {
        "semi_major_axis_m": el.semi_major_axis_m,
        "eccentricity": el.eccentricity,
        "inclination": el.inclination,
        "raan": el.raan,
        "arg_perigee": el.arg_perigee,
        "mean_anomaly": el.mean_anomaly,
    }
    for name, value in values.items():
        if not math.isfinite(value):
            fail(f"{name} must be finite, got {value!r}")
    if el.mu <= 0.0:
        fail(f"gravitational parameter must be positive, got {el.mu}")
    if not (0.0 <= el.eccentricity < 1.0):
        raise OrbitError(
            ErrorCode.FIELD_OUT_OF_RANGE,
            f"eccentricity must satisfy 0 <= e < 1 for bound orbits, got {el.eccentricity}",
            location="eccentricity",
        )
    if el.semi_major_axis_m <= 0.0:
        raise OrbitError(
            ErrorCode.FIELD_OUT_OF_RANGE,
            f"semi-major axis must be positive, got {el.semi_major_axis_m}",
            location="semi_major_axis_m",
        )
    if not (0.0 <= el.inclination <= math.pi):
        raise OrbitError(
            ErrorCode.FIELD_OUT_OF_RANGE,
            f"inclination must be in [0, pi] rad, got {el.inclination}",
            location="inclination_rad",
        )
    perigee_radius = el.semi_major_axis_m * (1.0 - el.eccentricity)
    if perigee_radius <= EARTH_RADIUS_M:
        fail(
            "perigee radius "
            f"{perigee_radius:.1f} m is at/below the Earth surface "
            f"({EARTH_RADIUS_M:.1f} m)"
        )
    period_s = TWO_PI * math.sqrt(el.semi_major_axis_m**3 / el.mu)
    if not (60.0 <= period_s <= 200 * 86400.0):
        fail(f"orbital period {period_s:.1f} s outside supported [60 s, 200 d]")


def elements_to_state(el: ClassicalElements) -> State:
    """Keplerian elements -> Cartesian state at the same epoch (two-body).

    Frame: J2000 (EME2000) geocentric inertial.
    """
    eccentric_anomaly = mean_to_eccentric_anomaly(
        el.mean_anomaly, el.eccentricity
    )
    cos_e = math.cos(eccentric_anomaly)
    sin_e = math.sin(eccentric_anomaly)
    # Position in the perifocal (p, q, w) frame.
    p = el.semi_major_axis_m * (1.0 - el.eccentricity**2)
    x_p = el.semi_major_axis_m * (cos_e - el.eccentricity)
    y_p = (
        el.semi_major_axis_m
        * math.sqrt(1.0 - el.eccentricity**2)
        * sin_e
    )
    radius = el.semi_major_axis_m * (1.0 - el.eccentricity * cos_e)
    n = math.sqrt(el.mu / el.semi_major_axis_m**3)
    common = n * el.semi_major_axis_m**2 / radius
    vx_p = -common * sin_e
    vy_p = common * math.sqrt(1.0 - el.eccentricity**2) * cos_e

    cos_o = math.cos(el.raan)
    sin_o = math.sin(el.raan)
    cos_i = math.cos(el.inclination)
    sin_i = math.sin(el.inclination)
    cos_w = math.cos(el.arg_perigee)
    sin_w = math.sin(el.arg_perigee)

    r11 = cos_o * cos_w - sin_o * sin_w * cos_i
    r12 = -cos_o * sin_w - sin_o * cos_w * cos_i
    r21 = sin_o * cos_w + cos_o * sin_w * cos_i
    r22 = -sin_o * sin_w + cos_o * cos_w * cos_i
    r31 = sin_w * sin_i
    r32 = cos_w * sin_i

    position = (
        r11 * x_p + r12 * y_p,
        r21 * x_p + r22 * y_p,
        r31 * x_p + r32 * y_p,
    )
    velocity = (
        r11 * vx_p + r12 * vy_p,
        r21 * vx_p + r22 * vy_p,
        r31 * vx_p + r32 * vy_p,
    )
    if not all(math.isfinite(c) for c in position + velocity):
        raise OrbitError(
            ErrorCode.NON_PHYSICAL_ORBIT,
            "non-finite Cartesian state produced (degenerate orbit?)",
        )
    return State(epoch=el.epoch, position=position, velocity=velocity)


def mean_to_eccentric_anomaly(mean_anomaly: float, eccentricity: float) -> float:
    """Solve Kepler's equation M = E - e sin E (Newton iteration).

    Works for circular orbits (e=0) and high eccentricities near 1.
    Raises ITERATION_NOT_CONVERGED if the tolerance is not reached.
    """
    if not (0.0 <= eccentricity < 1.0):
        raise OrbitError(
            ErrorCode.FIELD_OUT_OF_RANGE,
            f"eccentricity out of [0, 1): {eccentricity}",
            location="eccentricity",
        )
    # Reduce into the principal branch, but remember how many whole
    # revolutions the caller supplied so the returned E stays on that branch
    # (important for M > pi when propagating).
    revolutions = round(mean_anomaly / TWO_PI)
    m = mean_anomaly - revolutions * TWO_PI
    if m <= -math.pi:
        m += TWO_PI
        revolutions -= 1
    elif m > math.pi:
        m -= TWO_PI
        revolutions += 1
    # Initial guess: a third-order starter stays convergent for high e.
    eccentric_anomaly = (
        m
        + eccentricity * math.sin(m)
        * (
            1.0
            + eccentricity * math.cos(m)
            + eccentricity**2
            * (1.0 - 1.5 * math.sin(m) ** 2)
        )
    )
    for _ in range(80):
        f = eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly) - m
        fp = 1.0 - eccentricity * math.cos(eccentric_anomaly)
        step = f / fp
        eccentric_anomaly -= step
        if abs(step) < 1e-13:
            return eccentric_anomaly + revolutions * TWO_PI
    raise OrbitError(
        ErrorCode.ITERATION_NOT_CONVERGED,
        f"Kepler equation did not converge for e={eccentricity}, M={m}",
    )


def mean_motion_to_sma(mean_motion_rad_s: float, mu: float = EARTH_MU) -> float:
    """Mean motion (rad/s) -> semi-major axis (m)."""
    if mean_motion_rad_s <= 0.0:
        raise OrbitError(
            ErrorCode.FIELD_OUT_OF_RANGE,
            f"mean motion must be positive, got {mean_motion_rad_s}",
            location="mean_motion",
        )
    return (mu / mean_motion_rad_s**2) ** (1.0 / 3.0)


def sma_to_mean_motion(semi_major_axis_m: float, mu: float = EARTH_MU) -> float:
    """Semi-major axis (m) -> mean motion (rad/s)."""
    if semi_major_axis_m <= 0.0:
        raise OrbitError(
            ErrorCode.FIELD_OUT_OF_RANGE,
            f"semi-major axis must be positive, got {semi_major_axis_m}",
            location="semi_major_axis_m",
        )
    return math.sqrt(mu / semi_major_axis_m**3)


def wrap_angle(angle: float) -> float:
    """Wrap an angle into (-pi, pi]."""
    wrapped = (angle + math.pi) % TWO_PI - math.pi
    if wrapped == -math.pi:
        return math.pi
    return wrapped


def wrap_2pi(angle: float) -> float:
    """Wrap an angle into [0, 2*pi)."""
    return angle % TWO_PI
