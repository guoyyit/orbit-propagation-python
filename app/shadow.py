"""Earth-shadow (eclipse) models and illumination events.

Two geometric umbra models are supported (see docs/shadow.md):

* ``cylindrical``: the Sun is treated as infinitely distant with zero angular
  radius, so the umbra is a cylinder of radius ``R_E`` anti-parallel to the
  Sun direction.  Simple, but ignores the finite angular size of the Sun.
* ``conical``: the finite Sun (angular radius ~0.2666 deg) casts a conical
  umbra whose radius shrinks with distance behind the Earth.  More accurate;
  near LEO the two models differ by a few seconds per eclipse boundary.

The penumbra is geometrically detected but this implementation reports the
umbra boundary and treats penumbra as lit (documented limitation).
"""
from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from datetime import datetime

from .propagator import PropagateConfig, propagate
from .timeutil import datetime_to_jd

SUN_RADIUS_M = 6.957e8
SUN_DISTANCE_M = 1.495978707e11  # 1 AU, low-precision circular solar model
EARTH_RADIUS_M = 6_378_137.0


class ShadowModel(str, enum.Enum):
    CYLINDRICAL = "cylindrical"
    CONICAL = "conical"


@dataclass(frozen=True)
class IlluminationEvent:
    """Shadow boundary (``entering_shadow`` / ``leaving_shadow``)."""

    time: datetime
    kind: str


@dataclass(frozen=True)
class IlluminationResult:
    events: tuple[IlluminationEvent, ...]
    illuminated_seconds: float
    shadowed_seconds: float


def sun_position_j2000(when: datetime) -> tuple[float, float, float]:
    """Apparent geocentric Sun position in J2000 (m), low-precision model.

    Mean ecliptic longitude of the Sun (Vallado Alg. 29, ~0.01 deg accuracy),
    obliquity 23.439291 deg, circular 1-AU orbit.
    """
    jd = datetime_to_jd(when)
    t_ut1 = (jd - 2451545.0) / 36525.0
    # Mean longitude of the Sun (deg) and mean anomaly.
    mean_long = (280.460 + 36000.771 * t_ut1) % 360.0
    mean_anom = math.radians((357.5291092 + 35999.05034 * t_ut1) % 360.0)
    ecliptic_long = math.radians(
        mean_long + 1.914666471 * math.sin(mean_anom)
        + 0.019994643 * math.sin(2.0 * mean_anom)
    )
    obliquity = math.radians(23.439291)
    cos_obl = math.cos(obliquity)
    sin_obl = math.sin(obliquity)
    x = SUN_DISTANCE_M * math.cos(ecliptic_long)
    y = SUN_DISTANCE_M * cos_obl * math.sin(ecliptic_long)
    z = SUN_DISTANCE_M * sin_obl * math.sin(ecliptic_long)
    return (x, y, z)


def illumination_fraction(
    satellite_position: tuple[float, float, float],
    sun_position: tuple[float, float, float],
    model: ShadowModel,
) -> float:
    """Return 1.0 when fully sunlit, 0.0 when inside the umbra.

    Penumbra is treated as lit (returns 1.0); this function is the indicator
    used for event detection and illuminated-duration integration.
    """
    # Vector from Earth to satellite, and the unit Sun direction (Earth->Sun).
    sat_x, sat_y, sat_z = satellite_position
    sun_abs = math.sqrt(sum(c * c for c in sun_position))
    sun_hat = tuple(c / sun_abs for c in sun_position)
    sat_abs = math.sqrt(sat_x * sat_x + sat_y * sat_y + sat_z * sat_z)

    # Project the satellite along the Sun direction. Positive => sunlit side.
    along = (
        sat_x * sun_hat[0] + sat_y * sun_hat[1] + sat_z * sun_hat[2]
    )
    if along >= 0.0:
        return 1.0  # on the Sun-facing hemisphere, never in shadow

    # Perpendicular distance from the satellite to the shadow axis.
    perp_sq = sat_abs * sat_abs - along * along
    perp = math.sqrt(max(0.0, perp_sq))
    distance_behind = -along

    if model is ShadowModel.CYLINDRICAL:
        shadow_radius = EARTH_RADIUS_M
    else:
        # Conical umbra half-angle: tan(theta) = (R_sun - R_earth) / AU.
        cone_slope = (SUN_RADIUS_M - EARTH_RADIUS_M) / SUN_DISTANCE_M
        shadow_radius = EARTH_RADIUS_M - cone_slope * distance_behind
        if shadow_radius <= 0.0:
            return 1.0  # beyond the end of the umbral cone

    return 0.0 if perp <= shadow_radius else 1.0


def find_illumination(
    el,
    start: datetime,
    end: datetime,
    model: ShadowModel = ShadowModel.CONICAL,
    boundary_tol_s: float = 1.0,
    sample_step_s: float = 60.0,
    config: PropagateConfig | None = None,
) -> IlluminationResult:
    """Compute umbra entry/exit events and illuminated duration in [start,end]."""
    if end <= start:
        return IlluminationResult((), 0.0, 0.0)
    total = (end - start).total_seconds()
    n_steps = max(1, math.ceil(total / sample_step_s))

    def lit_at(t_s: float) -> float:
        when = start + _seconds(t_s)
        state = propagate(el, when, config)
        sun = sun_position_j2000(when)
        return illumination_fraction(state.position, sun, model)

    times = [min(total, k * sample_step_s) for k in range(n_steps + 1)]
    values = [lit_at(t) for t in times]

    def refine(lo: float, hi: float) -> float:
        flo = lit_at(lo)
        for _ in range(80):
            if abs(hi - lo) <= boundary_tol_s:
                break
            mid = 0.5 * (lo + hi)
            fm = lit_at(mid)
            if (flo <= 0.0) == (fm <= 0.0):
                lo = mid
                flo = fm
            else:
                hi = mid
        return 0.5 * (lo + hi)

    events: list[IlluminationEvent] = []
    shadowed = 0.0
    in_shadow = values[0] <= 0.0
    segment_start = 0.0
    for k in range(n_steps):
        t0, t1 = times[k], times[k + 1]
        v0, v1 = values[k], values[k + 1]
        if v0 >= 1.0 and v1 <= 0.0:  # entering shadow
            boundary = refine(t0, t1)
            events.append(
                IlluminationEvent(start + _seconds(boundary), "entering_shadow")
            )
            in_shadow = True
            segment_start = boundary
        elif v0 <= 0.0 and v1 >= 1.0:  # leaving shadow
            boundary = refine(t0, t1)
            events.append(
                IlluminationEvent(start + _seconds(boundary), "leaving_shadow")
            )
            shadowed += boundary - segment_start
            in_shadow = False
    if in_shadow:
        shadowed += total - segment_start

    illuminated = total - shadowed
    return IlluminationResult(
        events=tuple(events),
        illuminated_seconds=illuminated,
        shadowed_seconds=shadowed,
    )


def _seconds(value: float):
    from datetime import timedelta

    return timedelta(seconds=value)
