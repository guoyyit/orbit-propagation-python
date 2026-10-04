"""Ground-station visibility window computation.

Visibility windows are intervals within which the satellite elevation as seen
from a ground station is at/above ``min_elevation_deg``.  The time range is
scanned at ``sample_step_s``; every detected sign crossing of the elevation
mask is refined to ``boundary_tol_s`` with bisection.  Refining the tolerance
only shifts each boundary inward/outward monotonically, so windows never jump,
disappear, or reorder (see docs/accuracy.md and the convergence test).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime

from .elements import State
from .frames import ecef_to_enu, geodetic_to_ecef, j2000_to_ecef
from .propagator import PropagateConfig, propagate


@dataclass(frozen=True)
class GroundStation:
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float = 0.0


@dataclass(frozen=True)
class VisibilityWindow:
    start: datetime
    end: datetime
    max_elevation_deg: float


def station_ecef(station: GroundStation) -> tuple[float, float, float]:
    """ECEF position of the station (m)."""
    return geodetic_to_ecef(
        math.radians(station.latitude_deg),
        math.radians(station.longitude_deg),
        station.altitude_m,
    )


def elevation_at(station: GroundStation, state: State) -> float:
    """Elevation angle (rad) of the satellite as seen from the station."""
    sat_ecef = j2000_to_ecef(state.position, state.epoch)
    sta_ecef = station_ecef(station)
    delta = (
        sat_ecef[0] - sta_ecef[0],
        sat_ecef[1] - sta_ecef[1],
        sat_ecef[2] - sta_ecef[2],
    )
    enu = ecef_to_enu(
        delta,
        math.radians(station.latitude_deg),
        math.radians(station.longitude_deg),
    )
    horizontal = math.hypot(enu[0], enu[1])
    return math.atan2(enu[2], horizontal)


def find_windows(
    el,
    station: GroundStation,
    start: datetime,
    end: datetime,
    min_elevation_deg: float,
    boundary_tol_s: float = 1.0,
    config: PropagateConfig | None = None,
    sample_step_s: float = 60.0,
) -> list[VisibilityWindow]:
    """Find all visibility windows within [start, end], time-ordered.

    The scan step bounds the minimum detectable pass duration; boundary times
    are bisected to ``boundary_tol_s``.  Windows that only touch the interval
    endpoints are clipped to ``start``/``end``.
    """
    if end <= start:
        return []
    if not (-90.0 <= station.latitude_deg <= 90.0):
        from .errors import ErrorCode, OrbitError

        raise OrbitError(
            ErrorCode.FIELD_OUT_OF_RANGE,
            f"station latitude must be in [-90, 90], got {station.latitude_deg}",
            location="station.latitude_deg",
        )
    mask = math.radians(min_elevation_deg)
    total = (end - start).total_seconds()
    n_steps = max(1, math.ceil(total / sample_step_s))

    def elev_at(t_s: float) -> float:
        when = start + _seconds(t_s)
        state = propagate(el, when, config)
        return elevation_at(station, state)

    times = [
        min(total, k * sample_step_s) for k in range(n_steps + 1)
    ]
    values = [elev_at(t) - mask for t in times]

    def refine(lo: float, hi: float) -> float:
        """Bisect f(elev - mask) == 0 on [lo, hi] to boundary_tol_s."""
        flo = elev_at(lo) - mask
        for _ in range(80):
            if abs(hi - lo) <= boundary_tol_s:
                break
            mid = 0.5 * (lo + hi)
            fm = elev_at(mid) - mask
            if (flo <= 0.0) == (fm <= 0.0):
                lo = mid
                flo = fm
            else:
                hi = mid
        return 0.5 * (lo + hi)

    windows: list[VisibilityWindow] = []
    window_start: float | None = 0.0 if values[0] >= 0.0 else None
    for k in range(n_steps):
        t0, t1 = times[k], times[k + 1]
        v0, v1 = values[k], values[k + 1]
        if v0 < 0.0 <= v1:  # rising crossing
            window_start = refine(t0, t1)
        elif v0 >= 0.0 > v1:  # setting crossing
            window_end = refine(t0, t1)
            start_t = window_start if window_start is not None else t0
            windows.append(
                _build_window(el, station, start, start_t, window_end, config)
            )
            window_start = None
    if window_start is not None and values[-1] >= 0.0:
        windows.append(
            _build_window(el, station, start, window_start, total, config)
        )
    return windows


def _seconds(value: float):
    from datetime import timedelta

    return timedelta(seconds=value)


def _build_window(el, station, start, start_t, end_t, config) -> VisibilityWindow:
    win_start = start + _seconds(start_t)
    win_end = start + _seconds(end_t)
    # Maximum elevation: golden-section-free dense scan within the window at
    # the sampling resolution (elevation is unimodal within one pass).
    duration = end_t - start_t
    n = max(2, math.ceil(duration / 30.0))
    max_el = -math.inf
    for k in range(n + 1):
        t = start_t + duration * k / n
        state = propagate(el, start + _seconds(t), config)
        max_el = max(max_el, elevation_at(station, state))
    return VisibilityWindow(
        start=win_start,
        end=win_end,
        max_elevation_deg=math.degrees(max_el),
    )
