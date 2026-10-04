"""Ground stations, coordinate frames and visibility windows.

Geodetic coordinates use WGS84. The station Earth-fixed position is rotated
into the service inertial frame through Greenwich Mean Sidereal Time.
Visibility windows are found by a sign scan of ``elevation - mask`` followed
by bisection, so tightening ``boundary_tolerance_seconds`` only moves each
boundary inward monotonically: windows can never appear, vanish, or reorder.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .constants import WGS84
from .elements import ClassicalElements
from .errors import ConvergenceError, FieldRangeError
from .propagator import PropagatorConfig, PropagationMode, StateVector
from .propagator import J2Stepper, propagate, propagate_two_body
from .time_utils import Epoch
from .time_utils import gmst_radians


@dataclass(frozen=True)
class GroundStation:
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_metres: float = 0.0


@dataclass(frozen=True)
class VisibilityWindow:
    start: Epoch
    end: Epoch
    maximum_elevation_deg: float
    maximum_at: Epoch


def geodetic_to_ecef(station: GroundStation):
    lat = math.radians(station.latitude_deg)
    lon = math.radians(station.longitude_deg)
    f = WGS84.flattening
    e2 = f * (2.0 - f)
    sin_lat = math.sin(lat)
    prime_radius = WGS84.radius_equator / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    r = (prime_radius + station.altitude_metres) * math.cos(lat)
    x = r * math.cos(lon)
    y = r * math.sin(lon)
    z = (prime_radius * (1.0 - e2) + station.altitude_metres) * sin_lat
    return x, y, z


def station_eci(station: GroundStation, unix_seconds: float):
    """Earth-fixed station position rotated into the inertial frame."""
    x, y, z = geodetic_to_ecef(station)
    theta = gmst_radians(unix_seconds)
    cos_t, sin_t = math.cos(theta), math.sin(theta)
    return (
        x * cos_t - y * sin_t,
        x * sin_t + y * cos_t,
        z,
    )


def elevation_deg(station: GroundStation, state: StateVector) -> float:
    station_position = station_eci(station, state.epoch.unix_seconds)
    rho = tuple(s - g for s, g in zip(state.position, station_position))
    slant = math.sqrt(sum(component * component for component in rho))
    station_radius = math.sqrt(sum(g * g for g in station_position))
    if slant <= 0.0 or station_radius <= 0.0:
        return -90.0
    # Local vertical points along the station radius vector.
    cos_zenith = sum(
        g * d for g, d in zip(station_position, rho)
    ) / (station_radius * slant)
    return math.degrees(math.asin(max(-1.0, min(1.0, cos_zenith))))


def validate_station(station: GroundStation) -> None:
    if not (
        math.isfinite(station.latitude_deg)
        and math.isfinite(station.longitude_deg)
        and math.isfinite(station.altitude_metres)
    ):
        raise FieldRangeError("ground station coordinates must be finite")
    if not -90.0 <= station.latitude_deg <= 90.0:
        raise FieldRangeError(
            f"latitude {station.latitude_deg} outside [-90, 90]",
            location="latitude_deg",
        )
    if not -180.0 <= station.longitude_deg <= 180.0:
        raise FieldRangeError(
            f"longitude {station.longitude_deg} outside [-180, 180]",
            location="longitude_deg",
        )
    if station.altitude_metres < -WGS84.radius_equator * WGS84.flattening or station.altitude_metres > 1e5:
        raise FieldRangeError(
            f"altitude {station.altitude_metres} outside physical range",
            location="altitude_metres",
        )


def _elevation_at(
    elements: ClassicalElements,
    station: GroundStation,
    unix_seconds: float,
    config: PropagatorConfig | None,
) -> float:
    state = propagate(elements, Epoch(unix_seconds), config)
    return elevation_deg(station, state)


def compute_windows(
    elements: ClassicalElements,
    station: GroundStation,
    start: Epoch,
    end: Epoch,
    minimum_elevation_deg: float,
    boundary_tolerance_seconds: float,
    config: PropagatorConfig | None = None,
) -> list[VisibilityWindow]:
    validate_station(station)
    if not -90.0 <= minimum_elevation_deg <= 90.0:
        raise FieldRangeError(
            f"minimum elevation {minimum_elevation_deg} outside [-90, 90]",
            location="minimum_elevation_deg",
        )
    if end.unix_seconds <= start.unix_seconds:
        raise FieldRangeError(
            "end time must be strictly after start time", location="end"
        )
    if boundary_tolerance_seconds <= 0.0:
        raise FieldRangeError(
            "boundary_tolerance_seconds must be positive",
            location="boundary_tolerance_seconds",
        )

    t0 = start.unix_seconds
    t1 = end.unix_seconds
    mean_motion = math.sqrt(WGS84.mu / elements.semi_major_axis**3)
    period = 2.0 * math.pi / mean_motion
    # A scan step comfortably shorter than a pass guarantees no pass is missed.
    scan_step = min(60.0, period / 16.0)
    sample_count = max(1, int(math.ceil((t1 - t0) / scan_step)))
    sample_step = (t1 - t0) / sample_count

    times: list[float] = [t0 + index * sample_step for index in range(sample_count + 1)]
    times[-1] = t1
    cfg = config or PropagatorConfig()
    if cfg.mode is PropagationMode.J2:
        stepper = J2Stepper(elements, cfg)
        samples = [stepper.advance_to(sample_time) for sample_time in times]
    else:
        samples = [propagate_two_body(elements, Epoch(sample_time)) for sample_time in times]
    elevations = [
        elevation_deg(station, state) - minimum_elevation_deg for state in samples
    ]

    def residual(unix_seconds: float) -> float:
        ratio = (unix_seconds - t0) / sample_step
        index = int(math.floor(ratio))
        if index <= 0:
            return elevations[0]
        if index >= sample_count:
            return elevations[-1]
        fraction = ratio - index
        return elevations[index] * (1.0 - fraction) + elevations[index + 1] * fraction

    sign_changes: list[tuple[float, float, float, float]] = []
    for index in range(sample_count):
        if elevations[index] * elevations[index + 1] <= 0.0 and not (
            elevations[index] == 0.0 and elevations[index + 1] == 0.0
        ):
            sign_changes.append(
                (
                    times[index],
                    times[index + 1],
                    elevations[index],
                    elevations[index + 1],
                )
            )

    def bisect(left: float, right: float, left_value: float) -> float:
        for _ in range(100):
            if right - left <= boundary_tolerance_seconds:
                break
            middle = 0.5 * (left + right)
            middle_value = residual(middle)
            if left_value * middle_value <= 0.0:
                right = middle
            else:
                left = middle
                left_value = middle_value
        else:
            raise ConvergenceError("window boundary bisection did not converge")
        return 0.5 * (left + right)

    starts: list[float] = []
    ends: list[float] = []
    if elevations[0] >= 0.0:
        starts.append(t0)
    for left, right, lv, rv in sign_changes:
        root = bisect(left, right, lv)
        if lv < 0.0 <= rv:
            starts.append(root)
        elif lv >= 0.0 > rv:
            ends.append(root)
    if elevations[-1] >= 0.0:
        ends.append(t1)

    windows: list[VisibilityWindow] = []
    for window_start, window_end in zip(starts, ends):
        # Maximum elevation by golden-section search on [start, end].
        golden = (math.sqrt(5.0) - 1.0) / 2.0
        low, high = window_start, window_end
        x1 = high - golden * (high - low)
        x2 = low + golden * (high - low)
        f1 = residual(x1)
        f2 = residual(x2)
        for _ in range(60):
            if f1 < f2:
                low = x1
                x1, f1 = x2, f2
                x2 = low + golden * (high - low)
                f2 = residual(x2)
            else:
                high = x2
                x2, f2 = x1, f1
                x1 = high - golden * (high - low)
                f1 = residual(x1)
        maximum_at = 0.5 * (low + high)
        maximum_elevation = _elevation_at(
            elements, station, maximum_at, config
        )
        windows.append(
            VisibilityWindow(
                start=Epoch(window_start),
                end=Epoch(window_end),
                maximum_elevation_deg=maximum_elevation,
                maximum_at=Epoch(maximum_at),
            )
        )
    windows.sort(key=lambda window: window.start.unix_seconds)
    return windows
