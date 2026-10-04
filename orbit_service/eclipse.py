"""Solar illumination: cylindrical and conical Earth shadow models.

The Sun position uses a low-precision geometric ephemeris (mean longitude of
the Sun on the ecliptic, obliquity fixed at J2000) which is accurate to about
a degree and far below the event-time sensitivity of LEO shadow crossings.

* ``cylindrical``: Earth casts an infinite cylinder of radius R_E.
* ``conical``: Earth casts the true umbral cone using the Sun's finite angular
  radius (~0.2666 deg). The cone is shorter and narrower, so the conical model
  predicts slightly later entry / earlier exit and shorter eclipses than the
  cylindrical one.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from .constants import AU_METRES, WGS84
from .elements import ClassicalElements
from .errors import ConvergenceError, FieldRangeError
from .propagator import PropagatorConfig
from .propagator import J2Stepper, propagate_two_body, PropagationMode
from .time_utils import Epoch
from .time_utils import julian_day


class ShadowModel(str, Enum):
    CYLINDRICAL = "cylindrical"
    CONICAL = "conical"


@dataclass(frozen=True)
class EclipseEvent:
    time: Epoch
    kind: str  # "entry" (lit -> dark) or "exit" (dark -> lit)


@dataclass(frozen=True)
class IlluminationResult:
    events: list
    illuminated_seconds: float
    total_seconds: float


def sun_position_eci(unix_seconds: float):
    """Low-precision apparent Sun position in the inertial frame, metres."""
    days = julian_day(unix_seconds) - 2451545.0
    mean_longitude = math.radians((280.460 + 0.9856474 * days) % 360.0)
    mean_anomaly = math.radians((357.528 + 0.9856003 * days) % 360.0)
    longitude = mean_longitude + math.radians(1.915) * math.sin(mean_anomaly)
    longitude += math.radians(0.020) * math.sin(2.0 * mean_anomaly)
    obliquity = math.radians(23.4392911)
    distance = (
        AU_METRES
        * (1.0 - 0.01671 ** 2)
        / (1.0 + 0.01671 * math.cos(mean_anomaly))
    )
    x = distance * math.cos(longitude)
    y = distance * math.cos(obliquity) * math.sin(longitude)
    z = distance * math.sin(obliquity) * math.sin(longitude)
    return x, y, z


def is_illuminated(position, unix_seconds: float, model: ShadowModel) -> bool:
    sun = sun_position_eci(unix_seconds)
    sun_distance = math.sqrt(sum(component * component for component in sun))
    sun_hat = tuple(component / sun_distance for component in sun)
    rx, ry, rz = position
    projection = rx * sun_hat[0] + ry * sun_hat[1] + rz * sun_hat[2]
    perp2 = (
        rx * rx + ry * ry + rz * rz - projection * projection
    )
    # The shadow lies on the anti-solar side of Earth.
    if projection >= 0.0:
        return True

    if model is ShadowModel.CYLINDRICAL:
        return perp2 > WGS84.radius_equator**2

    # Conical umbra. Apparent half angles as seen from the body and the Earth.
    sun_radius = 6.957e8
    earth_radius = WGS84.radius_equator
    sin_theta_sun = (sun_radius - earth_radius) / sun_distance
    # Distance from Earth centre to the umbra tip along the anti-solar axis.
    cone_length = earth_radius / sin_theta_sun
    axial_depth = -projection
    if axial_depth >= cone_length:
        return True
    shadow_radius = earth_radius - axial_depth * sin_theta_sun
    return math.sqrt(max(perp2, 0.0)) > shadow_radius


def _sample_states(elements, start: Epoch, end: Epoch, sample_step: float, config):
    count = max(1, int(math.ceil((end.unix_seconds - start.unix_seconds) / sample_step)))
    step_span = (end.unix_seconds - start.unix_seconds) / count
    times = [start.unix_seconds + index * step_span for index in range(count + 1)]
    times[-1] = end.unix_seconds
    cfg = config or PropagatorConfig()
    if cfg.mode is PropagationMode.J2:
        stepper = J2Stepper(elements, cfg)
        states = [stepper.advance_to(instant) for instant in times]
    else:
        states = [propagate_two_body(elements, Epoch(instant)) for instant in times]
    return times, states


def compute_illumination(
    elements: ClassicalElements,
    start: Epoch,
    end: Epoch,
    model: ShadowModel,
    boundary_tolerance_seconds: float,
    config: PropagatorConfig | None = None,
) -> IlluminationResult:
    if end.unix_seconds <= start.unix_seconds:
        raise FieldRangeError(
            "end time must be strictly after start time", location="end"
        )
    if boundary_tolerance_seconds <= 0.0:
        raise FieldRangeError(
            "boundary_tolerance_seconds must be positive",
            location="boundary_tolerance_seconds",
        )

    mean_motion = math.sqrt(WGS84.mu / elements.semi_major_axis**3)
    period = 2.0 * math.pi / mean_motion
    sample_step = min(60.0, period / 16.0)
    times, states = _sample_states(
        elements, start, end, sample_step, config
    )
    lit_flags = [
        is_illuminated(state.position, instant, model)
        for instant, state in zip(times, states)
    ]

    cfg = config or PropagatorConfig()

    events: list[EclipseEvent] = []
    illuminated_seconds = 0.0

    def refine_root(
        left_time: float, right_time: float, left_lit: bool, right_lit: bool
    ) -> float:
        """Classic bisection inside one sub-step (at most the scan step)."""
        low, high = left_time, right_time
        j2 = cfg.mode is PropagationMode.J2
        anchor_state = states[index]

        def lit_value(instant: float) -> bool:
            if j2:
                probe = J2Stepper.from_state(anchor_state, cfg)
                state = probe.advance_to(instant)
            else:
                state = propagate_two_body(elements, Epoch(instant))
            return is_illuminated(state.position, instant, model)

        for _ in range(100):
            if high - low <= boundary_tolerance_seconds:
                return low
            mid = 0.5 * (low + high)
            mid_lit = lit_value(mid)
            if mid_lit == right_lit:
                high = mid
            else:
                low = mid
        raise ConvergenceError("eclipse boundary bisection did not converge")

    for index in range(len(times) - 1):
        left_lit, right_lit = lit_flags[index], lit_flags[index + 1]
        left_time, right_time = times[index], times[index + 1]
        if left_lit != right_lit:
            root = refine_root(left_time, right_time, left_lit, right_lit)
            events.append(
                EclipseEvent(
                    time=Epoch(root),
                    kind="entry" if right_lit is False else "exit",
                )
            )
            if left_lit:
                illuminated_seconds += root - left_time
            else:
                illuminated_seconds += right_time - root
        elif left_lit:
            illuminated_seconds += right_time - left_time

    return IlluminationResult(
        events=events,
        illuminated_seconds=illuminated_seconds,
        total_seconds=end.unix_seconds - start.unix_seconds,
    )
