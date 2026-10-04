"""Orbit propagation: analytic two-body and numerically integrated J2.

Inputs are treated as J2 mean elements in both modes. For J2 propagation the
mean elements are converted to osculating Cartesian state (first-order J2 short
period correction) and integrated with a fixed-step RK4 over an axis-aligned
Earth (J2 gravity only). Fixed steps make the integration bit-for-bit
deterministic; the step size is configurable via ``PropagatorConfig``.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math

from .constants import WGS84
from .elements import ClassicalElements
from .elements import coe_to_eci
from .errors import ConvergenceError
from .time_utils import Epoch


class PropagationMode(str, Enum):
    TWO_BODY = "two_body"
    J2 = "j2"


@dataclass(frozen=True)
class PropagatorConfig:
    mode: PropagationMode = PropagationMode.J2
    rtol: float = 1e-10
    atol_metres: float = 1e-6
    max_steps: int = 200_000
    step_seconds: float = 15.0


@dataclass(frozen=True)
class StateVector:
    position: tuple  # (x, y, z) metres, inertial frame
    velocity: tuple  # (vx, vy, vz) m/s
    epoch: Epoch


def propagate(
    elements: ClassicalElements,
    time: Epoch,
    config: PropagatorConfig | None = None,
) -> StateVector:
    cfg = config or PropagatorConfig()
    if cfg.mode is PropagationMode.TWO_BODY:
        return propagate_two_body(elements, time)
    return propagate_j2(elements, time, cfg)


def propagate_two_body(elements: ClassicalElements, time: Epoch) -> StateVector:
    mean_motion = math.sqrt(
        WGS84.mu / elements.semi_major_axis**3
    )
    shifted = ClassicalElements(
        semi_major_axis=elements.semi_major_axis,
        eccentricity=elements.eccentricity,
        inclination=elements.inclination,
        raan=elements.raan,
        argument_of_periapsis=elements.argument_of_periapsis,
        mean_anomaly=(
            elements.mean_anomaly
            + mean_motion * (time.unix_seconds - elements.epoch.unix_seconds)
        ),
        epoch=time,
    )
    position, velocity = coe_to_eci(shifted)
    return StateVector(position=position, velocity=velocity, epoch=time)


def _j2_acceleration(x: float, y: float, z: float):
    mu = WGS84.mu
    r_eq = WGS84.radius_equator
    j2 = WGS84.j2
    r2 = x * x + y * y + z * z
    r = math.sqrt(r2)
    factor = -1.5 * j2 * mu * r_eq * r_eq / (r2 * r2 * r)
    z2_over_r2 = 5.0 * z * z / r2
    ax = -mu * x / (r2 * r) + factor * x * (1.0 - z2_over_r2)
    ay = -mu * y / (r2 * r) + factor * y * (1.0 - z2_over_r2)
    az = -mu * z / (r2 * r) + factor * z * (3.0 - z2_over_r2)
    return ax, ay, az


def _rk4_step(state, dt: float):
    x, y, z, vx, vy, vz = state

    def derivative(s):
        ax, ay, az = _j2_acceleration(s[0], s[1], s[2])
        return (s[3], s[4], s[5], ax, ay, az)

    k1 = derivative(state)
    s2 = tuple(state[i] + 0.5 * dt * k1[i] for i in range(6))
    k2 = derivative(s2)
    s3 = tuple(state[i] + 0.5 * dt * k2[i] for i in range(6))
    k3 = derivative(s3)
    s4 = tuple(state[i] + dt * k3[i] for i in range(6))
    k4 = derivative(s4)
    return tuple(
        state[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        for i in range(6)
    )


def _mean_to_osc(position, velocity):
    """First-order J2 short-period correction (Kozai), mean -> osculating."""
    rx, ry, rz = position
    vx, vy, vz = velocity
    mu = WGS84.mu
    r_eq = WGS84.radius_equator
    j2 = WGS84.j2
    radius = math.sqrt(rx * rx + ry * ry + rz * rz)
    rxy = math.sqrt(rx * rx + ry * ry)
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    h = math.sqrt(hx * hx + hy * hy + hz * hz)
    p = h * h / mu
    if p <= 0.0 or rxy <= 0.0 or radius <= 0.0:
        return position, velocity
    node_x, node_y = -hy, hx
    node_mag = math.sqrt(node_x * node_x + node_y * node_y)
    if node_mag <= 0.0:
        return position, velocity

    n_hat = (node_x / node_mag, node_y / node_mag, 0.0)
    h_hat = (hx / h, hy / h, hz / h)
    m_hat = (
        h_hat[1] * n_hat[2] - h_hat[2] * n_hat[1],
        h_hat[2] * n_hat[0] - h_hat[0] * n_hat[2],
        h_hat[0] * n_hat[1] - h_hat[1] * n_hat[0],
    )
    radial_hat = (rx / radius, ry / radius, rz / radius)
    along_hat = (
        h_hat[1] * radial_hat[2] - h_hat[2] * radial_hat[1],
        h_hat[2] * radial_hat[0] - h_hat[0] * radial_hat[2],
        h_hat[0] * radial_hat[1] - h_hat[1] * radial_hat[0],
    )

    # Argument of latitude u and its trig functions, measured from the node.
    sin_u = (n_hat[0] * ry - n_hat[1] * rx) / rxy
    cos_u = (n_hat[0] * rx + n_hat[1] * ry) / rxy
    u2 = 2.0 * math.atan2(sin_u, cos_u)
    sin_2u = 2.0 * sin_u * cos_u
    cos_2u = cos_u * cos_u - sin_u * sin_u

    rv = rx * vx + ry * vy + rz * vz
    speed2 = vx * vx + vy * vy + vz * vz
    ex = (speed2 - mu / radius) * rx / mu - rv * vx / mu
    ey = (speed2 - mu / radius) * ry / mu - rv * vy / mu
    ez = (speed2 - mu / radius) * rz / mu - rv * vz / mu
    e_sin_u = ex * n_hat[0] + ey * n_hat[1]
    e_cos_u = (
        ex * m_hat[0] + ey * m_hat[1] + ez * m_hat[2]
    )

    sin_i = node_mag / h
    cos_i = hz / h
    si2 = sin_i * sin_i
    c = 0.75 * j2 * (r_eq / p) ** 2
    c2 = c * 2.0
    vp = math.sqrt(mu / p)

    # Radial / transverse / cross-track frame (RSW).
    r_hat = radial_hat
    s_hat = along_hat
    w_hat = (
        r_hat[1] * s_hat[2] - r_hat[2] * s_hat[1],
        r_hat[2] * s_hat[0] - r_hat[0] * s_hat[2],
        r_hat[0] * s_hat[1] - r_hat[1] * s_hat[0],
    )

    dr_r = c2 * p * si2 * cos_2u * 0.5
    dr_s = -c2 * p * si2 * sin_2u * 0.5
    dr_w = c2 * p * sin_i * cos_i * sin_u
    dv_r = -c2 * vp * si2 * sin_2u * 0.5
    dv_s = c * vp * (
        (3.0 * si2 - 2.0) * (e_sin_u * sin_u + e_cos_u * cos_u)
        + si2 * cos_2u
    )
    dv_w = 2.0 * c * vp * sin_i * cos_i * (e_cos_u - cos_u)

    osc_position = (
        rx + dr_r * r_hat[0] + dr_s * s_hat[0] + dr_w * w_hat[0],
        ry + dr_r * r_hat[1] + dr_s * s_hat[1] + dr_w * w_hat[1],
        rz + dr_r * r_hat[2] + dr_s * s_hat[2] + dr_w * w_hat[2],
    )
    osc_velocity = (
        vx + dv_r * r_hat[0] + dv_s * s_hat[0] + dv_w * w_hat[0],
        vy + dv_r * r_hat[1] + dv_s * s_hat[1] + dv_w * w_hat[1],
        vz + dv_r * r_hat[2] + dv_s * s_hat[2] + dv_w * w_hat[2],
    )
    return osc_position, osc_velocity


def propagate_j2(
    elements: ClassicalElements,
    time: Epoch,
    config: PropagatorConfig | None = None,
) -> StateVector:
    cfg = config or PropagatorConfig()
    dt_total = time.unix_seconds - elements.epoch.unix_seconds
    if dt_total == 0.0:
        position, velocity = coe_to_eci(elements)
        return StateVector(position=position, velocity=velocity, epoch=time)

    mean_position, mean_velocity = coe_to_eci(elements)
    position, velocity = _mean_to_osc(mean_position, mean_velocity)
    state = (
        position[0], position[1], position[2],
        velocity[0], velocity[1], velocity[2],
    )

    step = abs(cfg.step_seconds)
    if not math.isfinite(step) or step <= 0.0:
        raise ConvergenceError("step_seconds must be positive finite")
    direction = 1.0 if dt_total > 0.0 else -1.0
    remaining = abs(dt_total)
    steps = int(math.ceil(remaining / step - 1e-12))
    if steps > cfg.max_steps:
        raise ConvergenceError(
            f"propagation needs {steps} steps, exceeding max_steps={cfg.max_steps}"
        )
    base_dt = direction * remaining / steps
    for index in range(steps):
        state = _rk4_step(state, base_dt)
        if not all(math.isfinite(component) for component in state):
            raise ConvergenceError(
                f"non-finite state at integration step {index}", location="time"
            )
    return StateVector(
        position=(state[0], state[1], state[2]),
        velocity=(state[3], state[4], state[5]),
        epoch=time,
    )


class J2Stepper:
    """Stateful, monotonically advancing J2 integrator.

    Used by window/event scanners so each requested time continues from the
    previous one instead of restarting at the epoch. The sequence is purely a
    function of the requested times and the fixed step size, hence identical
    regardless of concurrency or call scheduling.
    """

    def __init__(
        self,
        elements: ClassicalElements,
        config: PropagatorConfig | None = None,
    ) -> None:
        cfg = config or PropagatorConfig()
        self._step = abs(cfg.step_seconds)
        if not math.isfinite(self._step) or self._step <= 0.0:
            raise ConvergenceError("step_seconds must be positive finite")
        self._max_steps = cfg.max_steps
        mean_position, mean_velocity = coe_to_eci(elements)
        position, velocity = _mean_to_osc(mean_position, mean_velocity)
        self._state = (
            position[0], position[1], position[2],
            velocity[0], velocity[1], velocity[2],
        )
        self._time = elements.epoch.unix_seconds
        self._step = abs(cfg.step_seconds)
        if not math.isfinite(self._step) or self._step <= 0.0:
            raise ConvergenceError("step_seconds must be positive finite")
        self._max_steps = cfg.max_steps

    @classmethod
    def from_state(
        cls,
        state: StateVector,
        config: PropagatorConfig | None = None,
    ) -> "J2Stepper":
        """Stepper initialized from an already-propagated osculating state."""
        self = cls.__new__(cls)
        cfg = config or PropagatorConfig()
        self._state = (
            state.position[0], state.position[1], state.position[2],
            state.velocity[0], state.velocity[1], state.velocity[2],
        )
        self._time = state.epoch.unix_seconds
        self._step = abs(cfg.step_seconds)
        if not math.isfinite(self._step) or self._step <= 0.0:
            raise ConvergenceError("step_seconds must be positive finite")
        self._max_steps = cfg.max_steps
        return self

    def advance_to(self, unix_seconds: float) -> StateVector:
        delta = unix_seconds - self._time
        if delta < 0.0:
            raise ConvergenceError(
                "J2Stepper only advances forward in time", location="time"
            )
        if delta > 0.0:
            steps = max(1, int(math.ceil(delta / self._step - 1e-12)))
            if steps > self._max_steps:
                raise ConvergenceError(
                    f"integration needs {steps} steps, "
                    f"exceeding max_steps={self._max_steps}"
                )
            dt = delta / steps
            for index in range(steps):
                self._state = _rk4_step(self._state, dt)
                if not all(math.isfinite(v) for v in self._state):
                    raise ConvergenceError(
                        f"non-finite state at integration step {index}",
                        location="time",
                    )
            self._time = unix_seconds
        return StateVector(
            position=(self._state[0], self._state[1], self._state[2]),
            velocity=(self._state[3], self._state[4], self._state[5]),
            epoch=Epoch(unix_seconds),
        )
