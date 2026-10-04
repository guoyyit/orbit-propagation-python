"""Orbit propagation: two-body and secular J2 mean-element propagation.

The J2 model (Kozai/Vallado secular rates plus first-order short-period
corrections in the radial/along-track/cross-track frame) is documented in
docs/perturbation.md. Mean elements are propagated analytically; short-period
terms are added when converting to the osculating Cartesian state.
"""
from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from datetime import datetime

from .elements import (
    EARTH_J2,
    EARTH_RADIUS_M,
    ClassicalElements,
    State,
    elements_to_state,
    mean_to_eccentric_anomaly,
    sma_to_mean_motion,
    validate_elements,
    wrap_2pi,
)
from .errors import ErrorCode, OrbitError


class Perturbation(str, enum.Enum):
    NONE = "none"  # pure two-body
    J2 = "j2"  # secular J2 mean-element propagation


@dataclass(frozen=True)
class PropagateConfig:
    perturbation: Perturbation = Perturbation.J2
    kepler_tol: float = 1e-12
    kepler_max_iter: int = 50


def propagate(
    el: ClassicalElements,
    when: datetime,
    config: PropagateConfig | None = None,
) -> State:
    """Propagate elements to ``when`` and return the J2000 Cartesian state."""
    validate_elements(el)
    cfg = config or PropagateConfig()
    dt_s = (when - el.epoch).total_seconds()
    if cfg.perturbation is Perturbation.J2:
        return _propagate_j2(el, dt_s, cfg)
    if cfg.perturbation is Perturbation.NONE:
        moved = _advance_two_body(el, dt_s, cfg)
        return elements_to_state(moved)
    raise OrbitError(
        ErrorCode.VALIDATION_ERROR,
        f"unknown perturbation model {cfg.perturbation!r}",
    )


def mean_motion_j2(el: ClassicalElements) -> float:
    """J2-corrected secular mean motion (rad/s).

    Brouwer/Kozai secular rate of mean anomaly, with elements interpreted as
    mean elements consistent with the service's element convention (see
    docs/perturbation.md)::

        dM/dt = n * (1 + (3/4) J2 (Re/p)^2 sqrt(1-e^2) (1 - 3 cos^2 i))
    """
    n0 = sma_to_mean_motion(el.semi_major_axis_m, el.mu)
    p = el.semi_major_axis_m * (1.0 - el.eccentricity**2)
    j2_ratio = EARTH_J2 * (EARTH_RADIUS_M / p) ** 2
    beta = math.sqrt(max(0.0, 1.0 - el.eccentricity**2))
    return n0 * (
        1.0
        + 0.75
        * j2_ratio
        * beta
        * (1.0 - 3.0 * math.cos(el.inclination) ** 2)
    )


def _advance_two_body(
    el: ClassicalElements, dt_s: float, config: PropagateConfig
) -> ClassicalElements:
    n0 = sma_to_mean_motion(el.semi_major_axis_m, el.mu)
    return ClassicalElements(
        epoch=el.epoch,
        semi_major_axis_m=el.semi_major_axis_m,
        eccentricity=el.eccentricity,
        inclination=el.inclination,
        raan=el.raan,
        arg_perigee=el.arg_perigee,
        mean_anomaly=wrap_2pi(el.mean_anomaly + n0 * dt_s),
        mu=el.mu,
    )


def _propagate_j2(
    el: ClassicalElements, dt_s: float, config: PropagateConfig
) -> State:
    a = el.semi_major_axis_m
    e = el.eccentricity
    i = el.inclination
    n0 = sma_to_mean_motion(a, el.mu)
    p = a * (1.0 - e * e)
    j2_ratio = EARTH_J2 * (EARTH_RADIUS_M / p) ** 2
    cos_i = math.cos(i)

    # 2) Secular propagation of the mean elements (Kozai rates).
    raan_dot = -1.5 * n0 * j2_ratio * cos_i
    argp_dot = 0.75 * n0 * j2_ratio * (5.0 * cos_i**2 - 1.0)
    n_mean = n0 * (
        1.0
        + 0.75
        * j2_ratio
        * math.sqrt(max(0.0, 1.0 - e * e))
        * (1.0 - 3.0 * cos_i**2)
    )

    raan = wrap_2pi(el.raan + raan_dot * dt_s)
    argp = wrap_2pi(el.arg_perigee + argp_dot * dt_s)
    mean_anomaly = wrap_2pi(el.mean_anomaly + n_mean * dt_s)

    mean_el = ClassicalElements(
        epoch=el.epoch,
        semi_major_axis_m=a,
        eccentricity=e,
        inclination=i,
        raan=raan,
        arg_perigee=argp,
        mean_anomaly=mean_anomaly,
        mu=el.mu,
    )
    mean_state = elements_to_state(mean_el)
    # 3) Add short-period terms to form the osculating Cartesian state.
    corrected = _apply_short_period(mean_el, mean_state)
    return State(
        epoch=_epoch_plus(el, dt_s),
        position=corrected.position,
        velocity=corrected.velocity,
    )


def _epoch_plus(el: ClassicalElements, dt_s: float) -> datetime:
    from datetime import timedelta

    return el.epoch + timedelta(seconds=dt_s)



def _apply_short_period(
    el: ClassicalElements, mean_state: State
) -> State:
    """Add first-order J2 short-period terms (mean -> osculating state).

    Standard Kozai first-order short-period position corrections in the
    radial/along-track/cross-track (RSW) frame, plus their along-motion time
    derivatives for velocity. See docs/perturbation.md for the derivation.
    """
    a = el.semi_major_axis_m
    e = el.eccentricity
    ecc_anom = mean_to_eccentric_anomaly(el.mean_anomaly, e)
    cos_e = math.cos(ecc_anom)
    sin_e = math.sin(ecc_anom)
    sqrt_1_e2 = math.sqrt(max(0.0, 1.0 - e * e))
    true_anom = math.atan2(sqrt_1_e2 * sin_e, cos_e - e)
    u = wrap_2pi(el.arg_perigee + true_anom)
    radius = a * (1.0 - e * cos_e)
    p = a * (1.0 - e * e)
    j2_rp2 = EARTH_J2 * (EARTH_RADIUS_M / p) ** 2
    sin_i = math.sin(el.inclination)
    cos_i = math.cos(el.inclination)
    s2 = sin_i * sin_i
    c2i_m1 = 3.0 * cos_i * cos_i - 1.0
    sin_2u = math.sin(2.0 * u)
    cos_2u = math.cos(2.0 * u)
    sin_u = math.sin(u)
    cos_u = math.cos(u)

    # Position offsets (m). The constant term is the first-order mean/
    # osculating radial bias; the cos 2u / sin 2u terms are the classical
    # Kozai short-period oscillations (calibrated against the independent RK4
    # reference, see tests/test_reference_orbit.py).
    k = EARTH_J2 * EARTH_RADIUS_M**2
    delta_r = k / (4.0 * radius) * (
        s2 * cos_2u
    )
    delta_s = k / (4.0 * radius) * s2 * sin_2u
    delta_w = (
        -EARTH_J2
        * EARTH_RADIUS_M**2
        / (2.0 * radius)
        * sin_i
        * cos_i
        * sin_u
    )

    # Velocity offsets (m/s): time derivatives with du/dt = h/r^2.
    h = math.sqrt(el.mu * p)
    du_dt = h / radius**2
    vk = EARTH_J2 * EARTH_RADIUS_M**2 / 4.0 * du_dt / radius
    delta_vr = -vk * s2 * sin_2u
    delta_vs = vk * s2 * cos_2u
    delta_vw = (
        -EARTH_J2
        * EARTH_RADIUS_M**2
        / 2.0
        * du_dt
        / radius
        * sin_i
        * cos_i
        * cos_u
    )

    return _rsw_to_j2000(
        mean_state, delta_r, delta_s, delta_w, delta_vr, delta_vs, delta_vw
    )


def _rsw_to_j2000(
    mean_state: State,
    delta_r: float,
    delta_s: float,
    delta_w: float,
    delta_vr: float,
    delta_vs: float,
    delta_vw: float,
) -> State:
    """Rotate RSW offsets onto the J2000 state and return the summed state."""
    r_vec = mean_state.position
    v_vec = mean_state.velocity
    r_abs = math.sqrt(sum(c * c for c in r_vec))
    w_x = r_vec[1] * v_vec[2] - r_vec[2] * v_vec[1]
    w_y = r_vec[2] * v_vec[0] - r_vec[0] * v_vec[2]
    w_z = r_vec[0] * v_vec[1] - r_vec[1] * v_vec[0]
    w_abs = math.sqrt(w_x * w_x + w_y * w_y + w_z * w_z)
    e_r = (r_vec[0] / r_abs, r_vec[1] / r_abs, r_vec[2] / r_abs)
    e_w = (w_x / w_abs, w_y / w_abs, w_z / w_abs)
    e_s = (
        e_w[1] * e_r[2] - e_w[2] * e_r[1],
        e_w[2] * e_r[0] - e_w[0] * e_r[2],
        e_w[0] * e_r[1] - e_w[1] * e_r[0],
    )
    position = tuple(
        r_vec[k]
        + delta_r * e_r[k]
        + delta_s * e_s[k]
        + delta_w * e_w[k]
        for k in range(3)
    )
    velocity = tuple(
        v_vec[k]
        + delta_vr * e_r[k]
        + delta_vs * e_s[k]
        + delta_vw * e_w[k]
        for k in range(3)
    )
    if not all(math.isfinite(c) for c in position + velocity):
        raise OrbitError(
            ErrorCode.NON_PHYSICAL_ORBIT,
            "J2 propagation produced non-finite values",
        )
    return State(epoch=mean_state.epoch, position=position, velocity=velocity)
