"""Independent high-fidelity reference propagator (RK4, J2 point-mass).

Used solely by tests as an independent comparison against the analytic
secular-J2 propagator, following the accuracy contract in docs/accuracy.md.

The equations of motion are integrated with fixed-step classical RK4 on the
J2000 inertial Cartesian state::

    r'' = -mu r / |r|^3 + a_J2(r)

where ``a_J2`` is the standard zonal acceleration. Initial conditions come
from the two-body element->state conversion, so this module shares *no*
propagation code with :mod:`app.propagator` - only constants.
"""
from __future__ import annotations

import math
from datetime import datetime

from .elements import (
    EARTH_J2,
    EARTH_MU,
    EARTH_RADIUS_M,
    ClassicalElements,
    State,
    elements_to_state,
    validate_elements,
)


def _j2_acceleration(
    position: tuple[float, float, float],
    mu: float = EARTH_MU,
    j2: float = EARTH_J2,
    radius_eq: float = EARTH_RADIUS_M,
) -> tuple[float, float, float]:
    x, y, z = position
    r = math.sqrt(x * x + y * y + z * z)
    factor = -mu / r**3
    zonal = 1.5 * j2 * (radius_eq / r) ** 2
    z_term = 5.0 * (z / r) ** 2
    ax = factor * x * (1.0 - zonal * (z_term - 1.0))
    ay = factor * y * (1.0 - zonal * (z_term - 1.0))
    az = factor * z * (1.0 - zonal * (z_term - 3.0))
    return ax, ay, az


def _derivative(state6: tuple[float, ...]) -> tuple[float, ...]:
    px, py, pz, vx, vy, vz = state6
    ax, ay, az = _j2_acceleration((px, py, pz))
    return (vx, vy, vz, ax, ay, az)


def _scale_add(
    s: tuple[float, ...], other: tuple[float, ...], scale: float
) -> tuple[float, ...]:
    return tuple(a + scale * b for a, b in zip(s, other))


def _rk4_step(state6: tuple[float, ...], h: float) -> tuple[float, ...]:
    k1 = _derivative(state6)
    k2 = _derivative(_scale_add(state6, k1, 0.5 * h))
    k3 = _derivative(_scale_add(state6, k2, 0.5 * h))
    k4 = _derivative(_scale_add(state6, k3, h))
    return tuple(
        s + h / 6.0 * (a + 2.0 * b + 2.0 * c + d)
        for s, a, b, c, d in zip(state6, k1, k2, k3, k4)
    )


def propagate_reference(
    el: ClassicalElements,
    when: datetime,
    step_s: float = 30.0,
) -> State:
    """Integrate J2 equations of motion with fixed-step RK4.

    The final fractional step (< ``step_s``) is integrated with RK4 too,
    so results are sampled exactly at ``when``. Backward propagation is
    supported via a negative step direction.
    """
    validate_elements(el)
    initial = elements_to_state(el)
    total = (when - el.epoch).total_seconds()
    if total == 0.0:
        return initial
    direction = 1.0 if total > 0 else -1.0
    remaining = abs(total)
    state6 = initial.position + initial.velocity
    while remaining > 1e-9:
        step = min(step_s, remaining) * direction
        state6 = _rk4_step(state6, step)
        remaining -= min(step_s, remaining)
    position = tuple(state6[0:3])
    velocity = tuple(state6[3:6])
    if not all(math.isfinite(c) for c in position + velocity):
        raise ArithmeticError("reference propagator produced non-finite values")
    return State(epoch=when, position=position, velocity=velocity)
