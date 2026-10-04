"""Tests for two-body conservation and J2 vs the independent RK4 reference."""
from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone

import pytest

from app.elements import EARTH_MU, ClassicalElements
from app.errors import ErrorCode, OrbitError
from app.propagator import Perturbation, PropagateConfig, propagate
from app.reference import propagate_reference

EPOCH = datetime(2024, 1, 1, tzinfo=timezone.utc)
LEO = ClassicalElements(
    EPOCH,
    7_000_000.0,
    0.001,
    math.radians(51.6),
    math.radians(120.0),
    math.radians(40.0),
    math.radians(10.0),
)

NONE = PropagateConfig(perturbation=Perturbation.NONE)
J2 = PropagateConfig(perturbation=Perturbation.J2)


def _energy(state) -> float:
    r = math.sqrt(sum(c * c for c in state.position))
    v2 = sum(c * c for c in state.velocity)
    return v2 / 2 - EARTH_MU / r


def _angular_momentum(state) -> float:
    p, v = state.position, state.velocity
    return math.sqrt(
        (p[1] * v[2] - p[2] * v[1]) ** 2
        + (p[2] * v[0] - p[0] * v[2]) ** 2
        + (p[0] * v[1] - p[1] * v[0]) ** 2
    )


def test_two_body_energy_angular_momentum_drift_bounded() -> None:
    initial = propagate(LEO, EPOCH, NONE)
    e0 = _energy(initial)
    h0 = _angular_momentum(initial)
    for hours in (1, 12, 24, 48, 72):
        state = propagate(LEO, EPOCH + timedelta(hours=hours), NONE)
        assert abs(_energy(state) - e0) < 1e-3  # J/kg
        assert abs(_angular_momentum(state) - h0) < 1e-3  # m^2/s


def test_reference_propagator_self_convergence() -> None:
    when = EPOCH + timedelta(hours=24)
    fine = propagate_reference(LEO, when, step_s=10.0)
    coarse = propagate_reference(LEO, when, step_s=30.0)
    diff = math.sqrt(
        sum((fine.position[k] - coarse.position[k]) ** 2 for k in range(3))
    )
    assert diff < 200.0  # metres over 24 h with 10 s vs 30 s RK4 steps


@pytest.mark.parametrize("hours", [1, 24, 72])
def test_j2_analytic_matches_reference_within_tolerance(hours: int) -> None:
    """First-order secular J2 accuracy contract (docs/accuracy.md).

    The analytic secular-J2 model agrees with the independent J2 RK4
    reference within a bounded, documented envelope for typical LEO.
    """
    when = EPOCH + timedelta(hours=hours)
    analytic = propagate(LEO, when, J2)
    reference = propagate_reference(LEO, when, step_s=10.0)
    diff = math.sqrt(
        sum(
            (analytic.position[k] - reference.position[k]) ** 2
            for k in range(3)
        )
    )
    # Bounded secular/short-period residual: grows sub-linearly, <=25 km @72h.
    assert diff < 25_000.0


def test_high_eccentricity_no_nan() -> None:
    el = ClassicalElements(
        EPOCH,
        65_000_000.0,
        0.85,
        math.radians(63.4),
        0.2,
        0.3,
        2.1,
    )
    for hours in (0, 6, 24):
        state = propagate(el, EPOCH + timedelta(hours=hours), J2)
        assert all(math.isfinite(c) for c in state.position + state.velocity)


def test_two_body_forward_backward_returns_to_epoch() -> None:
    """Negative-duration propagation is the exact inverse of positive."""
    middle = EPOCH + timedelta(hours=4)
    forward = propagate(LEO, middle, NONE)
    n = math.sqrt(EARTH_MU / LEO.semi_major_axis_m**3)
    # Elements at the middle epoch carry the advanced mean anomaly; propagate
    # backward by 4 h to recover the original epoch state.
    middle_el = ClassicalElements(
        middle,
        LEO.semi_major_axis_m,
        LEO.eccentricity,
        LEO.inclination,
        LEO.raan,
        LEO.arg_perigee,
        LEO.mean_anomaly + n * 4 * 3600,
    )
    back = propagate(middle_el, EPOCH, NONE)
    initial = propagate(LEO, EPOCH, NONE)
    for k in range(3):
        assert back.position[k] == pytest.approx(initial.position[k], abs=1e-5)
        assert back.velocity[k] == pytest.approx(initial.velocity[k], abs=1e-9)


def test_non_finite_elements_rejected() -> None:
    el = ClassicalElements(
        EPOCH, float("nan"), 0.001, 0.4, 0.1, 0.2, 0.3
    )
    with pytest.raises(OrbitError) as exc:
        from app.elements import validate_elements

        validate_elements(el)
    assert exc.value.code is ErrorCode.NON_PHYSICAL_ORBIT
