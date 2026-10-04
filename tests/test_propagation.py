"""Reference orbit agreement, conservation and degenerate-orbit robustness."""

from __future__ import annotations

import math

from orbit_service.constants import (
    J2_REFERENCE_METRES,
    TWOBODY_ANGMOM_DRIFT,
    TWOBODY_ENERGY_DRIFT,
)
from orbit_service.elements import (
    angular_momentum,
    coe_to_eci,
    specific_energy,
)
from orbit_service.propagator import (
    PropagatorConfig,
    PropagationMode,
    propagate,
    propagate_two_body,
)
from orbit_service.time_utils import Epoch

from .conftest import independent_rk4_j2


def test_j2_matches_independent_reference(iss_elements):
    duration = 2 * 86400
    target = Epoch(iss_elements.epoch.unix_seconds + duration)
    config = PropagatorConfig(step_seconds=15.0)
    state = propagate(iss_elements, target, config)
    reference = independent_rk4_j2(iss_elements, duration, 1.0)
    distance = math.sqrt(
        sum((state.position[i] - reference[i]) ** 2 for i in range(3))
    )
    assert distance < J2_REFERENCE_METRES


def test_iso_unix_time_equivalence_for_propagation(iss_elements):
    from orbit_service.time_utils import iso_from_unix, parse_time

    unix_target = iss_elements.epoch.unix_seconds + 12345.0
    via_unix = propagate(iss_elements, parse_time(unix_target))
    via_iso = propagate(iss_elements, parse_time(iso_from_unix(unix_target)))
    assert via_unix.position == via_iso.position


def test_two_body_energy_and_momentum_conserved(iss_elements):
    p0, v0 = coe_to_eci(iss_elements)
    e0 = specific_energy(p0, v0)
    h0 = angular_momentum(p0, v0)
    h0_mag = math.sqrt(sum(x * x for x in h0))
    target = Epoch(iss_elements.epoch.unix_seconds + 7 * 86400)
    state = propagate(iss_elements, target, PropagatorConfig(mode=PropagationMode.TWO_BODY))
    e1 = specific_energy(state.position, state.velocity)
    h1 = angular_momentum(state.position, state.velocity)
    h1_mag = math.sqrt(sum(x * x for x in h1))
    assert abs((e1 - e0) / e0) < TWOBODY_ENERGY_DRIFT
    assert abs((h1_mag - h0_mag) / h0_mag) < TWOBODY_ANGMOM_DRIFT


def test_degenerate_orbits_have_no_nan(iss_elements):
    from orbit_service.elements import ClassicalElements

    epoch = iss_elements.epoch
    cases = [
        (26_600_000.0, 0.7, math.radians(30)),
        (7_000_000.0, 1e-8, 0.0),
        (8_000_000.0, 1e-6, math.radians(179.9)),
    ]
    for axis, ecc, inc in cases:
        elements = ClassicalElements(axis, ecc, inc, 0.3, 0.7, 1.1, epoch)
        state = propagate(
            elements,
            Epoch(epoch.unix_seconds + 86400),
            PropagatorConfig(step_seconds=10.0),
        )
        assert all(math.isfinite(v) for v in state.position + state.velocity)
        assert math.sqrt(sum(x * x for x in state.position)) > 6_371_000


def test_propagation_is_bitwise_deterministic(iss_elements):
    target = Epoch(iss_elements.epoch.unix_seconds + 50000.0)
    first = propagate(iss_elements, target)
    second = propagate(iss_elements, target)
    assert first.position == second.position
    assert first.velocity == second.velocity
