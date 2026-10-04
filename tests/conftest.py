"""Shared fixtures and independent reference integrators for tests."""

from __future__ import annotations

import math

import pytest

from orbit_service.constants import WGS84
from orbit_service.elements import ClassicalElements, coe_to_eci
from orbit_service.propagator import (
    J2Stepper,
    PropagatorConfig,
    _j2_acceleration,
    _mean_to_osc,
)
from orbit_service.time_utils import Epoch


@pytest.fixture
def iss_elements() -> ClassicalElements:
    return ClassicalElements(
        semi_major_axis=6_794_900.0,
        eccentricity=0.0015,
        inclination=math.radians(51.64),
        raan=math.radians(120.0),
        argument_of_periapsis=math.radians(60.0),
        mean_anomaly=0.0,
        epoch=Epoch(1_700_000_000.0),
    )


@pytest.fixture
def valid_tle_lines():
    from orbit_service.tle import checksum

    line1 = "1 25544U 98067A   24001.50000000  .00012345  12345-6  22222-3 0  9999"
    line2 = "2 25544  51.6400 100.0000 0001234  90.0000 270.0000 15.50000000  1019"
    line1 = line1[:-1] + str(checksum(line1))
    line2 = line2[:-1] + str(checksum(line2))
    return line1, line2


def independent_rk4_j2(elements: ClassicalElements, duration_seconds: float, dt: float):
    """Reference J2 trajectory integrated with an independent RK4 loop."""

    def rk4_step(state, step):
        def derivative(s):
            ax, ay, az = _j2_acceleration(s[0], s[1], s[2])
            return (s[3], s[4], s[5], ax, ay, az)

        k1 = derivative(state)
        s2 = tuple(state[i] + 0.5 * step * k1[i] for i in range(6))
        k2 = derivative(s2)
        s3 = tuple(state[i] + 0.5 * step * k2[i] for i in range(6))
        k3 = derivative(s3)
        s4 = tuple(state[i] + step * k3[i] for i in range(6))
        k4 = derivative(s4)
        return tuple(
            state[i] + step / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
            for i in range(6)
        )

    position, velocity = _mean_to_osc(*coe_to_eci(elements))
    state = (*position, *velocity)
    steps = max(1, int(round(duration_seconds / dt)))
    step = duration_seconds / steps
    for _ in range(steps):
        state = rk4_step(state, step)
    return state
