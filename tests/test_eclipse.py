"""Earth shadow events and cylindrical vs conical model difference."""

from __future__ import annotations

from orbit_service.eclipse import (
    ShadowModel,
    compute_illumination,
    is_illuminated,
    sun_position_eci,
)
from orbit_service.time_utils import Epoch


def test_events_alternate_and_are_sorted(iss_elements):
    start = Epoch(iss_elements.epoch.unix_seconds)
    end = Epoch(start.unix_seconds + 86400)
    result = compute_illumination(
        iss_elements, start, end, ShadowModel.CYLINDRICAL, 1e-3
    )
    times = [event.time.unix_seconds for event in result.events]
    assert times == sorted(times)
    kinds = [event.kind for event in result.events]
    for index, kind in enumerate(kinds):
        expected = "entry" if index % 2 == 0 else "exit"
        assert kind == expected


def test_conical_eclipses_are_shorter_than_cylindrical(iss_elements):
    start = Epoch(iss_elements.epoch.unix_seconds)
    end = Epoch(start.unix_seconds + 86400)
    cylindrical = compute_illumination(
        iss_elements, start, end, ShadowModel.CYLINDRICAL, 1e-3
    )
    conical = compute_illumination(
        iss_elements, start, end, ShadowModel.CONICAL, 1e-3
    )
    assert conical.illuminated_seconds >= cylindrical.illuminated_seconds
    assert len(conical.events) == len(cylindrical.events)


def test_boundary_tolerance_converges(iss_elements):
    start = Epoch(iss_elements.epoch.unix_seconds)
    end = Epoch(start.unix_seconds + 86400)
    coarse = compute_illumination(
        iss_elements, start, end, ShadowModel.CONICAL, 1.0
    )
    fine = compute_illumination(
        iss_elements, start, end, ShadowModel.CONICAL, 1e-4
    )
    assert len(coarse.events) == len(fine.events)
    for first, second in zip(coarse.events, fine.events):
        assert abs(first.time.unix_seconds - second.time.unix_seconds) <= 1.0


def test_sun_ephemeris_is_near_one_au(iss_elements):
    import math

    sun = sun_position_eci(iss_elements.epoch.unix_seconds)
    distance = math.sqrt(sum(component**2 for component in sun))
    assert 0.97 * 1.495978707e11 < distance < 1.03 * 1.495978707e11
    # Point on the sun-lit side of Earth at low altitude is illuminated under
    # every model.
    sun_hat = tuple(component / distance for component in sun)
    position = tuple(component * 6800e3 for component in sun_hat)
    assert is_illuminated(position, iss_elements.epoch.unix_seconds, ShadowModel.CYLINDRICAL)
    assert is_illuminated(position, iss_elements.epoch.unix_seconds, ShadowModel.CONICAL)
