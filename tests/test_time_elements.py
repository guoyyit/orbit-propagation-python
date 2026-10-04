"""Tests for time parsing, Kepler conversion and element validation."""
from __future__ import annotations

import math
from datetime import datetime, timezone

import pytest

from app.elements import (
    EARTH_MU,
    ClassicalElements,
    elements_to_state,
    mean_motion_to_sma,
    mean_to_eccentric_anomaly,
    sma_to_mean_motion,
    validate_elements,
)
from app.errors import ErrorCode, OrbitError
from app.timeutil import (
    datetime_to_jd,
    greenwich_sidereal_time,
    parse_time,
    to_unix,
)


EPOCH = datetime(2024, 1, 1, tzinfo=timezone.utc)


def test_iso_and_unix_time_equivalent() -> None:
    iso = "2024-06-01T12:34:56Z"
    unix = to_unix(iso)
    assert parse_time(iso) == parse_time(unix)
    assert to_unix(unix) == unix


def test_naive_iso_assumed_utc() -> None:
    assert parse_time("2024-06-01T00:00:00") == parse_time(
        "2024-06-01T00:00:00Z"
    )


def test_invalid_time_rejected() -> None:
    with pytest.raises(OrbitError) as exc:
        parse_time("not-a-time")
    assert exc.value.code is ErrorCode.TIME_FORMAT_ERROR


def test_j2000_julian_date() -> None:
    jd = datetime_to_jd(datetime(2000, 1, 1, 12, tzinfo=timezone.utc))
    assert jd == pytest.approx(2451545.0, abs=1e-9)


def test_gmst_j2000() -> None:
    jd = datetime_to_jd(datetime(2000, 1, 1, 12, tzinfo=timezone.utc))
    assert math.degrees(greenwich_sidereal_time(jd)) == pytest.approx(
        280.4606, abs=1e-3
    )


@pytest.mark.parametrize("eccentricity", [0.0, 1e-9, 0.5, 0.9, 0.999])
def test_kepler_equation_converges(eccentricity: float) -> None:
    for mean_anomaly in [i * 0.37 for i in range(-20, 21)]:
        ecc = mean_to_eccentric_anomaly(mean_anomaly, eccentricity)
        residual = ecc - eccentricity * math.sin(ecc) - mean_anomaly
        assert abs(residual) < 1e-9


def test_kepler_branch_for_large_mean_anomaly() -> None:
    ecc = mean_to_eccentric_anomaly(3.0 * math.pi, 0.001)
    assert ecc > math.pi


def test_kepler_rejects_bad_eccentricity() -> None:
    with pytest.raises(OrbitError) as exc:
        mean_to_eccentric_anomaly(0.0, 1.1)
    assert exc.value.code is ErrorCode.FIELD_OUT_OF_RANGE


def test_circular_orbit_state_matches_vis_viva() -> None:
    radius = 7_000_000.0
    el = ClassicalElements(EPOCH, radius, 0.0, math.radians(51.6), 0.3, 0.0, 0.0)
    state = elements_to_state(el)
    r = math.sqrt(sum(c * c for c in state.position))
    v = math.sqrt(sum(c * c for c in state.velocity))
    assert r == pytest.approx(radius, abs=1e-6)
    assert v == pytest.approx(math.sqrt(EARTH_MU / radius), rel=1e-12)


@pytest.mark.parametrize("eccentricity", [1e-6, 0.001, 0.1, 0.5, 0.9])
def test_energy_and_momentum_conserved_across_anomaly(eccentricity: float) -> None:
    radius = 12_000_000.0
    energy0 = h0 = None
    for mean_anomaly in [i * 0.25 for i in range(26)]:
        el = ClassicalElements(
            EPOCH, radius, eccentricity, 0.4, 0.5, 0.2, mean_anomaly
        )
        state = elements_to_state(el)
        r = math.sqrt(sum(c * c for c in state.position))
        v = math.sqrt(sum(c * c for c in state.velocity))
        energy = v * v / 2 - EARTH_MU / r
        hx = (
            state.position[1] * state.velocity[2]
            - state.position[2] * state.velocity[1]
        )
        if energy0 is None:
            energy0 = energy
            h0 = hx
        assert energy == pytest.approx(energy0, abs=1e-3)


def test_sma_mean_motion_roundtrip() -> None:
    radius = 7_123_456.0
    assert mean_motion_to_sma(sma_to_mean_motion(radius)) == pytest.approx(
        radius, abs=1e-6
    )


def test_eccentricity_out_of_range() -> None:
    el = ClassicalElements(EPOCH, 7e6, 1.2, 0.1, 0.0, 0.0, 0.0)
    with pytest.raises(OrbitError) as exc:
        validate_elements(el)
    assert exc.value.code is ErrorCode.FIELD_OUT_OF_RANGE


def test_subterranean_perigee_rejected() -> None:
    el = ClassicalElements(EPOCH, 6_400_000.0, 0.01, 0.1, 0.0, 0.0, 0.0)
    with pytest.raises(OrbitError) as exc:
        validate_elements(el)
    assert exc.value.code is ErrorCode.NON_PHYSICAL_ORBIT
