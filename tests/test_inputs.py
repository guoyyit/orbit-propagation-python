"""Input validation, error-code discrimination and TLE/COE equivalence."""

from __future__ import annotations

import math

import pytest

from orbit_service.constants import TLE_COE_AGREEMENT_METRES, WGS84
from orbit_service.elements import coe_to_eci
from orbit_service.errors import (
    FieldRangeError,
    NonPhysicalOrbitError,
    TLEChecksumError,
    TLEFormatError,
    TimeFormatError,
)
from orbit_service.time_utils import parse_time
from orbit_service.tle import checksum, parse_tle, tle_to_elements


def test_iso_and_unix_time_agree():
    iso = parse_time("2024-01-01T00:00:00Z")
    unix = parse_time(1_704_067_200)
    assert abs(iso.unix_seconds - unix.unix_seconds) < 1e-6


@pytest.mark.parametrize("bad", ["nonsense", "", "2024-13-40T00:00:00Z", float("nan")])
def test_bad_times_rejected(bad):
    with pytest.raises(TimeFormatError):
        parse_time(bad)


def test_tle_checksum_location(valid_tle_lines):
    line1, line2 = valid_tle_lines
    broken = line1[:-1] + str((int(line1[-1]) + 1) % 10)
    with pytest.raises(TLEChecksumError) as info:
        parse_tle(broken, line2)
    assert info.value.location == "line1:69"


def test_tle_format_error_location(valid_tle_lines):
    _, line2 = valid_tle_lines
    with pytest.raises(TLEFormatError):
        parse_tle("too short", line2)


def test_tle_field_range_location(valid_tle_lines):
    line1, line2 = valid_tle_lines
    broken = list(line2)
    broken[8:16] = list("190.0000")
    broken = "".join(broken)
    broken = broken[:-1] + str(checksum(broken))
    with pytest.raises(FieldRangeError) as info:
        parse_tle(line1, broken)
    assert "line2:9-16" in info.value.location


def test_tle_and_coe_equivalent(valid_tle_lines):
    tle = parse_tle(*valid_tle_lines)
    elements = tle_to_elements(tle)
    # The same orbit expressed as COE (mean motion converted back to axis by the
    # TLE parser) must be identical within the documented tolerance.
    n = tle.mean_motion_rad_per_s
    axis = (WGS84.mu / (n * n)) ** (1 / 3)
    assert abs(axis - elements.semi_major_axis) < TLE_COE_AGREEMENT_METRES
    position, velocity = coe_to_eci(elements)
    assert all(math.isfinite(v) for v in position + velocity)


@pytest.mark.parametrize(
    "eccentricity",
    [-0.1, 1.0, 2.0],
)
def test_eccentricity_bounds(eccentricity):
    with pytest.raises(NonPhysicalOrbitError):
        parse_time_value = 0
        from orbit_service.elements import ClassicalElements

        elements = ClassicalElements(
            7_000_000.0, eccentricity, 0.1, 0.0, 0.0, 0.0, parse_time(0)
        )
        coe_to_eci(elements)


def test_perigee_below_earth_rejected():
    from orbit_service.elements import ClassicalElements

    elements = ClassicalElements(
        6_500_000.0, 0.5, 0.1, 0.0, 0.0, 0.0, parse_time(0)
    )
    with pytest.raises(NonPhysicalOrbitError):
        coe_to_eci(elements)
