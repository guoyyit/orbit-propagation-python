"""Tests for TLE parsing, checksum validation and element conversion."""
from __future__ import annotations

import math

import pytest

from app.errors import ErrorCode, OrbitError
from app.tle import parse_tle, tle_checksum, tle_to_elements

LINE0 = "ISS (ZARYA)"
LINE1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
LINE2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"


def test_checksum_value() -> None:
    assert tle_checksum(LINE1) == int(LINE1[-1])
    assert tle_checksum(LINE2) == int(LINE2[-1])


def test_parse_valid_tle() -> None:
    tle = parse_tle(LINE0, LINE1, LINE2)
    assert tle.satellite_number == 25544
    assert tle.classification == "U"
    assert tle.eccentricity == pytest.approx(0.0006703, abs=1e-12)
    assert tle.mean_motion_rev_per_day == pytest.approx(15.72125391)
    assert tle.bstar == pytest.approx(-1.1606e-5)
    assert tle.epoch.year == 2008


def test_tle_to_elements_semi_major_axis() -> None:
    tle = parse_tle(LINE0, LINE1, LINE2)
    el = tle_to_elements(tle)
    # ISS perigee altitude ~348 km (WGS-72).
    perigee = el.semi_major_axis_m * (1 - el.eccentricity)
    assert (perigee - 6_378_135.0) / 1000 == pytest.approx(348.3, abs=1.0)
    assert math.degrees(el.inclination) == pytest.approx(51.6416, abs=1e-6)


def test_bad_checksum_distinguished() -> None:
    bad_line2 = LINE2[:-1] + "0"
    with pytest.raises(OrbitError) as exc:
        parse_tle(LINE0, LINE1, bad_line2)
    assert exc.value.code is ErrorCode.TLE_CHECKSUM_ERROR
    assert "checksum" in exc.value.location


def test_bad_structure_distinguished() -> None:
    with pytest.raises(OrbitError) as exc:
        parse_tle(LINE0, LINE1[:-1], LINE2)
    assert exc.value.code is ErrorCode.TLE_FORMAT_ERROR


def test_line_number_mismatch() -> None:
    with pytest.raises(OrbitError) as exc:
        parse_tle(LINE0, "3" + LINE1[1:], LINE2)
    assert exc.value.code is ErrorCode.TLE_FORMAT_ERROR


def test_satellite_number_mismatch() -> None:
    bad_line2 = "2 25545  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"
    # rebuild correct checksum so only the number mismatch is tested
    bad_line2 = bad_line2[:68] + str(tle_checksum(bad_line2))
    with pytest.raises(OrbitError) as exc:
        parse_tle(LINE0, LINE1, bad_line2)
    assert exc.value.code is ErrorCode.TLE_FORMAT_ERROR


def test_out_of_range_eccentricity() -> None:
    # The implicit-decimal eccentricity field cannot encode >= 1 in 7
    # columns, so an out-of-range value is exercised via an invalid digit.
    broken = LINE2[:26] + "X000000" + LINE2[33:]
    broken = broken[:68] + str(tle_checksum(broken))
    with pytest.raises(OrbitError) as exc:
        parse_tle(LINE0, LINE1, broken)
    assert exc.value.code is ErrorCode.TLE_FORMAT_ERROR


def _build_tle(
    epoch_year: int,
    epoch_doy: float,
    inc_deg: float,
    raan_deg: float,
    eccentricity: float,
    argp_deg: float,
    mean_anom_deg: float,
    mean_motion_rev_day: float,
) -> tuple[str, str, str]:
    """Format a syntactically valid TLE from individual field values."""

    def checksum(line: str) -> str:
        return str(tle_checksum(line))

    line1 = (
        f"1 70000U 24001A   {epoch_year % 100:02d}{epoch_doy:012.8f} "
        f"{0.0:9.8f}  00000-0  00000-0 0  999"
    )
    line2 = (
        f"2 70000 {inc_deg:8.4f} {raan_deg:8.4f} {eccentricity * 1e7:07.0f}"
        f" {argp_deg:8.4f} {mean_anom_deg:8.4f} {mean_motion_rev_day:11.8f} 9999"
    )
    line1 = line1[:68] + str(tle_checksum(line1[:68] + "0"))
    line2 = line2[:68] + str(tle_checksum(line2[:68] + "0"))
    return ("SAT", line1, line2)


def test_tle_and_classical_elements_describe_same_orbit() -> None:
    """The same orbit given as TLE and classical elements converts to the
    same state (WGS-72 constants, same epoch)."""
    from datetime import datetime, timezone

    from app.elements import (
        WGS72_MU,
        ClassicalElements,
        elements_to_state,
        sma_to_mean_motion,
    )

    epoch = datetime(2024, 3, 20, 12, tzinfo=timezone.utc)
    a = 7_100_000.0
    e = 0.0123
    inc = math.radians(51.6)
    raan = math.radians(123.45)
    argp = math.radians(67.89)
    mean_anom = math.radians(234.56)

    n_rad_s = sma_to_mean_motion(a, WGS72_MU)
    n_rev_day = n_rad_s * 86400.0 / (2.0 * math.pi)

    lines = _build_tle(
        2024, 80.5, 51.6, 123.45, e, 67.89, 234.56, n_rev_day
    )
    tle = parse_tle(*lines)
    # Epoch matches day 80.5 of 2024 at the minute level the fields can express.
    assert tle.epoch == epoch

    from_tle = tle_to_elements(tle)
    classical = ClassicalElements(
        epoch=epoch,
        semi_major_axis_m=a,
        eccentricity=e,
        inclination=inc,
        raan=raan,
        arg_perigee=argp,
        mean_anomaly=mean_anom,
        mu=WGS72_MU,
    )
    s_tle = elements_to_state(from_tle)
    s_cls = elements_to_state(classical)
    for k in range(3):
        assert s_tle.position[k] == pytest.approx(s_cls.position[k], abs=1.0)
        assert s_tle.velocity[k] == pytest.approx(s_cls.velocity[k], abs=1e-4)
