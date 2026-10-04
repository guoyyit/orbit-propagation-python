"""Two-Line Element (TLE) parsing, validation and conversion.

Reference: NORAD/SPACETRACK REPORT NO. 3, models 1980/2004 ("TLE-3").

Validation failures are reported with distinguishable codes:

* ``TLE_FORMAT_ERROR``   - wrong line count/length/structure
* ``TLE_CHECKSUM_ERROR`` - NORAD modulo-10 checksum mismatch
* ``TLE_FIELD_ERROR``    - a numeric/field is out of its allowed range

Angles are converted to radians and rates to SI units.  Conversion to
classical elements uses the WGS-72 constants embedded in TLE mean elements;
the resulting semi-major axis is reported in metres and mean motion in rad/s.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone

from .elements import (
    WGS72_MU,
    ClassicalElements,
    mean_motion_to_sma,
)
from .errors import ErrorCode, OrbitError

TLE_LINE_LEN = 69
TLE_EPOCH_REF_YEAR = 2000  # two-digit year: >=57 -> 1900s, else 2000s


@dataclass(frozen=True)
class TLE:
    """Parsed TLE contents (angles rad, rates rad/s where applicable)."""

    name: str
    satellite_number: int
    classification: str
    launch_year: int
    launch_number: int
    launch_piece: str
    epoch_year: int
    epoch_doy: float
    inclination: float
    raan: float
    eccentricity: float
    arg_perigee: float
    mean_anomaly: float
    mean_motion_rev_per_day: float
    ndot_over_2: float
    nddot_over_6: float
    bstar: float
    ephemeris_type: int
    element_set_number: int

    @property
    def epoch(self) -> datetime:
        """UTC epoch of the TLE."""
        return _epoch_datetime(self.epoch_year, self.epoch_doy)


def tle_checksum(line: str) -> int:
    """NORAD checksum: digit sum with '-' counting as 1, modulo 10."""
    total = 0
    for char in line[:-1]:
        if char.isdigit():
            total += int(char)
        elif char == "-":
            total += 1
    return total % 10


def parse_tle(line0: str, line1: str, line2: str) -> TLE:
    """Parse and fully validate a 3-line TLE set."""
    l0 = line0.rstrip("\n")
    l1 = line1.rstrip("\n")
    l2 = line2.rstrip("\n")
    if l1[:1] != "1" or l2[:1] != "2":
        raise OrbitError(
            ErrorCode.TLE_FORMAT_ERROR,
            "line numbers must be '1' and '2' in column 1",
            location="line1/line2",
        )
    if len(l1) != TLE_LINE_LEN or len(l2) != TLE_LINE_LEN:
        raise OrbitError(
            ErrorCode.TLE_FORMAT_ERROR,
            f"data lines must be exactly {TLE_LINE_LEN} columns; "
            f"got {len(l1)} and {len(l2)}",
            location="line1/line2",
        )
    _check_checksum(l1, "line1")
    _check_checksum(l2, "line2")

    satnum = _parse_int(l1[2:7], "satellite_number", "line1")
    classification = l1[7]
    if classification not in ("U", "C", "S"):
        raise OrbitError(
            ErrorCode.TLE_FIELD_ERROR,
            f"classification must be U/C/S, got {classification!r}",
            location="line1.classification",
        )
    launch_year = _parse_int(l1[9:11], "launch_year", "line1")
    launch_number = _parse_int(l1[11:14], "launch_number", "line1")
    launch_piece = l1[14:17].strip()
    epoch_year = _parse_int(l1[18:20], "epoch_year", "line1")
    epoch_doy = _parse_float(l1[20:32], "epoch_day", "line1")
    if not (1.0 <= epoch_doy <= 366.99999999):
        raise OrbitError(
            ErrorCode.TLE_FIELD_ERROR,
            f"epoch day of year out of range: {epoch_doy}",
            location="line1.epoch",
        )
    ndot = _parse_float(l1[33:43], "first_derivative_mean_motion", "line1")
    nddot = _parse_decimal_assumption(l1[44:52], "second_derivative", "line1")
    bstar = _parse_decimal_assumption(l1[53:61], "bstar", "line1")
    ephemeris_type = _parse_int(l1[62:63], "ephemeris_type", "line1")
    element_set = _parse_int(l1[64:68], "element_set_number", "line1")

    satnum2 = _parse_int(l2[2:7], "satellite_number", "line2")
    if satnum2 != satnum:
        raise OrbitError(
            ErrorCode.TLE_FORMAT_ERROR,
            f"satellite number mismatch: {satnum} vs {satnum2}",
            location="line2.satellite_number",
        )
    inclination = math.radians(
        _parse_float(l2[8:16], "inclination", "line2")
    )
    raan = math.radians(_parse_float(l2[17:25], "raan", "line2"))
    eccentricity = _parse_implicit_decimal(l2[26:33], "eccentricity", "line2")
    if not (0.0 <= eccentricity < 1.0):
        raise OrbitError(
            ErrorCode.TLE_FIELD_ERROR,
            f"TLE eccentricity out of [0,1): {eccentricity}",
            location="line2.eccentricity",
        )
    arg_perigee = math.radians(
        _parse_float(l2[34:42], "arg_perigee", "line2")
    )
    mean_anomaly = math.radians(
        _parse_float(l2[43:51], "mean_anomaly", "line2")
    )
    mean_motion = _parse_float(l2[52:63], "mean_motion", "line2")
    if mean_motion <= 0.0:
        raise OrbitError(
            ErrorCode.TLE_FIELD_ERROR,
            f"mean motion must be positive: {mean_motion}",
            location="line2.mean_motion",
        )

    return TLE(
        name=l0.strip(),
        satellite_number=satnum,
        classification=classification,
        launch_year=1900 + launch_year if launch_year >= 57 else 2000 + launch_year,
        launch_number=launch_number,
        launch_piece=launch_piece,
        epoch_year=1900 + epoch_year if epoch_year >= 57 else 2000 + epoch_year,
        epoch_doy=epoch_doy,
        inclination=inclination,
        raan=raan,
        eccentricity=eccentricity,
        arg_perigee=arg_perigee,
        mean_anomaly=mean_anomaly,
        mean_motion_rev_per_day=mean_motion,
        ndot_over_2=ndot,
        nddot_over_6=nddot,
        bstar=bstar,
        ephemeris_type=ephemeris_type,
        element_set_number=element_set,
    )


def tle_to_elements(tle: TLE) -> ClassicalElements:
    """Convert TLE mean elements (WGS-72) to classical elements (SI units)."""
    mean_motion_rad_s = tle.mean_motion_rev_per_day * 2.0 * math.pi / 86400.0
    semi_major_axis = mean_motion_to_sma(mean_motion_rad_s, WGS72_MU)
    return ClassicalElements(
        epoch=tle.epoch,
        semi_major_axis_m=semi_major_axis,
        eccentricity=tle.eccentricity,
        inclination=tle.inclination,
        raan=tle.raan,
        arg_perigee=tle.arg_perigee,
        mean_anomaly=tle.mean_anomaly,
        mu=WGS72_MU,
    )


def _check_checksum(line: str, where: str) -> None:
    expected = _parse_int(line[68:69], "checksum", where)
    actual = tle_checksum(line)
    if actual != expected:
        raise OrbitError(
            ErrorCode.TLE_CHECKSUM_ERROR,
            f"checksum mismatch: expected {expected}, computed {actual}",
            location=f"{where}.checksum",
        )


def _parse_int(text: str, field: str, where: str) -> int:
    try:
        return int(text)
    except ValueError as exc:
        raise OrbitError(
            ErrorCode.TLE_FORMAT_ERROR,
            f"invalid integer field {field}: {text!r}",
            location=f"{where}.{field}",
        ) from exc


def _parse_float(text: str, field: str, where: str) -> float:
    try:
        return float(text)
    except ValueError as exc:
        raise OrbitError(
            ErrorCode.TLE_FORMAT_ERROR,
            f"invalid numeric field {field}: {text!r}",
            location=f"{where}.{field}",
        ) from exc


def _parse_implicit_decimal(text: str, field: str, where: str) -> float:
    """Parse an assumed leading decimal point (TLE eccentricity)."""
    try:
        value = float("0." + text.strip())
    except ValueError as exc:
        raise OrbitError(
            ErrorCode.TLE_FORMAT_ERROR,
            f"invalid decimal field {field}: {text!r}",
            location=f"{where}.{field}",
        ) from exc
    return value


def _parse_decimal_assumption(text: str, field: str, where: str) -> float:
    """Parse the NORAD 'assumed decimal + exponent' field (bstar, nddot).

    Format (8 columns)::

        [sign][digit][digit digit digit digit][sign][exponent]

    with an implicit decimal point after the first digit, e.g.
    ``-11606-4`` means ``-0.11606 e-4`` = -1.1606e-5 and ``12345-3``
    means ``0.12345 e-3``.
    """
    cleaned = text.strip()
    if cleaned in ("", "0", "00000", "00000-0"):
        return 0.0
    try:
        sign = 1.0
        if cleaned[0] in "+-":
            sign = -1.0 if cleaned[0] == "-" else 1.0
            cleaned = cleaned[1:]
        if len(cleaned) < 2:
            raise ValueError("missing exponent")
        exp_sign = 1.0
        if cleaned[-2] in "+-":
            exp_sign = -1.0 if cleaned[-2] == "-" else 1.0
            digits = cleaned[:-2]
            exp_digits = cleaned[-1]
        else:
            # No explicit exponent sign: older packing (rare); assume positive.
            digits = cleaned[:-1]
            exp_digits = cleaned[-1]
        mantissa = float("0." + digits.strip())
        exponent = exp_sign * int(exp_digits)
        return sign * mantissa * 10.0**exponent
    except ValueError as exc:
        raise OrbitError(
            ErrorCode.TLE_FORMAT_ERROR,
            f"invalid exponential field {field}: {text!r}",
            location=f"{where}.{field}",
        ) from exc


def _epoch_datetime(year: int, doy: float) -> datetime:
    start = datetime(year, 1, 1, tzinfo=timezone.utc)
    whole = int(math.floor(doy)) - 1
    fraction = doy - math.floor(doy)
    return start + timedelta(days=whole, seconds=fraction * 86400.0)
