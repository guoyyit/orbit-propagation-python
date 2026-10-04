"""Two-line element set parsing and conversion to classical elements.

Only the standard 69-column two-line form is accepted. Checksums are verified
with the modulo-10 rule and every fixed-width field is range checked; a failure
raises a distinguishable error carrying the offending line/column location.
"""

from __future__ import annotations

from datetime import datetime, timezone
import math
from dataclasses import dataclass

from .constants import WGS84
from .elements import ClassicalElements
from .errors import FieldRangeError, TLEChecksumError, TLEFormatError
from .time_utils import Epoch


@dataclass(frozen=True)
class TLE:
    line1: str
    line2: str
    satellite_number: str
    classification: str
    international_designator: str
    epoch: Epoch
    mean_motion_first_derivative: float
    mean_motion_second_derivative: float
    bstar: float
    ephemeris_type: int
    element_set: int
    inclination_rad: float
    raan_rad: float
    eccentricity: float
    argument_of_periapsis_rad: float
    mean_anomaly_rad: float
    mean_motion_rad_per_s: float


def checksum(line: str) -> int:
    """Modulo-10 TLE checksum over the line content before the check digit."""
    total = 0
    for char in line[:-1]:
        if char.isdigit():
            total += int(char)
        elif char == "-":
            total += 1
    return total % 10


def _check_line(line: str, expected_prefix: str, line_number: int) -> None:
    location = f"line{line_number}"
    if len(line) != 69:
        raise TLEFormatError(
            f"line {line_number} must be exactly 69 columns, got {len(line)}",
            location=location,
        )
    if line[0] != expected_prefix:
        raise TLEFormatError(
            f"line {line_number} must start with {expected_prefix!r}, got {line[0]!r}",
            location=f"{location}:1",
        )
    digit = line[-1]
    if not digit.isdigit():
        raise TLEFormatError(
            f"line {line_number} checksum column must be a digit",
            location=f"{location}:69",
        )
    if int(digit) != checksum(line):
        raise TLEChecksumError(
            f"line {line_number} checksum mismatch", location=f"{location}:69"
        )


def _read_float(
    text: str, location: str, *, allow_implicit_decimal: bool = False
) -> float:
    try:
        if allow_implicit_decimal:
            cleaned = text.strip()
            if not cleaned:
                raise ValueError("empty field")
            if cleaned[0] in "+-":
                sign = cleaned[0]
                rest = cleaned[1:]
            else:
                sign = ""
                rest = cleaned
            exponent = 0
            if len(rest) >= 2 and rest[-2] in "+-":
                exponent = int(rest[-1])
                rest = rest[:-2]
            digits = rest.replace(".", "")
            if not digits.isdigit() or not digits:
                raise ValueError(f"bad implicit-decimal field {text!r}")
            mantissa = float("0." + digits)
            return float(sign + "1") * mantissa * 10.0**exponent
        return float(text)
    except ValueError as exc:
        raise TLEFormatError(f"not a number: {text.strip()!r}", location=location) from exc


def _read_int(text: str, location: str) -> int:
    try:
        return int(text)
    except ValueError as exc:
        raise TLEFormatError(f"not an integer: {text.strip()!r}", location=location) from exc


def _epoch_from_year_day(year: int, day_of_year: float) -> Epoch:
    if not 0 <= year <= 99:
        raise FieldRangeError(
            f"TLE epoch year {year} outside 00..99", location="line1:19-20"
        )
    full_year = 2000 + year if year < 57 else 1900 + year
    if not 0.0 < day_of_year <= 366.0:
        raise FieldRangeError(
            f"TLE epoch day {day_of_year} outside (0, 366]", location="line1:21-32"
        )
    base = datetime(full_year, 1, 1, tzinfo=timezone.utc)
    return Epoch(base.timestamp() + (day_of_year - 1.0) * 86400.0)


def parse_tle(line1: str, line2: str) -> TLE:
    if not isinstance(line1, str) or not isinstance(line2, str):
        raise TLEFormatError("TLE lines must be strings")
    line1 = line1.rstrip("\n\r")
    line2 = line2.rstrip("\n\r")
    _check_line(line1, "1", 1)
    _check_line(line2, "2", 2)

    satellite_number_1 = line1[2:7].strip()
    satellite_number_2 = line2[2:7].strip()
    if satellite_number_1 != satellite_number_2 or not satellite_number_1.isdigit():
        raise TLEFormatError(
            "satellite numbers must match and be numeric",
            location="line1:3-7",
        )

    classification = line1[7]
    if classification not in ("U", "C", "S"):
        raise FieldRangeError(
            f"classification {classification!r} not in U/C/S",
            location="line1:8",
        )

    year = _read_int(line1[18:20], "line1:19-20")
    day_of_year = _read_float(line1[20:32], "line1:21-32")
    epoch = _epoch_from_year_day(year, day_of_year)

    first_derivative = _read_float(line1[33:43], "line1:34-43")
    second_derivative = _read_float(
        line1[44:52], "line1:45-52", allow_implicit_decimal=True
    )
    bstar = _read_float(line1[53:61], "line1:54-61", allow_implicit_decimal=True)
    ephemeris_type = _read_int(line1[62:63], "line1:63")
    element_set = _read_int(line1[64:68], "line1:65-68")

    inclination_deg = _read_float(line2[8:16], "line2:9-16")
    raan_deg = _read_float(line2[17:25], "line2:18-25")
    eccentricity_raw = line2[26:33]
    if not eccentricity_raw.strip().isdigit():
        raise TLEFormatError(
            "eccentricity field must be 7 digits", location="line2:27-33"
        )
    eccentricity = int(eccentricity_raw) / 1e7
    argp_deg = _read_float(line2[34:42], "line2:35-42")
    mean_anomaly_deg = _read_float(line2[43:51], "line2:44-51")
    mean_motion_rev_day = _read_float(line2[52:63], "line2:53-63")

    if not 0.0 <= inclination_deg <= 180.0:
        raise FieldRangeError(
            f"inclination {inclination_deg} deg outside [0, 180]",
            location="line2:9-16",
        )
    if not 0.0 <= raan_deg < 360.0:
        raise FieldRangeError(
            f"RAAN {raan_deg} deg outside [0, 360)", location="line2:18-25"
        )
    if not 0.0 <= eccentricity < 1.0:
        raise FieldRangeError(
            f"eccentricity {eccentricity} outside [0, 1)", location="line2:27-33"
        )
    if not 0.0 <= argp_deg < 360.0:
        raise FieldRangeError(
            f"argument of perigee {argp_deg} deg outside [0, 360)",
            location="line2:35-42",
        )
    if not 0.0 <= mean_anomaly_deg < 360.0:
        raise FieldRangeError(
            f"mean anomaly {mean_anomaly_deg} deg outside [0, 360)",
            location="line2:44-51",
        )
    if mean_motion_rev_day <= 0.0 or not math.isfinite(mean_motion_rev_day):
        raise FieldRangeError(
            f"mean motion {mean_motion_rev_day} rev/day must be positive",
            location="line2:53-63",
        )

    return TLE(
        line1=line1,
        line2=line2,
        satellite_number=satellite_number_1,
        classification=classification,
        international_designator=line1[9:17].strip(),
        epoch=epoch,
        mean_motion_first_derivative=first_derivative,
        mean_motion_second_derivative=second_derivative,
        bstar=bstar,
        ephemeris_type=ephemeris_type,
        element_set=element_set,
        inclination_rad=math.radians(inclination_deg),
        raan_rad=math.radians(raan_deg),
        eccentricity=eccentricity,
        argument_of_periapsis_rad=math.radians(argp_deg),
        mean_anomaly_rad=math.radians(mean_anomaly_deg),
        mean_motion_rad_per_s=mean_motion_rev_day * 2.0 * math.pi / 86400.0,
    )


def tle_to_elements(tle: TLE) -> ClassicalElements:
    n = tle.mean_motion_rad_per_s
    semi_major_axis = (WGS84.mu / (n * n)) ** (1.0 / 3.0)
    min_radius = semi_major_axis * (1.0 - tle.eccentricity)
    if min_radius < WGS84.radius_equator:
        raise FieldRangeError(
            f"perigee radius {min_radius:.1f} m is inside the Earth",
            location="line2:53-63",
        )
    return ClassicalElements(
        semi_major_axis=semi_major_axis,
        eccentricity=tle.eccentricity,
        inclination=tle.inclination_rad,
        raan=tle.raan_rad,
        argument_of_periapsis=tle.argument_of_periapsis_rad,
        mean_anomaly=tle.mean_anomaly_rad,
        epoch=tle.epoch,
    )
