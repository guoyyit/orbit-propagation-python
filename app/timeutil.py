"""Time parsing and conversions (UTC, ISO 8601, Unix seconds, Julian date).

Conventions (see docs/conventions.md):

* All civil input times are UTC. ISO 8601 strings may carry an explicit UTC
  offset (``Z`` or ``+HH:MM``); naive strings are assumed UTC.
* Unix timestamps are seconds since 1970-01-01T00:00:00Z (float, may be
  negative).
* UT1 is approximated by UTC. TT - UTC1 is the constant 69.184 s, matching
  the J2000 leap-second environment; this is accurate to a few ms for the
  1972-2100 service window and documented as a known limitation.
"""
from __future__ import annotations

import math
from datetime import datetime, timezone

from .errors import ErrorCode, OrbitError

J2000_JD = 2451545.0
TT_MINUS_UTC_S = 69.184


def parse_time(value: str | int | float | datetime) -> datetime:
    """Parse ISO 8601 string or Unix timestamp into an aware UTC datetime."""
    if isinstance(value, datetime):
        ts = value
        if ts.tzinfo is None:
            return ts.replace(tzinfo=timezone.utc)
        return ts.astimezone(timezone.utc)
    if isinstance(value, bool):
        raise OrbitError(
            ErrorCode.TIME_FORMAT_ERROR,
            "boolean is not a valid time; expected ISO 8601 string or Unix seconds",
        )
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise OrbitError(
                ErrorCode.TIME_FORMAT_ERROR,
                "Unix timestamp must be finite",
            )
        try:
            return datetime.fromtimestamp(float(value), tz=timezone.utc)
        except (OverflowError, OSError, ValueError) as exc:
            raise OrbitError(
                ErrorCode.TIME_FORMAT_ERROR, f"Unix timestamp out of range: {exc}"
            ) from exc
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise OrbitError(
                ErrorCode.TIME_FORMAT_ERROR, "empty time string"
            )
        try:
            parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
        except ValueError as exc:
            raise OrbitError(
                ErrorCode.TIME_FORMAT_ERROR,
                f"invalid ISO 8601 time {value!r}: {exc}",
            ) from exc
        if parsed.tzinfo is None:
            return parsed.replace(tzinfo=timezone.utc)
        return parsed.astimezone(timezone.utc)
    raise OrbitError(
        ErrorCode.TIME_FORMAT_ERROR,
        f"unsupported time type {type(value).__name__}",
    )


def to_unix(value: str | int | float | datetime) -> float:
    """Convert ISO 8601 / Unix timestamp / datetime to Unix seconds (UTC)."""
    if isinstance(value, bool):
        raise OrbitError(ErrorCode.TIME_FORMAT_ERROR, "boolean is not a valid time")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if not math.isfinite(float(value)):
            raise OrbitError(ErrorCode.TIME_FORMAT_ERROR, "timestamp must be finite")
        return float(value)
    return parse_time(value).timestamp()


def datetime_to_jd(ts: datetime) -> float:
    """Convert an aware UTC datetime to a Julian date (UT1 = UTC)."""
    if ts.tzinfo is None:
        ts = ts.replace(tzinfo=timezone.utc)
    ts = ts.astimezone(timezone.utc)
    year = ts.year
    month = ts.month
    if month <= 2:
        year -= 1
        month += 12
    a = year // 100
    b = 2 - a + a // 4
    day_fraction = (
        ts.day
        + (ts.hour + ts.minute / 60.0 + (ts.second + ts.microsecond / 1e6) / 3600.0)
        / 24.0
    )
    return (
        math.floor(365.25 * (year + 4716))
        + math.floor(30.6001 * (month + 1))
        + day_fraction
        + b
        - 1524.5
    )


def unix_to_jd(unix_seconds: float) -> float:
    """Convert Unix seconds to Julian date."""
    return datetime_to_jd(datetime.fromtimestamp(unix_seconds, tz=timezone.utc))


def julian_centuries_tt(jd_ut1: float) -> float:
    """Julian centuries of Terrestrial Time since J2000 (TT - UTC = 69.184 s)."""
    jd_tt = jd_ut1 + TT_MINUS_UTC_S / 86400.0
    return (jd_tt - J2000_JD) / 36525.0


def greenwich_sidereal_time(jd_ut1: float) -> float:
    """Greenwich Mean Sidereal Time in radians (0..2*pi) at the given JD.

    Vallado, *Fundamentals of Astrodynamics*, eq. 3-45 (IAU 1982 model).
    """
    t = (jd_ut1 - J2000_JD) / 36525.0
    gmst_seconds = (
        67310.54841
        + (876600.0 * 3600.0 + 8640184.812866) * t
        + 0.093104 * t * t
        - 6.2e-6 * t * t * t
    )
    gmst_deg = (gmst_seconds % 86400.0) / 240.0
    return math.radians(gmst_deg % 360.0)
