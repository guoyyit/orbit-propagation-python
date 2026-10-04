"""Time parsing and time-system helpers.

Every internal instant is a POSIX timestamp (seconds since
1970-01-01T00:00:00Z), a continuous float. UTC is treated as UT1 for the
sidereal-time conversion; leap seconds are ignored (POSIX time scale), which is
accurate well below the service's metre-level propagation tolerances.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timezone

from .errors import TimeFormatError


@dataclass(frozen=True)
class Epoch:
    """An instant expressed as Unix seconds (UTC, POSIX, leap-second agnostic)."""

    unix_seconds: float


def parse_time(value: str | int | float) -> Epoch:
    if isinstance(value, bool):
        raise TimeFormatError("boolean is not a valid time", location="time")
    if isinstance(value, (int, float)):
        seconds = float(value)
        if not math.isfinite(seconds):
            raise TimeFormatError("timestamp must be finite", location="time")
        return Epoch(seconds)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise TimeFormatError("empty time string", location="time")
        try:
            if text.endswith(("Z", "z")):
                text = text[:-1] + "+00:00"
            parsed = datetime.fromisoformat(text)
        except ValueError as exc:
            raise TimeFormatError(
                f"not an ISO 8601 instant: {value!r}", location="time"
            ) from exc
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        return Epoch(parsed.timestamp())
    raise TimeFormatError(
        f"time must be an ISO 8601 string or Unix seconds, got {type(value).__name__}",
        location="time",
    )


def iso_from_unix(unix_seconds: float) -> str:
    dt = datetime.fromtimestamp(unix_seconds, tz=timezone.utc)
    return dt.isoformat().replace("+00:00", "Z")


def julian_day(unix_seconds: float) -> float:
    """Julian Date (UT) for a POSIX timestamp."""
    return unix_seconds / 86400.0 + 2440587.5


def gmst_radians(unix_seconds: float) -> float:
    """Greenwich Mean Sidereal Time angle, radians in [0, 2*pi)."""
    jd = julian_day(unix_seconds)
    d = jd - 2451545.0
    centuries = d / 36525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * d
        + 0.000387933 * centuries * centuries
        - centuries**3 / 38_710_000.0
    )
    return math.fmod(math.radians(gmst_deg), 2.0 * math.pi) % (2.0 * math.pi)
