"""Distinguishable error types for the orbit service.

Every error exposed through the HTTP API carries a stable ``code`` string and,
when applicable, a JSON pointer style ``location`` describing which field was
rejected.
"""

from __future__ import annotations


class OrbitServiceError(Exception):
    """Base class for all domain errors."""

    code = "internal_error"
    http_status = 500

    def __init__(self, message: str = "", *, location: str | None = None) -> None:
        super().__init__(message)
        self.message = message or self.code
        self.location = location

    def to_dict(self) -> dict:
        payload: dict = {"code": self.code, "message": self.message}
        if self.location is not None:
            payload["location"] = self.location
        return payload


class TLEFormatError(OrbitServiceError):
    code = "tle_format_error"
    http_status = 422


class TLEChecksumError(OrbitServiceError):
    code = "tle_checksum_error"
    http_status = 422


class FieldRangeError(OrbitServiceError):
    code = "field_range_error"
    http_status = 422


class TimeFormatError(OrbitServiceError):
    code = "time_format_error"
    http_status = 422


class NonPhysicalOrbitError(OrbitServiceError):
    code = "non_physical_orbit"
    http_status = 422


class ConvergenceError(OrbitServiceError):
    code = "convergence_error"
    http_status = 422


class QueueOverflowError(OrbitServiceError):
    code = "queue_overflow"
    http_status = 503


class WaitTimeoutError(OrbitServiceError):
    code = "wait_timeout"
    http_status = 504


class SceneNotFoundError(OrbitServiceError):
    code = "scene_not_found"
    http_status = 404


class SceneConflictError(OrbitServiceError):
    code = "scene_conflict"
    http_status = 409
