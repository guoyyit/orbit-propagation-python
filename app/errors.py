"""Domain error types with distinguishable error codes."""
from __future__ import annotations

import enum
from dataclasses import dataclass


class ErrorCode(str, enum.Enum):
    """Stable, machine-readable error codes."""

    INTERNAL_ERROR = "INTERNAL_ERROR"
    VALIDATION_ERROR = "VALIDATION_ERROR"
    TLE_FORMAT_ERROR = "TLE_FORMAT_ERROR"
    TLE_CHECKSUM_ERROR = "TLE_CHECKSUM_ERROR"
    TLE_FIELD_ERROR = "TLE_FIELD_ERROR"
    ELEMENTS_FORMAT_ERROR = "ELEMENTS_FORMAT_ERROR"
    FIELD_OUT_OF_RANGE = "FIELD_OUT_OF_RANGE"
    NON_PHYSICAL_ORBIT = "NON_PHYSICAL_ORBIT"
    TIME_FORMAT_ERROR = "TIME_FORMAT_ERROR"
    ITERATION_NOT_CONVERGED = "ITERATION_NOT_CONVERGED"
    QUEUE_FULL = "QUEUE_FULL"
    QUEUE_TIMEOUT = "QUEUE_TIMEOUT"
    NOT_FOUND = "NOT_FOUND"
    STORAGE_ERROR = "STORAGE_ERROR"


@dataclass
class OrbitError(Exception):
    """Error carrying a stable code and a JSON-Pointer-like location."""

    code: ErrorCode
    message: str
    location: str | None = None

    def __str__(self) -> str:
        prefix = f"[{self.code.value}]"
        if self.location:
            prefix += f" {self.location}:"
        return f"{prefix} {self.message}"
