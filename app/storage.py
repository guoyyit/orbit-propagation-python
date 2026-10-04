"""Persistent registry of scenarios (JSON on local disk).

Writes are atomic (write a temp file, fsync, os.replace) so a crash never
leaves a half-written registry.  JSON keys are sorted and indentation fixed,
making the persisted bytes deterministic for identical inputs.  All data lives
locally; no external service is contacted.
"""
from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .errors import ErrorCode, OrbitError


@dataclass(frozen=True)
class StoredScenario:
    scenario_id: str
    payload: dict
    created_at_unix: float

    def to_dict(self) -> dict:
        return {
            "scenario_id": self.scenario_id,
            "created_at_unix": self.created_at_unix,
            "payload": self.payload,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "StoredScenario":
        return cls(
            scenario_id=data["scenario_id"],
            payload=data["payload"],
            created_at_unix=float(data["created_at_unix"]),
        )


class ScenarioStore:
    """File-backed, process-safe scenario storage (atomic writes, fsync)."""

    def __init__(self, path: str | Path = "data/scenarios.json") -> None:
        self.path = Path(path)

    def save(self, scenario: StoredScenario) -> None:
        """Insert or update one scenario; persisted (fsynced) before return."""
        data = self._load_all()
        data[scenario.scenario_id] = scenario.to_dict()
        self._write_all(data)

    def get(self, scenario_id: str) -> StoredScenario:
        """Fetch one scenario or raise NOT_FOUND."""
        data = self._load_all()
        if scenario_id not in data:
            raise OrbitError(
                ErrorCode.NOT_FOUND,
                f"scenario {scenario_id!r} not found",
                location="scenario_id",
            )
        return StoredScenario.from_dict(data[scenario_id])

    def list_ids(self) -> list[str]:
        """Return all stored scenario ids (sorted, deterministic)."""
        return sorted(self._load_all().keys())

    def delete(self, scenario_id: str) -> None:
        """Delete one scenario or raise NOT_FOUND."""
        data = self._load_all()
        if scenario_id not in data:
            raise OrbitError(
                ErrorCode.NOT_FOUND,
                f"scenario {scenario_id!r} not found",
                location="scenario_id",
            )
        del data[scenario_id]
        self._write_all(data)

    def _load_all(self) -> dict:
        if not self.path.exists():
            return {}
        try:
            with self.path.open("r", encoding="utf-8") as handle:
                data = json.load(handle)
        except (json.JSONDecodeError, OSError) as exc:
            raise OrbitError(
                ErrorCode.STORAGE_ERROR,
                f"failed to read scenario store {self.path}: {exc}",
            ) from exc
        if not isinstance(data, dict):
            raise OrbitError(
                ErrorCode.STORAGE_ERROR,
                f"scenario store {self.path} is not a JSON object",
            )
        return data

    def _write_all(self, data: dict) -> None:
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp_name = tempfile.mkstemp(
                dir=str(self.path.parent), suffix=".tmp", prefix=".scenarios-"
            )
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as handle:
                    json.dump(data, handle, sort_keys=True, indent=2,
                              separators=(",", ": "))
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(tmp_name, self.path)
                # Best-effort directory fsync for crash durability.
                try:
                    dir_fd = os.open(str(self.path.parent), os.O_RDONLY)
                    try:
                        os.fsync(dir_fd)
                    finally:
                        os.close(dir_fd)
                except OSError:
                    pass
            except BaseException:
                if os.path.exists(tmp_name):
                    os.unlink(tmp_name)
                raise
        except OSError as exc:
            raise OrbitError(
                ErrorCode.STORAGE_ERROR,
                f"failed to write scenario store {self.path}: {exc}",
            ) from exc
