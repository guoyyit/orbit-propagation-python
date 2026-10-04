"""Deterministic, atomic JSON persistence for registered scenes.

Scenes are plain JSON-serializable dicts. Writes go to a temporary file and
are atomically renamed, so a crash never leaves a partially written store.
Registered scenes survive process restarts and recompute identically because
the store only holds inputs, never computed state.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass
import json
import os

from .errors import SceneConflictError, SceneNotFoundError


@dataclass(frozen=True)
class StoreConfig:
    path: str = "data/scenes.json"


class SceneStore:
    def __init__(self, config: StoreConfig | None = None) -> None:
        self.config = config or StoreConfig()
        self._lock = threading.RLock()
        directory = os.path.dirname(self.config.path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        self._data: dict[str, dict] = {}
        if os.path.exists(self.config.path):
            with open(self.config.path, "r", encoding="utf-8") as handle:
                loaded = json.load(handle)
            if not isinstance(loaded, dict):
                raise ValueError("scene store root must be a JSON object")
            self._data = loaded

    def _flush(self) -> None:
        directory = os.path.dirname(self.config.path) or "."
        temporary = os.path.join(directory, f".{os.path.basename(self.config.path)}.tmp")
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump(self._data, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
        os.replace(temporary, self.config.path)

    def put(self, name: str, payload: dict, *, overwrite: bool = False) -> None:
        with self._lock:
            if not overwrite and name in self._data:
                raise SceneConflictError(
                    f"scene {name!r} already registered", location=name
                )
            self._data[name] = payload
            self._flush()

    def get(self, name: str) -> dict:
        with self._lock:
            if name not in self._data:
                raise SceneNotFoundError(
                    f"scene {name!r} is not registered", location=name
                )
            return json.loads(json.dumps(self._data[name]))

    def delete(self, name: str) -> None:
        with self._lock:
            if name not in self._data:
                raise SceneNotFoundError(
                    f"scene {name!r} is not registered", location=name
                )
            del self._data[name]
            self._flush()

    def names(self) -> list[str]:
        with self._lock:
            return sorted(self._data.keys())
