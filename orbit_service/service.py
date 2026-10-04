"""Application service: validates inputs and orchestrates computations."""

from __future__ import annotations

from dataclasses import dataclass
import math

from .constants import DEFAULT_ROOT_TOL_SECONDS
from .eclipse import IlluminationResult, ShadowModel, compute_illumination
from .persistence import SceneStore, StoreConfig
from .elements import ClassicalElements
from .errors import FieldRangeError, OrbitServiceError, TimeFormatError
from .groundstation import GroundStation, VisibilityWindow, compute_windows
from .propagator import PropagatorConfig, PropagationMode, StateVector, propagate
from .scheduler import Scheduler, SchedulerConfig
from .time_utils import Epoch, iso_from_unix, parse_time
from .tle import parse_tle, tle_to_elements


@dataclass(frozen=True)
class ServiceConfig:
    scheduler: SchedulerConfig = SchedulerConfig()
    store: StoreConfig = StoreConfig()


class OrbitService:
    def __init__(self, config: ServiceConfig | None = None) -> None:
        self.config = config or ServiceConfig()
        self.scheduler = Scheduler(self.config.scheduler)
        self.store = SceneStore(self.config.store)

    # ---- input parsing ---------------------------------------------------

    def _require(self, payload: dict, key: str, location: str | None = None):
        if not isinstance(payload, dict) or key not in payload:
            raise FieldRangeError(
                f"missing required field {key!r}", location=location or key
            )
        return payload[key]

    def parse_orbit(self, payload: dict, *, prefix: str = "orbit") -> ClassicalElements:
        if not isinstance(payload, dict):
            raise FieldRangeError("orbit must be an object", location=prefix)
        kind = payload.get("format") or payload.get("type")
        if kind is None:
            raise FieldRangeError(
                "orbit format must be 'tle' or 'coe'", location=f"{prefix}.format"
            )
        if kind == "tle":
            line1 = self._require(payload, "line1", f"{prefix}.line1")
            line2 = self._require(payload, "line2", f"{prefix}.line2")
            if not isinstance(line1, str) or not isinstance(line2, str):
                raise FieldRangeError(
                    "TLE lines must be strings", location=f"{prefix}.line1"
                )
            return tle_to_elements(parse_tle(line1, line2))
        if kind in ("coe", "classical", "elements"):
            angles = payload.get("angles", "radians")
            if angles not in ("radians", "degrees"):
                raise FieldRangeError(
                    "angles must be 'radians' or 'degrees'",
                    location=f"{prefix}.angles",
                )
            scale = math.pi / 180.0 if angles == "degrees" else 1.0

            def angle(key: str):
                return float(self._require(payload, key, f"{prefix}.{key}")) * scale

            epoch_value = self._require(payload, "epoch", f"{prefix}.epoch")
            elements = ClassicalElements(
                semi_major_axis=float(
                    self._require(payload, "semi_major_axis", f"{prefix}.semi_major_axis")
                ),
                eccentricity=float(
                    self._require(payload, "eccentricity", f"{prefix}.eccentricity")
                ),
                inclination=angle("inclination"),
                raan=angle("raan"),
                argument_of_periapsis=angle("argument_of_periapsis"),
                mean_anomaly=angle("mean_anomaly"),
                epoch=parse_time(epoch_value),
            )
            # Trigger physical validation through the conversion.
            from .elements import coe_to_eci

            coe_to_eci(elements)
            return elements
        raise FieldRangeError(
            f"unknown orbit format {kind!r}", location=f"{prefix}.format"
        )

    def parse_station(self, payload: dict, *, prefix: str = "station") -> GroundStation:
        if not isinstance(payload, dict):
            raise FieldRangeError("station must be an object", location=prefix)
        station = GroundStation(
            name=str(payload.get("name", prefix)),
            latitude_deg=float(self._require(payload, "latitude_deg", f"{prefix}.latitude_deg")),
            longitude_deg=float(self._require(payload, "longitude_deg", f"{prefix}.longitude_deg")),
            altitude_metres=float(payload.get("altitude_metres", 0.0)),
        )
        from .groundstation import validate_station

        validate_station(station)
        return station

    def _resolve_orbit(self, payload: dict, prefix: str = "orbit") -> ClassicalElements:
        if isinstance(payload, str):
            stored = self.store.get(f"orbit:{payload}")
            return self.parse_orbit(stored["payload"], prefix=payload)
        if isinstance(payload, dict) and "ref" in payload:
            stored = self.store.get(f"orbit:{payload['ref']}")
            return self.parse_orbit(stored["payload"], prefix=str(payload["ref"]))
        return self.parse_orbit(payload, prefix=prefix)

    def _resolve_station(self, payload, prefix: str = "station") -> GroundStation:
        if isinstance(payload, str):
            stored = self.store.get(f"station:{payload}")
            return self.parse_station(stored["payload"], prefix=payload)
        if isinstance(payload, dict) and "ref" in payload:
            stored = self.store.get(f"station:{payload['ref']}")
            return self.parse_station(stored["payload"], prefix=str(payload["ref"]))
        return self.parse_station(payload, prefix=prefix)

    def _propagator_config(self, payload: dict) -> PropagatorConfig:
        settings = payload.get("propagator", {}) if isinstance(payload, dict) else {}
        if not isinstance(settings, dict):
            raise FieldRangeError(
                "propagator must be an object", location="propagator"
            )
        mode_value = settings.get("mode", PropagationMode.J2.value)
        try:
            mode = PropagationMode(mode_value)
        except ValueError as exc:
            raise FieldRangeError(
                f"unknown propagation mode {mode_value!r}", location="propagator.mode"
            ) from exc
        return PropagatorConfig(
            mode=mode,
            step_seconds=float(settings.get("step_seconds", 15.0)),
            rtol=float(settings.get("rtol", 1e-10)),
            atol_metres=float(settings.get("atol_metres", 1e-6)),
            max_steps=int(settings.get("max_steps", 200_000)),
        )

    # ---- serialization --------------------------------------------------

    @staticmethod
    def _state_to_dict(state: StateVector) -> dict:
        return {
            "time": iso_from_unix(state.epoch.unix_seconds),
            "unix_seconds": state.epoch.unix_seconds,
            "position_metres": list(state.position),
            "velocity_mps": list(state.velocity),
            "frame": "ECI_J2000_LIKE",
        }

    @staticmethod
    def _window_to_dict(window: VisibilityWindow) -> dict:
        return {
            "start": iso_from_unix(window.start.unix_seconds),
            "start_unix": window.start.unix_seconds,
            "end": iso_from_unix(window.end.unix_seconds),
            "end_unix": window.end.unix_seconds,
            "duration_seconds": window.end.unix_seconds - window.start.unix_seconds,
            "maximum_elevation_deg": window.maximum_elevation_deg,
            "maximum_at": iso_from_unix(window.maximum_at.unix_seconds),
        }

    @staticmethod
    def _illumination_to_dict(result: IlluminationResult) -> dict:
        return {
            "events": [
                {
                    "time": iso_from_unix(event.time.unix_seconds),
                    "unix_seconds": event.time.unix_seconds,
                    "kind": event.kind,
                }
                for event in result.events
            ],
            "illuminated_seconds": result.illuminated_seconds,
            "eclipse_seconds": result.total_seconds - result.illuminated_seconds,
            "total_seconds": result.total_seconds,
        }

    # ---- operations ------------------------------------------------------

    def register_orbit(self, name: str, payload: dict, *, overwrite: bool = False) -> dict:
        self.parse_orbit(payload, prefix=name)
        self.store.put(f"orbit:{name}", {"kind": "orbit", "payload": payload}, overwrite=overwrite)
        return {"registered": name, "type": "orbit"}

    def register_station(self, name: str, payload: dict, *, overwrite: bool = False) -> dict:
        self.parse_station(payload, prefix=name)
        self.store.put(
            f"station:{name}", {"kind": "station", "payload": payload}, overwrite=overwrite
        )
        return {"registered": name, "type": "station"}

    def list_scenes(self) -> dict:
        orbits = [n[6:] for n in self.store.names() if n.startswith("orbit:")]
        stations = [n[8:] for n in self.store.names() if n.startswith("station:")]
        return {"orbits": sorted(orbits), "stations": sorted(stations)}

    def delete_scene(self, name: str) -> dict:
        removed = []
        for prefixed in (f"orbit:{name}", f"station:{name}"):
            try:
                self.store.delete(prefixed)
                removed.append(prefixed.split(":", 1)[0])
            except Exception:
                pass
        if not removed:
            from .errors import SceneNotFoundError

            raise SceneNotFoundError(f"scene {name!r} is not registered", location=name)
        return {"deleted": name, "types": removed}

    def propagate(self, payload: dict) -> dict:
        elements = self._resolve_orbit(self._require(payload, "orbit"), "orbit")
        target = parse_time(self._require(payload, "time", "time"))
        config = self._propagator_config(payload)
        state = propagate(elements, target, config)
        return self._state_to_dict(state)

    def visibility(self, payload: dict) -> dict:
        elements = self._resolve_orbit(self._require(payload, "orbit"), "orbit")
        station = self._resolve_station(self._require(payload, "station"), "station")
        start = parse_time(self._require(payload, "start", "start"))
        end = parse_time(self._require(payload, "end", "end"))
        minimum_elevation = float(payload.get("minimum_elevation_deg", 0.0))
        tolerance = float(
            payload.get("boundary_tolerance_seconds", DEFAULT_ROOT_TOL_SECONDS)
        )
        config = self._propagator_config(payload)
        windows = compute_windows(
            elements,
            station,
            start,
            end,
            minimum_elevation,
            tolerance,
            config,
        )
        return {
            "station": station.name,
            "minimum_elevation_deg": minimum_elevation,
            "window_count": len(windows),
            "windows": [self._window_to_dict(window) for window in windows],
        }

    def illumination(self, payload: dict) -> dict:
        elements = self._resolve_orbit(self._require(payload, "orbit"), "orbit")
        start = parse_time(self._require(payload, "start", "start"))
        end = parse_time(self._require(payload, "end", "end"))
        model_value = payload.get("shadow_model", ShadowModel.CYLINDRICAL.value)
        try:
            model = ShadowModel(model_value)
        except ValueError as exc:
            raise FieldRangeError(
                f"unknown shadow model {model_value!r}", location="shadow_model"
            ) from exc
        tolerance = float(
            payload.get("boundary_tolerance_seconds", DEFAULT_ROOT_TOL_SECONDS)
        )
        config = self._propagator_config(payload)
        result = compute_illumination(
            elements, start, end, model, tolerance, config
        )
        output = self._illumination_to_dict(result)
        output["shadow_model"] = model.value
        return output

    async def batch(self, payload: dict) -> dict:
        queries = self._require(payload, "queries", "queries")
        if not isinstance(queries, list) or not queries:
            raise FieldRangeError("queries must be a non-empty list", location="queries")

        def make_job(index: int, query: dict):
            operation = query.get("operation", "propagate")

            def job():
                if operation == "propagate":
                    return {"index": index, "result": self.propagate(query)}
                if operation == "visibility":
                    return {"index": index, "result": self.visibility(query)}
                if operation == "illumination":
                    return {"index": index, "result": self.illumination(query)}
                raise FieldRangeError(
                    f"unknown operation {operation!r}",
                    location=f"queries[{index}].operation",
                )

            return job

        jobs = [make_job(index, query) for index, query in enumerate(queries)]
        results = await self.scheduler.run_many(jobs)
        return {"results": results}
