"""Application service: orchestrates parsing, propagation, visibility, shadow."""
from __future__ import annotations

from dataclasses import dataclass

from .elements import ClassicalElements
from .errors import ErrorCode, OrbitError
from .models import (
    BatchRequest,
    PropagateRequest,
    VisibilityRequest,
)
from .propagator import Perturbation, PropagateConfig, propagate
from .scheduler import BoundedScheduler, SchedulerConfig
from .shadow import ShadowModel, find_illumination
from .station import GroundStation, find_windows
from .tle import parse_tle, tle_to_elements
from .timeutil import parse_time


@dataclass
class ServiceConfig:
    scheduler: SchedulerConfig | None = None
    storage_path: str = "data/scenarios.json"


class OrbitService:
    def __init__(self, config: ServiceConfig | None = None) -> None:
        self.config = config or ServiceConfig()
        self.scheduler = BoundedScheduler(self.config.scheduler)

    def resolve_orbit(
        self, orbit
    ) -> tuple[ClassicalElements, PropagateConfig]:
        """Validate orbit input; return (classical elements, prop config)."""
        if (orbit.tle is None) == (orbit.elements is None):
            raise OrbitError(
                ErrorCode.VALIDATION_ERROR,
                "exactly one of 'tle' or 'elements' must be provided",
                location=f"orbits[{orbit.orbit_id}]",
            )
        try:
            perturbation = Perturbation(orbit.perturbation)
        except ValueError as exc:
            raise OrbitError(
                ErrorCode.VALIDATION_ERROR,
                f"unknown perturbation {orbit.perturbation!r} "
                "(expected 'none' or 'j2')",
                location=f"orbits[{orbit.orbit_id}].perturbation",
            ) from exc
        if orbit.tle is not None:
            tle = parse_tle(
                orbit.tle.line0 or "",
                orbit.tle.line1,
                orbit.tle.line2,
            )
            el = tle_to_elements(tle)
        else:
            src = orbit.elements
            el = ClassicalElements(
                epoch=parse_time(src.epoch),
                semi_major_axis_m=src.semi_major_axis_m,
                eccentricity=src.eccentricity,
                inclination=src.inclination_rad,
                raan=src.raan_rad,
                arg_perigee=src.arg_perigee_rad,
                mean_anomaly=src.mean_anomaly_rad,
            )
        from .elements import validate_elements

        validate_elements(el)
        return el, PropagateConfig(perturbation=perturbation)

    def propagate(self, request: PropagateRequest) -> dict:
        el, config = self.resolve_orbit(request.orbit)
        points = []
        for value in request.times:
            when = parse_time(value)
            state = propagate(el, when, config)
            points.append(
                {
                    "time": when.isoformat().replace("+00:00", "Z"),
                    "unix": when.timestamp(),
                    "position_j2000_m": list(state.position),
                    "velocity_j2000_m_s": list(state.velocity),
                }
            )
        return {"orbit_id": request.orbit.orbit_id, "points": points}

    def visibility(self, request: VisibilityRequest) -> dict:
        el, config = self.resolve_orbit(request.orbit)
        station = GroundStation(
            name=request.station.name,
            latitude_deg=request.station.latitude_deg,
            longitude_deg=request.station.longitude_deg,
            altitude_m=request.station.altitude_m,
        )
        start = parse_time(request.start)
        end = parse_time(request.end)
        windows = find_windows(
            el,
            station,
            start,
            end,
            request.min_elevation_deg,
            boundary_tol_s=request.boundary_tol_s,
            config=config,
            sample_step_s=request.sample_step_s,
        )
        return {
            "orbit_id": request.orbit.orbit_id,
            "station": station.name,
            "windows": [
                {
                    "start": w.start.isoformat().replace("+00:00", "Z"),
                    "end": w.end.isoformat().replace("+00:00", "Z"),
                    "max_elevation_deg": w.max_elevation_deg,
                }
                for w in windows
            ],
        }

    def illumination(self, request) -> dict:
        el, config = self.resolve_orbit(request.orbit)
        try:
            model = ShadowModel(request.model)
        except ValueError as exc:
            raise OrbitError(
                ErrorCode.VALIDATION_ERROR,
                f"unknown shadow model {request.model!r}",
                location="model",
            ) from exc
        result = find_illumination(
            el,
            parse_time(request.start),
            parse_time(request.end),
            model=model,
            boundary_tol_s=request.boundary_tol_s,
            sample_step_s=request.sample_step_s,
            config=config,
        )
        return {
            "orbit_id": request.orbit.orbit_id,
            "model": model.value,
            "events": [
                {
                    "time": ev.time.isoformat().replace("+00:00", "Z"),
                    "kind": ev.kind,
                }
                for ev in result.events
            ],
            "illuminated_seconds": result.illuminated_seconds,
            "shadowed_seconds": result.shadowed_seconds,
        }

    async def batch(self, request: BatchRequest) -> dict:
        """Run propagation for all orbits under bounded concurrency.

        Visibility is computed for every orbit x station pair when stations
        are supplied.  Results are returned in the deterministic input order.
        """
        start = parse_time(request.start) if request.start is not None else None
        end = parse_time(request.end) if request.end is not None else None

        jobs = []
        for orbit in request.orbits:
            for station in request.stations:
                jobs.append(("visibility", orbit, station))
            if start is not None and end is not None:
                jobs.append(("propagate", orbit))

        async def run_job(job):
            kind = job[0]
            if kind == "propagate":
                orbit = job[1]
                el, config = self.resolve_orbit(orbit)
                points = []
                for value in (start, end):
                    state = propagate(el, value, config)
                    points.append(
                        {
                            "time": value.isoformat().replace("+00:00", "Z"),
                            "unix": value.timestamp(),
                            "position_j2000_m": list(state.position),
                            "velocity_j2000_m_s": list(state.velocity),
                        }
                    )
                return {"kind": "propagate", "orbit_id": orbit.orbit_id,
                        "points": points}
            orbit = job[1]
            station = job[2]
            return self.visibility(
                VisibilityRequest(
                    orbit=orbit,
                    station=station,
                    start=request.start,
                    end=request.end,
                    min_elevation_deg=request.min_elevation_deg,
                    boundary_tol_s=request.boundary_tol_s,
                    sample_step_s=request.sample_step_s,
                )
            )

        results = await self.scheduler.map_ordered(run_job, jobs)
        return {"results": results}
