"""Pydantic request/response models for the HTTP API."""
from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class TleInput(BaseModel):
    line0: str | None = None
    line1: str
    line2: str


class ElementsInput(BaseModel):
    epoch: str | int | float
    semi_major_axis_m: float
    eccentricity: float
    inclination_rad: float
    raan_rad: float = 0.0
    arg_perigee_rad: float = 0.0
    mean_anomaly_rad: float = 0.0


class OrbitInput(BaseModel):
    orbit_id: str
    tle: TleInput | None = None
    elements: ElementsInput | None = None
    perturbation: str = "j2"


class StationInput(BaseModel):
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float = 0.0


class TimeRange(BaseModel):
    start: str | int | float
    end: str | int | float


class PropagateRequest(BaseModel):
    orbit: OrbitInput
    times: list[str | int | float]


class VisibilityRequest(BaseModel):
    orbit: OrbitInput
    station: StationInput
    start: str | int | float
    end: str | int | float
    min_elevation_deg: float = 5.0
    boundary_tol_s: float = 1.0
    sample_step_s: float = 60.0


class IlluminationRequest(BaseModel):
    orbit: OrbitInput
    start: str | int | float
    end: str | int | float
    model: str = "conical"
    boundary_tol_s: float = 1.0
    sample_step_s: float = 60.0


class BatchRequest(BaseModel):
    orbits: list[OrbitInput]
    stations: list[StationInput] = Field(default_factory=list)
    start: str | int | float | None = None
    end: str | int | float | None = None
    min_elevation_deg: float = 5.0
    boundary_tol_s: float = 1.0
    sample_step_s: float = 60.0
    shadow_model: str = "conical"


class ErrorBody(BaseModel):
    code: str
    message: str
    location: str | None = None


class BatchItemError(BaseModel):
    index: int
    code: str
    message: str
    location: str | None = None


class ScenarioRequest(BaseModel):
    scenario_id: str
    request: dict[str, Any]
