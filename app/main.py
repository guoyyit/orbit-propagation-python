"""FastAPI application: HTTP surface for orbit propagation service."""
from __future__ import annotations

from datetime import datetime, timezone

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from .errors import ErrorCode, OrbitError
from .models import (
    BatchRequest,
    IlluminationRequest,
    PropagateRequest,
    ScenarioRequest,
    VisibilityRequest,
)
from .service import OrbitService
from .storage import ScenarioStore, StoredScenario
from .timeutil import to_unix


def create_app(service: OrbitService | None = None) -> FastAPI:
    app = FastAPI(title="Local Orbit Propagation Service", version="0.1.0")
    svc = service or OrbitService()
    app.state.service = svc
    app.state.store = ScenarioStore(svc.config.storage_path)
    store = app.state.store

    @app.exception_handler(OrbitError)
    async def orbit_error_handler(_: Request, exc: OrbitError) -> JSONResponse:
        status = 400
        if exc.code in (ErrorCode.QUEUE_FULL, ErrorCode.QUEUE_TIMEOUT):
            status = 503
        elif exc.code is ErrorCode.NOT_FOUND:
            status = 404
        return JSONResponse(
            status_code=status,
            content={
                "code": exc.code.value,
                "message": exc.message,
                "location": exc.location,
            },
        )

    @app.get("/health")
    async def health() -> dict:
        return {"status": "ok"}

    @app.post("/api/v1/propagate")
    async def propagate(request: PropagateRequest) -> dict:
        return svc.propagate(request)

    @app.post("/api/v1/visibility")
    async def visibility(request: VisibilityRequest) -> dict:
        return svc.visibility(request)

    @app.post("/api/v1/illumination")
    async def illumination(request: IlluminationRequest) -> dict:
        return svc.illumination(request)

    @app.post("/api/v1/batch")
    async def batch(request: BatchRequest) -> dict:
        return await svc.batch(request)

    @app.post("/api/v1/scenarios")
    async def save_scenario(request: ScenarioRequest) -> dict:
        if not request.scenario_id or not request.scenario_id.strip():
            raise OrbitError(
                ErrorCode.VALIDATION_ERROR,
                "scenario_id must be a non-empty string",
                location="scenario_id",
            )
        store.save(
            StoredScenario(
                scenario_id=request.scenario_id,
                payload=request.request,
                created_at_unix=to_unix(datetime.now(timezone.utc)),
            )
        )
        return {"scenario_id": request.scenario_id, "saved": True}

    @app.get("/api/v1/scenarios/{scenario_id}")
    async def get_scenario(scenario_id: str) -> dict:
        scenario = store.get(scenario_id)
        return {
            "scenario_id": scenario.scenario_id,
            "created_at_unix": scenario.created_at_unix,
            "request": scenario.payload,
        }

    return app


app = create_app()
