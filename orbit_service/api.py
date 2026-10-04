"""FastAPI HTTP layer."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from .service import OrbitService
from .errors import OrbitServiceError
from .persistence import StoreConfig
from .service import ServiceConfig


def create_app(service: OrbitService | None = None) -> FastAPI:
    service = service or OrbitService(
        ServiceConfig(store=StoreConfig(path="data/scenes.json"))
    )
    app = FastAPI(title="orbit-service", version="0.1.0")
    app.state.service = service

    @app.exception_handler(OrbitServiceError)
    async def handle_domain_error(
        request: Request, exc: OrbitServiceError
    ) -> JSONResponse:
        return JSONResponse(status_code=exc.http_status, content={"error": exc.to_dict()})

    @app.get("/health")
    async def health() -> dict:
        return {"status": "ok", "scenes": service.list_scenes()}

    @app.put("/orbits/{name}")
    async def put_orbit(name: str, payload: dict) -> dict:
        return service.register_orbit(name, payload, overwrite=True)

    @app.post("/orbits/{name}")
    async def post_orbit(name: str, payload: dict) -> dict:
        return service.register_orbit(name, payload, overwrite=False)

    @app.put("/stations/{name}")
    async def put_station(name: str, payload: dict) -> dict:
        return service.register_station(name, payload, overwrite=True)

    @app.post("/stations/{name}")
    async def post_station(name: str, payload: dict) -> dict:
        return service.register_station(name, payload, overwrite=False)

    @app.get("/scenes")
    async def list_scenes() -> dict:
        return service.list_scenes()

    @app.delete("/scenes/{name}")
    async def delete_scene(name: str) -> dict:
        return service.delete_scene(name)

    @app.post("/propagate")
    async def propagate(payload: dict) -> dict:
        return service.propagate(payload)

    @app.post("/visibility")
    async def visibility(payload: dict) -> dict:
        return service.visibility(payload)

    @app.post("/illumination")
    async def illumination(payload: dict) -> dict:
        return service.illumination(payload)

    @app.post("/batch")
    async def batch(payload: dict) -> dict:
        return await service.batch(payload)

    return app
