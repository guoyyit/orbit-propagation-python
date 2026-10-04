"""Tests for bounded concurrency, persistence, determinism and HTTP errors."""
from __future__ import annotations

import asyncio
import os

import pytest
from fastapi.testclient import TestClient

from app.errors import ErrorCode, OrbitError
from app.main import create_app
from app.models import ElementsInput, OrbitInput, PropagateRequest
from app.scheduler import BoundedScheduler, SchedulerConfig
from app.service import OrbitService, ServiceConfig
from app.storage import ScenarioStore, StoredScenario


ELEMENTS_BODY = {
    "epoch": "2024-03-20T00:00:00Z",
    "semi_major_axis_m": 7_000_000.0,
    "eccentricity": 0.001,
    "inclination_rad": 0.9013,
    "raan_rad": 4.318,
    "arg_perigee_rad": 0.698,
    "mean_anomaly_rad": 0.1745,
}


def test_concurrency_cap_respected() -> None:
    active = 0
    peak = 0

    async def job(_: int) -> int:
        nonlocal active, peak
        active += 1
        peak = max(peak, active)
        await asyncio.sleep(0.02)
        active -= 1
        return _

    scheduler = BoundedScheduler(
        SchedulerConfig(max_concurrency=2, max_queue_depth=8)
    )
    results = asyncio.run(scheduler.map_ordered(job, list(range(10))))
    assert results == list(range(10))
    assert peak == 2


def test_queue_full_raises_distinct_code() -> None:
    scheduler = BoundedScheduler(
        SchedulerConfig(max_concurrency=1, max_queue_depth=1, wait_timeout_s=5)
    )
    gate = asyncio.Event()

    async def slow() -> None:
        await gate.wait()

    async def scenario() -> None:
        t0 = asyncio.create_task(scheduler.run(slow))
        t1 = asyncio.create_task(scheduler.run(slow))
        await asyncio.sleep(0.01)
        with pytest.raises(OrbitError) as exc:
            await scheduler.run(slow)
        assert exc.value.code is ErrorCode.QUEUE_FULL
        gate.set()
        await asyncio.gather(t0, t1)

    asyncio.run(scenario())


def test_queue_timeout_raises_distinct_code() -> None:
    scheduler = BoundedScheduler(
        SchedulerConfig(max_concurrency=1, max_queue_depth=4, wait_timeout_s=0.02)
    )
    gate = asyncio.Event()

    async def hold() -> None:
        await gate.wait()

    async def scenario() -> None:
        holder = asyncio.create_task(scheduler.run(hold))
        await asyncio.sleep(0.01)
        with pytest.raises(OrbitError) as exc:
            await scheduler.run(lambda: asyncio.sleep(0))
        assert exc.value.code is ErrorCode.QUEUE_TIMEOUT
        gate.set()
        await holder

    asyncio.run(scenario())


def test_persistence_roundtrip_and_restart(tmp_path) -> None:
    path = tmp_path / "scenarios.json"
    store = ScenarioStore(path)
    store.save(StoredScenario("sc1", {"a": 1}, 123.0))

    reopened = ScenarioStore(path)
    scenario = reopened.get("sc1")
    assert scenario.payload == {"a": 1}
    assert scenario.created_at_unix == 123.0


def test_computation_deterministic_before_and_after_restart(tmp_path) -> None:
    path = tmp_path / "scenarios.json"
    svc1 = OrbitService(ServiceConfig(storage_path=str(path)))
    svc2 = OrbitService(ServiceConfig(storage_path=str(path)))

    request = PropagateRequest(
        orbit=OrbitInput(
            orbit_id="o", elements=ElementsInput(**ELEMENTS_BODY)
        ),
        times=["2024-03-20T01:00:00Z", 1710950400],
    )
    first = svc1.propagate(request)
    second = svc2.propagate(request)
    assert first == second


def test_batch_is_order_independent(tmp_path) -> None:
    async def run() -> list:
        from app.models import BatchRequest, StationInput

        svc = OrbitService()
        orbits = [
            OrbitInput(orbit_id=f"o{k}", elements=ElementsInput(**ELEMENTS_BODY))
            for k in range(8)
        ]
        stations = [
            StationInput(name=f"s{k}", latitude_deg=k * 5.0, longitude_deg=k * 7.0)
            for k in range(3)
        ]
        request = BatchRequest(
            orbits=orbits,
            stations=stations,
            start="2024-03-20T00:00:00Z",
            end="2024-03-20T02:00:00Z",
            sample_step_s=60.0,
        )
        return [await svc.batch(request) for _ in range(3)]

    runs = asyncio.run(run())
    assert runs[0] == runs[1] == runs[2]


def test_api_propagate_and_error_codes() -> None:
    client = TestClient(create_app())
    response = client.post(
        "/api/v1/propagate",
        json={"orbit": {"orbit_id": "o", "elements": ELEMENTS_BODY},
              "times": ["2024-03-20T00:30:00Z"]},
    )
    assert response.status_code == 200
    point = response.json()["points"][0]
    assert len(point["position_j2000_m"]) == 3

    bad = dict(ELEMENTS_BODY)
    bad["eccentricity"] = 1.5
    response = client.post(
        "/api/v1/propagate",
        json={"orbit": {"orbit_id": "bad", "elements": bad}, "times": [0]},
    )
    assert response.status_code == 400
    assert response.json()["code"] == "FIELD_OUT_OF_RANGE"


def test_api_iso_and_unix_times_agree() -> None:
    client = TestClient(create_app())
    iso = client.post(
        "/api/v1/propagate",
        json={"orbit": {"orbit_id": "o", "elements": ELEMENTS_BODY},
              "times": ["2024-03-20T00:30:00Z"]},
    ).json()
    unix = client.post(
        "/api/v1/propagate",
        json={"orbit": {"orbit_id": "o", "elements": ELEMENTS_BODY},
              "times": [1710894600]},
    ).json()
    assert iso["points"][0]["position_j2000_m"] == pytest.approx(
        unix["points"][0]["position_j2000_m"], abs=1e-6
    )


def test_api_requires_exactly_one_orbit_format() -> None:
    client = TestClient(create_app())
    body = {"orbit": {"orbit_id": "o"}, "times": [0]}
    response = client.post("/api/v1/propagate", json=body)
    assert response.status_code == 400
