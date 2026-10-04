"""Restart persistence and HTTP error contracts."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from orbit_service.api import create_app
from orbit_service.persistence import StoreConfig
from orbit_service.service import OrbitService, ServiceConfig


ORBIT_PAYLOAD = {
    "format": "coe",
    "semi_major_axis": 6_794_900.0,
    "eccentricity": 0.0015,
    "inclination": 0.901,
    "raan": 2.09,
    "argument_of_periapsis": 1.04,
    "mean_anomaly": 0.0,
    "epoch": 1_700_000_000.0,
}


def test_scene_survives_restart_and_recomputes_identically(tmp_path):
    path = str(tmp_path / "scenes.json")
    service = OrbitService(ServiceConfig(store=StoreConfig(path)))
    service.register_orbit("iss", ORBIT_PAYLOAD)
    before = service.propagate({"orbit": "iss", "time": 1_700_050_000})

    restarted = OrbitService(ServiceConfig(store=StoreConfig(path)))
    assert "iss" in restarted.list_scenes()["orbits"]
    after = restarted.propagate({"orbit": "iss", "time": 1_700_050_000})
    assert before == after


@pytest.fixture
def client(tmp_path):
    service = OrbitService(
        ServiceConfig(store=StoreConfig(str(tmp_path / "scenes.json")))
    )
    return TestClient(create_app(service))


def test_http_error_codes_are_distinguishable(client):
    # Non-physical eccentricity
    response = client.post(
        "/propagate",
        json={
            "orbit": {**ORBIT_PAYLOAD, "eccentricity": 1.5},
            "time": 1_700_000_000,
        },
    )
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "non_physical_orbit"

    # Missing field
    response = client.post("/propagate", json={"time": 1_700_000_000})
    assert response.json()["error"]["code"] == "field_range_error"

    # Unknown scene
    response = client.post(
        "/propagate", json={"orbit": "ghost", "time": 1_700_000_000}
    )
    assert response.status_code == 404
    assert response.json()["error"]["code"] == "scene_not_found"


def test_http_batch_and_visibility(client):
    client.put("/orbits/iss", json=ORBIT_PAYLOAD)
    client.put(
        "/stations/bj",
        json={"name": "beijing", "latitude_deg": 39.9, "longitude_deg": 116.4},
    )
    response = client.post(
        "/visibility",
        json={
            "orbit": "iss",
            "station": "bj",
            "start": 1_700_000_000,
            "end": 1_700_086_400,
            "minimum_elevation_deg": 10.0,
        },
    )
    assert response.status_code == 200
    assert response.json()["window_count"] >= 1

    response = client.post(
        "/batch",
        json={
            "queries": [
                {
                    "operation": "propagate",
                    "orbit": "iss",
                    "time": 1_700_000_000 + offset,
                }
                for offset in (0, 1000, 5000)
            ]
        },
    )
    assert response.status_code == 200
    assert [item["index"] for item in response.json()["results"]] == [0, 1, 2]
