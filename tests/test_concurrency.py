"""Bounded concurrency, queue limits and deterministic batch results."""

from __future__ import annotations

import asyncio
import threading
import time

import pytest

from orbit_service.errors import QueueOverflowError, WaitTimeoutError
from orbit_service.scheduler import Scheduler, SchedulerConfig
from orbit_service.service import OrbitService, ServiceConfig
from orbit_service.persistence import StoreConfig


async def _run(cfg, count):
    active = 0
    peak = 0
    lock = threading.Lock()

    def job(index):
        nonlocal active, peak
        with lock:
            active += 1
            peak = max(peak, active)
        time.sleep(0.02)
        with lock:
            active -= 1
        return index * index

    scheduler = Scheduler(cfg)
    values = await scheduler.run_many([lambda i=i: job(i) for i in range(count)])
    return values, peak


def test_concurrency_cap_and_order():
    values, peak = asyncio.run(
        _run(SchedulerConfig(max_concurrency=2, max_queue_depth=4), 6)
    )
    assert values == [i * i for i in range(6)]
    assert peak == 2


def test_queue_overflow_is_distinguishable():
    with pytest.raises(QueueOverflowError):
        asyncio.run(
            _run(SchedulerConfig(max_concurrency=1, max_queue_depth=1), 8)
        )


def test_wait_timeout_is_distinguishable():
    with pytest.raises(WaitTimeoutError):
        asyncio.run(
            _run(
                SchedulerConfig(
                    max_concurrency=1,
                    max_queue_depth=10,
                    wait_timeout_seconds=0.01,
                ),
                6,
            )
        )


def test_batch_results_independent_of_scheduling(tmp_path):
    store = StoreConfig(str(tmp_path / "scenes.json"))
    service = OrbitService(
        ServiceConfig(
            scheduler=SchedulerConfig(max_concurrency=1, max_queue_depth=20),
            store=store,
        )
    )
    service.register_orbit(
        "iss",
        {
            "format": "coe",
            "angles": "degrees",
            "semi_major_axis": 6_794_900.0,
            "eccentricity": 0.0015,
            "inclination": 51.64,
            "raan": 120.0,
            "argument_of_periapsis": 60.0,
            "mean_anomaly": 0.0,
            "epoch": 1_700_000_000.0,
        },
    )
    queries = [
        {
            "operation": "propagate",
            "orbit": "iss",
            "time": 1_700_000_000 + 1000 * i,
        }
        for i in range(12)
    ]
    first = asyncio.run(service.batch({"queries": queries}))
    second = asyncio.run(service.batch({"queries": list(reversed(queries))}))
    by_time_first = {
        queries[item["index"]]["time"]: item["result"]["position_metres"]
        for item in first["results"]
    }
    reversed_queries = list(reversed(queries))
    by_time_second = {
        reversed_queries[item["index"]]["time"]: item["result"]["position_metres"]
        for item in second["results"]
    }
    assert by_time_first == by_time_second
    assert [item["index"] for item in first["results"]] == list(range(12))
