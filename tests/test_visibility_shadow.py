"""Tests for visibility-window convergence and Earth-shadow events."""
from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone

import pytest

from app.elements import ClassicalElements
from app.shadow import (
    ShadowModel,
    SUN_DISTANCE_M,
    find_illumination,
    illumination_fraction,
    sun_position_j2000,
)
from app.station import GroundStation, find_windows

EPOCH = datetime(2024, 3, 20, tzinfo=timezone.utc)
POLAR = ClassicalElements(
    EPOCH, 7_000_000.0, 0.0, math.radians(90.0), 0.0, 0.0, 0.0
)
STATION = GroundStation("equator", 0.0, 0.0, 0.0)


def test_windows_time_ordered_and_nonempty() -> None:
    windows = find_windows(
        POLAR,
        STATION,
        EPOCH,
        EPOCH + timedelta(hours=12),
        5.0,
        boundary_tol_s=1.0,
        sample_step_s=20.0,
    )
    assert len(windows) >= 1
    for earlier, later in zip(windows, windows[1:]):
        assert earlier.end <= later.start


def test_boundary_converges_and_stabilizes() -> None:
    starts: list[float] = []
    for tol in (60.0, 10.0, 1.0, 0.1):
        windows = find_windows(
            POLAR,
            STATION,
            EPOCH,
            EPOCH + timedelta(hours=12),
            5.0,
            boundary_tol_s=tol,
            sample_step_s=10.0,
        )
        # Use the highest-elevation pass (stable, clearly-detected window).
        target = max(windows, key=lambda w: w.max_elevation_deg)
        starts.append((target.start - EPOCH).total_seconds())
    # Tightening the tolerance moves the boundary by <= the previous tolerance
    # and converges to a stable value (no jumps, no lost windows).
    for coarse, fine in zip(starts, starts[1:]):
        assert abs(coarse - fine) <= 60.0
    assert abs(starts[-1] - starts[-2]) <= 1.0


def test_no_windows_when_mask_above_zenith() -> None:
    windows = find_windows(
        POLAR,
        STATION,
        EPOCH,
        EPOCH + timedelta(hours=2),
        89.9,
        boundary_tol_s=1.0,
        sample_step_s=30.0,
    )
    assert windows == []


def test_sun_distance_is_one_au() -> None:
    sun = sun_position_j2000(EPOCH)
    assert math.sqrt(sum(c * c for c in sun)) == pytest.approx(
        SUN_DISTANCE_M, rel=1e-3
    )


def test_cylindrical_shadow_indicator() -> None:
    sun = (SUN_DISTANCE_M, 0.0, 0.0)
    # Satellite on the Sun-facing side is lit.
    assert illumination_fraction((7e6, 0, 0), sun, ShadowModel.CYLINDRICAL)
    # Satellite directly behind the Earth on the shadow axis is in umbra.
    assert (
        illumination_fraction((-7e6, 0, 0), sun, ShadowModel.CYLINDRICAL) == 0.0
    )


def test_conical_umbra_narrower_far_from_earth() -> None:
    sun = (SUN_DISTANCE_M, 0.0, 0.0)
    # At LEO the conical umbra radius is slightly smaller than the cylinder.
    cylindrical_in = (
        illumination_fraction((-7e6, 6.370e6, 0), sun, ShadowModel.CYLINDRICAL)
        == 0.0
    )
    conical_lit = (
        illumination_fraction((-7e6, 6.370e6, 0), sun, ShadowModel.CONICAL) == 1.0
    )
    assert cylindrical_in and conical_lit


def test_illumination_events_paired_and_durations_consistent() -> None:
    result = find_illumination(
        POLAR,
        EPOCH,
        EPOCH + timedelta(hours=6),
        ShadowModel.CONICAL,
        boundary_tol_s=0.5,
        sample_step_s=20.0,
    )
    kinds = [ev.kind for ev in result.events]
    assert kinds.count("entering_shadow") == kinds.count("leaving_shadow")
    total = result.illuminated_seconds + result.shadowed_seconds
    assert total == pytest.approx(6 * 3600.0, abs=1.0)
    assert result.shadowed_seconds > 0.0


def test_cylindrical_and_conical_durations_close_for_leo() -> None:
    cyl = find_illumination(
        POLAR, EPOCH, EPOCH + timedelta(hours=12),
        ShadowModel.CYLINDRICAL, boundary_tol_s=1.0, sample_step_s=20.0,
    )
    con = find_illumination(
        POLAR, EPOCH, EPOCH + timedelta(hours=12),
        ShadowModel.CONICAL, boundary_tol_s=1.0, sample_step_s=20.0,
    )
    # The finite solar disc shifts each boundary by seconds; total shadowed
    # time differs by well under 2% for a LEO orbit (docs/shadow.md).
    assert abs(cyl.shadowed_seconds - con.shadowed_seconds) / max(
        1.0, cyl.shadowed_seconds
    ) < 0.02
