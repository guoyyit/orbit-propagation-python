"""Visibility window correctness and boundary convergence."""

from __future__ import annotations

from orbit_service.groundstation import GroundStation, compute_windows, elevation_deg
from orbit_service.propagator import PropagatorConfig, propagate
from orbit_service.time_utils import Epoch


def test_window_boundaries_converge_without_jumps(iss_elements):
    station = GroundStation("beijing", 39.9, 116.4, 50.0)
    start = Epoch(iss_elements.epoch.unix_seconds)
    end = Epoch(iss_elements.epoch.unix_seconds + 86400)
    tolerances = [1.0, 0.1, 0.01, 0.001]
    reference = compute_windows(
        iss_elements, station, start, end, 10.0, tolerances[0]
    )
    previous = [(w.start.unix_seconds, w.end.unix_seconds) for w in reference]
    for tolerance in tolerances[1:]:
        windows = compute_windows(
            iss_elements, station, start, end, 10.0, tolerance
        )
        assert len(windows) == len(reference)
        for index, window in enumerate(windows):
            start_move = abs(window.start.unix_seconds - previous[index][0])
            end_move = abs(window.end.unix_seconds - previous[index][1])
            # Boundaries are found on a piecewise-linear interpolation of the
            # sampled elevation; tightening the bisection tolerance can move a
            # boundary by at most the previous tolerance plus interpolation
            # residual, bounded by a fraction of the coarser tolerance.
            assert start_move <= tolerances[tolerances.index(tolerance) - 1]
            assert end_move <= tolerances[tolerances.index(tolerance) - 1]
        previous = [(w.start.unix_seconds, w.end.unix_seconds) for w in windows]


def test_windows_are_sorted_and_non_overlapping(iss_elements):
    station = GroundStation("beijing", 39.9, 116.4, 50.0)
    start = Epoch(iss_elements.epoch.unix_seconds)
    end = Epoch(iss_elements.epoch.unix_seconds + 2 * 86400)
    windows = compute_windows(iss_elements, station, start, end, 5.0, 1e-3)
    for window in windows:
        assert window.end.unix_seconds > window.start.unix_seconds
        assert window.maximum_elevation_deg >= 5.0 - 1e-9
    for earlier, later in zip(windows, windows[1:]):
        assert earlier.end.unix_seconds <= later.start.unix_seconds


def test_reported_boundaries_bracket_mask_crossing(iss_elements):
    station = GroundStation("beijing", 39.9, 116.4, 50.0)
    start = Epoch(iss_elements.epoch.unix_seconds)
    end = Epoch(iss_elements.epoch.unix_seconds + 86400)
    config = PropagatorConfig()
    windows = compute_windows(
        iss_elements, station, start, end, 10.0, 1e-3, config
    )
    assert windows
    first = windows[0]
    before = elevation_deg(
        station,
        propagate(iss_elements, Epoch(first.start.unix_seconds - 5.0), config),
    )
    after = elevation_deg(
        station,
        propagate(iss_elements, Epoch(first.start.unix_seconds + 5.0), config),
    )
    assert before < 10.0 < after
