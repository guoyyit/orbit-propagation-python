"""Bounded-concurrency scheduler with a finite queue and wait timeout.

Jobs are plain zero-argument callables executed in worker threads. At most
``max_concurrency`` jobs run at once; up to ``max_queue_depth`` additional jobs
wait. Submitting beyond the queue raises ``QueueOverflowError`` immediately; a
job whose turn does not start within ``wait_timeout_seconds`` raises
``WaitTimeoutError``. Results are always returned in the original submission
order, independent of completion order, so scheduling cannot alter results.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import threading

from .errors import QueueOverflowError, WaitTimeoutError


@dataclass(frozen=True)
class SchedulerConfig:
    max_concurrency: int = 4
    max_queue_depth: int = 32
    wait_timeout_seconds: float = 30.0


class Scheduler:
    def __init__(self, config: SchedulerConfig | None = None) -> None:
        self.config = config or SchedulerConfig()
        if self.config.max_concurrency < 1:
            raise ValueError("max_concurrency must be >= 1")
        if self.config.max_queue_depth < 0:
            raise ValueError("max_queue_depth must be >= 0")
        if self.config.wait_timeout_seconds <= 0.0:
            raise ValueError("wait_timeout_seconds must be positive")
        self._slots = threading.BoundedSemaphore(self.config.max_concurrency)
        self._active = 0
        self._active_lock = threading.Lock()

    def _worker(self, job, box: dict) -> None:
        if not self._slots.acquire(timeout=self.config.wait_timeout_seconds):
            box["error"] = WaitTimeoutError(
                f"job waited longer than {self.config.wait_timeout_seconds}s "
                "for a worker slot"
            )
            return
        try:
            box["value"] = job()
        except BaseException as exc:  # captured, re-raised on the caller side
            box["error"] = exc
        finally:
            self._slots.release()

    async def run_many(self, jobs: list[Callable[[], Any]]) -> list[Any]:
        """Execute jobs with bounded concurrency; results preserve input order."""
        import asyncio

        loop = asyncio.get_running_loop()
        in_flight_limit = self.config.max_concurrency + self.config.max_queue_depth
        in_flight = threading.BoundedSemaphore(in_flight_limit)
        boxes: list[dict] = [None] * len(jobs)
        futures: list = []
        try:
            for index, job in enumerate(jobs):
                if not in_flight.acquire(blocking=False):
                    raise QueueOverflowError(
                        f"queue depth {self.config.max_queue_depth} exceeded"
                    )
                box: dict = {}

                def work(job=job, box=box, gate=in_flight) -> None:
                    try:
                        self._worker(job, box)
                    finally:
                        gate.release()

                boxes[index] = box
                futures.append(loop.run_in_executor(None, work))
        except BaseException:
            for future in futures:
                future.cancel()
            raise

        await asyncio.gather(*futures)
        results: list[Any] = []
        for box in boxes:
            if "error" in box:
                raise box["error"]
            results.append(box["value"])
        return results
