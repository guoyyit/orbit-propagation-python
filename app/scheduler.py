"""Bounded-concurrency scheduler with a finite FIFO queue and wait timeout."""
from __future__ import annotations

import asyncio
from dataclasses import dataclass
from typing import Any, Awaitable, Callable

from .errors import ErrorCode, OrbitError


@dataclass
class SchedulerConfig:
    max_concurrency: int = 4
    max_queue_depth: int = 128
    wait_timeout_s: float = 30.0

    def __post_init__(self) -> None:
        if self.max_concurrency < 1:
            raise OrbitError(
                ErrorCode.FIELD_OUT_OF_RANGE,
                "max_concurrency must be >= 1",
                location="max_concurrency",
            )
        if self.max_queue_depth < 0:
            raise OrbitError(
                ErrorCode.FIELD_OUT_OF_RANGE,
                "max_queue_depth must be >= 0",
                location="max_queue_depth",
            )
        if self.wait_timeout_s <= 0:
            raise OrbitError(
                ErrorCode.FIELD_OUT_OF_RANGE,
                "wait_timeout_s must be > 0",
                location="wait_timeout_s",
            )


class BoundedScheduler:
    """Runs coroutines under a concurrency cap; excess jobs queue (FIFO).

    Submission raises ``QUEUE_FULL`` immediately when the bounded queue is
    saturated and ``QUEUE_TIMEOUT`` if a queued job waits longer than the
    configured timeout.  Results are returned in input (submission) order, so
    scheduling timing never affects deterministic batch output.
    """

    def __init__(self, config: SchedulerConfig | None = None) -> None:
        self._config = config or SchedulerConfig()
        self._semaphore: asyncio.Semaphore | None = None
        self._inflight = 0

    def _sem(self) -> asyncio.Semaphore:
        if self._semaphore is None:
            self._semaphore = asyncio.Semaphore(self._config.max_concurrency)
        return self._semaphore

    async def run(
        self,
        fn: Callable[..., Awaitable[Any]],
        *args: Any,
        **kwargs: Any,
    ) -> Any:
        """Acquire a worker slot (queueing if needed) and run ``fn``."""
        sem = self._sem()
        capacity = self._config.max_concurrency + self._config.max_queue_depth
        if self._inflight >= capacity:
            raise OrbitError(
                ErrorCode.QUEUE_FULL,
                f"in-flight jobs reached capacity {capacity} "
                f"(concurrency {self._config.max_concurrency} + "
                f"queue {self._config.max_queue_depth})",
            )
        self._inflight += 1
        try:
            try:
                await asyncio.wait_for(
                    sem.acquire(), timeout=self._config.wait_timeout_s
                )
            except TimeoutError as exc:
                raise OrbitError(
                    ErrorCode.QUEUE_TIMEOUT,
                    "timed out after "
                    f"{self._config.wait_timeout_s}s waiting for a worker",
                ) from exc
            try:
                return await fn(*args, **kwargs)
            finally:
                sem.release()
        finally:
            self._inflight -= 1

    async def map_ordered(
        self,
        fn: Callable[[Any], Awaitable[Any]],
        items: list,
    ) -> list:
        """Run ``fn`` over items with bounded concurrency.

        All tasks are scheduled up to the concurrency/queue limits; results and
        exceptions are gathered in the original input order.
        """
        tasks = [
            asyncio.create_task(self.run(fn, item)) for item in items
        ]
        results = await asyncio.gather(*tasks, return_exceptions=True)
        ordered: list[Any] = []
        for result in results:
            if isinstance(result, BaseException):
                raise result
            ordered.append(result)
        return ordered

    @property
    def queue_depth(self) -> int:
        """Approximate number of jobs currently waiting for a worker slot."""
        return max(0, self._inflight - self._config.max_concurrency)
