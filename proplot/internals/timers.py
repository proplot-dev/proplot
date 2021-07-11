#!/usr/bin/env python3
"""
Utilities for timing ProPlot performance.
"""
import time

BENCHMARK = False  # toggle this to turn on benchmarking (see timers.py)


class _benchmark(object):
    """
    Context object for timing arbitrary blocks of code.
    """
    def __init__(self, message):
        self.message = message

    def __enter__(self):
        if BENCHMARK:
            self.time = time.perf_counter()

    def __exit__(self, *args):  # noqa: U100
        if BENCHMARK:
            print(f'{self.message}: {time.perf_counter() - self.time}s')
