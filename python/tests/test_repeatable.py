"""FEMM 5.0 must give the same numbers every time: the plan's "repeated runs bit-identical", and "worker counts bit-identical".

FEMM 5.0's own readme says its results are "bit-for-bit repeatable from run to run". These tests do not take that on trust:
they solve the default motor again and compare every number with ==, not with a tolerance, and compare the solution files
byte for byte. They start real FEMMs (a few seconds in all), so they are marked `femm`. Files go to P\\out\\tests\\.
Run from P:   uv run pytest tests/test_repeatable.py -q
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pytest

import femm_test_tools as tools

OUT = Path(__file__).resolve().parent.parent / "out" / "tests" / "repeatable"
POINTS = [(0.0, 0.0, 0.0), (0.0, 0.0, 10.0), (7.0, -6.0, 8.0), (13.0, 6.0, 8.0), (29.0, -4.7, 8.8), (41.0, 0.0, 10.0)]   # (theta_e, id, iq)


@pytest.mark.femm
def test_a_second_femm_gives_the_same_numbers_to_the_last_bit():
    first = tools.default_model_points(OUT, "first", POINTS)
    second = tools.default_model_points(OUT, "second", POINTS)
    assert first[0] == second[0]                       # every float of every point, compared with ==
    assert first[1] == second[1]                       # and the last solution file, byte for byte


@pytest.mark.femm
def test_the_fine_mesh_is_just_as_repeatable():
    first = tools.default_model_points(OUT, "fine_first", POINTS[:3], fine=True)
    second = tools.default_model_points(OUT, "fine_second", POINTS[:3], fine=True)
    assert first == second


@pytest.mark.femm
@pytest.mark.parametrize("workers", [2, 4])
def test_workers_that_run_at_the_same_time_get_the_serial_numbers(workers):
    """worker_scaling.py's claim, in small: the answer does not depend on how many FEMMs are busy next to each other."""
    serial = tools.default_model_points(OUT, "serial", POINTS)
    tasks = [(str(OUT), f"w{workers}_{k}", POINTS, False) for k in range(workers)]
    with ProcessPoolExecutor(workers) as pool:
        results = list(pool.map(tools.default_model_worker, tasks))
    assert len(results) == workers
    for numbers, ans_hash in results:
        assert numbers == serial[0]
        assert ans_hash == serial[1]
