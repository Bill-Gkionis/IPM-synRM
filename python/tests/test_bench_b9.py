"""Tests for B9 (bench/b9_rules.py): the rules themselves, shown with two bench runs at once. Everything that needs neither a solver nor a real bench run: the command line of a
child, starting and ending a child (with stand-in processes, and the rule that only B9's own children are ended), reading what a run printed (its pid, the waiting lines), the
completeness of a run's results, the comparison of two runs' numbers, the rules read from runs.json (priority, workers, threads), serialization on one root, the picture, the checks.

Run from P:   uv run pytest tests/test_bench_b9.py -q      (about 10 s: a few stand-in processes; no solver is run, no slot is taken)

In groups:
  1. The children: their command, starting them, waiting for them with a timeout, ending one that overran, and never ending anything else.
  2. Reading what a run said: the pid line, the waiting lines of the slot and of a root.
  3. Judging the runs: complete results, numbers that agree (exactly for B3, to a millionth for B6, a number that is not finite is a difference), the rules in runs.json,
     serialization on one root, who was left.
  4. The picture (a colour of its own for every kind of work) and the six checks.
  5. The case as a whole, run() with stand-in children (fix round 1 of Task 13): the watchers cover both pairs, the rules are read from the same-root run too, a child that cannot be
     started does not end the one that was, B3 is compared exactly and B6 is not.
"""
import ast
import collections
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from bench_test_tools import P, use_config


# ----------------------------------------------------------------------------------------------- 1. the children
def test_a_child_is_the_bench_started_on_a_root_of_its_own_with_the_cases_asked_for():
    from bench import b9_rules as b9
    command = b9.child_command(Path("C:/out/bench/B9/run1"), ("B3", "B6"))
    flat = [str(part) for part in command]
    assert flat[0] == sys.executable and flat[1] == str(P / "bench" / "run_bench.py")
    assert flat[flat.index("--out") + 1] == str(Path("C:/out/bench/B9/run1")) and flat[flat.index("--only") + 1] == "B3,B6"
    assert "B9" not in flat[flat.index("--only") + 1].split(",")                                         # a child never starts B9: it would start two more
    assert b9.CHILD_CASES == ("B3", "B6") and b9.SAME_ROOT_CASES == ("B6",) and "B9" not in b9.CHILD_CASES + b9.SAME_ROOT_CASES


def test_two_stand_in_children_are_started_at_once_below_normal_and_their_output_is_kept_in_files(tmp_path):
    import bench
    from bench import b9_rules as b9
    code = "import os, sys, time; print('pid', os.getpid(), flush=True); print('to stderr', file=sys.stderr, flush=True); time.sleep(1.0)"
    children = [b9.start_child([sys.executable, "-I", "-c", code], tmp_path / f"c{i}.out.txt", tmp_path / f"c{i}.err.txt") for i in (1, 2)]
    try:
        time.sleep(0.4)
        assert all(child.poll() is None for child in children)                                           # both are running at the same time
        classes = {bench.class_name(bench.priority_class(child.pid)) for child in children}
        assert classes == {"below_normal"}, classes                                                      # "The children run at below-normal priority"
        codes = b9.wait_for_children(children, timeout_s=60)
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
                child.wait(30)
    assert codes == [0, 0]
    for i in (1, 2):
        assert (tmp_path / f"c{i}.out.txt").read_text().startswith("pid ") and "to stderr" in (tmp_path / f"c{i}.err.txt").read_text()


def test_a_child_that_overruns_is_ended_on_the_timeout_and_the_other_is_left_alone(tmp_path):
    from bench import b9_rules as b9
    slow = b9.start_child([sys.executable, "-I", "-c", "import time; time.sleep(120)"], tmp_path / "slow.out.txt", tmp_path / "slow.err.txt")
    quick = b9.start_child([sys.executable, "-I", "-c", "import time; time.sleep(0.5)"], tmp_path / "quick.out.txt", tmp_path / "quick.err.txt")
    started = time.monotonic()
    try:
        codes = b9.wait_for_children([slow, quick], timeout_s=2.0)
    finally:
        for child in (slow, quick):
            if child.poll() is None:
                child.kill()
                child.wait(30)
    assert codes[0] is None and codes[1] == 0                                                            # None: it did not finish by itself, it was ended; the quick one finished
    assert 1.5 < time.monotonic() - started < 30 and slow.poll() is not None                             # ended at the timeout, not before, and it is gone


def kill_for_cleanup(pid: int) -> None:
    """The tests' own clean-up of a process the test started (not B9's code: tests are not scanned)."""
    import signal
    try:
        os.kill(pid, signal.SIGTERM)
    except OSError:
        pass


CHAIN = ("import subprocess, sys, time\np = subprocess.Popen([sys.executable, '-I', '-c', 'import time; time.sleep(120)'])\n"
         "print(p.pid, flush=True)\ntime.sleep(120)\n")


def test_a_pid_that_is_not_below_the_child_is_never_ended(tmp_path):
    """The venv's python.exe is a launcher in front of the interpreter that runs the bench, so the pid a child prints is its grandchild's. end_child ends that one too, but only
    when it is really below the child it was given: a pid that is not (a wrong line in the output, a number reused) is left alone, even with its right creation time."""
    from checks.matlab_check import process_created
    from bench import b9_rules as b9
    child = subprocess.Popen([sys.executable, "-I", "-c", CHAIN], stdout=subprocess.PIPE, text=True)
    stranger = subprocess.Popen([sys.executable, "-I", "-c", "import time; time.sleep(120)"])
    below = None
    try:
        below = int(child.stdout.readline())
        b9.end_child(child, real_pid=stranger.pid, real_created=process_created(stranger.pid))
        assert child.poll() is not None                                                                  # the child itself was ended
        assert stranger.poll() is None                                                                   # the stranger was not: it is not below the child
    finally:
        for process in (child, stranger):
            if process.poll() is None:
                process.kill()
            process.wait(30)
        if below is not None:
            kill_for_cleanup(below)                                                                      # the grandchild of the ended child, which end_child was not told about


def test_the_process_below_the_child_that_ran_its_code_is_ended_with_it(tmp_path):
    from checks.matlab_check import process_created, process_table
    from bench import b9_rules as b9
    child = subprocess.Popen([sys.executable, "-I", "-c", CHAIN], stdout=subprocess.PIPE, text=True)
    grandchild = None
    try:
        grandchild = int(child.stdout.readline())
        assert grandchild in process_table()
        b9.end_child(child, real_pid=grandchild, real_created=process_created(grandchild))             # the creation time the tracker recorded for it
        deadline = time.monotonic() + 10
        while grandchild in process_table() and time.monotonic() < deadline:
            time.sleep(0.1)
        assert child.poll() is not None and grandchild not in process_table()                           # both are gone
    finally:
        if child.poll() is None:
            child.kill()
        child.wait(30)
        if grandchild is not None:
            kill_for_cleanup(grandchild)


def test_the_process_below_the_child_is_ended_only_if_it_is_the_one_the_tracker_recorded_whatever_the_pid_says(tmp_path):
    """end_child opens ONE handle to the pid, compares the creation time it reads through that handle with the record (the creation time the tracker has for the process), and ends
    the process through that handle (review M1). A pid that has been given to another process since has another creation time, and a record that was not made cannot be compared:
    in both cases the process below the child is left alone, and the child itself is still ended."""
    from checks.matlab_check import process_created, process_table
    from bench import b9_rules as b9
    for label, created in (("a wrong creation time", lambda true: true + 10_000_000), ("no creation time at all", lambda true: None)):
        child = subprocess.Popen([sys.executable, "-I", "-c", CHAIN], stdout=subprocess.PIPE, text=True)
        grandchild = None
        try:
            grandchild = int(child.stdout.readline())
            true = process_created(grandchild)
            assert true is not None and grandchild in process_table()
            b9.end_child(child, real_pid=grandchild, real_created=created(true))
            assert child.poll() is not None, label                                                       # the child is ended whatever
            time.sleep(0.5)
            assert grandchild in process_table(), label                                                  # the process below it is not: it cannot be shown to be the recorded one
        finally:
            if child.poll() is None:
                child.kill()
            child.wait(30)
            if grandchild is not None:
                kill_for_cleanup(grandchild)


# ----------------------------------------------------------------------------------------------- 2. what a run said
RUN_OUTPUT = """bench: root D:\\Storage\\Motor Design\\Python to FEMM\\out\\bench\\B9\\run2 (pid 4321)
bench: 2 case(s), one at a time: B3, B6
=== B3
waiting for the solver slot: held by 'bench B6 ROSS' (pid 1234) since 14:02
[femm_pool bench B1 FEMM] starting 8 worker(s), pids [1, 2, 3], at below_normal priority, each with its own FEMM
=== B6
waiting for the bench root D:\\x: held by 'bench run B3,B6' (pid 777) since 09:05
"""


def test_the_pid_of_a_run_is_read_from_its_first_line_and_the_waiting_lines_name_their_holders():
    from bench import b9_rules as b9
    assert b9.run_pid(RUN_OUTPUT) == 4321 and b9.run_pid("nothing here") is None and b9.run_pid("bench: root x (pid 12a)") is None
    lines = b9.waiting_lines(RUN_OUTPUT)
    assert lines == [{"kind": "slot", "label": "bench B6 ROSS", "pid": 1234, "since": "14:02", "root": None, "line": "waiting for the solver slot: held by 'bench B6 ROSS' (pid 1234) since 14:02"},
                     {"kind": "root", "label": "bench run B3,B6", "pid": 777, "since": "09:05", "root": "D:\\x",
                      "line": "waiting for the bench root D:\\x: held by 'bench run B3,B6' (pid 777) since 09:05"}]
    assert b9.waiting_lines("a run whose note (holder.json) cannot be read\nwaiting for the solver slot: held by a run whose note (holder.json) cannot be read") == [
        {"kind": "slot", "label": None, "pid": None, "since": None, "root": None, "line": "waiting for the solver slot: held by a run whose note (holder.json) cannot be read"}]


def test_the_final_line_of_a_run_says_how_it_ended():
    from bench import b9_rules as b9
    assert b9.final_line("x\nbench: PASS, all 23 gating checks pass (47 s)\n") == "bench: PASS, all 23 gating checks pass (47 s)"
    assert b9.final_line("bench: FAIL, 1 of 22 gating checks fail: B3 (1) (47 s)\n\n").startswith("bench: FAIL")
    assert b9.final_line("no verdict") is None


# ----------------------------------------------------------------------------------------------- 3. judging the runs
def test_a_complete_case_folder_has_its_results_pictures_and_logs_and_a_damaged_one_says_what_is_missing(tmp_path):
    from bench import b9_rules as b9
    folder = tmp_path / "B6"
    folder.mkdir()
    checks = [{"name": f"c{i}", "verdict": "PASS"} for i in range(8)]
    (folder / "results.json").write_text(json.dumps({"case": "B6", "checks": checks}), encoding="utf-8")
    (folder / "figure.png").write_bytes(b"x" * 40_000)
    (folder / "run.log").write_text("log")
    (folder / "runs.json").write_text("[]")
    assert b9.case_problems(folder, 8) == []
    (folder / "figure.png").write_bytes(b"x" * 100)
    assert any("figure.png" in p for p in b9.case_problems(folder, 8))                                   # a picture that is not one
    (folder / "figure.png").write_bytes(b"x" * 40_000)
    assert any("8" in p and "checks" in p for p in b9.case_problems(folder, 9))                         # fewer checks than the case has
    (folder / "results.json").write_text(json.dumps({"case": "B6", "error": "boom", "checks": [{"name": "(did not finish)", "verdict": "FAIL"}] + checks}), encoding="utf-8")
    assert any("did not finish" in p for p in b9.case_problems(folder, 8))
    (folder / "run.log").unlink()
    (folder / "results.json").unlink()
    problems = b9.case_problems(folder, 8)
    assert any("run.log" in p for p in problems) and any("results.json" in p for p in problems)
    assert b9.case_problems(tmp_path / "missing", 8) and "folder" in b9.case_problems(tmp_path / "missing", 8)[0]


def test_two_roots_hold_the_same_files_or_the_difference_is_named(tmp_path):
    from bench import b9_rules as b9
    for root in ("one", "two"):
        (tmp_path / root / "B3" / "sub").mkdir(parents=True)
        (tmp_path / root / "B3" / "results.json").write_text("{}")
        (tmp_path / root / "B3" / "sub" / "x.sif").write_text("x")
    assert b9.tree_names(tmp_path / "one" / "B3") == ["results.json", "sub\\x.sif"] and b9.files_differ(tmp_path / "one" / "B3", tmp_path / "two" / "B3") == []
    (tmp_path / "two" / "B3" / "extra.txt").write_text("e")
    (tmp_path / "one" / "B3" / "sub" / "x.sif").unlink()
    differ = b9.files_differ(tmp_path / "one" / "B3", tmp_path / "two" / "B3")
    assert any("extra.txt" in d and "only" in d for d in differ) and any("x.sif" in d for d in differ)


def test_numbers_are_compared_leaf_by_leaf_and_what_belongs_to_one_run_only_is_left_out():
    """The two runs of a deterministic tool give the same numbers. The leaves that are THIS run's (seconds, process ids, how many looks the watcher took, the time of day) differ
    by nature and are left out; every other number in results.json is compared, and a changed one, or one that only a run has, is named."""
    from bench import b9_rules as b9
    a = {"psi": 16.62, "nodes": 15859, "ok": True, "name": "x", "seconds_total": 12.3, "pool": {"worker_pids": [1, 2], "seconds": 1.0, "time": "09:00", "workers": 8},
         "looks_at_the_processes": 500, "list": [{"v": 1.0}, {"v": 2.0}], "checks": [{"measured": 3.5, "difference_percent": 0.2}], "missing": 1.0}
    b = {"psi": 16.62, "nodes": 15859, "ok": True, "name": "x", "seconds_total": 99.9, "pool": {"worker_pids": [9, 8], "seconds": 7.0, "time": "09:05", "workers": 8},
         "looks_at_the_processes": 80, "list": [{"v": 1.0}, {"v": 2.0}], "checks": [{"measured": 3.5, "difference_percent": 0.2}]}
    flat = b9.flatten_numbers(a)
    assert "psi" in flat and "nodes" in flat and "pool.workers" in flat and "list[1].v" in flat and "checks[0].measured" in flat
    assert not any(k for k in flat if "seconds" in k or "pids" in k or k.endswith(".time") or "looks" in k) and "ok" not in flat and "name" not in flat     # bools and text are not numbers
    same = b9.compare_runs(a, b, 1e-6)
    assert same["different"] == [] and same["only_in_first"] == ["missing"] and same["only_in_second"] == [] and same["compared"] == len(flat) - 1
    assert same["exact"] == same["compared"] and same["worst_relative"] == 0.0
    c = json.loads(json.dumps(b))
    c["psi"] = 16.62 * (1 + 1e-5)
    c["list"][1]["v"] = 2.5
    changed = b9.compare_runs(a, c, 1e-6)
    assert [d[0] for d in changed["different"]] == ["list[1].v", "psi"] or [d[0] for d in changed["different"]] == ["psi", "list[1].v"]
    assert changed["worst_relative"] == pytest.approx(0.2, rel=1e-6)
    near = json.loads(json.dumps(b))
    near["psi"] = 16.62 * (1 + 1e-12)
    result = b9.compare_runs(a, near, 1e-6)
    assert result["different"] == [] and result["exact"] == result["compared"] - 1 and 0 < result["worst_relative"] <= 1e-9        # one ulp-sized difference: within the tolerance, not exact
    assert [d[0] for d in b9.compare_runs(a, near, 0.0)["different"]] == ["psi"]                           # and a difference at all when the tolerance is 0.0: exact
    assert [d[0] for d in b9.compare_runs(a, near)["different"]] == ["psi"]                                # exact is what a call that gives no tolerance means
    assert b9.compare_runs(a, c, tolerance=1e-4)["different"] == [("list[1].v", 2.0, 2.5)]                 # the tolerance is a share of the numbers themselves


@pytest.mark.parametrize("tolerance", [0.0, 1e-6])
def test_a_number_that_is_not_finite_is_a_difference_unless_both_runs_have_the_same_one(tolerance):
    """NaN and infinity passed compare_runs as agreement (review Important 2): relative = |x - y| / scale is NaN when either is, and 'NaN > tolerance' is False. A run that produced a
    NaN where the other produced a number is the nondeterministic failure the check exists to catch. Non-finite is a difference unless both are NaN or both are the same infinity."""
    from bench import b9_rules as b9
    nan, inf = float("nan"), float("inf")

    def compared(x, y):
        return b9.compare_runs({"v": x, "w": 1.0}, {"v": y, "w": 1.0}, tolerance)
    for x, y in ((nan, 5.0), (5.0, nan), (inf, 5.0), (5.0, inf), (-inf, 5.0), (5.0, -inf), (inf, -inf), (-inf, inf), (nan, inf), (inf, nan), (nan, 0.0), (0.0, nan), (inf, 0.0), (0.0, -inf)):
        result = compared(x, y)
        assert [d[0] for d in result["different"]] == ["v"] and result["worst_relative"] == inf, (x, y, result)
        assert result["compared"] == 2 and result["exact"] == 1                                              # the other number is the same in both
    for x, y in ((nan, nan), (inf, inf), (-inf, -inf)):
        result = compared(x, y)
        assert result["different"] == [] and result["exact"] == 2 and result["worst_relative"] == 0.0, (x, y, result)
    parsed = [json.loads(text) for text in ('{"v": NaN, "w": 1.0}', '{"v": 5.0, "w": 1.0}', '{"v": Infinity, "w": 1.0}')]       # the literals Python's json reads: they are numbers too
    assert [d[0] for d in b9.compare_runs(parsed[0], parsed[1], tolerance)["different"]] == ["v"] and [d[0] for d in b9.compare_runs(parsed[2], parsed[1], tolerance)["different"]] == ["v"]


def test_b3_is_compared_exactly_and_b6_to_a_millionth_and_each_case_has_its_own_limit():
    """FEMM's, Elmer's and Gmsh's numbers (B3) are identical to the last bit from one run to the next, so they are compared exactly. ROSS's (B6) are not: see the check's text."""
    from bench import b9_rules as b9
    assert b9.TOLERANCE_BY_CASE == {"B3": 0.0, "B6": 1e-6} and set(b9.TOLERANCE_BY_CASE) == set(b9.CHILD_CASES)
    assert not hasattr(b9, "NUMBER_TOLERANCE")                                                           # one limit for both cases was the unsupported claim


def records_good():
    return [{"kind": "femm_pool", "label": "bench B3 FEMM", "workers": 8, "workers_limit": 8, "workers_that_solved": 8, "worker_priority_classes_seen": {"below_normal": 8},
             "femm_priority_classes_seen": {"below_normal": 8}},
            {"kind": "run_solver", "label": "bench B3 ElmerSolver", "program": "ElmerSolver.exe", "threads": 4, "threads_limit": 8, "priority_classes_seen": {"elmersolver.exe": ["below_normal"]},
             "looks_at_the_processes": 120},
            {"kind": "run_solver", "label": "bench B3 ElmerGrid", "program": "ElmerGrid.exe", "threads": 4, "threads_limit": 8, "priority_classes_seen": {}, "looks_at_the_processes": 9},
            {"kind": "in_process", "label": "the bench's own Python", "threads": 1, "threads_limit": 8, "priority_class_seen": "below_normal"}]


def test_the_rules_are_read_from_a_runs_log_below_normal_and_at_most_eight_workers_or_threads():
    from bench import b9_rules as b9
    good = b9.runs_comply(records_good(), 8)
    assert good["ok"] and good["problems"] == [] and good["pool_records"] == 1 and good["solver_records"] == 2 and good["in_process_records"] == 1
    assert good["observed"] == 3 and good["max_workers"] == 8 and good["max_threads"] == 4                  # the ElmerGrid run was too short to be seen: noted, not a failure

    def broken(index, **changes):
        records = records_good()
        records[index] = {**records[index], **changes}
        return b9.runs_comply(records, 8)
    assert not broken(0, workers=9)["ok"] and not broken(0, workers_that_solved=9)["ok"]                 # a ninth worker
    assert not broken(0, worker_priority_classes_seen={"below_normal": 7, "normal": 1})["ok"]            # one worker at normal
    assert not broken(0, femm_priority_classes_seen={"above_normal": 8})["ok"]
    assert not broken(0, worker_priority_classes_seen={})["ok"]                                          # nothing observed: it cannot pass
    assert not broken(1, threads=9)["ok"] and not broken(1, threads=0)["ok"]                            # nine threads, or a thread count that is not one
    assert not broken(1, priority_classes_seen={"elmersolver.exe": ["below_normal", "normal"]})["ok"]
    assert not broken(3, priority_class_seen="normal")["ok"] and not broken(3, priority_class_seen="unreadable")["ok"]
    assert not broken(3, threads=12)["ok"]
    assert not b9.runs_comply([{"kind": "mystery"}], 8)["ok"]
    only_unseen = [r for r in records_good() if r.get("looks_at_the_processes") == 9]
    assert not b9.runs_comply(only_unseen, 8)["ok"]                                                      # every program unobserved: the case shows nothing about priority


def test_two_runs_on_one_root_are_serialized_when_the_second_did_its_first_work_after_the_first_ended():
    from bench import b9_rules as b9
    assert b9.serialized(first_ended=100.0, second_first_work=100.4) and b9.serialized(100.0, 140.0)
    assert b9.serialized(100.0, 99.5)                                                                    # within a second of the end: file times and the clock are not that exact
    assert not b9.serialized(100.0, 90.0)                                                                # started its work 10 s before the first ended: not serialized
    assert not b9.serialized(100.0, None)                                                                # no work seen: it cannot be said


def test_the_holder_samples_are_folded_into_segments_one_for_each_stretch_of_one_holder():
    from bench import b9_rules as b9
    samples = [(0.0, None), (0.25, ("bench B3 FEMM", 11)), (0.5, ("bench B3 FEMM", 11)), (0.75, ("bench B6 ROSS", 22)), (1.0, ("bench B6 ROSS", 22)), (1.25, None)]
    segments = b9.holder_segments(samples)
    assert [(s["holder"], s["start"], s["end"]) for s in segments] == [(None, 0.0, 0.25), (("bench B3 FEMM", 11), 0.25, 0.75), (("bench B6 ROSS", 22), 0.75, 1.25), (None, 1.25, 1.25)]
    assert b9.holder_segments([]) == []
    assert b9.holders_seen(samples) == [("bench B3 FEMM", 11), ("bench B6 ROSS", 22)]                      # in the order first seen, a free slot left out


# ----------------------------------------------------------------------------------------------- 4. the picture and the checks
def test_the_picture_shows_who_held_the_slot_and_who_held_the_root_and_has_no_blank_margin(tmp_path):
    import matplotlib.image as image
    from bench import b9_rules as b9
    slot = [(0.0, None), (1.0, ("bench B3 FEMM", 11)), (5.0, ("bench B6 ROSS", 22)), (9.0, ("bench B3 ROSS", 11)), (12.0, ("bench B6 ROSS", 22)), (15.0, None)]
    root = [(0.0, ("bench run B6", 33)), (6.0, ("bench run B6", 44)), (12.0, None)]
    target = tmp_path / "figure.png"
    b9.fig_timeline(slot, root, {11: "run 1", 22: "run 2", 33: "run C", 44: "run D"}, b9.figure_title(), target)
    picture = image.imread(str(target))
    height, width = picture.shape[:2]
    assert width <= 1.1 * 11.5 * 130 and height > 0.25 * width
    white = (picture[..., :3] > 0.97).all(axis=2)
    count = 0
    for line in white:
        if not line.all():
            break
        count += 1
    assert count <= 0.04 * height                                                                         # no blank band at the top


def test_the_figure_title_fits_the_figure_and_is_not_rewrapped():
    """Its second line was 104 characters against bench.TITLE_CHARS (92), so wrap_title left an orphan word on a line of its own (review M7). Every line fits now and the title is
    kept as it is written, as B7's is."""
    import bench
    from bench import b9_rules as b9
    title = b9.figure_title()
    assert all(len(line) <= bench.TITLE_CHARS for line in title.split("\n")) and bench.wrap_title(title) == title
    for needle in ("two bench runs", "solver slot", "lock", "B3", "B6"):
        assert needle in title, needle


def test_every_kind_of_work_has_a_colour_of_its_own():
    """The legend reused colours (review C1): B3's FEMM reference and Elmer at size 0.5 were both blue, B3's Elmer and 'run B6' both green, its laminated Elmer and B6's ROSS both red.
    Every label gets a colour of its own, as many labels as there are, and the colours are far enough apart to tell."""
    import math
    from bench import b9_rules as b9
    for n in range(1, 16):
        labels = [f"bench kind {i}" for i in range(n)]
        colours = b9.label_colours(labels)
        assert list(colours) == labels and len(set(colours.values())) == n, n
        rgb = [tuple(int(c[i:i + 2], 16) for i in (1, 3, 5)) for c in colours.values()]
        assert all(math.dist(a, b) >= 50 for i, a in enumerate(rgb) for b in rgb[i + 1:]), n
    assert b9.label_colours(["a", "b", "c"]) == b9.label_colours(["a", "b", "c"])                          # the same labels, the same colours


def test_the_picture_gives_each_kind_of_work_one_colour_in_both_panels_and_no_colour_to_two():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from bench import b9_rules as b9
    labels = ["bench B3 FEMM reference", "bench B3 Elmer", "bench B3 Elmer laminated", "bench B3 Elmer size 0.5", "bench B3 Elmer plain", "bench B6 ROSS", "bench run B6", "bench B3 gmsh",
              "bench B3 convert"]
    slot, t = [(0.0, None)], 1.0
    for i, label in enumerate(labels[:7]):
        slot.append((t, (label, 100 + i % 2)))                                                            # two runs (pids 100 and 101) take turns
        t += 8.0
    slot.append((t, None))
    root = [(0.0, (labels[6], 300)), (9.0, (labels[7], 301)), (18.0, (labels[8], 300)), (27.0, (labels[5], 301)), (36.0, None)]       # "bench B6 ROSS" is in both panels
    fig = b9.draw_timeline(slot, root, {100: "run 1", 101: "run 2", 300: "run C", 301: "run D"}, b9.figure_title())
    try:
        by_panel = [[tuple(np.round(c.get_facecolor()[0][:3], 3)) for c in ax.collections] for ax in fig.axes]
        bars = [colour for panel in by_panel for colour in panel]
        assert len(bars) == 7 + 4 and len(set(bars)) == len(labels)                                       # 9 kinds of work, 9 colours; a kind that is in both panels has one
        legend = fig.legends[0]
        swatches = [tuple(np.round(matplotlib.colors.to_rgb(h.get_facecolor()), 3)) for h in legend.legend_handles]
        texts = [t.get_text() for t in legend.get_texts()]
        assert len(swatches) == len(set(swatches)) == len(labels) and set(texts) == {label.replace("bench ", "") for label in labels}
        assert set(bars) == set(swatches)                                                                 # every bar is a colour of the legend, and the legend has no other
        assert by_panel[0][5] == by_panel[1][3]                                                           # B6 ROSS: the sixth bar above, the fourth below
    finally:
        plt.close(fig)


def test_the_legend_of_the_picture_does_not_cross_the_axis_label_below_the_lower_panel():
    """The legend's first row crossed the x-axis label of the lower panel ('seconds since the first of the two was started') when there are many kinds of work (the re-review's
    look at the real figure): the margin under the panels was a fixed 17 % of the figure, enough for two rows of legend, and the real run has 13 labels in four rows. The margin
    follows the number of rows now. Measured on the drawn figure: the legend's box meets no tick label, no axis label, neither panel, and lies whole inside the figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from bench import b9_rules as b9
    labels = [f"bench B3 kind of work number {i}" for i in range(13)]
    slot, t = [(0.0, None)], 1.0
    for i, label in enumerate(labels):
        slot.append((t, (label, 100 + i % 2)))
        t += 8.0
    slot.append((t, None))
    root = [(0.0, (labels[0], 300)), (9.0, (labels[1], 301)), (18.0, None)]
    for n in (13, 4, 1):                                                                                   # 13 labels in four rows, 4 in one row, and two (the slot's one and the second that the root's lock draws)
        kinds = labels[:n]
        fig = b9.draw_timeline([(s, h if h is None or h[0] in kinds else (kinds[0], h[1])) for s, h in slot], root, {100: "run 1", 101: "run 2", 300: "run C", 301: "run D"}, b9.figure_title())
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            box = fig.legends[0].get_window_extent(renderer)
            lower = fig.axes[1]
            for artist in (lower.xaxis.label, *lower.get_xticklabels(), *fig.axes):
                other = artist.get_window_extent(renderer)
                assert not box.overlaps(other), (n, getattr(artist, "get_text", lambda: "a panel")(), box, other)
            assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.x1 and box.y1 <= fig.bbox.y1, (n, box, fig.bbox)
        finally:
            plt.close(fig)


def good_context():
    """What make_checks needs, for runs in which everything is as the brief and the rulings say."""
    return {
        "limit": 8,
        "a": {"pids": [4001, 4002], "returncodes": [0, 0], "timed_out": [False, False], "seconds": 101.0,
              "waits": [[], [{"kind": "slot", "label": "bench B3 FEMM", "pid": 4001, "since": "09:05", "root": None, "line": "waiting for the solver slot: held by 'bench B3 FEMM' (pid 4001) since 09:05"}]],
              "final_lines": ["bench: PASS, all 23 gating checks pass (50 s)", "bench: PASS, all 23 gating checks pass (101 s)"]},
        "complete": {"problems": [], "cases_checked": 4, "checks_per_case": {"B3": 15, "B6": 8}},
        "numbers": {"compared": 600, "exact": 598, "worst_relative": 3.1e-8, "different": [], "only_in_first": [], "only_in_second": [],
                    "per_case": {"B3": {"compared": 350, "exact": 350, "worst_relative": 0.0, "tolerance": 0.0}, "B6": {"compared": 250, "exact": 248, "worst_relative": 3.1e-8, "tolerance": 1e-6}}},
        "left": {"alive_after": [], "members": 40, "failed_looks": 0, "looks": 2000, "slot_free": True, "names": ["python.exe", "femm_headlessw.exe", "fkn.exe", "elmersolver.exe"],
                 "slot_looks": 520, "slot_unreadable": 0, "lock_looks": 130, "lock_unreadable": 0},
        "rules": {"problems": [], "records": 14, "observed": 12, "max_workers": 8, "max_threads": 4, "classes_seen": {"python.exe": ["below_normal"], "femm_headlessw.exe": ["below_normal"]},
                  "peaks": {"femm_headlessw.exe": 8, "python.exe": 18}, "looks": 4000},
        "same": {"pids": [5001, 5002], "returncodes": [0, 0], "timed_out": [False, False], "c_ended_at": 200.0, "d_first_work_at": 200.6,
                 "waits": [{"kind": "root", "label": "bench run B6", "pid": 5001, "since": "09:10", "root": "D:\\x", "line": "waiting for the bench root D:\\x: held by 'bench run B6' (pid 5001) since 09:10"}],
                 "d_problems": [], "c_final": "bench: PASS, all 7 gating checks pass (24 s)", "d_final": "bench: PASS, all 7 gating checks pass (49 s)", "lock_free": True},
    }


def verdicts(context):
    from bench import b9_rules as b9
    return {c["name"]: c for c in b9.make_checks(context)}


def test_the_six_checks_of_the_ruling_are_there_in_order_and_all_pass_when_everything_is_as_ruled():
    from bench import b9_rules as b9
    checks = b9.make_checks(good_context())
    names = [c["name"] for c in checks]
    assert len(names) >= 6 and len(names) == len(set(names))
    for position, needle in enumerate(("waiting line", "complete", "agree", "solver process", "below-normal", "same root")):
        assert needle in names[position], (position, needle, names[position])
    assert all(c["verdict"] == "PASS" for c in checks if c.get("gating", True)), [n for n in names if verdicts(good_context())[n]["verdict"] != "PASS"]


def test_each_of_the_six_checks_turns_red_for_the_failure_it_is_about():
    def failed(**changes):
        context = good_context()
        for key, value in changes.items():
            context[key] = {**context[key], **value}
        return {n for n, c in verdicts(context).items() if c["verdict"] == "FAIL" and c.get("gating", True)}
    assert any("waiting line" in n for n in failed(a={"waits": [[], []]}))                                 # nobody waited
    assert any("waiting line" in n for n in failed(a={"waits": [[], [{"kind": "slot", "label": "x", "pid": 99, "since": "09:05", "root": None, "line": "..."}]]}))   # it named a pid that is not the other run's
    assert any("waiting line" in n for n in failed(a={"waits": [[], [{"kind": "root", "label": "x", "pid": 4001, "since": "09:05", "root": "y", "line": "..."}]]}))   # the root's line is not the slot's
    assert any("complete" in n for n in failed(complete={"problems": ["run2\\B6: figure.png is 100 bytes"]}))
    assert any("complete" in n for n in failed(a={"timed_out": [False, True]}))                            # a child had to be ended
    assert any("agree" in n for n in failed(numbers={"different": [("B3 psi", 1.0, 1.1)], "worst_relative": 0.1}))
    assert any("agree" in n for n in failed(numbers={"only_in_second": ["x"]}))
    assert any("agree" in n for n in failed(numbers={"compared": 0, "exact": 0}))                          # nothing compared: it cannot pass
    assert any("solver process" in n for n in failed(left={"alive_after": ["fkn.exe"]}))
    assert any("solver process" in n for n in failed(left={"slot_free": False}))
    assert any("solver process" in n for n in failed(left={"failed_looks": 5}))                            # the process list could not be read: it cannot pass
    assert any("below-normal" in n for n in failed(rules={"problems": ["B3: 9 workers"]}))
    assert any("below-normal" in n for n in failed(rules={"max_workers": 9}))
    assert any("below-normal" in n for n in failed(rules={"peaks": {"femm_headlessw.exe": 9}}))            # nine FEMMs alive at one look, whatever the logs say
    assert any("below-normal" in n for n in failed(rules={"classes_seen": {"python.exe": ["below_normal", "normal"]}}))
    assert any("below-normal" in n for n in failed(rules={"classes_seen": {}}))
    assert any("same root" in n for n in failed(same={"waits": []}))
    assert any("same root" in n for n in failed(same={"d_first_work_at": 150.0}))                          # the second did its first work while the first was still going
    assert any("same root" in n for n in failed(same={"d_problems": ["B6: results.json is not there"]}))
    assert any("same root" in n for n in failed(same={"timed_out": [True, False]}))
    assert any("same root" in n for n in failed(same={"lock_free": False}))


def test_the_check_of_the_numbers_says_what_is_known_about_ross_and_nothing_more():
    """The text claimed more than the evidence shows (review Important 4): '(alone, too...)', 'alone as well as side by side', 'it is ROSS's eigen-solver'. What is known: once, in one of
    the first concurrent runs, B6 differed from a solo B6 by up to 4.9e-8 in the modes near 450,000 rpm; two other concurrent pairs were identical and in a later one B6 differed again
    by up to 4.8e-8, in the same modes and with the same two values the eight solo runs at threads=2 gave (451,957.695 and 451,957.717 rpm for the highest); ROSS fixes its ARPACK start vector
    (v0 = ones), so a random start is not the cause; and 16 solo runs of B6's ROSS step (8 with threads=2, 8 with threads=1) took two values in each group: setting OMP_NUM_THREADS=1
    through run_solver did not remove the variation (the count was not read back from ROSS's BLAS, so that is all it shows); the cause is not isolated. The 1e-6 was set after that
    first run. The experiment is described by what it was, not by a name that means something only in the plan's workspace; and only numbers are compared, not labels (re-review)."""
    item = next(c for n, c in verdicts(good_context()).items() if "agree" in n)
    text = item["expected"] + " | " + item["measured"]
    for gone in ("(alone, too", "alone as well as side by side", "it is ROSS's eigen-solver", "not the two runs meeting", "so it is ROSS", "experiment E2", "fix round", "E2",
                 "The thread count is not the cause", "later concurrent pairs were identical"):
        assert gone not in text, gone
    for kept in ("exactly", "B3", "B6", "4.9e-8", "450,000", "once", "not isolated", "v0", "set after", "1e-06", "4.3e-8"):
        assert kept in text, kept
    for kept in ("16 solo runs of B6's ROSS step on 2026-11-06", "8 with threads=2, 8 with threads=1", "one after another through the slot",
                 "Setting OMP_NUM_THREADS=1 (run_solver's threads=1) did not remove the variation"):                       # what was done, said as what was done
        assert kept in text, kept
    for kept in ("two other concurrent pairs were identical", "in a later one B6 differed again by up to 4.8e-8, in the same modes",
                 "the same two values the eight solo runs at threads=2 gave (451,957.695 and 451,957.717 rpm for the highest; with threads=1 the solo runs gave 451,957.706 and 451,957.686)"):
        assert kept in text, kept                                                                                            # what the fourth real pair showed (T13-P3: threads=2 is B6's own)
    assert "45 times the 10,000 rpm B6 is about" in text and "this case is about" not in text and "the same two values the solo runs gave" not in text
    for kept in ("Only numbers are compared", "labels such as whirl directions are not", "Mixed in 12 of the 16, Forward in 4", "about 452,000 rpm", "gate on that noise"):
        assert kept in text, kept                                                                                            # why the labels are left out
    assert "B3 350 of 350 identical (worst 0.0e+00, limit 0)" in item["measured"] and "B6 248 of 250 identical (worst 3.1e-08, limit 1e-06)" in item["measured"]
    # B3 is exact: a difference of one part in 10^12 in B3 fails, in B6 it does not
    context = good_context()
    context["numbers"] = {**context["numbers"], "different": [("B3 psi", 16.6, 16.6 * (1 + 1e-12))], "worst_relative": 1e-12}
    assert next(c for n, c in verdicts(context).items() if "agree" in n)["verdict"] == "FAIL"


def test_the_rules_check_judges_every_case_of_b9s_runs_and_the_checks_report_what_they_looked_at():
    """Check 5's scope (a design decision): the cases of the runs B9 starts (B3 and B6, both pairs), not B1-B8's logs in the default root; the name says so. Check 4 says how
    often the slot and the root's lock were looked at and how many of those looks failed (HolderSampler.errors, which nothing read)."""
    checks = verdicts(good_context())
    rules = next(c for n, c in checks.items() if "below-normal" in n)
    assert "every case of B9's runs" in next(n for n in checks if "below-normal" in n) and "B3" in rules["expected"] and "same" in rules["expected"]
    left = next(c for n, c in checks.items() if "solver process" in n)
    assert "the slot was looked at 520 times (0 unreadable) and the root's lock 130 times (0 unreadable)" in left["measured"]
    context = good_context()
    context["left"] = {**context["left"], "slot_unreadable": 3, "lock_unreadable": 1}
    assert "the slot was looked at 520 times (3 unreadable) and the root's lock 130 times (1 unreadable)" in next(c for n, c in verdicts(context).items() if "solver process" in n)["measured"]


def test_a_child_that_exits_one_is_not_a_failure_of_the_rules_the_runs_are_judged_by_completeness():
    """B2 exits 1 by a standing ruling; any case can. The exit status and the final line are reported, never gating."""
    context = good_context()
    context["a"] = {**context["a"], "returncodes": [1, 1], "final_lines": ["bench: FAIL, 1 of 22 gating checks fail: B3 (1)", "bench: FAIL, 1 of 22 gating checks fail: B3 (1)"]}
    checks = verdicts(context)
    assert all(c["verdict"] == "PASS" for c in checks.values() if c.get("gating", True))
    assert any(not c.get("gating", True) and "exit" in n for n, c in checks.items())


def test_b9_holds_no_slot_starts_only_the_bench_and_ends_only_what_it_started():
    """B9 must not hold the solver slot while its children run (they need it: they would wait for it until the timeout): no run_solver, femm_pool or solver_slot call in the module.
    It starts one kind of program, the bench (run_bench.py) through one function, and ends a process in one function that wait_for_children calls on the timeout and nothing else."""
    from test_bench import queue_calls
    from test_check_tools import scopes_calling
    source = (P / "bench" / "b9_rules.py").read_text(encoding="utf-8")
    assert queue_calls(source) == (0, [])
    for door in ("solver_slot", "femm_slot"):
        assert scopes_calling(source, "femm_queue", door) == []
    assert scopes_calling(source, None, "start_child")                                                     # the bench is started in this module, through start_child
    assert scopes_calling(source, None, "end_child") == ["wait_for_children"]
    assert scopes_calling(source, None, "wait_for_children") and "timeout_s" in source
    tree = ast.parse(source)
    popen = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == "Popen"]
    assert len(popen) == 1                                                                                 # one place that starts a program
    assert "BELOW_NORMAL_PRIORITY_CLASS" in source                                                         # and it starts it below normal
    assert scopes_calling(source, "subprocess", "Popen") == ["start_child"]
    # the enders, counted (review M2): the scan of tests/test_bench.py keeps a SET of (function, name), so a second kill or TerminateProcess in end_child would not show there
    from test_check_tools import ENDERS
    end_child = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "end_child")
    counted = lambda node: collections.Counter(n.attr for n in ast.walk(node) if isinstance(n, ast.Attribute) and n.attr in ENDERS) + collections.Counter(
        n.id for n in ast.walk(node) if isinstance(n, ast.Name) and n.id in ENDERS)
    assert counted(end_child) == counted(tree) == collections.Counter({"kill": 1, "TerminateProcess": 1}), (counted(end_child), counted(tree))
    # and start_child is only ever given a command made by child_command(...): the one kind of program B9 starts is the bench
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "start_child"]
    assert calls and all(c.args and isinstance(c.args[0], ast.Call) and isinstance(c.args[0].func, ast.Name) and c.args[0].func.id == "child_command" for c in calls)
    named = [n for n in ast.walk(tree) if isinstance(n, ast.Name) and n.id == "start_child" and not any(n is c.func for c in calls)]
    assert named == []                                                                                     # never passed around as a value, so no other way in



# ----------------------------------------------------------------------------------------------- 5. the case as a whole: run() with stand-in children
# run() starts two bench runs on two roots and two on one, watches the whole tree from outside and judges what the runs left. These tests drive the REAL run() with stand-ins for the
# children (a few lines of Python that print the first line a bench run prints and leave a result tree; the same-root pair uses the real root lock, so the second really waits and
# says who holds it), because what they show is how run() is put together: which stretch of the run the watchers cover, what is read from where, and what happens to a child that
# cannot be started. The checks themselves are tried on hand-made contexts above.
STAND_IN = r'''
import json, os, shutil, subprocess, sys, time
from pathlib import Path

root, cases, mode, project, slot_dir = Path(sys.argv[1]), sys.argv[2].split(","), sys.argv[3], sys.argv[4], sys.argv[5]
print(f"bench: root {root} (pid {os.getpid()})", flush=True)
root.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, project)
same = root.name == "same"
if mode == "hang_pair_a" and not same:
    time.sleep(300)


def write_tree():
    for case in cases:
        folder = root / case
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True)
        jitter = 0.0
        if root.name == "run2" and ((mode == "jitter_b6" and case == "B6") or (mode.startswith("jitter_b3") and case == "B3")):
            jitter = 1e-8 if case == "B6" else 1e-12
        results = {"case": case, "checks": [{"name": f"c{i}", "verdict": "PASS"} for i in range({"B3": 15, "B6": 8}[case])], "psi": 16.6 * (1 + jitter), "seconds_total": time.time()}
        (folder / "results.json").write_text(json.dumps(results), encoding="utf-8")
        (folder / "figure.png").write_bytes(b"x" * 30000)
        (folder / "run.log").write_text("log", encoding="utf-8")
        record = {"kind": "in_process", "label": "the bench's own Python", "threads": 1, "threads_limit": 8,
                  "priority_class_seen": "normal" if (same and "bad_records" in mode) else "below_normal"}
        (folder / "runs.json").write_text(json.dumps([record]), encoding="utf-8")
    (root / "summary.json").write_text(json.dumps({"cases": cases}), encoding="utf-8")


if same:
    from bench import hold_root                                      # the real root lock: the second run waits, and says who holds it
    with hold_root("bench run " + ",".join(cases), folder=root):
        write_tree()
        if "straggler" in mode:                                      # something this run started that outlives it
            straggler = subprocess.Popen([sys.executable, "-I", "-c", "import time; time.sleep(300)"])
            print("straggler", straggler.pid, flush=True)
        time.sleep(0.8)
else:
    import femm_queue                                                # the real solver slot (in a lock folder of the test's own): the second run waits, and says who holds it
    if root.name == "run2":
        time.sleep(0.5)
    with femm_queue.solver_slot("bench stand-in " + root.name, lock_dir=slot_dir):
        write_tree()
        time.sleep(1.0)
'''


@pytest.fixture
def stand_in(tmp_path, monkeypatch):
    """B9's run() with stand-in children. state["mode"] says what the stand-ins do; state["fail_call"] = n makes the n-th command a program that does not exist (so that
    starting it fails); the constants that are minutes long are seconds."""
    from bench import b9_rules as b9
    use_config(monkeypatch, tmp_path)
    script = tmp_path / "fake_run_bench.py"
    script.write_text(STAND_IN, encoding="utf-8")
    state = {"mode": "plain", "calls": 0, "fail_call": None, "tmp": tmp_path, "slot_dir": tmp_path / "out" / ".femm"}        # the lock folder B9's own look at the slot reads

    def child_command(root_folder, cases):
        state["calls"] += 1
        if state["fail_call"] == state["calls"]:
            return [str(tmp_path / "no_such_program.exe")]
        return [sys.executable, "-I", str(script), str(root_folder), ",".join(cases), state["mode"], str(P), str(state["slot_dir"])]
    monkeypatch.setattr(b9, "child_command", child_command)
    for name, value in (("GRACE_S", 0.3), ("PAIR_TIMEOUT_S", 120), ("SAME_TIMEOUT_S", 120), ("LOCK_WAIT_S", 30), ("POLL_S", 0.05)):
        monkeypatch.setattr(b9, name, value)
    yield state
    for name in ("c.out.txt", "d.out.txt"):                                                              # the stand-ins' stragglers are the test's to end
        path = tmp_path / "B9" / "same" / name
        if path.is_file():
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("straggler "):
                    kill_for_cleanup(int(line.split()[1]))


def run_b9(tmp_path, monkeypatch=None, draw=True):
    """Run B9's run() in tmp_path/B9; (its results, {check name: check})."""
    import bench
    from bench import b9_rules as b9
    if not draw and monkeypatch is not None:
        monkeypatch.setattr(b9, "fig_timeline", lambda *args, **kwargs: None)                              # the picture is tried in the first of these runs only
    folder = tmp_path / "B9"
    folder.mkdir()
    results = b9.run(folder, bench.Log(folder))
    return results, {c["name"]: c for c in results["checks"]}


def failed(checks):
    return {n: c for n, c in checks.items() if c["verdict"] == "FAIL" and c.get("gating", True)}


def test_run_with_stand_in_children_passes_every_check_and_the_two_runs_of_b6_may_differ_by_a_hundredth_of_its_limit(stand_in, tmp_path):
    """The whole case on stand-ins that behave: six checks pass, the picture is drawn, B3's numbers are identical and B6's differ by 1e-8 (under its limit of 1e-6)."""
    stand_in["mode"] = "jitter_b6"
    results, checks = run_b9(tmp_path)
    assert not failed(checks), {n: c["measured"] for n, c in failed(checks).items()}
    assert (tmp_path / "B9" / "figure.png").stat().st_size > 5_000
    per_case = results["numbers"]["per_case"]
    assert per_case["B3"]["worst_relative"] == 0.0 and per_case["B3"]["tolerance"] == 0.0
    assert 0 < per_case["B6"]["worst_relative"] <= 1.1e-8 and per_case["B6"]["tolerance"] == 1e-6
    assert results["pair_a"]["returncodes"] == [0, 0] and results["pair_b"]["returncodes"] == [0, 0]
    assert results["rules"]["records"] == 5                                                                 # four case folders in the two roots, one in the same root


def test_a_leftover_of_the_same_root_pair_and_the_rules_of_its_log_are_seen_and_so_is_a_difference_in_b3(stand_in, tmp_path, monkeypatch):
    """Review Important 1: the watchers were stopped, and 'alive' and 'slot free' read, before the pair on one root ran, so a process left by that pair was never seen, and its runs.json
    was not read for the rules. Here the same-root runs each leave a process running and write a record at normal priority; B3's second run differs from its first by one part in 10^12
    (exact for B3). Three checks must fail, each for its own reason."""
    stand_in["mode"] = "jitter_b3_straggler_bad_records"
    results, checks = run_b9(tmp_path, monkeypatch, draw=False)
    bad = failed(checks)
    agree = next(c for n, c in bad.items() if "agree" in n)
    assert "B3 psi" in agree["measured"]                                                                  # B3 differs by 1e-12 and is compared exactly
    left = next(c for n, c in bad.items() if "solver process" in n)
    assert "python.exe" in left["measured"].split("alive after:")[1].split(";")[0], left["measured"]       # the stragglers of the second pair
    rules = next(c for n, c in bad.items() if "below-normal" in n)
    assert "same" in rules["measured"] and "priority normal" in rules["measured"], rules["measured"]       # the same-root run's runs.json was read
    assert len(bad) == 3, list(bad)


def test_if_the_second_child_of_the_pair_cannot_be_started_the_first_is_waited_for_not_ended_and_the_step_is_named(stand_in, tmp_path, monkeypatch):
    """Review M1: a start that raised left the first child running on its own and ended the whole case with a traceback. Now the first child is waited for, with the same timeout
    (it is ended only on that timeout, as before), and the case FAILS naming the step."""
    stand_in["fail_call"] = 2                                                                             # the second child of the pair on two roots
    results, checks = run_b9(tmp_path, monkeypatch, draw=False)
    assert results["pair_a"]["returncodes"][0] == 0 and results["pair_a"]["returncodes"][1] != 0            # the first finished by itself (not None: not ended); the second never started
    assert (tmp_path / "B9" / "run1" / "B3" / "results.json").is_file()                                   # and left its results
    complete = next(c for n, c in failed(checks).items() if "complete" in n)
    assert "starting run 2 failed" in complete["measured"] and "FileNotFoundError" in complete["measured"], complete["measured"]


def test_if_the_second_run_on_the_same_root_cannot_be_started_the_first_is_waited_for_and_the_step_is_named(stand_in, tmp_path, monkeypatch):
    stand_in["fail_call"] = 4                                                                             # the second child of the pair on one root
    results, checks = run_b9(tmp_path, monkeypatch, draw=False)
    assert results["pair_b"]["returncodes"][0] == 0 and results["pair_b"]["returncodes"][1] != 0
    same = next(c for n, c in failed(checks).items() if "same root" in n)
    assert "starting the second run failed" in same["measured"] and "FileNotFoundError" in same["measured"], same["measured"]
    assert not [n for n in failed(checks) if "same root" not in n]                                         # nothing else is wrong


def test_children_that_overrun_are_ended_with_the_interpreter_below_them_and_nothing_is_left(stand_in, tmp_path, monkeypatch):
    """The timeout path through run(): both children of the pair on two roots sleep, B9 ends them at its timeout (end_child: the launcher by its own handle, the interpreter below it
    through the handle it opens, if its creation time is the one the tracker recorded) and the case sees that nothing of the tree is left."""
    from bench import b9_rules as b9
    monkeypatch.setattr(b9, "PAIR_TIMEOUT_S", 6)
    stand_in["mode"] = "hang_pair_a"
    results, checks = run_b9(tmp_path, monkeypatch, draw=False)
    assert results["pair_a"]["returncodes"] == [None, None]                                               # ended on the timeout, not finished
    left = next(c for n, c in checks.items() if "solver process" in n)
    assert left["verdict"] == "PASS" and results["processes_below_b9"]["alive_after"] == [], left["measured"]       # the interpreters are gone too
    assert any("complete" in n for n in failed(checks))                                                   # and the timeout is a failure of completeness


def test_the_docstring_and_the_comment_name_the_solo_runs_that_gave_the_two_values():
    """T13-P3. 'the same two values the solo runs gave' was true of the eight solo runs at threads=2 (B6's own setting) and not of those at threads=1 (451,957.706 and 451,957.686 rpm).
    The docstring's numbers-to-expect and the comment above TOLERANCE_BY_CASE say which; and 'The six checks' names the rule, not a file of the plan's workspace."""
    from bench import b9_rules as b9
    source = (P / "bench" / "b9_rules.py").read_text(encoding="utf-8")
    assert source.count("the two values the eight solo runs at threads=2 gave") + source.count("the same two values the eight solo runs at threads=2 gave") == 3, "docstring, comment and check text"
    assert "the two values the solo runs gave" not in source
    assert "The six checks (the plan's B9 row" in b9.__doc__
