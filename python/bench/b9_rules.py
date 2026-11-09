"""B9: the rules themselves: two bench runs started at once, in two processes, on a machine that has one solver slot.

What it shows, and why
    R10 says: one heavy run at a time on this PC, every solver at below-normal priority, at most 8 workers or threads, a timeout on every run, nothing ended but what a run started.
    Task 3 built it (femm_queue) and Task 4 and Tasks 11 to 13 use it, one case after another. The rules were made for the day two runs meet, so B9 makes them meet: it starts two
    `run_bench.py` processes at once, on the same light cases, and looks at what happens. The worry was "overlapping or damaging the machine"; the case is sound when the second
    run waits for the slot and says who holds it, when both runs' results are complete and agree, when nothing of theirs is left running, when every log and every process shows
    below-normal priority and at most 8 workers or threads, and when two runs on the SAME root take turns instead of emptying each other's folders.

Method (the files and the protocol)
    * Two roots. `run_bench.py --out ROOT` (bench/run_bench.py, this task): every case writes under ROOT\\<case>\\, and fresh() refuses anything outside the active root. B9's two
      children get out\\bench\\B9\\run1 and run2, so that two runs never write into one folder, and `--only` rewrites the summary.json of their own root only. They run B3 (Gmsh,
      ElmerGrid, ElmerSolver and a FEMM pool: 26 s alone) and B6 (ROSS in its own environment: 21 s): light, they take the solver slot many times, they use a pool (8 workers), run_solver
      (threads) and the in-process Gmsh library, so every kind of record of the log is there. B2 exits 1 by a standing ruling and is left out; every run is judged by the
      completeness of its results, not by its exit status (which is reported). B9 is not among them: a child that started B9 would start two more.
    * The children are the one program B9 starts itself. run_solver would hold the solver slot for the length of the call, and the children need the slot: they would wait for it,
      and for each other, until the timeout. So B9 starts them directly (the named exception of tests/test_bench.py, one function: start_child, only ever given a command made by
      child_command), below normal, with their output in files, and waits for both with a timeout. It ends a child only on that timeout, and only the child it started (end_child:
      the child by the handle Popen holds, and the interpreter below it that ran its code, if it is really below the child -- the venv's python.exe is a launcher in front of the
      interpreter -- through ONE handle opened to its pid, after reading through that handle that it was created when the tracker recorded it: no check-then-kill by pid). A child
      that cannot be started does not end the one that was: that one is waited for, with the same timeout, and the case FAILS naming the step. A killed run frees everything by
      itself: its pool workers watch their owner, and run_solver's job closes with its caller.
    * Watched from outside while they run, from before the first child to after the last one has ended (both pairs, one start and one stop): bench.ChildWatcher reads the priority
      class of every process below B9 (and counts the most of each name alive at one look, which is how 'at most 8 workers' is read from outside), checks.matlab_check.JobTracker
      follows every process that descends from B9 (known by pid and creation time) and says at the end which are still alive, and a bench.HolderSampler looks at
      femm_queue.holder() every 0.25 s (and, for the pair on one root, at that root's lock): the picture is its record of who held them. What is alive, and whether the slot is free,
      are read once, at the end.
    * What the runs said: each child prints "bench: root <ROOT> (pid N)" first (the interpreter that runs the bench is the process BELOW the launcher B9 started, so only the run
      can say its pid) and, when it has to wait, "waiting for the solver slot: held by '<label>' (pid N) since HH:MM". Those lines are read from the children's output.
    * Then two runs on ONE root, out\\bench\\B9\\same: the first is started, B9 waits until the root's lock note (holder.json) names it, and starts the second. The second must print
      "waiting for the bench root <ROOT>: held by 'bench run B6' (pid N) ...", must not touch the root's folders until the first has ended (its first file is newer than the first
      run's summary.json, written just before the lock was let go), and must finish complete.

The six checks (the plan's B9 row, made checkable)
    The plan's B9 row: the second run waits for the slot and prints the holder, no solver process is left afterwards, every case's log shows below-normal priority and at most
    8 workers or threads.
    1. the second run printed the slot's waiting line, naming the holder (a pid that is the other run's);
    2. both runs' results are complete: every check of the chosen cases is present, each case has its results.json, its pictures and its logs, the two roots hold the same files, and
       neither child had to be ended;
    3. the two runs' numbers agree (every number of every results.json except this run's seconds, process ids and looks at the process list): B3's (FEMM, Elmer, Gmsh) exactly,
       B6's (ROSS) to 1e-6 of themselves (TOLERANCE_BY_CASE says what is known of ROSS's modes near 450,000 rpm, which differed by up to 4.9e-8: not isolated; a number that is not
       finite is a difference);
    4. no solver process is left: nothing that descended from B9, in either pair, is alive after the last child ended (and a process list that could not be read does not count as
       empty), and the slot is free;
    5. every case of B9's runs shows below-normal priority and at most 8 workers or threads (runs.json of every case of both pairs, the one on a single root included), and from
       outside: every process under B9 was seen below normal, and never more than 8 FEMMs alive at one look;
    6. two runs on the same root are serialized by the root lock.

Run (from the Python to FEMM folder; about 3 minutes: the two runs take turns on the slot, 2 x 50 s, and the same-root pair 2 x 25 s)
    uv run python bench/run_bench.py --only B9
Writes out\\bench\\B9\\: results.json, figure.png (who held the slot and the root's lock, over time), run1\\ and run2\\ (a bench run each: B3\\, B6\\, summary.json, child.out.txt,
child.err.txt), same\\ (the pair on one root: B6\\, summary.json, .lock\\, c.out.txt, d.out.txt), run.log, runs.json.

Numbers to expect (this PC, 2026-11-06; four real runs, the last inside the whole bench; the case takes 148 to 154 s: two roots 100 to 102 s, one root 49 s)
    two roots   both runs end complete and exit 0 (B3 and B6 in each: 15 gating and 8 reported checks, all passing); the first took the slot at once, the second waited (4 or 5
                waiting lines each, the first one naming the other's pid, e.g. "held by 'bench B3 FEMM reference' (pid N)"); 309 to 310 processes were below B9 over the whole run, both pairs
                (python.exe, femm_headlessw.exe, fkn.exe, ElmerSolver.exe, ElmerGrid.exe), none left, the slot free (2,420 looks at the process list, none failed; the slot looked at 597 times and the
                root's lock 192, none unreadable); every process seen below normal in 4,633 looks; at most 8 femm_headlessw.exe alive at one look; 30 records in the runs' runs.json (both pairs),
                at most 8 workers and 4 threads; 2,037 numbers compared (B3: 1,787, B6: 250), all identical in the second and the third real run; in the first and the fourth, B6's high modes
                (critical speeds of 450,000 rpm) differed by up to 4.9e-8 and 4.8e-8 between run 1 and run 2 (the fourth: 68 of B6's 250 numbers, the highest mode 451,957.695 against
                451,957.717 rpm, the two values the eight solo runs at threads=2 gave too; threads=1 gave 451,957.706 and 451,957.686), and 16 solo runs of B6's ROSS step (8 with
                threads=2, 8 with threads=1) took two values in each group, up to 4.9e-8 and 4.3e-8 apart: B3 is compared exactly (identical in all four), B6 to 1e-6
                (TOLERANCE_BY_CASE says what is known)
    one root    the second printed "waiting for the bench root ...: held by 'bench run B6' (pid N)", started its work 0.1 to 0.2 s after the first's last file, finished complete (24 s each); the lock
                was free afterwards
    the picture a Gantt chart of who held the slot (two runs taking turns, B3's FEMM reference, its Elmer runs, B6's ROSS) and the root's lock (two bars in a row)

Self-check (the answers follow).
  1. Why does B9 not take the solver slot itself, to be safe, while it starts the two runs?
  2. The second run's output says "held by 'bench B3 FEMM' (pid 4001)", and the first child B9 started has pid 3990. Is that a failure?
  3. Why are the seconds of a run left out when the two runs' numbers are compared?
Answers.
  1. The children need the slot. B9 holding it would make both wait for it until B9's timeout, which is a deadlock that the lock was built to make impossible: a process that waits for
     the slot while holding it for itself is refused (RuntimeError), and a process that holds it while waiting for its children's solvers waits for ever. B9 only watches.
  2. No. The process B9 started is the launcher in front of the interpreter (the virtual environment's python.exe); the holder is the interpreter's pid, which the child printed on
     its first line. The check compares the pid in the waiting line with the one the run printed, and, to end a child, requires that pid to be below the one B9 started, and the
     process on it to have the creation time the tracker recorded (read through the one handle that then ends it).
  3. They are the run's own: the second run waited for the first, so it took longer by exactly the wait. The numbers that describe the MOTOR (flux linkage, stress, nodes, indices) are what
     two runs of deterministic tools must repeat.
"""
from __future__ import annotations

import colorsys
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import femm_queue
import paths
from checks.matlab_check import JobTracker, Watcher, process_table

from . import OWN_LIMIT, HolderSampler, Log, check, root_lock_dir, watch_children, wrap_title

CASE = "B9"
TITLE = "The rules themselves: two bench runs started at once, on two roots (and two on one root), take turns on the solver slot"
CHILD_CASES = ("B3", "B6")                                  # the light cases both children run (B2 exits 1 by a standing ruling; B9 itself is never a child's case)
SAME_ROOT_CASES = ("B6",)                                   # the case of the two runs on one root
MIN_CHECKS = {"B3": 15, "B6": 8}                            # the checks each case has (the standalone runs of 11-06: out\bench\B3 and B6 results.json): fewer means a damaged run
PAIR_TIMEOUT_S = 900                                        # the two runs on two roots: 2 x 50 s of slot time, and the wait for each other
SAME_TIMEOUT_S = 400                                        # the two runs on one root: 2 x 25 s
LOCK_WAIT_S = 120                                           # how long B9 waits for the first run to hold its root before it starts the second
POLL_S = 0.25                                               # the holder is looked at this often
GRACE_S = 2.0                                               # after the children have ended: what they started may still be ending
# Two runs' numbers agree to this share of themselves, case by case. B3 (FEMM, Elmer, Gmsh): exactly, as they were in every real pair so far. B6 (ROSS): 1e-6, a limit SET AFTER THE
# FIRST REAL RUN (it was 1e-9 before). What is known: once, in one of the first concurrent runs, ROSS's modes near 450,000 rpm (critical speeds of 446,165, 447,547 and 451,958 rpm,
# no meaning for a motor) differed from a solo run's by up to 4.9e-8; two other concurrent pairs were identical, and in a later one B6 differed again by up to 4.8e-8, in the same
# modes and with the same two values the eight solo runs at threads=2 gave (451,957.695 and 451,957.717 rpm for the highest; with threads=1 the solo runs gave 451,957.706 and
# 451,957.686); and 16 solo runs of B6's ROSS step on 2026-11-06 (8 with threads=2, B6's own setting, and 8 with threads=1, one after another through the slot) took two values in
# each group, up to 4.9e-8 and 4.3e-8 apart. Setting OMP_NUM_THREADS=1 (run_solver's threads=1) did not remove the variation; the count was set through the environment and not
# read back from ROSS's BLAS, so that is all it shows. ROSS fixes its ARPACK start vector (ROSS 3.1.0, ross/rotor_assembly.py:2422: v0=np.ones(A.shape[0]), opened 2026-11-06), so
# a random start is not the cause either; the cause is not isolated. 1e-6 is 20 times above the largest difference seen and 100 times below the tightest physical tolerance of the
# bench (0.01 %). Only numbers are compared, not labels (check 3 says why).
TOLERANCE_BY_CASE = {"B3": 0.0, "B6": 1e-6}
NOT_STARTED = "not started"                                 # what a child's exit code says when it could not be started
LEGEND_COLUMNS = 4                                          # the picture's legend: columns of kinds of work, and (in inches) the height of one row of it and its padding,
LEGEND_ROW_IN, LEGEND_PAD_IN = 0.2, 0.12                    # the room under the lower panel for its tick labels and axis label, and the height of the two panels with the space between
BELOW_PANELS_IN, PANELS_IN = 0.75, 4.1
SLACK_S = 1.0                                               # file times and the clock are not exact to better than this
FIGURE_MIN_BYTES = 20_000                                   # a case's figure.png is a real picture
PICTURE_MIN_BYTES = 1_000                                   # any other picture of a case, at least
ONLY_RUN_BY_NAME = ("run1", "run2", "same")                 # the roots B9 makes under its own folder


# ============================================================================================== starting and ending the children
def child_command(root_folder, cases) -> list:
    """`python bench\\run_bench.py --out ROOT --only CASES`: a bench run on a root of its own. The interpreter is this one (the virtual environment's)."""
    return [sys.executable, paths.P / "bench" / "run_bench.py", "--out", Path(root_folder), "--only", ",".join(cases)]


def start_child(command, out_path, err_path) -> subprocess.Popen:
    """Start a bench run (the one kind of program B9 starts, named in tests/test_bench.py): below normal, no console input, its output in two files. Returns the Popen."""
    with open(out_path, "wb") as out, open(err_path, "wb") as err:
        return subprocess.Popen([str(part) for part in command], cwd=paths.P, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                                creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS)


def is_below(pid: int, ancestor: int) -> bool:
    """Is process `pid` a descendant of process `ancestor`, by the parent links of the process list now? (Both must be in the list.)"""
    table = process_table()
    current = pid
    for _ in range(16):
        if current not in table:
            return False
        current = table[current][0]
        if current == ancestor:
            return True
    return False


PROCESS_QUERY_LIMITED_INFORMATION = 0x1000                  # the access end_child asks for with its one handle: to read when the process was created ...
PROCESS_TERMINATE = 0x0001                                  # ... and to end it


def recorded_creation(tracker, pid: int) -> int | None:
    """The creation time (a Windows FILETIME) the JobTracker recorded for process `pid`: that of its newest record; None if it has none or could not read one."""
    records = [m for m in tracker.members() if m.pid == pid and m.created is not None]
    return max(records, key=lambda m: m.first_seen).created if records else None


def end_child(child, real_pid: int | None = None, real_created: int | None = None) -> None:
    """End a child B9 started (the only place in bench/ that ends a process, named in tests/test_bench.py), on its timeout: the child itself, by the handle Popen holds, and the
    interpreter below it that ran its code, if `real_pid` is really below this child. That one is ended through ONE handle opened to its pid: the creation time read through the
    handle has to be `real_created`, the one the tracker recorded for that process, and only then is it ended, through that same handle. A pid that Windows has given to another
    process since, and a process that has no record, are left alone, and there is no check-then-kill by pid. Whatever the interpreter had started ends by itself when it dies: its
    pool workers watch their owner and run_solver's job closes with its caller."""
    if real_pid is not None and real_created is not None and real_pid != child.pid and is_below(real_pid, child.pid):
        import ctypes
        from ctypes import wintypes
        k = ctypes.WinDLL("kernel32", use_last_error=True)
        k.OpenProcess.argtypes, k.OpenProcess.restype = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD], wintypes.HANDLE
        k.GetProcessTimes.argtypes, k.GetProcessTimes.restype = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4, wintypes.BOOL
        k.CloseHandle.argtypes, k.CloseHandle.restype = [wintypes.HANDLE], wintypes.BOOL
        handle = k.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION | PROCESS_TERMINATE, False, real_pid)
        if handle:
            try:
                created, exited, kernel, user = (wintypes.FILETIME() for _ in range(4))
                if k.GetProcessTimes(handle, ctypes.byref(created), ctypes.byref(exited), ctypes.byref(kernel), ctypes.byref(user)) \
                        and (created.dwHighDateTime << 32 | created.dwLowDateTime) == real_created:
                    k.TerminateProcess(ctypes.c_void_p(handle), 1)
            finally:
                k.CloseHandle(handle)
    child.kill()
    child.wait(30)


def wait_for_children(children, timeout_s: float, out_paths=None, on_exit=None, created_of=None) -> list[int | None]:
    """Wait for every child, at most `timeout_s` seconds. A child still running then is ended (end_child; its real pid is read from its output file in `out_paths`, and its creation
    time is `created_of(pid)`, the tracker's record of it) and its code is None; a child that finished has its exit code. `on_exit(index)` is called when a child is first seen to
    have finished (not for one that was ended)."""
    seen: set[int] = set()
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        for index, child in enumerate(children):
            if index not in seen and child.poll() is not None:
                seen.add(index)
                if on_exit is not None:
                    on_exit(index)
        if len(seen) == len(children):
            break
        time.sleep(0.2)
    codes: list[int | None] = []
    for index, child in enumerate(children):
        if child.poll() is None:
            pid = run_pid(Path(out_paths[index]).read_text(encoding="utf-8", errors="replace")) if out_paths else None
            end_child(child, pid, created_of(pid) if created_of is not None and pid is not None else None)
            codes.append(None)
        else:
            if index not in seen and on_exit is not None:
                on_exit(index)
            codes.append(child.returncode)
    return codes


# ============================================================================================== what a run said
def run_pid(text: str) -> int | None:
    """The pid a bench run printed on its first line ('bench: root <ROOT> (pid N)'); None if there is no such line."""
    found = re.search(r"^bench: root .* \(pid (\d+)\)\s*$", text, re.M)
    return int(found.group(1)) if found else None


SLOT_LINE = re.compile(r"^waiting for the solver slot: held by (?:'(?P<label>.*)' \(pid (?P<pid>\d+)\) since (?P<since>\S+)|a run whose note.*)$")
ROOT_LINE = re.compile(r"^waiting for the bench root (?P<root>.*): held by (?:'(?P<label>.*)' \(pid (?P<pid>\d+)\) since (?P<since>\S+)|a run whose note.*)$")


def waiting_lines(text: str) -> list[dict]:
    """Every line in which a run said it was waiting, in order: {"kind": "slot" or "root", "label", "pid", "since", "root", "line"} (the holder's label, pid and time are None if its
    note could not be read)."""
    out = []
    for line in text.splitlines():
        line = line.rstrip()
        slot, root = SLOT_LINE.match(line), ROOT_LINE.match(line)
        match = slot or root
        if match:
            out.append({"kind": "slot" if slot else "root", "label": match.group("label"), "pid": int(match.group("pid")) if match.group("pid") else None, "since": match.group("since"),
                        "root": root.group("root") if root else None, "line": line})
    return out


def final_line(text: str) -> str | None:
    """The line with which a bench run ended ('bench: PASS, ...' or 'bench: FAIL, ...'); None if it printed none (it did not finish)."""
    lines = [line.strip() for line in text.splitlines() if line.startswith(("bench: PASS", "bench: FAIL"))]
    return lines[-1] if lines else None


# ============================================================================================== judging the runs
def tree_names(folder) -> list[str]:
    """The files under `folder`, relative to it, sorted."""
    folder = Path(folder)
    return sorted(str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file())


def files_differ(first, second) -> list[str]:
    """How two case folders differ in the files they hold: '<name>: only in the first/second'. Empty if they hold the same."""
    a, b = set(tree_names(first)), set(tree_names(second))
    return [f"{name}: only in the first" for name in sorted(a - b)] + [f"{name}: only in the second" for name in sorted(b - a)]


def case_problems(folder, min_checks: int) -> list[str]:
    """What is missing or wrong in a case's folder after its run: results.json (with at least `min_checks` checks, none that says the case did not finish), figure.png (a real
    picture), every other picture, run.log and runs.json. Empty if the run is complete."""
    folder = Path(folder)
    if not folder.is_dir():
        return [f"the case folder {folder} is not there"]
    problems: list[str] = []
    results = folder / "results.json"
    if not results.is_file():
        problems.append("results.json is not there")
    else:
        try:
            data = json.loads(results.read_text(encoding="utf-8"))
        except ValueError:
            problems.append("results.json is not JSON")
        else:
            checks = data.get("checks") or []
            if "error" in data or any(c.get("name") == "(did not finish)" for c in checks):
                problems.append(f"the case did not finish: {data.get('error', 'a check says so')}")
            if len(checks) < min_checks:
                problems.append(f"{len(checks)} checks, the case has at least {min_checks}")
    for name in ("run.log", "runs.json"):
        if not (folder / name).is_file():
            problems.append(f"{name} is not there")
    figure = folder / "figure.png"
    if not figure.is_file():
        problems.append("figure.png is not there")
    elif figure.stat().st_size < FIGURE_MIN_BYTES:
        problems.append(f"figure.png is {figure.stat().st_size} bytes: not a picture")
    for png in sorted(folder.rglob("*.png")):
        if png != figure and png.stat().st_size < PICTURE_MIN_BYTES:
            problems.append(f"{png.relative_to(folder)} is {png.stat().st_size} bytes: not a picture")
    return problems


def _ignored(key: str) -> bool:
    """A key whose number belongs to THIS run and differs by nature: seconds of any kind, the time of day, process ids, looks at the process list."""
    k = key.lower()
    return "seconds" in k or k in ("time", "timestamp", "elapsed", "pid", "pids", "worker_pids", "femm_pids") or k.endswith(("_pid", "_pids")) or "looks" in k


def flatten_numbers(data, prefix: str = "") -> dict[str, float]:
    """{path: value} of every number in a JSON document (not the booleans, not the text), except the keys that belong to one run (seconds, pids, looks, the time of day)."""
    out: dict[str, float] = {}
    if isinstance(data, dict):
        for key, value in data.items():
            if not _ignored(str(key)):
                out.update(flatten_numbers(value, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(data, (list, tuple)):
        for i, value in enumerate(data):
            out.update(flatten_numbers(value, f"{prefix}[{i}]"))
    elif isinstance(data, bool):
        pass
    elif isinstance(data, (int, float)):
        out[prefix] = float(data)
    return out


def compare_runs(first, second, tolerance: float = 0.0) -> dict:
    """Two runs' results.json against each other: {"compared", "exact", "worst_relative", "different": [(path, first, second)], "only_in_first", "only_in_second"}. A number differs
    when it is further apart than `tolerance` of the larger of the two (0.0, the default: any difference at all). A number that is not finite is a difference whatever the tolerance,
    unless both are NaN or both are the same infinity; its worst_relative is infinity."""
    a, b = flatten_numbers(first), flatten_numbers(second)
    common = sorted(set(a) & set(b))
    different, exact, worst = [], 0, 0.0
    for key in common:
        x, y = a[key], b[key]
        if x == y or (math.isnan(x) and math.isnan(y)):
            exact += 1
            continue
        if not (math.isfinite(x) and math.isfinite(y)):                      # a NaN or an infinity on one side only, or two different infinities
            different.append((key, x, y))
            worst = math.inf
            continue
        relative = abs(x - y) / max(abs(x), abs(y))                         # x != y, so the larger is not 0
        worst = max(worst, relative)
        if relative > tolerance:
            different.append((key, x, y))
    return {"compared": len(common), "exact": exact, "worst_relative": worst, "different": different, "only_in_first": sorted(set(a) - set(b)), "only_in_second": sorted(set(b) - set(a))}


def runs_comply(records: list[dict], limit: int) -> dict:
    """The rules, read from one case's runs.json: every record below normal and at most `limit` workers or threads. {"ok", "problems", "pool_records", "solver_records",
    "in_process_records", "observed", "max_workers", "max_threads"}. A pool record has to say what it saw of its workers' and FEMMs' classes; a program too short to be seen is noted
    (not a failure), unless no program of the case was seen at all; a class that is not below normal, or a count above the limit, is a problem."""
    problems: list[str] = []
    counts = {"femm_pool": 0, "run_solver": 0, "in_process": 0}
    observed = unseen = 0
    max_workers = max_threads = 0
    for record in records:
        kind, label = record.get("kind"), record.get("label", "?")
        if kind == "femm_pool":
            counts[kind] += 1
            workers, solved = record.get("workers"), record.get("workers_that_solved")
            max_workers = max(max_workers, workers if isinstance(workers, int) else 0)
            if not isinstance(workers, int) or workers > limit:
                problems.append(f"{label}: {workers} workers, at most {limit}")
            if isinstance(solved, int) and solved > limit:
                problems.append(f"{label}: {solved} workers solved, at most {limit}")
            for key in ("worker_priority_classes_seen", "femm_priority_classes_seen"):
                classes = record.get(key) or {}
                if not classes:
                    problems.append(f"{label}: {key} is empty: nothing was seen")
                elif set(classes) != {"below_normal"}:
                    problems.append(f"{label}: {key} {classes}")
            if record.get("worker_priority_classes_seen"):
                observed += 1
        elif kind == "run_solver":
            counts[kind] += 1
            threads = record.get("threads")
            max_threads = max(max_threads, threads if isinstance(threads, int) else 0)
            if not isinstance(threads, int) or not 1 <= threads <= limit:
                problems.append(f"{label}: {threads} threads, between 1 and {limit}")
            classes = record.get("priority_classes_seen") or {}
            seen = {c for values in classes.values() for c in values}
            if seen - {"below_normal"}:
                problems.append(f"{label}: seen at {sorted(seen)}")
            if classes:
                observed += 1
            else:
                unseen += 1
        elif kind == "in_process":
            counts[kind] += 1
            threads = record.get("threads")
            max_threads = max(max_threads, threads if isinstance(threads, int) else 0)
            if record.get("priority_class_seen") != "below_normal":
                problems.append(f"{label}: priority {record.get('priority_class_seen')}")
            else:
                observed += 1
            if threads is not None and not (isinstance(threads, int) and 1 <= threads <= limit):
                problems.append(f"{label}: {threads} threads, between 1 and {limit}")
        else:
            problems.append(f"a record of kind {kind!r} that B9 does not know")
    if counts["run_solver"] and unseen == counts["run_solver"]:
        problems.append("no program of this case was seen while it ran: the log shows nothing about its priority")
    return {"ok": not problems, "problems": problems, "pool_records": counts["femm_pool"], "solver_records": counts["run_solver"], "in_process_records": counts["in_process"],
            "observed": observed, "max_workers": max_workers, "max_threads": max_threads}


def serialized(first_ended: float, second_first_work: float | None, slack: float = SLACK_S) -> bool:
    """Two runs on one root took turns: the second's first file is no older than the first's last (within `slack` seconds). No file of the second: it cannot be said."""
    return second_first_work is not None and second_first_work >= first_ended - slack


def holder_segments(samples) -> list[dict]:
    """Fold a HolderSampler's samples [(t, holder)] into stretches of one holder: [{"holder", "start", "end"}], a stretch ending at the first look that saw another (the last one at
    its own last look)."""
    segments: list[dict] = []
    for t, holder in samples:
        if segments and segments[-1]["holder"] == holder:
            segments[-1]["end"] = t
            continue
        if segments:
            segments[-1]["end"] = t
        segments.append({"holder": holder, "start": t, "end": t})
    return segments


def holders_seen(samples) -> list[tuple]:
    """The holders in a HolderSampler's samples, each once, in the order first seen; a free slot is not a holder."""
    out: list[tuple] = []
    for _, holder in samples:
        if holder is not None and holder not in out:
            out.append(holder)
    return out


# ============================================================================================== the picture
def figure_title() -> str:
    """The title of figure.png. Every line is at most bench.TITLE_CHARS long, so that wrap_title leaves it as written (the second line was 104 characters and left an orphan word)."""
    return ("B9: two bench runs started at once take turns\n"
            "above: two roots, one solver slot (B3, then B6, in each run)\n"
            "below: one root, one lock (B6 in each run)\n"
            "a bar is a stretch in which a run held it; the colour says what the run was doing")


def label_colours(labels) -> dict[str, str]:
    """A colour of its own for every label, as many as there are: the old palette had six colours and reused them from the seventh kind of work on, so that two kinds shared one.
    The hues are spaced evenly round the circle and neighbours alternate between a strong dark shade and a soft lighter one, so that any two are far apart (50 or more in RGB for up to
    15 labels). The same labels give the same colours. {label: '#rrggbb'} in the order of the labels (each once)."""
    labels = list(dict.fromkeys(labels))
    n = max(len(labels), 1)
    out = {}
    for i, label in enumerate(labels):
        r, g, b = colorsys.hsv_to_rgb(i / n, (0.9, 0.5)[i % 2], (0.58, 0.72)[i % 2])
        out[label] = "#%02x%02x%02x" % (round(255 * r), round(255 * g), round(255 * b))
    return out


def draw_timeline(slot_samples, root_samples, names: dict, title: str):
    """Who held what, over time, as a matplotlib figure (not saved): the solver slot while the two runs on two roots took turns (above), and the lock of the one root while the two runs
    on it took turns (below). One row for each run (named by `names`: {pid: text}), a bar for each stretch it held the lock, labelled with what it was doing; a free lock is a gap.
    Every kind of work (a label) has a colour of its own, the same in both panels."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = ((None, slot_samples, "the solver slot (two runs, two roots)", "seconds since the two runs were started"),
              (None, root_samples, "the lock of the root (two runs, one root)", "seconds since the first of the two was started"))
    kinds = [s["holder"][0] for _, samples, _, _ in panels for s in holder_segments(samples) if s["holder"] is not None]
    colours = label_colours(kinds)                                                  # what a run was doing -> its colour
    # The room under the panels is made to fit: the tick labels and the axis label of the lower panel (BELOW_PANELS_IN), then the legend, as high as its rows are (ROW_IN each).
    # A margin that was a fixed share of the figure held two rows of legend and the real run has four, so the first row crossed the axis label.
    text = wrap_title(title)
    lines = text.count("\n") + 1
    legend_in = (LEGEND_ROW_IN * -(-len(colours) // LEGEND_COLUMNS) + LEGEND_PAD_IN) if colours else 0.0
    top_in, bottom_in = 0.2 * lines + 0.45, BELOW_PANELS_IN + legend_in
    height = top_in + PANELS_IN + bottom_in
    fig, axes = plt.subplots(2, 1, figsize=(11, height))
    fig.subplots_adjust(left=0.1, right=0.985, bottom=bottom_in / height, top=1 - top_in / height, hspace=0.62)
    fig.text(0.005, 1 - 0.05 / height, text, ha="left", va="top", fontsize=10.5, color="#1d1d1f")
    panels = tuple((ax, samples, heading, x_label) for ax, (_, samples, heading, x_label) in zip(axes, panels))
    for ax, samples, heading, x_label in panels:
        segments = [s for s in holder_segments(samples) if s["holder"] is not None]
        pids = []
        for s in segments:
            if s["holder"][1] not in pids:
                pids.append(s["holder"][1])
        for s in segments:
            row = pids.index(s["holder"][1])
            label = s["holder"][0]
            ax.broken_barh([(s["start"], max(s["end"] - s["start"], 0.15))], (row - 0.36, 0.72), facecolors=colours[label])
            if s["end"] - s["start"] > 6:
                ax.text((s["start"] + s["end"]) / 2, row, label.replace("bench ", ""), ha="center", va="center", fontsize=8.5, color="white")
        ax.set_yticks(range(len(pids)))
        ax.set_yticklabels([names.get(p, f"pid {p}") for p in pids])
        ax.set_ylim(-0.7, max(len(pids) - 0.3, 0.7))
        ax.invert_yaxis()
        ax.set_xlabel(x_label)
        ax.set_title(heading, fontsize=9.5, loc="left")
        ax.grid(axis="x", alpha=0.3)
        if not segments:
            ax.text(0.5, 0.5, "nobody was seen holding it", transform=ax.transAxes, ha="center", va="center")
    if colours:
        from matplotlib.patches import Patch
        fig.legend(handles=[Patch(color=colour, label=label.replace("bench ", "")) for label, colour in colours.items()], loc="lower center", ncol=LEGEND_COLUMNS, fontsize=8,
                   frameon=False, bbox_to_anchor=(0.5, 0.0))
    return fig


def fig_timeline(slot_samples, root_samples, names: dict, title: str, target) -> None:
    """draw_timeline, saved with a tight box."""
    import matplotlib.pyplot as plt
    fig = draw_timeline(slot_samples, root_samples, names, title)
    fig.savefig(target, dpi=130, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


# ============================================================================================== the case
def read_text(path) -> str:
    try:
        return Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def run(folder: Path, log: Log) -> dict:
    """Run B9 in `folder` (already empty). Returns the results dict, whose "checks" list every verdict.

    The watchers (the tracker of every process below B9, the looks at the slot, the classes of every process) are started once, before the first child, and stopped once, after the
    last child has ended and its grace: the pair on one root is looked at as the pair on two roots is, and what is left alive, and whether the slot is free, are read once, at the end."""
    folder = Path(folder)
    started = time.perf_counter()
    limit = paths.setting("femm", "workers")
    log.say(f"B9: two runs on two roots (run1, run2: {', '.join(CHILD_CASES)}), then two on one root (same: {', '.join(SAME_ROOT_CASES)}); the solver slot is held by "
            f"{femm_queue.holder() or 'nobody'} now")
    roots = [folder / "run1", folder / "run2"]
    outs, errs = [r / "child.out.txt" for r in roots], [r / "child.err.txt" for r in roots]
    for r in roots:
        r.mkdir(parents=True)
    same = folder / "same"
    same.mkdir()
    c_out, c_err, d_out, d_err = same / "c.out.txt", same / "c.err.txt", same / "d.out.txt", same / "d.err.txt"
    lock = root_lock_dir(same)
    zero = time.perf_counter()
    tracker = JobTracker(os.getpid(), process_table())
    watcher, sampler = Watcher(tracker, zero), HolderSampler(zero, POLL_S)
    root_sampler = HolderSampler(zero, POLL_S, lock_dir=lock)
    problems: list[str] = []                                                       # what is wrong with the completeness of the runs, the steps that failed first
    ended: dict = {}
    watcher.start()
    sampler.start()
    with watch_children() as classes:
        try:
            # ---- the pair on two roots, started at once; a child that cannot be started does not end the one that was: that one is waited for, with the same timeout
            pair_started = time.perf_counter()
            children = []
            for number, (r, o, e) in enumerate(zip(roots, outs, errs), start=1):
                try:
                    children.append(start_child(child_command(r, CHILD_CASES), o, e))
                except (OSError, ValueError) as error:
                    problems.append(f"the pair on two roots: starting run {number} failed ({type(error).__name__}: {error}); the run started before it was waited for, not ended")
                    break
            codes_a = wait_for_children(children, PAIR_TIMEOUT_S, out_paths=outs, created_of=lambda pid: recorded_creation(tracker, pid)) + [NOT_STARTED] * (len(roots) - len(children))
            time.sleep(GRACE_S)
            seconds_a = time.perf_counter() - pair_started
            pair_a_ended = time.perf_counter() - zero                              # the slot's samples after this belong to the second pair
            # ---- the pair on one root: the first is started, B9 waits until the root's lock note names it, and starts the second
            root_sampler.zero = time.perf_counter()
            root_sampler.start()
            child_c = child_d = None
            problems_b: list[str] = []
            try:
                child_c = start_child(child_command(same, SAME_ROOT_CASES), c_out, c_err)
            except (OSError, ValueError) as error:
                problems_b.append(f"the pair on one root: starting the first run failed ({type(error).__name__}: {error})")
            if child_c is not None:
                deadline = time.monotonic() + LOCK_WAIT_S
                while time.monotonic() < deadline and child_c.poll() is None:
                    c_pid = run_pid(read_text(c_out))
                    note = read_json(lock / "holder.json")
                    if c_pid is not None and isinstance(note, dict) and note.get("pid") == c_pid:
                        break
                    time.sleep(0.1)
                try:
                    child_d = start_child(child_command(same, SAME_ROOT_CASES), d_out, d_err)
                except (OSError, ValueError) as error:
                    problems_b.append(f"the pair on one root: starting the second run failed ({type(error).__name__}: {error}); the first was waited for, not ended")

            def first_ended(index: int) -> None:
                if index == 0:                                                      # C ended: its summary.json is the last file it wrote, and D cannot have written its own yet
                    summary = same / "summary.json"
                    ended["c"] = summary.stat().st_mtime if summary.is_file() else None
            began = [(c, out) for c, out in ((child_c, c_out), (child_d, d_out)) if c is not None]
            codes_b = wait_for_children([c for c, _ in began], SAME_TIMEOUT_S, out_paths=[out for _, out in began], on_exit=first_ended,
                                        created_of=lambda pid: recorded_creation(tracker, pid)) + [NOT_STARTED] * (2 - len(began))
            root_sampler.stop()
            time.sleep(GRACE_S)
        finally:
            if root_sampler.is_alive():
                root_sampler.stop()
            watcher.stop()
            sampler.stop()
    # ---- what is left, read once, now that both pairs have ended: no process of the whole tree, the slot, the root's lock
    alive = [m.name for m in tracker.alive(process_table())]
    members = tracker.members()
    slot_free = femm_queue.holder() is None
    lock_free = femm_queue.holder(lock) is None
    slot_a = [sample for sample in sampler.samples if sample[0] <= pair_a_ended]
    texts_a = [read_text(o) for o in outs]
    pids_a = [run_pid(t) for t in texts_a]
    waits_a = [waiting_lines(t) for t in texts_a]
    log.say(f"B9: pair on two roots done in {seconds_a:.0f} s: exit codes {codes_a}, pids {pids_a}, waiting lines {[len(w) for w in waits_a]}")
    # ---- completeness, numbers and rules, case by case (the pair on two roots)
    numbers, rules_problems = {"compared": 0, "exact": 0, "worst_relative": 0.0, "different": [], "only_in_first": [], "only_in_second": [], "per_case": {}}, []
    record_total = observed = max_workers = max_threads = 0
    per_case = {}
    for case in CHILD_CASES:
        folders = [r / case for r in roots]
        for r, f in zip(roots, folders):
            problems += [f"{r.name}\\{case}: {p}" for p in case_problems(f, MIN_CHECKS[case])]
        if all(f.is_dir() for f in folders):
            problems += [f"{case}: {d}" for d in files_differ(*folders)]
            a, b = read_json(folders[0] / "results.json"), read_json(folders[1] / "results.json")
            if a is not None and b is not None:
                compared = compare_runs(a, b, TOLERANCE_BY_CASE[case])
                per_case[case] = {"compared": compared["compared"], "exact": compared["exact"], "worst_relative": compared["worst_relative"], "different": len(compared["different"]),
                                  "checks": [len(a.get("checks", [])), len(b.get("checks", []))], "tolerance": TOLERANCE_BY_CASE[case]}
                numbers["per_case"][case] = {"compared": compared["compared"], "exact": compared["exact"], "worst_relative": compared["worst_relative"], "tolerance": TOLERANCE_BY_CASE[case]}
                for key in ("compared", "exact"):
                    numbers[key] += compared[key]
                numbers["worst_relative"] = max(numbers["worst_relative"], compared["worst_relative"])
                numbers["different"] += [(f"{case} {k}", x, y) for k, x, y in compared["different"]]
                numbers["only_in_first"] += [f"{case} {k}" for k in compared["only_in_first"]]
                numbers["only_in_second"] += [f"{case} {k}" for k in compared["only_in_second"]]
        for r, f in zip(roots, folders):
            runs = read_json(f / "runs.json") or []
            judged = runs_comply(runs, limit)
            rules_problems += [f"{r.name}\\{case}: {p}" for p in judged["problems"]]
            record_total += len(runs)
            observed += judged["observed"]
            max_workers, max_threads = max(max_workers, judged["max_workers"]), max(max_threads, judged["max_threads"])
    for r in roots:
        summary = read_json(r / "summary.json")
        if summary is None or summary.get("cases") != list(CHILD_CASES):
            problems.append(f"{r.name}: summary.json is not there or does not list {list(CHILD_CASES)}")
    # ---- the pair on one root: its log is judged by the same rules (the second run's runs.json is what is left in the root: it emptied the first's folder)
    texts_b = [read_text(c_out), read_text(d_out)]
    pids_b = [run_pid(t) for t in texts_b]
    d_files = [p for p in (same / SAME_ROOT_CASES[0]).rglob("*") if p.is_file()] if (same / SAME_ROOT_CASES[0]).is_dir() else []
    d_first = min((p.stat().st_mtime for p in d_files), default=None)
    d_problems = problems_b + case_problems(same / SAME_ROOT_CASES[0], MIN_CHECKS[SAME_ROOT_CASES[0]])
    for case in SAME_ROOT_CASES:
        runs = read_json(same / case / "runs.json") or []
        judged = runs_comply(runs, limit)
        rules_problems += [f"same\\{case}: {p}" for p in judged["problems"]]
        record_total += len(runs)
        observed += judged["observed"]
        max_workers, max_threads = max(max_workers, judged["max_workers"]), max(max_threads, judged["max_threads"])
    log.say(f"B9: pair on one root done: exit codes {codes_b}, pids {pids_b}, the second's waiting lines {[w['line'] for w in waiting_lines(texts_b[1])][:1]}; "
            f"{len(members)} processes were below B9 over the whole run, {len(alive)} left")
    # ---- the picture, the numbers, the verdicts
    names = {pids_a[0]: "run 1", pids_a[1]: "run 2", pids_b[0]: "first run (same root)", pids_b[1]: "second run (same root)"}
    names = {pid: text for pid, text in names.items() if pid is not None}
    fig_timeline(slot_a, root_sampler.samples, names, figure_title(), folder / "figure.png")
    peaks = classes.peaks()
    classes_seen = classes.summary()
    context = {
        "limit": limit,
        "a": {"pids": pids_a, "returncodes": codes_a, "timed_out": [c is None for c in codes_a], "seconds": seconds_a, "waits": waits_a, "final_lines": [final_line(t) for t in texts_a]},
        "complete": {"problems": problems, "cases_checked": len(CHILD_CASES) * 2, "checks_per_case": {c: per_case.get(c, {}).get("checks") for c in CHILD_CASES}},
        "numbers": numbers,
        "left": {"alive_after": alive, "members": len(members), "failed_looks": watcher.failed, "looks": watcher.looks, "slot_free": slot_free,
                 "names": sorted({m.name.lower() for m in members}), "slot_looks": len(sampler.samples), "slot_unreadable": sampler.errors, "lock_looks": len(root_sampler.samples),
                 "lock_unreadable": root_sampler.errors},
        "rules": {"problems": rules_problems, "records": record_total, "observed": observed, "max_workers": max_workers, "max_threads": max_threads, "classes_seen": classes_seen,
                  "peaks": peaks, "looks": classes.looks},
        "same": {"pids": pids_b, "returncodes": codes_b, "timed_out": [c is None for c in codes_b], "c_ended_at": ended.get("c"), "d_first_work_at": d_first,
                 "waits": [w for w in waiting_lines(texts_b[1]) if w["kind"] == "root"], "d_problems": d_problems, "c_final": final_line(texts_b[0]), "d_final": final_line(texts_b[1]),
                 "lock_free": lock_free},
    }
    results = {"case": CASE, "title": TITLE, "cases": {"two_roots": list(CHILD_CASES), "one_root": list(SAME_ROOT_CASES)},
               "pair_a": {"roots": [str(r) for r in roots], "pids": pids_a, "returncodes": codes_a, "seconds": round(seconds_a, 1), "waiting_lines": [[w["line"] for w in ws] for ws in waits_a],
                          "final_lines": context["a"]["final_lines"], "slot_holders": [list(h) for h in holders_seen(slot_a)], "slot_looks": len(slot_a),
                          "per_case_numbers": per_case, "problems": problems},
               "processes_below_b9": {"members": [(m.name, round(m.first_seen, 1), round(m.last_seen, 1)) for m in members][:400], "count": len(members), "alive_after": alive,
                                      "names": context["left"]["names"], "peaks": peaks, "classes_seen": classes_seen, "failed_looks": watcher.failed},
               "looks": {"slot": len(sampler.samples), "slot_unreadable": sampler.errors, "root_lock": len(root_sampler.samples), "root_lock_unreadable": root_sampler.errors},
               "rules": {"problems": rules_problems, "records": record_total, "observed": observed, "max_workers": max_workers, "max_threads": max_threads, "limit": limit},
               "pair_b": {"root": str(same), "pids": pids_b, "returncodes": codes_b, "waiting_lines": [w["line"] for w in waiting_lines(texts_b[1])], "final_lines": [context["same"]["c_final"],
                          context["same"]["d_final"]], "c_summary_written_at": ended.get("c"), "d_first_file_at": d_first, "lock_holders": [list(h) for h in holders_seen(root_sampler.samples)],
                          "d_problems": d_problems},
               "numbers": {k: v for k, v in numbers.items() if k != "different"} | {"different": [list(d) for d in numbers["different"]][:20], "tolerance_by_case": TOLERANCE_BY_CASE}}
    results["checks"] = make_checks(context)
    results["seconds_total"] = round(time.perf_counter() - started, 1)
    return results


# ============================================================================================== the checks
def make_checks(context: dict) -> list[dict]:
    """The six checks of the ruling, in order, and what is reported beside them (the shape of `context` is good_context's in tests/test_bench_b9.py)."""
    out: list[dict] = []
    limit, a, same = context["limit"], context["a"], context["same"]
    # 1. the second run printed the slot's waiting line, naming the holder
    pids = a["pids"]
    named = [(i, w) for i, waits in enumerate(a["waits"]) for w in waits if w["kind"] == "slot" and w["pid"] is not None and w["pid"] in [p for j, p in enumerate(pids) if j != i and p is not None]]
    out.append(check("the second run printed the slot's waiting line, naming the holder", bool(named),
                     "the run that did not get the solver slot first prints 'waiting for the solver slot: held by '<label>' (pid N) since HH:MM', and N is the other run's pid (the brief: 'The "
                     "second waits for the slot and prints the holder')",
                     (f"run {named[0][0] + 1} said: {named[0][1]['line']} (the other run's pid is {[p for j, p in enumerate(pids) if j != named[0][0]][0]})" if named
                      else f"no waiting line that names the other run's pid; the lines: {[[w['line'] for w in ws] for ws in a['waits']]}, pids {pids}"),
                     "the two runs' output (child.out.txt)"))
    # 2. complete
    complete, problems = context["complete"], context["complete"]["problems"]
    out.append(check("both runs' results are complete: every check of the chosen cases, each case's results, pictures and logs, the same files in both roots", not problems and not any(a["timed_out"])
                     and complete["cases_checked"] >= 1,
                     "every case folder of both roots has results.json (with at least the checks the case has: " + ", ".join(f"{k} {v}" for k, v in MIN_CHECKS.items())
                     + "), figure.png and its other pictures, run.log and runs.json; the two roots hold the same files; neither run had to be ended on its timeout; judged by completeness, "
                     "not by the exit status (the ruling)", ("complete: " + f"{complete['cases_checked']} case folders, checks per case {complete['checks_per_case']}, "
                                                              f"exit codes {a['returncodes']} (reported below)" if not problems and not any(a["timed_out"])
                                                              else f"INCOMPLETE: {problems[:6]}; timed out: {a['timed_out']}"), "the case folders of out\\bench\\B9\\run1 and run2"))
    # 3. agree
    numbers = context["numbers"]
    agree = numbers["compared"] > 0 and not numbers["different"] and not numbers["only_in_first"] and not numbers["only_in_second"]
    out.append(check("the two runs' numbers agree: the same deterministic tools give the same numbers", agree,
                     f"every number of every results.json (except each run's own seconds, process ids and looks at the process list) is the same in both runs: exactly for B3 (FEMM, Elmer "
                     f"and Gmsh gave numbers identical to the last bit in every real pair so far), to {TOLERANCE_BY_CASE['B6']:g} of itself for B6 (ROSS) ({OWN_LIMIT}). ROSS is the "
                     f"exception because its modes near 450,000 rpm (critical speeds of 446,165, 447,547 and 451,958 rpm, which mean nothing for a motor) were seen to differ: once, in one of "
                     f"the first concurrent runs, a B6 differed from a solo B6 by up to 4.9e-8; two other concurrent pairs were identical, and in a later one B6 differed again by up to 4.8e-8, in "
                     f"the same modes and with the same two values the eight solo runs at threads=2 gave (451,957.695 and 451,957.717 rpm for the highest; with threads=1 the solo runs "
                     f"gave 451,957.706 and 451,957.686); and 16 solo runs of B6's ROSS step on 2026-11-06 (8 with "
                     f"threads=2, 8 with threads=1, one after another through the slot) took two values in each group, up to 4.9e-8 and 4.3e-8 apart. Setting OMP_NUM_THREADS=1 (run_solver's "
                     f"threads=1) did not remove the variation (the count was set through the environment, not read back from ROSS's BLAS), and ROSS fixes its ARPACK start vector (v0 = "
                     f"ones: ross/rotor_assembly.py:2422 of ROSS 3.1.0), so a random start is not the cause either; the cause is not isolated. The limit of {TOLERANCE_BY_CASE['B6']:g} was "
                     f"set after that first run (it was 1e-9 before): 20 times above the largest difference seen and 100 times below the bench's tightest physical tolerance (0.01 %). Only "
                     f"numbers are compared: labels such as whirl directions are not, because ROSS's label for its highest mode (about 452,000 rpm, 45 times the 10,000 rpm B6 is about) "
                     f"changes between solo runs (Mixed in 12 of the 16, Forward in 4), and a comparison of labels would gate on that noise",
                     f"{numbers['compared']} numbers compared, {numbers['exact']} identical, the worst relative difference {numbers['worst_relative']:.2e}"
                     + ("; by case: " + ", ".join(f"{k} {v['exact']} of {v['compared']} identical (worst {v['worst_relative']:.1e}, limit {v['tolerance']:g})"
                                                  for k, v in numbers.get("per_case", {}).items()) if numbers.get("per_case") else "")
                     + (f"; DIFFERENT: {numbers['different'][:4]}" if numbers["different"] else "") + (f"; only in the first: {numbers['only_in_first'][:4]}" if numbers["only_in_first"] else "")
                     + (f"; only in the second: {numbers['only_in_second'][:4]}" if numbers["only_in_second"] else ""), "compare_runs on the two roots' results.json"))
    # 4. no solver process left
    left = context["left"]
    out.append(check("no solver process is left afterwards, and the slot is free", not left["alive_after"] and left["slot_free"] and left["failed_looks"] == 0 and left["looks"] > 0
                     and left["members"] > 0,
                     "no process that descended from B9 (the runs of both pairs, their pool workers and FEMMs, ElmerSolver, ElmerGrid, Gmsh, ROSS's python) is alive after the last run ended, "
                     "watched from before the first child to after the last; every look at the process list worked (a list that could not be read is not an empty one); "
                     "femm_queue.holder() says the slot is free",
                     f"{left['members']} processes were below B9 over the whole run, both pairs ({', '.join(left['names'])}); alive after: {left['alive_after'] or 'none'}; failed looks "
                     f"{left['failed_looks']} of {left['looks']}; slot free: {left['slot_free']}; the slot was looked at {left['slot_looks']} times ({left['slot_unreadable']} unreadable) "
                     f"and the root's lock {left['lock_looks']} times ({left['lock_unreadable']} unreadable)", "checks.matlab_check.JobTracker, femm_queue.holder()"))
    # 5. below normal, at most 8
    rules = context["rules"]
    classes = rules["classes_seen"]
    seen = {c for values in classes.values() for c in values}
    femm_peak = rules["peaks"].get("femm_headlessw.exe", 0)
    out.append(check(f"every case of B9's runs shows below-normal priority and at most {limit} workers or threads", not rules["problems"] and rules["max_workers"] <= limit and rules["max_threads"] <= limit
                     and femm_peak <= limit and bool(classes) and seen == {"below_normal"} and rules["looks"] > 0 and rules["records"] > 0,
                     f"every record of every runs.json of the cases B9's runs made (B3 and B6 in each of the two roots, and B6 in the root of the pair on one root, same): pools of at most {limit} "
                     f"workers whose workers and FEMMs were seen below normal, programs with 1 to {limit} threads seen below normal, the in-process record below normal; and from outside, "
                     f"while both pairs went: every process below B9 seen below normal, and never more than {limit} femm_headlessw.exe alive at one look (the pools took turns). B1 to B8's "
                     f"logs in the default root are not B9's runs: each of those cases has its own priority checks",
                     f"{rules['records']} records of both runs ({rules['observed']} observed), at most {rules['max_workers']} workers and {rules['max_threads']} threads; from outside, "
                     f"{rules['looks']} looks: " + "; ".join(f"{n}: {', '.join(v)}" for n, v in sorted(classes.items())) + f"; the most alive at one look: femm_headlessw.exe {femm_peak}, "
                     f"python.exe {rules['peaks'].get('python.exe', 0)} (the runs and their pool workers)" + (f"; PROBLEMS: {rules['problems'][:4]}" if rules["problems"] else ""),
                     "runs.json of every case, and bench.ChildWatcher over the whole run"))
    # 6. the same root
    waits_root = [w for w in same["waits"] if w["kind"] == "root" and same["pids"][0] is not None and w["pid"] == same["pids"][0]]
    ordered = serialized(same["c_ended_at"], same["d_first_work_at"]) if same["c_ended_at"] is not None else False
    ok6 = bool(waits_root) and ordered and not same["d_problems"] and not any(same["timed_out"]) and same["lock_free"]
    out.append(check("two runs on the same root are serialized by the root lock", ok6,
                     "the second run prints 'waiting for the bench root <ROOT>: held by '<label>' (pid N) ...' with the first run's pid; its first file in the root is newer than the first run's "
                     "summary.json (it emptied nothing, and wrote nothing, until the first had ended); it finishes complete; the lock is free afterwards (the ruling: 'a run holds a lock on its "
                     "root for the whole run')",
                     (f"the second said: {waits_root[0]['line']}; its first file in the root came {same['d_first_work_at'] - same['c_ended_at']:+.2f} s after the first run's summary.json "
                      f"(the first run's last file); the second's folder: " + ("complete" if not same["d_problems"] else str(same["d_problems"]))
                      + f"; lock free afterwards: {same['lock_free']}") if waits_root and ordered else
                     f"NOT SERIALIZED or incomplete: waiting lines {[w['line'] for w in same['waits']]}, first run ended {same['c_ended_at']}, second's first file {same['d_first_work_at']}, "
                     f"problems {same['d_problems']}, timed out {same['timed_out']}, lock free {same['lock_free']}", "the two runs on out\\bench\\B9\\same"))
    # reported beside them
    out.append(check("exit status and final line of each run (reported: B2 exits 1 by a standing ruling, any case can; the runs are judged by completeness)", True,
                     "each run ends with 'bench: PASS' or 'bench: FAIL' and an exit code",
                     f"two roots: exit codes {a['returncodes']}, {a['final_lines']}; one root: exit codes {same['returncodes']}, first '{same['c_final']}', second '{same['d_final']}'",
                     "the runs' output", gating=False, reason="the ruling: judge each run by completeness, not by its exit status"))
    return out
