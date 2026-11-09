"""check_queue: the live check of femm_queue, with real FEMMs. Three runs start together and have to take turns.

What it does, and why
    femm_queue's tests (tests/test_femm_queue.py) use stand-in processes. This is the real thing, the check of Task 3 of the tools
    plan: two FEMM pool jobs (two workers each, each worker with a FEMM of its own, solving small magnet models) and one
    run_solver call, started within three seconds of each other, as if from three terminals. The solver slot lets one of them in;
    the other two must print their waiting line, wait, and run after it. Nothing may overlap, every worker and FEMM must be at
    below-normal priority, and when it is all over no FEMM of ours may be left running.

Method
    This script starts three processes of itself and prints their output together, each line with the seconds since the start:
        A  (at 0 s)   --job A    a femm_pool of 2 workers: 4 small solves, each worker holding the slot a little longer (--hold)
        B  (at 1 s)   --job B    the same, while A holds the slot
        C  (at 2 s)   --solver   run_solver of a stand-in program (a Python one-liner, NOT a real solver), while A and B want the slot
    Every job writes when it took the slot and when it was about to give it back (out/check_queue/<name>.json). The check
    then tests: (1) no two of those intervals overlap; (2) B and C printed a waiting line that names a holder; (3) every worker,
    FEMM and the stand-in solver ran at below-normal priority (the operating system's own answer, read back); (4) the solves
    are right (a magnet drives flux, no magnet drives none: nobody read another job's file); (5) every FEMM that a worker reported is
    gone, and the number of femm_headlessw.exe on the PC is what it was before (a FEMM window that somebody has open by hand is left alone, and
    counted on both sides); (6) the slot is free afterwards. It never ends a process itself: it only looks.
    The jobs use out/check_queue/<job>/ with a file name of its own per task: two FEMMs never share a file name.

Run (from this folder; about 9 s on 2026-11-05; nothing else heavy should be running, it takes the real slot)
    uv run python check_queue.py [--hold SECONDS]        --hold: how long each task keeps its worker, default 1.5 s
    Exit status 0 if every check passes, 1 if one fails. The internal options --job and --solver start one of the three runs.

Numbers to expect (run three times on 2026-11-05 with FEMM 5.0.20261101: 8.9 s each time; the pids and clock times differ every run)
    [  0.0 s] check_queue: the slot is free; FEMM programs on the PC before: 0
    [  0.3 s] A | [femm_pool check_queue job A] starting 2 worker(s), pids [21760, 14472], at below_normal priority, each with its own FEMM
    [  0.6 s] A | [femm_pool check_queue job A] all workers are ready
    [  1.3 s] B | waiting for the solver slot: held by 'check_queue job A' (pid 5432) since 20:05
    [  2.3 s] C | waiting for the solver slot: held by 'check_queue job A' (pid 5432) since 20:05
    [  3.4 s] B | waiting for the solver slot: held by 'check_queue job A' (pid 5432) since 20:05
    [  3.7 s] A | job A done: 4 solves, By = 0.1582, 0.0000, 0.1582, 0.0000 T
    [  5.4 s] C | solver done: exit 0, ran 1.1 s, priority class 0x4000, OMP_NUM_THREADS=2
    [  5.5 s] B | [femm_pool check_queue job B] starting 2 worker(s), pids [14056, 11256], at below_normal priority, ...
    [  8.9 s] B | job B done: 4 solves, By = 0.1582, 0.0000, 0.1582, 0.0000 T
    [  8.9 s] the order in which the runs held the slot: A (20:05:44 to 20:05:47) -> C (20:05:48 to 20:05:49) -> B (20:05:49 to 20:05:52)
               B printed 2 waiting line(s)
               C printed 1 waiting line(s)
               4 workers, 4 FEMMs; priority classes seen: ['0x4000']
               FEMMs reported by the workers still running: 0; FEMM programs on the PC after: 0 (before: 0)
    [  8.9 s] PASS: the three runs took turns, none overlapped, nothing is left running
    A is first because it started first. C and B follow in the order in which they look at the lock after A lets go (it was C
    first in all three runs). B's and C's lines name A while A holds the slot. By = 0.1582 T is the magnet's field at its centre,
    the same number learn/03_pyfemm_magnet.py gets for the same model; 0.0000 T is the block with no magnetisation.
    Afterwards out/.femm/ holds only the lock file solver.lock (no note: holder.json is removed when a run gives the slot back).

Self-check (the answers follow).
  1. B asked for the slot one second after A. Why does B print "waiting for the solver slot: held by 'check_queue job A' (pid N)"
     and not simply wait?  Where does the pid come from?
  2. After A gives the slot back, B and C both want it. Which one gets it, and what does the other print next?
  3. Why does a job write the time it is about to give the slot back, and not the time it has given it back?
  4. The check counts femm_headlessw.exe before and after, and does not just expect 0. Why?
Answers.
  1. A line per poll tells whoever watches that the run is queued and why. The name and pid come from holder.json, the note the holder
     wrote after it took the lock (femm_queue.py explains why the note, and not the lock, carries them).
  2. Whichever of them looks at the lock first after the release; there is no queue order. The other one keeps waiting and its next
     line names the new holder (that is why C's lines can change from A to B).
  3. Once it has given the slot back another run can take it at once, so a time written afterwards could be later than that run's
     start, and the check would cry "overlap" for runs that did take turns. Written just before, it can only be early, so the check
     never gives a false alarm. (The price is a gap of a few milliseconds in which it could not see a real overlap. The lock, not
     this clock, is what stops runs from overlapping; the clock is a second opinion.)
  4. A FEMM window that somebody has open by hand, if there is one, is also a femm_headlessw.exe or femm_gui.exe on the PC and is none of this script's business.
     It must be the same before and after; only the FEMMs the workers reported are expected to be gone.
"""
from __future__ import annotations

import argparse
import ctypes
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from ctypes import wintypes
from pathlib import Path

import femm_queue
import paths

OUT = paths.path("out") / "check_queue"
BELOW_NORMAL = 0x4000                                           # BELOW_NORMAL_PRIORITY_CLASS
N45_COERCIVITY = 1_030_272                                      # A/m, the motor's magnet (FEMM: mu_r 1.05)


# ============================================================================================== what the runs do
def priority_class(pid: int) -> int:
    """The Windows priority class of a running process (0x4000 is BELOW_NORMAL), asked of the operating system."""
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.GetPriorityClass.argtypes = [wintypes.HANDLE]
    kernel32.GetPriorityClass.restype = wintypes.DWORD
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    handle = kernel32.OpenProcess(0x1000, False, pid)           # PROCESS_QUERY_LIMITED_INFORMATION
    if not handle:
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        return kernel32.GetPriorityClass(handle)
    finally:
        kernel32.CloseHandle(handle)


def solve_tiny(args):
    """One small model on this worker's own FEMM: a 10 x 2 mm magnet in a 60 mm box of air (learn/03's model), solved; By at its centre.

    args = (job, index, hold_s). Even indexes get the N45 magnet, odd ones a block with no magnetisation, so a job that read
    another job's answer file would show it. The file name has the job and the task in it, so no two FEMMs ever share one.
    The sleep keeps this worker (and so the slot) a little longer, long enough for the other runs to be seen waiting.
    """
    import femm
    import femm_session
    job, index, hold_s = args
    folder = OUT / f"job_{job}"
    folder.mkdir(parents=True, exist_ok=True)
    fem = folder / f"task{index}.fem"
    fem.with_suffix(".ans").unlink(missing_ok=True)             # a finished .ans next to the model can only be this run's answer
    femm.newdocument(0)                                         # 0 = magnetics
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, 20, 30)   # units, 2D planar, 20 mm deep
    femm.mi_getmaterial("Air")
    femm.mi_addmaterial("magnet", 1.05, 1.05, N45_COERCIVITY if index % 2 == 0 else 0)
    femm.mi_drawrectangle(-30, -30, 30, 30)                     # the box of air
    femm.mi_drawrectangle(-5, -1, 5, 1)                         # the magnet
    for (x, y), material, direction in (((0, 0), "magnet", 90), ((0, 20), "Air", 0)):
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop(material, 1, 0, "", direction, 0, 0)
        femm.mi_clearselected()
    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # no flux leaves the box
    for x, y in ((0, -30), (30, 0), (0, 30), (-30, 0)):
        femm.mi_selectsegment(x, y)
    femm.mi_setsegmentprop("A=0", 0, 1, 0, 0)
    femm.mi_clearselected()
    femm.mi_saveas(str(fem))
    femm.mi_analyze(1)                                          # 1 = no progress window
    femm.mi_loadsolution()
    by = femm.mo_getb(0, 0)[1]
    time.sleep(hold_s)
    femm_pid = femm_session.pid(femm._default)
    return {"index": index, "worker": os.getpid(), "femm": femm_pid, "by": by,
            "worker_class": priority_class(os.getpid()), "femm_class": priority_class(femm_pid)}


STAND_IN = ("import ctypes, os, time\n"
            "k = ctypes.WinDLL('kernel32'); k.GetCurrentProcess.restype = ctypes.c_void_p; k.GetPriorityClass.argtypes = [ctypes.c_void_p]\n"
            "t0 = time.time(); time.sleep(1)\n"
            "print(t0, time.time(), k.GetPriorityClass(k.GetCurrentProcess()), os.environ.get('OMP_NUM_THREADS'))\n")


def run_job(name: str, hold_s: float) -> None:
    """A femm_pool job: two workers, four tasks. Writes out/check_queue/<name>.json with the times it held the slot."""
    OUT.mkdir(parents=True, exist_ok=True)
    with femm_queue.femm_pool(f"check_queue job {name}", workers=2, task_timeout_s=60, poll_s=2) as pool:
        entered = time.time()
        results = pool.map(solve_tiny, [(name, index, hold_s) for index in range(4)])
        left = time.time()                                      # inside the slot: the next run cannot have started yet
    (OUT / f"{name}.json").write_text(json.dumps({"name": name, "entered": entered, "left": left, "results": results}))
    print(f"job {name} done: {len(results)} solves, By = " + ", ".join(f"{r['by']:.4f}" for r in results) + " T", flush=True)


def run_stand_in_solver() -> None:
    """A run_solver call (of a stand-in program, not a solver): below normal, 2 threads, under a timeout, under the same slot."""
    OUT.mkdir(parents=True, exist_ok=True)
    run = femm_queue.run_solver([sys.executable, "-c", STAND_IN], label="check_queue stand-in solver", cwd=OUT, timeout_s=60,
                                threads=2, poll_s=2)
    started, ended, priority, threads = run.stdout.split()
    (OUT / "solver.json").write_text(json.dumps({"name": "C", "entered": float(started), "left": float(ended),
                                                 "class": int(priority), "threads": threads, "returncode": run.returncode}))
    print(f"solver done: exit {run.returncode}, ran {run.seconds:.1f} s, priority class {int(priority):#x}, OMP_NUM_THREADS={threads}", flush=True)


# ============================================================================================== the check
def femm_programs() -> int:
    """How many FEMM programs run on the PC now (tasklist, by image name: this only COUNTS, it never ends anything)."""
    out = subprocess.run(["tasklist", "/NH"], capture_output=True, text=True).stdout.lower()
    return sum(1 for line in out.splitlines() if line.startswith(("femm_headlessw.exe", "femm_headless.exe", "femm_gui.exe", "femm.exe")))


def alive(pid: int) -> bool:
    out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"], capture_output=True, text=True).stdout
    return str(pid) in out


def check(hold_s: float) -> int:
    if OUT.exists():
        shutil.rmtree(OUT)                                      # out/check_queue of the last run: the one folder this script writes
    OUT.mkdir(parents=True)
    started = time.monotonic()
    lines: dict[str, list[str]] = {"A": [], "B": [], "C": []}
    lock = threading.Lock()

    def say(who: str, text: str) -> None:
        with lock:
            lines[who].append(text)
            print(f"[{time.monotonic() - started:5.1f} s] {who} | {text}", flush=True)

    def launch(who: str, *arguments: str) -> tuple[subprocess.Popen, threading.Thread]:
        proc = subprocess.Popen([sys.executable, "-u", str(Path(__file__).resolve()), *arguments], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", cwd=Path(__file__).resolve().parent)
        reader = threading.Thread(target=lambda: [say(who, line.rstrip()) for line in proc.stdout], daemon=True)
        reader.start()
        return proc, reader

    before = femm_programs()
    held = femm_queue.holder()
    print(f"[  0.0 s] check_queue: the slot is " + ("free" if held is None else f"held by {held['label']!r}, so A waits its turn too")
          + f"; FEMM programs on the PC before: {before}", flush=True)
    runs = [launch("A", "--job", "A", "--hold", str(hold_s))]
    time.sleep(1)
    runs.append(launch("B", "--job", "B", "--hold", str(hold_s)))
    time.sleep(1)
    runs.append(launch("C", "--solver"))
    problems: list[str] = []
    for who, (proc, reader) in zip("ABC", runs):
        try:
            proc.wait(timeout=300)
        except subprocess.TimeoutExpired:                       # nothing is ended from here: each run has its own timeouts
            problems.append(f"run {who} (pid {proc.pid}) was still going after 300 s and is left alone")
        reader.join(10)
    seconds = time.monotonic() - started

    notes = {name: json.loads((OUT / f"{name}.json").read_text()) for name in ("A", "B", "solver") if (OUT / f"{name}.json").exists()}
    if len(notes) != 3 or any(proc.returncode != 0 for proc, _ in runs):
        problems.append(f"a run did not finish cleanly (exit codes {[proc.returncode for proc, _ in runs]}, notes {sorted(notes)})")
    else:
        spans = sorted((note["entered"], note["left"], note["name"]) for note in notes.values())
        for (_, left, first), (entered, _, second) in zip(spans, spans[1:]):
            if entered < left:
                problems.append(f"{first} and {second} overlapped by {left - entered:.2f} s")
        clock = lambda moment: time.strftime("%H:%M:%S", time.localtime(moment))
        order = " -> ".join(f"{name} ({clock(entered)} to {clock(left)})" for entered, left, name in spans)
        print(f"[{seconds:5.1f} s] the order in which the runs held the slot: {order}", flush=True)
    for who in ("B", "C"):
        waiting = [text for text in lines[who] if text.startswith("waiting for the solver slot: held by '")]
        print(f"           {who} printed {len(waiting)} waiting line(s)", flush=True)
        if not waiting:
            problems.append(f"{who} never printed a waiting line")
    reported = [r for name in ("A", "B") if name in notes for r in notes[name]["results"]]
    femm_pids = {r["femm"] for r in reported}
    worker_pids = {r["worker"] for r in reported}
    classes = {r["worker_class"] for r in reported} | {r["femm_class"] for r in reported} | ({notes["solver"]["class"]} if "solver" in notes else set())
    print(f"           {len(worker_pids)} workers, {len(femm_pids)} FEMMs; priority classes seen: {sorted(hex(c) for c in classes)}", flush=True)
    if classes != {BELOW_NORMAL}:
        problems.append(f"priority classes {sorted(hex(c) for c in classes)}, expected only {BELOW_NORMAL:#x}")
    if len(femm_pids) != len(worker_pids) or len(worker_pids) != 4:
        problems.append(f"expected 4 workers with a FEMM each, saw {len(worker_pids)} workers and {len(femm_pids)} FEMMs")
    for r in reported:
        if (r["by"] < 0.1) if r["index"] % 2 == 0 else (abs(r["by"]) > 1e-6):
            problems.append(f"wrong answer for task {r['index']}: By = {r['by']}")
    still = [pid for pid in femm_pids | worker_pids if alive(pid)]
    after = femm_programs()
    print(f"           FEMMs reported by the workers still running: {len(still)}; FEMM programs on the PC after: {after} (before: {before})", flush=True)
    if still or after != before:
        problems.append(f"left running: {still}; FEMM programs {before} before, {after} after")
    if femm_queue.holder() is not None:
        problems.append(f"the slot is still held: {femm_queue.holder()}")
    print(f"[{seconds:5.1f} s] " + ("PASS: the three runs took turns, none overlapped, nothing is left running" if not problems
                                  else "FAIL:\n  " + "\n  ".join(problems)), flush=True)
    return 1 if problems else 0


if __name__ == "__main__":                                      # a pool worker is a new Python that imports this file: the guard keeps it from re-running the check
    parser = argparse.ArgumentParser(description="Live check of femm_queue: three runs start together and take turns.")
    parser.add_argument("--hold", type=float, default=1.5, help="seconds each task keeps its worker (default 1.5)")
    parser.add_argument("--job", choices=("A", "B"), help="(internal) run one femm_pool job")
    parser.add_argument("--solver", action="store_true", help="(internal) run the stand-in solver")
    options = parser.parse_args()
    if options.job:
        run_job(options.job, options.hold)
    elif options.solver:
        run_stand_in_solver()
    else:
        sys.exit(check(options.hold))
