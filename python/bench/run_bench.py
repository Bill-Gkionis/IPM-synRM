"""run_bench.py: run the bench's cases, print every number against its expectation, and exit with the verdict.

Why it exists
    The cases (bench/b1_femm.py, b2_gmsh.py, b3_elmer.py, b4_ccx.py, b5_prepomax.py, b6_ross.py, b7_salib.py, b8_matlab.py, b9_rules.py) each try one chain of tools on the motor's own geometry. This script is the one
    command that runs them: it empties a case's folder, runs the case, writes what the case returned as results.json and what it said
    as run.log and runs.json, and prints each check as one entry (what was expected and where the expectation comes from, what was
    measured, the verdict). Read bench/__init__.py first: it says what the bench is and what every case leaves.

Run (from the Python to FEMM folder; the cases that run a solver take the machine's one solver slot and wait, printing who holds it, if
another run has it)
    uv run python bench/run_bench.py                       every case
    uv run python bench/run_bench.py --only B1,B3          these cases (commas or repeated --only, any case)
    uv run python bench/run_bench.py --out out\\bench\\mine --only B6   the same under another root (a folder inside out\\; relative: from the project folder). B9 starts
                                                           two of these at once, each in a root of its own
    uv run python bench/run_bench.py --only B1 --mesh-study   B1 also solves a finer mesh (every block's size halved): the evidence that a
                                                           difference between two geometries is the mesh (it runs by itself when a number misses)
    uv run python -m bench.run_bench --only B2             the same, started as a module
    uv run pytest tests/test_bench*.py -q                  the tests: no solver is run
Exit status: 0 if every GATING check of every case passed, 1 if one failed or a case did not finish, 2 for a wrong --only or --out, 130 after Ctrl-C
(what was done is saved). A failed expectation is a FAIL line with the difference, not a hidden one; the tolerances and expected values are the
brief's and are never adjusted by a run. A check that is not gating (B3's plain-curve run: the verdict is the laminated run's) is printed as
REPORTED (PASS) or REPORTED (FAIL), with its own verdict and the reason, and counted apart: it cannot change the exit status, and it cannot
hide. A check that cannot be made on this machine (a file it reads is not there) is printed as SKIPPED, with the reason.
Writes ROOT\\summary.json (ROOT is out\\bench, or the --out; the cases of this run, and their checks counted: `gating` {checks, passed, failed}, `reported` {checks,
passed, failed, skipped}, `failed` and `reported_failed` naming what did not pass, by case) and, in ROOT\\<case>\\, what bench/__init__.py lists.

The root and its lock (Task 13): a run prints "bench: root <ROOT> (pid N)" first, then takes the lock of its root and keeps it until the last case is done. A second run on the same
root waits, printing "waiting for the bench root <ROOT>: held by 'bench run B1,B2' (pid N) since HH:MM", and so cannot empty the first run's case folders. Two runs on two roots go
on side by side and meet only at the machine's solver slot. `--only` rewrites the root's summary.json with only the cases it ran, which is why B9's two runs use two roots.

A case that stops with an error is a failed check "(did not finish)" that says what it was, with the traceback in ROOT\\<case>\\traceback.txt;
the cases after it still run. A tool that is not installed is that error with the tool's name in it.

Numbers to expect (the whole bench, `run_bench.py --mesh-study` (B1 also solves its finer mesh), on this PC, 2026-11-06 14:09: 784 s)
    141 checks: 119 gating (116 pass, 3 FAIL) and 22 reported, not gating (18 pass, 3 fail, 1 skipped). The last line it printed:
        bench: FAIL, 3 of 119 gating checks fail: B2 (3); 22 reported, not gating (18 pass, 3 fail, 1 skipped) (784 s)
    Exit status 1 whenever B2 runs: its 3 FAILs are the ruled chord error, not a fault. A straight-edged mesh of a curved boundary always loses the circular segments, and path 1's
    mesh of the band-free sector misses motorgeo's areas by that: N45 -0.172 %, Air +0.319 % and Coil -0.237 % against the brief's 0.1 % (b2_gmsh.py explains each to the printed
    digit; the tolerance is not moved after a run, because that would hide that to_gmsh's mesh sizes are the cause, and a to_gmsh change cures it, which is for the Electromagnetics
    workstream). `--only B1,B3` exits 0; any run that includes B2 exits 1 until then. The 3 reported FAILs: B4's "the note's 15 MPa against this rotor's nominal stress" and
    "the pressure model's bore against Lame's disc", and B5's "the importer reads the whole deck" (the *DLOAD CENTRIF line and the section's thickness); the skipped one is B5's
    PrePoMax automation (MCP) session.
    Each case alone (checks, seconds in the whole run): B1 16, 119 s with --mesh-study; B2 27, 11 s; B3 15, 26 s; B4 22, 135 s; B5 13, 242 s; B6 8, 23 s; B7 11, 44 s; B8 22, 35 s; B9 7, 148 s.
    Each case's own docstring has its numbers in full.

Self-check (the answers follow).
  1. B2 runs for a minute and B1 is started in a second terminal. What happens?
  2. Why does a case that raises not stop the whole run?
Answers.
  1. B2 uses Gmsh as a library, inside the process: it takes no solver slot (a mesh is seconds of CPU), so B1 starts at once. B2's one
     program, the Gmsh app that draws its picture, does take it, for a few seconds. If B1 asks while B2 has it, B1 prints
     "waiting for the solver slot: held by '...' (pid N) since HH:MM" every 30 s and goes on when it is free.
  2. The cases do not depend on each other except B3's reference solve, which B3 makes itself. A bug in one chain should not hide the
     verdicts of the others; the exit status still says that something failed.
"""
from __future__ import annotations

import argparse
import importlib
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Callable

HERE = Path(__file__).resolve().parent
P = HERE.parent
if str(P) not in sys.path:                    # `python bench/run_bench.py` puts bench\ first on the path, not the project folder
    sys.path.insert(0, str(P))

import bench  # noqa: E402  (after the path is right)

CASES = {"B1": "b1_femm", "B2": "b2_gmsh", "B3": "b3_elmer", "B4": "b4_ccx", "B5": "b5_prepomax", "B6": "b6_ross",
         "B7": "b7_salib", "B8": "b8_matlab", "B9": "b9_rules"}       # the case -> its module in this folder, in the order they run


def registry() -> dict[str, Callable]:
    """{case: run(folder, log) -> results dict}. The module is imported when its case is run, so `--only B2` does not import FEMM's."""
    def lazy(module: str) -> Callable:
        def run(folder, log):
            return importlib.import_module(f"bench.{module}").run(folder, log)
        return run
    return {case: lazy(module) for case, module in CASES.items()}


_default_registry = registry                   # main() has a parameter of the same name


def chosen_cases(known, only: list[str] | None) -> list[str] | None:
    """The cases to run, in the registry's order; None (after saying why on stderr) if --only names one that is not a case."""
    if not only:
        return list(known)
    wanted = {part.strip().upper() for text in only for part in text.split(",") if part.strip()}
    unknown = sorted(wanted - set(known))
    if unknown or not wanted:
        print(f"bench: --only {', '.join(unknown) or '(nothing)'} is not a case. The cases are: {', '.join(known)}", file=sys.stderr)
        return None
    return [case for case in known if case in wanted]


def verdict_word(item: dict) -> str:
    """What the table shows before a check's name: PASS or FAIL for a gating check; REPORTED (PASS), REPORTED (FAIL) or SKIPPED for one that does not gate."""
    if bench.is_gating(item):
        return item["verdict"]
    return "SKIPPED" if item["verdict"] == "SKIPPED" else f"REPORTED ({item['verdict']})"


def print_checks(case: str, checks: list[dict]) -> None:
    for item in checks:
        expected, measured = bench.describe(item)
        print(f"  {verdict_word(item)}  {case}: {item['name']}", flush=True)
        print(f"        expected: {expected}\n        measured: {measured}", flush=True)


def tally(done: dict[str, list[dict]]) -> dict:
    """What a run counted, by the rule of the exit status: the gating checks (passed, failed: any verdict but PASS fails), the reported ones, which do not
    gate (passed, failed, skipped), and the names that did not pass, by case. {"checks", "gating", "reported", "failed", "reported_failed"}."""
    everything = [c for checks in done.values() for c in checks]
    gating = [c for c in everything if bench.is_gating(c)]
    reported = [c for c in everything if not bench.is_gating(c)]
    passed = lambda items: sum(1 for c in items if c["verdict"] == "PASS")
    failed = {case: [c["name"] for c in checks if bench.is_gating(c) and c["verdict"] != "PASS"] for case, checks in done.items()}
    reported_failed = {case: [c["name"] for c in checks if not bench.is_gating(c) and c["verdict"] == "FAIL"] for case, checks in done.items()}
    return {"checks": len(everything),
            "gating": {"checks": len(gating), "passed": passed(gating), "failed": len(gating) - passed(gating)},
            "reported": {"checks": len(reported), "passed": passed(reported), "failed": sum(1 for c in reported if c["verdict"] == "FAIL"),
                         "skipped": sum(1 for c in reported if c["verdict"] == "SKIPPED")},
            "failed": {case: items for case, items in failed.items() if items},
            "reported_failed": {case: items for case, items in reported_failed.items() if items}}


def main(argv: list[str] | None = None, *, registry: dict | None = None) -> int:
    """`run_bench.py [--only CASES] [--out ROOT]`: run the cases and print the verdicts. 0 if every check passed, 1 if not, 2 for a wrong --only or --out, 130 after Ctrl-C."""
    parser = argparse.ArgumentParser(prog="run_bench.py", description="Run the bench's cases on the motor's geometry and print every check.")
    parser.add_argument("--only", action="append", metavar="CASES", help="run only these cases: comma-separated names, or the option repeated (B1, B2, B3)")
    parser.add_argument("--mesh-study", action="store_true", help="B1: also solve a finer mesh (every block's mesh size x 0.5), whether or not a number missed")
    parser.add_argument("--out", metavar="ROOT", help="write under ROOT instead of out\\bench: a folder inside out\\ (a relative one starts at the project folder); the run keeps a lock on it")
    args = parser.parse_args(argv)
    import femm_session
    femm_session.lower_priority()              # R10: the bench's own Python (Gmsh's library, numpy, matplotlib) runs below normal too; never raised, never idle -> normal
    bench.OPTIONS.clear()
    bench.OPTIONS["mesh_study"] = bool(args.mesh_study)
    for stream in (sys.stdout, sys.stderr):    # a Windows console cannot print every character a tool may say: never fail on one
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(errors="replace")
    known = _default_registry() if registry is None else registry
    names = chosen_cases(known, args.only)
    if names is None:
        return 2
    if args.out:
        try:
            chosen = bench.root_allowed(args.out)
        except ValueError as error:
            print(f"bench: --out {args.out}: {error}", file=sys.stderr)
            return 2
        with bench.use_root(chosen):
            return run_in_root(names, known)
    return run_in_root(names, known)


def run_in_root(names: list[str], known: dict) -> int:
    """Say where this run writes and who it is, take the lock of its root and run the cases while holding it. Returns the exit status.

    The first line ("bench: root <ROOT> (pid N)") is what B9 reads the pid of its children from: the interpreter that runs this is not always the process that was started (a
    virtual environment's launcher stands in front of it), so only the process itself can say."""
    print(f"bench: root {bench.root().resolve()} (pid {os.getpid()})", flush=True)
    try:
        with bench.hold_root(f"bench run {','.join(names)}"):
            return run_locked(names, known)
    except KeyboardInterrupt:                  # only the wait for the root can end here: run_locked answers its own Ctrl-C
        print("bench: interrupted while waiting for the root", file=sys.stderr)
        return 130


def run_locked(names: list[str], known: dict) -> int:
    """The cases, one at a time, while this run holds the lock of its root. Returns the exit status."""
    started = time.perf_counter()
    done: dict[str, list[dict]] = {}
    interrupted = False
    print(f"bench: {len(names)} case(s), one at a time: {', '.join(names)}", flush=True)
    try:
        for name in names:
            print(f"=== {name}", flush=True)
            folder = bench.fresh(bench.case_folder(name))
            log = bench.Log(folder)
            case_started = time.perf_counter()
            try:
                results = known[name](folder, log)
                checks = list(results.get("checks", []))
                results["seconds_total"] = round(time.perf_counter() - case_started, 1)
            except Exception as error:         # a bug, or a tool that did something unexpected: one line here, the rest in a file
                (folder / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
                checks = [bench.check("(did not finish)", False, "the case runs to its end", f"error: {type(error).__name__}: {error}", "run_bench.py")]
                results = {"case": name, "error": f"{type(error).__name__}: {error}", "checks": checks,
                           "seconds_total": round(time.perf_counter() - case_started, 1)}
                print(f"      the traceback is in {folder / 'traceback.txt'}", flush=True)
            bench.write_json(folder / "results.json", results)
            log.add_run(bench.process_record())
            log.save()
            print_checks(name, checks)
            done[name] = checks
    except KeyboardInterrupt:                  # run_solver and the pool have ended whatever they started; say what was done
        interrupted = True
        print("bench: interrupted", file=sys.stderr)
    seconds = time.perf_counter() - started
    counts = tally(done)
    failed, gating, reported = counts["failed"], counts["gating"], counts["reported"]
    bench.root().mkdir(parents=True, exist_ok=True)
    bench.write_json(bench.root() / "summary.json", {"cases": list(done), **counts, "interrupted": interrupted, "seconds": round(seconds, 1),
                                                       "time": time.strftime("%Y-%m-%d %H:%M:%S")})
    beside = (f"; {reported['checks']} reported, not gating ({reported['passed']} pass, {reported['failed']} fail, {reported['skipped']} skipped)"
              if reported["checks"] else "")
    print()
    if interrupted or failed or not done:
        who = ", ".join(f"{name} ({len(items)})" for name, items in failed.items()) or "none finished"
        print(f"bench: FAIL, {gating['failed']} of {gating['checks']} gating checks fail: {who}" + (" (interrupted)" if interrupted else "") + f"{beside} ({seconds:.0f} s)")
        return 130 if interrupted else 1
    print(f"bench: PASS, all {gating['checks']} gating checks pass{beside} ({seconds:.0f} s)")
    return 0


if __name__ == "__main__":                     # a pool worker is a new Python that imports this file: the guard keeps it from running the bench again
    sys.exit(main())
