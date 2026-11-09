"""check_tools.py: every tool's round trip, in one command.

Why it exists
    A tool that starts and prints its version is not yet a tool that works. Before the test bench uses FEMM, Gmsh, CalculiX,
    Elmer, ROSS, SALib, pymoo or MATLAB on the motor's geometry, each of them has to give ONE known answer, driven from Python the
    way the pipeline will drive it. This script does that for all nine (eleven rows, because Gmsh and PrePoMax have two cases
    each), prints a table, and exits 0 only when every row passes. Run it again after any tool is installed, updated or moved:
    a row that fails names the tool, and what it measured against what it expected.

What each row proves (the cases are in checks\\, one module each; the first line of each module's docstring says what it does)
    FEMM        check_femm.py's 10 A solve: the torque from the air-gap element and from the weighted stress tensor agree to 0.1 %
    Gmsh        motorgeo.to_gmsh meshes the default rotor (> 1,000 nodes), and the app prints 4.15.2
    CalculiX    a ring that spins at 10,000 rpm: the hoop stress at the bore is within 0.5 % of Lame's 14.098 MPa
    PrePoMax    the program's version resource is 3.0.0, and the solver its settings.bin names is CalculiX 2.23
    ROSS        a 300 mm shaft of 15 mm steel on two very stiff bearings: its first natural frequency is within 2 % of the formula
    Elmer       a round conductor of 10 mm radius carrying 100 A: B is within 2 % of 2.0 mT at r = R and of 1.0 mT at r = R/2
    SALib       the Sobol first-order indices of the Ishigami function are within 0.05 of 0.314, 0.442 and 0.000
    pymoo       NSGA-II on ZDT1 (40 individuals, 250 generations) ends within 0.05 (IGD) of the true front
    MATLAB      matlab -batch loads a copy of Improved_plant.slx without simulating it, then simulates a small model it builds

How it runs the tools (R10 of the tools plan, femm_queue.py)
    One at a time. Every program (ElmerSolver, ccx, matlab, the python of ROSS's environment, even `gmsh -version`) is started
    through femm_queue.run_solver: it takes the machine-wide solver slot, runs at below-normal priority under a timeout, in a
    Windows job that the system ends with this script, and nothing is ever ended by its process number (no taskkill): the job
    ends what it holds, at the timeout and when this script dies, however it dies. If another session holds the slot, this script
    prints who and waits. Python that runs in-process (Gmsh's API, SALib, pymoo) takes no slot, and runs below normal priority all the
    same: the script lowers its own process when it starts, as run_bench.py does (it is seconds of CPU, but R10 says below normal).
    Every path comes from config.toml; every case works in its own folder, out\\check_tools\\<name>\\, emptied at the start of
    the case (the scratch of other tasks in out\\check_tools is never touched).
    One run at a time: a case empties its folder BEFORE it asks for the solver slot, so a second run started during the first would
    delete files the first is using. A run therefore holds the lock of out\\check_tools (the folder .lock inside it; femm_queue's own
    operating-system lock, as run_bench.py's is for a bench root) from its first row to its last, and a second run prints "waiting
    for another check_tools run: held by '<label>' (pid N) since HH:MM" and starts when the first is done (Ctrl-C while it waits
    ends it, status 130). The system frees the lock when a run dies, however it dies.

Run (from the Python to FEMM folder)
    uv run python check_tools.py                       all rows
    uv run python check_tools.py --only gmsh,salib     only these rows (the names: femm, gmsh, calculix, prepomax, ross, elmer,
                                                       salib, pymoo, matlab; commas or repeated --only)
    uv run pytest tests/test_check_tools.py -q         the tests of this script: no tool is run
Writes out\\check_tools.csv with the columns  tool, case, expected, measured, verdict, seconds  (the rows of this run only: a run
with --only replaces the file). Exit status 0 if every row passed, 1 if one failed or a tool is missing, 2 for a wrong --only,
130 after Ctrl-C (the rows made so far are in the CSV).
`seconds` is the time of the row's own part, including any wait for the solver slot (a wait is printed as it happens). What each
program printed is kept in its case folder (<program>.out.txt and .err.txt) for the look afterwards.

A tool that is not there is a row, not a crash: `missing: ccx (D:\\...\\ccx_static.exe)` in the measured column, exit status 1.
A case that stops with an error is a failed row `error: <what>`; its traceback is in out\\check_tools\\<name>\\traceback.txt.

Numbers to expect (2026-11-06, this PC; a full run takes about 30 s on a quiet PC; the MATLAB row is the one that stretches when
other sessions are busy)
    FEMM        T_age 0.8087 and T_wst 0.8085 N.m, 0.025 % apart                          1.6 s
    Gmsh        9,176 nodes (MSH 2.2); the app prints 4.15.2                              0.3 s, 0.1 s
    CalculiX    bore hoop stress 14.1060 MPa against Lame's 14.0977 (+0.06 %)            3 s
    PrePoMax    version 3.0.0; its settings.bin names ccx_static.exe, which says 2.23     0.0 s, 0.0 s
    ROSS        f1 = 339.17 Hz against the formula's 340.18 Hz (-0.30 %)                  7 s
    Elmer       |B(R)| 1.9904 mT (-0.48 %), |B(R/2)| 1.0019 mT (+0.19 %)                  1.4 s
    SALib       S1 = 0.316, 0.438, 0.001 against 0.314, 0.442, 0.000                      0.1 s
    pymoo       IGD 0.0134 against the limit 0.05 (40 individuals x 250 generations, seed 1)       0.4 s
                (the brief's 40 x 100 never reached 0.05, on any of 100 seeds: a design decision set 250 generations on
                2026-11-06; the evidence is in checks\\pymoo_check.py)
    MATLAB      exit code 0; the job held matlab.exe, MATLAB.exe and MathWorksServiceHost.exe client processes, and was empty
                within 0.05 s of the last matlab.exe                                      11 s
    With every row passing the verdict line reads "check_tools: PASS, all 11 rows pass" and the exit status is 0.

Self-check (the answers follow).
  1. ccx_static.exe -v exits with 201, and ElmerSolver with no case exits with 1. Why do the version checks look at the printed
     text and not at the exit code?
  2. A row passes with `measured` 14.106 MPa against Lame's 14.098. Why is that a pass, and what would a 14.5 be?
  3. The CSV is overwritten by a run with --only. Why is that better than merging the new rows into the old file?
Answers.
  1. Many programs use the exit code for something else than "all is well": ccx returns 201 for "I only printed my version" and
     ElmerSolver returns 1 because no case was given. What proves that the right program of the right version ran is the text it
     printed, so that is what is compared.
  2. 0.06 % is far inside the 0.5 % allowed, so the mesh is fine enough and the solver agrees with the exact stress. 14.5 is 2.9 %
     high: that would say the solver or the model is wrong, which is exactly what the row is there to catch.
  3. A merged file would mix rows of different days and different tool versions, and could show a pass for a tool that has since
     been moved. The CSV is "what this run measured", nothing older.
"""
from __future__ import annotations

import argparse
import csv
import sys
import textwrap
import time
import traceback

import checks
import paths
from checks import COLUMNS, FAIL, PASS, Missing, Row

CSV_NAME = "check_tools.csv"
WIDTHS = {"tool": 10, "case": 40, "expected": 30, "measured": 46, "verdict": 7, "seconds": 7}      # the table's columns, in characters


def format_table(rows: list[Row]) -> str:
    """The rows as a plain table; a long text is wrapped over lines inside its column and never cut."""
    lines: list[str] = []

    def emit(cells: list[str]) -> None:
        wrapped = [textwrap.wrap(cell, WIDTHS[column], break_long_words=True, break_on_hyphens=False) or [""]
                   for cell, column in zip(cells, COLUMNS)]
        for i in range(max(len(part) for part in wrapped)):
            lines.append("  ".join((part[i] if i < len(part) else "").ljust(WIDTHS[column])
                                   for part, column in zip(wrapped, COLUMNS)).rstrip())

    emit(list(COLUMNS))
    lines.append("  ".join("-" * WIDTHS[column] for column in COLUMNS))
    for row in rows:
        emit([row.tool, row.case, row.expected, row.measured, row.verdict, f"{row.seconds:.1f}"])
    return "\n".join(lines)


def write_csv(rows: list[Row]):
    """out\\check_tools.csv: one line per row, the six columns, UTF-8."""
    out = paths.path("out")
    out.mkdir(parents=True, exist_ok=True)
    target = out / CSV_NAME
    with open(target, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow([row.tool, row.case, row.expected, row.measured, row.verdict, f"{row.seconds:.2f}"])
    return target


def chosen_names(registry: dict, only: list[str] | None) -> list[str] | None:
    """The row names to run, in the registry's order; None (after saying why on stderr) if --only names one that is not a row."""
    if not only:
        return list(registry)
    wanted = {part.strip().lower() for text in only for part in text.split(",") if part.strip()}
    unknown = sorted(wanted - set(registry))
    if unknown or not wanted:
        print(f"check_tools: --only {', '.join(unknown) or '(nothing)'} is not a row. The rows are: {', '.join(registry)}", file=sys.stderr)
        return None
    return [name for name in registry if name in wanted]


def run_one(name: str, check: checks.Check, rows: list[Row]) -> None:
    """Run one check and append its rows (printing each as it is made). A missing tool or an error becomes a failed row."""
    folder = checks.case_folder(name)
    (folder / "traceback.txt").unlink(missing_ok=True)      # one left by an earlier failed run must not pass for this run's (some folders are not emptied)
    started = time.perf_counter()
    made = 0

    def keep(row: Row) -> None:
        nonlocal made
        made += 1
        rows.append(row)
        print(f"  {row.verdict}  {row.tool}: {row.case}  ({row.seconds:.1f} s)", flush=True)
        print(f"        expected: {row.expected}\n        measured: {row.measured}", flush=True)       # whole, not wrapped like the table

    try:
        with checks.intact_path():                  # a library (Gmsh's) may cut this process's PATH short: the next check must not suffer
            for row in check.run(folder):
                keep(row)
    except Missing as missing:                      # one row for each tool that is not there (a case may need several)
        for item in (missing, *missing.also):
            keep(Row(check.label, item.case, "the tool is installed", f"missing: {item.name} ({item.where})", FAIL, 0.0))
    except Exception as error:                      # a bug, or a tool that did something unexpected: one line here, the rest in a file
        trace = folder / "traceback.txt"
        trace.write_text(traceback.format_exc(), encoding="utf-8")
        keep(Row(check.label, "(did not finish)", "the case runs to its end", f"error: {type(error).__name__}: {error}", FAIL,
                 round(time.perf_counter() - started, 2)))
        print(f"      the traceback is in {trace}", flush=True)
    if not made:
        keep(Row(check.label, "(no result)", "at least one row", "the check produced no row", FAIL, 0.0))


def main(argv: list[str] | None = None, *, registry: dict | None = None) -> int:
    """`check_tools.py [--only NAMES]`: run the rows, print the table, write the CSV. 0 if every row passed, 1 if not, 2 for a wrong --only.

    `registry` is {name: Check}; the tests hand in their own, so that no real tool is run."""
    parser = argparse.ArgumentParser(prog="check_tools.py", description="Run every tool's round trip and print the table.")
    parser.add_argument("--only", action="append", metavar="NAMES",
                        help="run only these rows: comma-separated names, or the option repeated (femm, gmsh, calculix, prepomax, "
                             "ross, elmer, salib, pymoo, matlab)")
    args = parser.parse_args(argv)
    import femm_session
    femm_session.lower_priority()                   # R10: this script's own Python (Gmsh's library, SALib, pymoo) runs below normal too; never raised
    for stream in (sys.stdout, sys.stderr):         # a Windows console cannot print every character a tool may say: never fail on one
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(errors="replace")
    registry = checks.registry() if registry is None else registry
    names = chosen_names(registry, args.only)
    if names is None:
        return 2
    try:
        with checks.hold_the_run(f"check_tools {','.join(names)}"):
            return run_locked(names, registry)
    except KeyboardInterrupt:                       # only the wait for the lock can end here: run_locked answers its own Ctrl-C
        print("check_tools: interrupted while waiting for another run", file=sys.stderr)
        return 130


def run_locked(names: list[str], registry: dict) -> int:
    """The rows, one at a time, while this run holds the lock of out\\check_tools. Returns the exit status."""
    rows: list[Row] = []
    started = time.perf_counter()
    print(f"check_tools: {len(names)} check(s), one at a time: {', '.join(names)}", flush=True)
    interrupted = False
    try:
        for name in names:
            print(f"- {name}", flush=True)
            run_one(name, registry[name], rows)
    except KeyboardInterrupt:                       # run_solver has ended whatever it started; say what was done and what was not
        interrupted = True
        print("check_tools: interrupted", file=sys.stderr)
    target = write_csv(rows)
    print()
    print(format_table(rows))
    print()
    print(f"csv: {target}")
    failed = [row for row in rows if row.verdict != PASS]
    seconds = time.perf_counter() - started
    if interrupted or failed or not rows:
        who = ", ".join(dict.fromkeys(row.tool for row in failed)) or "none finished"
        print(f"check_tools: FAIL, {len(failed)} of {len(rows)} rows fail: {who}" + (" (interrupted)" if interrupted else "") +
              f" ({seconds:.0f} s)")
        return 130 if interrupted else 1
    print(f"check_tools: PASS, all {len(rows)} rows pass ({seconds:.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
