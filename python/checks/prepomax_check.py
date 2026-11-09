"""prepomax_check: PrePoMax is the version we installed, and the solver it would run is CalculiX 2.23.

What it does, and why
    PrePoMax is the pre- and post-processor of the structural work: the window in which a CalculiX model is drawn, meshed and
    looked at. It is a window, so it cannot be driven from a script the way the solvers are; what a script can check is that the
    program is the version we installed and that the solver it is set to call is the one we expect. Task 1 pointed PrePoMax at
    CalculiX 2.23 (its own folder ships 2.22): this row is the check that this is still so.
    Two rows:
      * the version resource of PrePoMax.exe (the block of text inside a Windows program that Explorer shows under Properties >
        Details) says 3.0.0. Reading it starts nothing: femm_session.file_version does it with three calls into version.dll;
      * the solver. PrePoMax keeps its settings in settings.bin, next to PrePoMax.exe, as a .NET BinaryFormatter stream; the key
        `_calculixExe` of the object PrePoMax.SolverSettings is the solver it runs. nrbf.py reads that file without starting
        PrePoMax. A path inside PrePoMax's own folder is stored with a leading "#" (the shipped value is
        #\\Solver\\ccx_dynamic.exe, meaning <PrePoMax folder>\\Solver\\ccx_dynamic.exe); any other path is stored as it is.
        The program found there is then run with -v through femm_queue.run_solver and must say version 2.23.

If the setting cannot be read reliably (the file is not a stream this reader understands completely, or has no such key), the row
FAILS and says why. It is not a pass on the strength of a guess.

Numbers to expect (2026-11-06): PrePoMax.exe 3.0.0; settings.bin names ccx_static.exe in D:\\Apps\\CalculiX\\2.23, which prints
"This is Version 2.23" (and exits with 201: the exit code is not looked at, the text is).

Self-check (the answers follow).
  1. Why does the row not start PrePoMax and ask it?
  2. PrePoMax is updated and ships its own CalculiX 2.22 again, and the setting goes back to #\\Solver\\ccx_dynamic.exe. What does
     this row show?
Answers.
  1. PrePoMax opens a window and waits for a person; there is no "print the solver" option to run it with. The settings file is
     where it keeps the answer, and reading a file starts nothing.
  2. The second row fails: the path resolves to PrePoMax's own Solver folder, whose ccx_dynamic.exe prints "This is Version 2.22",
     and the row says "expected 2.23" next to what was measured.
"""
from __future__ import annotations

import re
import time
from pathlib import Path

import femm_queue

from . import FAIL, PASS, Missing, Row, fresh, need_file, need_tool, nrbf, save_log

EXPECT_VERSION = "3.0.0"
EXPECT_SOLVER = "2.23"
CLASS = "PrePoMax.SolverSettings"
KEY = "_calculixExe"


class SettingsError(Exception):
    """settings.bin could not be read reliably, or does not name a solver."""


def solver_path(settings: Path, program: Path) -> Path:
    """The solver PrePoMax would run, from its settings.bin. A path inside PrePoMax's folder is stored as "#\\relative\\path"."""
    try:
        members = nrbf.read_class_members(Path(settings).read_bytes(), CLASS)
    except (OSError, nrbf.NrbfError) as error:
        raise SettingsError(f"cannot read {Path(settings).name}: {error}") from error
    value = members.get(KEY)
    if not isinstance(value, str) or not value.strip():
        raise SettingsError(f"{Path(settings).name} has no {KEY} (it holds {value!r}; the members are {', '.join(members)})")
    if value.startswith("#"):
        return Path(program).parent / value[1:].lstrip("\\/")
    return Path(value)


def run(folder: Path):
    import femm_session

    program = need_tool("prepomax")
    started = time.perf_counter()
    version = femm_session.file_version(program) or "(none)"
    yield Row("PrePoMax", "PrePoMax.exe: the program's version resource", f"version {EXPECT_VERSION}", f"version {version}",
              PASS if version == EXPECT_VERSION else FAIL, round(time.perf_counter() - started, 2))

    case = "the solver that PrePoMax's settings.bin names, run with -v"
    expected = f"CalculiX {EXPECT_SOLVER}"
    started = time.perf_counter()
    settings = need_file("PrePoMax settings.bin", program.parent / "settings.bin", case)
    try:
        solver = solver_path(settings, program)
    except SettingsError as error:
        yield Row("PrePoMax", case, expected, f"cannot tell: {error}", FAIL, round(time.perf_counter() - started, 2))
        return
    if not solver.is_file():
        raise Missing(f"the solver of PrePoMax's {KEY}", solver, case)
    fresh(folder)
    result = femm_queue.run_solver([solver, "-v"], label="check_tools PrePoMax solver -v", cwd=folder, timeout_s=60, threads=1)
    save_log(folder, "solver_v", result)
    seconds = round(time.perf_counter() - started, 2)
    text = (result.stdout + result.stderr).strip()
    found = re.search(r"Version\s+(\d+\.\d+(?:\.\d+)?)", text)
    got = found.group(1) if found else ("timed out" if result.timed_out else f"no version in: {text[:80]!r}")
    yield Row("PrePoMax", case, expected, f"CalculiX {got} ({solver})", PASS if got == EXPECT_SOLVER else FAIL, seconds)
