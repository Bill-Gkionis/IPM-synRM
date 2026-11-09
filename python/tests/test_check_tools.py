"""Tests for check_tools.py and the cases under checks/: the framework and the pure logic of each case. No tool is run.

Run from P:   uv run pytest tests/test_check_tools.py -q      (a few seconds; no FEMM, no solver, no MATLAB is started)

The real round trips (a solve in every tool) are what `uv run python check_tools.py` does; a test that did the same would only
repeat it. These tests pin what can go wrong around it, in four groups:
  1. The three behaviours of the brief: a tool that is missing gives a row that names it and exit code 1 (and no traceback), the
     CSV has the six columns, and --only runs only the rows it names. Then the rest of the framework: a case that fails with an
     error, rows made before a missing tool, a case folder that is only ever emptied inside out/check_tools.
  2. The known answers each case compares with, worked out here from the formulas (Lame's spinning ring, the Euler-Bernoulli
     shaft, Ishigami's Sobol indices, Ampere's law for a round conductor), so a wrong "expected" value in a case shows up as a
     failing test and not as a tool that seems to be broken.
  3. The two readers that have to be exactly right: the .NET BinaryFormatter reader that finds PrePoMax's solver in settings.bin,
     and the tracker that says which processes a MATLAB run left in its job (a process is known by its pid AND its creation time:
     the three ways a reused pid fooled it are tests), then the rows that cannot start (a broken settings.bin).
  4. The rule behind every case (R10 of the tools plan): a program is started only through femm_queue.run_solver, with no exception,
     and nothing ends a process by its number. The scans are tried on small sources first, so that a scan that finds nothing is
     known to be able to find something.
"""
import ast
import csv
import os
import struct
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

import check_tools
import checks
import femm_queue
import paths
from checks import FAIL, PASS, Check, Missing, Row
import queue_test_tools as qt
from bench_test_tools import the_test_run_is_not_lowered      # noqa: F401  (autouse, as in the bench's tests: check_tools.main() lowers its own process, and a test must not lower the whole pytest run)
from conftest import cut_process_path, process_path

P = Path(__file__).resolve().parent.parent
COLUMNS = ["tool", "case", "expected", "measured", "verdict", "seconds"]
NAMES = ["femm", "gmsh", "calculix", "prepomax", "ross", "elmer", "salib", "pymoo", "matlab"]


def use_config(monkeypatch, tmp_path, tools: str = "", extra_paths: str = "") -> Path:
    """A config.toml of the test's own (so that nothing is written under the real out\\, and a tool can be 'missing'). Returns its out."""
    out = tmp_path / "out"
    config = tmp_path / "config.toml"
    config.write_text(f"[tools]\n{tools}\n[paths]\nout = '{out}'\n{extra_paths}\n[femm]\nworkers = 8\npriority = 'below_normal'\n", encoding="utf-8")
    monkeypatch.setattr(paths, "CONFIG", config)
    return out


@pytest.fixture
def nothing_may_start(monkeypatch):
    """Any attempt to start a program, to open a pool or to take the machine's solver slot fails the test: these tests run no tool. check_tools.main() does take one lock, the one
    of its own run (out\\check_tools\\.lock, a lock folder of its own, which in these tests is in the test's temporary out): that is allowed, and only that."""
    def refuse(*args, **kwargs):
        raise AssertionError("a test of check_tools started a program or took the solver slot")
    real_slot = femm_queue.solver_slot

    def only_a_lock_folder_of_its_own(*args, lock_dir=None, **kwargs):
        if lock_dir is None:
            refuse()
        return real_slot(*args, lock_dir=lock_dir, **kwargs)
    monkeypatch.setattr(femm_queue, "run_solver", refuse)
    monkeypatch.setattr(femm_queue, "femm_pool", refuse)
    monkeypatch.setattr(femm_queue, "solver_slot", only_a_lock_folder_of_its_own)


def fake(label: str, rows=(), *, then=None, calls=None, name=""):
    """A Check that yields `rows`, then raises `then` if given; it notes its name in `calls` when it is started."""
    def run(folder):
        if calls is not None:
            calls.append(name)
        yield from rows
        if then is not None:
            raise then
    return Check(label, run)


def row(tool="Fake", case="a case", verdict=PASS, seconds=1.25, measured="1.0", expected="1.0") -> Row:
    return Row(tool, case, expected, measured, verdict, seconds)


def read_csv(out: Path) -> list[list[str]]:
    with open(out / "check_tools.csv", newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle))


# ----------------------------------------------------------------------------------------------- 1. the framework
def test_a_missing_tool_gives_a_row_that_names_it_and_exit_code_1(monkeypatch, tmp_path, capsys, nothing_may_start):
    where = "Z:\\no such folder\\ccx_static.exe"
    out = use_config(monkeypatch, tmp_path, f"ccx = '{where}'")
    assert check_tools.main(["--only", "calculix"]) == 1
    table = read_csv(out)
    assert table[0] == COLUMNS and len(table) == 2
    tool, case, expected, measured, verdict, seconds = table[1]
    assert tool == "CalculiX" and verdict == FAIL
    assert measured == f"missing: ccx ({where})"
    shown = capsys.readouterr().out
    assert f"missing: ccx ({where})" in shown and "Traceback" not in shown


def test_a_tool_the_config_does_not_name_or_leaves_empty_is_missing_in_the_same_words(monkeypatch, tmp_path, nothing_may_start):
    out = use_config(monkeypatch, tmp_path, "gmsh = ''")             # no ccx line at all
    assert check_tools.main(["--only", "calculix"]) == 1
    assert read_csv(out)[1][3] == "missing: ccx (not named in config.toml)"
    out = use_config(monkeypatch, tmp_path, "ccx = ''")
    assert check_tools.main(["--only", "calculix"]) == 1
    assert read_csv(out)[1][3] == "missing: ccx (empty in config.toml)"


def test_both_tools_of_a_case_are_named_when_both_are_missing_and_the_python_environment_of_ross_too(monkeypatch, tmp_path, nothing_may_start):
    """Elmer needs ElmerSolver and ElmerGrid: each missing one is a row of its own. ROSS lives in its own environment, whose
    python is built from paths.P (the project folder), so a project folder without envs\\ross is a missing tool named 'ross'."""
    out = use_config(monkeypatch, tmp_path, "elmer = 'Z:\\a\\ElmerSolver.exe'\nelmergrid = 'Z:\\a\\ElmerGrid.exe'")
    assert check_tools.main(["--only", "elmer"]) == 1
    assert [r[3] for r in read_csv(out)[1:]] == ["missing: elmer (Z:\\a\\ElmerSolver.exe)", "missing: elmergrid (Z:\\a\\ElmerGrid.exe)"]
    monkeypatch.setattr(paths, "P", tmp_path / "no project here")
    assert check_tools.main(["--only", "ross"]) == 1
    wanted = tmp_path / "no project here" / "envs" / "ross" / ".venv" / "Scripts" / "python.exe"
    assert read_csv(out)[1][3] == f"missing: ross ({wanted})"


def test_the_csv_has_the_six_columns_and_one_line_per_row(monkeypatch, tmp_path, nothing_may_start):
    out = use_config(monkeypatch, tmp_path)
    registry = {"one": fake("One", [row("One", "first case", seconds=1.254), row("One", "second case", FAIL, 0.5, "0.9", "1.0")])}
    assert check_tools.main([], registry=registry) == 1
    table = read_csv(out)
    assert table[0] == COLUMNS
    assert table[1:] == [["One", "first case", "1.0", "1.0", "PASS", "1.25"], ["One", "second case", "1.0", "0.9", "FAIL", "0.50"]]
    assert all(len(line) == 6 for line in table)


def test_a_text_with_commas_quotes_and_non_ascii_survives_the_csv(monkeypatch, tmp_path, nothing_may_start):
    out = use_config(monkeypatch, tmp_path)
    awkward = 'B = 2.0 mT, "within" 2 %, \u00b1 0.05 \u03c0'
    check_tools.main([], registry={"one": fake("One", [row(measured=awkward)])})
    assert read_csv(out)[1][3] == awkward


def test_exit_code_0_only_when_every_row_passes(monkeypatch, tmp_path, capsys, nothing_may_start):
    use_config(monkeypatch, tmp_path)
    good = {"a": fake("A", [row("A"), row("A", "again")]), "b": fake("B", [row("B")])}
    assert check_tools.main([], registry=good) == 0
    assert "PASS" in capsys.readouterr().out.strip().splitlines()[-1]
    bad = {"a": fake("A", [row("A")]), "b": fake("B", [row("B", verdict=FAIL)])}
    assert check_tools.main([], registry=bad) == 1
    last = capsys.readouterr().out.strip().splitlines()[-1]
    assert "FAIL" in last and "1 of 2" in last and "B" in last          # the verdict line says how many, and which tool


def test_only_runs_only_the_named_rows(monkeypatch, tmp_path, capsys, nothing_may_start):
    out = use_config(monkeypatch, tmp_path)
    calls: list[str] = []
    registry = {name: fake(name.upper(), [row(name.upper())], calls=calls, name=name) for name in ("a", "b", "c")}
    assert check_tools.main(["--only", "b"], registry=registry) == 0
    assert calls == ["b"] and [r[0] for r in read_csv(out)[1:]] == ["B"]
    calls.clear()
    assert check_tools.main(["--only", "C,a", "--only", "b"], registry=registry) == 0         # commas and repeats; the registry's order wins
    assert calls == ["a", "b", "c"]
    calls.clear()
    assert check_tools.main(["--only", "a,c"], registry=registry) == 0
    assert calls == ["a", "c"] and [r[0] for r in read_csv(out)[1:]] == ["A", "C"]


def test_a_name_that_is_not_a_row_is_refused_and_the_valid_names_are_listed(monkeypatch, tmp_path, capsys, nothing_may_start):
    use_config(monkeypatch, tmp_path)
    calls: list[str] = []
    registry = {name: fake(name, [row()], calls=calls, name=name) for name in ("alpha", "beta")}
    assert check_tools.main(["--only", "alpha,gamma"], registry=registry) == 2
    assert calls == []                                                  # nothing ran, not even the valid one
    shown = capsys.readouterr()
    assert "gamma" in shown.err + shown.out and "alpha" in shown.err + shown.out and "beta" in shown.err + shown.out


def test_the_real_registry_has_the_nine_rows_in_the_order_of_the_brief():
    assert list(checks.registry()) == NAMES
    assert [check.label for check in checks.registry().values()] == [
        "FEMM", "Gmsh", "CalculiX", "PrePoMax", "ROSS", "Elmer", "SALib", "pymoo", "MATLAB"]


def test_a_case_that_fails_with_an_error_is_a_failed_row_with_its_traceback_in_a_file(monkeypatch, tmp_path, capsys, nothing_may_start):
    out = use_config(monkeypatch, tmp_path)
    registry = {"broken": fake("Broken", [row("Broken", "made before")], then=ValueError("boom"))}
    assert check_tools.main([], registry=registry) == 1
    table = read_csv(out)
    assert table[1][4] == PASS and table[2][4] == FAIL                  # the row made before the error is kept
    assert table[2][3] == "error: ValueError: boom"
    shown = capsys.readouterr().out
    assert "Traceback" not in shown and "traceback.txt" in shown
    assert "ValueError: boom" in (out / "check_tools" / "broken" / "traceback.txt").read_text(encoding="utf-8")


def test_a_traceback_of_an_earlier_failed_run_is_removed_when_the_check_runs_again(monkeypatch, tmp_path, nothing_may_start):
    """A folder of earlier work (out\\check_tools\\matlab_mcp, say) is never emptied, so a traceback.txt would stay for ever."""
    out = use_config(monkeypatch, tmp_path)
    folder = out / "check_tools" / "once"
    folder.mkdir(parents=True)
    (folder / "traceback.txt").write_text("an old failure")
    (folder / "task 9's file.txt").write_text("keep")
    assert check_tools.main([], registry={"once": fake("Once", [row()])}) == 0
    assert not (folder / "traceback.txt").exists() and (folder / "task 9's file.txt").read_text() == "keep"


def test_rows_made_before_a_missing_tool_are_kept_and_the_missing_one_follows(monkeypatch, tmp_path, nothing_may_start):
    out = use_config(monkeypatch, tmp_path)
    registry = {"two": fake("Two", [row("Two", "python part")], then=Missing("gmsh", "Z:\\x\\gmsh.exe"))}
    assert check_tools.main([], registry=registry) == 1
    table = read_csv(out)
    assert [(line[1], line[4]) for line in table[1:]] == [("python part", PASS), ("(not started)", FAIL)]
    assert table[2][3] == "missing: gmsh (Z:\\x\\gmsh.exe)"


def test_a_check_that_yields_no_row_is_a_failure_not_a_pass(monkeypatch, tmp_path, nothing_may_start):
    """A case that quietly did nothing must not make the run look good."""
    out = use_config(monkeypatch, tmp_path)
    assert check_tools.main([], registry={"empty": fake("Empty", [])}) == 1
    assert read_csv(out)[1][4] == FAIL


def test_the_table_shows_every_column_and_does_not_cut_a_long_text():
    words = [f"word{i}" for i in range(60)]
    shown = check_tools.format_table([row("ToolX", "caseY", measured=" ".join(words)), row("ToolZ", verdict=FAIL, seconds=12.5)])
    for word in ("tool", "case", "expected", "measured", "verdict", "seconds", "ToolX", "caseY", "ToolZ", "PASS", "FAIL", "12.5"):
        assert word in shown
    assert all(word in shown.split() for word in words)                      # wrapped over several lines inside its column, not cut
    assert len(shown.splitlines()) > 6


def test_a_case_folder_is_emptied_only_inside_out_check_tools(monkeypatch, tmp_path):
    out = use_config(monkeypatch, tmp_path)
    folder = checks.case_folder("demo")
    assert folder == out / "check_tools" / "demo" and folder.is_dir()
    (folder / "sub").mkdir()
    (folder / "sub" / "old.txt").write_text("old")
    (folder / "old.txt").write_text("old")
    assert checks.fresh(folder) == folder and folder.is_dir() and list(folder.iterdir()) == []
    for refused in (tmp_path, out, out / "check_tools", out / "check_tools" / "demo" / ".." / ".." / "..", Path(paths.P)):
        with pytest.raises(ValueError):
            checks.fresh(refused)
    assert tmp_path.is_dir() and out.is_dir()


def test_the_folders_of_earlier_tasks_are_refused_by_name_and_only_matlab_mcp_run_may_be_emptied(monkeypatch, tmp_path):
    """Review M5: asked, fresh() emptied out\\check_tools\\matlab_mcp and verify (Task 9's test and Task 10's scratch)."""
    out = use_config(monkeypatch, tmp_path)
    root = out / "check_tools"
    for name in ("verify", "matlab_mcp"):
        (root / name / "sub").mkdir(parents=True)
        (root / name / "theirs.txt").write_text("keep")
    (root / "matlab_mcp" / "work").mkdir()
    for refused in (root / "verify", root / "verify" / "sub", root / "VERIFY", root / "matlab_mcp", root / "matlab_mcp" / "work",
                    root / "matlab_mcp" / "sub", root / "Matlab_MCP"):
        with pytest.raises(ValueError, match="earlier task"):
            checks.fresh(refused)
    assert (root / "verify" / "theirs.txt").read_text() == "keep" and (root / "matlab_mcp" / "theirs.txt").read_text() == "keep"
    run = root / "matlab_mcp" / "run"
    (run / "results").mkdir(parents=True)
    (run / "results" / "old.txt").write_text("old")
    assert checks.fresh(run) == run and list(run.iterdir()) == []                  # the one sub-folder of matlab_mcp that may be emptied
    assert (root / "matlab_mcp" / "theirs.txt").read_text() == "keep" and (root / "matlab_mcp" / "work").is_dir()
    assert checks.fresh(root / "ross") == root / "ross"                              # and an ordinary case folder, as before


def test_the_framework_makes_only_the_folders_of_the_rows_it_runs(monkeypatch, tmp_path, nothing_may_start):
    """out\\check_tools also holds the scratch of earlier tasks (matlab_mcp\\, verify\\); a run must leave all of it alone."""
    out = use_config(monkeypatch, tmp_path)
    (out / "check_tools" / "verify").mkdir(parents=True)
    (out / "check_tools" / "verify" / "keep.txt").write_text("keep")
    check_tools.main([], registry={"one": fake("One", [row()])})
    assert (out / "check_tools" / "verify" / "keep.txt").read_text() == "keep"
    assert sorted(p.name for p in (out / "check_tools").iterdir()) == [".lock", "one", "verify"]          # (the run's own lock folder is the one thing it adds)


# ----------------------------------------------------------------------------------------------- 1b. the run's own priority and lock (final review: B5 and N2)
def test_check_tools_lowers_its_own_process_as_the_bench_does(monkeypatch, tmp_path, nothing_may_start):
    """B5 of the final review. The rows that run in this very process (Gmsh's library, SALib, pymoo) are seconds of CPU, and they ran at normal priority while
    run_bench.py lowers its own process at the start (R10: below normal). check_tools.main() does the same now, through femm_session.lower_priority. (The autouse fixture
    imported at the top of this file stands in for that call in the other tests, as it does for the bench's, so that no test lowers the whole pytest run.)"""
    import femm_session
    use_config(monkeypatch, tmp_path)
    lowered = []
    monkeypatch.setattr(femm_session, "lower_priority", lambda: lowered.append(True))
    assert check_tools.main([], registry={"one": fake("One", [row()])}) == 0
    assert lowered == [True]


def test_a_run_holds_the_lock_of_out_check_tools_for_its_whole_length_and_lets_go_of_it(monkeypatch, tmp_path, nothing_may_start):
    """N2 of the final review (and the ledger's concern 3): every row empties its case folder, so two runs at once would empty each other's folders. One lock for the whole
    run, taken with femm_queue's operating-system lock on out\\check_tools\\.lock (as bench.hold_root does for a bench root, and not the machine's solver slot): held while a row
    runs, naming this run, and free afterwards, however the run ends."""
    out = use_config(monkeypatch, tmp_path)
    lock_dir = out / "check_tools" / ".lock"
    seen = []

    def run(folder):
        seen.append(femm_queue.holder(lock_dir))
        yield row()

    assert check_tools.main([], registry={"one": Check("One", run), "two": Check("Two", run)}) == 0
    assert [(note["label"], note["pid"]) for note in seen] == [("check_tools one,two", os.getpid())] * 2
    assert femm_queue.holder(lock_dir) is None
    failing = {"boom": fake("Boom", [row()], then=ValueError("boom"))}
    assert check_tools.main([], registry=failing) == 1
    assert femm_queue.holder(lock_dir) is None


def test_a_second_run_waits_for_the_first_and_says_who_holds_the_lock_instead_of_emptying_its_folders(monkeypatch, tmp_path, capsys, nothing_may_start):
    """N2: another check_tools run holds the lock (a stand-in process that takes the same lock folder). The second run prints which run it waits for and starts no row; when
    the first lets go, it runs. Without the lock it would have started at once and emptied the folders the first is using."""
    out = use_config(monkeypatch, tmp_path)
    lock_dir = out / "check_tools" / ".lock"
    calls: list[str] = []
    result = {}
    with qt.Holder(lock_dir, "check_tools femm,gmsh") as first:
        waiter = threading.Thread(target=lambda: result.update(code=check_tools.main([], registry={"one": fake("One", [row()], calls=calls, name="one")})), daemon=True)
        waiter.start()
        time.sleep(1.5)
        assert waiter.is_alive() and calls == [], "the second run started while the first held the lock"
        text = capsys.readouterr().out
        assert "waiting for another check_tools run" in text and "'check_tools femm,gmsh'" in text and f"(pid {first.pid})" in text, text
        first.release()
        waiter.join(30)
        assert not waiter.is_alive(), "the second run never got the lock after the first let go"
    assert result["code"] == 0 and calls == ["one"]
    assert femm_queue.holder(lock_dir) is None


def test_ctrl_c_while_a_run_waits_for_the_lock_ends_it_with_130_and_takes_nothing(monkeypatch, tmp_path):
    """N2: the wait for the lock can be interrupted (a stand-in process gets the SIGINT a Ctrl-C raises while it waits behind a holder). It returns 130, runs no row, and leaves
    the holder's lock alone."""
    out = use_config(monkeypatch, tmp_path)
    lock_dir = out / "check_tools" / ".lock"
    config = paths.CONFIG
    code = ("import signal, sys, threading\n"
            f"sys.path[:0] = [{str(P)!r}, {str(P / 'tests')!r}]\n"
            "from pathlib import Path\n"
            "import paths\n"
            f"paths.CONFIG = Path({str(config)!r})\n"
            "import check_tools\n"
            "from checks import Check\n"
            "def never(folder):\n"
            "    print('A ROW RAN')\n"
            "    yield\n"
            "threading.Timer(1.0, signal.raise_signal, args=(signal.SIGINT,)).start()\n"
            "print('RETURNED', check_tools.main([], registry={'one': Check('One', never)}))\n")
    with qt.Holder(lock_dir, "check_tools the first") as first:
        done = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, timeout=120)
        assert "RETURNED 130" in done.stdout and "A ROW RAN" not in done.stdout, (done.stdout, done.stderr)
        assert femm_queue.holder(lock_dir)["pid"] == first.pid, "the interrupted run disturbed the holder's lock"
    assert femm_queue.holder(lock_dir) is None


# ----------------------------------------------------------------------------------------------- 2. the known answers
def test_lame_gives_the_ring_its_hoop_stress():
    from checks import ccx_check
    assert ccx_check.lame_hoop(ccx_check.BORE) == pytest.approx(14.0977, abs=5e-4)           # the brief's 14.098 MPa at the bore
    assert ccx_check.lame_hoop(ccx_check.OUTER) == pytest.approx(3.362, abs=5e-4)            # the rim, four times smaller
    assert ccx_check.OMEGA == pytest.approx(1047.1976, abs=1e-4)                              # 10,000 rpm


def test_the_euler_bernoulli_shaft_is_340_hz():
    from checks import ross_check
    assert ross_check.euler_bernoulli_f1(211e9, 7810, 0.015, 0.300) == pytest.approx(340.2, abs=0.1)     # "about 340 Hz" in the brief
    assert ross_check.euler_bernoulli_f1(211e9, 7810, 0.015, 0.150) == pytest.approx(4 * 340.2, abs=0.4)   # half as long: four times as stiff


def test_the_pymoo_case_is_the_one_the_controller_ruled_on_2026_11_06():
    """The brief's 40 x 100 missed its own limit on all of 100 seeds (IGD 0.0514 to 0.1846), so the ruling is 40 individuals x 250
    generations, seed 1, limit 0.05 unchanged (the evidence is in checks/pymoo_check.py). A change here is a new ruling."""
    from checks import pymoo_check
    assert (pymoo_check.POPULATION, pymoo_check.GENERATIONS, pymoo_check.SEED, pymoo_check.LIMIT) == (40, 250, 1, 0.05)


def test_ishigami_first_order_indices():
    from checks import salib_check
    s1, s2, s3 = salib_check.ishigami_first_order(7.0, 0.1)
    assert (s1, s2, s3) == pytest.approx((0.3139, 0.4424, 0.0), abs=5e-5)                    # the brief's 0.314, 0.442, 0.000


def test_ampere_law_for_a_round_conductor():
    from checks import elmer_check
    assert elmer_check.field_of_conductor(100.0, 0.010, 0.010) == pytest.approx(2.0e-3, rel=1e-12)    # at r = R: 2.0 mT
    assert elmer_check.field_of_conductor(100.0, 0.010, 0.005) == pytest.approx(1.0e-3, rel=1e-12)    # at r = R/2: 1.0 mT
    assert elmer_check.field_of_conductor(100.0, 0.010, 0.020) == pytest.approx(1.0e-3, rel=1e-12)    # outside it falls as 1/r
    assert elmer_check.current_density(100.0, 0.010) == pytest.approx(318309.886, rel=1e-8)           # A/m2: 100 A over pi R^2


def test_the_elmer_input_file_carries_the_numbers_it_was_made_from():
    from checks import elmer_check
    sif = elmer_check.case_sif(current=100.0, radius=0.010)
    assert "MagnetoDynamics2D" in sif and "MagnetoDynamicsCalcFields" in sif and "SaveLine" in sif
    assert "Current Density = Real 318309.88" in sif


def test_the_matlab_script_loads_the_plant_without_simulating_it():
    """The brief: load a scratch copy of Improved_plant.slx WITHOUT simulating it, then simulate a tiny model the script builds."""
    from checks import matlab_check
    plant, tiny = matlab_check.PLANT_PART, matlab_check.TINY_PART
    assert "load_system" in plant and "close_system" in plant
    for word in ("sim(", "SimulationCommand", "'update'"):
        assert word not in plant, f"the plant part of the script must not simulate or update the model: {word}"
    assert "sim(" in tiny and "new_system" in tiny
    assert matlab_check.script_text().index(plant.strip()[:30]) < matlab_check.script_text().index(tiny.strip()[:30])


def test_the_two_reference_files_the_matlab_row_copies_are_named_in_the_config():
    plant, params = paths.path("simulink_plant"), paths.path("motor_structure_mat")
    assert plant == P / "reference" / "Improved_plant.slx" and plant.is_file()
    assert params == P / "reference" / "motor structure.mat" and params.is_file()


# ----------------------------------------------------------------------------------------------- 2b. what Gmsh does to the PATH
# Found by the first full run (2026-11-06): after gmsh.initialize() a tool row could not start powershell.exe by name.
# The Python module cuts the PATH of the running PROCESS short (os.environ, Python's own copy, still shows all of it).
# process_path(), cut_process_path() and the path_is_put_back fixture live in tests/conftest.py (the bench's tests use them too).


def test_a_path_that_a_library_cut_short_is_put_back_after_the_block(path_is_put_back):
    whole = process_path()
    with checks.intact_path():
        cut_process_path("C:\\only this")
        assert process_path() == "C:\\only this"
    assert process_path() == whole == path_is_put_back
    with pytest.raises(RuntimeError):                                   # and also when the block ends in an error
        with checks.intact_path():
            cut_process_path("C:\\only that")
            raise RuntimeError("stop")
    assert process_path() == whole


def test_gmsh_initialize_leaves_the_path_whole_inside_intact_path(path_is_put_back):
    """The real cause. Self-check: without intact_path, gmsh 4.15.2 leaves the process PATH short (1,864 characters became
    313 on 2026-11-06). This test does not insist on that (a later Gmsh may behave); it insists that inside intact_path it is whole."""
    import gmsh
    whole = process_path()
    with checks.intact_path():
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.finalize()
    assert process_path() == whole


def test_the_framework_puts_the_path_back_after_every_check(path_is_put_back, monkeypatch, tmp_path, nothing_may_start):
    use_config(monkeypatch, tmp_path)
    whole = process_path()

    def cuts_the_path(folder):
        cut_process_path("C:\\only this")
        yield row()

    assert check_tools.main([], registry={"cuts": Check("Cuts", cuts_the_path)}) == 0
    assert process_path() == whole


# ----------------------------------------------------------------------------------------------- 3a. PrePoMax's settings.bin
def lps(text: str) -> bytes:
    """A length-prefixed string of the .NET BinaryFormatter: the length as a 7-bit-per-byte number, then the UTF-8 bytes."""
    data = text.encode("utf-8")
    n, head = len(data), b""
    while True:
        low = n & 0x7F
        n >>= 7
        head += bytes([low | (0x80 if n else 0)])
        if not n:
            return head + data


def i32(n: int) -> bytes:
    return struct.pack("<i", n)


def settings_stream(exe: str, *, extra: bytes = b"", end: bool = True, library_name: str = "PrePoMax") -> bytes:
    """A small BinaryFormatter stream shaped like PrePoMax's: header, library, one class with a string, an Int32, a Boolean and a
    second string that is only a reference to the first (the form a repeated string takes), then MessageEnd."""
    data = b"\x00" + i32(1) + i32(-1) + i32(1) + i32(0)                           # SerializedStreamHeader
    data += b"\x0c" + i32(2) + lps(library_name)                                  # BinaryLibrary 2
    names = ["_calculixExe", "_calculixDefaultSolverType", "_started", "_again"]
    data += b"\x05" + i32(1) + lps("PrePoMax.SolverSettings") + i32(len(names)) + b"".join(lps(n) for n in names)
    data += bytes([1, 0, 0, 1]) + bytes([8, 1])                                   # types: String, Primitive(Int32), Primitive(Boolean), String
    data += i32(2)                                                                # the library id
    data += b"\x06" + i32(3) + lps(exe)                                           # _calculixExe: a BinaryObjectString, object 3
    data += i32(3)                                                                # _calculixDefaultSolverType = 3
    data += b"\x01"                                                               # _started = True
    data += b"\x09" + i32(3)                                                      # _again: MemberReference to object 3
    data += extra
    return data + (b"\x0b" if end else b"")


def test_the_reader_finds_a_string_member_of_a_class_in_a_binary_formatter_stream():
    from checks import nrbf
    members = nrbf.read_class_members(settings_stream("#\\Solver\\ccx_dynamic.exe"), "PrePoMax.SolverSettings")
    assert members["_calculixExe"] == "#\\Solver\\ccx_dynamic.exe"
    assert members["_calculixDefaultSolverType"] == 3 and members["_started"] is True
    assert members["_again"] == "#\\Solver\\ccx_dynamic.exe"                      # a reference is followed to the string it names


def test_the_reader_reads_a_long_string_and_non_ascii_text():
    from checks import nrbf
    long = "D:\\Apps\\" + "folder with a long name \u00e9\\" * 12 + "ccx_static.exe"            # longer than 127 bytes: a two-byte length
    assert nrbf.read_class_members(settings_stream(long), "PrePoMax.SolverSettings")["_calculixExe"] == long


def test_the_reader_refuses_what_it_cannot_read_with_a_reason():
    """The brief: if the value cannot be read reliably, the row says so. Every way a stream can be wrong is an NrbfError."""
    from checks import nrbf
    good = settings_stream("x")
    for broken, why in ((good[:-1], "ends before MessageEnd"), (good[:40], "cut in the middle"), (good + b"\x00", "bytes after the end"),
                        (b"", "empty"), (b"\x07" + good[1:], "does not start with a header"),
                        (settings_stream("x", end=False) + b"\xff\x0b", "unknown record type")):
        with pytest.raises(nrbf.NrbfError):
            nrbf.read_class_members(broken, "PrePoMax.SolverSettings")
    with pytest.raises(nrbf.NrbfError, match="PrePoMax.Other"):                    # a well-formed stream without the class asked for
        nrbf.read_class_members(good, "PrePoMax.Other")


def test_a_path_inside_the_prepomax_folder_is_stored_with_a_hash_and_is_resolved_against_the_program(tmp_path):
    from checks import prepomax_check
    program = tmp_path / "PrePoMax v3.0.0" / "PrePoMax.exe"
    settings = tmp_path / "settings.bin"
    settings.write_bytes(settings_stream("#\\Solver\\ccx_dynamic.exe"))
    assert prepomax_check.solver_path(settings, program) == program.parent / "Solver" / "ccx_dynamic.exe"
    settings.write_bytes(settings_stream("D:\\somewhere else\\ccx_static.exe"))
    assert prepomax_check.solver_path(settings, program) == Path("D:\\somewhere else\\ccx_static.exe")


def test_a_settings_file_with_no_solver_path_is_an_error_that_says_so(tmp_path):
    from checks import prepomax_check
    settings = tmp_path / "settings.bin"
    settings.write_bytes(settings_stream(""))
    with pytest.raises(prepomax_check.SettingsError, match="_calculixExe"):
        prepomax_check.solver_path(settings, tmp_path / "PrePoMax.exe")
    settings.write_bytes(b"not a formatter stream")
    with pytest.raises(prepomax_check.SettingsError, match="settings.bin"):
        prepomax_check.solver_path(settings, tmp_path / "PrePoMax.exe")


def test_the_reader_reads_the_real_settings_bin_to_the_last_byte():
    """Self-check: the synthetic streams above are my reading of the format; this is PrePoMax's own file. If the reader
    mis-reads any of its 105 objects it raises (it must end exactly at MessageEnd), and if the value were wrong the next test
    of the row (the solver's -v) would show it."""
    from checks import prepomax_check
    try:
        program = paths.tool("prepomax")
    except FileNotFoundError:
        pytest.skip("PrePoMax is not installed here")
    settings = program.parent / "settings.bin"
    if not settings.is_file():
        pytest.skip("PrePoMax has not written its settings.bin yet")
    solver = prepomax_check.solver_path(settings, program)
    assert solver.suffix.lower() == ".exe" and solver.is_file()


# ----------------------------------------------------------------------------------------------- 3b. who is in a MATLAB run's job
def clock(made: dict | None = None):
    """A stand-in for process_created: the creation time of a fake pid is 1000 + pid unless the test sets another in `made`, which
    it may change between two looks (that is how a pid is given to a new process). Returns (made, the function)."""
    made = {} if made is None else made
    return made, (lambda pid: made.get(pid, 1000 + pid))


def test_the_process_list_sees_this_process_and_its_parent():
    from checks import matlab_check
    table = matlab_check.process_table()
    parent, name = table[os.getpid()]                       # the table maps a pid to (parent pid, program name)
    assert name.lower().startswith("python") and parent == os.getppid()


def test_a_process_creation_time_is_read_through_a_handle_that_can_only_read():
    from checks import matlab_check
    created, now = matlab_check.process_created(os.getpid()), matlab_check.filetime_now()
    assert created is not None and 0 <= now - created < 10 ** 7 * 24 * 3600      # this test process was created within the last day
    assert matlab_check.process_created(2 ** 31 - 4) is None                      # no such process: no handle, so no time
    assert matlab_check.PROCESS_QUERY_LIMITED_INFORMATION == 0x1000               # the right that reads a process and cannot end it


def test_the_tracker_follows_the_children_of_the_program_and_nothing_else():
    """me = 10 starts the program 100 (matlab.exe), which starts 101, which starts 102 and then ends, and 103 appears under 102.
    Not members: 50 (was a child of mine before the run), 999 (a stranger), 7 (a service that 101 only talks to)."""
    from checks import matlab_check
    me = 10
    before = {50: (me, "old.exe"), 7: (1, "service.exe"), me: (1, "python.exe")}
    tracker = matlab_check.JobTracker(me, before, created=clock()[1])
    tracker.update(before, now=0.0)
    assert tracker.members() == []
    tracker.update({**before, 100: (me, "matlab.exe"), 999: (5, "stranger.exe")}, now=1.0)
    tracker.update({**before, 100: (me, "matlab.exe"), 101: (100, "MATLAB.exe"), 999: (5, "stranger.exe")}, now=2.0)
    tracker.update({**before, 101: (100, "MATLAB.exe"), 102: (101, "helper.exe"), 999: (5, "stranger.exe")}, now=3.0)       # 100 has ended
    tracker.update({**before, 102: (101, "helper.exe"), 103: (102, "grandchild.exe")}, now=4.0)                             # 101 has ended: orphans
    tracker.update({**before, 103: (102, "grandchild.exe")}, now=5.0)                                                       # 102 has ended too
    got = [(m.pid, m.name, m.first_seen, m.last_seen) for m in tracker.members()]
    assert got == [(100, "matlab.exe", 1.0, 2.0), (101, "MATLAB.exe", 2.0, 3.0), (102, "helper.exe", 3.0, 4.0), (103, "grandchild.exe", 4.0, 5.0)]
    assert [m.ppid for m in tracker.members()] == [me, 100, 101, 102]
    assert tracker.last_seen("MATLAB.exe") == 3.0 and tracker.last_seen("nothing.exe") is None
    assert tracker.still_present(now=5.0) == ["grandchild.exe"]


def test_the_tracker_never_counts_a_process_that_existed_before_the_run_even_if_it_has_the_right_parent():
    from checks import matlab_check
    tracker = matlab_check.JobTracker(10, {20: (10, "earlier child.exe")}, created=clock()[1])
    tracker.update({20: (10, "earlier child.exe"), 21: (20, "its child.exe")}, now=1.0)
    assert tracker.members() == []


def test_a_new_process_of_the_job_that_gets_a_pid_from_the_baseline_is_still_a_member():
    """Reuse 1 (review M1): pid 50 was an old child of mine and has ended; Windows gives 50 to the program's new helper. By pid
    alone it was taken for the old one and left out."""
    from checks import matlab_check
    me = 10
    made, created = clock({50: 500})
    tracker = matlab_check.JobTracker(me, {50: (me, "old.exe"), me: (1, "python.exe")}, created=created)
    tracker.update({50: (me, "old.exe")}, now=0.0, wall=600)
    assert tracker.members() == []                                        # the old one stays out
    made[50] = 2000                                                       # the old one ended; the new process is created at 2000
    tracker.update({50: (me, "helper.exe")}, now=1.0, wall=2100)
    assert [(m.pid, m.name) for m in tracker.members()] == [(50, "helper.exe")]


def test_a_strangers_child_is_not_a_member_because_its_parent_has_the_reused_pid_of_a_member_that_ended():
    """Reuse 2 (M1): the program 100 ended; a stranger got the pid 100 and started 300. By pid alone, 300 was a child of a member."""
    from checks import matlab_check
    me = 10
    made, created = clock({100: 1000})
    tracker = matlab_check.JobTracker(me, {me: (1, "python.exe")}, created=created)
    tracker.update({100: (me, "matlab.exe")}, now=1.0, wall=1100)
    tracker.update({}, now=2.0, wall=2000)                                # 100 is gone (found gone at wall 2000)
    made.update({100: 3000, 300: 3100})                                   # a stranger takes pid 100 and starts 300, both after 2000
    tracker.update({100: (5, "stranger.exe"), 300: (100, "child of the stranger.exe")}, now=3.0, wall=3200)
    assert [m.name for m in tracker.members()] == ["matlab.exe"]


def test_a_stranger_on_the_reused_pid_of_a_member_does_not_extend_its_last_seen():
    """Reuse 3 (M1): a desktop MATLAB.exe gets the pid of the job's MATLAB.exe after it ended: 'empty N s after the last matlab.exe'
    would have been wrong by pid and name alone."""
    from checks import matlab_check
    me = 10
    made, created = clock({101: 1000})
    tracker = matlab_check.JobTracker(me, {me: (1, "python.exe")}, created=created)
    tracker.update({101: (me, "MATLAB.exe")}, now=1.0, wall=1100)
    tracker.update({101: (me, "MATLAB.exe")}, now=2.0, wall=1200)
    made[101] = 5000                                                      # the member ended; somebody's own MATLAB.exe is created with its pid
    tracker.update({101: (5, "MATLAB.exe")}, now=3.0, wall=5100)
    tracker.update({101: (5, "MATLAB.exe")}, now=4.0, wall=5200)
    [member] = tracker.members()
    assert (member.pid, member.first_seen, member.last_seen) == (101, 1.0, 2.0)


def test_a_child_born_while_its_parent_lived_is_a_member_even_if_first_seen_after_the_parent_ended():
    from checks import matlab_check
    me = 10
    made, created = clock({100: 1000, 101: 1500, 102: 2500})
    tracker = matlab_check.JobTracker(me, {me: (1, "python.exe")}, created=created)
    tracker.update({100: (me, "launcher.exe")}, now=1.0, wall=1100)
    tracker.update({101: (100, "late child.exe"), 102: (100, "too late.exe")}, now=2.0, wall=2000)    # 100 ended at wall 2000
    assert [m.name for m in tracker.members()] == ["launcher.exe", "late child.exe"]      # 101 was created at 1500, 102 at 2500: after the end


def test_a_process_whose_creation_time_cannot_be_read_is_a_member_that_says_so():
    from checks import matlab_check
    me = 10
    tracker = matlab_check.JobTracker(me, {me: (1, "python.exe")}, created=lambda pid: None)       # every handle refused
    tracker.update({100: (me, "protected.exe")}, now=1.0)
    tracker.update({100: (me, "protected.exe")}, now=2.0)
    [member] = tracker.members()
    assert member.created is None and member.last_seen == 2.0                                      # known by pid and name
    assert "creation time unreadable" in tracker.summary(returned_at=3.0, program="protected.exe")


def test_the_job_summary_names_the_members_and_how_long_the_job_took_to_empty():
    from checks import matlab_check
    tracker = matlab_check.JobTracker(10, {}, created=clock()[1])
    tracker.update({100: (10, "matlab.exe")}, now=0.1)
    tracker.update({100: (10, "matlab.exe"), 101: (100, "MATLAB.exe"), 102: (101, "MathWorksServiceHost.exe")}, now=0.2)
    tracker.update({101: (100, "MATLAB.exe"), 102: (101, "MathWorksServiceHost.exe")}, now=60.0)
    tracker.update({102: (101, "MathWorksServiceHost.exe")}, now=61.0)
    text = tracker.summary(returned_at=62.5, program="matlab.exe")        # Windows names are not case-sensitive: the launcher and the engine are both "matlab.exe"
    assert "matlab.exe" in text and "MATLAB.exe" in text and "MathWorksServiceHost.exe" in text
    assert "2.50 s later" in text                           # the last matlab.exe was seen at 60.0; run_solver returned at 62.5 (the job was empty)
    assert "after MATLAB: MathWorksServiceHost.exe" in text  # the helper that was still there when MATLAB's own processes were gone


def test_the_tracker_tells_which_members_are_alive_by_creation_time():
    from checks import matlab_check
    made, created = clock()
    tracker = matlab_check.JobTracker(10, {}, created=created)
    tracker.update({100: (10, "matlab.exe"), 101: (100, "MATLAB.exe")}, now=1.0)
    assert [m.pid for m in tracker.alive({100: (10, "matlab.exe"), 101: (100, "MATLAB.exe")})] == [100, 101]
    made[101] = 99999                                                     # pid 101 is now another process
    assert [m.pid for m in tracker.alive({100: (10, "matlab.exe"), 101: (100, "MATLAB.exe")})] == [100]


def test_a_process_list_that_cannot_be_read_is_reported_and_not_taken_for_an_empty_job(monkeypatch):
    """Review M2: a persistent OSError was swallowed one look at a time and the summary said 'the program ended before the first look'."""
    from checks import matlab_check
    monkeypatch.setattr(matlab_check, "process_table", lambda: (_ for _ in ()).throw(OSError("the list is gone")))
    tracker = matlab_check.JobTracker(10, {}, created=clock()[1])
    watcher = matlab_check.Watcher(tracker, time.perf_counter())
    watcher.start()
    time.sleep(0.3)
    watcher.stop()
    assert watcher.looks >= 2 and watcher.failed == watcher.looks and "the list is gone" in str(watcher.last_error) and watcher.error is None
    text = tracker.summary(1.0, looks=watcher.looks, failed=watcher.failed, last_error=watcher.last_error)
    assert "UNKNOWN" in text and "FAILED" in text and "the list is gone" in text and "ended before the first look" not in text
    some = tracker.summary(1.0, looks=10, failed=2)                       # a few failed looks beside a job that was seen: said, not hidden
    assert "2 of 10 looks at the process list FAILED" in some and "may be incomplete" in some
    quiet = tracker.summary(1.0, looks=10, failed=0)
    assert "ended before the first look" in quiet and "FAILED" not in quiet


def test_a_scratch_copy_that_vanishes_during_the_matlab_run_fails_the_row_and_says_which(monkeypatch, tmp_path):
    """Review M9: a copy that was gone counted as 'unchanged'. The run is faked (run_solver is replaced), so no MATLAB starts."""
    files = {name: tmp_path / name for name in ("matlab.exe", "plant.slx", "params.mat")}
    for file in files.values():
        file.write_bytes(b"x" + file.name.encode())
    out = use_config(monkeypatch, tmp_path, f"matlab = '{files['matlab.exe']}'",
                     f"simulink_plant = '{files['plant.slx']}'\nmotor_structure_mat = '{files['params.mat']}'")
    deleted = []

    def fake_run_solver(command, *, label, cwd, timeout_s, threads=None, **rest):
        folder = Path(cwd)
        (folder / "result.txt").write_text("OK=1\nPLANT_LOADED=1\nRELEASE=2025b\nTINY_ERR_GAIN=0.000e+00\nMAXNUMCOMPTHREADS=1\n")
        for name in deleted:
            (folder / name).unlink()
        return femm_queue.SolverRun(0, "", "", 1.5)

    monkeypatch.setattr(femm_queue, "run_solver", fake_run_solver)
    assert check_tools.main(["--only", "matlab"]) == 0                    # the control: everything in place, the faked run passes
    assert "references unchanged" in read_csv(out)[1][3]
    deleted.append("plant.slx")
    assert check_tools.main(["--only", "matlab"]) == 1
    line = read_csv(out)[1]
    assert line[4] == FAIL and "a scratch copy is MISSING (plant.slx)" in line[3]


# ----------------------------------------------------------------------------------------------- 3b. MATLAB's threads (final review: B1)
def matlab_row_setup(monkeypatch, tmp_path, result_text: str):
    """The MATLAB row with dummy files and a run_solver that only writes `result_text` into result.txt and notes the command. Returns (out, seen)."""
    files = {name: tmp_path / name for name in ("matlab.exe", "plant.slx", "params.mat")}
    for file in files.values():
        file.write_bytes(b"x" + file.name.encode())
    out = use_config(monkeypatch, tmp_path, f"matlab = '{files['matlab.exe']}'",
                     f"simulink_plant = '{files['plant.slx']}'\nmotor_structure_mat = '{files['params.mat']}'")
    seen = {}

    def fake_run_solver(command, *, label, cwd, timeout_s, threads=None, **rest):
        seen.update(command=[str(part) for part in command], threads=threads)
        (Path(cwd) / "result.txt").write_text(result_text)
        return femm_queue.SolverRun(0, "", "", 1.5)

    monkeypatch.setattr(femm_queue, "run_solver", fake_run_solver)
    return out, seen


GOOD_RESULT = "OK=1\nPLANT_LOADED=1\nRELEASE=2025b\nTINY_ERR_GAIN=0.000e+00\nMAXNUMCOMPTHREADS=1\n"


def test_matlab_is_started_with_one_computational_thread_and_the_row_says_how_many_it_reports(monkeypatch, tmp_path):
    """B1 of the final review. R10 caps solver threads at 8 and OMP_NUM_THREADS does not cap MATLAB's own computational threads (its default is the number of physical cores, 12
    here). The command carries -singleCompThread, as B8's does, and the script writes maxNumCompThreads into result.txt; the row shows it."""
    from checks import matlab_check
    out, seen = matlab_row_setup(monkeypatch, tmp_path, GOOD_RESULT)
    assert check_tools.main(["--only", "matlab"]) == 0
    assert "-singleCompThread" in seen["command"] and seen["command"].index("-singleCompThread") < seen["command"].index("-batch"), seen["command"]
    line = read_csv(out)[1]
    assert line[4] == PASS and "1 computational thread" in line[3], line[3]
    assert "maxNumCompThreads" in matlab_check.script_text() and "MAXNUMCOMPTHREADS" in matlab_check.script_text()


def test_the_matlab_row_fails_when_matlab_reports_more_threads_than_the_cap_or_none(monkeypatch, tmp_path):
    out, _ = matlab_row_setup(monkeypatch, tmp_path, GOOD_RESULT.replace("MAXNUMCOMPTHREADS=1", "MAXNUMCOMPTHREADS=12"))
    assert check_tools.main(["--only", "matlab"]) == 1
    line = read_csv(out)[1]
    assert line[4] == FAIL and "12 computational threads" in line[3] and "at most 8" in line[3], line[3]
    out, _ = matlab_row_setup(monkeypatch, tmp_path, GOOD_RESULT.replace("MAXNUMCOMPTHREADS=1\n", ""))
    assert check_tools.main(["--only", "matlab"]) == 1
    line = read_csv(out)[1]
    assert line[4] == FAIL and "did not report" in line[3], line[3]


# ----------------------------------------------------------------------------------------------- 3c. the rows that cannot start
def test_every_way_settings_bin_can_be_wrong_is_a_failed_row_with_its_reason(monkeypatch, tmp_path, nothing_may_start):
    """Review M8: the PrePoMax FAIL path, through check_tools itself and not only through solver_path."""
    program = tmp_path / "PrePoMax v3" / "PrePoMax.exe"
    program.parent.mkdir()
    program.write_bytes(b"not a program")                                # a file that is there, with no version resource
    out = use_config(monkeypatch, tmp_path, f"prepomax = '{program}'")
    settings = program.parent / "settings.bin"
    good = settings_stream("#\\Solver\\ccx_dynamic.exe")
    inside = program.parent / "Solver" / "ccx_dynamic.exe"
    variants = [
        ("not a formatter stream at all", b"not a formatter stream", "cannot tell: cannot read settings.bin"),
        ("cut short by one byte", good[:-1], "cannot tell: cannot read settings.bin"),
        ("a stray byte after the end", good + b"\x00", "cannot tell: cannot read settings.bin"),
        ("no solver path in it", settings_stream(""), "cannot tell: settings.bin has no _calculixExe"),
        ("a solver that is not there", settings_stream("Z:\\gone\\ccx.exe"), "missing: the solver of PrePoMax's _calculixExe (Z:\\gone\\ccx.exe)"),
        ("a path inside PrePoMax that is not there", good, f"missing: the solver of PrePoMax's _calculixExe ({inside})"),
    ]
    for what, data, reason in variants:
        settings.write_bytes(data)
        assert check_tools.main(["--only", "prepomax"]) == 1, what
        version, solver = read_csv(out)[1:3]
        assert version[4] == FAIL and "(none)" in version[3], what             # the dummy program states no version: a failed row, too
        assert solver[4] == FAIL and reason in solver[3], f"{what}: {solver[3]}"
    settings.unlink()
    assert check_tools.main(["--only", "prepomax"]) == 1
    assert read_csv(out)[2][3] == f"missing: PrePoMax settings.bin ({settings})"


# ----------------------------------------------------------------------------------------------- 4. R10: only run_solver starts a program
STARTER_NAMES = {"subprocess": {"Popen", "run", "call", "check_call", "check_output", "getoutput", "getstatusoutput"},
                 "os": {"system", "popen", "startfile"}, "multiprocessing": {"Process"}}
STARTER_PREFIXES = {"os": ("spawn", "exec"), "asyncio": ("create_subprocess_",)}
STARTER_MODULES = ("subprocess", "os", "asyncio", "multiprocessing")
ENDERS = {"kill", "killpg", "terminate", "TerminateProcess", "TerminateJobObject"}


def is_starter(module: str, name: str) -> bool:
    return name in STARTER_NAMES.get(module, ()) or name.startswith(STARTER_PREFIXES.get(module, ()))


class ScopeWalker(ast.NodeVisitor):
    """Walks a module and keeps the name of the function or class it is in (Class.method), as the scope of what it finds."""

    def __init__(self):
        self.scope: list[str] = []

    def here(self) -> str:
        return ".".join(self.scope) or "<module>"

    def visit_ClassDef(self, node):
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    visit_FunctionDef = visit_AsyncFunctionDef = visit_ClassDef


def programs_started(source: str) -> set[tuple[str, str]]:
    """(scope, what) for every way the source starts a program or a process other than through run_solver: calls and uses of
    subprocess.* / os.system, os.popen, os.spawn*, os.exec*, asyncio.create_subprocess_*, multiprocessing.Process, under whatever
    name the module was imported (`import subprocess as sp`), `from subprocess import run` (the import itself is found), a star
    import from one of those modules, and getattr(<one of those modules>, ...)."""
    bound: dict[str, str] = {}
    found: set[tuple[str, str]] = set()

    class Walker(ScopeWalker):
        def visit_Import(self, node):
            for alias in node.names:
                top = alias.name.split(".")[0]
                if top in STARTER_MODULES:
                    bound[alias.asname or top] = top

        def visit_ImportFrom(self, node):
            module = (node.module or "").split(".")[0]
            if module in STARTER_MODULES:
                for alias in node.names:
                    if alias.name == "*" or is_starter(module, alias.name):
                        found.add((self.here(), f"from {module} import {alias.name}"))

        def visit_Attribute(self, node):
            if isinstance(node.value, ast.Name) and node.value.id in bound and is_starter(bound[node.value.id], node.attr):
                found.add((self.here(), f"{bound[node.value.id]}.{node.attr}"))
            self.generic_visit(node)

        def visit_Call(self, node):
            if isinstance(node.func, ast.Name) and node.func.id == "getattr" and node.args and isinstance(node.args[0], ast.Name) \
                    and node.args[0].id in bound:
                found.add((self.here(), f"getattr({bound[node.args[0].id]}, ...)"))
            self.generic_visit(node)

    Walker().visit(ast.parse(source))
    return found


def scopes_calling(source: str, owner: str | None, name: str) -> list[str]:
    """The scopes that call `owner.name(...)` (or `name(...)` when owner is None), sorted, once for each call."""
    found: list[str] = []

    class Walker(ScopeWalker):
        def visit_Call(self, node):
            func = node.func
            if (owner is None and isinstance(func, ast.Name) and func.id == name) or (
                    owner is not None and isinstance(func, ast.Attribute) and func.attr == name
                    and isinstance(func.value, ast.Name) and func.value.id == owner):
                found.append(self.here())
            self.generic_visit(node)

    Walker().visit(ast.parse(source))
    return sorted(found)


def process_enders(source: str) -> set[tuple[str, str]]:
    """What in the source could end a process: `taskkill` in any string that is not a docstring, and the names kill, killpg, terminate,
    TerminateProcess and TerminateJobObject (as an attribute or a name). Only run_solver's job ends processes."""
    tree = ast.parse(source)
    docstrings = {id(node.body[0].value) for node in ast.walk(tree) if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
                  and node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant)
                  and isinstance(node.body[0].value.value, str)}
    found: set[tuple[str, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings and "taskkill" in node.value.lower():
            found.add(("string", "taskkill"))
        if isinstance(node, ast.Attribute) and node.attr in ENDERS:
            found.add(("attribute", node.attr))
        if isinstance(node, ast.Name) and node.id in ENDERS:
            found.add(("name", node.id))
    return found


def test_the_scan_for_programs_started_outside_run_solver_finds_them_however_they_are_written():
    """Self-check: the scan has to find each of these before 'it found nothing' means anything (review M4)."""
    samples = {
        "import subprocess\nsubprocess.run(['x'])\n": {("<module>", "subprocess.run")},
        "import subprocess as sp\ndef f():\n    sp.Popen('x')\n": {("f", "subprocess.Popen")},
        "from subprocess import check_output as co\nco('x')\n": {("<module>", "from subprocess import check_output")},
        "from subprocess import Popen, PIPE\n": {("<module>", "from subprocess import Popen")},
        "from subprocess import *\n": {("<module>", "from subprocess import *")},
        "import os\nos.system('y')\nos.spawnlp(0, 'x')\nos.execvp('x', [])\nos.popen('z')\n":
            {("<module>", "os.system"), ("<module>", "os.spawnlp"), ("<module>", "os.execvp"), ("<module>", "os.popen")},
        "import asyncio\nasync def g():\n    await asyncio.create_subprocess_exec('x')\n": {("g", "asyncio.create_subprocess_exec")},
        "import subprocess\nrun = getattr(subprocess, 'run')\n": {("<module>", "getattr(subprocess, ...)")},
        "import multiprocessing as mp\nclass C:\n    def m(self):\n        mp.Process(target=f)\n": {("C.m", "multiprocessing.Process")},
        "import subprocess\nsubprocess.PIPE\nsubprocess.TimeoutExpired\nsubprocess.CREATE_NO_WINDOW\n": set(),     # constants and exceptions start nothing
        "import os\nos.path.join('a', 'b')\nos.environ\nos.getpid()\n": set(),
    }
    for source, expected in samples.items():
        assert programs_started(source) == expected, source


def test_the_scan_for_something_that_ends_a_process_finds_it():
    assert process_enders("import subprocess\nsubprocess.run(['taskkill', '/F'])\n") == {("string", "taskkill")}
    assert process_enders("import os\nos.kill(1, 9)\nos.killpg(1, 9)\n") == {("attribute", "kill"), ("attribute", "killpg")}
    assert process_enders("proc.terminate()\nproc.kill()\nctypes.windll.kernel32.TerminateProcess(h, 1)\n") == {
        ("attribute", "terminate"), ("attribute", "kill"), ("attribute", "TerminateProcess")}
    assert process_enders('def f():\n    """no taskkill here: a docstring says so."""\n    return 1\n') == set()


MODULES_THAT_RUN_A_PROGRAM = ("femm_check", "gmsh_check", "ccx_check", "prepomax_check", "ross_check", "elmer_check", "matlab_check")


def test_nothing_in_checks_starts_a_program_at_all():
    """R10: one heavy run at a time, below normal priority, under a timeout, in a job that dies with its caller. A case that called
    subprocess itself would be none of these, and there is no exception: no module of checks/ starts a program or a process in any
    way but through femm_queue.run_solver. check_tools.py is scanned too."""
    offenders = []
    for file in sorted((P / "checks").glob("*.py")) + [P / "check_tools.py"]:
        found = programs_started(file.read_text(encoding="utf-8"))
        if found:
            offenders.append(f"{file.name}: starts {sorted(found)}")
    assert not offenders, "a program is started without run_solver:\n" + "\n".join(offenders)


def test_nothing_in_checks_ends_a_process_by_its_number_or_at_all():
    """Fix 1: whatever a run started is in run_solver's job, and the job ends it. No taskkill, no os.kill, no Popen.kill/terminate."""
    offenders = []
    for file in sorted((P / "checks").glob("*.py")) + [P / "check_tools.py"]:
        found = process_enders(file.read_text(encoding="utf-8"))
        if found:
            offenders.append(f"{file.name}: {sorted(found)}")
    assert not offenders, "something ends a process itself:\n" + "\n".join(offenders)


def test_the_modules_that_run_a_program_call_run_solver():
    for module in MODULES_THAT_RUN_A_PROGRAM:
        source = (P / "checks" / f"{module}.py").read_text(encoding="utf-8")
        assert scopes_calling(source, "femm_queue", "run_solver"), f"{module} runs a program, so it has to call femm_queue.run_solver"


def test_matlab_check_uses_the_one_sha256_helper_of_the_package(tmp_path):
    """Review M3: the helper was written twice, verbatim; it is written once now, in checks/__init__.py."""
    import hashlib
    from checks import matlab_check
    file = tmp_path / "x.bin"
    file.write_bytes(b"abc" * 400_000)
    assert checks.sha256(file) == hashlib.sha256(b"abc" * 400_000).hexdigest()
    assert matlab_check.sha256 is checks.sha256


def test_the_docstring_lists_every_exit_status():
    """Review M10: 130 (Ctrl-C) was missing."""
    assert "130" in check_tools.__doc__
    for text in ("Exit status 0", "1 if one failed", "2 for a wrong --only"):
        assert text in check_tools.__doc__
    assert "130" in (P / "README.md").read_text(encoding="utf-8")


def test_motorgeo_to_gmsh_puts_the_process_path_back_by_itself(path_is_put_back, tmp_path):
    """159 of the final review. gmsh.initialize() cuts this process's own PATH short and finalize() does not put it back (the check above shows it); the guard used to sit around
    the callers (check_tools, the bench), so any other caller (check_geometry.py, example.py, learn/tools_figures.py) cut the PATH of its own process. motorgeo.to_gmsh saves
    the PATH before gmsh.initialize() and restores it after finalize(), in a finally. The PATH read from Windows, not Python's copy (which Gmsh never cuts), is the same before and after."""
    import motorgeo as mg
    p = mg.Params(barrier="none")
    regions = mg.build(p, gap="sliding")
    before = process_path()
    nodes, groups = mg.to_gmsh(p, regions, str(tmp_path / "t.msh"), size_max=3.0)
    assert nodes > 100 and groups, "the mesh was not made"
    assert process_path() == before, f"to_gmsh left the process's PATH at {len(process_path())} characters, from {len(before)}"


def test_motorgeo_to_gmsh_puts_the_process_path_back_when_it_fails_too(path_is_put_back, tmp_path):
    """The same in the finally: a call that raises after gmsh.initialize() (no regions at all) has cut the PATH, and it comes back. (Gmsh is finalized by the test: a failed
    to_gmsh leaves it initialized, as it always did.)"""
    import gmsh

    import motorgeo as mg
    p = mg.Params(barrier="none")
    before = process_path()
    try:
        with pytest.raises(IndexError):
            mg.to_gmsh(p, [], str(tmp_path / "t.msh"))
        assert process_path() == before
    finally:
        if gmsh.isInitialized():
            gmsh.finalize()
        os.environ["PATH"] = before
