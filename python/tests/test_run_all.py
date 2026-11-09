"""Tests for validation/run_all.py: the runner that re-runs every script on FEMM 5.0.

The runner is only useful if it can be trusted for a long run without anybody watching, so these tests check the parts
that a long run depends on, without running a single study and without starting FEMM:

1. The table of steps: unique names, every script is really there, every "needs" points at an EARLIER step, every
   step is listed in the module's docstring under its number, and the order puts barrier_study before reference_open.
2. The command line: --list, --only, --from, unknown names, the time limits.
3. The table of results (run_all.csv): rows are replaced, kept, renumbered and dropped as the docstring says, and a
   run marks its steps `not-run` first, so a run that stops early cannot leave an old "ok" behind.
4. Process identity and which FEMM programs count as "left running", on made-up process lists (a pool worker's FEMM,
   an orphan, a FEMM of somebody's own window, a FEMM that was there before, an id that was handed to another process).
5. Real child processes, with tiny scripts written to a temp folder: a good step, a failing one, a step that prints
   non-ASCII, a step that runs too long (everything seen below it must go), a FEMM the step leaves behind (ended, if it
   was seen below the step; only reported, if not), an unrelated orphan FEMM that must survive a timeout, and one that
   was already there. The "FEMM" is a copy of cmd.exe under a name of its own (the runner is told to count that name as
   a FEMM program), so these tests never touch a real FEMM and a real runner active at the same time cannot count it.
6. The job object, the crash-window setting that children inherit, and the lock.
7. main() with a fake step runner: order, exit codes, the console in cp1252, Ctrl-C.

Run from P:   uv run pytest tests/test_run_all.py -q       (about 30 s; nothing here starts FEMM)
"""
import ast
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

import paths
import run_all as ra

P = Path(__file__).resolve().parent.parent
CMD_EXE = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "cmd.exe"


def end_tree(pid):
    """Test clean-up only: end a process and everything below it with Windows' own taskkill."""
    exe = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "taskkill.exe"
    subprocess.run([str(exe), "/PID", str(pid), "/F", "/T"], capture_output=True, check=False)


def end_tree_if_same(pid, identity):
    """Test clean-up for a process this test did not start itself and holds no handle to (its number came out of a log): end it, and everything below it, only if it is still
    the process that was seen then, the same program with the same creation time (ra.process_identity). `identity` is None when nothing was seen: nothing is ended."""
    if identity is not None and ra.process_identity(pid) == identity:
        end_tree(pid)


def first_number(text):
    """The first line of a log that is just a number (the step scripts print the pid they want the test to look at)."""
    return int(next(ln for ln in text.splitlines() if ln.strip().isdigit()))


# ----------------------------------------------------------------------------- 1. the table of steps
def test_step_names_are_unique_and_safe_in_a_file_name():
    names = [s.name for s in ra.STEPS]
    assert len(names) == len(set(names))
    for name in names:
        assert name.replace("_", "").isalnum(), name


def test_every_script_of_the_table_exists_where_the_step_runs_it():
    for s in ra.STEPS:
        assert (P / s.folder).is_dir(), s.name
        assert (P / s.folder / s.script).is_file(), f"{s.name}: {s.folder}/{s.script}"


def test_every_input_comes_from_an_earlier_step():
    position = {s.name: i for i, s in enumerate(ra.STEPS)}
    for s in ra.STEPS:
        for provider, what in s.needs:
            assert provider in position, f"{s.name} needs {what!r} from {provider!r}, which is no step"
            assert position[provider] < position[s.name], f"{s.name} needs {what!r} from {provider}, which runs later"


def test_the_inputs_the_controller_named_are_in_the_table():
    needs = {s.name: {p for p, _ in s.needs} for s in ra.STEPS}
    assert "barrier_study" in needs["materials_fields"]                       # fields.py reads barrier_study's .ans files
    assert {"materials_fields", "materials_current"} <= needs["materials_magnetic"]
    assert {"barrier_study", "reference_open", "flux_explorer_sweep"} <= needs["pf_harmonics"]
    assert "validate_fix_zeros" in needs["barrier_plots"] or "validate_full" in needs["barrier_plots"]
    assert needs["hand_reference_check"] == {"validate_full"}
    assert [s.name for s in ra.STEPS].index("hand_reference_check") == [s.name for s in ra.STEPS].index("validate_full") + 1


def test_every_step_is_explained_in_the_docstring():
    """The docstring's section "What each step proves" lists every step as "<its number> <its name>".

    The number is the step's position in the table, which is also the NN of its log file (NN-name.log). An earlier version
    of this test looked for a stem such as "current" or "thermal" anywhere in the whole docstring, and those words occur
    all over it, so the six materials_* steps passed without being listed at all. Now only that section is read, the
    whole name must stand there as a word of its own ("check_femm" is not found inside "check_femm_fine"), and the number
    in front of it must be the step's position.
    """
    doc = ra.__doc__
    assert "What each step proves" in doc and "Self-check" in doc, "the docstring's sections have been renamed"
    section = doc.split("What each step proves", 1)[1].split("Self-check", 1)[0]
    for number, s in enumerate(ra.STEPS, start=1):
        entry = rf"(?<![\w-]){number}\s+{re.escape(s.name)}(?!\w)"
        assert re.search(entry, section), f"step {number} ({s.name}) is not listed as '{number} {s.name}' in the docstring"


def test_the_materials_steps_follow_build_md_and_the_learn_steps_run_in_learn():
    names = [s.name for s in ra.STEPS]
    first = names.index("materials_fields")
    assert names[first:first + 6] == ["materials_fields", "materials_current", "materials_magnetic",
                                      "materials_thermal", "materials_structural", "materials_crosssection"]
    for s in ra.STEPS:
        assert (s.folder == "learn") == s.name.startswith("learn_"), s.name
    figure_steps = [s for s in ra.STEPS if s.name in ("learn_geometry_figure", "learn_toy_front_figure")]
    for s in figure_steps:
        assert s.mkdirs == (r"out\learn",)                                    # the figure's folder must exist
        assert s.args[0].startswith(r"..\out\learn" + "\\")                   # and is inside P\out, never the project notes


def test_no_argument_points_into_the_docs_figures_folder():
    """The runner keeps the materials figures out of the repository's docs figures folder (MOTOR_FIGURES_DIR is out/materials/figures), so no step
    is handed that folder, or a path inside it, either. An argument counts as a path read from the step's own folder, and from P."""
    docs = Path(os.path.normpath(paths.path("docs_figures")))
    for s in ra.STEPS:
        for base in (P / s.folder, P):
            for arg in s.args:
                target = Path(os.path.normpath(base / arg))
                assert docs != target and docs not in target.parents, s.name


def test_barrier_study_runs_before_reference_open_which_needs_its_folder():
    names = [s.name for s in ra.STEPS]
    assert names.index("barrier_study") < names.index("reference_open")
    assert "fem" in dict(ra.STEPS[names.index("reference_open")].needs)["barrier_study"]


# ----------------------------------------------------------------------------- 2. the command line
def test_list_prints_every_step_and_what_it_needs(capsys):
    assert ra.main(["--list"]) == 0
    out = capsys.readouterr().out
    for s in ra.STEPS:
        assert s.name in out and ra.command_text(s) in out
        for provider, what in s.needs:
            assert what in out and provider in out
    assert "FEMM 4.2" in out and "FEMM 5.0" in out and "time limit" in out


def test_list_respects_only_and_from(capsys):
    ra.main(["--list", "--only", "check_femm"])
    lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.strip()]
    assert len(lines) == 2 and "check_femm " in lines[1] and "check_femm_fine" not in lines[1]


def test_only_runs_in_table_order_and_takes_exact_names():
    chosen = ra.select_steps(["check_femm,check_geometry", "check_vs_fem"], None)
    assert [s.name for s in chosen] == ["check_geometry", "check_vs_fem", "check_femm"]
    assert [s.name for s in ra.select_steps(["check_femm"], None)] == ["check_femm"]       # not check_femm_fine


def test_from_runs_that_step_and_all_later_ones():
    chosen = ra.select_steps(None, "learn_05_pymoo_toy")
    assert [s.name for s in chosen] == ["learn_05_pymoo_toy", "learn_geometry_figure", "learn_toy_front_figure"]
    assert len(ra.select_steps(None, None)) == len(ra.STEPS)


def test_unknown_names_are_refused_with_a_hint(capsys):
    with pytest.raises(ValueError, match="check_femm"):
        ra.select_steps(["check_fem"], None)
    with pytest.raises(SystemExit) as stop:
        ra.main(["--only", "nonsense"])
    assert stop.value.code == 2
    assert "nonsense" in capsys.readouterr().err
    with pytest.raises(SystemExit) as stop:
        ra.main(["--only", "example", "--from", "evaluate"])
    assert stop.value.code == 2


def test_every_step_has_the_time_it_took_in_the_first_femm_5_sweep():
    for s in ra.STEPS:
        assert isinstance(s.measured, int) and s.measured >= 1, s.name


def test_time_limits_are_five_times_the_measured_time_but_at_least_five_minutes():
    by_name = {s.name: s for s in ra.STEPS}
    assert ra.timeout_for(by_name["check_femm"]) == 300                       # 5 x 2 s is under the minimum
    assert ra.timeout_for(by_name["validate_full"]) == 300                    # 5 x 27 s is too
    assert ra.timeout_for(by_name["barrier_study"]) == 5 * 131                # 11 minutes, not the 3 hours of 3 x FEMM 4.2
    assert ra.timeout_for(by_name["barrier_study"], factor=10) == 10 * 131
    assert ra.timeout_for(by_name["barrier_study"], factor=0) == 300          # and never below the minimum
    assert max(ra.timeout_for(s) for s in ra.STEPS) < 15 * 60                 # nothing in this table may hide a hang for hours


def test_the_runner_imports_no_study_module():
    """A study that is broken must not be able to break the runner: it only launches them."""
    tree = ast.parse(Path(ra.__file__).read_text(encoding="utf-8"))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0:
            imported.add(node.module.split(".")[0])
    studies = {p.stem for folder in (P, P / "materials", P / "learn", P / "validation") for p in folder.glob("*.py")}
    assert imported & studies == {"femm_session"}                             # the one module of ours it uses, for femm_dir()


def test_the_children_get_the_environment_the_docstring_promises(monkeypatch):
    """Needs FEMM 5.0 to be installed (femm_dir() finds it), but starts no FEMM."""
    monkeypatch.delenv("FEMM_DIR", raising=False)                             # a stale terminal: the registry must answer
    try:
        env = ra.child_environment()
    except FileNotFoundError:
        pytest.skip("FEMM 5.0 is not installed here")
    assert (Path(env["FEMM_DIR"]) / "femm_headlessw.exe").is_file()
    assert Path(env["MOTOR_FIGURES_DIR"]) == P / "out" / "materials" / "figures"
    assert env["PYTHONUNBUFFERED"] == "1" and env["PYTHONIOENCODING"] == "utf-8"


def test_durations_are_written_for_people():
    assert ra.fmt_duration(1.8) == "1.8 s" and ra.fmt_duration(2) == "2 s" and ra.fmt_duration(None) == "-"
    assert ra.fmt_duration(300) == "5 min 00 s" and ra.fmt_duration(3 * 3600 + 120) == "3 h 02 min"


# ----------------------------------------------------------------------------- 3. the table of results
def result(name, status="ok", code=0, left=0):
    step = next(s for s in ra.STEPS if s.name == name)
    return ra.Result(1, step, code, 1.25, left, status, Path("x.log"), "2026-11-04 20:00:00")


def rows_of(*results):
    return [ra.ledger_row(r) for r in results]


def read_ledger(path):
    import csv
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_ledger_has_the_columns_the_docstring_describes(tmp_path):
    path = tmp_path / "run_all.csv"
    ra.update_ledger(path, rows_of(result("check_femm")))
    rows = read_ledger(path)
    assert list(rows[0]) == list(ra.COLUMNS)
    assert list(ra.COLUMNS)[:6] == ["step", "name", "command", "exit_code", "seconds", "femm_left_running"]
    assert rows[0]["step"] == "3" and rows[0]["name"] == "check_femm" and rows[0]["command"] == "python check_femm.py"
    assert rows[0]["exit_code"] == "0" and rows[0]["femm_left_running"] == "0" and rows[0]["status"] == "ok"


def test_a_partial_run_replaces_only_its_own_rows_and_keeps_table_order(tmp_path):
    path = tmp_path / "run_all.csv"
    ra.update_ledger(path, rows_of(result("example"), result("check_geometry", "failed", 1)))
    ra.update_ledger(path, rows_of(result("check_geometry")))                 # the fixed step, run alone
    rows = read_ledger(path)
    assert [(r["name"], r["status"]) for r in rows] == [("check_geometry", "ok"), ("example", "ok")]
    assert [r["step"] for r in rows] == ["1", "5"]


def test_the_ledger_drops_steps_that_are_gone_and_survives_a_broken_file(tmp_path):
    path = tmp_path / "run_all.csv"
    path.write_text("step,name\n1,a_step_that_was_renamed\n2,check_femm\n", encoding="utf-8")
    ra.update_ledger(path, rows_of(result("example")))
    assert [r["name"] for r in read_ledger(path)] == ["check_femm", "example"]
    path.write_bytes(b"\xff\xfe not a csv at all \x00")
    ra.update_ledger(path, rows_of(result("example")))
    assert [r["name"] for r in read_ledger(path)] == ["example"]


def test_a_pending_row_says_not_run_and_has_nothing_that_looks_like_a_result(tmp_path):
    path = tmp_path / "run_all.csv"
    ra.update_ledger(path, rows_of(result("example")))                        # an old "ok"
    ra.update_ledger(path, [ra.pending_row(next(s for s in ra.STEPS if s.name == "example"), 5)])
    row = read_ledger(path)[0]
    assert row["status"] == "not-run" and row["exit_code"] == "" and row["seconds"] == "" and row["finished"] == ""
    assert row["femm_left_running"] == "" and row["command"] == "python example.py"


# ----------------------------------------------------------------------------- 4. identity, and which FEMMs count
#   process table: {pid: (parent pid, program)}
TABLE = {
    10: (1, "python.exe"),            # the step's python.exe (launcher)
    11: (10, "python.exe"),           # its interpreter
    12: (11, "python.exe"),           # a pool worker
    13: (12, "femm_headlessw.exe"),   # the worker's FEMM
    14: (13, "fkn.exe"),              # ... and its solver
    20: (19, "femm_headlessw.exe"),   # an orphan: its parent (19) has gone
    30: (31, "fkn.exe"),              # a solver of somebody's own FEMM window ...
    31: (2, "femm_gui.exe"),          # ... whose parent is running outside the step
    40: (41, "femm_headlessw.exe"),   # a FEMM that was there before the step
    41: (2, "explorer.exe"),
    2: (0, "explorer.exe"),
    1: (0, "uv.exe"),
}
BORN = {pid: 1000 + pid for pid in TABLE}                                     # creation times: a made-up, consistent order


def identity_of(pid):
    return (TABLE[pid][1], BORN[pid]) if pid in TABLE else None


def seen(*pids):
    """`own` as the runner builds it: the identities of the processes seen below the step."""
    return {pid: identity_of(pid) for pid in pids}


def test_a_femm_seen_below_the_step_is_ours_and_the_others_are_only_reported():
    own = seen(10, 11, 12, 13)                                                # the solver 14 was not seen (born between looks)
    ours, others = ra.leaked_femms(set(), TABLE, own, identity_of)
    assert set(ours) == {13}                                                  # seen, still the same process: the only one to end
    assert others == {14, 20}                                                 # the unseen solver below the step; the orphan
    # 30 is a solver of somebody's own window (its parent 31 is alive and outside the step) and 40 hangs below explorer.exe:
    # neither is ever counted or touched, even though nobody listed them as "before".


def test_a_femm_that_was_there_before_is_never_counted_even_below_the_step():
    ours, others = ra.leaked_femms({13, 14, 20}, TABLE, seen(10, 11, 12, 13, 14), identity_of)
    assert ours == {} and others == set()


def test_an_id_that_was_handed_to_another_process_is_not_the_one_that_was_seen():
    reused = {**seen(10, 11, 12, 13), 13: ("femm_headlessw.exe", 5)}          # seen as one process; now another has the id
    ours, others = ra.leaked_femms(set(), TABLE, reused, identity_of)
    assert 13 not in ours                                                     # never ended: it is not the process that was seen
    assert 13 in others                                                       # but it hangs below the worker, so it is reported


def test_a_worker_that_outlived_its_parent_still_owns_its_femm():
    table = {12: (11, "python.exe"), 13: (12, "femm_headlessw.exe"), 99: (0, "x.exe")}    # 11, the interpreter, has gone
    ident = lambda pid: (table[pid][1], 1000 + pid) if pid in table else None
    ours, others = ra.leaked_femms(set(), table, {12: ident(12)}, ident)
    assert ours == {} and others == {13}                                      # below a seen worker: reported (never seen itself)
    ours, others = ra.leaked_femms(set(), table, {12: ident(12), 13: ident(13)}, ident)
    assert set(ours) == {13} and others == set()
    ours, others = ra.leaked_femms(set(), table, {}, ident)                   # nothing seen: 13's parent 12 is alive, outside
    assert ours == {} and others == set()


def test_descendants_follow_the_tree_and_survive_a_loop():
    assert ra.descendants(TABLE, 10) == {11, 12, 13, 14}
    assert ra.descendants({1: (2, "a.exe"), 2: (1, "b.exe")}, 1) == {2}


def test_a_child_older_than_the_parent_it_names_is_a_stale_link():
    """Windows keeps naming a parent that has ended, and gives its id to a new process. Such a link must not be followed."""
    table = {10: (1, "a.exe"), 11: (10, "b.exe"), 12: (11, "c.exe")}
    born = {10: 100, 11: 200, 12: 150}                                        # 12 names parent 11 but is older than it
    assert ra.descendants(table, 10, born.get) == {11}
    assert ra.descendants(table, 10) == {11, 12}                              # without creation times the link is believed
    assert ra.descendants(table, 10, lambda pid: None) == set()               # a process that cannot be asked is not claimed


def test_the_process_list_sees_this_python():
    table = ra.running_processes()
    assert table[os.getpid()][1].startswith("python")
    assert ra.femm_pids(table) == {p for p, (_, name) in table.items() if name in ra.FEMM_IMAGES}


def test_a_process_identity_is_its_program_and_creation_time():
    name, born = ra.process_identity(os.getpid())
    assert name.startswith("python") and born > 0
    assert ra.process_identity(os.getpid()) == (name, born)                   # the same twice
    assert ra.process_identity(0x7FFFFFF0) is None                            # no such process


def test_a_process_is_ended_only_if_it_is_still_the_one_that_was_seen():
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        time.sleep(0.5)
        name, born = ra.process_identity(sleeper.pid)
        assert not ra.terminate_if_same(sleeper.pid, (name, born + 1))        # another creation time: the id was reused
        assert not ra.terminate_if_same(sleeper.pid, ("fkn.exe", born))       # another program
        assert not ra.terminate_if_same(0x7FFFFFF0, (name, born))             # no such process
        time.sleep(0.3)
        assert sleeper.poll() is None                                         # none of that touched it
        assert ra.terminate_if_same(sleeper.pid, (name, born))
        assert sleeper.wait(timeout=10) is not None
    finally:
        end_tree(sleeper.pid)


# ----------------------------------------------------------------------------- 5. real child processes
@pytest.fixture
def child_env():
    return {**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"}


def run(tmp_path, body, child_env, timeout=60, grace=1.0, name="tiny"):
    """Write a script, run it as a step, return the Result."""
    script = tmp_path / f"{name}.py"
    script.write_text(body, encoding="utf-8")
    step = ra.Step(name, str(script))                                         # absolute: the folder is not used
    return ra.run_step(step, 7, child_env, tmp_path / "logs", timeout, grace=grace)


def test_a_good_step_is_ok_and_its_log_is_utf8_with_a_header_and_a_footer(tmp_path, child_env):
    r = run(tmp_path, "print('torque 0.81 N\u00b7m, \u03c8 12 mWb, 30\u00b0')\n", child_env)
    assert (r.status, r.exit_code, r.femm_left) == ("ok", 0, 0)
    assert r.log.name == "07-tiny.log"
    text = r.log.read_bytes().decode("utf-8")                                 # strict: must be valid UTF-8
    assert "torque 0.81 N\u00b7m, \u03c8 12 mWb, 30\u00b0" in text
    assert text.startswith("# run_all step 7: tiny\n# command: python ")
    assert "exit code 0" in text.splitlines()[-1] and "status ok" in text.splitlines()[-1]


def test_a_failing_step_keeps_its_exit_code_and_its_traceback(tmp_path, child_env):
    r = run(tmp_path, "import sys\nprint('before')\nraise SystemExit(3)\n", child_env)
    assert (r.status, r.exit_code) == ("failed", 3)
    r = run(tmp_path, "raise ValueError('boom')\n", child_env, name="boom")
    assert (r.status, r.exit_code) == ("failed", 1)
    assert "ValueError: boom" in r.log.read_text(encoding="utf-8")            # stderr goes to the same log


def test_a_step_run_again_keeps_the_log_of_the_run_before(tmp_path, child_env):
    first = run(tmp_path, "print('the first run')\n", child_env)
    second = run(tmp_path, "print('the second run')\n", child_env)
    assert first.log == second.log
    assert "the second run" in second.log.read_text(encoding="utf-8")
    previous = second.log.with_name(second.log.name + ".prev")
    assert "the first run" in previous.read_text(encoding="utf-8") and "the second" not in previous.read_text(encoding="utf-8")
    run(tmp_path, "print('the third run')\n", child_env)
    assert "the second run" in previous.read_text(encoding="utf-8")           # one generation: the older one goes


def test_the_log_grows_while_the_step_runs(tmp_path, child_env):
    script = tmp_path / "slow.py"
    script.write_text("import time\nprint('first line', flush=True)\ntime.sleep(4)\nprint('second line')\n")
    step = ra.Step("slow", str(script))
    box = []
    worker = threading.Thread(target=lambda: box.append(ra.run_step(step, 8, child_env, tmp_path / "logs", 60)))
    worker.start()
    log = tmp_path / "logs" / "08-slow.log"
    deadline = time.time() + 3
    seen_line = False
    while time.time() < deadline and not seen_line:
        seen_line = log.exists() and "first line" in log.read_text(encoding="utf-8", errors="replace")
        time.sleep(0.1)
    still_running = worker.is_alive()
    worker.join()
    assert seen_line and still_running                                        # read it while the step was still going
    assert box[0].status == "ok"


def test_a_step_over_its_time_limit_is_ended_with_everything_seen_below_it(tmp_path, child_env):
    body = ("import subprocess, sys, time\n"
            "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)'])\n"
            "print(child.pid, flush=True)\n"
            "time.sleep(120)\n")
    started = time.time()
    r = run(tmp_path, body, child_env, timeout=2)
    assert r.status == "timeout" and r.exit_code != 0
    assert time.time() - started < 40
    grandchild = first_number(r.log.read_text(encoding="utf-8"))
    assert grandchild not in ra.running_processes()                           # the whole tree went
    assert "status timeout" in r.log.read_text(encoding="utf-8").splitlines()[-1]


@pytest.fixture
def fake_femm(tmp_path, monkeypatch):
    """A program with a name of its own (a copy of cmd.exe) that the runner is told to count as a FEMM program.

    Not fkn.exe: a real runner active at the same time (a sweep in the background) would count and might end this one.
    Returns a function that starts one; it also has .exe and .name.
    """
    if not CMD_EXE.is_file():
        pytest.skip("no cmd.exe to copy")
    name = f"fakefemm_runall_test_{os.getpid()}.exe"
    exe = tmp_path / "fakefemm" / name
    exe.parent.mkdir()
    shutil.copy(CMD_EXE, exe)
    monkeypatch.setattr(ra, "FEMM_IMAGES", (name,))
    started = []                                # the Popen objects, not the numbers: while one is alive its handle is open, and the number is not given to another process

    def start():
        proc = subprocess.Popen([str(exe), "/c", "ping -n 120 127.0.0.1 >nul"], stdin=subprocess.DEVNULL,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        started.append(proc)
        time.sleep(0.5)
        assert ra.running_processes().get(proc.pid, (0, ""))[1] == name, "the fake FEMM did not start"
        return proc.pid
    start.exe, start.name, start.started = exe, name, started
    yield start
    for proc in started:
        end_tree(proc.pid)                      # certainly the process this fixture started: its handle is still open
        try:
            proc.wait(10)
        except subprocess.TimeoutExpired:
            proc.kill()


def fake_launching_script(fake_femm, after, detached=False):
    """A step script: start a fake FEMM, print its pid, then do `after` (python lines)."""
    flags = ", creationflags=0x00000008 | 0x00000200" if detached else ""            # DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
    return ("import subprocess, time\n"
            f"p = subprocess.Popen([r'{fake_femm.exe}', '/c', 'ping -n 120 127.0.0.1 >nul'], stdin=subprocess.DEVNULL,\n"
            f"    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL{flags})\n"
            "print(p.pid, flush=True)\n" + after)


def test_a_femm_seen_below_the_step_that_outlives_it_is_counted_and_ended(tmp_path, child_env, fake_femm):
    already_there = fake_femm()
    body = fake_launching_script(fake_femm, "time.sleep(1.6)\n")              # long enough to be seen; then it exits, leaving it
    r = run(tmp_path, body, child_env, grace=1.0)
    leaked = first_number(r.log.read_text(encoding="utf-8"))
    assert (r.status, r.exit_code, r.femm_left) == ("leak", 0, 1)
    assert f"still running: {fake_femm.name} (pid {leaked}), ended by the runner" in r.log.read_text(encoding="utf-8")
    assert leaked not in ra.running_processes()                               # one leak must not disturb the next step
    assert already_there in ra.running_processes()                            # the one that was there before: never touched


def test_a_femm_the_runner_never_saw_is_counted_and_reported_but_left_alone(tmp_path, child_env, fake_femm, monkeypatch):
    monkeypatch.setattr(ra, "look_below", lambda own, table, root: None)     # the runner sees nothing below the step
    already_there = fake_femm()
    r = run(tmp_path, fake_launching_script(fake_femm, "", detached=True), child_env, grace=1.0)
    leaked = first_number(r.log.read_text(encoding="utf-8"))
    seen = ra.process_identity(leaked)                                        # what the number was when it was read: this test holds no handle to it
    try:
        assert seen is not None and seen[0] == fake_femm.name
        assert (r.status, r.exit_code, r.femm_left) == ("leak", 0, 1)         # an orphan: counted, so the run does not pass
        assert f"still running: {fake_femm.name} (pid {leaked}), not seen below the step: left alone" in \
            r.log.read_text(encoding="utf-8")
        assert leaked in ra.running_processes() and already_there in ra.running_processes()
    finally:
        end_tree_if_same(leaked, seen)


def test_a_timeout_ends_the_femm_of_the_step_but_not_an_unrelated_orphan_or_the_one_before(tmp_path, child_env, fake_femm):
    """The case of somebody else's FEMM from a debugging session that has since exited, born while the step ran."""
    already_there = fake_femm()
    unrelated = {}

    def somebody_elses_session():
        """A python that starts a FEMM program detached and exits at once: that FEMM is an orphan, and not the step's."""
        code = (f"import subprocess; p = subprocess.Popen([r'{fake_femm.exe}', '/c', 'ping -n 120 127.0.0.1 >nul'], "
                "stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, creationflags=0x8 | 0x200); "
                "print(p.pid, flush=True)")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=30)
        unrelated["pid"] = int(out.stdout.split()[0])
    timer = threading.Timer(1.2, somebody_elses_session)                      # while the step is running
    timer.start()
    body = fake_launching_script(fake_femm, "time.sleep(120)\n")
    r = run(tmp_path, body, child_env, timeout=4, grace=1.0)
    timer.join()
    own = first_number(r.log.read_text(encoding="utf-8"))
    seen = ra.process_identity(unrelated["pid"]) if "pid" in unrelated else None      # what the orphan's number was when it was read (no handle to it)
    try:
        table = ra.running_processes()
        assert r.status == "timeout"
        assert own not in table                                               # the step's own FEMM was ended with it
        assert already_there in table                                         # the one that existed before was not
        assert unrelated["pid"] in table                                      # and neither was the orphan that is not the step's
        assert r.femm_left == 1                                               # ... which is reported, so nobody overlooks it
        assert "not seen below the step: left alone" in r.log.read_text(encoding="utf-8")
    finally:
        if "pid" in unrelated:
            end_tree_if_same(unrelated["pid"], seen)


def test_ctrl_c_ends_the_step_and_is_recorded_as_interrupted(tmp_path, child_env, monkeypatch):
    real_wait = subprocess.Popen.wait
    calls = []

    def interrupted_once(self, timeout=None):
        calls.append(timeout)
        if len(calls) == 2:
            raise KeyboardInterrupt                   # Ctrl-C arrives while the runner waits, half a second into the step
        return real_wait(self, timeout)
    monkeypatch.setattr(subprocess.Popen, "wait", interrupted_once)
    r = run(tmp_path, "import os, time\nprint(os.getpid(), flush=True)\ntime.sleep(120)\n", child_env)
    assert (r.status, r.femm_left) == ("interrupted", 0)
    assert first_number(r.log.read_text(encoding="utf-8")) not in ra.running_processes()


def test_ctrl_c_while_the_runner_waits_for_femm_to_go_still_gives_the_step_its_row(tmp_path, child_env, monkeypatch):
    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt
    monkeypatch.setattr(ra, "leaked_femms", interrupted)                      # the wait after the step: Ctrl-C arrives here
    r = run(tmp_path, "print('done')\n", child_env)
    assert (r.status, r.exit_code) == ("interrupted", 0)
    assert "status interrupted" in r.log.read_text(encoding="utf-8").splitlines()[-1]


def test_a_fault_in_the_runner_while_a_step_runs_ends_the_step_and_is_raised(tmp_path, child_env, monkeypatch):
    real = ra.look_below
    calls = []

    def bug_once(own, table, root):
        calls.append(root)
        if len(calls) == 1:
            raise RuntimeError("a bug in the runner")
        return real(own, table, root)
    monkeypatch.setattr(ra, "look_below", bug_once)
    with pytest.raises(RuntimeError, match="a bug in the runner"):
        run(tmp_path, "import os, time\nprint(os.getpid(), flush=True)\ntime.sleep(120)\n", child_env)
    log = tmp_path / "logs" / "07-tiny.log"
    assert first_number(log.read_text(encoding="utf-8")) not in ra.running_processes()    # never left behind a crashed runner


# ----------------------------------------------------------------------------- 7. main() with a fake step runner
@pytest.fixture
def fake_run(tmp_path, monkeypatch):
    """main() with its files in tmp_path and run_step replaced by a table of outcomes. Returns (calls, outcomes)."""
    import femm_session
    monkeypatch.setattr(ra, "LEDGER", tmp_path / "run_all.csv")
    monkeypatch.setattr(ra, "LOGS", tmp_path / "logs")
    monkeypatch.setattr(ra, "LOCK", tmp_path / "run_all.lock")
    monkeypatch.setattr(ra, "child_environment", lambda: {"FEMM_DIR": "x", "MOTOR_FIGURES_DIR": "y"})
    monkeypatch.setattr(ra, "join_kill_on_exit_job", lambda: None)
    monkeypatch.setattr(ra, "suppress_crash_dialogs", lambda: None)           # this is the test process: leave it as it is
    monkeypatch.setattr(femm_session, "lower_priority", lambda: None)
    calls, outcomes = [], {}

    def fake(step, number, env, logs, timeout, grace=ra.GRACE):
        calls.append(step.name)
        what = outcomes.get(step.name, ("ok", 0, 0))
        if isinstance(what, Exception):
            raise what
        status, code, left = what[:3]
        log = logs / f"{number:02d}-{step.name}.log"
        if len(what) > 3:                                                     # a log with this text in it
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(what[3], encoding="utf-8")
        return ra.Result(number, step, code, 0.1, left, status, log, "2026-11-04 20:00:00")
    monkeypatch.setattr(ra, "run_step", fake)
    return calls, outcomes


def test_main_runs_the_steps_in_table_order_records_each_and_exits_1_if_one_failed(fake_run, tmp_path, capsys):
    calls, outcomes = fake_run
    outcomes["check_vs_fem"] = ("failed", 1, 0)
    code = ra.main(["--only", "example,check_vs_fem,check_geometry"])
    out = capsys.readouterr().out
    assert calls == ["check_geometry", "check_vs_fem", "example"]
    assert code == 1
    assert [(r["name"], r["status"]) for r in read_ledger(tmp_path / "run_all.csv")] == [
        ("check_geometry", "ok"), ("check_vs_fem", "failed"), ("example", "ok")]
    assert "NOT OK: step 2 check_vs_fem: failed, exit code 1" in out
    assert "3 step(s) run, 2 ok, 1 not ok." in out
    assert "3 of 33 steps recorded: 2 ok, 1 not ok, 0 not run." in out
    assert not (tmp_path / "run_all.lock").exists()                           # the lock is always given back


def test_main_exits_0_when_every_step_is_ok(fake_run, capsys):
    assert ra.main(["--only", "check_geometry,check_femm"]) == 0
    assert "2 step(s) run, 2 ok, 0 not ok." in capsys.readouterr().out


def test_a_full_run_that_passes_leaves_every_step_ok(fake_run, tmp_path, capsys):
    calls, outcomes = fake_run
    assert ra.main([]) == 0
    assert calls == [s.name for s in ra.STEPS]
    assert {r["status"] for r in read_ledger(tmp_path / "run_all.csv")} == {"ok"}
    assert "33 of 33 steps recorded: 33 ok, 0 not ok, 0 not run." in capsys.readouterr().out


def test_a_leak_makes_the_run_fail_even_though_the_script_exited_0(fake_run, tmp_path, capsys):
    calls, outcomes = fake_run
    outcomes["check_femm"] = ("leak", 0, 2)
    assert ra.main(["--only", "check_femm"]) == 1
    row = read_ledger(tmp_path / "run_all.csv")[0]
    assert (row["exit_code"], row["femm_left_running"], row["status"]) == ("0", "2", "leak")
    assert "2 FEMM program(s) left running" in capsys.readouterr().out


def test_the_summary_blames_the_step_whose_files_a_failed_step_needed(fake_run, capsys):
    calls, outcomes = fake_run
    outcomes["barrier_study"] = ("failed", 1, 0)
    outcomes["barrier_plots"] = ("failed", 1, 0)
    ra.main(["--only", "barrier_study,barrier_plots,optimize_demo"])
    out = capsys.readouterr().out
    assert "its input comes from barrier_study, which also failed" in out
    assert calls == ["barrier_study", "barrier_plots", "optimize_demo"]       # the run goes on after a failure


def test_a_bug_in_the_runner_is_recorded_and_the_next_step_still_runs(fake_run, tmp_path, capsys):
    calls, outcomes = fake_run
    outcomes["check_vs_fem"] = RuntimeError("bug")
    code = ra.main(["--only", "check_geometry,check_vs_fem,check_femm"])
    assert code == 1 and calls == ["check_geometry", "check_vs_fem", "check_femm"]
    rows = {r["name"]: r["status"] for r in read_ledger(tmp_path / "run_all.csv")}
    assert rows == {"check_geometry": "ok", "check_vs_fem": "runner-error", "check_femm": "ok"}


def test_ctrl_c_stops_the_run_and_the_steps_not_run_say_not_run_not_an_old_ok(fake_run, tmp_path, capsys):
    calls, outcomes = fake_run
    path = tmp_path / "run_all.csv"
    ra.update_ledger(path, rows_of(result("check_geometry"), result("check_vs_fem"), result("check_femm")))   # an old run: all ok
    outcomes["check_vs_fem"] = ("interrupted", 1, 0)
    assert ra.main(["--only", "check_geometry,check_vs_fem,check_femm"]) == 1
    assert calls == ["check_geometry", "check_vs_fem"]                        # check_femm was not started
    rows = {r["name"]: r["status"] for r in read_ledger(path)}
    assert rows == {"check_geometry": "ok", "check_vs_fem": "interrupted", "check_femm": "not-run"}
    out = capsys.readouterr().out
    assert "1 not run (the run was stopped)" in out and "3 of 33 steps recorded: 1 ok, 1 not ok, 1 not run." in out


def test_a_full_run_stopped_half_way_does_not_claim_the_rest_is_ok(fake_run, tmp_path, capsys):
    calls, outcomes = fake_run
    assert ra.main([]) == 0                                                   # an older, complete run: 33 ok
    outcomes["example"] = ("interrupted", 1, 0)
    ra.main([])
    rows = read_ledger(tmp_path / "run_all.csv")
    by_status = {}
    for r in rows:
        by_status.setdefault(r["status"], []).append(r["name"])
    assert by_status["interrupted"] == ["example"] and len(by_status["ok"]) == 4
    assert len(by_status["not-run"]) == len(ra.STEPS) - 5                     # everything after the interruption
    assert "33 of 33 steps recorded: 4 ok, 1 not ok, 28 not run." in capsys.readouterr().out


def test_a_console_that_cannot_show_a_character_does_not_end_the_run(fake_run, monkeypatch):
    """A redirected stdout is often cp1252. A failed step's last log lines are printed; one letter outside cp1252 used to
    raise UnicodeEncodeError inside print_summary, which runs in a finally."""
    import io
    calls, outcomes = fake_run
    outcomes["check_femm"] = ("failed", 1, 0, "the solver said \u03c8 = 12 mWb \u2192 \u00b0\n")    # psi, an arrow: not in cp1252
    raw = io.BytesIO()
    console = io.TextIOWrapper(raw, encoding="cp1252", errors="strict", write_through=True)
    monkeypatch.setattr(sys, "stdout", console)
    assert ra.main(["--only", "check_femm"]) == 1                             # no UnicodeEncodeError
    console.flush()
    text = raw.getvalue().decode("cp1252")
    assert "NOT OK: step 3 check_femm" in text and "the solver said ? = 12 mWb ? \u00b0" in text


def test_a_second_runner_is_refused_and_does_not_touch_the_first_ones_lock(fake_run, tmp_path, capsys):
    calls, outcomes = fake_run
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        time.sleep(0.5)
        (tmp_path / "run_all.lock").write_text(str(sleeper.pid))
        assert ra.main(["--only", "check_geometry"]) == 2
        assert calls == [] and (tmp_path / "run_all.lock").read_text() == str(sleeper.pid)
        assert "another run_all.py is running" in capsys.readouterr().out
    finally:
        end_tree(sleeper.pid)
        sleeper.wait()


def test_without_femm_dir_the_runner_says_so_and_gives_the_lock_back(fake_run, tmp_path, monkeypatch, capsys):
    calls, outcomes = fake_run

    def no_femm():
        raise FileNotFoundError("FEMM_DIR is not set")
    monkeypatch.setattr(ra, "child_environment", no_femm)
    assert ra.main(["--only", "check_geometry"]) == 2
    assert "cannot start: FEMM_DIR is not set" in capsys.readouterr().out
    assert calls == [] and not (tmp_path / "run_all.lock").exists()


# ----------------------------------------------------------------------------- 6. the job object, the crash window, the lock
def run_python(tmp_path, source):
    """Run a script with run_all importable; return its output."""
    script = tmp_path / "helper.py"
    script.write_text(f"import sys\nsys.path.insert(0, r'{Path(ra.__file__).parent}')\nimport run_all\n" + source,
                      encoding="utf-8")
    out = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    return out.stdout


def test_a_process_that_ends_without_cleaning_up_takes_its_children_with_it(tmp_path):
    out = run_python(tmp_path, (
        "import os, subprocess\n"
        "run_all.join_kill_on_exit_job()\n"
        "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)'], stdin=subprocess.DEVNULL,\n"
        "                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)\n"
        "print(child.pid, flush=True)\n"
        "os._exit(0)                       # no cleanup of any kind\n"))
    pid = int(out.split()[0])
    deadline = time.time() + 10
    while time.time() < deadline and pid in ra.running_processes():
        time.sleep(0.2)
    assert pid not in ra.running_processes()


def test_the_children_inherit_the_switched_off_crash_window(tmp_path):
    out = run_python(tmp_path, (
        "import ctypes, subprocess\n"
        "get = ctypes.windll.kernel32.GetErrorMode\n"
        "before = get()\n"
        "run_all.suppress_crash_dialogs()\n"
        "child = subprocess.run([sys.executable, '-c', 'import ctypes; print(ctypes.windll.kernel32.GetErrorMode())'],\n"
        "                       capture_output=True, text=True)\n"
        "print(before, get(), child.stdout.strip())\n"))
    before, mine, child = (int(v) for v in out.split())
    wanted = ra._SEM_FAILCRITICALERRORS | ra._SEM_NOGPFAULTERRORBOX
    assert mine & wanted == wanted and child & wanted == wanted               # this process, and the one it started
    assert mine & before == before                                            # nothing that was set before was taken away


def test_the_lock_refuses_a_second_runner_and_replaces_a_stale_lock(tmp_path):
    lock = tmp_path / "run_all.lock"
    assert ra.take_lock(lock) is None and lock.read_text() == str(os.getpid())
    ra.drop_lock(lock)
    assert not lock.exists()
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        time.sleep(0.5)
        lock.write_text(str(sleeper.pid))                                     # a live python holds it
        holder = ra.take_lock(lock)
        assert holder is not None and holder != os.getpid()
        assert lock.read_text() == str(sleeper.pid)                           # left as it was
    finally:
        end_tree(sleeper.pid)                                                 # the launcher and the interpreter below it
        sleeper.wait()
    assert ra.take_lock(lock) is None                                         # now the holder has gone: stale, replaced
    ra.drop_lock(lock)


def test_a_config_toml_that_is_not_valid_is_cannot_start_with_exit_2_and_not_a_traceback(fake_run, tmp_path, monkeypatch, capsys):
    """B4 of the final review (Task 2's minor): child_environment() asks femm_session.femm_dir(), which reads config.toml; a file that is not valid TOML raises ValueError
    (paths.py names the file and the line), and main() caught only FileNotFoundError and ImportError, so the runner ended with a traceback instead of 'cannot start'."""
    calls, outcomes = fake_run

    def broken_config():
        raise ValueError(f"{tmp_path / 'config.toml'} is not valid TOML: Invalid statement (at line 3, column 5)")
    monkeypatch.setattr(ra, "child_environment", broken_config)
    assert ra.main(["--only", "check_geometry"]) == 2
    shown = capsys.readouterr().out
    assert "cannot start: " in shown and "config.toml is not valid TOML" in shown and "Traceback" not in shown
    assert calls == [] and not (tmp_path / "run_all.lock").exists()


def test_child_environment_raises_a_value_error_for_a_config_that_is_not_toml(monkeypatch, tmp_path):
    """The premise of the test above, with the real child_environment and a real broken file."""
    bad = tmp_path / "config.toml"
    bad.write_text("[tools\nfemm = \n", encoding="utf-8")
    monkeypatch.setattr(paths, "CONFIG", bad)
    monkeypatch.delenv("FEMM_DIR", raising=False)
    with pytest.raises(ValueError, match="config.toml"):
        ra.child_environment()


# ----------------------------------------------------------------------------- 131: a test ends only a process it started, or one it has seen
def test_the_fake_femm_fixture_holds_a_handle_on_every_program_it_starts(fake_femm):
    """131 of the final review. The fixture used to keep the bare pids of the programs it started and end them with taskkill at the end; a pid names a process only while it
    lives, and Windows may give a dead process's number to a stranger. A Popen object holds the process's handle, and while a handle is open the number is not given away: the
    fixture keeps the Popen objects, so what it ends in the clean-up is certainly what it started."""
    first = fake_femm()
    assert [process.pid for process in fake_femm.started] == [first]
    assert all(isinstance(process, subprocess.Popen) for process in fake_femm.started)
    second = fake_femm()
    assert [process.pid for process in fake_femm.started] == [first, second]


def test_a_process_the_test_holds_no_handle_to_is_ended_only_if_it_is_still_the_one_that_was_seen(tmp_path):
    """131: the two tests that read a FEMM's pid from a step's log (nothing but a number) end it only if ra.process_identity (program and creation time) is still what it was when
    the number was read. A number that has been given to another process since is left alone."""
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        time.sleep(0.5)
        seen = ra.process_identity(sleeper.pid)
        name, born = seen
        end_tree_if_same(sleeper.pid, (name, born + 1))                          # another creation time: the number was reused
        end_tree_if_same(sleeper.pid, ("fkn.exe", born))                         # another program
        end_tree_if_same(sleeper.pid, None)                                      # nothing was seen: nothing is ended
        time.sleep(0.3)
        assert sleeper.poll() is None                                            # none of that touched it
        end_tree_if_same(sleeper.pid, seen)
        assert sleeper.wait(timeout=10) is not None
    finally:
        sleeper.kill()
        sleeper.wait()
