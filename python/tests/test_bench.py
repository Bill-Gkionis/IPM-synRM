"""Tests for the bench's shared part (bench/__init__.py, bench/run_bench.py) and for the scans that read every module of bench/. No solver.

Run from P:   uv run pytest tests/test_bench.py -q      (a few seconds; no FEMM, no Elmer, no solver slot is taken)
All of them:  uv run pytest tests/test_bench*.py -q      (about 30 s: Gmsh runs, in this process)

The bench has one test file for its shared part and one for each case (tests/test_bench_b1.py ... test_bench_b9.py). The cases
themselves (`uv run python bench/run_bench.py --only B1,B2,B3`) run the real tools on the motor and compare
their numbers with the expectations of the brief. A test that did the same would only repeat them. These tests pin what can go wrong around
them. Here, in groups:
  1. The shared part (bench/__init__.py): where a case writes, what may be emptied, how a number is held against its expectation, whether a
     check gates the verdict (a check that does not is reported, not hidden; one that cannot be made is skipped), the one copy of the motor's
     parameters and of the field picture (its title wrapped), the JSON that carries the results (numpy numbers included), the priority classes
     the log must name, the watcher that reads them.
  2. R10 for the bench: nothing in bench/ starts a program except through femm_queue, nothing ends a process, and every call has its time
     limit, however the call is spelled (the scan is tried on small sources first).
  3. The PATH that Gmsh cuts: a scan that every call that starts Gmsh is inside `with intact_path()`. (Each function that starts Gmsh has a test
     of its own, in the case's file, that reads the PATH of the process from Windows.)
  4. run_bench.py: --only, the exit codes (a check that does not gate leaves the status alone, one that gates makes it 1), the counts in
     summary.json, a case that raises, a case folder that is emptied only for the case that runs.
  5. The bench's own process: lowered, and the record of it that B9 reads, with the threads read back from Gmsh.
  6. The roots and the lock (Task 13, B9's rulings): `run_bench.py --out ROOT`, the root every case writes under and fresh() obeys, the lock a run holds on its
     root for the whole run (a second run on the same root waits and names the holder), and the watcher's peak count.
"""
import ast
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import paths
from bench_test_tools import P, the_test_run_is_not_lowered, use_config      # noqa: F401  (the fixture is autouse)


# ----------------------------------------------------------------------------------------------- 1. the shared part: the shape of a check, gating, the results file, the priority classes, the log
def test_a_case_writes_under_out_bench_and_the_folder_is_made(monkeypatch, tmp_path):
    import bench
    out = use_config(monkeypatch, tmp_path)
    assert bench.root() == out / "bench"
    folder = bench.case_folder("B1")
    assert folder == out / "bench" / "B1" and folder.is_dir()


def test_only_a_case_folder_inside_out_bench_may_be_emptied(monkeypatch, tmp_path):
    import bench
    out = use_config(monkeypatch, tmp_path)
    folder = bench.case_folder("B2")
    (folder / "sub").mkdir()
    (folder / "sub" / "old.txt").write_text("old")
    (folder / "figure.png").write_bytes(b"png")
    assert bench.fresh(folder) == folder
    assert folder.is_dir() and list(folder.iterdir()) == []
    # the root itself, a folder outside it, a name that is not a case, and a path that climbs out with ".." are all refused
    for bad in (bench.root(), out, tmp_path, bench.root() / "notes", bench.root() / "B1" / ".." / ".." / "B1", bench.root() / "B1" / "fem"):
        with pytest.raises(ValueError):
            bench.fresh(bad)
    assert folder.is_dir()


def test_a_number_held_against_its_expectation_says_how_far_off_it_is_and_where_the_expectation_comes_from():
    import bench
    ok = bench.compare("psi_pm", 16.65, 16.61, 0.5, "README row", unit="mWb")
    assert ok["verdict"] == "PASS" and ok["name"] == "psi_pm" and ok["source"] == "README row" and ok["unit"] == "mWb"
    assert ok["difference_percent"] == pytest.approx(100 * (16.65 / 16.61 - 1))
    bad = bench.compare("T", 1.14, 1.127, 0.5, "README row")
    assert bad["verdict"] == "FAIL" and bad["difference_percent"] == pytest.approx(100 * (1.14 / 1.127 - 1))
    edge = bench.compare("x", 1.005, 1.0, 0.5, "s")                      # exactly on the tolerance passes: "within 0.5 %"
    assert edge["verdict"] == "PASS"
    assert bench.compare("neg", -1.0, -1.0, 0.1, "s")["verdict"] == "PASS" and bench.compare("neg", -1.01, -1.0, 0.5, "s")["verdict"] == "FAIL"
    with pytest.raises(ValueError):                                      # a relative difference from zero means nothing
        bench.compare("zero", 1.0, 0.0, 1.0, "s")


def test_a_check_that_is_not_a_number_has_the_same_shape():
    import bench
    yes = bench.check("every group is there", True, "35 + 5 groups", "35 + 5 groups", "motorgeo.build")
    no = bench.check("every group is there", False, "35 + 5 groups", "34 + 5 groups, missing coil_0_0", "motorgeo.build")
    assert yes["verdict"] == "PASS" and no["verdict"] == "FAIL"
    for item in (yes, no):
        assert set(item) >= {"name", "expected", "measured", "source", "verdict"}


def test_every_check_says_whether_it_gates_the_verdict_and_one_that_does_not_says_why():
    """I2 of the review of Task 11: B3's laminated run carries the verdict and the plain run is reported beside it. A check that does not gate the verdict
    is still made, still shown with its own pass or fail, and still in results.json: nothing is hidden, it only does not decide the exit status."""
    import bench
    assert bench.check("x", True, "e", "m", "s")["gating"] is True and bench.compare("n", 1.0, 1.0, 1.0, "s")["gating"] is True      # gating unless said otherwise
    why = "the verdict is the laminated run's"
    words = "reported, not gating (the verdict is the laminated run's)"
    no = bench.check("x", False, "the same as FEMM's", "17", "s", gating=False, reason=why)
    assert no["gating"] is False and no["verdict"] == "FAIL" and no["expected"] == f"{words}: the same as FEMM's"      # its own verdict is kept, and its expected text says so
    number = bench.compare("n", 18.0, 16.61, 3.0, "FEMM", "mWb", gating=False, reason=why)
    assert number["gating"] is False and number["verdict"] == "FAIL" and number["expected"] == 16.61 and number["note"] == words    # the number stays a number
    assert bench.describe(number)[0] == f"{words}: 16.61 mWb within 3 %  [FEMM]"
    assert bench.describe(no)[0] == f"{words}: the same as FEMM's  [s]"                                                          # not said twice
    assert bench.describe(bench.compare("n", 18.0, 16.61, 3.0, "FEMM", "mWb"))[0] == "16.61 mWb within 3 %  [FEMM]"          # a gating check reads as before
    with pytest.raises(ValueError):                                                                                            # a check that does not gate has to say why
        bench.check("x", True, "e", "m", "s", gating=False)
    with pytest.raises(ValueError):
        bench.compare("n", 1.0, 1.0, 1.0, "s", gating=False)


def test_a_check_that_cannot_be_made_here_is_skipped_with_its_reason_and_never_gates():
    import bench
    item = bench.skipped("steel: psi_PM against the study's own", "the study's number, within 0.01 %", "out\\barrier\\summary.csv is not there: barrier_study.py makes it", "barrier_study.py")
    assert item["verdict"] == "SKIPPED" and item["gating"] is False and item["measured"].startswith("out\\barrier\\summary.csv is not there")
    assert set(item) >= {"name", "expected", "measured", "source", "verdict", "gating"}
    assert bench.describe(item) == ("the study's number, within 0.01 %  [barrier_study.py]", item["measured"])
    assert bench.is_gating(item) is False and bench.is_gating(bench.check("x", True, "e", "m", "s")) is True and bench.is_gating({"verdict": "FAIL"}) is True


def test_the_parameters_of_the_motor_go_into_results_json_through_one_public_helper():
    """M8 of the review of Task 11: _params was written twice (B1, B2) and B3 imported B2's private copy. There is one, in bench/__init__.py."""
    import bench
    import motorgeo as mg
    got = bench.params(mg.Params(barrier="none"))
    assert got["barrier"] == "none" and got["slots"] == 36 and got["pole_pairs"] == 4 and got["slot_fillet"] == 0.3
    assert set(got) == set(mg.Params.__dataclass_fields__)
    json.dumps(got)                                                              # and it is what results.json can hold


CASE_MODULES = ("b1_femm", "b2_gmsh", "b3_elmer", "b4_ccx", "b5_prepomax", "b6_ross", "b7_salib", "b8_matlab", "b9_rules")


def case_module_names(tree, cases=CASE_MODULES) -> set:
    """The names by which a file knows the case modules: their own, and every alias it gives them (`from . import b4_ccx as b4`, `from bench import b4_ccx as b4`,
    `import bench.b4_ccx as b4`)."""
    known = set(cases)
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "") in ("", "bench"):
            known |= {a.asname or a.name for a in node.names if a.name in cases}
        elif isinstance(node, ast.Import):
            known |= {a.asname for a in node.names if a.asname and a.name.split(".")[-1] in cases}
    return known


def private_name_offenders(file_name: str, tree, cases=CASE_MODULES) -> list[str]:
    """What a file takes from another case's private names: `from .b2_gmsh import _params`, and `b2_gmsh._params` (or `b2._params` through an alias) used as an attribute."""
    known = case_module_names(tree, cases)
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").split(".")[-1] in tuple(cases) + ("run_bench",):
            offenders += [f"{file_name}:{node.lineno} imports {a.name} from {node.module}" for a in node.names if a.name.startswith("_")]
        if isinstance(node, ast.Attribute) and node.attr.startswith("_") and isinstance(node.value, ast.Name) and node.value.id in known:
            offenders.append(f"{file_name}:{node.lineno} uses {node.value.id}.{node.attr}")
    return offenders


def test_no_case_defines_its_own_copy_of_a_shared_helper_or_reaches_into_another_case_s_private_names():
    """The cases share what they share through bench/__init__.py (params, fig_field, legend_and_model_figure, check, compare ...). A name with a leading underscore belongs to its
    module: another module that imports it (`from .b2_gmsh import _params`) or reaches for it (`b2_gmsh._params(p)`, or through an alias: `from . import b4_ccx as b4` and `b4._x`)
    is a copy that was not made."""
    offenders = []
    for file in sorted((P / "bench").glob("*.py")):
        tree = ast.parse(file.read_text(encoding="utf-8"))
        offenders += private_name_offenders(file.name, tree)
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name in ("_params", "params", "fig_field", "legend_and_model_figure") and file.name != "__init__.py":
                offenders.append(f"{file.name}:{node.lineno} defines {node.name}, which bench/__init__.py has")
    assert not offenders, "\n".join(offenders)
    init = ast.parse((P / "bench" / "__init__.py").read_text(encoding="utf-8"))
    assert {n.name for n in init.body if isinstance(n, ast.FunctionDef)} >= {"params", "fig_field", "legend_and_model_figure"}


def test_the_scan_for_private_names_sees_a_case_module_through_an_alias():
    """Fix round 1, item 11: B5 imports B4 as `from . import b4_ccx as b4`, and a scan that keyed on the full module names never saw `b4._x`."""
    def offenders(source):
        return private_name_offenders("b5_prepomax.py", ast.parse(source))
    assert offenders("from . import b4_ccx as b4\nvalue = b4._hidden\n") == ["b5_prepomax.py:2 uses b4._hidden"]
    assert offenders("from bench import b4_ccx as other\nother._hidden()\n") == ["b5_prepomax.py:2 uses other._hidden"]
    assert offenders("import bench.b4_ccx as b4\nb4._hidden\n") == ["b5_prepomax.py:2 uses b4._hidden"]
    assert offenders("from . import b4_ccx\nb4_ccx._hidden\n") == ["b5_prepomax.py:2 uses b4_ccx._hidden"]               # without an alias: seen before the fix too
    assert offenders("from .b4_ccx import _hidden\n") == ["b5_prepomax.py:1 imports _hidden from b4_ccx"]
    assert offenders("from . import b4_ccx as b4\nb4.public()\nb4.CONSTANT\n") == []                                       # public names are what the cases share
    assert offenders("from . import something_else as b4\nb4._hidden\n") == []                                              # an alias of something that is not a case module


def test_the_field_picture_wraps_a_long_title_so_that_the_picture_has_no_blank_margin(tmp_path):
    """M11 of the review of Task 11: B1's one-line title was wider than the plot, and the picture saved with a tight box was 45 % blank. The one function that
    draws the field wraps every line of the title to the width of the figure, so no case can bring that back."""
    import matplotlib.image as image
    import bench
    xy = np.array([[10.0, 0.0], [20.0, 0.0], [10.0, 10.0], [20.0, 10.0]])
    triangles = np.array([[0, 1, 2], [1, 3, 2]])
    potential = np.array([0.0, 1e-5, 2e-5, 3e-5])
    flux_density = np.array([0.5, 1.0])
    title = ("FEMM 5.0 (femm_pool): the band-free rotor, barrier=\"none\", no load at the drawn position\nflux density |B| and flux lines (contours of A); "
             "fine mesh, 15,859 nodes in the 90-degree sector; psi_PM 16.63 mWb, MTPA 1.1291 N.m, cogging 50.3 mN.m peak to peak")
    bench.fig_field(xy, triangles, potential, flux_density, title, tmp_path / "field.png")
    height, width = image.imread(str(tmp_path / "field.png")).shape[:2]
    assert width <= 1.1 * 9 * 150, f"{width} px wide: the title is wider than the 9-inch figure"
    assert height > 0.4 * width                                                      # and the cross-section, not a strip of title, is most of it


def test_a_tool_s_legend_and_model_are_put_together_under_a_caption_that_does_not_widen_the_picture(tmp_path):
    """Fix round 1, item 12: the function that sets a tool's own picture (cgx's frames in B4, PrePoMax's viewport in B5) together, a legend at the left and the model at the right under
    a caption, lives next to fig_field in bench/__init__.py. Both images keep their own pixels, the caption's lines are wrapped to TITLE_CHARS, and the picture is not mostly blank."""
    import matplotlib.image as image
    import bench
    legend = np.zeros((300, 60, 3), np.uint8)
    legend[:, :, 0] = 200
    model = np.zeros((300, 540, 3), np.uint8)
    model[:, :, 2] = 180
    caption = "PrePoMax: a caption of two lines\n" + "and a long second line " * 8
    bench.legend_and_model_figure(legend, model, caption, tmp_path / "both.png")
    picture = image.imread(str(tmp_path / "both.png"))
    height, width = picture.shape[:2]
    assert width <= 1.1 * 11 * 130, f"{width} px wide: the caption widened the 11-inch figure"
    assert height > 0.45 * width                                                     # the two images, not a strip of caption, are most of it
    red = (picture[..., 0] > 0.7) & (picture[..., 1] < 0.2) & (picture[..., 2] < 0.2)        # the legend's pixels (not the white around them)
    blue = (picture[..., 2] > 0.6) & (picture[..., 0] < 0.2) & (picture[..., 1] < 0.2)       # the model's
    assert red.any() and blue.any() and 5 * red.sum() < blue.sum() < 12 * red.sum()          # 60 and 540 pixels wide: each keeps its own width (a ratio of 9)
    columns = np.flatnonzero(red.any(axis=0))
    assert columns.max() < np.flatnonzero(blue.any(axis=0)).min()                            # the legend is at the left of the model


def test_a_picture_title_is_wrapped_line_by_line_and_a_short_one_is_left_alone():
    """bench.wrap_title is what fig_field and B2's bands.png wrap their titles with: each line of the title to TITLE_CHARS characters, the breaks the case wrote kept."""
    import bench
    short = "FEMM 5.0: a short first line\nand a second"
    assert bench.wrap_title(short) == short
    wrapped = bench.wrap_title("one line that is far too long for a figure of nine inches " * 5 + "\nand a short one")
    lines = wrapped.split("\n")
    assert len(lines) >= 4 and all(len(line) <= bench.TITLE_CHARS for line in lines) and lines[-1] == "and a short one"


def test_the_results_file_holds_numpy_numbers_paths_and_non_ascii_text(tmp_path):
    import bench
    data = {"a": np.float64(1.5), "b": np.arange(3), "c": np.int64(7), "d": Path("x") / "y", "e": "ψ and µ", "f": (1, 2), "g": np.bool_(True),
            "h": float("nan")}
    target = bench.write_json(tmp_path / "results.json", data)
    back = json.loads(target.read_text(encoding="utf-8"))
    assert back["a"] == 1.5 and back["b"] == [0, 1, 2] and back["c"] == 7 and back["e"] == "ψ and µ" and back["f"] == [1, 2] and back["g"] is True
    assert back["d"] in ("x\\y", "x/y")
    assert back["h"] is None                                             # NaN is not JSON: it becomes null, never the word NaN


def test_the_priority_class_of_this_process_is_read_back_and_named():
    import bench
    value = bench.priority_class(os.getpid())
    assert isinstance(value, int) and bench.class_name(value) in ("normal", "below_normal", "idle", "above_normal", "high", "realtime")
    assert bench.class_name(0x4000) == "below_normal" and bench.class_name(0x20) == "normal" and bench.class_name(0x40) == "idle"
    assert bench.class_name(0x1234) == "0x1234"                          # an unknown value is shown, not hidden
    assert bench.priority_class(2 ** 31 - 2) is None                     # no such process: None, not an exception


def test_the_watcher_sees_the_priority_class_of_a_child_this_process_starts():
    """The log has to name the priority class of the program a run_solver call started. The watcher reads it from outside while the
    program runs: here a stand-in child born below normal."""
    import bench
    child = None
    with bench.watch_children() as watcher:
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(1.5)"], creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        try:
            child.wait(timeout=30)
        finally:
            if child.poll() is None:
                child.kill()
    seen = watcher.seen()
    names = {name.lower() for name in seen}
    assert "python.exe" in names, seen
    assert {c for name, classes in seen.items() if name.lower() == "python.exe" for c in classes} == {0x4000}
    assert watcher.summary()["python.exe"] == ["below_normal"]


def test_the_log_keeps_what_was_said_and_the_runs_and_writes_both(tmp_path, capsys):
    import bench
    log = bench.Log(tmp_path)
    log.say("first line")
    log.add_run({"kind": "femm_pool", "label": "bench B1", "workers": 8, "priority_classes": {"below_normal": 8}})
    with log.tee():
        print("[femm_pool bench B1] starting 8 worker(s)")                 # what a pool prints is kept too
    log.save()
    text = (tmp_path / "run.log").read_text(encoding="utf-8")
    assert "first line" in text and "[femm_pool bench B1] starting 8 worker(s)" in text
    runs = json.loads((tmp_path / "runs.json").read_text(encoding="utf-8"))
    assert runs == [{"kind": "femm_pool", "label": "bench B1", "workers": 8, "priority_classes": {"below_normal": 8}}]
    shown = capsys.readouterr().out
    assert "first line" in shown and "starting 8 worker(s)" in shown     # and it was shown as it happened


# ----------------------------------------------------------------------------------------------- 2. R10 for the bench: no program started outside the queue, no process ended, every call has its time limit
# B9 shows the rules with TWO bench runs at once, so it has to start two run_bench.py processes (Python, never a solver) as its children. Starting them through run_solver
# would hold the solver slot that the children need, and they would wait for it, and for each other, until the timeout. That is the ONE exception to "nothing in bench/ starts
# a program", named here with its reason (R10 and B9's rule: the two runs must be two real processes that meet only at the solver slot). The only process B9 may end is its own
# child, on its own timeout (end_child, and only it).
# Anything else that starts a program, or ends a process, in any module of bench/ (B9's included) is still an offender.
NAMED_EXCEPTIONS = {
    "b9_rules.py": {
        "starts": {("start_child", "subprocess.Popen")},
        "ends": {("end_child", "kill"), ("end_child", "TerminateProcess")},
        "reason": "B9 starts its two run_bench.py children itself (Python, never a solver): run_solver would hold the slot they need and deadlock them until the timeout; "
                  "it ends only those children (the child by its own handle, and the interpreter below it through one handle it opens and checks against the tracker's record), "
                  "on its own timeout.",
    },
}


def scoped_enders(source: str) -> set[tuple[str, str]]:
    """(the function it is in, the name) for everything in the source that could end a process: kill, killpg, terminate, TerminateProcess, TerminateJobObject, as an attribute
    or a bare name. The scope is what lets the named exception be one function, not a whole file."""
    from test_check_tools import ENDERS, ScopeWalker
    found: set[tuple[str, str]] = set()

    class Walker(ScopeWalker):
        def visit_Attribute(self, node):
            if node.attr in ENDERS:
                found.add((self.here(), node.attr))
            self.generic_visit(node)

        def visit_Name(self, node):
            if node.id in ENDERS:
                found.add((self.here(), node.id))

    Walker().visit(ast.parse(source))
    return found


def program_and_ender_offenders(file_name: str, source: str) -> list[str]:
    """What in one module of bench/ starts a program or ends a process beyond what NAMED_EXCEPTIONS allows for that file."""
    from test_check_tools import process_enders, programs_started
    allowed = NAMED_EXCEPTIONS.get(file_name, {})
    offenders = []
    started = programs_started(source)
    if started != allowed.get("starts", set()):
        offenders.append(f"{file_name} starts {sorted(started)}, may start {sorted(allowed.get('starts', set()))}")
    ended = scoped_enders(source)
    if ended != allowed.get("ends", set()):
        offenders.append(f"{file_name} ends a process in {sorted(ended)}, may end in {sorted(allowed.get('ends', set()))}")
    if ("string", "taskkill") in process_enders(source):
        offenders.append(f"{file_name} runs taskkill (by name or number, it is no exception)")
    return offenders


def test_nothing_in_bench_starts_a_program_except_through_femm_queue_and_nothing_ends_a_process():
    """R10 for the bench, with its one named exception: B9's own two children (see NAMED_EXCEPTIONS)."""
    files = sorted((P / "bench").glob("*.py"))
    assert files, "bench/ has no module"
    offenders = []
    for file in files:
        offenders += program_and_ender_offenders(file.name, file.read_text(encoding="utf-8"))
    assert not offenders, "\n".join(offenders)


def test_the_named_exception_is_b9s_start_child_and_end_child_and_nothing_wider():
    """The exception is one function that starts and one that ends: the scan on B9's own source finds exactly those, and on the same source with one more place that starts a
    program, or one more that ends a process, or taskkill, it fails (the scan is tried on small sources first, so that 'found nothing' means something)."""
    assert set(NAMED_EXCEPTIONS) == {"b9_rules.py"} and NAMED_EXCEPTIONS["b9_rules.py"]["reason"]
    assert NAMED_EXCEPTIONS["b9_rules.py"]["starts"] == {("start_child", "subprocess.Popen")}
    assert NAMED_EXCEPTIONS["b9_rules.py"]["ends"] == {("end_child", "kill"), ("end_child", "TerminateProcess")}
    b9 = (P / "bench" / "b9_rules.py").read_text(encoding="utf-8")
    assert program_and_ender_offenders("b9_rules.py", b9) == []
    base = "import subprocess\ndef start_child(c):\n    return subprocess.Popen(c)\ndef end_child(p, k, h):\n    p.kill()\n    k.TerminateProcess(h, 1)\n"
    assert program_and_ender_offenders("b9_rules.py", base) == []
    assert program_and_ender_offenders("b9_rules.py", base + "def sneaky(c):\n    subprocess.run(c)\n")                               # a second place that starts
    assert program_and_ender_offenders("b9_rules.py", base + "def sneaky(p):\n    p.terminate()\n")                                    # a second way to end
    assert program_and_ender_offenders("b9_rules.py", base + "def sneaky(p):\n    p.kill()\n")                                         # the same ender in another function
    assert program_and_ender_offenders("b9_rules.py", base + "def sneaky(k, h):\n    k.TerminateProcess(h, 1)\n")                       # the handle ender in another function
    assert program_and_ender_offenders("b9_rules.py", base.replace("    k.TerminateProcess(h, 1)\n", ""))                                  # the exception is for BOTH enders, or the file is wrong
    assert program_and_ender_offenders("b9_rules.py", base + "def sneaky():\n    subprocess.call('taskkill /F /PID 1')\n")           # taskkill
    assert program_and_ender_offenders("b9_rules.py", "import subprocess\ndef start_child(c):\n    return subprocess.Popen(c)\n")      # the exception is for BOTH, or the file is wrong
    assert program_and_ender_offenders("b7_salib.py", base)                                                                           # no other file has the exception
    assert program_and_ender_offenders("b7_salib.py", "import numpy\nx = 1\n") == []
    assert scoped_enders("import os\ndef f(pid):\n    os.kill(pid, 9)\nclass C:\n    def m(self, p):\n        p.terminate()\n") == {("f", "kill"), ("C.m", "terminate")}
    assert scoped_enders("def f(k, h):\n    k.TerminateProcess(h, 1)\n") == {("f", "TerminateProcess")}


# ----------------------------------------------------------------------------------------------- 2b. no import that nothing uses (Task 13's modules; review M7)
TASK_13_MODULES = ("femm_points.py", "b7_salib.py", "b8_matlab.py", "b9_rules.py")


def unused_imports(source: str) -> list[str]:
    """The names a module imports and never uses (a name is used when it appears anywhere as a Name, which includes the base of an attribute like np.array, and annotations count,
    since the modules use `from __future__ import annotations` and a parsed annotation is still a Name). `from __future__ import ...` is not an import of a name."""
    tree = ast.parse(source)
    imported: dict[str, int] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imported[(alias.asname or alias.name).split(".")[0]] = node.lineno
        elif isinstance(node, ast.ImportFrom) and node.module != "__future__":
            for alias in node.names:
                imported[alias.asname or alias.name] = node.lineno
    used = {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
    return sorted(f"{name} (line {line})" for name, line in imported.items() if name not in used)


def test_the_scan_for_unused_imports_finds_them_and_the_modules_of_task_13_have_none():
    """B7 imported `paths` and B9 `numpy` and used neither (review M7)."""
    assert unused_imports("import os\nimport numpy as np\nfrom pathlib import Path\nfrom . import check, skipped\nprint(os.getcwd(), Path('.'), check)\n") == ["np (line 2)", "skipped (line 4)"]
    assert unused_imports("from __future__ import annotations\nimport scipy.io\nx = scipy.io.loadmat\ndef f(a: Missing) -> None:\n    pass\n") == []
    assert unused_imports("import json\nimport numpy as np\n\ndef f(a: np.ndarray):\n    return a\n") == ["json (line 1)"]
    for name in TASK_13_MODULES:
        assert unused_imports((P / "bench" / name).read_text(encoding="utf-8")) == [], name


QUEUE_LIMITS = {"femm_pool": "task_timeout_s", "run_solver": "timeout_s"}            # femm_queue's two doors, and the keyword that gives each its time limit


def queue_calls(source: str) -> tuple[int, list[str]]:
    """(how many calls of femm_queue's femm_pool and run_solver the source makes, those among them that do not pass their time limit).

    A call is seen however it is spelled: `femm_queue.run_solver(...)` or `fq.run_solver(...)` (an attribute), `run_solver(...)` after `from femm_queue import
    run_solver` (the bare name) and `go(...)` after `from femm_queue import run_solver as go` (an alias). A bare `run_solver` or `femm_pool` is taken for the
    queue's whether or not the import is seen, so a name that is imported some other way is not a way round the scan."""
    tree = ast.parse(source)
    named = {name: name for name in QUEUE_LIMITS}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "femm_queue":
            named.update({a.asname: a.name for a in node.names if a.asname and a.name in QUEUE_LIMITS})
    found, missing = 0, []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        target = f.attr if isinstance(f, ast.Attribute) and f.attr in QUEUE_LIMITS else named.get(f.id) if isinstance(f, ast.Name) else None
        if target:
            found += 1
            if QUEUE_LIMITS[target] not in {k.arg for k in node.keywords}:
                missing.append(f"{target}(...) at line {node.lineno} without {QUEUE_LIMITS[target]}")
    return found, missing


def test_the_scan_for_the_queue_s_time_limits_sees_every_way_a_call_can_be_spelled():
    """M10 of the review of Task 11: the scan saw only attribute calls, so `from femm_queue import run_solver` and then `run_solver(cmd)` went past it."""
    assert queue_calls("femm_queue.run_solver(cmd, label='x', timeout_s=5)\nfemm_queue.femm_pool('x', task_timeout_s=9)\n") == (2, [])
    assert queue_calls("femm_queue.run_solver(cmd, label='x')\n") == (1, ["run_solver(...) at line 1 without timeout_s"])
    assert queue_calls("import femm_queue as fq\nwith fq.femm_pool('x') as pool:\n    pass\n") == (1, ["femm_pool(...) at line 2 without task_timeout_s"])
    bare = "from femm_queue import run_solver, femm_pool\nrun_solver(cmd, label='x')\nwith femm_pool('x') as pool:\n    pass\n"
    assert queue_calls(bare) == (2, ["run_solver(...) at line 2 without timeout_s", "femm_pool(...) at line 3 without task_timeout_s"])
    assert queue_calls("from femm_queue import run_solver as go\ngo(cmd, timeout_s=3)\ngo(cmd)\n") == (2, ["run_solver(...) at line 3 without timeout_s"])
    assert queue_calls("from femm_queue import run_solver\nrun_solver(cmd, timeout_s=1)\n") == (1, [])
    assert queue_calls("subprocess_free = run_something(cmd)\nother.run(cmd)\n") == (0, [])                  # other names are not the queue's


def test_every_case_that_runs_a_tool_goes_through_the_queue_and_every_wait_has_a_limit():
    from test_check_tools import scopes_calling
    b1 = (P / "bench" / "b1_femm.py").read_text(encoding="utf-8")
    assert scopes_calling(b1, "femm_queue", "femm_pool"), "B1 runs FEMM: it has to open a femm_pool"
    seen, offenders = 0, []
    for file in sorted((P / "bench").glob("*.py")):
        found, missing = queue_calls(file.read_text(encoding="utf-8"))
        seen += found
        offenders += [f"{file.name}: {m}" for m in missing]
    assert not offenders, "\n".join(offenders)
    assert seen >= 3, "the scan should see B1's two pools and run_program's run_solver"


# ----------------------------------------------------------------------------------------------- 3. the PATH that Gmsh cuts: every call that starts Gmsh is inside intact_path
# gmsh.initialize() cuts the PATH of the running process (Task 4 found it: 1,864 characters became 313, and os.environ never shows it). The
# fix is `with checks.intact_path():` around every call that starts Gmsh. Each function that does so has a test that reads the process's own
# PATH from Windows (tests/test_bench_b2.py and test_bench_b3.py); this scan is the net under them: a Gmsh call added to any module of
# bench/ without the guard fails here, whoever forgot it and wherever it is.
def gmsh_calls(source: str) -> tuple[list[str], list[str]]:
    """The calls in `source` that start Gmsh (`gmsh.initialize()` itself, and `to_gmsh(...)`, which calls it), as (guarded, unguarded).

    Guarded means: inside the body of a `with intact_path():` or `with checks.intact_path():` statement. A function defined inside such a
    block is not inside it (it runs when it is called, not where it is written). Each item is "<call> at line N"."""
    guarded, unguarded = [], []

    def is_guard(node) -> bool:
        if not isinstance(node, (ast.With, ast.AsyncWith)):
            return False
        for item in node.items:
            call = item.context_expr
            if isinstance(call, ast.Call):
                name = call.func.attr if isinstance(call.func, ast.Attribute) else getattr(call.func, "id", "")
                if name == "intact_path":
                    return True
        return False

    def starts_gmsh(node) -> str:
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Attribute) and f.attr == "initialize" and isinstance(f.value, ast.Name) and f.value.id == "gmsh":
                return "gmsh.initialize()"
            if (isinstance(f, ast.Attribute) and f.attr == "to_gmsh") or (isinstance(f, ast.Name) and f.id == "to_gmsh"):
                return "to_gmsh()"
        return ""

    def visit(node, inside: bool) -> None:
        found = starts_gmsh(node)
        if found:
            (guarded if inside else unguarded).append(f"{found} at line {node.lineno}")
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
            inside = False
        for field, value in ast.iter_fields(node):
            for child in (value if isinstance(value, list) else [value]):
                if isinstance(child, ast.AST):
                    visit(child, inside or (field == "body" and is_guard(node)))

    visit(ast.parse(source), False)
    return guarded, unguarded


def test_the_scan_for_gmsh_calls_without_the_path_guard_finds_what_it_should():
    """The scan is tried on small sources first, so that a scan that finds nothing is known to be able to find something."""
    both = "import gmsh\nwith intact_path():\n    gmsh.initialize()\n    try:\n        n = mg.to_gmsh(p, r, 'x.msh')\n    finally:\n        gmsh.finalize()\n"
    assert gmsh_calls(both) == (["gmsh.initialize() at line 3", "to_gmsh() at line 5"], [])
    assert gmsh_calls("with checks.intact_path():\n    gmsh.initialize()\n") == (["gmsh.initialize() at line 2"], [])
    assert gmsh_calls("gmsh.initialize()\n") == ([], ["gmsh.initialize() at line 1"])
    assert gmsh_calls("def f():\n    return to_gmsh(a, b, c)\n") == ([], ["to_gmsh() at line 2"])                    # by the bare name
    assert gmsh_calls("with intact_path():\n    pass\nmg.to_gmsh(a, b, c)\n") == ([], ["to_gmsh() at line 3"])         # after the block, not in it
    assert gmsh_calls("with open('f') as f:\n    gmsh.initialize()\n") == ([], ["gmsh.initialize() at line 2"])       # another kind of block
    assert gmsh_calls("with intact_path():\n    def later():\n        gmsh.initialize()\n") == ([], ["gmsh.initialize() at line 3"])   # written there, run later
    assert gmsh_calls("gmsh.finalize()\nother.initialize()\ngmsh.option.setNumber('a', 1)\n") == ([], [])             # not calls that start Gmsh


def test_every_call_in_bench_that_starts_gmsh_is_inside_the_path_guard():
    """Self-check: take the `with intact_path():` line off one of them (b2_gmsh.py has three, b3_elmer.py one) and this fails, naming the line."""
    guarded_total, offenders = 0, []
    for file in sorted((P / "bench").glob("*.py")):
        guarded, unguarded = gmsh_calls(file.read_text(encoding="utf-8"))
        guarded_total += len(guarded)
        offenders += [f"{file.name}: {call}" for call in unguarded]
    assert not offenders, "Gmsh cuts the PATH of the process; wrap these in `with intact_path():`\n" + "\n".join(offenders)
    assert guarded_total >= 4, "the scan should see at least inspect_mesh, make_path1_mesh, mesh_faces and interface_facts"


# ----------------------------------------------------------------------------------------------- 4. run_bench.py: --only, the exit codes (gating checks only), a case that raises, the folders
def fake_case(checks=None, then=None, calls=None, name="B1"):
    """A case that makes no solver run: it returns the given checks (or raises `then`) and notes that it was called."""
    import bench

    def run(folder, log):
        if calls is not None:
            calls.append(name)
        log.say(f"{name} says hello")
        if then is not None:
            raise then
        (Path(folder) / "figure.png").write_bytes(b"png")
        return {"case": name, "checks": checks if checks is not None else [bench.check("a thing", True, "x", "x", "nowhere")]}
    return run


def test_only_takes_case_names_in_any_case_with_commas_or_repeats_and_keeps_the_order_of_the_registry():
    from bench import run_bench
    registry = {"B1": None, "B2": None, "B3": None}
    assert run_bench.chosen_cases(registry, None) == ["B1", "B2", "B3"]
    assert run_bench.chosen_cases(registry, ["b3,B1"]) == ["B1", "B3"]
    assert run_bench.chosen_cases(registry, ["B2", "b3"]) == ["B2", "B3"]
    assert run_bench.chosen_cases(registry, ["B7", "B1"]) is None and run_bench.chosen_cases(registry, [""]) is None


def test_a_wrong_name_is_exit_code_2_and_names_the_cases(capsys, monkeypatch, tmp_path):
    from bench import run_bench
    use_config(monkeypatch, tmp_path)
    assert run_bench.main(["--only", "B7"], registry={"B1": fake_case()}) == 2
    err = capsys.readouterr().err
    assert "B7" in err and "B1" in err


def test_every_check_passing_is_exit_code_0_and_each_case_leaves_its_files(monkeypatch, tmp_path, capsys):
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    calls = []
    assert run_bench.main([], registry={"B1": fake_case(calls=calls), "B2": fake_case(calls=calls, name="B2")}) == 0
    assert calls == ["B1", "B2"]
    for name in ("B1", "B2"):
        folder = out / "bench" / name
        assert (folder / "results.json").is_file() and (folder / "run.log").is_file() and (folder / "runs.json").is_file()
        assert f"{name} says hello" in (folder / "run.log").read_text(encoding="utf-8")
        assert json.loads((folder / "results.json").read_text(encoding="utf-8"))["case"] == name
    shown = capsys.readouterr().out
    assert "PASS" in shown and "bench: PASS" in shown
    assert json.loads((out / "bench" / "summary.json").read_text(encoding="utf-8"))["cases"] == ["B1", "B2"]


def test_a_failed_expectation_is_exit_code_1_and_is_reported_not_hidden(monkeypatch, tmp_path, capsys):
    import bench
    from bench import run_bench
    use_config(monkeypatch, tmp_path)
    bad = [bench.compare("psi_pm", 17.0, 16.61, 0.5, "README row", "mWb")]
    assert run_bench.main([], registry={"B1": fake_case(bad)}) == 1
    shown = capsys.readouterr().out
    assert "FAIL" in shown and "psi_pm" in shown and "README row" in shown and "bench: FAIL" in shown


def reported_and_gating(passing=True):
    """One gating check and, beside it, one that is reported and does not gate (as B3's plain run is beside its laminated run)."""
    import bench
    return [bench.check("laminated: the pattern is alike", passing, "alike", "alike" if passing else "not alike", "FEMM"),
            bench.check("plain: the pattern is alike", False, "alike", "not alike", "FEMM", gating=False, reason="the verdict is the laminated run's")]


def test_a_reported_check_that_fails_leaves_the_exit_status_at_0_and_is_still_shown_and_counted(monkeypatch, tmp_path, capsys):
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    assert run_bench.main([], registry={"B3": fake_case(reported_and_gating(), name="B3")}) == 0
    shown = capsys.readouterr().out
    assert "  PASS  B3: laminated: the pattern is alike" in shown
    assert "  REPORTED (FAIL)  B3: plain: the pattern is alike" in shown                                       # its own verdict, shown
    assert "reported, not gating (the verdict is the laminated run's): alike" in shown
    assert "bench: PASS, all 1 gating checks pass" in shown and "1 reported, not gating (0 pass, 1 fail, 0 skipped)" in shown
    summary = json.loads((out / "bench" / "summary.json").read_text(encoding="utf-8"))
    assert summary["checks"] == 2 and summary["failed"] == {}
    assert summary["gating"] == {"checks": 1, "passed": 1, "failed": 0}
    assert summary["reported"] == {"checks": 1, "passed": 0, "failed": 1, "skipped": 0} and summary["reported_failed"] == {"B3": ["plain: the pattern is alike"]}
    results = json.loads((out / "bench" / "B3" / "results.json").read_text(encoding="utf-8"))
    assert [c["gating"] for c in results["checks"]] == [True, False] and results["checks"][1]["verdict"] == "FAIL"      # nothing is hidden in the file either


def test_a_gating_check_that_fails_is_exit_code_1_whatever_the_reported_ones_say(monkeypatch, tmp_path, capsys):
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    assert run_bench.main([], registry={"B3": fake_case(reported_and_gating(passing=False), name="B3")}) == 1
    shown = capsys.readouterr().out
    assert "  FAIL  B3: laminated: the pattern is alike" in shown and "bench: FAIL, 1 of 1 gating checks fail: B3 (1)" in shown
    summary = json.loads((out / "bench" / "summary.json").read_text(encoding="utf-8"))
    assert summary["failed"] == {"B3": ["laminated: the pattern is alike"]} and summary["gating"] == {"checks": 1, "passed": 0, "failed": 1}


def test_a_skipped_check_is_shown_counted_apart_and_never_fails_the_run(monkeypatch, tmp_path, capsys):
    import bench
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    checks = [bench.check("a thing", True, "x", "x", "nowhere"), bench.skipped("steel: psi_PM, the study's own", "the study's number", "out\\barrier\\summary.csv is not there", "barrier_study.py")]
    assert run_bench.main([], registry={"B1": fake_case(checks)}) == 0
    shown = capsys.readouterr().out
    assert "  SKIPPED  B1: steel: psi_PM, the study's own" in shown and "measured: out\\barrier\\summary.csv is not there" in shown
    summary = json.loads((out / "bench" / "summary.json").read_text(encoding="utf-8"))
    assert summary["reported"] == {"checks": 1, "passed": 0, "failed": 0, "skipped": 1} and summary["reported_failed"] == {}
    assert summary["gating"] == {"checks": 1, "passed": 1, "failed": 0}


def test_a_case_that_raises_is_a_failed_case_with_its_traceback_in_a_file_and_the_next_case_still_runs(monkeypatch, tmp_path, capsys):
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    calls = []
    registry = {"B1": fake_case(then=RuntimeError("the pool lost a worker"), calls=calls), "B2": fake_case(calls=calls, name="B2")}
    assert run_bench.main([], registry=registry) == 1
    assert calls == ["B1", "B2"]
    assert "the pool lost a worker" in (out / "bench" / "B1" / "traceback.txt").read_text(encoding="utf-8")
    shown = capsys.readouterr().out
    assert "error: RuntimeError: the pool lost a worker" in shown and "Traceback" not in shown


def test_a_case_folder_is_emptied_before_its_run_and_only_that_case_s(monkeypatch, tmp_path):
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    (out / "bench" / "B1").mkdir(parents=True)
    (out / "bench" / "B1" / "old.txt").write_text("from the last run")
    (out / "bench" / "B3").mkdir()
    (out / "bench" / "B3" / "other.txt").write_text("another case's")
    assert run_bench.main(["--only", "B1"], registry={"B1": fake_case(), "B3": fake_case(name="B3")}) == 0
    assert not (out / "bench" / "B1" / "old.txt").exists()
    assert (out / "bench" / "B3" / "other.txt").read_text() == "another case's"


def test_the_real_registry_has_the_cases_of_tasks_11_to_13_and_runs_them_in_order():
    from bench import run_bench
    assert list(run_bench.registry()) == ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9"]
    assert run_bench.CASES["B7"] == "b7_salib" and run_bench.CASES["B8"] == "b8_matlab" and run_bench.CASES["B9"] == "b9_rules"


# ----------------------------------------------------------------------------------------------- 5. the bench's own process: lowered, and its record (the threads read back from Gmsh)
def test_the_fixture_that_keeps_run_bench_from_lowering_the_pytest_run_is_in_force_in_this_file():
    """run_bench.main() lowers its own process (R10). The fixture that stops a test from lowering the whole pytest run lives in bench_test_tools.py and is imported
    at the top of this file. This asserts WHAT stands in the place of femm_session.lower_priority: the fixture's own function. It does not read the priority class
    from the operating system, because in a whole-file run the tests before this one (they call run_bench.main()) would already have lowered the process, and
    the test could not go red. Take the import away and this fails whatever ran first (P1 of the re-review of Task 11)."""
    import bench_test_tools
    import femm_session
    assert femm_session.lower_priority is bench_test_tools.lower_priority_is_off


def test_the_bench_process_is_lowered_and_each_case_s_runs_json_ends_with_the_class_the_operating_system_reports(monkeypatch, tmp_path):
    """R10 says each run is below normal. The bench's own Python (Gmsh's library, numpy, matplotlib) is a run too: it is lowered when the bench starts,
    and the record that B9 reads names the class read back from the operating system, not the class that was asked for."""
    import femm_session
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    lowered = []
    monkeypatch.setattr(femm_session, "lower_priority", lambda: lowered.append(True))
    assert run_bench.main([], registry={"B1": fake_case()}) == 0
    assert lowered == [True]
    runs = json.loads((out / "bench" / "B1" / "runs.json").read_text(encoding="utf-8"))
    record = runs[-1]
    assert record["kind"] == "in_process" and record["pid"] == os.getpid() and record["threads"] >= 1
    import bench
    assert record["priority_class_seen"] == bench.class_name(bench.priority_class(os.getpid()))
    assert "Gmsh" in record["label"]


def test_the_process_record_names_the_class_it_reads_back():
    import bench
    record = bench.process_record()
    assert record["kind"] == "in_process" and record["pid"] == os.getpid() and record["threads"] >= 1
    assert record["priority_class_seen"] in ("normal", "below_normal", "idle", "above_normal", "high", "realtime")


class FakeGmsh:
    """A stand-in for the gmsh module that answers its options with the numbers it was given, and notes the calls: so that a test can tell a number that is READ
    from Gmsh from one that is written into the code."""

    def __init__(self, options):
        import types
        self.calls, self.on = [], False
        self.option = types.SimpleNamespace(getNumber=lambda name: options[name], setNumber=lambda name, value: self.calls.append(("set", name, value)))

    def initialize(self, *args, **kwargs):
        self.calls.append("initialize")
        self.on = True

    def finalize(self):
        self.calls.append("finalize")
        self.on = False

    def isInitialized(self):
        return int(self.on)


def test_the_threads_of_the_process_record_are_read_back_from_gmsh_and_not_written_into_the_code(monkeypatch):
    """M9 of the review of Task 11: `"threads": 1` was a constant, so B9's 'at most 8 threads' could not fail for the work done in this process. The record
    now asks Gmsh (General.NumThreads, and Mesh.MaxNumThreads2D, which takes over for the 2-D mesher when it is not 0) and says what it was told."""
    import bench
    options = {"General.NumThreads": 6.0, "Mesh.MaxNumThreads1D": 0.0, "Mesh.MaxNumThreads2D": 0.0, "Mesh.MaxNumThreads3D": 0.0}
    fake = FakeGmsh(options)
    monkeypatch.setitem(sys.modules, "gmsh", fake)
    record = bench.process_record()
    assert record["threads"] == 6 and record["gmsh_options"]["General.NumThreads"] == 6 and record["gmsh_options"]["Mesh.MaxNumThreads2D"] == 0
    assert fake.calls[0] == "initialize" and fake.calls[-1] == "finalize" and not fake.on            # a session opened for the question, and closed
    options["Mesh.MaxNumThreads2D"] = 2.0                                                            # the 2-D mesher's own limit wins when it is set
    assert bench.process_record()["threads"] == 2
    options["General.NumThreads"] = 12.0                                                             # and a number above the 8 of config.toml is shown, not hidden
    options["Mesh.MaxNumThreads2D"] = 0.0
    record = bench.process_record()
    assert record["threads"] == 12 and record["threads_limit"] == paths.setting("femm", "workers")


def test_a_gmsh_that_is_already_running_is_asked_and_left_running(monkeypatch):
    import bench
    fake = FakeGmsh({"General.NumThreads": 3.0, "Mesh.MaxNumThreads1D": 0.0, "Mesh.MaxNumThreads2D": 0.0, "Mesh.MaxNumThreads3D": 0.0})
    fake.on = True
    monkeypatch.setitem(sys.modules, "gmsh", fake)
    assert bench.process_record()["threads"] == 3 and fake.on and "initialize" not in fake.calls and "finalize" not in fake.calls


def test_the_real_gmsh_is_asked_and_the_answer_is_in_the_record(path_is_put_back):
    import gmsh

    import bench
    import checks
    with checks.intact_path():
        gmsh.initialize()
        try:
            general, two = int(gmsh.option.getNumber("General.NumThreads")), int(gmsh.option.getNumber("Mesh.MaxNumThreads2D"))
        finally:
            gmsh.finalize()
    record = bench.process_record()
    assert record["gmsh_options"]["General.NumThreads"] == general and record["gmsh_options"]["Mesh.MaxNumThreads2D"] == two
    assert record["threads"] == (two or general) and record["threads"] >= 1
    from conftest import process_path
    assert process_path() == path_is_put_back                                                       # and the PATH of the process is whole afterwards


# ----------------------------------------------------------------------------------------------- 6. the roots and the lock (Task 13, B9's rulings)
# B9 starts two bench runs at once. Two runs that write into the same out\bench\<case>\ would empty each other's folders (fresh()), so (1) a run can be told another root
# (`run_bench.py --out ROOT`; every case writes under ROOT\<case>\, and fresh() refuses anything outside it), and (2) a run holds a lock on its root for the whole run,
# so that two runs on the SAME root take turns. The lock is femm_queue's own (an operating-system lock on a file), taken with solver_slot(lock_dir=...): reused, not copied.
def test_the_root_is_out_bench_until_a_run_is_told_another_and_is_put_back(monkeypatch, tmp_path):
    import bench
    out = use_config(monkeypatch, tmp_path)
    assert bench.root() == out / "bench" == bench.default_root()
    other = out / "bench" / "B9" / "run1"
    with bench.use_root(other) as active:
        assert bench.root() == active == other.resolve()
        assert bench.case_folder("B1") == other.resolve() / "B1" and (other / "B1").is_dir()
    assert bench.root() == out / "bench"
    with pytest.raises(RuntimeError):
        with bench.use_root(other):
            raise RuntimeError("a case failed")
    assert bench.root() == out / "bench"                                  # put back after an error too


def test_a_root_has_to_be_inside_out_and_is_not_out_itself(monkeypatch, tmp_path):
    """fresh() empties ROOT\\B<n>. A root anywhere would let a typo aim it at somebody's folder named B1: it has to be inside out\\ (config.toml [paths] out), where the bench's own
    files are, and not out\\ itself. A relative one starts at the project folder, as in paths.py."""
    import bench
    out = use_config(monkeypatch, tmp_path)
    monkeypatch.setattr(paths, "P", tmp_path)
    for good in (out / "bench", out / "bench" / "B9" / "run1", out / "elsewhere", "out/bench2"):
        assert bench.root_allowed(good).is_absolute() and out.resolve() in bench.root_allowed(good).parents
    for bad in (out, tmp_path, tmp_path / "reference", out / ".." / "reference", out.parent / "out2", "reference", "..", "C:/"):
        with pytest.raises(ValueError, match="out"):
            bench.root_allowed(bad)


def make_junction(link: Path, target: Path):
    """A directory junction (mklink /J: it needs no administrator rights) at `link` pointing at `target`; the tests' own, so skipped with the reason where it cannot be made."""
    made = subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(target)], capture_output=True, text=True)
    if made.returncode != 0 or not link.exists():
        pytest.skip(f"a junction could not be made here: {made.stdout.strip()} {made.stderr.strip()}")


def test_a_junction_inside_out_is_judged_by_where_it_leads_and_fresh_cannot_be_led_out_by_one(monkeypatch, tmp_path):
    """root_allowed and fresh resolve the path first, so a junction is judged by where it leads (review M9: no case tried one). A junction in out\\bench that leads to somebody's
    folder is refused as a root; one that leads to another place inside out\\ is allowed, as that place; and a case folder B1 that is a junction to somebody's B1 is not emptied."""
    import bench
    out = use_config(monkeypatch, tmp_path)
    monkeypatch.setattr(paths, "P", tmp_path)
    (out / "bench" / "real").mkdir(parents=True)
    somebody = tmp_path / "somebody"
    (somebody / "B1").mkdir(parents=True)
    (somebody / "B1" / "keep.txt").write_text("somebody's file")
    made = []
    try:
        away, alias, case = out / "bench" / "away", out / "bench" / "alias", out / "bench" / "real" / "B1"
        for link, target in ((away, somebody), (alias, out / "bench" / "real"), (case, somebody / "B1")):
            make_junction(link, target)
            made.append(link)
        with pytest.raises(ValueError, match="out"):
            bench.root_allowed(away)                                                                       # in out\ by its name, in somebody's folder by where it leads
        with pytest.raises(ValueError, match="out"):
            bench.root_allowed(away / "B1")
        assert bench.root_allowed(alias) == (out / "bench" / "real").resolve()                             # inside out\ by where it leads: allowed, as that place
        with bench.use_root(out / "bench" / "real"):
            with pytest.raises(ValueError):
                bench.fresh(case)                                                                           # B1 of the root is a junction to somebody's B1: it resolves outside the root
        with pytest.raises(ValueError):
            with bench.use_root(away):
                pass
    finally:
        for link in made:
            os.rmdir(link)                                                                                  # a junction is removed as itself, whatever it points at
    assert (somebody / "B1" / "keep.txt").read_text() == "somebody's file"


def test_fresh_refuses_anything_outside_the_active_root(monkeypatch, tmp_path):
    import bench
    out = use_config(monkeypatch, tmp_path)
    default_b1 = bench.case_folder("B1")                                  # out\bench\B1: the default root's
    (default_b1 / "keep.txt").write_text("the other root's")
    with bench.use_root(out / "bench" / "B9" / "run1") as active:
        mine = bench.case_folder("B1")
        assert mine == active / "B1" and mine != default_b1
        (mine / "old.txt").write_text("old")
        assert bench.fresh(mine) == mine and list(mine.iterdir()) == []
        for bad in (default_b1, active, active / "B1" / "fem", out / "bench" / "B9" / "run1" / ".." / ".." / "B1", out / "bench" / "B9"):
            with pytest.raises(ValueError):
                bench.fresh(bad)
    assert (default_b1 / "keep.txt").read_text() == "the other root's"


def test_main_writes_under_the_root_it_is_given_and_says_where_and_who_it_is(monkeypatch, tmp_path, capsys):
    import bench
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    root = out / "bench" / "B9" / "run1"
    assert run_bench.main(["--out", str(root), "--only", "B1"], registry={"B1": fake_case(), "B2": fake_case(name="B2")}) == 0
    assert (root / "B1" / "results.json").is_file() and (root / "B1" / "run.log").is_file() and not (root / "B2").exists()
    assert json.loads((root / "summary.json").read_text(encoding="utf-8"))["cases"] == ["B1"]
    assert not (out / "bench" / "B1").exists() and not (out / "bench" / "summary.json").exists()       # the default root is untouched
    assert bench.root() == out / "bench"                                                                # and is the root again afterwards
    shown = capsys.readouterr().out
    assert f"bench: root {root.resolve()} (pid {os.getpid()})" in shown                                 # B9 reads the pid of its children from this line


def test_a_root_outside_out_is_refused_with_exit_code_2_and_nothing_is_made(monkeypatch, tmp_path, capsys):
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    elsewhere = tmp_path / "elsewhere"
    assert run_bench.main(["--out", str(elsewhere)], registry={"B1": fake_case()}) == 2
    err = capsys.readouterr().err
    assert "--out" in err and str(out.resolve()) in err
    assert not elsewhere.exists() and not out.exists()


def test_a_run_holds_the_lock_of_its_root_for_the_whole_run_and_lets_go_after_it(monkeypatch, tmp_path):
    import bench
    import femm_queue
    from bench import run_bench
    use_config(monkeypatch, tmp_path)
    seen = []

    def case(folder, log):
        seen.append(femm_queue.holder(bench.root_lock_dir()))
        return {"case": "B1", "checks": [bench.check("a thing", True, "x", "x", "nowhere")]}
    assert run_bench.main([], registry={"B1": case, "B2": case}) == 0
    assert len(seen) == 2 and all(h and h["pid"] == os.getpid() and h["label"].startswith("bench run") for h in seen)
    assert "B1" in seen[0]["label"] and "B2" in seen[0]["label"]                                        # the label names the cases
    assert femm_queue.holder(bench.root_lock_dir()) is None                                             # let go after the run
    assert (bench.root_lock_dir() / "solver.lock").exists()                                             # and the lock file is never deleted (femm_queue's fourth fact)


def test_a_second_run_on_the_same_root_waits_and_prints_who_holds_it(tmp_path):
    import threading
    import time

    import bench
    first_in, let_go = threading.Event(), threading.Event()

    def first():
        with bench.hold_root("bench run first", folder=tmp_path):
            first_in.set()
            let_go.wait(20)
    holder_thread = threading.Thread(target=first)
    holder_thread.start()
    assert first_in.wait(10)
    said, got = [], []

    def second():
        with bench.hold_root("bench run second", folder=tmp_path, say=said.append, poll_s=0.05, look_s=0.01):
            got.append(time.monotonic())
    waiter = threading.Thread(target=second)
    waiter.start()
    time.sleep(0.5)
    assert not got and said                                                                              # still waiting, and it said so
    assert "waiting for the bench root" in said[0] and "'bench run first'" in said[0] and f"(pid {os.getpid()})" in said[0] and str(tmp_path) in said[0]
    released = time.monotonic()
    let_go.set()
    holder_thread.join(20)
    waiter.join(20)
    assert got and got[0] >= released                                                                    # it went on only after the first let go


def test_a_second_run_in_another_process_waits_for_the_first_and_names_its_pid(tmp_path):
    """The real case: two run_bench.py processes. A stand-in process holds the root for two seconds; this one asks and has to wait for it."""
    import time

    import bench
    code = ("import os, sys, time\nsys.path.insert(0, sys.argv[1])\nimport bench\n"
            "with bench.hold_root('bench run stand-in', folder=sys.argv[2]):\n    print('held', os.getpid(), flush=True)\n    time.sleep(2.0)\n")
    child = subprocess.Popen([sys.executable, "-I", "-c", code, str(P), str(tmp_path)], stdout=subprocess.PIPE, text=True)
    try:
        line = child.stdout.readline().split()
        assert line[0] == "held"
        said, started = [], time.monotonic()
        with bench.hold_root("bench run second", folder=tmp_path, say=said.append, poll_s=0.05, look_s=0.02):
            waited = time.monotonic() - started
        assert 0.5 < waited < 20 and said and f"(pid {line[1]})" in said[0] and "'bench run stand-in'" in said[0]
    finally:
        if child.poll() is None:
            child.kill()
        child.wait(30)


def test_ctrl_c_while_waiting_for_the_root_leaves_no_lock_behind(tmp_path):
    import threading

    import bench
    first_in, let_go = threading.Event(), threading.Event()

    def first():
        with bench.hold_root("bench run first", folder=tmp_path):
            first_in.set()
            let_go.wait(20)
    holder_thread = threading.Thread(target=first)
    holder_thread.start()
    assert first_in.wait(10)

    def interrupted(text):
        raise KeyboardInterrupt
    with pytest.raises(KeyboardInterrupt):
        with bench.hold_root("bench run second", folder=tmp_path, say=interrupted, poll_s=0.05, look_s=0.01):
            pytest.fail("the root was not free")
    let_go.set()
    holder_thread.join(20)
    with bench.hold_root("bench run third", folder=tmp_path, say=lambda text: pytest.fail("it had to wait: " + text), look_s=0.01):
        pass                                                                                             # free at once: the interrupted wait took nothing and left nothing


def test_a_run_waits_for_the_lock_of_its_root_before_it_empties_any_case_folder(monkeypatch, tmp_path, capsys):
    """The danger the lock is for: a second run on the same root would empty the first run's case folders (fresh()) while the first is still working in them.
    Another thread plays the first run: it holds the root for a second; the run under test (in this thread, as a bench run is the main thread of its process) has to wait,
    and while it waits the first run's file must still be there and no case may have run."""
    import threading
    import time

    import bench
    from bench import run_bench
    out = use_config(monkeypatch, tmp_path)
    (out / "bench" / "B1").mkdir(parents=True)
    old = out / "bench" / "B1" / "old.txt"
    old.write_text("the first run's")
    calls, while_waiting = [], []
    holder_in = threading.Event()

    def first_run():
        with bench.hold_root("bench run first", folder=out / "bench"):
            holder_in.set()
            time.sleep(1.0)
            while_waiting.append((old.exists(), list(calls)))                                           # what the waiting run has done meanwhile: nothing
    thread = threading.Thread(target=first_run)
    thread.start()
    assert holder_in.wait(10)
    started = time.monotonic()
    assert run_bench.main(["--only", "B1"], registry={"B1": fake_case(calls=calls)}) == 0
    thread.join(30)
    assert while_waiting == [(True, [])]                                                                 # waiting: the old folder was still there and no case had run
    assert time.monotonic() - started >= 0.9 and calls == ["B1"] and not old.exists()                    # it went on after the first let go, and then emptied the folder
    assert "waiting for the bench root" in capsys.readouterr().out


def test_the_root_lock_is_the_queues_own_and_not_a_copy_of_it():
    """'Through the same operating-system mechanism as femm_queue's slot (reuse it, don't copy it)': bench/ takes the lock with femm_queue.solver_slot(lock_dir=...) in one place
    and has no locking code of its own."""
    from test_check_tools import scopes_calling
    init = (P / "bench" / "__init__.py").read_text(encoding="utf-8")
    assert scopes_calling(init, "femm_queue", "solver_slot") == ["hold_root"]
    for file in sorted((P / "bench").glob("*.py")):
        text = file.read_text(encoding="utf-8")
        assert "msvcrt" not in text and "fcntl" not in text and "LK_NBLCK" not in text, f"{file.name} has locking code of its own"


def test_the_watcher_counts_the_most_children_of_one_name_alive_at_one_look():
    """B9's 'at most 8 workers' is read from outside: the most processes of one name that were alive, under the bench, at any one look."""
    import time

    import bench
    kids = []
    with bench.watch_children() as watcher:
        for _ in range(3):
            kids.append(subprocess.Popen([sys.executable, "-I", "-c", "import time; time.sleep(2.0)"], creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS))
        try:
            time.sleep(1.0)
        finally:
            for kid in kids:
                kid.wait(30)
    peaks = watcher.peaks()
    assert peaks["python.exe"] >= 3, peaks                       # three children (each may be a launcher and the interpreter below it: so at least three)
    assert all(isinstance(v, int) and v >= 1 for v in peaks.values())


# ----------------------------------------------------------------------------------------------- 4. the entry point states what to expect (final review, Important 3)
def test_run_bench_states_the_numbers_to_expect_and_that_a_run_with_b2_exits_with_status_1():
    """A3 of the final review. The plan's rule is that every new script says the numbers to expect; run_bench.py, the bench's one entry point, did not, and its first
    whole run ends 'bench: FAIL, 3 of 119 gating checks fail: B2 (3)' with exit status 1, by ruling (the chord error of path 1's mesh). The README's Files table has a row
    for it too: the command, --only, --out, --mesh-study, summary.json and the same exit status."""
    from bench import run_bench
    doc = run_bench.__doc__
    assert "Numbers to expect" in doc and "run_bench.py --mesh-study" in doc and "2026-11-06" in doc
    for text in ("141 checks", "119 gating", "116 pass", "22 reported", "18 pass, 3 fail, 1 skipped", "784 s",
                 "bench: FAIL, 3 of 119 gating checks fail: B2 (3); 22 reported, not gating (18 pass, 3 fail, 1 skipped) (784 s)",
                 "Exit status 1 whenever B2 runs", "chord error"):
        assert text in doc, text
    for case in ("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9"):
        assert re.search(rf"\b{case}\b.*\bs\b", doc.split("Numbers to expect", 1)[1]), case               # each case's checks and seconds are in the table
    readme = (P / "README.md").read_text(encoding="utf-8")
    rows = [line for line in readme.splitlines() if line.startswith("| `bench/run_bench.py`") or line.startswith("| `bench/`")]
    assert len(rows) == 1, rows
    for text in ("uv run python bench/run_bench.py", "--only", "--out", "--mesh-study", "summary.json", "exit status 1", "B2", "tests/test_bench"):
        assert text in rows[0], text


def test_no_file_of_the_bench_or_its_tests_points_at_the_interfaces_file_of_the_plans_workspace():
    """T13-P5. The interfaces file of Task 13, a file of the plan's own workspace (on this PC only, and not in git), was named as the source of a rule in three places, so a
    reader of the code was sent to a file that is not in the repository. The rule is written where it applies."""
    name = "task-13-" + "interfaces"                                                          # (built from two pieces: this test does not name it itself)
    offenders = []
    for folder in ("bench", "tests", "checks"):
        for file in sorted((P / folder).glob("*.py")):
            for number, line in enumerate(file.read_text(encoding="utf-8").splitlines(), start=1):
                if name in line:
                    offenders.append(f"{folder}/{file.name}:{number}: {line.strip()[:100]}")
    assert not offenders, "a file of the plan's workspace is named in the code:\n" + "\n".join(offenders)
