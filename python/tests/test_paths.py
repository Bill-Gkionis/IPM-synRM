"""Tests for paths.py and config.toml: every tool path and every machine folder comes from one file.

Run from P:   uv run pytest tests/test_paths.py -q      (a few seconds; no FEMM is started)

Four groups of tests:
  1. paths.py on its own: where a relative path points, what a name that is missing or empty says, what a broken file says.
  2. config.toml as it is on this PC: every tool it names is really there (this is the test that notices a program that was
     moved or deleted), and the values other scripts read from it.
  3. Who reads it. femm_session.femm_dir() takes [tools] femm first and FEMM_DIR only if the config names no FEMM (R4 of the
     tools plan). The scripts that read the reference model, and materials/common.py, take their paths from the config, and
     MOTOR_FIGURES_DIR still wins over vault_figures.
  4. R5, the rule behind all of it: no script types a drive path of its own. The scan (drive_path_literals) is tried on a
     small source first, so that a scan which finds nothing is known to be able to find something ("Self-check").
"""
import ast
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

import femm_session
import paths

P = Path(__file__).resolve().parent.parent
TOOLS = {"femm", "ccx", "cgx", "gmsh", "prepomax", "elmer", "elmergrid", "matlab"}      # what config.toml must name


def write_config(folder: Path, text: str) -> Path:
    """A config.toml of the test's own, so a wrong one can be tried without touching the real one."""
    config = folder / "config.toml"
    config.write_text(text, encoding="utf-8")
    return config


def load_fresh(name: str, file: Path):
    """The file as a module object of its own: not the one other tests may have imported, or changed (hand_reference_check
    sets validate_full.FEM; check_vs_fem reads its model from pytest's command line when it is imported under pytest)."""
    spec = importlib.util.spec_from_file_location(name, file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ----------------------------------------------------------------------------------------------- 1. paths.py on its own
def test_relative_path_resolves_against_project():
    assert paths.path("out") == P / "out"
    assert paths.path("femm42_out") == P / "femm42-out"
    assert paths.path("reference_fem") == P / "reference" / "IPM synRM analysis.FEM"


def test_the_default_config_is_config_toml_next_to_paths_py():
    assert paths.P == P and paths.CONFIG == P / "config.toml" and paths.CONFIG.is_file()
    assert Path(paths.__file__).resolve().parent == P


def test_a_relative_path_starts_at_the_project_whatever_the_current_folder_is(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    assert paths.path("out") == P / "out"


def test_an_absolute_path_is_returned_as_written_and_need_not_exist(tmp_path):
    far = tmp_path / "elsewhere" / "not made yet"
    config = write_config(tmp_path, f"[paths]\nfar = '{far}'\nnear = 'sub\\dir'\n")
    assert paths.path("far", config) == far
    assert paths.path("near", config) == P / "sub" / "dir"            # relative entries start at P even in another config
    assert not far.exists()


def test_setting_returns_the_value_with_the_type_toml_gave_it():
    assert paths.setting("femm", "workers") == 8 and isinstance(paths.setting("femm", "workers"), int)
    assert paths.setting("femm", "priority") == "below_normal"


def test_each_call_reads_the_config_it_is_given(tmp_path):
    here = tmp_path / "program.exe"
    here.write_text("x")
    config = write_config(tmp_path, f"[tools]\nmine = '{here}'\n[paths]\nout = 'other out'\n[femm]\nworkers = 3\n")
    assert paths.tool("mine", config) == here
    assert paths.path("out", config) == P / "other out"
    assert paths.setting("femm", "workers", config) == 3


def test_missing_tool_names_itself(tmp_path):
    config = write_config(tmp_path, "[tools]\nelmer = 'Z:\\nowhere.exe'\n")
    with pytest.raises(FileNotFoundError, match="elmer") as caught:
        paths.tool("elmer", config)
    assert "Z:\\nowhere.exe" in str(caught.value) and "config.toml" in str(caught.value)      # and where it looked, and in which file


def test_a_tool_left_empty_in_the_config_is_reported_as_missing(tmp_path):
    config = write_config(tmp_path, "[tools]\nccx = ''\n")
    with pytest.raises(FileNotFoundError, match="ccx"):
        paths.tool("ccx", config)


def test_a_tool_the_config_does_not_name_is_reported_as_missing_and_the_names_it_has_are_listed(tmp_path):
    """One kind of error for every way a tool can be missing (no line, an empty line, a line that leads nowhere), so a
    caller that goes through the tools needs one `except FileNotFoundError` and never shows a traceback."""
    config = write_config(tmp_path, "[tools]\nccx = ''\n")
    with pytest.raises(FileNotFoundError, match="ccz.*ccx"):            # a typo: the message names it and what could have been meant
        paths.tool("ccz", config)
    with pytest.raises(FileNotFoundError, match="ccz"):                 # and a config with no [tools] section at all
        paths.tool("ccz", write_config(tmp_path, "[paths]\nout = 'out'\n"))


def test_a_path_or_setting_the_file_does_not_have_is_a_key_error_that_lists_the_names_it_has(tmp_path):
    config = write_config(tmp_path, "[tools]\nccx = ''\n[paths]\nout = 'out'\n[femm]\nworkers = 8\n")
    with pytest.raises(KeyError, match="nowhere.*out"):
        paths.path("nowhere", config)
    with pytest.raises(KeyError, match="threads.*workers"):
        paths.setting("femm", "threads", config)
    with pytest.raises(KeyError, match="nosection"):
        paths.setting("nosection", "x", config)


def test_an_empty_path_entry_is_refused_not_taken_for_the_project_folder(tmp_path):
    config = write_config(tmp_path, "[paths]\nout = ''\n")
    with pytest.raises(ValueError, match="out"):
        paths.path("out", config)


def test_a_file_that_is_not_toml_is_reported_with_its_name(tmp_path):
    config = write_config(tmp_path, "[tools\nccx = \n")
    with pytest.raises(ValueError, match="config.toml"):
        paths.setting("femm", "workers", config)


def test_a_config_file_that_is_not_there_is_a_file_not_found_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="nothing.toml"):
        paths.path("out", tmp_path / "nothing.toml")


def test_running_paths_py_lists_the_names_and_exits_1_when_a_tool_is_missing(tmp_path, capsys):
    assert paths.main([]) == 0
    assert "8 of 8 tools found" in capsys.readouterr().out
    config = write_config(tmp_path, "[tools]\nelmer = 'Z:\\nowhere.exe'\ngmsh = ''\n")
    assert paths.main([str(config)]) == 1
    shown = capsys.readouterr().out
    assert "elmer" in shown and "gmsh" in shown and "MISSING" in shown and "0 of 2 tools found" in shown


# ----------------------------------------------------------------------------------------------- 2. config.toml on this PC
def test_configured_tools_exist():
    tools = tomllib.loads(paths.CONFIG.read_text(encoding="utf-8"))["tools"]
    assert set(tools) == TOOLS, "a tool was added to or dropped from config.toml: change TOOLS here too"
    for name in tools:
        assert paths.tool(name).exists(), name


def test_the_other_files_the_config_names_are_there():
    assert paths.path("femm_matlib").is_file()
    assert paths.path("vault_figures").is_dir()
    assert paths.path("reference_fem").is_file()
    assert paths.path("femm42_out").is_dir()


# ----------------------------------------------------------------------------------------------- 3. who reads the config
def test_femm_dir_comes_from_config(monkeypatch, tmp_path):
    """R4: with FEMM_DIR set to a wrong folder, femm_dir() still returns [tools] femm."""
    for wrong in (tmp_path, tmp_path / "no such folder"):               # a folder with no FEMM in it, and one that is not there
        monkeypatch.setenv("FEMM_DIR", str(wrong))
        assert femm_session.femm_dir() == paths.tool("femm")
    monkeypatch.delenv("FEMM_DIR")
    assert femm_session.femm_dir() == paths.tool("femm")


def test_femm_dir_uses_femm_dir_only_when_the_config_names_no_femm(monkeypatch, tmp_path):
    real = paths.tool("femm")
    for text in ("[tools]\nfemm = ''\n", "[tools]\nccx = ''\n", "[paths]\nout = 'out'\n"):      # an empty line, no line, no section
        monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, text))
        monkeypatch.setenv("FEMM_DIR", str(real))
        assert femm_session.femm_dir() == real                          # as before the config existed
        monkeypatch.setenv("FEMM_DIR", str(tmp_path))                   # a wrong FEMM_DIR is still refused, and named
        with pytest.raises(FileNotFoundError, match="FEMM_DIR"):
            femm_session.femm_dir()


def test_a_wrong_femm_line_in_the_config_is_an_error_and_not_a_fall_back_to_femm_dir(monkeypatch, tmp_path):
    """A good FEMM_DIR must not hide a stale config line: the error names the config, so the line gets fixed."""
    monkeypatch.setenv("FEMM_DIR", str(paths.tool("femm")))
    for wrong in (tmp_path / "no such folder", tmp_path):               # not there at all, and there but not FEMM's
        monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, f"[tools]\nfemm = '{wrong}'\n"))
        with pytest.raises(FileNotFoundError, match="config.toml") as caught:
            femm_session.femm_dir()
        assert "femm" in str(caught.value)


def test_the_three_scripts_that_read_the_reference_model_read_the_configured_one(monkeypatch):
    """M7 of the final review: check_vs_fem.py, validate_full.py and materials/common.py used to name a file in the
    IPM-synRM clone on the Desktop. Now they take the copy in P\\reference from config.toml."""
    monkeypatch.setattr(sys, "argv", ["check_vs_fem.py"])               # check_vs_fem takes its model from the command line, if given
    reference = paths.path("reference_fem")
    assert reference.is_file()
    assert Path(load_fresh("check_vs_fem_under_test", P / "check_vs_fem.py").FEM) == reference
    assert Path(load_fresh("validate_full_under_test", P / "validate_full.py").FEM) == reference
    assert load_fresh("common_under_test", P / "materials" / "common.py").REFERENCE_FEM == reference


def test_the_materials_scripts_take_the_library_and_the_notes_folder_from_the_config():
    cm = load_fresh("common_under_test", P / "materials" / "common.py")
    assert cm.MATLIB == paths.path("femm_matlib")
    assert cm.FIGURES == paths.path("vault_figures")


def test_a_figure_is_copied_into_the_vault_only_when_asked(monkeypatch, tmp_path):
    """138 of the final review. save() copies a figure to the folder MOTOR_FIGURES_DIR names, to the config's vault_figures if the variable says `vault`, and nowhere if it is
    unset (the rule: no run writes into the project notes; figures are copied by hand). config.toml's vault_figures stays: it is the place used when asked.

    The copy is replaced by a recorder, so that this test cannot write into the real notes' folder even if it fails."""
    cm = load_fresh("common_under_test", P / "materials" / "common.py")
    monkeypatch.setattr(cm, "OUT", tmp_path / "out")
    (tmp_path / "out").mkdir()
    copies = []
    monkeypatch.setattr(cm.shutil, "copy2", lambda source, destination: copies.append(Path(destination)))
    plt = cm.style()
    fig, ax = plt.subplots(figsize=(1, 1))
    ax.plot([0, 1], [0, 1])
    try:
        monkeypatch.setenv("MOTOR_FIGURES_DIR", str(tmp_path / "mine"))
        cm.save(fig, "t.png")
        monkeypatch.delenv("MOTOR_FIGURES_DIR")
        cm.save(fig, "t.png")                                          # unset: no copy
        monkeypatch.setenv("MOTOR_FIGURES_DIR", "")
        cm.save(fig, "t.png")                                          # empty: no copy
        monkeypatch.setenv("MOTOR_FIGURES_DIR", "vault")
        cm.save(fig, "t.png")                                          # asked for: the config's vault_figures
    finally:
        plt.close(fig)
    assert copies == [tmp_path / "mine" / "t.png", paths.path("vault_figures") / "t.png"]


# ----------------------------------------------------------------------------------------------- 4. R5: no drive path in a script
# A drive path is a single letter, a colon and a slash: C:\Users, D:/Apps. The letter must not be the end of a longer word,
# or "http://" and a regular expression like r"workers:\s+" would count (the "p:/" of http, the "s:\" of workers).
DRIVE_PATH = re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]")
RETIRED = re.compile("|".join(("Desk" + "top", "IPM-" + "synRM", "python-to-" + "femm", r"scripts[\\/]auto" + "mation")))
NOT_SCANNED = ("learn/", "tests/", "envs/")           # R5: teaching examples and tests may show machine paths; scripts may not


def literals_matching(source: str, pattern: re.Pattern) -> list[tuple[int, str]]:
    """(line, text) of every string literal in `source` that `pattern` finds, except docstrings.

    A docstring is the first statement of a module, class or function when it is a string: it is documentation, and may
    show an example path. Everything else is code: a path inside an f-string counts, a comment is not in the tree at all.
    """
    tree = ast.parse(source)
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            first = node.body[0] if node.body else None
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                docstrings.add(id(first.value))
    return [(node.lineno, node.value) for node in ast.walk(tree)
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings
            and pattern.search(node.value)]


def drive_path_literals(source: str) -> list[tuple[int, str]]:
    return literals_matching(source, DRIVE_PATH)


def tracked_python_files(folder: Path = P) -> list[str]:
    """The .py files git tracks in `folder` (P: so scratch scripts under out\\ never count), as paths relative to it."""
    if shutil.which("git") is None:
        pytest.skip("git is not installed, so the tracked files cannot be listed")
    try:
        # -z ends each name with a NUL instead of a newline, so a name with unusual characters is never quoted away; git writes the bytes of the name as they are
        # (UTF-8), so they are decoded as UTF-8 and not with the locale's code page, which would turn a name with a non-ASCII letter into another name and drop it
        listing = subprocess.run(["git", "ls-files", "-z", "*.py"], cwd=folder, capture_output=True, encoding="utf-8", check=True).stdout
    except subprocess.CalledProcessError:
        pytest.skip("this folder is not inside a git repository, so there is no list of tracked files")
    return [name for name in listing.split("\0") if name and (folder / name).is_file()]


SAMPLE = r'''
r"""A module docstring may show a path: D:\Apps\femm-qt\bin."""
import re


def f():
    """So may a function's: C:/Users/user/x."""
    return r"C:\Users\user\x.FEM"
BAD_SLASH = 'D:/Apps/tool.exe'
BAD_F = f"{BAD_SLASH}E:\\folder"
URL = "http://example.org/a"
WORKERS = re.compile(r"workers:\s+([\d.]+) s")
NOTE = "Note: see the next line"
'''


def test_the_scan_finds_drive_paths_in_code_and_spares_docstrings_urls_and_patterns():
    """Self-check: which literals of SAMPLE does the scan flag? Three: the raw string with C:\\ in f(), the string with
    D:/ and the f-string's literal part with E:\\. Not flagged: the module's and the function's docstrings (documentation),
    "http://..." (its p:/ is the end of a word, not a drive letter), the pattern with workers:\\s (same reason), and "Note: see"
    (no slash after the colon)."""
    found = [text for _, text in drive_path_literals(SAMPLE)]
    assert sorted(found) == sorted(["C:\\Users\\user\\x.FEM", "D:/Apps/tool.exe", "E:\\folder"])


def test_no_hard_coded_drive_paths():
    """R5: not one string literal that is not a docstring, in a tracked .py file outside learn/, tests/ and envs/, holds a
    drive path. config.toml is the one place for machine paths."""
    found = []
    for name in tracked_python_files():
        if name.startswith(NOT_SCANNED):
            continue
        for line, text in drive_path_literals((P / name).read_text(encoding="utf-8")):
            found.append(f"{name}:{line}: {text[:90]!r}")
    assert not found, "a drive path is typed into a script; it belongs in config.toml (paths.path / paths.tool):\n" + "\n".join(found)


def test_no_script_names_a_place_that_is_retired():
    """The IPM-synRM clone on the Desktop is restructured, and D:\\dev\\scripts\\automation and D:\\dev\\python-to-femm are gone.
    No tracked .py file (tests included) may read from them: they appear in docstrings as history only."""
    found = []
    for name in tracked_python_files():
        for line, text in literals_matching((P / name).read_text(encoding="utf-8"), RETIRED):
            found.append(f"{name}:{line}: {text[:90]!r}")
    assert not found, "a script reads from a place that no longer exists or is being restructured:\n" + "\n".join(found)


def test_a_tracked_file_with_a_non_ascii_name_is_in_the_list_the_drive_path_scan_reads(tmp_path):
    """B4 of the final review (Task 2's minor): the listing was decoded with the locale's code page, so a .py file whose name has a letter outside ASCII came back as another
    name, was not found on disk, and dropped out of the R5 scan without a word. Shown in a repository made for the test."""
    if shutil.which("git") is None:
        pytest.skip("git is not installed")
    name = "pr\u00fcfung.py"
    (tmp_path / name).write_text("x = 1\n", encoding="utf-8")
    (tmp_path / "plain.py").write_text("y = 2\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "."], cwd=tmp_path, check=True)
    assert sorted(tracked_python_files(tmp_path)) == sorted([name, "plain.py"])


def test_motorgeos_docstring_names_the_reference_model_where_it_is_now():
    """B4 of the final review (Task 2's minor): the restructure moved FEMM/models/ to femm/ in the IPM-synRM repository, and motorgeo.py's first lines still named the old place.
    The byte-exact copy that the code reads is reference/IPM synRM analysis.FEM."""
    import motorgeo
    assert "FEMM/models" not in motorgeo.__doc__ and "reference/IPM synRM analysis.FEM" in motorgeo.__doc__


def test_pytests_temporary_folders_are_under_out_dot_pytest_unless_the_user_chose_a_place(monkeypatch, tmp_path):
    """149 of the final review. pytest keeps ~300 files (1.2 MB) of every run in the user's Temp on C:, three runs at a time. conftest.py sets PYTEST_DEBUG_TEMPROOT, the root pytest
    takes its numbered temp folders from, to out\\.pytest (git-ignored, on D:), with setdefault, so a user's own setting wins. Not --basetemp: that empties its folder at every start
    and would clobber a second pytest run going at the same time; the numbered folders of the temp root do not. The folder has to exist (pytest makes pytest-of-<user> inside it, but not it).

    The name starts with a dot on purpose: the stale scan (rulings.results, which stale.py runs on out\\) does not look into a folder whose name starts with a dot, at any depth,
    as it does not look into out\\.femm. A folder called out\\pytest was scanned, and the tests' own csv and json files (about 300 for the three runs pytest keeps) were listed as
    unstamped results next to the studies'. The last lines show that."""
    import conftest
    import rulings
    config = tmp_path / "config.toml"
    config.write_text(f"[paths]\nout = '{tmp_path / 'out'}'\n", encoding="utf-8")
    monkeypatch.setattr(paths, "CONFIG", config)
    monkeypatch.delenv("PYTEST_DEBUG_TEMPROOT", raising=False)
    root = conftest.temp_root_under_out()
    assert root == tmp_path / "out" / ".pytest" and root.is_dir() and os.environ["PYTEST_DEBUG_TEMPROOT"] == str(root)
    monkeypatch.setenv("PYTEST_DEBUG_TEMPROOT", str(tmp_path / "mine"))
    assert conftest.temp_root_under_out() is None and os.environ["PYTEST_DEBUG_TEMPROOT"] == str(tmp_path / "mine"), "a user's own setting has to win"
    leftovers = root / "pytest-of-someone" / "pytest-1"
    leftovers.mkdir(parents=True)
    (leftovers / "result.csv").write_text("a\n1\n", encoding="utf-8")
    (tmp_path / "out" / "study.csv").write_text("a\n1\n", encoding="utf-8")
    assert [path.name for path in rulings.results(tmp_path / "out")] == ["study.csv"], "the stale scan looked into pytest's temp root"


def test_this_tests_own_temporary_folder_is_under_the_temp_root_conftest_chose(tmp_path):
    root = Path(os.environ["PYTEST_DEBUG_TEMPROOT"]).resolve()
    assert root in tmp_path.resolve().parents, f"{tmp_path} is not under {root}"
