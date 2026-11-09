"""Tests for materials/common.py: where the materials scripts read FEMM's library and where they write figures.

  - Figures: every materials script saves through common.save(), which used to copy each figure into the notes' folder
    unless MOTOR_FIGURES_DIR said otherwise. Since the final review of the tools plan the copy is opt-in (the rule: no
    run writes into the project notes; figures are copied by hand): unset or empty, a figure stays in out/materials; MOTOR_FIGURES_DIR=<folder>
    copies it there too; MOTOR_FIGURES_DIR=vault copies it into the notes' folder. These tests point the notes' folder at a temporary
    one too, so even a failing test cannot write into the real notes' folder.
  - Library: FEMM 5.0's material library is the one the scripts read now. Its "M-19 Steel" is gone (FEMM 4.2's block is
    kept in motorgeo.M19_STEEL_FEMM42), and the other blocks the notes use are the same in both libraries.

Run from P:   uv run pytest tests/test_materials_common.py -q
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

import paths

P = Path(__file__).resolve().parent.parent
LIBRARY_42 = Path(r"D:\Apps\femm42\bin\matlib.dat")        # FEMM 4.2's, the old program: not one of config.toml's tools, kept only to compare with
LIBRARY_50 = paths.path("femm_matlib")                      # FEMM 5.0's, from config.toml; the test below checks its blocks, so it is FEMM 5.0's


@pytest.fixture(scope="module")
def cm():
    """materials/common.py as a module. (It is a script helper, imported as `common` from its own folder.)"""
    spec = importlib.util.spec_from_file_location("materials_common", P / "materials" / "common.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def sandbox(cm, tmp_path, monkeypatch):
    """Point both of save()'s destinations at temporary folders: nothing in a test can reach the real notes' folder."""
    monkeypatch.setattr(cm, "OUT", tmp_path / "out")
    monkeypatch.setattr(cm, "FIGURES", tmp_path / "vault")
    (tmp_path / "out").mkdir()
    return tmp_path


def a_figure(cm):
    plt = cm.style()
    fig, ax = plt.subplots(figsize=(1, 1))
    ax.plot([0, 1], [0, 1])
    return fig, plt


def test_save_copies_to_motor_figures_dir_and_not_to_the_vault(cm, sandbox, monkeypatch):
    monkeypatch.setenv("MOTOR_FIGURES_DIR", str(sandbox / "mine"))
    fig, plt = a_figure(cm)
    cm.save(fig, "t.png")
    plt.close(fig)
    assert (sandbox / "out" / "t.png").is_file()              # the numbers' folder always gets it
    assert (sandbox / "mine" / "t.png").is_file()             # the requested folder gets the copy (and is made if missing)
    assert not (sandbox / "vault").exists(), "the notes' folder must not even be created"


@pytest.mark.parametrize("value", [None, ""])
def test_save_copies_nowhere_but_out_when_the_variable_is_unset_or_empty(cm, sandbox, monkeypatch, value):
    """138 of the final review: unset means no copy into the notes' folder (it used to). The figure is in the numbers' folder, and nothing else is made."""
    if value is None:
        monkeypatch.delenv("MOTOR_FIGURES_DIR", raising=False)
    else:
        monkeypatch.setenv("MOTOR_FIGURES_DIR", value)
    fig, plt = a_figure(cm)
    cm.save(fig, "t.png")
    plt.close(fig)
    assert (sandbox / "out" / "t.png").is_file()
    assert not (sandbox / "vault").exists(), "a figure was copied into the notes' folder although nobody asked"
    assert sorted(p.name for p in sandbox.iterdir()) == ["out"]


@pytest.mark.parametrize("word", ["vault", "Vault", "VAULT"])
def test_save_copies_into_the_notes_folder_only_when_the_variable_says_vault(cm, sandbox, monkeypatch, word):
    """The one way a run writes into the notes' folder, on purpose (regenerating the guide's figures): MOTOR_FIGURES_DIR=vault names the notes' folder (config.toml, vault_figures),
    whatever the case of the word. A folder that happens to be called vault is written ./vault."""
    monkeypatch.setenv("MOTOR_FIGURES_DIR", word)
    fig, plt = a_figure(cm)
    cm.save(fig, "t.png")
    plt.close(fig)
    assert (sandbox / "vault" / "t.png").is_file() and (sandbox / "out" / "t.png").is_file()


def test_a_folder_called_vault_is_asked_for_with_a_path(cm, sandbox, monkeypatch):
    elsewhere = sandbox / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    monkeypatch.setenv("MOTOR_FIGURES_DIR", "./vault")
    fig, plt = a_figure(cm)
    cm.save(fig, "t.png")
    plt.close(fig)
    assert (elsewhere / "vault" / "t.png").is_file() and not (sandbox / "vault").exists()


def test_the_variable_is_read_when_save_is_called(cm, sandbox, monkeypatch):
    """Set after the module was imported, it still counts."""
    fig, plt = a_figure(cm)
    monkeypatch.setenv("MOTOR_FIGURES_DIR", str(sandbox / "late"))
    cm.save(fig, "t.png")
    plt.close(fig)
    assert (sandbox / "late" / "t.png").is_file() and not (sandbox / "vault").exists()


@pytest.mark.parametrize("spelling", ["out", "out/../out"])
def test_the_figures_folder_may_be_the_numbers_folder_itself(cm, sandbox, monkeypatch, spelling):
    """MOTOR_FIGURES_DIR = out/materials names the folder save() has just written the figure into: copying a file onto
    itself used to raise shutil.SameFileError. However the folder is spelled, the figure is simply left where it is."""
    monkeypatch.setenv("MOTOR_FIGURES_DIR", str(sandbox / spelling))
    fig, plt = a_figure(cm)
    path = cm.save(fig, "t.png")
    plt.close(fig)
    assert path == sandbox / "out" / "t.png" and path.is_file() and path.stat().st_size > 0
    assert sorted(p.name for p in (sandbox / "out").iterdir()) == ["t.png"]
    assert not (sandbox / "vault").exists()


@pytest.mark.skipif(not LIBRARY_50.is_file(), reason="FEMM 5.0 is not installed")
def test_matlib_is_femm5s_library_and_has_no_m19_steel(cm):
    assert Path(cm.MATLIB) == LIBRARY_50
    names = {b["BlockName"] for b in cm.read_matlib()}
    assert {"Air", "N45", "M-19 14 mil"} <= names
    assert "M-19 Steel" not in names                          # why motorgeo carries FEMM 4.2's block itself


@pytest.mark.skipif(not (LIBRARY_50.is_file() and LIBRARY_42.is_file()), reason="needs both FEMM libraries")
def test_the_blocks_the_notes_use_are_the_same_in_both_libraries(cm):
    old, new = cm.read_matlib(LIBRARY_42), cm.read_matlib(LIBRARY_50)
    asm = lambda lib: next(b for b in lib if b["BlockName"] == "M-19" and "Metals Handbook" in b["folder"])
    for pick in (lambda lib: next(b for b in lib if b["BlockName"] == "Air"),
                 lambda lib: next(b for b in lib if b["BlockName"] == "N45"), asm):
        a, b = pick(old), pick(new)
        assert {k: v for k, v in a.items() if k not in ("bh", "folder")} == {k: v for k, v in b.items() if k not in ("bh", "folder")}
        assert np.array_equal(a["bh"], b["bh"]), a["BlockName"]
    assert len(asm(new)["bh"]) > 20                           # a real curve, so the comparison above compared something


@pytest.mark.skipif(not LIBRARY_42.is_file(), reason="FEMM 4.2 is not installed")
def test_m19_steel_block_has_the_shape_read_matlib_returns(cm):
    from_library = next(b for b in cm.read_matlib(LIBRARY_42) if b["BlockName"] == "M-19 Steel")
    built = cm.m19_steel_block()
    assert set(built) == set(from_library)
    assert built["folder"] == from_library["folder"] == "Soft Magnetic Materials / Silicon Iron"
    for key in from_library:
        if key == "bh":
            assert np.array_equal(built["bh"], from_library["bh"]) and built["bh"].shape == (47, 2)
        else:
            assert built[key] == from_library[key], key
