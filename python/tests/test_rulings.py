"""Tests for rulings.py and stale.py: the rulings file's loader, the stamp a result carries, and the list of what is stale.

Run from P:   uv run pytest tests/test_rulings.py -q      (a few seconds; no solver is started; git is asked a few dozen times)

The idea they pin (the idea of orchestration): a result says which rulings it was computed from, so that when one ruling
changes, the results that used it can be listed and re-run. Five groups of tests:
  1. The brief's six (T5 of the Tools workstream), in the brief's words.
  2. load() and value(): the format, the keys flattened by dots, and what a wrong file says (a bad status, a missing field, anything
     that is not a ruling, an empty table, a file that is not UTF-8 without a BOM). Where the file is when nobody says (config.toml's
     [paths] rulings), and what a missing one says.
  3. stamp(): the fields of the file it writes, the commit of the code (R6 of the tools plan), `dirty` for the code folder and not for
     the rest of the repository, what it says when git cannot tell, and its refusals (one string for the keys, a result that is not
     there, a key the rulings do not have).
  4. stale(): a fresh result is not listed, every kind of result is looked at (folders whose names start with a dot excepted), a stamp
     is not a result, a stamp that cannot be read is no stamp, and True is not 1.
  5. stale.py, the command line around stale(): what it prints and the exit status (0 nothing stale, 1 something is, 2 it could not look).
They use tests/fixtures/rulings_min.toml (three quantities) and a copy of it in a temporary folder, which a test edits.
The real rulings.toml is written by the orchestration session (O-1) and does not exist yet; no test needs it.
A test that compares with the commit of the code needs P to be inside a git checkout that has a commit. When it is not (a zip export of
the code, say) that test is skipped, with the reason, and not failed: a stamp then has no commit by design.
"""
import dataclasses
import json
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import pytest

import paths
import rulings
import stale as stale_cli

P = Path(__file__).resolve().parent.parent
FIX = Path(__file__).with_name("fixtures") / "rulings_min.toml"


# ----------------------------------------------------------------------------------------------- helpers
def copy_of_fixture(folder: Path) -> Path:
    """A rulings file of the test's own: the fixture copied. The test may change it; the fixture never changes."""
    copy = folder / "rulings.toml"
    shutil.copy(FIX, copy)
    return copy


def edit(path: Path, old: str, new: str) -> None:
    """Replace `old` by `new` in a text file, which has to hold it exactly once (a test that edits the wrong line must fail, not pass)."""
    text = path.read_text(encoding="utf-8")
    assert text.count(old) == 1, f"{old!r} is in {path} {text.count(old)} times"
    path.write_text(text.replace(old, new), encoding="utf-8")


def drop_table(path: Path, name: str) -> None:
    """Delete the table [name] (its header and its lines, up to the next blank line) from a rulings file."""
    text = path.read_text(encoding="utf-8")
    cut = re.sub(rf"^\[{re.escape(name)}\]\n.*?\n\n", "", text, count=1, flags=re.S | re.M)
    assert cut != text, f"there is no table [{name}] in {path}"
    path.write_text(cut, encoding="utf-8")


def drop_line(path: Path, field: str) -> None:
    """Delete the first line that starts with `field =` (in the fixture that is the bridge's: it comes first)."""
    text = path.read_text(encoding="utf-8")
    cut = re.sub(rf"^{field}\s*=.*\n", "", text, count=1, flags=re.M)
    assert cut != text, f"there is no line for {field!r} in {path}"
    path.write_text(cut, encoding="utf-8")


def make_result(folder: Path, name: str = "a.csv") -> Path:
    """A result file, as a study would write it (the stamp comes after it, never before)."""
    result = folder / name
    result.parent.mkdir(parents=True, exist_ok=True)
    result.write_text("x\n1\n", encoding="utf-8")
    return result


def write_config(folder: Path, text: str) -> Path:
    """A config.toml of the test's own, to point paths.path("rulings") and paths.path("out") somewhere else."""
    config = folder / "config.toml"
    config.write_text(text, encoding="utf-8")
    return config


def require_git() -> None:
    if shutil.which("git") is None:
        pytest.skip("git is not installed, so there is no commit to compare with")


def head_of_the_code() -> str:
    """`git rev-parse HEAD` run in P: what R6 says the stamp's `commit` is. Skips the test, with the reason, when P is not in a git
    checkout with a commit (a zip export of the code, say): there is then no commit to compare with, which is no failure of the stamp."""
    require_git()
    done = subprocess.run(["git", "rev-parse", "HEAD"], cwd=paths.P, capture_output=True, text=True)
    if done.returncode != 0:
        why = done.stderr.strip().splitlines()[0] if done.stderr.strip() else "no message"
        pytest.skip(f"{paths.P} is not inside a git checkout that has a commit (git said: {why}), so there is no commit to compare with")
    return done.stdout.strip()


def make_repo(folder: Path) -> str:
    """A git repository of the test's own, with one commit; returns that commit. (Nothing about the real repository is touched.)

    It has the shape of Storage in small: a folder `code` (this project's) and a folder `other` (anything else the repository holds), a
    file in each, and a .gitignore that ignores *.log. The tests ask about `code`."""
    def git(*args: str) -> str:
        return subprocess.run(["git", "-c", "user.name=test", "-c", "user.email=test@example.org", "-c", "commit.gpgsign=false", *args],
                              cwd=folder, capture_output=True, text=True, check=True).stdout.strip()
    git("init", "-q")
    (folder / "code").mkdir()
    (folder / "other").mkdir()
    (folder / "code" / "a.txt").write_text("a", encoding="utf-8")
    (folder / "other" / "b.txt").write_text("b", encoding="utf-8")
    (folder / ".gitignore").write_text("*.log\n", encoding="utf-8")
    git("add", ".")
    git("commit", "-q", "-m", "first")
    return git("rev-parse", "HEAD")


def stamped_then_changed(folder: Path) -> tuple[Path, Path]:
    """a.csv stamped with the bridge ruling (step EM-1) and the ruling then changed from 1.4 to 1.2; b.csv never stamped. Returns (rulings file, a.csv)."""
    rules = copy_of_fixture(folder)
    first = make_result(folder, "a.csv")
    make_result(folder, "b.csv")
    rulings.stamp(first, ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules)
    edit(rules, "value  = 1.4", "value  = 1.2")
    return rules, first


# ----------------------------------------------------------------------------------------------- 1. the brief's six
def test_load_flattens_keys():
    """[rotor.bridge_outer] is the key "rotor.bridge_outer": the dots of a table's name are the nesting. The fixture has two quantities
    under [rotor...] and one under [stator...], so the flattening has to walk into a group that holds more than one."""
    found = rulings.load(FIX)
    assert set(found) == {"rotor.bridge_outer", "stator.slots", "rotor.steel"}
    assert found["rotor.bridge_outer"].value == 1.4
    assert found["rotor.bridge_outer"].status == "assumed"


def test_bad_status_is_rejected(tmp_path):
    bad = copy_of_fixture(tmp_path)
    edit(bad, 'status = "assumed"', 'status = "maybe"')
    with pytest.raises(ValueError, match="status") as caught:
        rulings.load(bad)
    assert "rotor.bridge_outer" in str(caught.value) and "maybe" in str(caught.value)      # which ruling, and what it said


def test_changed_ruling_marks_result_stale(tmp_path):
    rules = copy_of_fixture(tmp_path)
    result = make_result(tmp_path, "a.csv")
    rulings.stamp(result, ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules)
    assert rulings.stale(tmp_path, rulings_path=rules) == []            # nothing has changed yet: the result is fresh
    edit(rules, "value  = 1.4", "value  = 1.2")
    found = rulings.stale(tmp_path, rulings_path=rules)
    assert len(found) == 1
    assert found[0].result == result
    assert found[0].reason == "changed"
    assert found[0].changed == {"rotor.bridge_outer": (1.4, 1.2)}
    assert found[0].step == "EM-1"                                       # the step to re-run is the one the stamp names


def test_unstamped_result_is_stale_not_fresh(tmp_path):
    """A result that no stamp vouches for is listed: nothing says which rulings it used, so nothing can say they still hold."""
    rules = copy_of_fixture(tmp_path)
    result = make_result(tmp_path, "b.csv")
    found = rulings.stale(tmp_path, rulings_path=rules)
    assert [(s.result, s.reason, s.step, s.changed) for s in found] == [(result, "unstamped", None, {})]


def test_removed_key_is_reported(tmp_path):
    """A stamped key that the rulings file no longer has: reason "missing-key", and the key is listed with None for its new value.
    The other stamped key did not change, so it is not listed."""
    rules = copy_of_fixture(tmp_path)
    result = make_result(tmp_path)
    rulings.stamp(result, ["rotor.bridge_outer", "stator.slots"], step="EM-1", script="study.py", rulings_path=rules)
    drop_table(rules, "rotor.bridge_outer")
    found = rulings.stale(tmp_path, rulings_path=rules)
    assert len(found) == 1 and found[0].reason == "missing-key"
    assert found[0].changed == {"rotor.bridge_outer": (1.4, None)}


def test_stamp_records_the_commit(tmp_path):
    """R6: the stamp's field is `commit` (not the brief's proc_commit), and it is the Storage repository's HEAD, `git rev-parse HEAD` in P."""
    head = head_of_the_code()
    written = rulings.stamp(make_result(tmp_path), ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=FIX)
    recorded = json.loads(written.read_text(encoding="utf-8"))
    assert recorded["commit"] == head
    assert "proc_commit" not in recorded


# ----------------------------------------------------------------------------------------------- 2. load(), value(), where the file is
def test_load_reads_every_field_of_a_ruling_with_the_type_toml_gave_the_value():
    found = rulings.load(FIX)
    assert found["rotor.bridge_outer"] == rulings.Ruling(
        key="rotor.bridge_outer", value=1.4, unit="mm", status="assumed", source="OnShape Bridge_Thickness; motorgeo Params.bridge",
        step="OP-1", names={"onshape": "Bridge_Thickness", "motorgeo": "bridge"})
    assert found["stator.slots"].value == 36 and isinstance(found["stator.slots"].value, int) and found["stator.slots"].status == "ruled"
    assert found["rotor.steel"].value == "M-19" and found["rotor.steel"].names == {}


def test_a_ruling_cannot_be_changed_after_it_was_loaded():
    ruling = rulings.load(FIX)["rotor.bridge_outer"]
    with pytest.raises(dataclasses.FrozenInstanceError):
        ruling.value = 2.0


@pytest.mark.parametrize("field", ["value", "unit", "status", "source", "step", "names"])
def test_a_missing_field_is_rejected_and_named(tmp_path, field):
    bad = copy_of_fixture(tmp_path)
    drop_line(bad, field)
    with pytest.raises(ValueError, match=field) as caught:
        rulings.load(bad)
    assert "rotor.bridge_outer" in str(caught.value)


GOOD = 'value = 1\nunit = "mm"\nstatus = "ruled"\nsource = "s"\nstep = "T"\nnames = {}\n'      # the body of one valid ruling
BAD_FILES = {
    "an extra field is not dropped silently": (f"[a.b]\n{GOOD}note = 'x'\n", "note"),
    "a stray number beside the tables": (f"version = 1\n[a.b]\n{GOOD}", "version"),
    "the same key from two tables": (f"[a.b]\n{GOOD}[\"a.b\"]\n{GOOD}", "a.b"),
    "a list is not a value": (f"[a.b]\n{GOOD.replace('value = 1', 'value = [1, 2]')}", "value"),
    "a date is not a value": (f"[a.b]\n{GOOD.replace('value = 1', 'value = 2026-11-06')}", "value"),
    "nan is not a value": (f"[a.b]\n{GOOD.replace('value = 1', 'value = nan')}", "nan"),
    "a unit that is not text": (f"[a.b]\n{GOOD.replace('unit = \"mm\"', 'unit = 3')}", "unit"),
    "names that are not text": (f"[a.b]\n{GOOD.replace('names = {}', 'names = { onshape = 3 }')}", "names"),
    "a file that is not TOML": ("[a.b\nvalue =", "valid TOML"),
    # a table is a group of rulings only if it holds nothing but tables; every other table is a ruling, and is checked as one
    "an empty table is not a group": ("[a.b]\n", r"'a\.b' lacks value, unit, status, source, step, names"),
    "a table with only names": ("[a.b]\nnames = { onshape = 'x' }\n", r"'a\.b' lacks value, unit, status, source, step\."),
    "a table with only keys a ruling does not have": ("[a.b]\nfoo = 1\n", r"'a\.b' lacks value, unit, status, source, step, names and has 'foo'"),
    "a group that also holds a plain line": (f"[a]\nx = 1\n[a.b]\n{GOOD}", r"'a' lacks .* and has 'x', 'b'"),
}


@pytest.mark.parametrize("text, said", list(BAD_FILES.values()), ids=list(BAD_FILES))
def test_a_file_that_is_not_in_the_format_is_a_value_error_that_says_what_and_where(tmp_path, text, said):
    """The format is checked as the file is read, so that a typo is an error naming the ruling, and never a ruling quietly missing or
    quietly changed. The good body is accepted: only what each case changes is wrong."""
    good = tmp_path / "good.toml"
    good.write_text(f"[a.b]\n{GOOD}", encoding="utf-8")
    assert rulings.load(good)["a.b"].value == 1
    bad = tmp_path / "bad.toml"
    bad.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match=said) as caught:
        rulings.load(bad)
    assert "bad.toml" in str(caught.value)


def test_a_file_with_a_byte_order_mark_is_refused_and_the_message_says_what_to_do(tmp_path):
    """Some editors save UTF-8 with a BOM, the three bytes EF BB BF in front. tomllib would say only "Invalid statement (at line 1, column 1)";
    the message names the file, says what the mark is, and says what the file has to be: UTF-8 without a BOM."""
    bom = tmp_path / "bom.toml"
    bom.write_bytes(b"\xef\xbb\xbf" + f"[a.b]\n{GOOD}".encode("utf-8"))
    with pytest.raises(ValueError, match="UTF-8 without a BOM") as caught:
        rulings.load(bom)
    assert "bom.toml" in str(caught.value) and "byte-order mark" in str(caught.value)


@pytest.mark.parametrize("encoding", ["latin-1", "utf-16"])
def test_a_file_in_another_encoding_is_refused_and_the_message_names_the_file(tmp_path, encoding):
    """UnicodeDecodeError does not name the file, and a learner who gets it does not know which file to open again. The same text in
    UTF-8 loads, a non-ASCII letter included (the micro sign is one byte in latin-1 and two in UTF-8)."""
    text = "[a.b]\n" + GOOD.replace('unit = "mm"', 'unit = "µm"')
    utf8 = tmp_path / "utf8.toml"
    utf8.write_bytes(text.encode("utf-8"))
    assert rulings.load(utf8)["a.b"].unit == "µm"
    other = tmp_path / f"{encoding}.toml"
    other.write_bytes(text.encode(encoding))
    with pytest.raises(ValueError, match="UTF-8 without a BOM") as caught:
        rulings.load(other)
    assert f"{encoding}.toml" in str(caught.value)


def test_value_gives_a_rulings_value_and_a_key_it_does_not_have_is_a_key_error_that_suggests_the_nearest():
    assert rulings.value("rotor.bridge_outer", FIX) == 1.4
    assert rulings.value("stator.slots", FIX) == 36
    with pytest.raises(KeyError, match="rotor.bridge_outr") as caught:      # a typo: the message names it, and what could have been meant
        rulings.value("rotor.bridge_outr", FIX)
    assert "rotor.bridge_outer" in str(caught.value)


def test_the_default_rulings_file_is_the_one_config_toml_names(tmp_path, monkeypatch):
    """load(None) and value(key) read [paths] rulings of config.toml, like every other path. In the real config it is rulings.toml next to
    paths.py (O-1 writes it); the test points the config at the fixture instead."""
    assert paths.path("rulings") == P / "rulings.toml"
    monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, f"[paths]\nrulings = '{FIX}'\n"))
    assert rulings.load()["rotor.bridge_outer"].value == 1.4
    assert rulings.value("stator.slots") == 36


def test_a_missing_rulings_file_names_the_path_and_the_session_that_writes_it(tmp_path, monkeypatch):
    nothing = tmp_path / "nowhere" / "rulings.toml"
    monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, f"[paths]\nrulings = '{nothing}'\n"))
    for call in (lambda: rulings.load(), lambda: rulings.load(nothing), lambda: rulings.value("rotor.bridge_outer"),
                 lambda: rulings.stamp(make_result(tmp_path), ["rotor.bridge_outer"], step="EM-1", script="study.py"),
                 lambda: rulings.stale(tmp_path)):
        with pytest.raises(FileNotFoundError) as caught:
            call()
        assert str(nothing) in str(caught.value) and "O-1" in str(caught.value)


# ----------------------------------------------------------------------------------------------- 3. stamp()
def test_a_stamp_has_the_documented_fields_and_sits_next_to_its_result(tmp_path):
    head_of_the_code()      # skips this test when P is not in a git checkout: the stamp then has no commit and no dirty flag, by design
    result = make_result(tmp_path, "a.csv")
    written = rulings.stamp(result, ["rotor.bridge_outer", "stator.slots", "rotor.steel"], step="EM-1", script="study.py", rulings_path=FIX)
    assert written == tmp_path / "a.csv.stamp.json"                      # <result>.stamp.json: the whole name of the result, then the suffix
    assert b"\r" not in written.read_bytes()                             # the same bytes on every machine: lines end in \n here too, not \r\n
    recorded = json.loads(written.read_text(encoding="utf-8"))
    assert set(recorded) == {"step", "script", "commit", "dirty", "written", "rulings"}
    assert recorded["step"] == "EM-1" and recorded["script"] == "study.py"
    assert recorded["rulings"] == {"rotor.bridge_outer": 1.4, "stator.slots": 36, "rotor.steel": "M-19"}
    assert isinstance(recorded["dirty"], bool)
    assert datetime.fromisoformat(recorded["written"]).tzinfo is not None      # the time of the stamp, with its offset


def test_stamping_again_after_the_rerun_replaces_the_stamp_and_makes_the_result_fresh(tmp_path):
    """The loop the stamps are for: a ruling changes, the result is listed, the step is run again, its result is stamped again, and the list is empty."""
    rules, result = stamped_then_changed(tmp_path)
    assert [s.result.name for s in rulings.stale(tmp_path, rulings_path=rules)] == ["a.csv", "b.csv"]
    again = rulings.stamp(result, ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules)
    assert json.loads(again.read_text(encoding="utf-8"))["rulings"] == {"rotor.bridge_outer": 1.2}
    assert [s.result.name for s in rulings.stale(tmp_path, rulings_path=rules)] == ["b.csv"]         # b.csv was never stamped: still listed


def test_stamp_refuses_one_string_for_the_keys(tmp_path):
    """A string is a sequence of one-letter strings in Python: stamp(p, "rotor.bridge_outer") would look for the keys "r", "o", "t"..."""
    result = make_result(tmp_path)
    with pytest.raises(TypeError, match="list"):
        rulings.stamp(result, "rotor.bridge_outer", step="EM-1", script="study.py", rulings_path=FIX)
    assert not (tmp_path / "a.csv.stamp.json").exists()


def test_stamp_refuses_a_result_that_is_not_there(tmp_path):
    with pytest.raises(FileNotFoundError, match="missing.csv"):
        rulings.stamp(tmp_path / "missing.csv", ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=FIX)
    assert not (tmp_path / "missing.csv.stamp.json").exists()


def test_stamp_names_a_key_the_rulings_do_not_have_and_writes_nothing(tmp_path):
    result = make_result(tmp_path)
    with pytest.raises(KeyError, match="rotor.nothing"):
        rulings.stamp(result, ["rotor.bridge_outer", "rotor.nothing"], step="EM-1", script="study.py", rulings_path=FIX)
    assert not (tmp_path / "a.csv.stamp.json").exists()                  # a stamp that lists only some of the keys would be a false promise


def test_dirty_describes_the_code_folder_and_not_the_rest_of_the_repository(tmp_path):
    """`dirty` answers "was this code that commit?", and the code is all in the folder asked about (P). Git lists what differs under it:
    a tracked file that was changed, or an untracked file that is not ignored. A change elsewhere in the repository (here `other`, and the
    top of the repository) does not count, and an ignored file does not either."""
    require_git()
    head = make_repo(tmp_path)
    code, other = tmp_path / "code", tmp_path / "other"
    assert rulings._git_state(code) == (head, False, None)
    (other / "b.txt").write_text("changed", encoding="utf-8")            # tracked and changed, but outside the folder
    (other / "new.txt").write_text("x", encoding="utf-8")                # untracked, outside
    (tmp_path / "top.txt").write_text("x", encoding="utf-8")             # untracked, at the top of the repository
    assert rulings._git_state(code) == (head, False, None)               # edits elsewhere are not this code's business
    (code / "x.log").write_text("x", encoding="utf-8")                   # untracked, inside, but ignored
    assert rulings._git_state(code) == (head, False, None)
    (code / "new.txt").write_text("x", encoding="utf-8")                 # untracked, inside
    assert rulings._git_state(code) == (head, True, None)
    (code / "new.txt").unlink()
    (code / "sub").mkdir()
    (code / "sub" / "deep.txt").write_text("x", encoding="utf-8")        # untracked, inside, in a new sub-folder
    assert rulings._git_state(code) == (head, True, None)
    (code / "sub" / "deep.txt").unlink()
    (code / "sub").rmdir()
    (code / "a.txt").write_text("changed", encoding="utf-8")             # tracked, inside, changed
    assert rulings._git_state(code) == (head, True, None)
    (code / "a.txt").write_text("a", encoding="utf-8")                   # put back: clean again
    assert rulings._git_state(code) == (head, False, None)


def test_untracked_files_count_whatever_the_users_git_config_says(tmp_path):
    """status.showUntrackedFiles = no is a setting people have, and plain `git status` would then never list a new file. The stamp asks for
    untracked files explicitly (--untracked-files=normal), so a new file in the code folder is always `dirty`."""
    require_git()
    head = make_repo(tmp_path)
    subprocess.run(["git", "config", "status.showUntrackedFiles", "no"], cwd=tmp_path, capture_output=True, text=True, check=True)
    (tmp_path / "code" / "new.txt").write_text("x", encoding="utf-8")
    assert rulings._git_state(tmp_path / "code") == (head, True, None)


def test_a_stamp_asks_git_about_the_code_folder(tmp_path, monkeypatch):
    """The folder git is asked about is P, whatever folder the study was started from: the stamp says whether THIS code was that commit."""
    asked = []
    monkeypatch.setattr(rulings, "_git_state", lambda folder: asked.append(folder) or ("abc123", False, None))
    monkeypatch.chdir(tmp_path)
    written = rulings.stamp(make_result(tmp_path), ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=FIX)
    assert asked == [paths.P]
    recorded = json.loads(written.read_text(encoding="utf-8"))
    assert recorded["commit"] == "abc123" and recorded["dirty"] is False and "commit_note" not in recorded


def test_a_folder_that_is_not_a_repository_has_no_commit_and_the_stamp_says_why(tmp_path, monkeypatch):
    require_git()
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path.parent))  # git does not look above the test's folder for a repository
    commit, dirty, why = rulings._git_state(tmp_path)
    assert commit is None and dirty is None
    assert "not a git repository" in why


def test_without_git_the_commit_is_none_and_the_stamp_says_why(tmp_path, monkeypatch):
    monkeypatch.setattr(rulings.shutil, "which", lambda name: None)
    written = rulings.stamp(make_result(tmp_path), ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=FIX)
    recorded = json.loads(written.read_text(encoding="utf-8"))
    assert recorded["commit"] is None and recorded["dirty"] is None
    assert "git" in recorded["commit_note"]
    assert recorded["rulings"] == {"rotor.bridge_outer": 1.4}            # the rest of the stamp is complete


@pytest.mark.parametrize("trouble", [OSError("cannot start"), subprocess.TimeoutExpired("git", 30)], ids=["cannot start", "does not answer"])
def test_a_git_that_cannot_start_or_does_not_answer_is_no_reason_to_lose_the_stamp(tmp_path, monkeypatch, trouble):
    require_git()

    def broken(*args, **kwargs):
        raise trouble
    monkeypatch.setattr(rulings.subprocess, "run", broken)
    written = rulings.stamp(make_result(tmp_path), ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=FIX)
    recorded = json.loads(written.read_text(encoding="utf-8"))
    assert recorded["commit"] is None and recorded["dirty"] is None and recorded["commit_note"]
    assert recorded["rulings"] == {"rotor.bridge_outer": 1.4}


# ----------------------------------------------------------------------------------------------- 4. stale()
def test_a_result_whose_rulings_have_not_changed_is_fresh_even_when_others_have(tmp_path):
    """Only the keys a result was stamped with count: the steel ruling changes, and a result that used only the bridge stays fresh."""
    rules = copy_of_fixture(tmp_path)
    rulings.stamp(make_result(tmp_path), ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules)
    edit(rules, 'value  = "M-19"', 'value  = "M-15"')
    edit(rules, "value  = 36", "value  = 24")
    assert rulings.stale(tmp_path, rulings_path=rules) == []


def test_every_kind_of_result_is_looked_at_and_a_stamp_is_not_a_result(tmp_path):
    """csv, npz, mat and json, in sub-folders too, and nothing else. The two stamped files are fresh; their stamps (which are .json files)
    are not listed as unstamped results; notes.txt is not a result."""
    rules = copy_of_fixture(tmp_path)
    for name in ("a.csv", "b.npz", "c.mat", "d.json", "sub/deep/e.csv"):
        make_result(tmp_path, name)
    make_result(tmp_path, "notes.txt")
    for name in ("a.csv", "d.json"):
        rulings.stamp(tmp_path / name, ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules)
    found = rulings.stale(tmp_path, rulings_path=rules)
    assert [s.result.relative_to(tmp_path).as_posix() for s in found] == ["b.npz", "c.mat", "sub/deep/e.csv"]
    assert [p.relative_to(tmp_path).as_posix() for p in rulings.results(tmp_path)] == ["a.csv", "b.npz", "c.mat", "d.json", "sub/deep/e.csv"]


def test_folders_whose_names_start_with_a_dot_are_skipped_but_a_root_that_is_one_is_looked_at(tmp_path):
    """The queue's lock folder, out\\.femm, holds holder.json: a note about who has the slot, not a result (femm_queue._lock_folder is
    paths.path("out") / ".femm"). No study writes into a folder whose name starts with a dot, so none is looked into, at any depth, and the slot's own holder.json can
    never make `stale.py` exit 1 (a run's unstamped results elsewhere under out\\ still do). A dot INSIDE a name does not make it hidden (T5-P6): the folder mesh_0.5 and the file
    a.b.csv are looked at. The root itself may be named anything."""
    rules = copy_of_fixture(tmp_path)
    for name in ("result.csv", "sub/ok.csv", "a.b.csv", "mesh_0.5/m.csv", ".femm/holder.json", ".cache/x.json", "sub/.deep/y.csv"):
        make_result(tmp_path, name)
    kept = ["a.b.csv", "mesh_0.5/m.csv", "result.csv", "sub/ok.csv"]
    assert [p.relative_to(tmp_path).as_posix() for p in rulings.results(tmp_path)] == kept
    assert [s.result.relative_to(tmp_path).as_posix() for s in rulings.stale(tmp_path, rulings_path=rules)] == kept
    assert [p.name for p in rulings.results(tmp_path / ".femm")] == ["holder.json"]         # asked for by name, it is looked at


def test_a_stamp_that_cannot_be_read_counts_as_no_stamp(tmp_path):
    """A half-written stamp, or a file that is not a stamp, vouches for nothing."""
    rules = copy_of_fixture(tmp_path)
    for name, damage in (("a.csv", "{not json"), ("b.csv", "[1, 2]"), ("c.csv", '{"step": "EM-1"}'), ("d.csv", '{"rulings": 3}')):
        result = make_result(tmp_path, name)
        rulings.stamp(result, ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules).write_text(damage, encoding="utf-8")
    assert [(s.result.name, s.reason) for s in rulings.stale(tmp_path, rulings_path=rules)] == [(n, "unstamped") for n in ("a.csv", "b.csv", "c.csv", "d.csv")]


def test_true_is_not_1_and_1_point_0_is_1(tmp_path):
    """Python says True == 1. A ruling that goes from 1 to true is a change all the same; from 1 to 1.0 it is not (the same number)."""
    rules = tmp_path / "rulings.toml"
    rules.write_text(f"[a.b]\n{GOOD}", encoding="utf-8")
    rulings.stamp(make_result(tmp_path), ["a.b"], step="T", script="s.py", rulings_path=rules)
    rules.write_text(f"[a.b]\n{GOOD.replace('value = 1', 'value = 1.0')}", encoding="utf-8")
    assert rulings.stale(tmp_path, rulings_path=rules) == []
    rules.write_text(f"[a.b]\n{GOOD.replace('value = 1', 'value = true')}", encoding="utf-8")
    assert [s.changed for s in rulings.stale(tmp_path, rulings_path=rules)] == [{"a.b": (1, True)}]


def test_a_root_that_is_not_a_folder_is_an_error_and_not_an_empty_list(tmp_path):
    """"Nothing is stale" must never come from looking in the wrong place."""
    with pytest.raises(FileNotFoundError, match="nope"):
        rulings.stale(tmp_path / "nope", rulings_path=FIX)


# ----------------------------------------------------------------------------------------------- 5. stale.py, the command line
def test_stale_py_prints_each_stale_result_its_reason_its_changed_keys_and_its_step_and_exits_1(tmp_path, capsys):
    rules, _ = stamped_then_changed(tmp_path)
    assert stale_cli.main([str(tmp_path), "--rulings", str(rules)]) == 1
    shown = capsys.readouterr().out
    assert "2 of 2 results" in shown
    assert "a.csv: changed, re-run step EM-1" in shown
    assert "rotor.bridge_outer: 1.4 -> 1.2" in shown
    assert "b.csv: unstamped" in shown and "not known" in shown          # no stamp, so no step to name


def test_stale_py_says_which_ruling_is_gone(tmp_path, capsys):
    rules = copy_of_fixture(tmp_path)
    rulings.stamp(make_result(tmp_path), ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules)
    drop_table(rules, "rotor.bridge_outer")
    assert stale_cli.main([str(tmp_path), "--rulings", str(rules)]) == 1
    shown = capsys.readouterr().out
    assert "a.csv: missing-key, re-run step EM-1" in shown
    assert "rotor.bridge_outer: 1.4 -> (not in the rulings file)" in shown


def test_stale_py_exits_0_and_says_how_much_it_checked_when_nothing_is_stale(tmp_path, capsys):
    rules = copy_of_fixture(tmp_path)
    rulings.stamp(make_result(tmp_path, "a.csv"), ["rotor.bridge_outer"], step="EM-1", script="study.py", rulings_path=rules)
    assert stale_cli.main([str(tmp_path), "--rulings", str(rules)]) == 0
    shown = capsys.readouterr().out
    assert "nothing is stale" in shown and "1 result " in shown          # one result, looked at and found fresh


def test_stale_py_with_no_results_says_nothing_was_checked(tmp_path, capsys):
    """Exit status 0 (nothing is stale), but the line says why: there was nothing to look at."""
    assert stale_cli.main([str(tmp_path), "--rulings", str(FIX)]) == 0
    assert "no results" in capsys.readouterr().out


def test_stale_py_without_a_root_looks_under_out_and_without_rulings_reads_the_configured_file(tmp_path, capsys, monkeypatch):
    out = tmp_path / "out"
    make_result(out, "b.csv")
    monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, f"[paths]\nout = '{out}'\nrulings = '{FIX}'\n"))
    assert stale_cli.main([]) == 1
    assert "b.csv: unstamped" in capsys.readouterr().out


def test_stale_py_help_says_that_folders_starting_with_a_dot_are_skipped(capsys):
    with pytest.raises(SystemExit) as caught:
        stale_cli.main(["--help"])
    assert caught.value.code == 0
    shown = " ".join(capsys.readouterr().out.split())                    # argparse wraps its lines: compare the words, not the line breaks
    assert "start with a dot" in shown and ".femm" in shown


def test_stale_py_exits_2_with_a_message_when_it_cannot_look(tmp_path, capsys):
    """Not a traceback, and not the 1 that means "something is stale": a root that is not there, a rulings file that is not there, a rulings file that is wrong."""
    assert stale_cli.main([str(tmp_path / "nope"), "--rulings", str(FIX)]) == 2
    assert "nope" in capsys.readouterr().err
    assert stale_cli.main([str(tmp_path), "--rulings", str(tmp_path / "no rulings.toml")]) == 2
    assert "O-1" in capsys.readouterr().err
    bad = copy_of_fixture(tmp_path)
    edit(bad, 'status = "assumed"', 'status = "maybe"')
    assert stale_cli.main([str(tmp_path), "--rulings", str(bad)]) == 2
    assert "status" in capsys.readouterr().err


def test_a_folder_that_cannot_be_read_is_an_error_and_not_a_silent_gap(tmp_path, monkeypatch, capsys):
    """T5-P3. os.walk skips a folder it cannot read without a word, as rglob did, so a result in it was never looked at and the answer read 'nothing is stale' for a place
    that was not looked at. results() now stops with the error (the walk's onerror re-raises it), like every other thing here that cannot be looked at; stale.py already
    maps an OSError to exit status 2 and says it on the error stream."""
    make_result(tmp_path, "ok.csv")
    make_result(tmp_path, "locked/inside.csv")
    real_scandir = os.scandir

    def scandir_that_refuses_one_folder(path=".", *args, **kwargs):
        if Path(path).name == "locked":
            raise PermissionError(13, "Access is denied", str(path))
        return real_scandir(path, *args, **kwargs)
    monkeypatch.setattr(os, "scandir", scandir_that_refuses_one_folder)
    with pytest.raises(PermissionError, match="locked"):
        rulings.results(tmp_path)
    with pytest.raises(PermissionError, match="locked"):
        rulings.stale(tmp_path, rulings_path=FIX)
    assert stale_cli.main([str(tmp_path), "--rulings", str(FIX)]) == 2
    assert "locked" in capsys.readouterr().err


def test_the_docstrings_say_what_is_true_of_the_exit_status_of_the_old_proc_folder_and_of_the_reserved_word_names():
    """T5-P4, P7 and P8, wording that the code decides. (P4) Only the slot's own holder.json cannot make stale.py exit 1: a run's unstamped results still do. (P7) The brief's
    proc_commit was named after the old proc folder, which was never a repository, and `stale.py out` exits 2 until rulings.toml is written: the 266 files were read with the
    fixture. (P8) `names` is reserved, and it is said where somebody who writes a rulings.toml looks: the README's row, and the fixture's header."""
    flat = lambda text: " ".join(text.split())                                                   # (a docstring's line breaks are not part of what it says)
    doc, rulings_doc = flat(stale_cli.__doc__), flat(rulings.__doc__)
    assert "a run that holds the slot cannot make the exit status 1" not in doc and "the slot's own holder.json can never make the exit status 1" in doc
    assert "exits 2 until the orchestration session (O-1) has written rulings.toml" in doc and "rulings_min.toml" in doc.split("Numbers to expect", 1)[1].split("Self-check", 1)[0]
    assert "went with the old proc repository" not in rulings_doc and "never a repository of its own" in rulings_doc
    readme = (P / "README.md").read_text(encoding="utf-8")
    row = next(line for line in readme.splitlines() if line.startswith("| `rulings.py` + `stale.py`"))
    assert "`names` is reserved" in row and "[group.names]" in row
    header = FIX.read_text(encoding="utf-8").split("\n[rotor.bridge_outer]\n", 1)[0]               # (the comments name the table too: cut at the table's own line)
    assert "`names` is reserved" in header and "[rotor.names]" in header
