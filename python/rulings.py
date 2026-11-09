r"""rulings: the rulings file's loader, the stamp that a result carries, and the list of what is stale.

Why it exists
    A design is a pile of numbers that somebody had to decide: the outer bridge is 1.4 mm, the lamination steel is M-19, the machine
    has 36 slots. Some are only ASSUMED until a step of the plan rules them, and then they are RULED. Every result the pipeline
    writes (a CSV of torques, an .npz of fields, a .mat for Simulink) was computed from some of them. When a ruling changes, the
    results that used the old value are out of date, and without a record nobody can tell which they are: a careful person re-runs
    everything and a hurried one re-runs nothing. This module is that record (the idea of orchestration):
        rulings.toml   every quantity once: its value, unit, status, source, and the step that rules it
        stamp()        a result says which rulings it used, and which code made it
        stale()        compares the stamps with the file as it is NOW, and lists the results that are out of date, each with the
                       step that makes it again
    The picture is a tin of paint with a sticker: "mixed from batch 14". When batch 14 is recalled nobody opens every tin: the
    stickers are read. A tin with no sticker is the one not to trust.

The file (rulings.toml)
    TOML, like config.toml (paths.py explains the format), but here every ruling is a TABLE. One table per quantity. The dots in
    the table's name are the nesting, and load() turns them back into one key, "rotor.bridge_outer" here:
        [rotor.bridge_outer]
        value  = 1.4                   a number, a string or true/false (not a list, a date or nan)
        unit   = "mm"                  text; "" for a number that has no unit
        status = "assumed"             "assumed" or "ruled"
        source = "OnShape Bridge_Thickness; motorgeo Params.bridge"      where the number comes from
        step   = "OP-1"                the step of the plan that rules it
        names  = { onshape = "Bridge_Thickness", motorgeo = "bridge" }   what the other tools call it; {} if none
    Every field is required and no other is allowed. A table that holds nothing but other tables only GROUPS rulings ([rotor] above
    [rotor.bridge_outer]); every other table is a ruling and is checked as one, an empty table and a table with a typo for a field
    included (`names` is reserved: a table that has a `names` key is a ruling, whatever else it lacks). A typo is therefore an error
    that names the ruling and the file, never a ruling that is quietly missing or quietly changed. The file is UTF-8 without a
    byte-order mark (BOM); some editors add one, and load() says so. The file is written by the orchestration session (O-1); the tests
    use tests/fixtures/rulings_min.toml. Where it is when nobody says: config.toml, [paths] rulings = 'rulings.toml', next to paths.py
    like every relative path. Until O-1 has written it, load() raises FileNotFoundError, and the message says so.

The calls (each takes the rulings file as an optional `path` or `rulings_path`; None means the configured one)
    load(path=None)             {key: Ruling}. A Ruling is frozen: key, value, unit, status, source, step, names. ValueError, naming
                                the file (and the ruling), for a file that is not UTF-8 without a BOM, a status that is not assumed
                                or ruled, a missing field, or anything in the file that is not a ruling
    value(key, path=None)       one ruling's value: 1.4, 36 or 'M-19'. KeyError, with the nearest keys, for a key the file lacks
    stamp(result, keys, *, step, script, rulings_path=None)
                                write <result>.stamp.json beside the result (below) and return its path
    results(root)               the result files under `root`: *.csv, *.npz, *.mat and *.json, but not *.stamp.json; sorted. Folders
                                whose names start with a dot are not looked into (the queue's out\.femm holds holder.json, a note and
                                not a result; no study writes into such a folder). The root itself may be named anything
    stale(root, *, rulings_path=None)
                                a list of Stale(result, step, reason, changed), one for each result under `root` that is out of date

The stamp (what stamp() writes, here for a result a.csv)
    {"step": "EM-1",                   the step that made the result: the one to run again when a ruling changes
     "script": "barrier_study.py",     the script that made it
     "commit": "1ddc0f5...",           `git rev-parse HEAD` run in P: the Storage repository's HEAD (R6 of the tools plan; the brief's
                                       name for the field, proc_commit, was named after the old proc folder, which was never a repository
                                       of its own)
     "dirty": false,                   true when git lists anything UNDER P: a tracked file that differs from the commit, or a file
                                       that is neither tracked nor ignored (git status --porcelain --untracked-files=normal -- .).
                                       A change elsewhere in the repository does not count: the question is "was this code that
                                       commit?", and the code is all in P
     "written": "2026-11-06T15:20:41+02:00",
     "rulings": {"rotor.bridge_outer": 1.4}}      each key the study named, with the value the file had when stamp() ran
    If git is missing, does not answer, or P is not inside a repository that has a commit, "commit" and "dirty" are null and a further
    field, "commit_note", says why (if git names the commit and then cannot say whether the tree is clean, only "dirty" is null). The
    stamp is still written: what matters most in it is "rulings". Git is a quick tool, not a solver, so asking it takes no solver slot.

What stale() says, for each result that is out of date (a result that is up to date is not listed)
    "unstamped"    no stamp, or a stamp that cannot be read (a half-written file, or a file that is not a stamp). Nothing vouches
                   for it, so it is listed and not trusted: step is None, because nothing says which step made it
    "changed"      a stamped key has another value now. `changed` maps each such key to (the stamped value, the value now)
    "missing-key"  a stamped key is no longer in the rulings file at all. It is in `changed` too, as (the stamped value, None):
                   a ruling's value is never None, so None means "gone". It outranks "changed" when a result has both
    For "changed" and "missing-key", `step` is the step the stamp names. Comparing: 1 and 1.0 are the same number, but 1 and true are
    not the same value (Python says True == 1, so the comparison asks for the type of a flag first).

Traps
    * A stamp covers only the keys it is given. A study that reads rotor.mag_t but stamps only rotor.bridge_outer is not listed when
      mag_t changes. Stamp every ruling the result depends on.
    * stamp() reads the rulings file when it is CALLED, not when the study started. If the file changes during a long run, the stamp
      records the new values for a result that was computed from the old ones. Do not edit the rulings file while a run that will
      stamp is going; the solver slot (femm_queue.py) makes the runs take turns, but nothing guards the file.
    * Stamp after the result is written: stamp() refuses a result that is not there. A study that stops between writing its result
      and stamping it leaves no NEW stamp. On a first run that means no stamp, so the result is listed. On a re-run the OLD stamp
      stays: the result is listed if a stamped ruling has moved since, and reads fresh if none has. Either way it is never called
      fresh unless the rulings its stamp names still hold.
    * The list is only as good as the folder it is given: stale() raises FileNotFoundError for a root that is not a folder, instead
      of returning the empty list that would read as "nothing is stale". It never defaults to the whole of P.

Run (from the Python to FEMM folder)
    uv run pytest tests/test_rulings.py -q        the tests: no solver is started, git is asked a few dozen times
    uv run python stale.py                        the command line around stale(): its docstring says what it prints

Self-check (the answers follow).
  1. A study stamps its result with ["rotor.bridge_outer"] but also reads rotor.steel. The steel ruling changes. Is the result listed?
  2. Why is a result with no stamp listed as stale and not taken as fresh?
  3. A stamp says "commit": null. What does that tell, and was the result's stamp still any use?
  4. A ruling goes from 1.4 to 1.2, and the study is run again. What must the study do for stale() to stop listing the result?
  5. Why does stamp() refuse a key the rulings file does not have, and write nothing, instead of stamping the keys it does have?
Answers.
  1. No. A stamp covers only the keys it was given, and stale() compares only those. The fix is in the study: stamp every key it used.
  2. Fresh would mean "checked and still right", and nothing was checked. The safe answer to "can I trust this result?" when no record
     exists is no, and the list is exactly the list of results to run again or to stamp.
  3. Git could not name the code's commit (not installed, no answer, or not in a repository); "commit_note" says which. The
     "rulings" part is intact, so stale() works as always. Only the question "which code made this?" has no answer.
  4. Write a new stamp, with stamp(), after the new result is written. It replaces the old stamp file and records 1.2; then the
     stamp and the file agree. Running the study without stamping leaves the old stamp, which still says 1.4.
  5. A stamp that lists only some of the keys would promise less than the study meant to, and a changed ruling among the missing
     ones would never show. A typo in a key is an error to fix now, not a silent hole in the record.
"""
from __future__ import annotations

import codecs
import json
import math
import os
import shutil
import subprocess
import tomllib
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from difflib import get_close_matches
from pathlib import Path

import paths                                    # config.toml: [paths] rulings is where the rulings file is when nobody says

__all__ = ["Ruling", "load", "value", "stamp", "Stale", "stale", "results"]

STATUSES = ("assumed", "ruled")                 # what a ruling's status may say
FIELDS = ("value", "unit", "status", "source", "step", "names")     # what a ruling has: all of them, and nothing else
STAMP_SUFFIX = ".stamp.json"                    # <result>.stamp.json sits beside the result
SCANNED = ("*.csv", "*.npz", "*.mat", "*.json")  # the files stale() takes for results (a *.stamp.json is a stamp, not a result)
_GIT_TIMEOUT_S = 30.0                           # git rev-parse and git status answer in a few hundredths of a second; this is the net


# ============================================================================================== the file
@dataclass(frozen=True)
class Ruling:
    """One quantity of the design, as the rulings file has it. Frozen: none of its fields can be assigned again once it is loaded.

    key     the dotted name, "rotor.bridge_outer"
    value   the number, string or flag that is assumed or ruled
    unit    "mm", "A/mm2", ... or "" if there is none
    status  "assumed" (a working value) or "ruled" (a step of the plan decided it)
    source  where the value comes from
    step    the step of the plan that rules it
    names   what other tools call it, e.g. {"onshape": "Bridge_Thickness", "motorgeo": "bridge"}
    """
    key: str
    value: float | int | str | bool
    unit: str
    status: str
    source: str
    step: str
    names: dict[str, str]


def _file(path: Path | None) -> Path:
    """The rulings file to read: `path` if one was given, else the one config.toml's [paths] rulings names."""
    return Path(path) if path is not None else paths.path("rulings")


def _is_group(table: dict) -> bool:
    """A table only GROUPS rulings if it is not empty, holds nothing but tables, and has no `names` key. Every other table is checked as a
    ruling, so that an empty table, a table with a typo for a field and a table with only `names` are errors that say what they lack,
    and never rulings that are quietly missing. (`names` is the one reserved word: it marks a ruling even when nothing else is there.)"""
    return bool(table) and "names" not in table and all(isinstance(item, dict) for item in table.values())


def _ruling(key: str, table: dict, where: Path) -> Ruling:
    """The Ruling that `table` (the lines under [key]) describes. ValueError, naming `where` and the key, if it is not one."""
    missing = [name for name in FIELDS if name not in table]
    unknown = [name for name in table if name not in FIELDS]
    if missing or unknown:
        said = []
        if missing:
            said.append(f"lacks {', '.join(missing)}")
        if unknown:
            said.append(f"has {', '.join(map(repr, unknown))}, which a ruling does not have")
        raise ValueError(f"{where}: ruling {key!r} {' and '.join(said)}. A ruling has {', '.join(FIELDS)} and nothing else (write "
                         "names = {} if no other tool has a name for it); a table that only groups rulings holds nothing but tables")
    number = table["value"]
    if not isinstance(number, (bool, int, float, str)):
        raise ValueError(f"{where}: ruling {key!r}: value = {number!r} must be a number, a string or true/false, not a {type(number).__name__}")
    if isinstance(number, float) and math.isnan(number):
        raise ValueError(f"{where}: ruling {key!r}: value = nan is not a value. nan is never equal to itself, so every result "
                         "stamped with it would be stale for ever")
    for name in ("unit", "source", "step"):
        if not isinstance(table[name], str):
            raise ValueError(f"{where}: ruling {key!r}: {name} = {table[name]!r} must be text")
    if table["status"] not in STATUSES:
        raise ValueError(f"{where}: ruling {key!r}: status = {table['status']!r}, but a ruling's status is 'assumed' or 'ruled'")
    names = table["names"]
    if not isinstance(names, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in names.items()):
        raise ValueError(f"{where}: ruling {key!r}: names = {names!r} must be a table of text, such as {{ onshape = \"Bridge_Thickness\" }}, "
                         "or {} if there is none")
    return Ruling(key=key, value=number, unit=table["unit"], status=table["status"], source=table["source"], step=table["step"],
                  names=dict(names))


def _walk(table: dict, prefix: str, where: Path, found: dict[str, Ruling]) -> None:
    """Add every ruling below `table` to `found`, under its dotted key. `prefix` is the key of `table` itself ("" at the top)."""
    for name, item in table.items():
        key = f"{prefix}.{name}" if prefix else name
        if not isinstance(item, dict):                          # (only at the top of the file: a group holds nothing but tables)
            raise ValueError(f"{where}: {key} = {item!r} is not a ruling. A ruling is a table, [{key}], with {', '.join(FIELDS)}; "
                             "a plain line can only sit inside one")
        if _is_group(item):
            _walk(item, key, where, found)                      # a group of rulings: [rotor] only holds [rotor.bridge_outer] and the like
            continue
        if key in found:
            raise ValueError(f"{where}: two tables make the key {key!r} (a name with a dot in quotes, and the same name as nested tables)")
        found[key] = _ruling(key, item, where)


def load(path: Path | None = None) -> dict[str, Ruling]:
    """Every ruling of the file, by dotted key: load()["rotor.bridge_outer"].value is 1.4. `path` None means config.toml's [paths] rulings.

    FileNotFoundError, naming the file and the session that writes it (O-1), if it is not there. ValueError, naming the file (and the
    ruling), for a file that is not UTF-8 without a BOM, is not TOML, or has a ruling that is not in the format (module docstring): a
    status that is not assumed or ruled, a field that is missing or not allowed, an empty table, a value that is not a number, a
    string or a flag, a key that two tables make."""
    where = _file(path)
    try:
        raw = where.read_bytes()                                # bytes first, so that the encoding can be checked and named
    except FileNotFoundError:
        raise FileNotFoundError(f"the rulings file {where} is not there. The real one is written by the orchestration session (O-1), so "
                                "until then pass another rulings file (the tests use tests/fixtures/rulings_min.toml)") from None
    if raw.startswith(codecs.BOM_UTF8):                         # tomllib would only say "Invalid statement (at line 1, column 1)"
        raise ValueError(f"{where} starts with a byte-order mark (BOM, the bytes EF BB BF that some editors put in front of UTF-8 text). "
                         "A rulings file is UTF-8 without a BOM: save it again that way")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:                         # its own message does not say which file
        raise ValueError(f"{where} is not UTF-8 (the byte at position {error.start} is not valid UTF-8). A rulings file is UTF-8 "
                         "without a BOM: save it again that way") from error
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:                    # its message says "line 3, column 5" but not which file
        raise ValueError(f"{where} is not valid TOML: {error}") from error
    found: dict[str, Ruling] = {}
    _walk(data, "", where, found)
    return found


def _lookup(found: dict[str, Ruling], key: str, where: Path) -> Ruling:
    """found[key]; KeyError, naming the key, the file and the nearest keys, if there is none (a typo is the usual reason)."""
    if key not in found:
        near = get_close_matches(key, found, n=3)
        hint = f"; closest names: {', '.join(map(repr, near))}" if near else f" (it has {len(found)} rulings)"
        raise KeyError(f"{where} has no ruling {key!r}{hint}")
    return found[key]


def value(key: str, path: Path | None = None):
    """The value of one ruling: 1.4, 36 or 'M-19' (the type TOML gave it). KeyError, with the nearest keys, if the file has no such ruling."""
    return _lookup(load(path), key, _file(path)).value


# ============================================================================================== the stamp
def _first_line(text: str) -> str:
    text = text.strip()
    return text.splitlines()[0] if text else "it printed no message"


def _git_state(folder: Path) -> tuple[str | None, bool | None, str | None]:
    """(commit, dirty, why not) for the code in `folder` (P): `git rev-parse HEAD` and `git status --porcelain --untracked-files=normal -- .`
    run there.

    Normally (the commit, whether the code differs from it, None). When git cannot say, the last item is the reason and what git could
    not say is None: git is not installed, it did not start or did not answer in 30 s, `folder` is not inside a repository, or the
    repository has no commit yet. `dirty` is True when the status lists anything UNDER `folder`: a tracked file that differs from the
    commit, or an untracked file that is not ignored. A change elsewhere in the repository does not count: the question the stamp answers
    is "was this code that commit?", and the code is all in `folder`. --untracked-files=normal asks for the untracked files whatever the
    user's own status.showUntrackedFiles says. --no-optional-locks keeps this read-only call from taking the repository's index lock,
    so asking never makes somebody else's `git commit` fail."""
    git = shutil.which("git")
    if git is None:
        return None, None, "git is not installed (there is no git on the PATH), so the commit of the code is not known"
    try:
        head = subprocess.run([git, "rev-parse", "HEAD"], cwd=folder, capture_output=True, text=True, errors="replace",
                              stdin=subprocess.DEVNULL, timeout=_GIT_TIMEOUT_S)
        if head.returncode != 0:
            return None, None, f"git could not name the commit in {folder}: {_first_line(head.stderr)}"
        status = subprocess.run([git, "--no-optional-locks", "status", "--porcelain", "--untracked-files=normal", "--", "."], cwd=folder,
                                capture_output=True, text=True, errors="replace", stdin=subprocess.DEVNULL, timeout=_GIT_TIMEOUT_S)
    except (OSError, subprocess.TimeoutExpired) as error:
        return None, None, f"git did not answer in {folder}: {error}"
    commit = head.stdout.strip()
    if status.returncode != 0:
        return commit, None, f"git named the commit but could not say whether the tree is clean in {folder}: {_first_line(status.stderr)}"
    return commit, bool(status.stdout.strip()), None


def _stamp_path(result: Path) -> Path:
    """Where the stamp of `result` is: the result's whole name and then the suffix, so a.csv has a.csv.stamp.json beside it."""
    return result.with_name(result.name + STAMP_SUFFIX)


def stamp(result: Path, keys: Iterable[str], *, step: str, script: str, rulings_path: Path | None = None) -> Path:
    """Write <result>.stamp.json: which rulings `result` was computed from, which step and script made it, and which commit of the code.

    result        the result file, which must already be there (stamp after writing it, never before)
    keys          the rulings the result depends on, as a list of dotted keys. Their values are read from the rulings file NOW
    step, script  the step that makes this result (stale.py names it when a ruling changes) and the script that does
    rulings_path  the rulings file; None means config.toml's [paths] rulings
    Returns the stamp's path. A stamp that is already there is replaced: re-running a study and stamping again is how a result becomes
    fresh again. TypeError for one string instead of a list of keys; FileNotFoundError if the result or the rulings file is not there;
    KeyError for a key the file lacks (nothing is written then). The fields are in the module docstring."""
    if isinstance(keys, (str, bytes)):
        raise TypeError(f"keys must be a list of ruling keys such as ['rotor.bridge_outer'], not the one string {keys!r} "
                        "(a string would be taken letter by letter)")
    result = Path(result)
    if not result.is_file():
        raise FileNotFoundError(f"cannot stamp {result}: it is not there. Write the result first, and stamp it after")
    where = _file(rulings_path)
    found = load(where)
    used = {key: _lookup(found, key, where).value for key in keys}
    commit, dirty, note = _git_state(paths.P)
    record = {"step": step, "script": script, "commit": commit, "dirty": dirty,
              "written": datetime.now().astimezone().isoformat(timespec="seconds"), "rulings": used}
    if note:
        record["commit_note"] = note
    target = _stamp_path(result)
    target.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8", newline="\n")      # \n on every system, not \r\n on Windows
    return target


# ============================================================================================== what is stale
@dataclass(frozen=True)
class Stale:
    """A result that is out of date, and why.

    result   the result file
    step     the step to run again: the one its stamp names; None if it has no stamp (nothing says which step made it)
    reason   "changed" (a stamped ruling has another value now), "unstamped" (no usable stamp) or "missing-key" (a stamped ruling is gone)
    changed  each stamped key that differs: key -> (the stamped value, the value now). A key that is gone has None for the value now
    """
    result: Path
    step: str | None
    reason: str
    changed: dict[str, tuple]


def _unreadable(error: OSError) -> None:
    """os.walk's onerror: the walk would skip a folder it cannot read in silence, so say it instead (results() lets it out as it is)."""
    raise error


def results(root: Path) -> list[Path]:
    """The result files under `root`, sorted: *.csv, *.npz, *.mat and *.json in the folder and below it, except the stamps (*.stamp.json).

    A folder whose name starts with a dot is not looked into, at any depth: the queue's lock folder out\\.femm is one, and its holder.json
    is a note about who has the solver slot, not a result. No study writes into such a folder. The root itself may be named anything.
    FileNotFoundError if `root` is not a folder: a folder that is not there has no results, but "nothing is stale" must never
    come from looking in the wrong place. For the same reason a folder under `root` that cannot be read (PermissionError, say) is an
    error, and not a gap that the walk steps over: its results would never be looked at, and the answer would read "nothing is stale"
    for a place that was not looked at."""
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"{root} is not a folder, so there are no results to look at there")
    suffixes = tuple(pattern[1:] for pattern in SCANNED)        # "*.csv" -> ".csv": the case is ignored, as Windows ignores it in a name
    found = []
    for folder, subfolders, names in os.walk(root, onerror=_unreadable):
        subfolders[:] = [name for name in subfolders if not name.startswith(".")]      # os.walk goes only where this list says
        for name in names:
            lower = name.lower()
            if lower.endswith(suffixes) and not lower.endswith(STAMP_SUFFIX):
                found.append(Path(folder) / name)
    return sorted(found)


def _read_stamp(result: Path) -> dict | None:
    """The stamp of `result` as a dict, or None if there is none or it cannot be used (not JSON, not an object, no "rulings" object)."""
    try:
        recorded = json.loads(_stamp_path(result).read_text(encoding="utf-8"))
    except (OSError, ValueError):                               # no file, a file cut short, or not UTF-8: it vouches for nothing
        return None
    return recorded if isinstance(recorded, dict) and isinstance(recorded.get("rulings"), dict) else None


def _same(stamped, now) -> bool:
    """Is the stamped value the ruling's value now? 1 and 1.0 are the same number; 1 and true are not the same value, although
    Python says True == 1."""
    if isinstance(stamped, bool) != isinstance(now, bool):
        return False
    return stamped == now


def stale(root: Path, *, rulings_path: Path | None = None) -> list[Stale]:
    """The results under `root` that are out of date with the rulings, sorted by path. A result that is up to date is not in the list.

    For each result (results(root), which does not look into folders whose names start with a dot, such as the queue's out\\.femm): no
    usable stamp -> "unstamped"; otherwise each stamped key is looked up in the rulings file as it is now, and a key with another value
    makes the result "changed", a key that is gone "missing-key" (module docstring). The rulings file is read once. FileNotFoundError
    if it is not there (naming O-1) or if `root` is not a folder; ValueError if it is not in the format."""
    current = load(rulings_path)
    listed = []
    for result in results(root):
        recorded = _read_stamp(result)
        if recorded is None:
            listed.append(Stale(result, None, "unstamped", {}))
            continue
        changed = {}
        for key, stamped in recorded["rulings"].items():
            if key not in current:
                changed[key] = (stamped, None)
            elif not _same(stamped, current[key].value):
                changed[key] = (stamped, current[key].value)
        if changed:
            step = recorded.get("step")
            gone = any(key not in current for key in changed)
            listed.append(Stale(result, step if isinstance(step, str) else None, "missing-key" if gone else "changed", changed))
    return listed
