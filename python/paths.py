r"""paths.py: every tool path and every machine folder comes from config.toml, and from nowhere else.

Why this exists
    A script that holds a line such as   MATLIB = r"D:\Apps\femm-qt\bin\matlib.dat"   works on one PC on one day. Move
    the program, install it again, or hand the project to somebody else, and every script that typed the place has to be
    found and edited. The one that is forgotten fails or, worse, goes on reading an old copy. So the places are written
    once, in config.toml (next to this file), and a script asks for one by its name. The rule, checked by
    tests/test_paths.py (R5 of the tools plan): no script holds a drive path of its own.

The file, config.toml
    TOML is a small text format for settings: [sections], `key = value` lines and # for comments. Python reads it with
    its standard library (tomllib), so nothing has to be installed. One thing to watch: write a Windows path in SINGLE
    quotes, 'D:\Apps\x'. Single quotes make a literal string, in which a backslash is just a backslash. In "double
    quotes" the same text would start escape sequences, and the file would not load.
        [tools]   the programs: femm (FEMM 5.0's folder), ccx, cgx, gmsh, prepomax, elmer, elmergrid, matlab
        [paths]   folders and files: out, rulings, femm42_out, reference_fem, simulink_plant, motor_structure_mat, dxf_rotor,
                  dxf_rotor_with_band, dxf_stator, femm_matlib, vault_figures
        [femm]    settings of the FEMM runs: workers, priority

The three calls (each takes an optional `config`, another file to read instead; the tests use it to try wrong ones)
    tool(name)              the [tools] entry as a Path. FileNotFoundError, naming the tool, if the file has no such line,
                            the line is empty, or the file or folder it names is not there: one kind of error for
                            every way a program can be missing, so a caller that lists the tools needs one `except`.
    path(name)              the [paths] entry as a Path. A relative one starts at P, the folder of this file, so 'out'
                            is P\out whichever folder the script was started from. An absolute one is returned as
                            written. The path need not exist yet (a script makes out\ the first time it writes there).
    setting(section, key)   the value as TOML read it: 8 is an int, 'below_normal' a str.
    A [paths] name or a setting the file does not have at all (a typo) is a KeyError that lists the names it does have.
    A file that is not valid TOML is a ValueError that names the file and the line.

Who else decides
    femm_session.femm_dir() takes [tools] femm first, and FEMM_DIR (environment, then registry) only if the config names
    no FEMM. MOTOR_FIGURES_DIR, if it is set, is still used instead of vault_figures (materials/common.save()).

Run (from the Python to FEMM folder):
    uv run python paths.py [other_config.toml]
Lists every name in the file with the path it gives, marks a tool it cannot find MISSING, and ends with
"N of N tools found". Exit status 0 if every tool is there, 1 if one is not. It starts no program.
Expect with the config.toml of 2026-11-05: "8 of 8 tools found", exit status 0.

Self-check (the answers follow).
  1. config.toml says  out = 'out'.  A script is started from P\validation.  Which folder does path("out") give?
  2. Why does the file write  'D:\Apps\femm-qt\bin'  and not  "D:\Apps\femm-qt\bin" ?
  3. A program is moved, and its line in config.toml is not changed. What does tool() do, and why is that better than
     handing back the path anyway?
Answers.
  1. P\out, the same as from any folder: a relative path starts at P, the folder that holds paths.py, not at the folder
     the script was started from.
  2. In double quotes TOML reads a backslash as the start of an escape sequence (\A is not one, so the file would not
     load). Single quotes make a literal string.
  3. It raises FileNotFoundError with the tool's name, the path it looked at and the line to correct. A path that is not
     there would only fail later, inside the program that was asked to start, with a message about that program.
"""
from __future__ import annotations

import sys
import tomllib
from pathlib import Path

P = Path(__file__).resolve().parent     # the project folder, "Python to FEMM": a relative path in config.toml starts here
CONFIG = P / "config.toml"              # the file every call reads unless it is handed another


def _load(config: Path | None) -> tuple[dict, Path]:
    """(the parsed file, its path). CONFIG is looked up now, not when this module was loaded, so a test can point it elsewhere."""
    where = CONFIG if config is None else Path(config)
    with open(where, "rb") as handle:               # tomllib wants a file opened in binary mode: it decodes the UTF-8 itself
        try:
            return tomllib.load(handle), where
        except tomllib.TOMLDecodeError as error:    # its message says "line 3, column 5" but not which file
            raise ValueError(f"{where} is not valid TOML: {error}") from error


def _entry(section: str, key: str, config: Path | None):
    """(the value of [section] key, the path of the config it came from). KeyError, naming what the file has, if it has none."""
    data, where = _load(config)
    if section not in data:
        raise KeyError(f"{where} has no [{section}] section (it has: {', '.join(data) or 'nothing'})")
    if key not in data[section]:
        raise KeyError(f"{where} has no {key!r} in [{section}] (it has: {', '.join(data[section]) or 'nothing'})")
    return data[section][key], where


def _resolve(value: str) -> Path:
    """P / value. A relative value starts at P; an absolute one wins over P in a join, so it comes back as written."""
    return P / value


def tool(name: str, config: Path | None = None) -> Path:
    """The program (or, for femm, the folder) that [tools] `name` names. FileNotFoundError, naming the tool, if it is not there."""
    try:
        value, where = _entry("tools", name, config)
    except KeyError as error:                        # no such line: still "that tool cannot be found", and in the same error type
        raise FileNotFoundError(f"tool {name!r} is not named in the config: {error.args[0]}") from None
    if not value:
        raise FileNotFoundError(f"tool {name!r} is empty in {where}: put the path of the program on the line [tools] {name}")
    found = _resolve(value)
    if not found.exists():
        raise FileNotFoundError(f"tool {name!r} is not at {found}: it was moved or never installed. "
                                f"Correct the line [tools] {name} in {where}")
    return found


def path(name: str, config: Path | None = None) -> Path:
    """The folder or file that [paths] `name` names; a relative entry starts at P. It need not exist."""
    value, where = _entry("paths", name, config)
    if not value:
        raise ValueError(f"[paths] {name} is empty in {where}: an empty path would mean the project folder itself")
    return _resolve(value)


def setting(section: str, key: str, config: Path | None = None):
    """The value of `key` in `[section]`, with the type TOML gave it (8 is an int, 'below_normal' a str)."""
    return _entry(section, key, config)[0]


def main(argv: list[str] | None = None) -> int:
    """`python paths.py [config.toml]`: list every name in the file and what it gives. Returns 1 if a tool is not found."""
    args = sys.argv[1:] if argv is None else argv
    config = Path(args[0]) if args else None
    data, where = _load(config)
    print(f"config {where}")
    found = problems = 0
    print("[tools]")
    for name in data.get("tools", {}):
        try:
            print(f"  {name:10s} {tool(name, config)}")
            found += 1
        except FileNotFoundError as error:
            problems += 1
            print(f"  {name:10s} MISSING: {error}")
    print("[paths]")
    for name in data.get("paths", {}):
        try:
            there = path(name, config)
            print(f"  {name:14s} {there}" + ("" if there.exists() else "   (not there yet)"))
        except ValueError as error:
            problems += 1
            print(f"  {name:14s} ERROR: {error}")
    print("[femm]")
    for key, value in data.get("femm", {}).items():
        print(f"  {key:10s} {value!r}")
    print(f"{found} of {len(data.get('tools', {}))} tools found")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
