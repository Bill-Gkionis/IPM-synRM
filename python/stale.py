r"""stale.py: list the results that are out of date with the rulings, and the step that makes each one again.

Why it exists
    rulings.py lets a result carry a stamp: which rulings it was computed from. When a different value is ruled, this is the
    question that follows: "which of my results did that change?" This script asks it of a folder of results and prints the answer,
    so that the step to run again is read off a list and not remembered. It starts no solver and changes nothing; it only reads.

Run (from the Python to FEMM folder)
    uv run python stale.py                          every result under out\ (config.toml, [paths] out), against rulings.toml
    uv run python stale.py out\bench                only the results under that folder
    uv run python stale.py --rulings my_rulings.toml    against another rulings file (the tests use tests\fixtures\rulings_min.toml)
Exit status: 0 nothing is stale, 1 something is, 2 it could not look (the folder is not there, the rulings file is not there or is not
in the format, a name in config.toml is missing). The message for a 2 goes to the error stream and names what is wrong.
Folders whose names start with a dot are not looked into, at any depth: the queue's lock folder out\.femm holds holder.json, a note about
who has the solver slot and not a result, and no study writes into a dot-folder. So the slot's own holder.json can never make the exit
status 1 (the unstamped results that a run writes elsewhere under out\ still do). A folder named explicitly as the one to look through
may be called anything.

What it prints. A header with the counts, then one block for each result that is out of date, paths relative to the folder looked at:
    stale.py: 2 of 3 results under D:\...\out are stale (rulings: D:\...\rulings.toml)
    a.csv: changed, re-run step EM-1
        rotor.bridge_outer: 1.4 -> 1.2
    b.csv: unstamped, no stamp says which rulings it used, so the step to re-run is not known
The reasons are those of rulings.stale() (its docstring and the module docstring of rulings.py have them): "changed", "unstamped",
"missing-key". For a key that is gone the arrow ends in "(not in the rulings file)". Nothing stale gives one line, with the number of
results that were looked at, and a folder without any results says so ("no results ... so nothing was checked"): "nothing is stale"
is never printed for a place where nothing was looked at.
Numbers to expect. Until the studies stamp their results, everything is listed as "unstamped": on 2026-11-06 `stale.py out --rulings
tests\fixtures\rulings_min.toml` (the fixture: the real rulings.toml is not written yet) read the 266 result files under out\ in 0.2 s, found
no stamp among them, and exited 1. Without --rulings it exits 2 until the orchestration session (O-1) has written rulings.toml, and says
the rulings file is not there. The scan takes every *.csv, *.npz, *.mat and *.json under
the folder (a *.stamp.json excepted), logs and tables included, because a script that has not adopted rulings.stamp() yet cannot be
told from one that never will. The list shrinks as the studies adopt it (R8 of the tools plan: the Electromagnetics workstream does).

Self-check (the answers follow).
  1. stale.py prints "no results under <folder>, so nothing was checked" and exits 0. Is the project up to date?
  2. A result is listed as "unstamped". Which step does the list say to run again?
  3. Why does a rulings file that is not there give the exit status 2, and not 1?
Answers.
  1. Not known. 0 means "nothing was listed", and here nothing was looked at: a wrong folder looks exactly like a clean one, which is
     why the line says what happened instead of "nothing is stale".
  2. None. Nothing says which step made it: only a stamp names the step. Find the script that wrote the file, run it, and have it
     stamp the result; from then on a changed ruling will name the step.
  3. 1 means "something is stale", and a reader of that status goes to re-run steps. A rulings file that is missing says nothing about
     any result: nothing was compared. 2 is "could not look", with the reason printed.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import paths
import rulings


def _noun(count: int) -> str:
    return "result" if count == 1 else "results"


def _describe(item: rulings.Stale, root: Path) -> list[str]:
    """The lines that tell one stale result's story: its name, why, the step to run again, and each key that moved."""
    name = item.result.relative_to(root)
    if item.reason == "unstamped":
        lines = [f"{name}: unstamped, no stamp says which rulings it used, so the step to re-run is not known"]
    else:
        again = f"re-run step {item.step}" if item.step else "its stamp names no step to re-run"
        lines = [f"{name}: {item.reason}, {again}"]
    for key, (old, new) in item.changed.items():
        lines.append(f"    {key}: {old!r} -> {'(not in the rulings file)' if new is None else repr(new)}")
    return lines


def main(argv: list[str] | None = None) -> int:
    """`python stale.py [folder] [--rulings PATH]`: print what is stale. Returns 0 (nothing is), 1 (something is) or 2 (could not look)."""
    parser = argparse.ArgumentParser(prog="stale.py", description="List the results that are out of date with the rulings.",
                                     epilog="Exit status: 0 nothing is stale, 1 something is, 2 it could not look.")
    parser.add_argument("root", nargs="?", type=Path, default=None,
                        help="the folder to look through, sub-folders too, except those whose names start with a dot, such as the queue's "
                             "lock folder .femm (default: [paths] out of config.toml)")
    parser.add_argument("--rulings", type=Path, default=None, metavar="PATH",
                        help="the rulings file (default: [paths] rulings of config.toml)")
    args = parser.parse_args(argv)
    try:
        root = args.root if args.root is not None else paths.path("out")
        where = args.rulings if args.rulings is not None else paths.path("rulings")
        listed = rulings.stale(root, rulings_path=where)
        looked_at = len(rulings.results(root))
    except KeyError as error:                                   # config.toml lacks a name: paths.py's message is the first argument
        print(f"stale.py: {error.args[0]}", file=sys.stderr)
        return 2
    except (OSError, ValueError) as error:                      # a folder or a rulings file that is not there, a file not in the format
        print(f"stale.py: {error}", file=sys.stderr)
        return 2
    if not looked_at:
        print(f"stale.py: no results under {root} (it looks for {', '.join(rulings.SCANNED)}), so nothing was checked")
        return 0
    if not listed:
        print(f"stale.py: nothing is stale ({looked_at} {_noun(looked_at)} under {root} checked against {where})")
        return 0
    verb = "is" if len(listed) == 1 else "are"
    print(f"stale.py: {len(listed)} of {looked_at} {_noun(looked_at)} under {root} {verb} stale (rulings: {where})")
    for item in listed:
        print("\n".join(_describe(item, Path(root))))
    return 1


if __name__ == "__main__":
    sys.exit(main())
