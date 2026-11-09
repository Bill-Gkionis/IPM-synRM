"""compare_materials.py: do two runs of the materials scripts give the same numbers, key by key?

The materials scripts (materials/*.py) write every number a note quotes into out/materials/*.json (plus one CSV) and
draw figures. Move the scripts to another computer, another Python or another FEMM, run them again, and this script
says which numbers moved. It reads the JSON files of the old run, finds the same file in the new run, and compares
every value:  the same keys, the same list lengths, and each number to a stated tolerance (default: exactly equal;
the scripts write numbers rounded to 6 significant digits, so equal inputs give equal files). Then the CSV files, line
by line. With --figures it also compares the PNG files pixel by pixel.

Run (from the Python to FEMM folder):
    uv run python validation/compare_materials.py OLD_DIR NEW_DIR [--rtol 0] [--figures]
for example
    uv run python validation/compare_materials.py femm42-out\\materials out\\materials --figures
(femm42-out is FEMM 4.2's run of 2026-10-28: a checked copy of the `out` folder that it left in
D:/dev/scripts/automation/proc/out. Its materials/ folder holds the numbers of the old run.)
With --figures the two folders may hold nothing but pictures, so the figures can be compared on their own. That is how the
report's comparison of FEMM 5.0's figures with the project notes' was made (the notes' folder is only read):
    uv run python validation/compare_materials.py "<notes folder>/docs/figures" out/materials/figures --figures
and it prints, for each of the 16 pictures, "equal" or "N of M pixels differ" (report, section 2.9 and open question 3).

Prints one line per file ("equal (123 values)", or each difference with its size) and then a total. A file that only
the old run has ("MISSING from NEW") and a file that only the new run has ("only in NEW") are differences. A CSV row
with another number of cells is a difference too. Exit status 0 if every compared value is equal within the
tolerance, 1 if anything differs, 2 if there was nothing to compare (a folder that is missing, or has no JSON or CSV
files, or, with --figures, no PNG files either) or the options were wrong.

What to expect when only FEMM changes: fields.json holds numbers measured on solved FEMM fields (the weakest triangle
of the magnets, the share of magnet area below a flux density), and magnetic.json's F3 and F4 sections are built from
them, so those are the values that can move. current.json, thermal.json and structural.json read no FEMM result and
should be equal.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


def _same(a: float, b: float, rtol: float) -> bool:
    """Equal, or both "not a number" (a blank in a table, such as mu_r where H = 0), or within rtol of each other."""
    if a == b or (math.isnan(a) and math.isnan(b)):
        return True
    return math.isclose(a, b, rel_tol=rtol, abs_tol=0.0)


def _walk(old, new, where, rtol, found: list, count: list) -> None:
    """Compare two parsed JSON values. Appends (where, text, relative_delta_or_None) to `found` for each difference."""
    if isinstance(old, dict) and isinstance(new, dict):
        for key in old:
            if key not in new:
                found.append((f"{where}/{key}", "is only in OLD", None))
            else:
                _walk(old[key], new[key], f"{where}/{key}", rtol, found, count)
        for key in new:
            if key not in old:
                found.append((f"{where}/{key}", "is only in NEW", None))
    elif isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            found.append((where, f"has {len(old)} items in OLD and {len(new)} in NEW", None))
        for k, (a, b) in enumerate(zip(old, new)):
            _walk(a, b, f"{where}[{k}]", rtol, found, count)
    elif isinstance(old, bool) or isinstance(new, bool) or isinstance(old, str) or isinstance(new, str) or old is None or new is None:
        count[0] += 1
        if old != new:
            found.append((where, f"is {old!r} in OLD and {new!r} in NEW", None))
    elif isinstance(old, (int, float)) and isinstance(new, (int, float)):
        count[0] += 1
        if not _same(old, new, rtol):
            delta = (new - old) / abs(old) if old else None
            found.append((where, f"is {old!r} in OLD and {new!r} in NEW", delta))
    else:
        count[0] += 1
        found.append((where, f"is {type(old).__name__} in OLD and {type(new).__name__} in NEW", None))


def compare_json(old_path: Path, new_path: Path, rtol: float) -> tuple[list, int]:
    found, count = [], [0]
    _walk(json.loads(old_path.read_text(encoding="utf-8")), json.loads(new_path.read_text(encoding="utf-8")), "", rtol, found, count)
    return found, count[0]


def compare_csv(old_path: Path, new_path: Path, rtol: float) -> tuple[list, int]:
    """A CSV with one header line: the same header, the same number of rows, the same number of cells in every row,
    then every cell (a number to rtol, anything that is not a number as the same text)."""
    a = old_path.read_text(encoding="utf-8").strip().splitlines()
    b = new_path.read_text(encoding="utf-8").strip().splitlines()
    found, count = [], 0
    if not a or not b:
        return ([] if not a and not b else [("file", f"is empty in {'OLD' if not a else 'NEW'}", None)]), 0
    if a[0] != b[0]:
        found.append(("header", f"is {a[0]!r} in OLD and {b[0]!r} in NEW", None))
    if len(a) != len(b):
        found.append(("row count", f"is {len(a) - 1} in OLD and {len(b) - 1} in NEW", None))
    for r, (la, lb) in enumerate(zip(a[1:], b[1:]), start=1):
        cells_a, cells_b = la.split(","), lb.split(",")
        if len(cells_a) != len(cells_b):                    # zipping the cells alone would hide the ones that are missing
            found.append((f"row {r}", f"has {len(cells_a)} columns in OLD and {len(cells_b)} in NEW", None))
        for c, (x, y) in enumerate(zip(cells_a, cells_b)):
            count += 1
            try:
                fx, fy = float(x), float(y)
            except ValueError:                              # a text cell
                if x != y:
                    found.append((f"row {r}, column {c}", f"is {x!r} in OLD and {y!r} in NEW", None))
                continue
            if not _same(fx, fy, rtol):
                found.append((f"row {r}, column {c}", f"is {x} in OLD and {y} in NEW", (fy - fx) / abs(fx) if fx else None))
    return found, count


def compare_png(old_path: Path, new_path: Path) -> tuple[list, int]:
    """Same size and same pixels? Returns (differences, number of pixels compared)."""
    import numpy as np
    from PIL import Image
    a = np.asarray(Image.open(old_path).convert("RGBA"))
    b = np.asarray(Image.open(new_path).convert("RGBA"))
    if a.shape != b.shape:
        return [("size", f"is {a.shape[1]}x{a.shape[0]} in OLD and {b.shape[1]}x{b.shape[0]} in NEW", None)], 0
    n = int((a != b).any(axis=2).sum())
    return ([("pixels", f"{n} of {a.shape[0] * a.shape[1]} pixels differ", None)] if n else []), a.shape[0] * a.shape[1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Compare two runs of the materials scripts, key by key.")
    ap.add_argument("old", help="the folder of the first run, for example the one FEMM 4.2 and the old scripts made")
    ap.add_argument("new", help="the folder of the second run, for example out/materials")
    ap.add_argument("--rtol", type=float, default=0.0, help="relative tolerance for numbers (default 0: exactly equal)")
    ap.add_argument("--figures", action="store_true", help="also compare the PNG files, pixel by pixel")
    args = ap.parse_args(argv)
    if not args.rtol >= 0:
        ap.error("--rtol must be zero or positive")
    old_dir, new_dir = Path(args.old), Path(args.new)
    for label, folder in (("OLD", old_dir), ("NEW", new_dir)):
        if not folder.is_dir():
            print(f"{label} is not a folder: {folder}", file=sys.stderr)
            return 2
    kinds = [("*.json", compare_json), ("*.csv", compare_csv)]
    if args.figures:
        kinds.append(("*.png", lambda a, b, rtol: compare_png(a, b)))
    if all(next(old_dir.glob(pattern), None) is None for pattern, _ in kinds):
        print(f"no {'JSON, CSV or PNG' if args.figures else 'JSON or CSV'} files in OLD ({old_dir}): there is nothing to compare", file=sys.stderr)
        return 2
    total_files = total_values = total_different = 0
    for pattern, compare in kinds:
        for old_path in sorted(old_dir.glob(pattern)):
            new_path = new_dir / old_path.name
            total_files += 1
            if not new_path.is_file():
                print(f"{old_path.name}: MISSING from NEW")
                total_different += 1
                continue
            found, count = compare(old_path, new_path, args.rtol)
            total_values += count
            total_different += len(found)
            if not found:
                print(f"{old_path.name}: equal ({count} values)")
                continue
            print(f"{old_path.name}: {len(found)} differences in {count} values")
            for where, text, delta in found:
                print(f"    {where} {text}" + (f"  ({delta * 100:+.3g} %)" if delta is not None else ""))
        for new_path in sorted(new_dir.glob(pattern)):          # a file that only the new run has is a difference too
            if not (old_dir / new_path.name).is_file():
                print(f"{new_path.name}: only in NEW")
                total_files += 1
                total_different += 1
    print(f"\n{total_files} files, {total_values} values compared, {total_different} differences")
    return 0 if total_different == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
