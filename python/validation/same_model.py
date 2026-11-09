"""same_model.py: are two FEMM problem files (.fem) the same model?

Why this exists: the motor model was built and solved with FEMM 4.2, and now it is built and solved with FEMM 5.0.
Before comparing any solved number, we want proof that both FEMMs were given, and kept, the same model. FEMM writes
a model to a plain-text file, so "the same model" can be checked without FEMM: read both files into numbers and
names, then compare the numbers and the names. This script does that and prints every difference it finds.

Run (from the Python to FEMM folder):
    uv run python validation/same_model.py OLD.fem NEW.fem [--rtol 1e-9]
Exit status: 0 = the same model, 1 = a difference in the model, 2 = a file could not be read (it is missing, it is not
a FEMM .fem, it stops in the middle of a section) or an option was wrong.

What a .fem file holds (FEMM 4.2 and 5.0 write the same layout, "Format = 4.0"; FEMM 5.0's own description of it is
D:\\Apps\\femm-qt\\docs\\formats\\magnetics-file-formats.pdf):

  [Key] = value lines      the problem definition: frequency, precision, units, planar or axisymmetric, depth, ...
  [PointProps]             then [BdryProps], [BlockProps] and [CircuitProps]: a count, then that many blocks of
                           <Name> = value lines between <BeginX> and <EndX>. These are the boundary conditions, the
                           materials (a material ends with its table of B-H points) and the circuits.
  [NumPoints] = n          n lines "x y point-property group [local-size]": the nodes (corners)
  [NumSegments] = n        "node0 node1 max-side boundary hidden group": straight edges, nodes counted from 0
  [NumArcSegments] = n     "node0 node1 sweep-degrees max-side-degrees boundary hidden group computed-side": circular
                           edges, counter-clockwise from node0 to node1
  [NumHoles] = n           "x y group": regions that are not meshed ("<No Mesh>" in motorgeo)
  [NumBlockLabels] = n     "x y material mesh-size circuit magnet-direction group turns external ["function"]":
                           one label per region
  The material, boundary, circuit and point-property numbers in those lines count from 1 in the order the blocks
  are listed; 0 means none. Here they are turned back into names, so two files that list the same things in a
  different order still compare equal. Everything from a line that says [Solution] on (a .ans file) is cut off:
  it is a result, not the model.

What "the same" means here:
  * Numbers are compared as numbers, not as text: 1e-008 (FEMM 4.2 on Windows) and 1e-08 (FEMM 5.0) are one number.
    Two numbers are the same when they agree to --rtol (default 1e-9) of the larger one. There is no absolute floor:
    0 and 1e-20 are different numbers, and --rtol 0 asks for numbers that are exactly equal.
  * Positions are the one exception, because turning a part leaves noise where a coordinate should be 0 (3e-15 for
    0). Nodes and labels are not compared by their coordinates but matched by them: two are the same point if they lie
    within rtol times the size of the model (its largest coordinate) of each other, and never nearer than 1e-12 file
    units (the neighbour search cannot use a distance of 0: it would find nothing, not even a point itself).
  * Nodes, segments, arcs and labels are matched by where they are, not by their number in the file, so the same
    corners stored in a different order are still the same model.
  * What is not the model is listed apart and never decides the verdict: how a number is written, the [Comment]
    line, the line endings, the order in which things are stored, and the one value in the file that FEMM computes
    itself instead of being given (the last number of an arc line, "mySideLength": FEMM 4.2 stores the side length it
    worked out, FEMM 5.0 stores 0 for "let the mesher choose").
  * What the tool does not understand it says out loud. A section it does not know (a name that starts with Num or ends
    in Props, such as a future [NumFoo] = 3 with three lines under it) is a WARNING: its lines were not compared, and
    the result then says "as far as the parts it understands go". A file with no nodes, no labels and no materials, or
    one that stops in the middle of a section, is refused (exit status 2): it is not called "the same model".

Numbers to expect (2026-11-04). The default model as check_femm.py writes it on FEMM 5.0 (out\\ipm.fem, coarse mesh)
against FEMM 4.2's out\\ipm.fem of 2026-10-29:
    compared: 13 header fields, 7 boundaries, 4 materials, 47 B-H points, 3 circuits, 167 nodes, 136 segments,
              62 arcs, 40 labels (holes included)
    MODEL DIFFERENCES: 0
    NOT THE MODEL: 3 lines. [Precision] is written 1e-008 and 1e-08, [Comment] is "Add comments here." and "", and
    mySideLength differs on 62 of 62 arcs
    RESULT: the same model
The model that check_femm.py --fine writes has one difference from FEMM 4.2's file: [DoSmartMesh] is 0 and 1.

The functions are importable (parse_fem, compare), so another check can read a model, change it on purpose and
compare again: validation/hand_reference_check.py does that.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

PROP_SECTIONS = ("PointProps", "BdryProps", "BlockProps", "CircuitProps")      # a count of <Begin..>..<End..> blocks
ROW_SECTIONS = ("NumPoints", "NumSegments", "NumArcSegments", "NumHoles", "NumBlockLabels")   # a count of text lines
MIN_VALUES = {"NumPoints": 4, "NumSegments": 6, "NumArcSegments": 7, "NumHoles": 3, "NumBlockLabels": 9}   # per line
PROP_WORDS = {"PointProps": ("point property", "point properties"), "BdryProps": ("boundary", "boundaries"),
              "BlockProps": ("material", "materials"), "CircuitProps": ("circuit", "circuits")}   # (one, many)
NOT_THE_MODEL = ("Comment",)       # header lines that are a note to the reader, not part of the problem
DERIVED = {"myside": "mySideLength, the side length FEMM computes for itself (FEMM 5.0 writes 0 = choose it)"}
NO_MESH = "<No Mesh>"
MIN_POSITION_TOL = 1e-12           # file units: two positions nearer than this are one, however small --rtol is

_SECTION = re.compile(r"^\[(\w+)\]\s*=\s*(.*)$")        # [Precision]   =  1e-08
_TAG = re.compile(r"^<(\w+)>\s*=\s*(.*)$")              # <Mu_x> = 4416
_SOLUTION = re.compile(rb"^\[Solution\]", re.MULTILINE)  # a .ans file: the input, then this bare line, then the result


# ----------------------------------------------------------------------------- 1. reading a file
def _value(text: str):
    """One value from a file: a float if it is a number, else the text without its quotation marks."""
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        return text[1:-1]
    try:
        return float(text)
    except ValueError:
        return text


def _read_blocks(lines: list[str], i: int, count: int, section: str):
    """`count` blocks starting at line i. Returns ([(fields, rows)], next line).

    fields: {"Mu_x": 4416.0, "BlockName": "M-19 Steel", ...}; rows: the untagged lines, e.g. the B-H points.
    Raises ValueError if the file ends before the blocks that the count announces.
    """
    blocks = []
    for n in range(count):
        while i < len(lines) and not lines[i].strip().startswith("<Begin"):
            i += 1
        if i >= len(lines):
            raise ValueError(f"the file ends before block {n + 1} of the {count} that [{section}] announces")
        i += 1
        fields, rows = {}, []
        while i < len(lines) and not lines[i].strip().startswith("<End"):
            text = lines[i].strip()
            i += 1
            tagged = _TAG.match(text)
            if tagged:
                fields[tagged.group(1)] = _value(tagged.group(2))
            elif text:
                rows.append([float(v) for v in text.split()])
        if i >= len(lines):
            raise ValueError(f"the file ends inside block {n + 1} of [{section}]")
        blocks.append((fields, rows))
        i += 1
    return blocks, i


def _name_of(fields: dict, fallback: str) -> str:
    """A block's name: its <BlockName>, <BdryName>, <CircuitName> or <PointName> field."""
    return next((str(v) for k, v in fields.items() if k.endswith("Name")), fallback)


def _pick(names: list[str], index: int) -> str:
    """The name behind a 1-based index in a file's list of properties; '' for 0 (none)."""
    return names[index - 1] if 0 < index <= len(names) else ""


def _parse(raw: bytes) -> dict:
    """The model part of a file (everything before [Solution]) as Python data. See parse_fem for the shape."""
    lines = raw.decode("latin-1").splitlines()
    header, props, order, rows, warnings = {}, {}, {}, {}, []
    i = 0
    while i < len(lines):
        section = _SECTION.match(lines[i].strip())
        i += 1
        if not section:
            continue
        key, text = section.group(1), section.group(2).strip()
        if key in PROP_SECTIONS:
            blocks, i = _read_blocks(lines, i, int(text), key)
            props[key], order[key] = {}, []
            for k, (fields, block_rows) in enumerate(blocks):
                name = _name_of(fields, f"#{k + 1}")
                props[key][name] = {"fields": fields, "rows": block_rows}
                order[key].append(name)
        elif key in ROW_SECTIONS:
            n = int(text)
            chunk = lines[i:i + n]              # the text of each line; split below, because a label may end in a quoted string
            if len(chunk) < n:
                raise ValueError(f"[{key}] announces {n} lines and the file ends after {len(chunk)}")
            for number, line in enumerate(chunk, start=1):
                if line.lstrip().startswith("[") or len(line.split()) < MIN_VALUES[key]:
                    raise ValueError(f"line {number} of [{key}] is not a line of at least {MIN_VALUES[key]} "
                                     f"numbers: {line.strip()[:60]!r}")
            rows[key] = chunk
            i += n
        else:
            if key.startswith("Num") or key.endswith("Props"):
                warnings.append(f"[{key}] is a section this tool does not know: whatever lines belong to it were not "
                                f"compared (has FEMM's file format changed?)")
            header[key] = text[1:-1] if len(text) >= 2 and text[0] == text[-1] == '"' else text
    point_names, bdry_names = order.get("PointProps", []), order.get("BdryProps", [])
    block_names, circuit_names = order.get("BlockProps", []), order.get("CircuitProps", [])
    model = {"header": header, "props": props, "order": order, "crlf": b"\r\n" in raw, "warnings": warnings}
    model["nodes"], model["segments"], model["arcs"], model["labels"] = [], [], [], []
    for r in (line.split() for line in rows.get("NumPoints", [])):
        model["nodes"].append({"x": float(r[0]), "y": float(r[1]), "prop": _pick(point_names, int(float(r[2]))),
                               "group": int(float(r[3])), "lfs": float(r[4]) if len(r) > 4 else -1.0})   # -1 = no limit
    for r in (line.split() for line in rows.get("NumSegments", [])):
        model["segments"].append({"n0": int(r[0]), "n1": int(r[1]), "maxside": float(r[2]),
                                  "boundary": _pick(bdry_names, int(float(r[3]))), "hidden": int(float(r[4])),
                                  "group": int(float(r[5]))})
    for r in (line.split() for line in rows.get("NumArcSegments", [])):
        model["arcs"].append({"n0": int(r[0]), "n1": int(r[1]), "sweep": float(r[2]), "maxside": float(r[3]),
                              "boundary": _pick(bdry_names, int(float(r[4]))), "hidden": int(float(r[5])),
                              "group": int(float(r[6])), "myside": float(r[7]) if len(r) > 7 else 0.0})
    for r in (line.split(None, 9) for line in rows.get("NumBlockLabels", [])):
        model["labels"].append({"x": float(r[0]), "y": float(r[1]),
                                "material": _pick(block_names, int(float(r[2]))) or NO_MESH, "meshsize": float(r[3]),
                                "circuit": _pick(circuit_names, int(float(r[4]))), "magdir": _value(r[5]),
                                "group": int(float(r[6])), "turns": int(float(r[7])), "external": int(float(r[8])),
                                "magdirfctn": _value(r[9]) if len(r) > 9 else ""})
    for r in (line.split() for line in rows.get("NumHoles", [])):      # a hole is a label with no material
        model["labels"].append({"x": float(r[0]), "y": float(r[1]), "material": NO_MESH, "group": int(float(r[2]))})
    for item in (*model["segments"], *model["arcs"]):
        if not (0 <= item["n0"] < len(model["nodes"]) and 0 <= item["n1"] < len(model["nodes"])):
            raise ValueError(f"an edge names node {item['n0']} or {item['n1']}, and the file has {len(model['nodes'])} nodes")
    return model


def parse_fem(path) -> dict:
    """Read a .fem into plain Python data, ready to compare (and, for another check, to edit and compare again).

    Returns a dict:
      "header"   {"Precision": "1e-08", "Depth": "20", ...}: the text after each [Key] =
      "props"    {"BlockProps": {"M-19 Steel": {"fields": {...}, "rows": [[b, h], ...]}, ...}, "BdryProps": ..., ...}
      "order"    {"BlockProps": [names in the order the file lists them], ...}
      "nodes"    [{"x", "y", "prop", "group", "lfs"}]                               in file order
      "segments" [{"n0", "n1", "maxside", "boundary", "hidden", "group"}]           n0, n1 are indexes into "nodes"
      "arcs"     [{"n0", "n1", "sweep", "maxside", "boundary", "hidden", "group", "myside"}]
      "labels"   [{"x", "y", "material", "meshsize", "circuit", "magdir", "group", "turns", "external", "magdirfctn"}]
                 holes are labels too: material "<No Mesh>", with only x, y and group
      "crlf"     True if the file uses Windows line endings
      "warnings" what the tool met and did not understand (see the module docstring)
    Everything from a [Solution] line on (a .ans file) is cut off before reading: it is a result, not the model.

    Raises OSError if the file cannot be opened, and ValueError, saying why, if it is not a FEMM .fem this tool can
    read: it has no nodes, no labels and no materials (so it is not a FEMM file, or it is an empty model), or it stops
    in the middle of a section, or a line holds text where numbers belong.
    """
    raw = Path(path).read_bytes()
    cut = _SOLUTION.search(raw)
    try:
        model = _parse(raw[:cut.start()] if cut else raw)
    except (IndexError, KeyError, ValueError) as exc:
        raise ValueError(f"{path} is not a FEMM .fem file that this tool can read ({exc})") from exc
    if not (model["nodes"] or model["labels"] or model["props"].get("BlockProps")):
        raise ValueError(f"{path} has no nodes, no block labels and no materials: it is not a FEMM .fem file, "
                         f"or it holds an empty model")
    return model


# ----------------------------------------------------------------------------- 2. comparing numbers
def _close(a, b, rtol: float) -> bool:
    """Equal if both are text and the same, or both are numbers that agree to rtol of the larger one (no absolute floor)."""
    if isinstance(a, str) or isinstance(b, str):
        return str(a) == str(b)
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=0.0)


def _show(v) -> str:
    return repr(v) if isinstance(v, str) else f"{v:.12g}"


def _compare_fields(where: str, a: dict, b: dict, rtol: float, out: list[str]) -> None:
    """Compare two {field: value} dicts; add one line to `out` for each field that differs or is missing."""
    for key in sorted(set(a) | set(b)):
        if key not in a or key not in b:
            out.append(f"{where}: field {key} is only in {'OLD' if key in a else 'NEW'}")
        elif not _close(a[key], b[key], rtol):
            out.append(f"{where}: {key} is {_show(a[key])} in OLD and {_show(b[key])} in NEW")


# ----------------------------------------------------------------------------- 3. the problem definition and the properties
def _compare_header(old: dict, new: dict, rtol: float, differences: list[str], formatting: list[str]) -> int:
    keys = list(dict.fromkeys([*old["header"], *new["header"]]))
    for key in keys:
        a, b = old["header"].get(key), new["header"].get(key)
        if a is None or b is None:
            differences.append(f"problem definition: [{key}] is only in {'OLD' if b is None else 'NEW'}")
        elif key in NOT_THE_MODEL:
            if a != b:
                formatting.append(f'problem definition: [{key}] is "{a}" in OLD and "{b}" in NEW (a note, not the model)')
        elif not _close(_value(a), _value(b), rtol):
            differences.append(f"problem definition: [{key}] is {a} in OLD and {b} in NEW")
        elif a != b:
            formatting.append(f"problem definition: [{key}] is written {a} in OLD and {b} in NEW (the same number)")
    return len(keys)


def _compare_props(old: dict, new: dict, rtol: float, differences: list[str], formatting: list[str]) -> dict:
    """Boundaries, materials (with their B-H points), circuits and point properties, matched by name."""
    counts = {}
    for section in PROP_SECTIONS:
        a, b = old["props"].get(section, {}), new["props"].get(section, {})
        one, many = PROP_WORDS[section]
        counts[many] = len(a)
        for name in a:
            if name not in b:
                differences.append(f'{one} "{name}" is only in OLD')
        for name in b:
            if name not in a:
                differences.append(f'{one} "{name}" is only in NEW')
        if list(a) != list(b) and set(a) == set(b):
            formatting.append(f"{many}: listed in a different order: OLD {list(a)}, NEW {list(b)}")
        rows_are = "B-H points" if section == "BlockProps" else "data lines"
        for name in a:
            if name not in b:
                continue
            where = f'{one} "{name}"'
            _compare_fields(where, a[name]["fields"], b[name]["fields"], rtol, differences)
            ra, rb = a[name]["rows"], b[name]["rows"]
            if section == "BlockProps":
                counts["B-H points"] = counts.get("B-H points", 0) + len(ra)
            if len(ra) != len(rb):
                differences.append(f"{where}: {len(ra)} {rows_are} in OLD and {len(rb)} in NEW")
            for k, (pa, pb) in enumerate(zip(ra, rb)):
                if len(pa) != len(pb) or not all(_close(u, v, rtol) for u, v in zip(pa, pb)):
                    differences.append(f"{where}: {rows_are[:-1]} {k} is ({', '.join(map(_show, pa))}) in OLD and "
                                       f"({', '.join(map(_show, pb))}) in NEW")
    return counts


# ----------------------------------------------------------------------------- 4. the geometry
def _match(a_xy: np.ndarray, b_xy: np.ndarray, tol: float) -> list[int]:
    """For each point of a_xy, the index of the point of b_xy that lies within tol of it, else -1."""
    if len(a_xy) == 0 or len(b_xy) == 0:
        return [-1] * len(a_xy)
    distance, index = cKDTree(b_xy).query(a_xy, distance_upper_bound=tol)
    return [int(j) if np.isfinite(d) else -1 for d, j in zip(distance, index)]


def _xy(items: list[dict]) -> np.ndarray:
    return np.array([[it["x"], it["y"]] for it in items], dtype=float).reshape(-1, 2)


def _at(item: dict) -> str:
    return f"({item['x']:.9g}, {item['y']:.9g})"


def _compare_points(kind: str, old_items: list[dict], new_items: list[dict], tol: float, rtol: float, keys: tuple,
                    differences: list[str], formatting: list[str]) -> list[int]:
    """Nodes or labels: match by position, then compare the properties in `keys`. Returns old index -> new index."""
    to_new = _match(_xy(old_items), _xy(new_items), tol)
    for k, j in enumerate(to_new):
        if j < 0:
            differences.append(f"{kind} {_at(old_items[k])} is only in OLD")
    for j in sorted(set(range(len(new_items))) - set(to_new)):
        differences.append(f"{kind} {_at(new_items[j])} is only in NEW")
    for j, n in Counter(j for j in to_new if j >= 0).items():
        if n > 1:
            differences.append(f"{n} {kind}s of OLD lie on the one in NEW at {_at(new_items[j])}")
    for k, j in enumerate(to_new):
        if j >= 0:
            a = {c: old_items[k][c] for c in keys if c in old_items[k]}
            b = {c: new_items[j][c] for c in keys if c in new_items[j]}
            _compare_fields(f"{kind} {_at(old_items[k])}", a, b, rtol, differences)
    moved = sum(j != k for k, j in enumerate(to_new) if j >= 0)
    if moved:
        formatting.append(f"{kind}s: stored in a different order ({moved} of {len(old_items)} have another number)")
    return to_new


def _compare_edges(kind: str, old: dict, new: dict, to_new: list[int], rtol: float, keys: tuple,
                   differences: list[str], formatting: list[str]) -> None:
    """Segments or arcs: an edge is its two end nodes (found in NEW through the node matching) plus its properties."""
    old_edges, new_edges = old[kind + "s"], new[kind + "s"]
    ordered = kind == "arc"                     # an arc runs counter-clockwise from n0 to n1; a segment has no direction

    def key_of(n0, n1):
        return (n0, n1) if ordered else frozenset((n0, n1))

    def where(model, e):
        a, b = model["nodes"][e["n0"]], model["nodes"][e["n1"]]
        return f"{kind} {_at(a)} -> {_at(b)}" if ordered else f"{kind} {_at(a)} - {_at(b)}"

    in_new = {}
    for j, e in enumerate(new_edges):
        in_new.setdefault(key_of(e["n0"], e["n1"]), []).append(j)
    used, moved, computed = set(), 0, {name: [] for name in DERIVED}
    for k, e in enumerate(old_edges):
        n0, n1 = to_new[e["n0"]], to_new[e["n1"]]
        candidates = [j for j in in_new.get(key_of(n0, n1), []) if j not in used] if n0 >= 0 and n1 >= 0 else []
        if not candidates:
            differences.append(f"{where(old, e)} is only in OLD")
            continue
        j = candidates[0]
        used.add(j)
        moved += j != k
        _compare_fields(where(old, e), {c: e[c] for c in keys}, {c: new_edges[j][c] for c in keys}, rtol, differences)
        for name in DERIVED:                              # values FEMM works out itself: noted, never a difference
            if name in e and not _close(e[name], new_edges[j][name], rtol):
                computed[name].append((where(old, e), e[name], new_edges[j][name]))
    for j, e in enumerate(new_edges):
        if j not in used:
            differences.append(f"{where(new, e)} is only in NEW")
    for name, found in computed.items():
        if found:
            first = found[0]
            formatting.append(f"{kind}s: {DERIVED[name]} differs on {len(found)} of {len(old_edges)}; "
                              f"first: {first[0]} is {_show(first[1])} in OLD and {_show(first[2])} in NEW")
    if moved:
        formatting.append(f"{kind}s: stored in a different order ({moved} of {len(old_edges)} have another number)")


def compare(old: dict, new: dict, rtol: float = 1e-9) -> dict:
    """Compare two parsed models (parse_fem). Returns {"differences", "formatting", "warnings", "compared"}.

    differences: one line for each way the model differs (an empty list means: the same model).
    formatting:  lines about what is not the model: how numbers are written, comments, order, line endings, and the
                 values FEMM computes for itself.
    warnings:    what either file held that this tool did not understand, so did not compare (see the module docstring).
    compared:    how many of each kind of thing were compared, so a clean result can be seen to cover something.
    rtol must be zero or positive; 0 asks for exactly equal numbers (positions are still matched to 1e-12).
    """
    if not rtol >= 0:
        raise ValueError(f"rtol must be zero or positive, not {rtol}")
    differences, formatting = [], []
    compared = {"header fields": _compare_header(old, new, rtol, differences, formatting)}
    compared.update(_compare_props(old, new, rtol, differences, formatting))
    coordinates = [abs(v) for model in (old, new) for item in (*model["nodes"], *model["labels"])
                   for v in (item["x"], item["y"])]
    size = max(coordinates) if coordinates else 1.0
    tol = max(rtol * size, MIN_POSITION_TOL)      # how close two positions must be to be one position
    to_new = _compare_points("node", old["nodes"], new["nodes"], tol, rtol, ("prop", "group", "lfs"), differences, formatting)
    _compare_edges("segment", old, new, to_new, rtol, ("maxside", "boundary", "hidden", "group"), differences, formatting)
    _compare_edges("arc", old, new, to_new, rtol, ("sweep", "maxside", "boundary", "hidden", "group"), differences, formatting)
    label_keys = ("material", "meshsize", "circuit", "magdir", "group", "turns", "external", "magdirfctn")
    _compare_points("label", old["labels"], new["labels"], tol, rtol, label_keys, differences, formatting)
    if old["crlf"] != new["crlf"]:
        ends = lambda m: "CRLF" if m["crlf"] else "LF"
        formatting.append(f"line endings: OLD uses {ends(old)}, NEW uses {ends(new)}")
    compared.update(nodes=len(old["nodes"]), segments=len(old["segments"]), arcs=len(old["arcs"]),
                    **{"labels (holes included)": len(old["labels"])})
    warnings = [f"{name}: {w}" for name, model in (("OLD", old), ("NEW", new)) for w in model.get("warnings", [])]
    return {"differences": differences, "formatting": formatting, "warnings": warnings, "compared": compared}


# ----------------------------------------------------------------------------- 5. the command line
def report(old_path, new_path, result: dict) -> str:
    """The text the command prints."""
    lines = [f"OLD  {old_path}", f"NEW  {new_path}",
             "compared: " + ", ".join(f"{n} {name}" for name, n in result["compared"].items() if n), "",
             f"MODEL DIFFERENCES: {len(result['differences'])}"]
    lines += [f"  {d}" for d in result["differences"]]
    lines += ["", f"NOT THE MODEL (formatting, notes, order, values FEMM computes): {len(result['formatting'])}"]
    lines += [f"  {d}" for d in result["formatting"]]
    warnings = result.get("warnings", [])
    if warnings:
        lines += ["", f"WARNINGS (parts of a file this tool did not compare): {len(warnings)}"]
        lines += [f"  {w}" for w in warnings]
    verdict = "the same model" if not result["differences"] else "NOT the same model"
    if warnings and not result["differences"]:
        verdict += " (as far as the parts it understands go: see the WARNINGS)"
    lines += ["", "RESULT: " + verdict]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Are two FEMM .fem files the same model? Prints every difference.")
    ap.add_argument("old", help="the first .fem (for example the one FEMM 4.2 wrote)")
    ap.add_argument("new", help="the second .fem (for example the one FEMM 5.0 wrote)")
    ap.add_argument("--rtol", type=float, default=1e-9,
                    help="two numbers are the same if they agree to this fraction of the larger one (default 1e-9; 0 = exactly equal)")
    args = ap.parse_args(argv)
    if not args.rtol >= 0:
        ap.error("--rtol must be zero or positive")
    models = []
    for label, path in (("OLD", args.old), ("NEW", args.new)):
        try:
            models.append(parse_fem(path))
        except (OSError, ValueError) as exc:
            print(f"could not read {label}: {exc}", file=sys.stderr)
            return 2
    result = compare(models[0], models[1], args.rtol)
    print(report(args.old, args.new, result))
    return 0 if not result["differences"] else 1


if __name__ == "__main__":
    sys.exit(main())
