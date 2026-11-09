"""Read FEMM's text files without FEMM: the geometry of a .fem, and the mesh + vector potential of a .ans.

A .fem is plain text: [NumPoints] then one "x y ..." line per node, [NumSegments] "n0 n1 ...",
[NumArcSegments] "n0 n1 sweep_deg ...", [NumBlockLabels] "x y material# meshsize circuit# magdir group turns ...".
A .ans is the same file plus a [Solution] section: one "x y A ..." line per mesh node, then one
"n0 n1 n2 label# ..." line per triangle. First-order triangles, so B is constant inside each one.
"""
import math

import numpy as np


def _lines(path):
    with open(path, encoding="latin-1") as f:
        return f.read().splitlines()


def _section(lines, tag):
    for i, ln in enumerate(lines):
        if ln.startswith(tag):
            n = int(ln.split("=")[1])
            return [ln2.split() for ln2 in lines[i + 1:i + 1 + n]]
    return []


def read_fem(path):
    """Geometry and labels of a .fem (or of the input part of a .ans)."""
    lines = _lines(path)
    names = lambda tag: [ln.split("=", 1)[1].strip().strip('"') for ln in lines if tag in ln]
    pts = [(float(r[0]), float(r[1])) for r in _section(lines, "[NumPoints]")]
    segs = [(int(r[0]), int(r[1]), int(r[5]) if len(r) > 5 else 0) for r in _section(lines, "[NumSegments]")]
    arcs = [(int(r[0]), int(r[1]), float(r[2]), int(r[6]) if len(r) > 6 else 0) for r in _section(lines, "[NumArcSegments]")]
    materials, circuits = names("<BlockName>"), names("<CircuitName>")
    labels = []                         # in file order: a triangle's label number in the .ans indexes this list
    for r in _section(lines, "[NumBlockLabels]"):
        typ, circ = int(r[2]), int(r[4])
        labels.append({"xy": (float(r[0]), float(r[1])),
                       "material": materials[typ - 1] if typ > 0 else "<No Mesh>",
                       "circuit": circuits[circ - 1] if circ > 0 else "", "magdir": float(r[5]),
                       "group": int(r[6]), "turns": int(float(r[7]))})
    labels += [{"xy": (float(r[0]), float(r[1])), "material": "<No Mesh>", "circuit": "", "magdir": 0.0,
                "group": int(r[2]), "turns": 0} for r in _section(lines, "[NumHoles]")]   # FEMM saves <No Mesh> as holes
    return {"points": pts, "segments": segs, "arcs": arcs, "labels": labels,
            "materials": materials, "circuits": circuits}


def arc_points(a, b, sweep_deg, step_deg=0.5):
    """Points along a FEMM arc: from a to b, counter-clockwise, through sweep_deg degrees."""
    ang = math.radians(sweep_deg)
    (ax, ay), (bx, by) = a, b
    d = math.dist(a, b)
    h = d / 2 / math.tan(ang / 2)                     # centre sits h to the left of the chord's midpoint
    cx, cy = (ax + bx) / 2 - (by - ay) / d * h, (ay + by) / 2 + (bx - ax) / d * h
    r, t0 = math.dist(a, (cx, cy)), math.atan2(ay - cy, ax - cx)
    n = max(2, math.ceil(sweep_deg / step_deg))
    return [a] + [(cx + r * math.cos(t0 + ang * i / n), cy + r * math.sin(t0 + ang * i / n)) for i in range(1, n)] + [b]


def edges_as_polylines(fem, step_deg=0.5):
    """Every segment and arc as a list of (x, y) points: what shapely needs to rebuild the faces."""
    P = fem["points"]
    out = [[P[a], P[b]] for a, b, _ in fem["segments"]]
    out += [arc_points(P[a], P[b], sw, step_deg) for a, b, sw, _ in fem["arcs"]]
    return out


def read_ans(path):
    """Mesh and solution of a magnetostatic .ans: node xy [mm], nodal A [Wb/m], triangles, triangle -> label index."""
    lines = _lines(path)
    i = next(k for k, ln in enumerate(lines) if ln.startswith("[Solution]"))
    nn = int(lines[i + 1])
    nodes = np.array([[float(v) for v in ln.split()[:3]] for ln in lines[i + 2:i + 2 + nn]])
    ne = int(lines[i + 2 + nn])
    el = np.array([[int(v) for v in ln.split()[:4]] for ln in lines[i + 3 + nn:i + 3 + nn + ne]])
    return {"xy": nodes[:, :2], "A": nodes[:, 2], "tri": el[:, :3], "label": el[:, 3], **read_fem(path)}


def element_b(xy, A, tri):
    """Flux density per triangle [T] from the nodal A of first-order elements: B = (dA/dy, -dA/dx), mm -> m."""
    x, y = xy[tri, 0], xy[tri, 1]
    a = A[tri]
    det = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    dadx = ((a[:, 1] - a[:, 0]) * (y[:, 2] - y[:, 0]) - (a[:, 2] - a[:, 0]) * (y[:, 1] - y[:, 0])) / det
    dady = ((a[:, 2] - a[:, 0]) * (x[:, 1] - x[:, 0]) - (a[:, 1] - a[:, 0]) * (x[:, 2] - x[:, 0])) / det
    return 1e3 * dady, -1e3 * dadx
