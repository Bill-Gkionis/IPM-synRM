"""optimization.thermal_fem: the heat solvers, calibrated on exact answers, and the thermal network checked against a 2D finite-element slice of the stator.

What it is, and why
    A network is only as good as the numbers it was built from, and a finite-element model is only as good as the solver under it. So this file
    does two things, in this order (the order is the rule: calibrate, then trust):
      1. CALIBRATION. FEMM 5.0's heat-flow solver and CalculiX 2.23's heat transfer solve cases whose answer is known on paper, and are held to
         tolerances that were written down BEFORE the first run (the constants below; they are not changed after seeing a result).
      2. THE CHECK. FEMM draws one slot pitch of the stator (thermal_geometry.py): the winding as a block of the equivalent conductivity with the
         slot's copper loss as a volume heat source, the liner, the teeth and the yoke with their iron loss, the contact to the housing as a thin
         layer, the housing, and the cooling's heat transfer coefficient on the outer surface. The same loads go into the network's reduced
         form (ThermalModel.radial_2d: no end turns, no rotor, nothing leaves axially), and the temperatures are compared.

The exact answers (each built from the one before)
    A. A long solid cylinder, radius R, conductivity k, uniform heat generation q [W/m^3], a convective surface (h, T_inf). All the heat made leaves
       through the surface: 2 pi R Q' = q pi R^2, so the surface is  T_s = T_inf + q R / (2 h). Inside, heat flowing out through radius r is
       q pi r^2 per metre, and Fourier's law  -k dT/dr = q r / 2  integrates to  T(r) = T_s + q (R^2 - r^2) / (4 k): the centre is q R^2 / (4 k)
       above the surface. This is the case the brief names. With R = 40 mm, k = 1, q = 1e5, h = 20: T_s = 100 K, T_centre = 140 K.
    B. The same idea in layers, the slot's own sequence: winding (conductivity k_w, generation), liner (thin, low k), iron (k_fe, generation), contact
       (thin, low k), aluminium, and a convective surface; an adiabatic hole at the centre (so that no mesh has to reach r = 0). The heat made inside radius
       r is Q'(r) = sum of q pi (r_out^2 - r_in^2) over the layers below r, and  T(r) = T(r_out of the layer) + the integral of Q'(s) / (2 pi k s) ds  from r to
       the layer's outer radius: exact by quadrature. It tests what the slot model needs: a block with a source, a thin low-conductivity layer, a jump in
       conductivity by a factor of 150, and the convective boundary.
    FEMM solves A and B; CalculiX solves B (a structured mesh of 8-node plane elements, one 10 degree wedge, the radial edges untouched, the outer edge
    under *FILM). The element is CPS8, the plane-stress type, which a heat transfer step uses as a 2D heat element: the first try with the heat element's own
    name, DC2D8, made CalculiX 2.23 lose the last three elements of the deck (a warning for each, "value N in set L4 > ne") and crash with an access violation;
    the same deck with CPS8 ran and converged. Each case is solved on three meshes: an answer that does not move toward the exact one as the mesh is refined
    is a solver that cannot be trusted, however well it agrees once.

The check against the network (what is compared, and the tolerance, fixed in advance)
    Loads (W per metre of stack): the copper loss in the slots, and the iron loss of the teeth and the yoke (thermal.py's placeholder at 200 Hz). Outer h.
    For each of several variants (a poorly and a well impregnated winding, a lower fill, a thicker liner, a liquid-cooled h, copper only, iron only) the
    slice and the network get the same numbers, and the comparison is made on
        T_hot    the hottest point of the winding block (FEMM: the largest nodal value inside it)
        T_mean   the winding block's area average
        T_teeth, T_yoke   the area averages of the two iron regions
    each as a rise above the housing's OUTER surface (the internal drops: what the slot, liner, iron and contact do) and above the ambient (which adds
    the outside's convection, identical in both by construction: it dilutes a disagreement, so the internal drop is the demanding number).
    Tolerance (fixed before the slice was first solved): the internal drops of T_hot and T_mean within 10 %, and the rises above the ambient within 5 %. The
    5 % applies to every variant in which the outside carries most of the rise (h of 15 and 60 W/(m^2 K)); in the liquid-cooled variant (h = 2,000) the
    outside is a fifth of the rise, "above the ambient" is nearly the same number as "above the outer surface", and only the 10 % applies. The variants
    without copper loss (iron only) are judged on T_teeth and T_yoke instead of T_hot and T_mean. The reasons the two could differ are listed in the
    report with each result: the equivalent rectangle, the tooth as a 1D fin, the yoke as a 1D ring, the uniform wall temperature the network assumes.

Run (from the Python to FEMM folder; each solver run takes the solver slot, waits for it, and stays within about 15 minutes)
    uv run python -m optimization.thermal_fem --only femm_cal        FEMM on cases A and B (three meshes each)
    uv run python -m optimization.thermal_fem --only ccx_cal         CalculiX on case B (three meshes)
    uv run python -m optimization.thermal_fem --only sector          the slice against the network, every variant (FEMM, one pool session)
    uv run python -m optimization.thermal_fem --only sector_mesh     the base variant on three mesh settings
    uv run python -m optimization.thermal_fem --only sector_split    the same with the teeth in three blocks (like-for-like diagnostics)
    uv run python -m optimization.thermal_fem                        all of them (the FEMM parts share one pool session)
    uv run python -m optimization.thermal_fem --reevaluate           apply the tolerances to the saved slice rows again (no solver)
    uv run python -m optimization.thermal_fem --rederive             recompute the network's side of the saved slice rows with the present thermal.py (no solver)
    (optimization/thermal_ccx.py solves the same slice in CalculiX.)
Writes out\\optimization\\thermal\\fem_*.json and fem_*.csv. Exit status 0 only if every calibration and every check is inside its tolerance.

Numbers to expect (this PC, 2026-11-09, FEMM 5.0 build 1144, CalculiX 2.23; the sessions took about a minute each)
    FEMM, case A (cylinder, exact 140.000 K at the centre and 100.000 K at the surface): 139.9969 and 99.9990 on the 4 and 2 mm settings (the same mesh: the arcs' 0.5 degree
        segments set the size, 5,075 nodes), 139.9990 and 99.9998 on 1 mm (7,169 nodes): errors 0.002 % and 0.0007 %, against a tolerance of 0.5 %.
    FEMM, case B (layered annulus, exact: hole wall 142.6453 K, winding mean 138.1562, liner drop 2.0768, contact drop 0.8339): the hole's wall 142.6498 / 142.6468 /
        142.6447 on 28,649 / 55,528 / 164,399 nodes (+0.0031 %, +0.0011 %, -0.0004 %), the winding mean -0.0075 / -0.0032 / -0.0016 %, the liner's drop -0.125 / +0.009 /
        -0.013 %, the contact layer's -0.014 / -0.010 / -0.004 %; 0.4 to 2.8 s a solve. Tolerances 1 % and 3 %: every row inside them, and the error falls with the mesh.
    CalculiX, case B (CPS8, 188 / 373 / 743 nodes): the hole's wall -0.0002 %, the winding mean -0.0067 / -0.0017 / -0.0005 %, the liner's drop +0.011 %, the contact layer's
        +0.016 % (the same on all three: the linear flow in a thin layer is exact in one quadratic element); 0.1 to 0.2 s a solve.
    The slice (66,000 to 75,000 nodes, 4 to 5 s a solve with FEMM, three workers): T_hot and T_mean as drops above the outer surface, network against finite elements:
        base +0.02 % and +0.74 %; poor impregnation +1.5 and +1.1; good impregnation +0.02 and +0.81; low fill +1.6 and +1.1; thick liner -0.8 and +0.3; liquid h +0.02 and +0.74
        (above the ambient +0.02 and +0.59); natural h +1.9 and +3.5; copper only -1.1 and -0.9. All inside the tolerance. The iron-only variant FAILS its criterion: T_teeth's
        drop is 18 % high (0.96 K against 0.81 K), because the FEMM region "teeth" is the whole tooth (shoe, fin and root) and the network's node is the fin's mean; with the teeth
        in three blocks the fin's mean is +7.5 % (0.07 K) and the tip +6.4 % (and +0.8 % and +0.9 % in the base variant). The three mesh settings move T_mean by 0.002 K.
    Flank / bottom split of the slot's heat: FEMM 13.6 / 6.0 W per metre (69 / 31 %) against the network's 15.1 / 4.9 W (75 / 25 %): the round bottom takes more than the rectangle's.
    One refinement was tried and reverted: sending the root's share (30 %) of the tooth iron loss into the yoke node improved the region average and made the like-for-like fin and
    hot-spot comparisons worse (hot spot +0.02 to -2.7 %, fin +0.8 to -4.0 % in the base variant).

Self-check (the answers follow).
  1. Why is a solid cylinder with an insulated centre hole still an exact test of the solver?
  2. Both programs agree with each other to 0.01 % and disagree with the exact solution by 2 %. What is wrong, and which program is to blame?
  3. The network's internal drop is 12 % smaller than the finite-element one. Why do both numbers have to be quoted, not only the rise above ambient?
Answers.
  1. The hole's wall carries no heat (nothing to enter or leave), which is the same as the symmetry of the solid cylinder about its centre, so the heat
     flow through any circle is the generation inside it and the exact solution above applies unchanged between the hole and the surface.
  2. Neither solver: a shared error points at the problem statement (the exact solution, the units, the geometry the two programs were given). The
     mesh refinement shows it: if both move together toward the exact answer, it is a mesh effect; if both stay, the statement is wrong.
  3. The outside is most of the rise and is the same in both, so a 12 % error in the small internal drop is 2 % of the total. It would hide a wrong
     slot model behind a right convection coefficient.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.spatial import cKDTree
from shapely.geometry import LineString
from shapely.ops import unary_union

import femm_queue
import motorgeo as mg
import paths
from optimization import thermal, thermal_geometry

# ---------------------------------------------------------------------------------------------------- tolerances: written down before any run
TOL_CYLINDER_PERCENT = 0.5            # case A: the centre and the surface rise above T_inf, and the heat flow out, within 0.5 %
TOL_COMPOSITE_PERCENT = 1.0           # case B: the temperature at the hole's wall and the mean over the winding layer, within 1 % of their rise
TOL_THIN_LAYER_PERCENT = 3.0          # case B: the drop across each thin layer, within 3 %
TOL_SECTOR_INTERNAL_PERCENT = 10.0    # slice: T_hot and T_mean above the housing's outer surface, network against finite elements
TOL_SECTOR_AMBIENT_PERCENT = 5.0      # slice: the same two above the ambient
POOL_WORKERS = 3
TASK_TIMEOUT_S = 300
OUT = lambda: paths.path("out") / "optimization" / "thermal"


# ---------------------------------------------------------------------------------------------------- the exact solutions
@dataclass(frozen=True)
class Layer:
    name: str
    r_in: float            # m
    r_out: float
    k: float               # W/(m K)
    q: float = 0.0         # W/m^3


def composite_exact(layers: list[Layer], h: float, t_inf: float) -> dict:
    """Exact 1D radial solution of the layered cylinder with an ADIABATIC inner boundary (case B). Returns a function T(r) [K, relative to the same
    scale as t_inf], the heat flow per metre, and helper means. Quadrature (scipy.integrate.quad), accurate to 1e-10."""
    def q_inside(r: float) -> float:
        tot = 0.0
        for L in layers:
            if r > L.r_in:
                tot += L.q * math.pi * (min(r, L.r_out) ** 2 - L.r_in ** 2)
        return tot
    q_total = q_inside(layers[-1].r_out)
    t_surface = t_inf + q_total / (2 * math.pi * layers[-1].r_out * h)

    def T(r: float) -> float:
        t = t_surface
        for L in reversed(layers):
            lo = max(r, L.r_in)
            if r >= L.r_out:
                continue
            val, _ = integrate.quad(lambda s, L=L: q_inside(s) / (2 * math.pi * L.k * s), lo, L.r_out, epsabs=1e-13, epsrel=1e-13)
            t += val
        return t

    def mean(layer: Layer) -> float:
        val, _ = integrate.quad(lambda r: T(r) * 2 * math.pi * r, layer.r_in, layer.r_out, epsabs=1e-12, epsrel=1e-11)
        return val / (math.pi * (layer.r_out ** 2 - layer.r_in ** 2))
    return {"T": T, "q_total": q_total, "t_surface": t_surface, "mean": mean}


def cylinder_exact(radius: float, k: float, q: float, h: float, t_inf: float) -> dict:
    """Case A: T_s = T_inf + q R / (2 h), T(r) = T_s + q (R^2 - r^2) / (4 k), heat flow q pi R^2 per metre."""
    ts = t_inf + q * radius / (2 * h)
    return {"T_surface": ts, "T_centre": ts + q * radius ** 2 / (4 * k), "q_total": q * math.pi * radius ** 2}


# the cases
CYL = dict(radius=0.040, k=1.0, q=1.0e5, h=20.0, t_inf=0.0)
COMPOSITE = [Layer("winding", 0.002, 0.012, 0.66, 2.0e5), Layer("liner", 0.012, 0.01225, 0.139), Layer("iron", 0.01225, 0.040, 21.9, 1.0e5),
             Layer("contact", 0.040, 0.0401, 0.259), Layer("aluminium", 0.0401, 0.045, 154.0)]
COMPOSITE_H, COMPOSITE_TINF = 15.0, 0.0


# ---------------------------------------------------------------------------------------------------- FEMM tasks (run in a pool worker)
def _first(x):
    return x[0] if isinstance(x, (list, tuple)) else x


def femm_cylinder_task(task: dict) -> dict:
    """Case A in FEMM: a solid cylinder, generation q, convection h. Returns the centre and surface temperature and the node count."""
    import femm
    started = time.perf_counter()
    c = task["case"]
    R = c["radius"] * 1e3
    mesh = task["mesh_mm"]
    path = str(Path(task["folder"]) / f"{task['tag']}.feh")
    femm.newdocument(2)
    femm.hi_probdef("millimeters", "planar", 1e-10, 1000.0, 30)
    femm.hi_addmaterial("mat", c["k"], c["k"], c["q"], 0)
    femm.hi_addboundprop("conv", 2, 0, 0, c["t_inf"], c["h"], 0)
    femm.hi_addnode(R, 0)
    femm.hi_addnode(-R, 0)
    femm.hi_addarc(R, 0, -R, 0, 180, 0.5)
    femm.hi_addarc(-R, 0, R, 0, 180, 0.5)
    femm.hi_addblocklabel(R / 2, 0)
    femm.hi_selectlabel(R / 2, 0)
    femm.hi_setblockprop("mat", 0, mesh, 0)
    femm.hi_clearselected()
    for y in (R, -R):
        femm.hi_selectarcsegment(0, y)
    femm.hi_setarcsegmentprop(0.5, "conv", 0, 0, "<None>")
    femm.hi_clearselected()
    femm.hi_saveas(path)
    femm.hi_analyze(1)
    femm.hi_loadsolution()
    t_c = _first(femm.ho_getpointvalues(0.0, 0.0))
    t_s = _first(femm.ho_getpointvalues(R - 1e-6, 0.0))
    nodes = int(femm.ho_numnodes())
    femm.ho_selectblock(R / 2, 0)
    t_mean = _first(femm.ho_blockintegral(0))
    femm.ho_clearblock()
    femm.ho_close()
    femm.hi_close()
    for ext in (".feh", ".anh"):
        try:
            os.remove(path[:-4] + ext)
        except OSError:
            pass
    return {"status": "ok", "T_centre": float(t_c), "T_surface": float(t_s), "T_mean": float(t_mean), "nodes": nodes, "mesh_mm": mesh,
            "seconds": time.perf_counter() - started}


def femm_composite_task(task: dict) -> dict:
    """Case B in FEMM: concentric rings (winding, liner, iron, contact, aluminium), an empty hole at the centre, convection outside."""
    import femm
    started = time.perf_counter()
    layers = [thermal_fem_layer(d) for d in task["layers"]]
    scale = task["mesh_scale"]
    path = str(Path(task["folder"]) / f"{task['tag']}.feh")
    femm.newdocument(2)
    femm.hi_probdef("millimeters", "planar", 1e-10, 1000.0, 30)
    femm.hi_addboundprop("conv", 2, 0, 0, task["t_inf"], task["h"], 0)
    radii = [layers[0].r_in] + [L.r_out for L in layers]
    for r in radii:
        r_mm = r * 1e3
        femm.hi_addnode(r_mm, 0)
        femm.hi_addnode(-r_mm, 0)
        femm.hi_addarc(r_mm, 0, -r_mm, 0, 180, 0.5)
        femm.hi_addarc(-r_mm, 0, r_mm, 0, 180, 0.5)
    base = {"winding": 0.6, "liner": 0.08, "iron": 1.2, "contact": 0.05, "aluminium": 0.6}
    femm.hi_addblocklabel(0, 0)                           # the hole at the centre is a closed region: FEMM wants it labelled "no mesh"
    femm.hi_selectlabel(0, 0)
    femm.hi_setblockprop("<No Mesh>", 0, 0, 0)
    femm.hi_clearselected()
    for L in layers:
        femm.hi_addmaterial(L.name, L.k, L.k, L.q, 0)
        x = (L.r_in + L.r_out) / 2 * 1e3
        femm.hi_addblocklabel(x, 0)
        femm.hi_selectlabel(x, 0)
        femm.hi_setblockprop(L.name, 0, base[L.name] * scale, 0)
        femm.hi_clearselected()
    ro = layers[-1].r_out * 1e3
    for y in (ro, -ro):
        femm.hi_selectarcsegment(0, y)
    femm.hi_setarcsegmentprop(0.5, "conv", 0, 0, "<None>")
    femm.hi_clearselected()
    femm.hi_saveas(path)
    femm.hi_analyze(1)
    femm.hi_loadsolution()
    out = {"status": "ok", "nodes": int(femm.ho_numnodes()), "mesh_scale": scale}
    out["T_hole_wall"] = float(_first(femm.ho_getpointvalues(layers[0].r_in * 1e3 + 1e-9, 0.0)))
    out["T_winding_out"] = float(_first(femm.ho_getpointvalues(layers[0].r_out * 1e3 - 1e-9, 0.0)))
    out["T_liner_out"] = float(_first(femm.ho_getpointvalues(layers[1].r_out * 1e3 + 1e-9, 0.0)))
    out["T_iron_out"] = float(_first(femm.ho_getpointvalues(layers[2].r_out * 1e3 - 1e-9, 0.0)))
    out["T_contact_out"] = float(_first(femm.ho_getpointvalues(layers[3].r_out * 1e3 + 1e-9, 0.0)))
    out["T_surface"] = float(_first(femm.ho_getpointvalues(ro - 1e-9, 0.0)))
    x = (layers[0].r_in + layers[0].r_out) / 2 * 1e3
    femm.ho_selectblock(x, 0)
    out["T_winding_mean"] = float(_first(femm.ho_blockintegral(0)))
    femm.ho_clearblock()
    femm.ho_close()
    femm.hi_close()
    for ext in (".feh", ".anh"):
        try:
            os.remove(path[:-4] + ext)
        except OSError:
            pass
    out["seconds"] = time.perf_counter() - started
    return out


def thermal_fem_layer(d: dict) -> Layer:
    return Layer(d["name"], d["r_in"], d["r_out"], d["k"], d["q"])


def _snap(points: list, tol: float = 1e-5) -> dict:
    """Map every vertex to the first vertex within tol (so that nearly equal corners of two polygons become one node)."""
    arr = np.array(points)
    tree = cKDTree(arr)
    rep = {}
    for i, p in enumerate(points):
        j = min(tree.query_ball_point(p, tol))
        rep[i] = j
    return rep


def noded_segments(regions: dict) -> tuple[list, list]:
    """All the regions' boundaries as one set of segments, noded at every shared vertex and every crossing: returns (nodes, segments (i, j))."""
    pts = []
    owner = []
    for name, poly in regions.items():
        for part in getattr(poly, "geoms", [poly]):
            for ring in [part.exterior, *part.interiors]:
                for pt in list(ring.coords)[:-1]:
                    pts.append(tuple(pt))
                    owner.append(name)
    rep = _snap(pts)
    snapped = [pts[rep[i]] for i in range(len(pts))]
    lines = []
    k = 0
    for name, poly in regions.items():
        for part in getattr(poly, "geoms", [poly]):
            for ring in [part.exterior, *part.interiors]:
                n = len(ring.coords) - 1
                ring_pts = snapped[k:k + n]
                k += n
                for a, b in zip(ring_pts, ring_pts[1:] + ring_pts[:1]):
                    if a != b:
                        lines.append(LineString([a, b]))
    merged = unary_union(lines)
    segs = []
    for g in getattr(merged, "geoms", [merged]):
        c = list(g.coords)
        segs += list(zip(c[:-1], c[1:]))
    # nodes
    nodes = {}
    out_segs = []
    for a, b in segs:
        ia = nodes.setdefault((round(a[0], 6), round(a[1], 6)), len(nodes))
        ib = nodes.setdefault((round(b[0], 6), round(b[1], 6)), len(nodes))
        if ia != ib:
            out_segs.append((min(ia, ib), max(ia, ib)))
    out_segs = sorted(set(out_segs))
    node_list = [None] * len(nodes)
    for xy, i in nodes.items():
        node_list[i] = xy
    return node_list, out_segs


def femm_sector_task(task: dict) -> dict:
    """The slice in FEMM: regions (polygons as vertex lists with material, heat generation and mesh size), the outer convection. Returns the
    temperatures the comparison needs."""
    import femm
    import shapely
    from shapely.geometry import MultiPolygon, Polygon
    started = time.perf_counter()
    regions = {}
    for name, d in task["regions"].items():
        parts = [Polygon(ring) for ring in d["parts"]]
        regions[name] = parts[0] if len(parts) == 1 else MultiPolygon(parts)
    nodes, segs = noded_segments(regions)
    path = str(Path(task["folder"]) / f"{task['tag']}.feh")
    femm.newdocument(2)
    femm.hi_probdef("millimeters", "planar", 1e-10, 1000.0, 30)
    for name, d in task["regions"].items():
        femm.hi_addmaterial(name, d["k"], d["k"], d["q"], 0)
    femm.hi_addboundprop("conv", 2, 0, 0, task["t_inf"], task["h"], 0)
    for x, y in nodes:
        femm.hi_addnode(x, y)
    for i, j in segs:
        femm.hi_addsegment(nodes[i][0], nodes[i][1], nodes[j][0], nodes[j][1])
    labels = {}
    for name, d in task["regions"].items():
        for part in getattr(regions[name], "geoms", [regions[name]]):
            pt = part.representative_point()
            femm.hi_addblocklabel(pt.x, pt.y)
            femm.hi_selectlabel(pt.x, pt.y)
            femm.hi_setblockprop(name, 0, d["mesh"], 0)
            femm.hi_clearselected()
            labels.setdefault(name, []).append((pt.x, pt.y))
    r_out = task["r_outer"]
    n_outer = 0
    for i, j in segs:
        (x1, y1), (x2, y2) = nodes[i], nodes[j]
        if abs(math.hypot(x1, y1) - r_out) < 1e-4 and abs(math.hypot(x2, y2) - r_out) < 1e-4:
            femm.hi_selectsegment((x1 + x2) / 2, (y1 + y2) / 2)
            femm.hi_setsegmentprop("conv", 0, 1, 0, 0, "<None>")
            femm.hi_clearselected()
            n_outer += 1
    femm.hi_saveas(path)
    femm.hi_analyze(1)
    femm.hi_loadsolution()
    nn = int(femm.ho_numnodes())
    xy = np.array(femm.ho_getallnodes(), float)
    tt = np.array(femm.ho_getallt(), float)
    out = {"status": "ok", "nodes": nn, "outer_segments": n_outer}
    means, areas = {}, {}
    for name, pts in labels.items():
        femm.ho_clearblock()
        for x, y in pts:
            femm.ho_selectblock(x, y)
        means[name] = float(_first(femm.ho_blockintegral(0)))
        areas[name] = float(_first(femm.ho_blockintegral(1)))
        femm.ho_clearblock()
    out["mean"] = means
    out["area"] = areas
    # the hottest node in each region, and the outer surface's mean temperature
    for name in ("winding", "teeth", "teeth_shoe", "teeth_fin", "teeth_root", "yoke"):
        if name not in regions:
            continue
        mask = shapely.contains_xy(regions[name].buffer(1e-6), xy[:, 0], xy[:, 1])
        out.setdefault("max", {})[name] = float(tt[mask].max()) if mask.any() else float("nan")
    r_bore = task.get("r_bore")
    if r_bore:                                                      # the bore surface (the tooth tips): the mean over the nodes on it
        on_bore = np.abs(np.hypot(xy[:, 0], xy[:, 1]) - r_bore) < 1e-4
        if on_bore.any():
            out["T_bore"] = float(tt[on_bore].mean())
    r = np.hypot(xy[:, 0], xy[:, 1])
    ring = np.abs(r - r_out) < 1e-4
    ang = np.degrees(np.arctan2(xy[ring, 1], xy[ring, 0]))
    order = np.argsort(ang)
    a_s, t_s = ang[order], tt[ring][order]
    out["T_outer"] = float(np.trapezoid(t_s, a_s) / (a_s[-1] - a_s[0]))
    # heat through the slot's wall: both flanks (straight contours) and the bottom (a polyline), at the slot's outline
    out["flux"] = _slot_fluxes(femm, task["slot_outline_parts"])
    femm.ho_close()
    femm.hi_close()
    for ext in (".feh", ".anh"):
        try:
            os.remove(path[:-4] + ext)
        except OSError:
            pass
    out["seconds"] = time.perf_counter() - started
    return out


def _slot_fluxes(femm, parts: dict) -> dict:
    """Heat flow [W per metre] through the slot's flanks and bottom, from line integrals of the flux along the slot outline's pieces."""
    res = {}
    for name, pts in parts.items():
        femm.ho_clearcontour()
        femm.ho_seteditmode("contour")
        femm.ho_addcontour(*pts[0])
        for p in pts[1:]:
            femm.ho_addcontour(*p)
        v = femm.ho_lineintegral(1)
        res[name] = float(abs(_first(v)))
        femm.ho_clearcontour()
    return res


# ---------------------------------------------------------------------------------------------------- CalculiX: the layered annulus
def ccx_number(value: float) -> str:
    """A number as it goes into a deck: 12 significant digits (CalculiX reads 20 characters of a field and cuts the rest)."""
    return f"{value:.12g}"


def write_ccx_annulus(path: Path, layers: list[Layer], rows: list[int], h: float, t_inf: float, wedge_deg: float = 10.0) -> dict:
    """A structured mesh of one wedge of the layered annulus (units mm, tonne, s, so k in mW/(mm K) = W/(m K) numerically, h in mW/(mm^2 K) = 1e-3
    x W/(m^2 K), heat generation in mW/mm^3 = 1e-6 x W/m^3). Elements CPS8 (8-node plane elements, used as 2D heat elements), `rows[i]` radial elements in layer i and one
    element across the wedge; the radial edges are untouched (adiabatic), the outer edge carries *FILM. Returns the radial nodes (r_mm, node id) along
    the wedge's first edge."""
    n_t = 1
    radii_edges = []                                     # (layer index, r0, r1) of every element row
    for li, (L, nr) in enumerate(zip(layers, rows)):
        # a geometric grading would sharpen the thin layers, but uniform rows keep the check simple: the thin layers have many rows of their own
        edges = np.linspace(L.r_in * 1e3, L.r_out * 1e3, nr + 1)
        for a, b in zip(edges[:-1], edges[1:]):
            radii_edges.append((li, a, b))
    n_r = len(radii_edges)
    theta = math.radians(wedge_deg)
    nid = {}
    coords = {}

    def node(i: int, j: int) -> int:
        key = (i, j)
        if key not in nid:
            nid[key] = len(nid) + 1
        return nid[key]
    elems = []
    r_grid = []                                          # radius of grid index i (0 ... 2 n_r)
    for ie, (li, a, b) in enumerate(radii_edges):
        r_grid_vals = (a, (a + b) / 2, b)
        i0 = 2 * ie
        quad = []
        for (di, dj) in ((0, 0), (2, 0), (2, 2), (0, 2), (1, 0), (2, 1), (1, 2), (0, 1)):
            i, j = i0 + di, dj
            r = r_grid_vals[di]
            ang = theta * j / (2 * n_t)
            n = node(i, j)
            coords[n] = (r * math.cos(ang), r * math.sin(ang))
            quad.append(n)
        elems.append((li, quad))
    lines = ["*HEADING", "layered annulus, one wedge, CPS8 used as a 2D heat element", "*NODE"]
    for n in sorted(coords):
        lines.append(f"{n}, {ccx_number(coords[n][0])}, {ccx_number(coords[n][1])}, 0.0")
    by_layer = {}
    for ei, (li, quad) in enumerate(elems, start=1):
        by_layer.setdefault(li, []).append(ei)
    lines.append("*ELEMENT, TYPE=CPS8, ELSET=EALL")
    for ei, (li, quad) in enumerate(elems, start=1):
        lines.append(f"{ei}, " + ", ".join(str(n) for n in quad))
    for li, ids in by_layer.items():
        lines.append(f"*ELSET, ELSET=L{li}")
        lines += [", ".join(str(x) for x in ids[k:k + 10]) for k in range(0, len(ids), 10)]
    for li, L in enumerate(layers):
        lines += [f"*MATERIAL, NAME=M{li}", "*CONDUCTIVITY", ccx_number(L.k), f"*SOLID SECTION, ELSET=L{li}, MATERIAL=M{li}", "1000."]
    outer = by_layer[len(layers) - 1]
    lines += ["*INITIAL CONDITIONS, TYPE=TEMPERATURE", f"NALL, {ccx_number(t_inf)}"]
    lines += ["*STEP", "*HEAT TRANSFER, STEADY STATE", "1., 1."]
    for li, L in enumerate(layers):
        if L.q:
            lines += ["*DFLUX", f"L{li}, BF, {ccx_number(L.q * 1e-6)}"]
    lines.append("*FILM")
    for ei in outer[-n_t:]:                              # the last row of elements: its outer edge (face 2) is the surface
        lines.append(f"{ei}, F2, {ccx_number(t_inf)}, {ccx_number(h * 1e-3)}")
    lines += ["*NODE FILE", "NT", "*END STEP"]
    # the node set for the initial condition
    insert = lines.index("*MATERIAL, NAME=M0")
    lines[insert:insert] = ["*NSET, NSET=NALL, GENERATE", f"1, {len(coords)}, 1"]
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")       # LF endings, as the bench writes its decks
    return {"n_nodes": len(coords), "n_elements": len(elems), "edge_nodes": {i: n for (i, j), n in nid.items() if j == 0}, "grid_r": radii_edges}


def read_frd_temperature(path: Path) -> dict[int, float]:
    """Nodal temperatures from CalculiX's ASCII .frd (the block that follows ' -4  NDTEMP' or ' -4  T')."""
    out, block = {}, None
    with open(path, encoding="ascii", errors="replace") as f:
        for line in f:
            if line.startswith(" -4"):
                block = line.split()[1]
            elif line.startswith(" -3"):
                block = None
            elif block and block.upper() in ("NDTEMP", "T", "TEMP") and line.startswith(" -1"):
                node = int(line[3:13])
                out[node] = float(line[13:25])
    return out


def run_ccx_annulus(folder: Path, tag: str, rows: list[int], layers=None, h=COMPOSITE_H, t_inf=COMPOSITE_TINF) -> dict:
    """Case B in CalculiX through femm_queue.run_solver: write the deck, solve, read the temperatures along the wedge's first edge, return the same
    quantities as the FEMM task."""
    layers = layers or COMPOSITE
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    meta = write_ccx_annulus(folder / f"{tag}.inp", layers, rows, h, t_inf)
    run = femm_queue.run_solver([paths.tool("ccx"), "-i", tag], label=f"thermal_fem ccx {tag}", cwd=folder, timeout_s=300, threads=2)
    if run.timed_out or run.returncode not in (0, None) or not (folder / f"{tag}.frd").exists():
        return {"status": "error", "returncode": run.returncode, "stderr": run.stderr[-500:], "stdout": run.stdout[-500:]}
    temps = read_frd_temperature(folder / f"{tag}.frd")
    prof = sorted((layers_r, temps[n]) for layers_r, n in _edge_profile(meta).items())
    r_arr = np.array([r for r, _ in prof]) * 1e-3
    t_arr = np.array([t for _, t in prof])
    out = {"status": "ok", "nodes": meta["n_nodes"], "elements": meta["n_elements"], "seconds": run.seconds}

    def at(r):                                          # temperature at a radius that is a grid node (to 1e-9)
        i = int(np.argmin(np.abs(r_arr - r)))
        assert abs(r_arr[i] - r) < 1e-9, (r, r_arr[i])
        return float(t_arr[i])
    out["T_hole_wall"] = at(layers[0].r_in)
    out["T_winding_out"] = at(layers[0].r_out)
    out["T_liner_out"] = at(layers[1].r_out)
    out["T_iron_out"] = at(layers[2].r_out)
    out["T_contact_out"] = at(layers[3].r_out)
    out["T_surface"] = at(layers[-1].r_out)
    sel = (r_arr >= layers[0].r_in - 1e-12) & (r_arr <= layers[0].r_out + 1e-12)
    out["T_winding_mean"] = float(np.trapezoid(t_arr[sel] * r_arr[sel], r_arr[sel]) * 2 / (layers[0].r_out ** 2 - layers[0].r_in ** 2))
    return out


def _edge_profile(meta: dict) -> dict:
    """{radius [mm]: node id} along the wedge's first edge, from the grid."""
    prof = {}
    for ie, (li, a, b) in enumerate(meta["grid_r"]):
        for di, r in ((0, a), (1, (a + b) / 2), (2, b)):
            prof[round(r, 9)] = meta["edge_nodes"][2 * ie + di]
    return prof


# ---------------------------------------------------------------------------------------------------- the calibration drivers
def _pct(value: float, exact: float, scale: float | None = None) -> float:
    """Error as a percentage of `scale` (default: the exact value)."""
    return 100.0 * (value - exact) / (exact if scale is None else scale)


def calibrate_femm(pool) -> dict:
    """Cases A and B in FEMM, coarse and fine; returns the rows and whether every one is inside its tolerance."""
    folder = OUT() / "fem_cal"
    folder.mkdir(parents=True, exist_ok=True)
    tasks = []
    for mesh in (4.0, 2.0, 1.0):
        tasks.append({"kind": "A", "case": CYL, "mesh_mm": mesh, "folder": str(folder), "tag": f"cylA_{mesh:g}"})
    layers = [dict(name=L.name, r_in=L.r_in, r_out=L.r_out, k=L.k, q=L.q) for L in COMPOSITE]
    for scale in (2.0, 1.0, 0.5):
        tasks.append({"kind": "B", "layers": layers, "mesh_scale": scale, "t_inf": COMPOSITE_TINF, "h": COMPOSITE_H, "folder": str(folder), "tag": f"compB_{scale:g}"})
    results = pool.map(femm_calibration_dispatch, tasks)
    rows, ok = [], True
    ex_a = cylinder_exact(**CYL)
    ex_b = composite_exact(COMPOSITE, COMPOSITE_H, COMPOSITE_TINF)
    for t, r in zip(tasks, results):
        if r.get("status") != "ok":
            rows.append({"case": t["tag"], "status": r})
            continue
        if t["kind"] == "A":
            row = {"case": t["tag"], "mesh": t["mesh_mm"], "nodes": r["nodes"], "seconds": r["seconds"],
                   "T_centre": r["T_centre"], "T_centre_exact": ex_a["T_centre"], "err_centre_%": _pct(r["T_centre"] - CYL["t_inf"], ex_a["T_centre"] - CYL["t_inf"]),
                   "T_surface": r["T_surface"], "T_surface_exact": ex_a["T_surface"], "err_surface_%": _pct(r["T_surface"] - CYL["t_inf"], ex_a["T_surface"] - CYL["t_inf"]),
                   "T_mean": r["T_mean"]}
            ex_mean = ex_a["T_surface"] + CYL["q"] * CYL["radius"] ** 2 / (8 * CYL["k"])
            row["err_mean_%"] = _pct(r["T_mean"], ex_mean)
            row["pass"] = abs(row["err_centre_%"]) <= TOL_CYLINDER_PERCENT and abs(row["err_surface_%"]) <= TOL_CYLINDER_PERCENT
        else:
            row = _composite_row(t["tag"], r, ex_b)
            row["mesh"] = t["mesh_scale"]
        rows.append(row)
    # the finest mesh of each case must pass its tolerance; the coarse ones are there to show the convergence
    finest_a = [r for r in rows if r["case"].startswith("cylA")][-1]
    finest_b = [r for r in rows if r["case"].startswith("compB")][-1]
    ok = bool(finest_a["pass"] and finest_b["pass"])
    return {"rows": rows, "ok": ok, "tolerances": {"cylinder_%": TOL_CYLINDER_PERCENT, "composite_%": TOL_COMPOSITE_PERCENT, "thin_layer_%": TOL_THIN_LAYER_PERCENT}}


def femm_calibration_dispatch(task: dict) -> dict:
    return femm_cylinder_task(task) if task["kind"] == "A" else femm_composite_task(task)


def _composite_row(tag: str, r: dict, ex: dict) -> dict:
    T = ex["T"]
    L = COMPOSITE
    t_hole, t_w_out, t_l_out, t_i_out, t_c_out, t_s = (T(L[0].r_in), T(L[0].r_out), T(L[1].r_out), T(L[2].r_out), T(L[3].r_out), T(L[-1].r_out))
    rise = t_hole - COMPOSITE_TINF
    mean = ex["mean"](L[0])
    row = {"case": tag, "nodes": r["nodes"], "seconds": r["seconds"],
           "T_hole_wall": r["T_hole_wall"], "T_hole_wall_exact": t_hole, "err_hole_wall_%": _pct(r["T_hole_wall"] - COMPOSITE_TINF, rise),
           "T_winding_mean": r["T_winding_mean"], "T_winding_mean_exact": mean, "err_winding_mean_%": _pct(r["T_winding_mean"] - COMPOSITE_TINF, mean - COMPOSITE_TINF),
           "liner_drop": r["T_winding_out"] - r["T_liner_out"], "liner_drop_exact": t_w_out - t_l_out,
           "contact_drop": r["T_iron_out"] - r["T_contact_out"], "contact_drop_exact": t_i_out - t_c_out}
    row["err_liner_drop_%"] = _pct(row["liner_drop"], row["liner_drop_exact"])
    row["err_contact_drop_%"] = _pct(row["contact_drop"], row["contact_drop_exact"])
    row["pass"] = (abs(row["err_hole_wall_%"]) <= TOL_COMPOSITE_PERCENT and abs(row["err_winding_mean_%"]) <= TOL_COMPOSITE_PERCENT
                   and abs(row["err_liner_drop_%"]) <= TOL_THIN_LAYER_PERCENT and abs(row["err_contact_drop_%"]) <= TOL_THIN_LAYER_PERCENT)
    return row


def calibrate_ccx() -> dict:
    """Case B in CalculiX, coarse and fine (each solve takes the solver slot)."""
    folder = OUT() / "ccx_cal"
    ex = composite_exact(COMPOSITE, COMPOSITE_H, COMPOSITE_TINF)
    rows = []
    for tag, rows_per_layer in (("coarse", [10, 2, 20, 1, 4]), ("medium", [20, 4, 40, 2, 8]), ("fine", [40, 8, 80, 4, 16])):
        r = run_ccx_annulus(folder / tag, tag, rows_per_layer)
        if r.get("status") != "ok":
            rows.append({"case": tag, "status": r})
            continue
        row = _composite_row(tag, r, ex)
        row["rows_per_layer"] = rows_per_layer
        rows.append(row)
    finest = rows[-1]
    return {"rows": rows, "ok": bool(finest.get("pass", False)), "tolerances": {"composite_%": TOL_COMPOSITE_PERCENT, "thin_layer_%": TOL_THIN_LAYER_PERCENT}}


# ---------------------------------------------------------------------------------------------------- the slice against the network
def sector_variant_inputs(variant: dict, p: mg.Params) -> tuple[thermal.Assumptions, dict]:
    """(assumptions, loads) of a variant: the loads are W per metre of stack for ALL slots, teeth and yoke."""
    a = thermal.Assumptions(**variant.get("assume", {}))
    return a, variant


def build_sector_task(variant: dict, p: mg.Params, mesh_scale: float, folder: Path, tag: str, split_teeth: bool = False) -> tuple[dict, dict]:
    """The pool task of the slice for a variant, and the network's matching inputs."""
    a = thermal.Assumptions(**variant.get("assume", {}))
    fill = variant.get("fill", 0.45)
    model = thermal.ThermalModel(stack_mm=1000.0, cooling=thermal.Cooling.natural(), assume=a, winding=thermal.WindingSpec(fill))
    s = model._slot_branch(fill)
    t_liner = a.liner_mm
    r_contact_area = a.contact_factor * a.contact_nrel_mm2kw * 1e-6                 # m^2 K/W
    t_contact = 0.1
    k_contact = t_contact * 1e-3 / r_contact_area
    sl = thermal_geometry.build_slice(p, t_liner_mm=t_liner, t_contact_mm=t_contact, r_housing_out_mm=model.r_ho * 1e3, split_teeth=split_teeth)
    q_slot = variant["p_slot"]                                                      # W per metre in ONE slot
    area_block_m2 = sl.regions["winding"].area * 1e-6
    q_w = q_slot / area_block_m2                                                    # W/m^3 in the block
    f_t, f_y = variant.get("iron_teeth_w_m3", 0.0), variant.get("iron_yoke_w_m3", 0.0)
    mesh = {"winding": 0.25, "liner": 0.1, "teeth": 0.4, "teeth_shoe": 0.4, "teeth_fin": 0.4, "teeth_root": 0.4, "yoke": 0.5, "contact": 0.1, "housing": 0.8}
    regions = {}
    k_fe = thermal.MATERIALS["m19"]["k_in"]
    mat = {"winding": (s["k_eq"], q_w), "liner": (thermal.MATERIALS["liner"]["k"], 0.0), "teeth": (k_fe, f_t), "teeth_shoe": (k_fe, f_t),
           "teeth_fin": (k_fe, f_t), "teeth_root": (k_fe, f_t), "yoke": (k_fe, f_y), "contact": (k_contact, 0.0), "housing": (thermal.MATERIALS["aluminium"]["k"], 0.0)}
    for name, poly in sl.regions.items():
        parts = [list(part.exterior.coords)[:-1] for part in getattr(poly, "geoms", [poly])]
        regions[name] = {"parts": parts, "k": mat[name][0], "q": mat[name][1], "mesh": mesh[name] * mesh_scale}
    # the slot's outline pieces for the heat-flow contours: the two flanks and the bottom arc
    sd = p.slot_half_deg
    a_s = math.radians(sd)
    r_b = p.r_wall * math.tan(a_s)
    c_b = p.r_wall / math.cos(a_s)
    flank_lo = [thermal_geometry.pol(p.r_shoe, -sd), thermal_geometry.pol(p.r_wall, -sd)]
    flank_hi = [thermal_geometry.pol(p.r_shoe, sd), thermal_geometry.pol(p.r_wall, sd)]
    arc = thermal_geometry.arc_points(r_b, -(90.0 + sd), 90.0 + sd, step=5.0, cx=c_b)
    task = {"regions": regions, "h": variant["h"], "t_inf": 0.0, "r_outer": model.r_ho * 1e3, "r_bore": p.r_si, "folder": str(folder), "tag": tag,
            "slot_outline": list(sl.slot.exterior.coords), "slot_outline_parts": {"flank_lo": flank_lo, "flank_hi": flank_hi, "bottom": arc}}
    areas = {k: v.area * 1e-6 for k, v in sl.regions.items()}
    if split_teeth:
        areas["teeth"] = areas["teeth_shoe"] + areas["teeth_fin"] + areas["teeth_root"]
    net = {"model": model, "fill": fill, "areas_m2": areas, "q_w": q_w, "t_contact": t_contact, "k_contact": k_contact}
    return task, net


def network_for_variant(variant: dict, net: dict, p: mg.Params) -> dict:
    """The network's radial_2d answer for the same loads: all 36 slots' worth of power per metre."""
    model = net["model"]
    n = p.slots
    areas = net["areas_m2"]
    p_slot = variant["p_slot"] * n
    p_teeth = variant.get("iron_teeth_w_m3", 0.0) * areas["teeth"] * n
    p_yoke = variant.get("iron_yoke_w_m3", 0.0) * areas["yoke"] * n
    out = model.radial_2d(p_slot, p_teeth, p_yoke, variant["h"], fill=net["fill"], t_amb=0.0)
    q_total = p_slot + p_teeth + p_yoke
    out["T_outer"] = q_total / (variant["h"] * 2 * math.pi * model.r_ho)                      # the outer surface, per metre
    out["q_total"] = q_total
    return out


SECTOR_VARIANTS = {
    # loads per metre of stack: 20 W in each slot (all 36: 720 W/m, about the copper loss of the forced-air case's 60 mm stack), the iron-loss placeholder at 200 Hz as
    # volume densities (teeth 1.19e5 W/m^3 = 15.6 W/kg, yoke 7.5e4 = 9.9 W/kg); h = 60 W/(m^2 K) on the housing's outer cylinder: the whole housing's conductance
    # (end bells, fins) divided by the cylinder's area, in the range the three coolings give
    "base": {"p_slot": 20.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 60.0, "fill": 0.45},
    "poor_impregnation": {"p_slot": 20.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 60.0, "fill": 0.45, "assume": {"k_matrix": 0.15}},
    "good_impregnation": {"p_slot": 20.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 60.0, "fill": 0.45, "assume": {"k_matrix": 0.8}},
    "low_fill": {"p_slot": 20.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 60.0, "fill": 0.30},
    "thick_liner": {"p_slot": 20.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 60.0, "fill": 0.45, "assume": {"liner_mm": 0.40}},
    "liquid_h": {"p_slot": 20.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 2000.0, "fill": 0.45},
    "natural_h": {"p_slot": 5.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 15.0, "fill": 0.45},
    "copper_only": {"p_slot": 20.0, "iron_teeth_w_m3": 0.0, "iron_yoke_w_m3": 0.0, "h": 60.0, "fill": 0.45},
    "iron_only": {"p_slot": 0.0, "iron_teeth_w_m3": 1.19e5, "iron_yoke_w_m3": 7.5e4, "h": 60.0, "fill": 0.45},
}


def compare_sector(pool, p: mg.Params, variants: dict | None = None, mesh_scale: float = 1.0, split_teeth: bool = False) -> dict:
    """Solve the slice for every variant and compare with the network. Returns the rows, and whether the base and every variant are inside tolerance."""
    variants = variants or SECTOR_VARIANTS
    folder = OUT() / "fem_sector"
    folder.mkdir(parents=True, exist_ok=True)
    tasks, nets = [], []
    for name, v in variants.items():
        t, n = build_sector_task(v, p, mesh_scale, folder, f"sector_{name}_{mesh_scale:g}{'_split' if split_teeth else ''}", split_teeth=split_teeth)
        tasks.append(t)
        nets.append(n)
    results = pool.map(femm_sector_task, tasks)
    rows = []
    for (name, v), t, n, r in zip(variants.items(), tasks, nets, results):
        if r.get("status") != "ok":
            rows.append({"variant": name, "status": r})
            continue
        net = network_for_variant(v, n, p)
        row = {"variant": name, "mesh_scale": mesh_scale, "nodes": r["nodes"], "seconds": r["seconds"], "k_eq": net["k_eq"], "f_loc": net["f_loc"]}
        if split_teeth:
            a_s, a_f, a_r = n["areas_m2"]["teeth_shoe"], n["areas_m2"]["teeth_fin"], n["areas_m2"]["teeth_root"]
            t_teeth = (r["mean"]["teeth_shoe"] * a_s + r["mean"]["teeth_fin"] * a_f + r["mean"]["teeth_root"] * a_r) / (a_s + a_f + a_r)
        else:
            t_teeth = r["mean"]["teeth"]
        fe = {"T_hot": r["max"]["winding"], "T_mean": r["mean"]["winding"], "T_teeth": t_teeth, "T_yoke": r["mean"]["yoke"], "T_outer": r["T_outer"]}
        nw = {"T_hot": net["T_hot"], "T_mean": net["T_W"], "T_teeth": net["T_T"], "T_yoke": net["T_Y"], "T_outer": net["T_outer"]}
        keys_cmp = ["T_hot", "T_mean", "T_teeth", "T_yoke"]
        if split_teeth:                              # like for like: the network's tooth node is the fin's mean, and its tip is the bore's end of the fin
            fe["T_fin"], nw["T_fin"] = r["mean"]["teeth_fin"], net["T_T"]
            keys_cmp.append("T_fin")
            if "T_bore" in r:
                fe["T_tip"], nw["T_tip"] = r["T_bore"], net["T_tip"]
                keys_cmp.append("T_tip")
        for key in keys_cmp:
            fe_int, nw_int = fe[key] - fe["T_outer"], nw[key] - nw["T_outer"]
            row[f"{key}_fe"], row[f"{key}_net"] = fe[key], nw[key]
            row[f"{key}_drop_fe"], row[f"{key}_drop_net"] = fe_int, nw_int
            row[f"{key}_drop_err_%"] = _pct(nw_int, fe_int) if abs(fe_int) > 1e-9 else float("nan")
            row[f"{key}_amb_err_%"] = _pct(nw[key], fe[key]) if abs(fe[key]) > 1e-9 else float("nan")
        row["T_outer_fe"], row["T_outer_net"] = fe["T_outer"], nw["T_outer"]
        row["q_flank_fe"], row["q_bottom_fe"] = r["flux"]["flank_lo"] + r["flux"]["flank_hi"], r["flux"]["bottom"]
        row["q_flank_net_W_per_slot"], row["q_bottom_net_W_per_slot"] = net["q_flank"] / p.slots, net["q_bottom"] / p.slots
        row["heat_made_slice_W_per_m"] = sum(d["q"] * n["areas_m2"][name_] for name_, d in t["regions"].items())        # every region's source times its area
        row["has_copper"] = v["p_slot"] > 0
        row["h_outer"] = v["h"]
        rows.append(row)
    return evaluate_sector_rows(rows)


def rederive_network(name: str, split: bool) -> dict:
    """Recompute the network's side of a saved slice file with the present thermal.py (the solvers' temperatures stay as they were solved), and judge it again. The file as
    it was is kept as <name>_before_rederive.json the first time. Used when the network changes after the finite-element runs (it did once: the root's share of the tooth
    iron loss now goes to the yoke node, see thermal.py)."""
    path = OUT() / f"{name}.json"
    backup = OUT() / f"{name}_before_rederive.json"
    if not backup.exists():
        backup.write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
    saved = json.loads(path.read_text(encoding="utf-8"))
    p = thermal.design_params()
    for row in saved["rows"]:
        if row.get("status", "ok") != "ok" or "variant" not in row:
            continue
        v = SECTOR_VARIANTS[row["variant"]]
        _, n = build_sector_task(v, p, 1.0, OUT(), "x", split_teeth=split)
        net = network_for_variant(v, n, p)
        suffix = "fe" if "T_hot_fe" in row else "ccx"
        nw = {"T_hot": net["T_hot"], "T_mean": net["T_W"], "T_teeth": net["T_T"], "T_yoke": net["T_Y"], "T_fin": net["T_T"], "T_tip": net["T_tip"]}
        t_out_fe = row[f"T_outer_{suffix}"]
        for key, value in nw.items():
            if f"{key}_{suffix}" not in row:
                continue
            fe_v = row[f"{key}_{suffix}"]
            fe_d, nw_d = fe_v - t_out_fe, value - net["T_outer"]
            row[f"{key}_net"], row[f"{key}_drop_net"] = value, nw_d
            row[f"{key}_drop_{suffix}" if suffix == "ccx" else f"{key}_drop_fe"] = fe_d
            row[f"{key}_drop_err_%"] = _pct(nw_d, fe_d) if abs(fe_d) > 1e-9 else float("nan")
            row[f"{key}_amb_err_%"] = _pct(value, fe_v) if abs(fe_v) > 1e-9 else float("nan")
        row["T_outer_net"] = net["T_outer"]
        if "q_flank_net_W_per_slot" in row:
            row["q_flank_net_W_per_slot"], row["q_bottom_net_W_per_slot"] = net["q_flank"] / p.slots, net["q_bottom"] / p.slots
        row.setdefault("has_copper", v["p_slot"] > 0)
        row.setdefault("h_outer", v["h"])
    res = evaluate_sector_rows(saved["rows"])
    _write(name, res)
    return res


def evaluate_sector_rows(rows: list[dict]) -> dict:
    """Apply the tolerances (module docstring) to the rows: T_hot and T_mean (T_teeth and T_yoke where there is no copper loss) within
    TOL_SECTOR_INTERNAL_PERCENT as drops above the outer surface, and within TOL_SECTOR_AMBIENT_PERCENT as rises above the ambient, except where the
    outside is a small part of the rise (h_outer above 1,000: only the internal tolerance applies). Sets row["pass"]; returns the summary."""
    for row in rows:
        if "T_hot_drop_err_%" not in row:
            row["pass"] = False
            continue
        keys = ("T_hot", "T_mean") if row.get("has_copper", True) else ("T_teeth", "T_yoke")
        outside_small = row.get("h_outer", 0.0) > 1000.0
        row["pass"] = all(abs(row[f"{k}_drop_err_%"]) <= TOL_SECTOR_INTERNAL_PERCENT
                          and (outside_small or abs(row[f"{k}_amb_err_%"]) <= TOL_SECTOR_AMBIENT_PERCENT) for k in keys)
    return {"rows": rows, "ok": all(r.get("pass", False) for r in rows), "tolerances": {"internal_%": TOL_SECTOR_INTERNAL_PERCENT, "ambient_%": TOL_SECTOR_AMBIENT_PERCENT}}


# ---------------------------------------------------------------------------------------------------- main
def _write(name: str, data: dict) -> None:
    OUT().mkdir(parents=True, exist_ok=True)
    (OUT() / f"{name}.json").write_text(json.dumps(data, indent=1, default=float), encoding="utf-8")
    rows = [r for r in data.get("rows", []) if isinstance(r, dict)]
    if rows:
        keys = sorted({k for r in rows for k in r})
        with open(OUT() / f"{name}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows([{k: r.get(k, "") for k in keys} for r in rows])


def _print_rows(rows: list, digits: int = 3) -> None:
    for r in rows:
        print("  ", {k: (round(v, digits) if isinstance(v, float) else v) for k, v in r.items()})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").strip().splitlines()[0])
    ap.add_argument("--only", choices=["femm_cal", "ccx_cal", "sector", "sector_mesh", "sector_split"], action="append")
    ap.add_argument("--reevaluate", action="store_true", help="apply the sector tolerances to the saved fem_sector.json again (no solver)")
    ap.add_argument("--rederive", action="store_true",
                    help="recompute the NETWORK's side of the saved slice results (fem_sector, fem_sector_split, ccx_sector) with the present thermal.py, keeping the solvers' numbers (no solver)")
    args = ap.parse_args(argv)
    if args.rederive:
        ok_all = True
        for name, split in (("fem_sector", False), ("fem_sector_split", True), ("ccx_sector", True)):
            if (OUT() / f"{name}.json").exists():
                res = rederive_network(name, split)
                ok_all &= res["ok"]
                print(f"{name}: network side recomputed:", "PASS" if res["ok"] else "FAIL")
                _print_rows(res["rows"])
        return 0 if ok_all else 1
    if args.reevaluate:
        saved = json.loads((OUT() / "fem_sector.json").read_text(encoding="utf-8"))
        for row in saved["rows"]:
            row.setdefault("has_copper", row["variant"] != "iron_only")
            row.setdefault("h_outer", SECTOR_VARIANTS[row["variant"]]["h"])
        res = evaluate_sector_rows(saved["rows"])
        _write("fem_sector", res)
        print("slice against network (re-evaluated):", "PASS" if res["ok"] else "FAIL")
        _print_rows(res["rows"])
        return 0 if res["ok"] else 1
    todo = args.only or ["femm_cal", "ccx_cal", "sector", "sector_mesh", "sector_split"]
    ok = True
    p = thermal.design_params()
    if any(t in todo for t in ("femm_cal", "sector", "sector_mesh", "sector_split")):
        # ONE pool session for all the FEMM work: the solver slot is taken once, so other jobs waiting for the slot run in between sessions
        with femm_queue.femm_pool("thermal_fem: FEMM heat flow (calibration, slice)", workers=POOL_WORKERS, task_timeout_s=TASK_TIMEOUT_S) as pool:
            if "femm_cal" in todo:
                res = calibrate_femm(pool)
                _write("fem_calibration_femm", res)
                ok &= res["ok"]
                print("FEMM calibration:", "PASS" if res["ok"] else "FAIL")
                _print_rows(res["rows"], 4)
            if "sector" in todo:
                res = compare_sector(pool, p)
                _write("fem_sector", res)
                ok &= res["ok"]
                print("slice against network:", "PASS" if res["ok"] else "FAIL")
                _print_rows(res["rows"])
            if "sector_split" in todo:
                res = compare_sector(pool, p, split_teeth=True)
                _write("fem_sector_split", res)
                print("slice against network, teeth in three blocks (like-for-like diagnostics):", "PASS" if res["ok"] else "FAIL")
                _print_rows(res["rows"])
            if "sector_mesh" in todo:
                rows = []
                for scale in (2.0, 1.0, 0.5):
                    rows += compare_sector(pool, p, {"base": SECTOR_VARIANTS["base"]}, mesh_scale=scale)["rows"]
                _write("fem_sector_mesh", {"rows": rows})
                print("slice mesh pair:")
                _print_rows(rows)
    if "ccx_cal" in todo:
        res = calibrate_ccx()
        _write("fem_calibration_ccx", res)
        ok &= res["ok"]
        print("CalculiX calibration:", "PASS" if res["ok"] else "FAIL")
        _print_rows(res["rows"], 4)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
