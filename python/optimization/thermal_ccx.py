"""optimization.thermal_ccx: CalculiX 2.23 solves the same slot-pitch slice as FEMM did, so that two independent solvers can be held against the thermal network.

What it is, and why
    thermal_fem.py checked the network against FEMM's heat-flow solver on a slice of the stator (one slot pitch). The project's plan (FEA.md section 3) names
    CalculiX as the thermal checker and FEMM as the second opinion; this file runs the same slice in CalculiX, so that the hot spot of the slot is agreed by two
    programs that share nothing but the geometry and the loads. It reuses thermal_geometry.build_slice (the regions, in millimetres) and thermal_fem's variants,
    network and verdict, and touches no FEMM.

Method
    1. The regions (winding, liner, teeth in three blocks, yoke, contact layer, housing) are polygons of 0.5 degree segments. Gmsh's OpenCASCADE kernel makes each one a
       surface and `fragment`s them, which makes their shared walls one set of curves (a conformal mesh); a Distance and a Threshold field refine the mesh around the
       liner and the contact layer (0.1 mm at the walls, 0.6 mm in the bulk); second-order triangles (6 nodes). The Python gmsh module is run inside intact_path()
       (it cuts the PATH of the process, checks/__init__.py).
    2. The deck is plain text in units of mm, tonne, s: conductivity in mW/(mm K) is numerically W/(m K); the heat generation in mW/mm^3 is 1e-6 x W/m^3; the film
       coefficient 1e-3 x W/(m^2 K); thickness 1000 mm (one metre of stack). The elements are CPS6 (plane-stress type, used as 2D heat elements by a heat transfer step:
       thermal_fem.py's note says why DC2D6 and DC2D8 are not used). Every region is an element set with its material and its *DFLUX body heat source; the outer edge of
       the housing carries *FILM with the variant's h and a sink at 0 (so the temperatures are rises above the ambient). The cut lines through the tooth centres carry
       nothing: untouched edges are adiabatic, which is the symmetry of the slice.
    3. ccx_static.exe runs through femm_queue.run_solver (below normal priority, 2 threads, a timeout, one at a time on the PC). The temperatures come back from the .frd.
    4. The region means are integrals over the curved 6-node triangles (a 6-point Gauss rule of degree 4: exact for what a quadratic element can hold), the hot spot is the
       largest nodal value inside the winding, the outer surface and the bore are the means over the nodes on them.

Run (from the Python to FEMM folder; one ccx solve per variant, about a second each, each under the solver slot)
    uv run python -m optimization.thermal_ccx                  every variant of thermal_fem.SECTOR_VARIANTS
    uv run python -m optimization.thermal_ccx --only base --only iron_only
Writes out\\optimization\\thermal\\ccx_sector.json and .csv; each variant's deck, result and log are in out\\optimization\\thermal\\ccx_sector\\<variant>\\.
Exit status 0 if the stated tolerances (thermal_fem.py) hold for the CalculiX numbers against the network.

Numbers to expect (this PC, 2026-11-09, CalculiX 2.23, Gmsh 4.15.2): the mesh has 29,319 nodes and 14,554 six-node triangles (15,711 with the 0.40 mm liner), made in about 2 s; a ccx
    solve takes about 2.2 s. CalculiX and FEMM agree on every variant: the winding's hot spot to 0.008 %, its mean to 0.033 %, the teeth's mean to 0.009 %, the yoke's to 0.0003 %
    (for the base variant, 47.653 / 46.459 / 44.241 / 42.790 C above the ambient in CalculiX against 47.652 / 46.457 / 44.240 / 42.790 in FEMM). So the network's errors against
    CalculiX are the errors against FEMM (base: hot spot +0.007 %, mean +0.70 %; the same 8 of 9 variants inside the tolerance, and the iron-only variant failing the same
    T_teeth criterion, 18.3 %, for the same reason: see thermal_fem.py).

Self-check (the answers follow).
  1. Why is the winding's mean from an integral over the elements, and not the mean of the nodal values?
  2. The two solvers agree to 0.1 % and the network is 2 % away. Whose number is wrong?
Answers.
  1. The nodes crowd where the mesh is fine (the liner, the walls), so a mean of nodal values weights those places more than their area deserves. The integral weights every
     bit of area equally.
  2. Neither, necessarily: the finite-element models solve the geometry they were given exactly, and the network is an approximation of it (an equivalent rectangle, a one-dimensional
     fin). A 2 % difference is the size of those approximations, which is what the check is for.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

import femm_queue
import paths
from checks import intact_path
from optimization import thermal, thermal_fem, thermal_geometry

# the 6-point Gauss rule of degree 4 on a triangle (Dunavant): exact for a quadratic function over a 6-node triangle
_A1, _A2 = 0.445948490915965, 0.091576213509771
_W1, _W2 = 0.223381589678011, 0.109951743655322
GAUSS = [((_A1, _A1, 1 - 2 * _A1), _W1), ((_A1, 1 - 2 * _A1, _A1), _W1), ((1 - 2 * _A1, _A1, _A1), _W1),
         ((_A2, _A2, 1 - 2 * _A2), _W2), ((_A2, 1 - 2 * _A2, _A2), _W2), ((1 - 2 * _A2, _A2, _A2), _W2)]


def shape6(l1: float, l2: float, l3: float):
    """The six shape functions of a quadratic triangle at the area coordinates (l1, l2, l3) and their derivatives with respect to (l1, l2)."""
    n = np.array([l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1), 4 * l1 * l2, 4 * l2 * l3, 4 * l3 * l1])
    d1 = np.array([4 * l1 - 1, 0.0, -(4 * l3 - 1), 4 * l2, -4 * l2, 4 * l3 - 4 * l1])
    d2 = np.array([0.0, 4 * l2 - 1, -(4 * l3 - 1), 4 * l1, 4 * l3 - 4 * l2, -4 * l1])
    return n, d1, d2


RULE = [(shape6(*point), weight) for point, weight in GAUSS]


def tri6_area_and_integral(xy: np.ndarray, values: np.ndarray) -> tuple[float, float]:
    """(area, the integral of the nodal field over the area) of a curved 6-node triangle: xy (6, 2), values (6,)."""
    area = integral = 0.0
    for (n, d1, d2), w in RULE:
        j = 0.5 * w * ((d1 @ xy[:, 0]) * (d2 @ xy[:, 1]) - (d1 @ xy[:, 1]) * (d2 @ xy[:, 0]))
        area += j
        integral += j * (n @ values)
    return area, integral


def mesh_slice(regions: dict, outer_r: float, bore_r: float, h_fine: float = 0.1, h_coarse: float = 0.6) -> dict:
    """Mesh the regions (shapely polygons, mm) in second-order triangles. Returns {"nodes": {id: (x, y)}, "tris": {region: [(6 node ids), ...]}, "outer": [(a, b, mid)],
    "bore": [node ids on the bore]}."""
    import gmsh
    with intact_path():
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("slice")
            occ = gmsh.model.occ
            points: dict = {}

            def point(xy):
                key = (round(xy[0], 6), round(xy[1], 6))
                if key not in points:
                    points[key] = occ.addPoint(key[0], key[1], 0)
                return points[key]

            order, surfaces = [], []
            for name, poly in regions.items():
                for part in getattr(poly, "geoms", [poly]):
                    ring = list(part.exterior.coords)[:-1]
                    tags = [point(c) for c in ring]
                    lines = [occ.addLine(a, b) for a, b in zip(tags, tags[1:] + tags[:1]) if a != b]
                    surfaces.append(occ.addPlaneSurface([occ.addCurveLoop(lines)]))
                    order.append(name)
            _, owners = occ.fragment([(2, surfaces[0])], [(2, s) for s in surfaces[1:]])
            owner = {}
            for name, dimtags in zip(order, owners):
                for _, tag in dimtags:
                    owner[tag] = name
            occ.synchronize()
            groups: dict = {}
            for tag, name in owner.items():
                groups.setdefault(name, []).append(tag)
            # fine near the thin layers: a Distance field from the liner's and the contact layer's boundary curves
            thin = [name for name in groups if name in ("liner", "contact")]
            curves = sorted({abs(c) for name in thin for _, c in gmsh.model.getBoundary([(2, t) for t in groups[name]], oriented=False)})
            field = gmsh.model.mesh.field
            distance = field.add("Distance")
            field.setNumbers(distance, "CurvesList", curves)
            field.setNumber(distance, "Sampling", 40)
            threshold = field.add("Threshold")
            field.setNumber(threshold, "InField", distance)
            field.setNumber(threshold, "SizeMin", h_fine)
            field.setNumber(threshold, "SizeMax", h_coarse)
            field.setNumber(threshold, "DistMin", 0.15)
            field.setNumber(threshold, "DistMax", 2.0)
            field.setAsBackgroundMesh(threshold)
            for option in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature", "Mesh.MeshSizeExtendFromBoundary"):
                gmsh.option.setNumber(option, 0)
            gmsh.option.setNumber("Mesh.ElementOrder", 2)
            gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
            for name, tags in groups.items():
                gmsh.model.addPhysicalGroup(2, tags, name=name)
            gmsh.model.mesh.generate(2)
            tags, coords, _ = gmsh.model.mesh.getNodes()
            xyz = np.asarray(coords).reshape(-1, 3)
            nodes = {int(t): (float(x), float(y)) for t, (x, y, _) in zip(tags, xyz)}
            tris = {}
            for name, surface_tags in groups.items():
                rows = []
                for tag in surface_tags:
                    types, _, connectivity = gmsh.model.mesh.getElements(2, tag)
                    for kind, conn in zip(types, connectivity):
                        assert kind == 9, f"expected 6-node triangles (type 9), got {kind}"
                        for row in np.asarray(conn).reshape(-1, 6):
                            t = tuple(int(n) for n in row)
                            (x1, y1), (x2, y2), (x3, y3) = nodes[t[0]], nodes[t[1]], nodes[t[2]]
                            if (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) < 0:        # make it counter-clockwise: corners 2 and 3 (and their mid nodes) exchanged
                                t = (t[0], t[2], t[1], t[5], t[4], t[3])
                            rows.append(t)
                tris[name] = rows
            outer, bore = [], set()
            for _, c in gmsh.model.getEntities(1):
                types, _, connectivity = gmsh.model.mesh.getElements(1, c)
                for kind, conn in zip(types, connectivity):
                    if kind != 8:
                        continue
                    for a, b, m in np.asarray(conn).reshape(-1, 3):
                        a, b, m = int(a), int(b), int(m)
                        if all(abs(math.hypot(*nodes[n]) - outer_r) < 0.02 for n in (a, b, m)):
                            outer.append((a, b, m))
                        for n in (a, b, m):
                            if abs(math.hypot(*nodes[n]) - bore_r) < 1e-3:
                                bore.add(n)
        finally:
            gmsh.finalize()
    used = {n for rows in tris.values() for t in rows for n in t}
    return {"nodes": {n: nodes[n] for n in sorted(used)}, "tris": tris, "outer": outer, "bore": sorted(bore & used)}


def write_deck(path: Path, mesh: dict, mats: dict, h: float) -> dict:
    """The deck: nodes, CPS6 elements per region (an element set each), materials, body heat sources, film on the outer edge. mats: {region: (k W/(m K), q W/m^3)}.
    Returns {"element_ids": {region: [ids]}, "elements": {id: nodes}}."""
    nodes = mesh["nodes"]
    lines = ["*HEADING", "slot pitch of the stator, CPS6 used as a 2D heat element", "*NODE"]
    for n, (x, y) in nodes.items():
        lines.append(f"{n}, {thermal_fem.ccx_number(x)}, {thermal_fem.ccx_number(y)}, 0.0")
    eid, ids, elements = 0, {}, {}
    lines.append("*ELEMENT, TYPE=CPS6, ELSET=EALL")
    for name, rows in mesh["tris"].items():
        for t in rows:
            eid += 1
            ids.setdefault(name, []).append(eid)
            elements[eid] = t
            lines.append(f"{eid}, " + ", ".join(str(n) for n in t))
    for name, items in ids.items():
        lines.append(f"*ELSET, ELSET=R_{name}")
        lines += [", ".join(str(x) for x in items[k:k + 12]) for k in range(0, len(items), 12)]
    edge_face = {}
    for e, t in elements.items():
        for face, (a, b) in enumerate(((t[0], t[1]), (t[1], t[2]), (t[2], t[0])), start=1):
            edge_face[frozenset((a, b))] = (e, face)
    for name, (k, q) in mats.items():
        lines += [f"*MATERIAL, NAME=M_{name}", "*CONDUCTIVITY", thermal_fem.ccx_number(k), f"*SOLID SECTION, ELSET=R_{name}, MATERIAL=M_{name}", "1000."]
    lines += ["*STEP", "*HEAT TRANSFER, STEADY STATE", "1., 1."]
    for name, (k, q) in mats.items():
        if q:
            lines += ["*DFLUX", f"R_{name}, BF, {thermal_fem.ccx_number(q * 1e-6)}"]
    lines.append("*FILM")
    for a, b, m in mesh["outer"]:
        e, face = edge_face[frozenset((a, b))]
        lines.append(f"{e}, F{face}, 0, {thermal_fem.ccx_number(h * 1e-3)}")
    lines += ["*NODE FILE", "NT", "*END STEP"]
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    return {"element_ids": ids, "elements": elements}


def run_variant(name: str, variant: dict, p, folder: Path, mesh_cache: dict) -> dict:
    """Mesh (once), write the deck, solve with ccx, and return the same numbers as thermal_fem's slice rows."""
    split = True
    task, net = thermal_fem.build_sector_task(variant, p, 1.0, folder, f"ccx_{name}", split_teeth=split)
    a = thermal.Assumptions(**variant.get("assume", {}))
    key = (a.liner_mm, round(net["model"].r_ho, 6))
    if key not in mesh_cache:
        from shapely.geometry import MultiPolygon, Polygon
        regions = {}
        for rname, d in task["regions"].items():                       # the vertices FEMM was given, so that both solve the same drawing
            parts = [Polygon(ring) for ring in d["parts"]]
            regions[rname] = parts[0] if len(parts) == 1 else MultiPolygon(parts)
        t0 = time.perf_counter()
        mesh_cache[key] = mesh_slice(regions, net["model"].r_ho * 1e3, p.r_si)
        mesh_cache[key]["mesh_seconds"] = time.perf_counter() - t0
    mesh = mesh_cache[key]
    job = Path(folder) / name
    job.mkdir(parents=True, exist_ok=True)
    mats = {r: (d["k"], d["q"]) for r, d in task["regions"].items()}
    deck = write_deck(job / "slice.inp", mesh, mats, variant["h"])
    run = femm_queue.run_solver([paths.tool("ccx"), "-i", "slice"], label=f"thermal_ccx {name}", cwd=job, timeout_s=300, threads=2)
    if run.timed_out or not (job / "slice.frd").exists():
        return {"variant": name, "status": "error", "returncode": run.returncode, "stdout": run.stdout[-600:]}
    temps = thermal_fem.read_frd_temperature(job / "slice.frd")
    nodes = mesh["nodes"]
    out = {"variant": name, "status": "ok", "nodes": len(nodes), "elements": len(deck["elements"]), "seconds": run.seconds,
           "mesh_seconds": mesh.get("mesh_seconds", 0.0)}
    means, areas, maxima = {}, {}, {}
    for region, items in deck["element_ids"].items():
        a_tot = i_tot = 0.0
        node_set = set()
        for e in items:
            t = deck["elements"][e]
            xy = np.array([nodes[n] for n in t])
            v = np.array([temps[n] for n in t])
            ar, it = tri6_area_and_integral(xy, v)
            a_tot += ar
            i_tot += it
            node_set.update(t)
        means[region], areas[region] = i_tot / a_tot, a_tot
        maxima[region] = max(temps[n] for n in node_set)
    out["mean"], out["area_mm2"], out["max"] = means, areas, maxima
    # the outer surface: the mean over the boundary edges (Simpson along each 3-node edge), and the bore: the mean over its nodes
    length = total = 0.0
    for a_, b_, m_ in mesh["outer"]:
        seg = math.dist(nodes[a_], nodes[b_])
        length += seg
        total += seg * (temps[a_] + 4 * temps[m_] + temps[b_]) / 6.0
    out["T_outer"] = total / length
    out["T_bore"] = float(np.mean([temps[n] for n in mesh["bore"]]))
    return out


def compare(rows_ccx: list[dict], p) -> dict:
    """The same comparison as thermal_fem.compare_sector, with CalculiX's numbers in place of FEMM's, plus the difference between the two solvers."""
    femm_split = {r["variant"]: r for r in json.loads((thermal_fem.OUT() / "fem_sector_split.json").read_text(encoding="utf-8"))["rows"]}
    rows = []
    for r in rows_ccx:
        if r.get("status") != "ok":
            rows.append(r)
            continue
        v = thermal_fem.SECTOR_VARIANTS[r["variant"]]
        task, net_inputs = thermal_fem.build_sector_task(v, p, 1.0, thermal_fem.OUT() / "ccx_sector", "x", split_teeth=True)
        net = thermal_fem.network_for_variant(v, net_inputs, p)
        ar = net_inputs["areas_m2"]
        a_s, a_f, a_r = ar["teeth_shoe"], ar["teeth_fin"], ar["teeth_root"]
        t_teeth = (r["mean"]["teeth_shoe"] * a_s + r["mean"]["teeth_fin"] * a_f + r["mean"]["teeth_root"] * a_r) / (a_s + a_f + a_r)
        fe = {"T_hot": r["max"]["winding"], "T_mean": r["mean"]["winding"], "T_teeth": t_teeth, "T_yoke": r["mean"]["yoke"], "T_fin": r["mean"]["teeth_fin"],
              "T_tip": r["T_bore"], "T_outer": r["T_outer"]}
        nw = {"T_hot": net["T_hot"], "T_mean": net["T_W"], "T_teeth": net["T_T"], "T_yoke": net["T_Y"], "T_fin": net["T_T"], "T_tip": net["T_tip"], "T_outer": net["T_outer"]}
        row = {"variant": r["variant"], "h_outer": v["h"], "has_copper": v["p_slot"] > 0, "nodes": r["nodes"], "elements": r["elements"], "seconds": r["seconds"]}
        fm = femm_split.get(r["variant"], {})
        for key in ("T_hot", "T_mean", "T_teeth", "T_yoke", "T_fin", "T_tip"):
            d_fe, d_nw = fe[key] - fe["T_outer"], nw[key] - nw["T_outer"]
            row[f"{key}_ccx"], row[f"{key}_net"] = fe[key], nw[key]
            row[f"{key}_drop_ccx"], row[f"{key}_drop_net"] = d_fe, d_nw
            row[f"{key}_drop_err_%"] = 100.0 * (d_nw - d_fe) / d_fe if abs(d_fe) > 1e-9 else float("nan")
            row[f"{key}_amb_err_%"] = 100.0 * (nw[key] - fe[key]) / fe[key] if abs(fe[key]) > 1e-9 else float("nan")
            if f"{key}_fe" in fm:
                row[f"{key}_femm"] = fm[f"{key}_fe"]
                row[f"{key}_ccx_vs_femm_%"] = 100.0 * (fe[key] - fm[f"{key}_fe"]) / fm[f"{key}_fe"] if abs(fm[f"{key}_fe"]) > 1e-9 else float("nan")
        row["T_outer_ccx"], row["T_outer_net"] = fe["T_outer"], nw["T_outer"]
        rows.append(row)
    return thermal_fem.evaluate_sector_rows(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").strip().splitlines()[0])
    ap.add_argument("--only", action="append", choices=list(thermal_fem.SECTOR_VARIANTS))
    args = ap.parse_args(argv)
    p = thermal.design_params()
    folder = thermal_fem.OUT() / "ccx_sector"
    folder.mkdir(parents=True, exist_ok=True)
    names = args.only or list(thermal_fem.SECTOR_VARIANTS)
    cache: dict = {}
    rows = [run_variant(n, thermal_fem.SECTOR_VARIANTS[n], p, folder, cache) for n in names]
    res = compare(rows, p)
    thermal_fem._write("ccx_sector", res)
    print("slice, CalculiX against the network:", "PASS" if res["ok"] else "FAIL")
    thermal_fem._print_rows(res["rows"])
    return 0 if res["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
