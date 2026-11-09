r"""optimization.compare_designs: one figure of the starting design and the optimized design side by side, for a README: what the optimization changed, and what
the field does at the study's comparison point.

What it is, and why
    A reader who has not seen the study needs one picture that says what changed and what it bought. The study compared two motors at the same copper loss, each
    at its own maximum-torque-per-ampere (MTPA) current angle, on the fine mesh: the STARTING design (design.params({}): the band-free rotor with the original stator)
    and the OPTIMIZED design (design.params(final.json's "overrides")). This module draws both from the parameters the study scored, solves each at its comparison point,
    and holds every number it prints against the study's own record (final.json) and against the values the study's write-up quotes. The model is FEMM 5.0, 2D, a
    20 mm slice, the fine mesh (smart mesh on), the 10-turn placeholder winding; currents are A peak.

The figure (reference_vs_optimized.png, 1800 px wide: two columns, "Starting design" and "Optimized design")
    top row      the geometry of one 90-degree sector (two of the eight poles; the motor repeats four times), both at the same scale: laminated steel in light gray, air white,
                 the N45 magnets red with an arrow along their magnetisation, the coil layers tinted by phase (A blue, B orange, C aqua), the shaft outline. A bar on the
                 drawing marks where each changed dimension is measured, and a leader brings its value out to a label: the centre rib, the outer bridge, the V angle (the angle
                 between the two magnet axes of a pole), the slot opening, the narrowest tooth and the yoke. On the optimized design the starting design's magnets are a dashed outline.
    bottom row   the magnetic field at the comparison point: |B| as a colour on ONE scale for both panels (0 to 2.2 T, with a colour bar), the flux lines over it (contours of the
                 vector potential A; one line for each 20 uWb through the 20 mm slice, the same step in both), and beside each panel the current and its angle, the mean
                 torque and the ripple, with the change against the starting design on the optimized one.
    The title and a one-line caption name the basis; four short lines at the foot say what the numbers are. The text of the figure is checked for overlaps and for leaving the
    page (Canvas.overlaps), and the drawn leaders are checked against the magnet slots; either problem is printed as a WARNING and makes the exit status 1.
Other outputs, in out\optimization\compare\ (git-ignored)
    reference_vs_optimized_geometry.png   the top row alone (with its own legend and notes)
    reference_vs_optimized_wide.png       a single wide row (2000 px): each design's geometry and field side by side, the numbers under them, the changes as one line
    reference_vs_optimized_numbers.json   every number printed in the figures, where it comes from, this script's own solves, and every check with its result
    field_reference.npz, field_optimized.npz   the solved fields (mesh, A, triangles, materials), so that the pictures can be drawn again without FEMM (--reuse)

The method, step by step
    1. The geometry. design.params(overrides) is a motorgeo.Params; motorgeo.build gives the regions of one sector and motorgeo.visible paints later regions over earlier
       ones (the magnets over the rotor steel, the slots over the stator steel); each visible region is drawn as a polygon with its holes. Every dimension in the figure is
       computed from these, never typed: the outer bridge by design.bridge_of (the rotor's radius less design.slot_reach, the slot's farthest point), the centre rib, the
       narrowest tooth and the yoke by motorgeo.check (distances between the drawn polygons), the V angle as 2 v_half_deg, the slot opening as slot_open. Each bar is found on the
       same polygons (the nearest points of two shapes, the farthest vertex of a slot) and its length is checked against the number it labels. The magnets are drawn as the
       model has them: each magnet slot (9.28 x 1.77 mm) is filled with N45.
    2. The solves. One task per design runs in a femm_queue.femm_pool worker (the solver slot: the run waits while another run holds it; 2 workers at below-normal
       priority, a task timeout of design.TASK_TIMEOUT_S; the pool starts each worker's own FEMM, and this module never starts or ends one). The task writes the fine-mesh
       model under a name of its own, removes it afterwards, and solves with motorgeo.FemmModel: (a) the no-load flux linkage at the two rotor positions of
       design.ANGLE_POSITIONS (the magnets' flux linkage, psi_pm); (b) the torque at the 15 rotor positions of design.RIPPLE_POSITIONS at the study's current and angle,
       which is design.score's own scoring step, so the mean torque and the ripple come out the study's; (c) one more position, whose solution file is read back with
       femfile.read_ans: the nodes, the vector potential A at each node, the triangles and their materials. That position is the electrical angle at which the rotor stands
       where it was drawn (view_angle: FemmModel.solve turns the rotor side of the air-gap element by theta_e / pole_pairs + axis_A - theta_d0, and this makes it 0), so the picture
       needs no turned rotor mesh. It is the same electrical angle in both designs.
    3. A to flux lines. The mesh is first-order triangles, so A is linear inside each one, and tricontour draws the lines of the solution exactly. A line of constant A is a flux
       line (the flux density is tangent to it); the flux between two lines is their difference in A times the stack length, so a fixed step in A is a fixed amount of flux: the
       smallest step of FLUX_STEPS_UWB that leaves 14 lines or fewer each side of zero on the field with the larger A (20 uWb here), the same in both panels, so that the line
       density can be compared and the larger flux of the optimized design shows.
    4. B per triangle to the colour map. B = (dA/dy, -dA/dx) is constant on a first-order triangle (femfile.element_b), so every triangle is one flat colour of |B|, and the
       steel's edges stay sharp (|B| jumps across a material boundary). The colour map (b_cmap) is one warm ramp that darkens steadily from warm white to deep wine, with its stops
       closer together between 0.3 and 1.5 T, where the steel works; values above 2.2 T take the last colour (the colour bar's arrow). FEMM leaves the middle third of the air
       gap (its air-gap element) unmeshed: the two arcs of nodes that bound it are joined by a strip of triangles (build_field), so that colour and flux lines run through the gap.
    5. The numbers. The current, angle, mean torque and ripple printed in the figure are the study's fine-mesh numbers (final.json, "compare"), not re-derived. The solves of step 2
       repeat the study's scoring and are held against them: the mean torque within 1 %, the ripple within 0.2 points, psi_pm within 0.01 mWb; they agree to rounding, because the code,
       the mesh and the current are the same. The dimensions are held against final.json (the optimized design's bridge, V angle and slot opening) and against the values the
       study quotes (bridge and rib within 0.01 mm, the V angle within 0.1 degree). A row of PASS or FAIL lines is printed, and the exit status is 1 if any fails.

Run (from the Python to FEMM folder)
    uv run python -m optimization.compare_designs            solve through the solver slot, then draw
    uv run python -m optimization.compare_designs --reuse    draw again from the saved fields: no FEMM, no solver slot
    options: --workers 2   (FEMM workers; one per design is all there is to run)

Numbers to expect (FEMM 5.0, this PC, 2026-11-09)
    two tasks in parallel, 10 and 13 s each, about 23 s of wall time in all (the three PNGs take most of the rest); the mesh 15,859 nodes and 30,610 triangles (starting) and 15,693 and 30,290 (optimized), plus
    530 strip triangles in the air-gap element. Mean torque 1.1288 and 1.8073 N.m, ripple 18.45 and 4.84 %, psi_pm 16.64 and 24.20 mWb, currents 10.0000 and 11.8627 A at
    22.46 and 22.47 degrees (all equal to the study's); the torque at the picture's position (theta_e = -170.1 = 189.9 degrees) 1.086 and 1.780 N.m. |A| reaches 185 and 228 uWb
    (through the slice), so 9 and 11 lines each side of zero at 20 uWb. The peak |B| in a few triangles at the slot corners is 2.85 and 3.18 T; 0.26 % and 0.37 % of the steel's
    area is above the 2.2 T scale. The PNGs: 1800 x 1778 px (about 1 MB), 1800 x 1067 px (0.3 MB) and 2000 x 1175 px (0.7 MB); with --reuse the run takes about 8 s. 53 checks are
    printed, all PASS, exit status 0.

Questions (the answers follow).
  1. Why is the field drawn at the electrical angle -170.1 degrees and not at 0?
  2. Why is one flux step used for both panels, and what would a step chosen for each panel hide?
  3. The colour bar ends in an arrow. What does it say, and how much of the picture does it concern?
Answers.
  1. At 0 degrees FemmModel.solve turns the rotor side of the air-gap element by 42.525 degrees (the axis of phase A is at 65.025 degrees and the first d-axis at 22.5), so the
     rotor's mesh would stand 42.525 degrees from the place it was drawn and the picture would have to turn it, and cut it at the sector's edges. At -170.1 degrees (equal to
     189.9 degrees) the turn is 0 degrees, and rotor and stator stand as drawn: one picture of both meshes, with the geometry row directly above it.
  2. A flux line is a step of the flux. With a step for each panel both would show about the same number of lines, and the optimized design's larger flux (|A| reaches 228
     uWb against 185) would not be seen; with one step the lines are as dense as the flux is.
  3. That every |B| above 2.2 T has the colour of 2.2 T, so the scale does not say how much above. It concerns 0.26 % and 0.37 % of the steel's area, a few triangles at the magnet
     slots' corners, where the solver's flux density is highest and least meaningful (a corner is a singularity of the field, and the mesh there is coarse).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
import traceback
import uuid
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import femfile
import femm_queue
import motorgeo as mg
import paths
from optimization import design

# ============================================================================================== where things are
OUT = paths.path("out") / "optimization" / "compare"
FINAL_JSON = paths.path("out") / "optimization" / "final" / "final.json"
PNG_MAIN = "reference_vs_optimized.png"
PNG_GEOMETRY = "reference_vs_optimized_geometry.png"
PNG_WIDE = "reference_vs_optimized_wide.png"
JSON_NUMBERS = "reference_vs_optimized_numbers.json"

# ============================================================================================== the look
SURFACE = "#ffffff"                  # the README page is white
INK, INK2, MUTED, AXIS = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"   # the dataviz reference ink tokens: primary, secondary, muted, hairline
STEEL = "#d9dde2"                    # laminated steel: light gray
EDGE = "#4a4f57"                     # the outline of every body
MAGNET = "#d33a39"                   # N45 magnets (the arrow shows the direction)
# the three phases: categorical slots 1 - 3 of the reference palette (blue, orange, aqua) taken 80 % toward the page white; as a set they
# pass the palette validator on a white surface (worst colour-blind distance 8.2, worst normal-vision distance 19.0)
PHASE = {"A": "#5593de", "B": "#ef865d", "C": "#49bf95"}
FLUX_LINE = "#14181d"                # flux lines
CMAP_NAME = "b_heat"                 # |B|: one warm sequential ramp, warm white (air) to deep wine (saturated steel); see b_cmap
FONT = ["Segoe UI", "DejaVu Sans"]
DPI = 150
W_PX = 1800                          # the figure's width in pixels (about 900 px when GitHub shrinks it)
HEAD_PT, VAL_PT = 10.5, 13.5         # a callout: the feature's name, then its value
TITLE_PT, CAPTION_PT, COLUMN_PT, LEGEND_PT, NOTE_PT = 19.0, 11.5, 17.0, 10.5, 10.0

# ============================================================================================== the study's operating point
B_MAX = 2.2                          # T, the top of the shared |B| scale
STACK_M = 0.02                       # the model is a 20 mm slice
POLE0, POLE1 = 0, 1                  # the two poles of the sector: 22.5 and 67.5 degrees (their magnets point out and in)
FLUX_STEPS_UWB = (10, 20, 25, 40, 50, 80, 100, 125, 200, 250, 400, 500)   # candidates for the flux between two neighbouring lines [uWb]

# ============================================================================================== where the labels go (mm, in the drawing)
XLIM = (-1.5, 96.0)                  # the drawing's window: the quarter disc [0, 64] x [0, 64] and a margin for labels on the right
YLIM_GEOMETRY = (-13.0, 66.0)        # the geometry panels have a margin below for the rotor's labels
YLIM_FIELD = (-3.0, 66.0)
SHELF_X = 69.0                       # a stator leader runs out along its radius to R_ELBOW, then level to this x; its label starts 1 mm further
R_ELBOW = 66.0
STATS_X = 65.0                       # the left edge of the numbers beside a field panel
LABEL_Y = -2.6                       # the rotor's leaders end on this line; their labels hang below it
RIB_DOGLEG_MM = 12.0                 # the rib's leader runs this far along the pole's d-axis toward the shaft, then straight down
V_LABEL_R = 11.0                     # the V's label sits this far from the apex, in the direction of the shaft, and its arc has radius V_ARC_R
V_ARC_R = 6.0
V_LINE_L = 11.0                      # the V's axes run this far toward the shaft from the apex


# ============================================================================================== the two designs
@dataclass
class Design:
    key: str                         # "reference" or "optimized"
    title: str                       # the column header
    overrides: dict                  # motorgeo.Params overrides (design.params)
    study: dict                      # the study's fine-mesh row for it (final.json "compare")
    p: mg.Params = None
    regions: list = field(default_factory=list)
    visible: dict = field(default_factory=dict)
    metrics: dict = field(default_factory=dict)
    dims: dict = field(default_factory=dict)
    marks: dict = field(default_factory=dict)

    def __post_init__(self):
        self.p = design.params(self.overrides)
        self.regions = mg.build(self.p)
        self.visible = mg.visible(self.regions)
        self.metrics, _ = mg.check(self.p, self.regions)
        self.dims = dimensions(self)
        self.marks = markers(self)


def dimensions(d: Design) -> dict:
    """Every dimension the figure prints, from the parameters with the project's own helpers (never typed in)."""
    p, m = d.p, d.metrics
    return {
        "bridge_mm": design.bridge_of(p),                 # the rotor's radius less the magnet slot's reach
        "slot_reach_mm": design.slot_reach(p),
        "bridge_check_mm": m["outer_bridge"],             # the same, read off the drawn polygons by motorgeo.check
        "rib_mm": m["centre_rib"],                        # the steel between the two slots of a pole, at its thinnest
        "v_angle_deg": 2.0 * p.v_half_deg,                # the angle between the two magnets of a pole
        "slot_open_mm": p.slot_open,
        "tooth_mm": m["tooth_min"],                       # the narrowest tooth: under the shoe
        "yoke_mm": m["yoke"],                             # the stator's outer radius less the slots' deepest point
        "h_open_mm": p.h_open,
        "coil_layer_area_mm2": m["coil_layer_area_mm2"],
    }


def _farthest(poly) -> tuple[float, float]:
    xy = np.asarray(poly.exterior.coords)
    return tuple(xy[int(np.argmax(np.hypot(xy[:, 0], xy[:, 1])))])


def markers(d: Design) -> dict:
    """Where each dimension sits on the drawing: {name: {"p0", "p1" (mm, the two ends of its bar), "mm" (its length)}}, found on the drawn polygons."""
    from shapely.ops import nearest_points, unary_union
    p = d.p
    poly = {r.name: r.polygon() for r in d.regions}
    out = {}
    a, b = nearest_points(poly[f"magnet_{POLE0}_up"], poly[f"magnet_{POLE0}_dn"])                 # centre rib: the slots' closest approach
    out["rib"] = {"p0": (a.x, a.y), "p1": (b.x, b.y), "mm": a.distance(b)}
    tip = np.array(_farthest(poly[f"magnet_{POLE0}_dn"]))                                          # outer bridge: the slot's farthest point,
    out["bridge"] = {"p0": tuple(tip), "p1": tuple(tip * p.r_ro / math.hypot(*tip)), "mm": p.r_ro - math.hypot(*tip)}   # straight out to the rim
    slot = [unary_union([poly[f"opening_{k}"], poly[f"coil_{k}_0"], poly[f"coil_{k}_1"]]) for k in range(round(p.slots / p.n_sectors))]
    k_tooth = 2                                                                                    # tooth: between slots 2 and 3, at its narrowest
    a, b = nearest_points(slot[k_tooth], slot[k_tooth + 1])
    out["tooth"] = {"p0": (a.x, a.y), "p1": (b.x, b.y), "mm": a.distance(b), "slots": (k_tooth, k_tooth + 1)}
    k_open = 1                                                                                     # slot opening, across the mouth at the bore
    phi = (k_open + 0.5) * 360 / p.slots
    a_o = math.degrees(math.asin(p.slot_open / 2 / p.r_si))
    q0, q1 = mg.pol(p.r_si, phi - a_o), mg.pol(p.r_si, phi + a_o)
    out["opening"] = {"p0": q0, "p1": q1, "mm": math.dist(q0, q1), "slot": k_open, "angle": phi}
    k_yoke = 4                                                                                     # yoke: from the slot's deepest point to the rim
    deep = np.array(_farthest(slot[k_yoke]))
    out["yoke"] = {"p0": tuple(deep), "p1": tuple(deep * p.r_so / math.hypot(*deep)), "mm": p.r_so - math.hypot(*deep), "slot": k_yoke,
                   "angle": math.degrees(math.atan2(deep[1], deep[0]))}
    b_half = math.radians(p.v_half_deg)                                                            # V angle: pole 1's two magnet axes meet on its d-axis
    theta_d = 180 / p.pole_pairs * (POLE1 + 0.5)
    x_apex = p.mag_x - p.mag_y / math.tan(b_half)
    out["v"] = {"apex": mg.rot((x_apex, 0.0), theta_d), "theta_d": theta_d, "half": p.v_half_deg, "deg": 2 * p.v_half_deg}
    return out


# ============================================================================================== the FEMM side (runs in a femm_pool worker)
def view_angle(p: mg.Params) -> float:
    """The electrical angle [deg] at which the rotor stands exactly where it was drawn: FemmModel.solve turns the rotor side of the air-gap element
    by theta_e / pole_pairs + axis_A - theta_d0, and this is the theta_e that makes that zero. The picture then needs no turned rotor mesh."""
    axes, _ = mg.phase_axes(p)
    return p.pole_pairs * (90.0 / p.pole_pairs - axes["A"])


def field_task(task: dict) -> dict:
    """One design at its comparison point, on the FINE mesh (a femm_pool worker runs this; the model and solution files are its own and are removed).

    task = {"tag", "key", "overrides", "fem_dir", "amps" [A peak], "gamma" [deg from the q-axis toward -d], "theta_view" [electrical deg]}.
    Solves, in this order, on one model: the no-load d-axis flux linkage at the two positions of design.ANGLE_POSITIONS (psi_pm); the torque at the 15
    positions of design.RIPPLE_POSITIONS (the study's own scoring: the same code, the same mesh); and the picture's position, whose solution is read back
    (femfile.read_ans) into the mesh, the nodal A and B of every triangle."""
    import femm
    started = time.perf_counter()
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")
    ans, kept = fem[:-4] + ".ans", fem[:-4] + "_view.ans"
    try:
        p = design.params(task["overrides"])
        regions = mg.build(p)
        mg.to_femm(p, regions, fem, fine=True)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        i_d, i_q = design.currents(float(task["amps"]), float(task["gamma"]))
        no_load = [model.solve(th, 0.0, 0.0)["psi_d"] for th in design.ANGLE_POSITIONS]
        torque = [model.solve(th, i_d, i_q)["T_age"] for th in design.RIPPLE_POSITIONS]
        view = model.solve(float(task["theta_view"]), i_d, i_q, keep_ans=kept)
        sol = femfile.read_ans(kept)
        material = np.array([sol["labels"][k]["material"] for k in sol["label"]])
        result = {"status": "ok", "key": task["key"], "xy": sol["xy"], "A": sol["A"], "tri": sol["tri"].astype(np.int32), "material": material,
                  "torque": [float(t) for t in torque], "psi_pm": float(np.mean(no_load)), "amps": float(task["amps"]), "gamma": float(task["gamma"]),
                  "view": {"theta_e": float(task["theta_view"]), "rotor_deg": float(view["rotor_deg"]), "T_age": float(view["T_age"]),
                           "T_wst": float(view["T_wst"]), "psi_d": float(view["psi_d"]), "psi_q": float(view["psi_q"])},
                  "nodes": int(len(sol["xy"])), "triangles": int(len(sol["tri"])), "seconds": time.perf_counter() - started}
    except Exception as error:                      # one failed solve must not leave files behind or end the pool
        result = {"status": "error", "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4),
                  "seconds": time.perf_counter() - started}
    finally:
        try:
            femm.mi_close()
        except Exception:
            pass
        for path in (fem, ans, kept):
            try:
                os.remove(path)
            except OSError:
                pass
    return result


def solve_designs(designs: list[Design], workers: int = 2) -> dict[str, dict]:
    """Run field_task for every design through femm_queue.femm_pool (the solver slot: this waits while another run holds it), and keep the fields."""
    from optimization import compare_designs as me       # by its module name, so that a worker can import the task function
    fem_dir = OUT / "fem"
    tasks = [{"tag": f"cmp_{d.key}_{uuid.uuid4().hex[:10]}", "key": d.key, "overrides": d.overrides, "fem_dir": str(fem_dir), "amps": d.study["amps"],
              "gamma": d.study["gamma"], "theta_view": view_angle(d.p)} for d in designs]
    with femm_queue.femm_pool("optimization: compare designs (fine mesh)", workers=min(workers, len(tasks)), task_timeout_s=design.TASK_TIMEOUT_S) as pool:
        rows = pool.map(me.field_task, tasks)
    try:
        fem_dir.rmdir()                                   # empty when every task removed its files
    except OSError:
        pass
    results = {}
    for task, row in zip(tasks, rows):
        if not (isinstance(row, dict) and row.get("status") == "ok"):
            raise RuntimeError(f"the solve of {task['key']} failed: {row.get('error') if isinstance(row, dict) else row}")
        results[task["key"]] = row
        save_result(row)
    return results


def save_result(row: dict) -> None:
    """Keep a solved field, so that the pictures can be drawn again without FEMM (--reuse)."""
    OUT.mkdir(parents=True, exist_ok=True)
    scalars = {k: v for k, v in row.items() if k not in ("xy", "A", "tri", "material")}
    np.savez_compressed(OUT / f"field_{row['key']}.npz", xy=row["xy"], A=row["A"], tri=row["tri"], material=row["material"], scalars=json.dumps(scalars))


def load_result(key: str) -> dict:
    with np.load(OUT / f"field_{key}.npz", allow_pickle=False) as z:
        row = json.loads(str(z["scalars"]))
        row.update({"xy": z["xy"], "A": z["A"], "tri": z["tri"], "material": z["material"]})
    return row


# ============================================================================================== the solved field as something to draw
@dataclass
class Field:
    key: str
    tri: object                      # matplotlib.tri.Triangulation: the mesh and the strip that fills the air-gap element
    bmag: np.ndarray                 # |B| of every triangle [T]
    A: np.ndarray                    # nodal A [Wb/m]
    n_mesh: int                      # the first n_mesh triangles are the solver's; the rest fill the gap element
    res: dict                        # the solve's numbers


def _zipper(ang_r, idx_r, ang_s, idx_s) -> np.ndarray:
    """Triangles that fill the strip between two arcs of nodes (both sorted by angle): one new triangle per step along either arc."""
    out, i, j = [], 0, 0
    while i < len(idx_r) - 1 or j < len(idx_s) - 1:
        if j == len(idx_s) - 1 or (i < len(idx_r) - 1 and ang_r[i + 1] <= ang_s[j + 1]):
            out.append((idx_r[i], idx_r[i + 1], idx_s[j]))
            i += 1
        else:
            out.append((idx_r[i], idx_s[j + 1], idx_s[j]))
            j += 1
    return np.array(out, dtype=np.int32)


def build_field(res: dict, p: mg.Params) -> Field:
    """Triangulation + |B| per triangle. FEMM leaves the middle third of the air gap (its air-gap element) unmeshed; two arcs of nodes bound it, so
    the strip between them is filled with triangles of those nodes (nodal A is known there), which lets the flux lines and the colour run through the gap."""
    import matplotlib.tri as mtri
    xy, tri, A = res["xy"], res["tri"], res["A"]
    r = np.hypot(xy[:, 0], xy[:, 1])
    ang = np.degrees(np.arctan2(xy[:, 1], xy[:, 0]))
    r_a, r_b = p.r_ro + p.airgap / 3, p.r_ro + 2 * p.airgap / 3        # motorgeo.build: the unmeshed ring lies between these radii
    ia, ib = np.flatnonzero(np.abs(r - r_a) < 1e-3), np.flatnonzero(np.abs(r - r_b) < 1e-3)
    ia, ib = ia[np.argsort(ang[ia])], ib[np.argsort(ang[ib])]
    strip = _zipper(ang[ia], ia, ang[ib], ib)
    every = np.vstack([tri, strip])
    bx, by = femfile.element_b(xy, A, every)
    return Field(res["key"], mtri.Triangulation(xy[:, 0], xy[:, 1], every), np.hypot(bx, by), A, len(tri), res)


def flux_step_uwb(fields: list[Field]) -> float:
    """The flux between two neighbouring lines [uWb], one value for every panel: the smallest of FLUX_STEPS_UWB that leaves at most 14 lines each side of zero
    on the field with the largest |A| (flux through a 20 mm slice is A x 0.02 m)."""
    a_max = max(float(np.max(np.abs(f.A))) for f in fields)
    for step in FLUX_STEPS_UWB:
        if a_max * STACK_M * 1e6 / step <= 14:
            return float(step)
    return float(FLUX_STEPS_UWB[-1])


def flux_levels(fields: list[Field], step_uwb: float) -> np.ndarray:
    dA = step_uwb * 1e-6 / STACK_M
    a_max = max(float(np.max(np.abs(f.A))) for f in fields)
    k = int(math.floor(a_max / dA))
    return dA * np.arange(-k, k + 1)


# ============================================================================================== drawing helpers
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": FONT, "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "text.color": INK,
                         "axes.unicode_minus": True, "path.simplify": False})
    return plt


def b_cmap():
    """The |B| colour map: one warm ramp that darkens steadily, warm white (air) to deep wine (saturated steel), with the stops closer together where the steel
    works (0.3 to 1.5 T) so that its loading can be told apart. Values above the top of the scale take the last colour."""
    from matplotlib.colors import LinearSegmentedColormap
    stops = [(0.0, "#fffdf8"), (0.3, "#fde9c4"), (0.7, "#fbc47a"), (1.1, "#f59a52"), (1.5, "#e5603f"), (1.8, "#c2282f"), (2.0, "#8f1230"), (B_MAX, "#4b0c27")]
    cmap = LinearSegmentedColormap.from_list("b_heat", [(v / B_MAX, c) for v, c in stops])
    cmap.set_over(stops[-1][1])
    return cmap


def _path_of(geom):
    """A shapely polygon (or several) as one matplotlib path, holes included (exterior counter-clockwise, holes clockwise)."""
    from matplotlib.path import Path as MPath
    from shapely.geometry.polygon import orient
    verts, codes = [], []
    for pg in getattr(geom, "geoms", [geom]):
        if pg.is_empty:
            continue
        pg = orient(pg, 1.0)
        for ring in (pg.exterior, *pg.interiors):
            xy = np.asarray(ring.coords)
            verts.extend(xy.tolist())
            codes.extend([MPath.MOVETO] + [MPath.LINETO] * (len(xy) - 2) + [MPath.CLOSEPOLY])
    return MPath(verts, codes)


def _halo(width: float):
    import matplotlib.patheffects as pe
    return [pe.Stroke(linewidth=width, foreground=SURFACE), pe.Normal()]


class Canvas:
    """A figure laid out in pixels from its top-left corner."""

    def __init__(self, width_px: int, height_px: int):
        plt = _plt()
        self.W, self.H = int(width_px), int(height_px)
        self.fig = plt.figure(figsize=(self.W / DPI, self.H / DPI), dpi=DPI, facecolor=SURFACE)

    def axes(self, x, y_top, w, h):
        return self.fig.add_axes([x / self.W, 1 - (y_top + h) / self.H, w / self.W, h / self.H])

    def text(self, x, y, s, **kw):
        kw.setdefault("va", "baseline")
        return self.fig.text(x / self.W, 1 - y / self.H, s, **kw)

    def line(self, x0, y0, x1, y1, **kw):
        from matplotlib.lines import Line2D
        self.fig.add_artist(Line2D([x0 / self.W, x1 / self.W], [1 - y0 / self.H, 1 - y1 / self.H], transform=self.fig.transFigure, **kw))

    def rect(self, x, y_top, w, h, **kw):
        from matplotlib.patches import Rectangle
        self.fig.add_artist(Rectangle((x / self.W, 1 - (y_top + h) / self.H), w / self.W, h / self.H, transform=self.fig.transFigure, **kw))

    def overlaps(self) -> list[str]:
        """Pairs of texts whose boxes overlap, and texts that leave the figure: the check behind 'no overlapping text'."""
        import matplotlib.text as mtext
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        boxes = []
        for t in self.fig.findobj(mtext.Text):
            s = t.get_text().strip()
            if s and t.get_visible():
                b = t.get_window_extent(renderer)
                if b.width > 1 and b.height > 1:
                    boxes.append((s, b))
        found = [f"text {s!r} leaves the figure" for s, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > self.W + 0.5 or b.y1 > self.H + 0.5]
        for i, (s1, a) in enumerate(boxes):
            for s2, b in boxes[i + 1:]:
                if a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1:
                    found.append(f"text {s1!r} overlaps {s2!r}")
        return found

    def save(self, path: Path) -> list[str]:
        problems = self.overlaps()
        self.fig.savefig(path, dpi=DPI, facecolor=SURFACE)
        import matplotlib.pyplot as plt
        plt.close(self.fig)
        return [f"{Path(path).name}: {p}" for p in problems]


MARGIN, GAP = 44, 36
COL_W = (W_PX - 2 * MARGIN - GAP) // 2
SCALE = COL_W / (XLIM[1] - XLIM[0])                       # pixels per millimetre, the same in every panel of the two-column figures
COL_X = (MARGIN, MARGIN + COL_W + GAP)


def panel(cv: Canvas, x: float, y_top: float, w_px: float, xlim: tuple[float, float], ylim: tuple[float, float]):
    """An axes whose data units are millimetres, drawn at w_px / (xlim span) pixels each. Returns it and its height in pixels."""
    h = (ylim[1] - ylim[0]) * w_px / (xlim[1] - xlim[0])
    ax = cv.axes(x, y_top, w_px, h)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    return ax, h


# ---------------------------------------------------------------------------------------------- the bodies
def draw_bodies(ax, d: Design, shaft_label: bool = True) -> None:
    """Steel, magnets and coil layers from the visible polygons of the regions (later regions are painted over earlier ones, as in motorgeo.visible)."""
    from matplotlib.patches import PathPatch
    for r in d.regions:
        if r.material in ("Air", "<No Mesh>"):
            continue
        g = d.visible[r.name]
        if g.is_empty:
            continue
        face = STEEL if r.material == "M-19 Steel" else MAGNET if r.material == "N45" else PHASE[r.circuit]
        ax.add_patch(PathPatch(_path_of(g), fc=face, ec=EDGE, lw=0.6, joinstyle="round", zorder=3 if r.material == "N45" else 2))
    if shaft_label:
        ax.text(3.4, 2.6, "shaft", ha="center", va="center", fontsize=9.5, color=MUTED, zorder=3)


def draw_magnet_arrows(ax, d: Design, length: float = 5.0) -> None:
    from matplotlib.patches import FancyArrowPatch
    for r in d.regions:
        if r.material != "N45":
            continue
        c, u = np.array(r.label), np.array(mg.pol(1.0, r.magdir))
        arrow = FancyArrowPatch(tuple(c - u * length / 2), tuple(c + u * length / 2), arrowstyle="-|>", mutation_scale=11, lw=1.35, color=INK,
                                shrinkA=0, shrinkB=0, zorder=4)
        arrow.set_path_effects(_halo(2.8))
        ax.add_patch(arrow)


def draw_ghost(ax, other: Design) -> None:
    """The other design's magnets as a dashed outline, so that what moved can be seen against what is there now."""
    from matplotlib.patches import PathPatch
    for r in other.regions:
        if r.material == "N45":
            ax.add_patch(PathPatch(_path_of(other.visible[r.name]), fc="none", ec=INK, lw=0.9, ls=(0, (3.2, 2.2)), zorder=3.5))


def draw_outlines(ax, d: Design, alpha: float = 0.5) -> None:
    """The bodies' outlines over a colour map, in a translucent dark line, so that the field can be read against the parts."""
    from matplotlib.patches import PathPatch
    for r in d.regions:
        if r.material in ("Air", "<No Mesh>"):
            continue
        g = d.visible[r.name]
        if not g.is_empty:
            ax.add_patch(PathPatch(_path_of(g), fc="none", ec=(0.05, 0.06, 0.08, alpha), lw=0.45, joinstyle="round", zorder=5))


# ---------------------------------------------------------------------------------------------- callouts
def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def stroke(ax, xs, ys, lw: float, **kw) -> None:
    """A black line with a white halo UNDER it (its own z-order), so that the halo of one annotation never cuts another annotation's line."""
    ax.plot(xs, ys, color=SURFACE, lw=lw + 1.8, zorder=5, solid_capstyle=kw.get("solid_capstyle", "round"), solid_joinstyle="round")
    ax.plot(xs, ys, color=INK, lw=lw, zorder=6, **kw)


def bar(ax, p0, p1, tick: float = 1.5, lw: float = 1.6) -> None:
    """A dimension bar: the line between the two ends, and a short tick across each end."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    t = _unit(p1 - p0)
    n = np.array([-t[1], t[0]])
    for a, b in ((p0, p1), (p0 - n * tick / 2, p0 + n * tick / 2), (p1 - n * tick / 2, p1 + n * tick / 2)):
        stroke(ax, [a[0], b[0]], [a[1], b[1]], lw, solid_capstyle="butt")


def leader(ax, pts, lw: float = 0.8) -> None:
    xs, ys = zip(*pts)
    stroke(ax, xs, ys, lw, solid_joinstyle="round", solid_capstyle="round")


def callout_text(ax, x, y, heading: str, value: str, ha: str = "left", va: str = "center"):
    """The feature's name (small, gray) over its value (large, bold), centred on y."""
    t1 = ax.text(x, y + 1.75, heading, ha=ha, va=va, fontsize=HEAD_PT, color=INK2, zorder=8)
    t2 = ax.text(x, y - 1.95, value, ha=ha, va=va, fontsize=VAL_PT, color=INK, fontweight="bold", zorder=8)
    return t1, t2


def _mid(a, b):
    return (np.asarray(a, float) + np.asarray(b, float)) / 2


def draw_callouts(ax, d: Design) -> list[str]:
    """The changed features (centre rib, outer bridge, V angle, slot opening, tooth, yoke): a bar on the feature, a leader, and its name and value. Returns
    warnings (a leader that crosses a magnet slot)."""
    from matplotlib.patches import Arc
    from shapely.geometry import LineString
    from shapely.ops import unary_union
    m, v = d.marks, d.dims
    warnings = []
    magnets = unary_union([r.polygon().buffer(0.15) for r in d.regions if r.material == "N45"])

    # centre rib and outer bridge: bars on the rotor, leaders through the plain steel to labels below the drawing. The rib's leader first follows the pole's
    # d-axis toward the shaft (it cannot touch either slot there), then drops straight down
    bar(ax, m["rib"]["p0"], m["rib"]["p1"])
    rib_mid = _mid(m["rib"]["p0"], m["rib"]["p1"])
    corner = rib_mid + RIB_DOGLEG_MM * np.array(mg.pol(1.0, 180 / d.p.pole_pairs * (POLE0 + 0.5) + 180.0))
    end = (corner[0], LABEL_Y + 0.3)
    for a, b in ((rib_mid, corner), (corner, end)):
        if LineString([tuple(a), tuple(b)]).intersects(magnets):
            warnings.append("rib: the leader crosses a magnet slot")
    leader(ax, [tuple(rib_mid), tuple(corner), end])
    callout_text(ax, corner[0], LABEL_Y - 4.4, "Centre rib", f"{v['rib_mm']:.2f} mm", ha="center")
    bar(ax, m["bridge"]["p0"], m["bridge"]["p1"])
    br_mid = _mid(m["bridge"]["p0"], m["bridge"]["p1"])
    end = (br_mid[0], LABEL_Y + 0.3)
    if LineString([tuple(br_mid), end]).intersects(magnets):
        warnings.append("bridge: the leader crosses a magnet slot")
    leader(ax, [tuple(br_mid), end])
    callout_text(ax, br_mid[0], LABEL_Y - 4.4, "Outer bridge", f"{v['bridge_mm']:.2f} mm", ha="center")

    # the V: the two magnet axes of pole 1, drawn through the apex where they meet (above the steel, under the magnets), the angle between them on the
    # shaft side, and its value beside the arc
    apex, th, b = np.array(m["v"]["apex"]), m["v"]["theta_d"], m["v"]["half"]
    reach = d.p.mag_y / math.sin(math.radians(b)) + d.p.slot_w / 2 + 0.9          # the axis ends a little beyond the magnet slot's far end
    for sgn in (+1, -1):
        u = np.array(mg.pol(1.0, th + sgn * b))
        p_a, p_b = apex - V_LINE_L * u, apex + reach * u
        ax.plot([p_a[0], p_b[0]], [p_a[1], p_b[1]], color=INK, lw=0.8, zorder=2.5)
    ax.add_patch(Arc(tuple(apex), 2 * V_ARC_R, 2 * V_ARC_R, theta1=th + 180.0 - b, theta2=th + 180.0 + b, color=INK, lw=1.4, zorder=6, path_effects=_halo(3.0)))
    where = apex + V_LABEL_R * np.array(mg.pol(1.0, th + 180.0))
    callout_text(ax, where[0], where[1], "V angle", f"{v['v_angle_deg']:.1f}°", ha="center")

    # stator features: a bar on the feature, out along its radius to R_ELBOW, level to SHELF_X, and the label there
    def stator_leader(start, heading, value):
        ang = math.atan2(start[1], start[0])
        elbow = np.array([R_ELBOW * math.cos(ang), R_ELBOW * math.sin(ang)])
        leader(ax, [tuple(start), tuple(elbow), (SHELF_X, elbow[1])])
        callout_text(ax, SHELF_X + 1.0, elbow[1], heading, value, ha="left")

    bar(ax, m["tooth"]["p0"], m["tooth"]["p1"])
    stator_leader(_mid(m["tooth"]["p0"], m["tooth"]["p1"]), "Narrowest tooth", f"{v['tooth_mm']:.2f} mm")
    bar(ax, m["opening"]["p0"], m["opening"]["p1"])
    stator_leader(_mid(m["opening"]["p0"], m["opening"]["p1"]), "Slot opening", f"{v['slot_open_mm']:.2f} mm")
    bar(ax, m["yoke"]["p0"], m["yoke"]["p1"])
    stator_leader(np.array(m["yoke"]["p1"]), "Yoke thickness", f"{v['yoke_mm']:.2f} mm")
    return warnings


# ---------------------------------------------------------------------------------------------- legend
def legend_row(cv: Canvas, y: float, items: list[tuple[str, str, str]]) -> None:
    """One row of swatches centred on the figure. items: (kind, colour, label); kind "swatch", "magnet" (red, with its arrow) or "ghost" (a dashed outline)."""
    from matplotlib.patches import FancyArrowPatch
    renderer = cv.fig.canvas.get_renderer()
    sw_w, sw_h, pad, gap = 30, 18, 9, 34
    widths = []
    for _, _, label in items:
        t = cv.text(0, 0, label, fontsize=LEGEND_PT, color=INK2, va="center")
        widths.append(t.get_window_extent(renderer).width)
        t.remove()
    x = (cv.W - (sum(sw_w + pad + w for w in widths) + gap * (len(items) - 1))) / 2
    for (kind, colour, label), w in zip(items, widths):
        if kind == "ghost":
            cv.rect(x, y - sw_h / 2, sw_w, sw_h, fc="none", ec=INK, lw=0.9, ls=(0, (3.2, 2.2)))
        else:
            cv.rect(x, y - sw_h / 2, sw_w, sw_h, fc=colour, ec=EDGE, lw=0.8)
        if kind == "magnet":
            cx = x + sw_w / 2
            cv.fig.add_artist(FancyArrowPatch(((cx - 8) / cv.W, 1 - y / cv.H), ((cx + 8) / cv.W, 1 - y / cv.H), transform=cv.fig.transFigure,
                                              arrowstyle="-|>", mutation_scale=8, lw=1.1, color=INK, shrinkA=0, shrinkB=0))
        cv.text(x + sw_w + pad, y, label, fontsize=LEGEND_PT, color=INK2, va="center")
        x += sw_w + pad + w + gap


GEOMETRY_LEGEND = [("swatch", STEEL, "Laminated steel"), ("magnet", MAGNET, "N45 magnet, arrow = magnetisation"), ("ghost", "", "dashed: the starting design's magnets"),
                   ("swatch", PHASE["A"], "Phase A"), ("swatch", PHASE["B"], "Phase B"), ("swatch", PHASE["C"], "Phase C"), ("swatch", SURFACE, "Air")]


# ============================================================================================== the panels
CAPTION = "FEMM 5.0 · 2D, 20 mm slice · fine mesh · same copper loss · each design at its own MTPA angle (maximum torque per ampere)"
TITLE = "IPM synRM (36 slots, 8 poles): starting design vs optimized design"


def geometry_panel(cv: Canvas, d: Design, other: Design, col: int, y_top: float) -> tuple[float, list[str]]:
    ax, h = panel(cv, COL_X[col], y_top, COL_W, XLIM, YLIM_GEOMETRY)
    draw_bodies(ax, d)
    if d.key == "optimized":
        draw_ghost(ax, other)
    draw_magnet_arrows(ax, d)
    warnings = draw_callouts(ax, d)
    ax.text(64.0, 63.5, "Geometry", ha="left", va="center", fontsize=HEAD_PT + 1.5, color=INK, fontweight="bold")
    ax.text(64.0, 59.6, "90° sector, 2 of 8 poles", ha="left", va="center", fontsize=HEAD_PT, color=MUTED)
    return h, warnings


def draw_field(ax, d: Design, f: Field, cmap, norm, levels) -> None:
    ax.tripcolor(f.tri, facecolors=f.bmag, cmap=cmap, norm=norm, edgecolors="face", linewidth=0.3, zorder=1)
    ax.tricontour(f.tri, f.A, levels=levels, colors=FLUX_LINE, linewidths=0.5, linestyles="solid", zorder=3)
    draw_outlines(ax, d)


def stats_rows(d: Design, other: Design) -> list[tuple[str, str, float]]:
    """(label, value, value size) of the operating point and what it gives: the study's own numbers (final.json), with the change against the other design
    on the optimized one."""
    s, so = d.study, other.study

    def change(key: str) -> str:
        return f"  ({100 * (s[key] / so[key] - 1):+.0f} %)".replace("-", "−") if d.key == "optimized" else ""

    return [("Current (peak), MTPA angle", f"{s['amps']:.2f} A at {s['gamma']:.1f}°", 13.5),
            ("Mean torque" + change("T_mean"), f"{s['T_mean']:.3f} N·m", 17.0),
            ("Torque ripple" + change("ripple_pct"), f"{s['ripple_pct']:.2f} %", 17.0),
            ("Magnet flux linkage" + change("psi_pm_mWb"), f"{s['psi_pm_mWb']:.2f} mWb", 17.0)]


def field_panel(cv: Canvas, d: Design, other: Design, f: Field, col: int, y_top: float, norm, cmap, levels) -> float:
    ax, h = panel(cv, COL_X[col], y_top, COL_W, XLIM, YLIM_FIELD)
    draw_field(ax, d, f, cmap, norm, levels)
    ax.text(64.0, 63.5, "Magnetic field", ha="left", va="center", fontsize=HEAD_PT + 1.5, color=INK, fontweight="bold")
    ax.text(64.0, 59.6, "at the operating point", ha="left", va="center", fontsize=HEAD_PT, color=MUTED)
    y = 50.0
    for label, value, size in stats_rows(d, other):
        ax.text(STATS_X, y + 2.2, label, ha="left", va="center", fontsize=HEAD_PT, color=INK2)
        ax.text(STATS_X, y - 2.0, value, ha="left", va="center", fontsize=size, color=INK, fontweight="bold")
        y -= 11.0
    return h


def header(cv: Canvas, title: str, caption: str) -> None:
    cv.text(MARGIN, 66, title, fontsize=TITLE_PT, color=INK, fontweight="bold")
    cv.text(MARGIN, 104, caption, fontsize=CAPTION_PT, color=INK2)


def column_headers(cv: Canvas, designs: list[Design], y: float) -> None:
    for col, d in enumerate(designs):
        cv.text(COL_X[col], y, d.title, fontsize=COLUMN_PT, color=INK, fontweight="bold")
        cv.line(COL_X[col], y + 14, COL_X[col] + COL_W, y + 14, color=AXIS, lw=1.0)


def notes(cv: Canvas, y: float, lines: list[str]) -> float:
    for k, text in enumerate(lines):
        cv.text(MARGIN, y + k * 25, text, fontsize=NOTE_PT, color=INK2)
    return y + 25 * len(lines)


def colour_bar(cv: Canvas, y: float, norm, cmap, step_uwb: float, x_center: float | None = None, width: float = 640) -> None:
    """The shared |B| scale with its label on the left and what a flux line is on the right."""
    from matplotlib.cm import ScalarMappable
    x0 = ((cv.W if x_center is None else 2 * x_center) - width) / 2
    cax = cv.axes(x0, y, width, 20)
    cb = cv.fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", extend="max", extendfrac=0.025, ticks=[0, 0.5, 1.0, 1.5, 2.0, B_MAX])
    cb.outline.set_edgecolor(AXIS)
    cb.ax.tick_params(labelsize=LEGEND_PT, colors=INK2, length=4, color=AXIS)
    cb.ax.set_xticklabels(["0", "0.5", "1.0", "1.5", "2.0", "2.2"])
    cv.text(x0 - 18, y + 16, "Flux density |B| (T)", fontsize=LEGEND_PT + 0.5, color=INK2, ha="right")
    cv.text(x0 + width + 30, y + 16, f"Flux lines: one per {step_uwb:g} µWb", fontsize=LEGEND_PT + 0.5, color=INK2, ha="left")


# ============================================================================================== the figures
def figure_geometry(designs: list[Design], path: Path) -> list[str]:
    """The top row alone: both geometries at the same scale, the changed features called out."""
    y_cols = 176
    y_row = y_cols + 28
    h_row = (YLIM_GEOMETRY[1] - YLIM_GEOMETRY[0]) * SCALE
    y_leg = y_row + h_row + 30
    y_note = y_leg + 52
    note_lines = geometry_notes(designs)
    cv = Canvas(W_PX, y_note + 25 * len(note_lines) + 28)
    header(cv, TITLE, CAPTION)
    column_headers(cv, designs, y_cols)
    warnings = []
    for col, d in enumerate(designs):
        _, w = geometry_panel(cv, d, designs[1 - col], col, y_row)
        warnings += [f"{d.key}: {x}" for x in w]
    legend_row(cv, y_leg, GEOMETRY_LEGEND)
    notes(cv, y_note, note_lines)
    return warnings + cv.save(path)


def geometry_notes(designs: list[Design]) -> list[str]:
    a, b = designs
    more = 100 * (b.dims["coil_layer_area_mm2"] / a.dims["coil_layer_area_mm2"] - 1)
    return ["Both sectors are drawn at the same scale; each bar marks where a dimension is measured. The V angle is the angle between the two magnets of a pole.",
            f"The optimized design's slots have {more:.0f} % more coil area per layer ({b.dims['coil_layer_area_mm2']:.1f} against {a.dims['coil_layer_area_mm2']:.1f} mm²), "
            f"so at the same copper loss it carries {b.study['amps']:.2f} A against {a.study['amps']:.2f} A.",
            "Both designs use the same 16 N45 magnets, drawn as modelled: each magnet slot is filled with N45."]


def figure_main(designs: list[Design], fields: dict[str, Field], path: Path, step_uwb: float, levels: np.ndarray, own: dict[str, dict]) -> list[str]:
    """The whole figure: geometry on top, the field under it, one colour scale."""
    from matplotlib.colors import Normalize
    cmap, norm = b_cmap(), Normalize(vmin=0.0, vmax=B_MAX)
    y_cols = 176
    y_geo = y_cols + 28
    y_leg = y_geo + (YLIM_GEOMETRY[1] - YLIM_GEOMETRY[0]) * SCALE + 26
    y_fld = y_leg + 30
    y_bar = y_fld + (YLIM_FIELD[1] - YLIM_FIELD[0]) * SCALE + 30
    y_note = y_bar + 92
    note_lines = main_notes(designs, step_uwb, own)
    cv = Canvas(W_PX, y_note + 25 * len(note_lines) + 24)
    header(cv, TITLE, CAPTION)
    column_headers(cv, designs, y_cols)
    warnings = []
    for col, d in enumerate(designs):
        _, w = geometry_panel(cv, d, designs[1 - col], col, y_geo)
        warnings += [f"{d.key}: {x}" for x in w]
    legend_row(cv, y_leg, GEOMETRY_LEGEND)
    for col, d in enumerate(designs):
        field_panel(cv, d, designs[1 - col], fields[d.key], col, y_fld, norm, cmap, levels)
    colour_bar(cv, y_bar, norm, cmap, step_uwb)
    notes(cv, y_note, note_lines)
    return warnings + cv.save(path)


def main_notes(designs: list[Design], step_uwb: float, own: dict[str, dict]) -> list[str]:
    a, b = designs
    return ["Torque and ripple: the study's fine-mesh scoring (mean torque; ripple = peak-to-peak as a share of the mean; 15 rotor positions over one ripple period).",
            "Magnet flux linkage: the no-load flux linkage on the d-axis (peak, per phase), from the same N45 magnets. Currents are A peak at the 10-turn placeholder winding.",
            "The optimized design's larger slots carry more current at the same copper loss. Current angle: from the q-axis toward −d (d = the magnet axis).",
            f"The field is one rotor position (rotor as drawn, the same electrical angle in both; torque there {own[a.key]['T_view']:.3f} and {own[b.key]['T_view']:.3f} N·m). Colours above 2.2 T take the last colour."]


def figure_wide(designs: list[Design], fields: dict[str, Field], path: Path, step_uwb: float, levels: np.ndarray, own: dict[str, dict]) -> list[str]:
    """A wide single row: for each design its geometry and its field side by side, the numbers under them, one colour scale at the foot."""
    from matplotlib.colors import Normalize
    cmap, norm = b_cmap(), Normalize(vmin=0.0, vmax=B_MAX)
    W = 2000
    margin, gap_in, gap_out = 40, 14, 56
    w_panel = (W - 2 * margin - 2 * gap_in - gap_out) / 4
    xlim, ylim = (-2.0, 66.0), (-2.0, 66.0)
    h_panel = (ylim[1] - ylim[0]) * w_panel / (xlim[1] - xlim[0])
    y_cols, y_cap, y_row = 168, 218, 232
    y_stats = y_row + h_panel + 44
    y_change = y_stats + 82 + 126
    y_bar = y_change + 62
    y_note = y_bar + 100
    note = ["Same basis in both: FEMM 5.0, 2D, 20 mm slice, fine mesh, the same copper loss, each design at its own MTPA current angle.",
            "Currents are A peak at the 10-turn placeholder winding. Colours above 2.2 T take the last colour."]
    cv = Canvas(W, y_note + 25 * len(note) + 20)
    cv.text(margin, 66, TITLE, fontsize=TITLE_PT, color=INK, fontweight="bold")
    cv.text(margin, 104, CAPTION, fontsize=CAPTION_PT, color=INK2)
    xs = [margin + k * (w_panel + gap_in) + (gap_out - gap_in) * (k // 2) for k in range(4)]
    for k, d in enumerate(designs):
        x0, x1 = xs[2 * k], xs[2 * k + 1] + w_panel
        cv.text(x0, y_cols, d.title, fontsize=COLUMN_PT, color=INK, fontweight="bold")
        cv.line(x0, y_cols + 14, x1, y_cols + 14, color=AXIS, lw=1.0)
        cv.text(xs[2 * k], y_cap, "Geometry", fontsize=HEAD_PT + 0.5, color=INK2, fontweight="bold")
        cv.text(xs[2 * k + 1], y_cap, "Magnetic field", fontsize=HEAD_PT + 0.5, color=INK2, fontweight="bold")
        ax, _ = panel(cv, xs[2 * k], y_row, w_panel, xlim, ylim)
        draw_bodies(ax, d, shaft_label=False)
        draw_magnet_arrows(ax, d, length=4.4)
        ax2, _ = panel(cv, xs[2 * k + 1], y_row, w_panel, xlim, ylim)
        draw_field(ax2, d, fields[d.key], cmap, norm, levels)
        for j, (label, value, size) in enumerate(stats_rows(d, designs[1 - k])):
            cx, cy = x0 + (j % 2) * (x1 - x0) / 2, y_stats + (j // 2) * 82
            cv.text(cx, cy, label, fontsize=HEAD_PT, color=INK2)
            cv.text(cx, cy + 40, value, fontsize=size + 1, color=INK, fontweight="bold")
    a, b = designs
    changes = (f"Changed:   centre rib {a.dims['rib_mm']:.2f} → {b.dims['rib_mm']:.2f} mm   outer bridge {a.dims['bridge_mm']:.2f} → {b.dims['bridge_mm']:.2f} mm   "
               f"V angle {a.dims['v_angle_deg']:.1f} → {b.dims['v_angle_deg']:.1f}°   slot opening {a.dims['slot_open_mm']:.2f} → {b.dims['slot_open_mm']:.2f} mm   "
               f"tooth width {a.dims['tooth_mm']:.2f} → {b.dims['tooth_mm']:.2f} mm   yoke {a.dims['yoke_mm']:.2f} → {b.dims['yoke_mm']:.2f} mm")
    cv.text(margin, y_change, changes, fontsize=LEGEND_PT + 1.0, color=INK)
    colour_bar(cv, y_bar, norm, cmap, step_uwb, x_center=W / 2, width=700)
    for k, text in enumerate(note):
        cv.text(margin, y_note + k * 25, text, fontsize=NOTE_PT, color=INK2)
    return cv.save(path)


# ============================================================================================== checks and the numbers file
# The values the study's write-up quotes, to the digits it prints them. They are NOT what the figure prints (that comes from the parameters); the check holds the one against the other.
QUOTED = {"reference": {"bridge_mm": 1.558, "rib_mm": 3.298, "v_angle_deg": 120.0, "slot_open_mm": 2.0, "tooth_mm": 3.97, "yoke_mm": 5.70, "psi_pm_mWb": 16.64,
                        "amps": 10.000, "T_mean": 1.1288, "ripple_pct": 18.45},
          "optimized": {"bridge_mm": 1.400, "rib_mm": 1.000, "v_angle_deg": 111.0, "slot_open_mm": 1.5, "tooth_mm": 3.00, "yoke_mm": 4.16, "psi_pm_mWb": 24.20,
                        "amps": 11.8627, "T_mean": 1.8073, "ripple_pct": 4.84}}
TOL_MM, TOL_DEG, TOL_PSI_MWB = 0.01, 0.1, 0.01
TOL_TORQUE_REL = 0.01                # this script's own torque must be within 1 % of the study's fine-mesh torque
TOL_RIPPLE_PT = 0.2                  # ripple in percentage points


def own_numbers(res: dict) -> dict:
    """Mean torque, ripple and psi_pm from this script's own solves (the same scoring as design.score, at the study's current and angle)."""
    t = np.asarray(res["torque"], float)
    mean = float(t.mean())
    pp = design.fourier_pkpk(t)
    return {"T_mean": mean, "ripple_pct": 100.0 * pp / mean, "ripple_pp": pp, "psi_pm_mWb": 1e3 * res["psi_pm"], "T_view": res["view"]["T_age"]}


def verify(designs: list[Design], results: dict[str, dict], fields: dict[str, Field], final: dict) -> list[dict]:
    checks = []

    def add(name, got, want, tol, kind="abs"):
        err = abs(got - want) if kind == "abs" else abs(got / want - 1.0)
        checks.append({"name": name, "got": float(got), "want": float(want), "tol": tol, "kind": kind, "error": float(err), "ok": bool(err <= tol)})

    for d in designs:
        q, dm, own = QUOTED[d.key], d.dims, own_numbers(results[d.key])
        for k, tol in (("bridge_mm", TOL_MM), ("rib_mm", TOL_MM), ("v_angle_deg", TOL_DEG), ("slot_open_mm", TOL_MM), ("tooth_mm", TOL_MM), ("yoke_mm", TOL_MM)):
            add(f"{d.key}: {k} against the quoted value", dm[k], q[k], tol)
        add(f"{d.key}: outer bridge, design.bridge_of against motorgeo.check", dm["bridge_mm"], dm["bridge_check_mm"], 1e-4)
        add(f"{d.key}: outer bridge against the study's record", dm["bridge_mm"], d.study["bridge_mm"], TOL_MM)
        for name, key in (("rib", "rib_mm"), ("bridge", "bridge_mm"), ("tooth", "tooth_mm"), ("yoke", "yoke_mm"), ("opening", "slot_open_mm")):
            add(f"{d.key}: the drawn {name} bar against its dimension", d.marks[name]["mm"], dm[key], 2e-3)
        add(f"{d.key}: current at the same copper loss (design.scoring_amps, scoring P) against the study's", design.scoring_amps(d.metrics, "P"), d.study["amps"], 1e-6)
        add(f"{d.key}: current against the quoted value", d.study["amps"], q["amps"], 5e-4)
        add(f"{d.key}: own mean torque against the study's fine-mesh torque", own["T_mean"], d.study["T_mean"], TOL_TORQUE_REL, "rel")
        add(f"{d.key}: own ripple against the study's", own["ripple_pct"], d.study["ripple_pct"], TOL_RIPPLE_PT)
        add(f"{d.key}: own magnet flux linkage against the study's", own["psi_pm_mWb"], d.study["psi_pm_mWb"], TOL_PSI_MWB)
        add(f"{d.key}: magnet flux linkage against the quoted value", d.study["psi_pm_mWb"], q["psi_pm_mWb"], TOL_PSI_MWB)
        add(f"{d.key}: study's torque against the quoted value", d.study["T_mean"], q["T_mean"], 5e-4)
        add(f"{d.key}: study's ripple against the quoted value", d.study["ripple_pct"], q["ripple_pct"], 0.01)
        add(f"{d.key}: the picture's rotor stands as drawn (rotor angle, deg)", results[d.key]["view"]["rotor_deg"], 0.0, 1e-6)
        magnet_mesh = float(fields[d.key].res["n45_area_mm2"])
        magnet_poly = sum(r.polygon().area for r in d.regions if r.material == "N45")
        add(f"{d.key}: N45 area of the mesh against the drawn magnets (mm2)", magnet_mesh, magnet_poly, 5e-3, "rel")
    opt = designs[1]
    add("optimized: outer bridge against final.json", opt.dims["bridge_mm"], final["bridge_mm"], TOL_MM)
    add("optimized: V angle against final.json (2 x v_half_deg)", opt.dims["v_angle_deg"], 2 * final["params"]["v_half_deg"], TOL_DEG)
    add("optimized: slot opening against final.json", opt.dims["slot_open_mm"], final["params"]["slot_open"], 1e-9)
    add("optimized: study's torque against final.json's fine block", opt.study["T_mean"], final["fine"]["T_mean"], 1e-9)
    add("optimized: study's ripple against final.json's fine block", opt.study["ripple_pct"], final["fine"]["ripple_pct"], 1e-9)
    add("optimized: study's angle against final.json's fine block", opt.study["gamma"], final["fine"]["gamma"], 1e-9)
    add("optimized: study's current against final.json's fine block", opt.study["amps"], final["fine"]["amps"], 1e-9)
    return checks


def write_numbers(designs, results, fields, final, checks, step_uwb, levels, warnings) -> Path:
    doc = {"what": "the numbers printed in reference_vs_optimized*.png, each with where it comes from",
           "basis": {"solver": "FEMM 5.0 (femm_queue.femm_pool)", "model": "2D planar, 20 mm slice, 90 degree sector, fine mesh (smart mesh on)", "winding": "10 turns per coil side (placeholder)",
                     "comparison": "the study's: the same copper loss (scoring P), each design at its own MTPA current angle; currents in A peak",
                     "field_position": "one rotor position: the rotor as drawn (rotor angle 0), the same electrical angle in both designs",
                     "torque_and_ripple": "the study's fine-mesh numbers (final.json compare); this script's own solves of the same 15 positions are held against them in checks",
                     "flux_lines": {"flux_per_line_uWb": step_uwb, "lines": int(len(levels)), "stack_m": STACK_M},
                     "colour_scale_T": [0.0, B_MAX], "colour_map": CMAP_NAME},
           "designs": {}, "checks": checks, "warnings": warnings, "all_checks_ok": all(c["ok"] for c in checks)}
    for d in designs:
        res, own = results[d.key], own_numbers(results[d.key])
        doc["designs"][d.key] = {
            "overrides": d.overrides,
            "printed": {"dimensions": {k: d.dims[k] for k in ("bridge_mm", "rib_mm", "v_angle_deg", "slot_open_mm", "tooth_mm", "yoke_mm")},
                        "operating_point": {"amps_peak": d.study["amps"], "gamma_deg": d.study["gamma"], "mean_torque_Nm": d.study["T_mean"], "ripple_pct": d.study["ripple_pct"]}},
            "sources": {"dimensions": "design.bridge_of, motorgeo.check and the parameters (Params), computed here",
                        "operating_point": "final.json compare row (the study's fine-mesh scoring)"},
            "other_dimensions": {"h_open_mm": d.dims["h_open_mm"], "slot_reach_mm": d.dims["slot_reach_mm"], "coil_layer_area_mm2": d.dims["coil_layer_area_mm2"],
                                 "magnet_flux_linkage_mWb_study": d.study["psi_pm_mWb"]},
            "own_solves": {"T_mean_Nm": own["T_mean"], "ripple_pct": own["ripple_pct"], "psi_pm_mWb": own["psi_pm_mWb"], "torque_at_the_15_positions_Nm": res["torque"],
                           "torque_at_the_pictures_position_Nm": own["T_view"], "picture_theta_e_deg": res["view"]["theta_e"], "picture_rotor_angle_deg": res["view"]["rotor_deg"],
                           "mesh_nodes": res["nodes"], "mesh_triangles": res["triangles"], "seconds": res["seconds"]},
            "max_flux_density_T": {"steel": fields[d.key].res["b_steel_max"], "everywhere": fields[d.key].res["b_max"],
                                   "share_of_steel_area_above_scale": fields[d.key].res["steel_above_scale"]},
        }
    path = OUT / JSON_NUMBERS
    path.write_text(json.dumps(doc, indent=1, default=float), encoding="utf-8")
    return path


def annotate_field(f: Field, d: Design) -> None:
    """Numbers about the mesh that the checks and the JSON use: the magnets' area in the mesh, and how much steel is above the colour scale."""
    r = f.res
    corners = r["xy"][r["tri"]]
    area = 0.5 * np.abs((corners[:, 1, 0] - corners[:, 0, 0]) * (corners[:, 2, 1] - corners[:, 0, 1]) - (corners[:, 2, 0] - corners[:, 0, 0]) * (corners[:, 1, 1] - corners[:, 0, 1]))
    material = np.asarray(r["material"])
    b = f.bmag[:f.n_mesh]
    steel = material == "M-19 Steel"
    r["n45_area_mm2"] = float(area[material == "N45"].sum())
    r["b_steel_max"] = float(b[steel].max())
    r["b_max"] = float(b.max())
    r["steel_above_scale"] = float(area[steel & (b > B_MAX)].sum() / area[steel].sum())


# ============================================================================================== main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reuse", action="store_true", help="draw from the fields saved by an earlier run (no FEMM, no solver slot)")
    ap.add_argument("--workers", type=int, default=2, help="FEMM workers (one per design; at most 8)")
    args = ap.parse_args(argv)
    if not FINAL_JSON.exists():
        raise SystemExit(f"{FINAL_JSON} is not there: the study's last step (optimization.finalize) writes it")
    final = json.loads(FINAL_JSON.read_text(encoding="utf-8"))
    rows = {r["design"]: r for r in final["compare"]}
    designs = [Design("reference", "Starting design", {}, rows["reference"]), Design("optimized", "Optimized design", final["overrides"], rows["chosen (fine knee)"])]
    OUT.mkdir(parents=True, exist_ok=True)
    if args.reuse and not all((OUT / f"field_{d.key}.npz").exists() for d in designs):
        raise SystemExit(f"--reuse needs the fields of an earlier run in {OUT}; run once without it")
    results = {d.key: load_result(d.key) for d in designs} if args.reuse else solve_designs(designs, args.workers)
    fields = {d.key: build_field(results[d.key], d.p) for d in designs}
    for d in designs:
        annotate_field(fields[d.key], d)
    step = flux_step_uwb(list(fields.values()))
    levels = flux_levels(list(fields.values()), step)
    own = {d.key: own_numbers(results[d.key]) for d in designs}
    warnings = []
    warnings += figure_main(designs, fields, OUT / PNG_MAIN, step, levels, own)
    warnings += figure_geometry(designs, OUT / PNG_GEOMETRY)
    warnings += figure_wide(designs, fields, OUT / PNG_WIDE, step, levels, own)
    checks = verify(designs, results, fields, final)
    path = write_numbers(designs, results, fields, final, checks, step, levels, warnings)
    for c in checks:
        print(f"{'PASS' if c['ok'] else 'FAIL'}  {c['name']}: got {c['got']:.6g}, want {c['want']:.6g} (error {c['error']:.3g}, tolerance {c['tol']:g}{' relative' if c['kind'] == 'rel' else ''})")
    for w in warnings:
        print("WARNING", w)
    print(f"flux lines: one per {step:g} uWb ({len(levels)} levels); wrote {OUT / PNG_MAIN}, {OUT / PNG_GEOMETRY}, {OUT / PNG_WIDE} and {path}")
    return 0 if all(c["ok"] for c in checks) and not warnings else 1


if __name__ == "__main__":
    raise SystemExit(main())
