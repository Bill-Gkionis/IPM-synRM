"""optimization.sweeps_stator: the stator tier looked at one variable at a time, so that the optimizer is chosen after LOOKING at the landscape.

What it is, and why
    An optimizer is a way of walking a landscape, and which one suits depends on the landscape: smooth or rough, one valley or several, a best
    point inside the box or on its wall, and how much of what the black box says is the motor and how much is the mesh. So the study looks
    first. This script moves each of the five stator variables of optimization/design.py alone (the others at the reference, the band-free
    rotor Params(barrier="none")) across the whole of its feasible interval, scores every step with design.Evaluator, and writes the rows
    for sweeps_stator_figures.py to read. It starts FEMM only through femm_queue.femm_pool (one heavy run at a time, R10).
    The variables (motorgeo.Params fields): slot_half_deg [deg], half the angle of the radial-walled slot body (the tooth body spans
    10 - 2 x this); slot_open [mm], the opening at the bore; h_open [mm], the straight radial part of the opening (the tooth tip's thinnest
    height); r_shoe [mm], where the shoe ends and the walls start; r_wall [mm], where the walls end and the round bottom starts (the slot's
    depth, so the yoke). A word about "fillets": the slot bottom is the only rounded part motorgeo draws. It is tangent to both walls and
    has the radius R = r_wall x tan(slot_half_deg) (2.50 mm on the reference), so it moves with r_wall and slot_half_deg and there is no
    fillet variable. The shoe and opening corners are sharp in the model; a stamped lamination rounds them (not modelled).

Method, and the maths
    screen   Each variable alone is scanned far past any limit that could be wanted (SCAN: slot_half_deg 0..5.5 deg, slot_open 0..8 mm,
             h_open 0..3 mm, r_shoe 46..56 mm, r_wall 46..64 mm; 3,855 values in all) with design.screen, which needs no FEMM. The end of
             the run that holds the reference is then bisected between the last feasible and the first infeasible value until they are
             BISECT_TOL = 1e-5 apart (log2(0.01 / 1e-5) = 10 halvings), and the constraints that are broken just beyond it are named.
             One of them is NOT in design.G_NAMES: motorgeo._slot cuts the slot into two coil layers of equal area at the radius r_mid, and
             that boundary has to lie inside the straight walls. With S = (r_wall^2 - r_shoe^2) sin(2a) / 2 the straight part's area and
             B = R^2 (pi + 2a + sin 2a) / 2 the round bottom's (a = slot_half_deg, R = r_wall tan a), r_mid <= r_wall is the same as S >= B.
             At S = B the lower layer turns inside out ("invalid shapes: coil_k_1") while every G is still <= 0: a hidden constraint with no
             number to tell an optimizer how far away it is. layer_split() computes r_mid - r_wall so that it can be seen (it is 0.00000
             at both of the ends it limits: S = B = 10.395 mm2 at r_shoe 53.62215, 7.9186 mm2 at r_wall 48.65294).
    sweeps   About 20 values per variable: 13 even steps across the interval and 3 more at 1.5, 3.5 and 5.5 % of its length from each end,
             because an optimizer held by a constraint spends its time at the end. The ends sit MARGIN = 0.001 inside their constraint
             (rounding of the cache key, 6 decimals, cannot flip them); an end that no constraint limits (a zero-width opening, a knife-edge
             tip) starts at SOFT_LOW = 0.01 mm, 10 micrometres, absurdly small on purpose. The reference row {} is in every sweep.
    J, AT    scoring "J" gives each candidate the reference slot's current density, I = 10 A x A / A_ref (A = its coil layer's area), so
             J_rel = (I / A) / (10 A / A_ref) = 1 and the copper loss is Pcu_rel = (I / 10 A)^2 x A_ref / A = A / A_ref. Scoring "AT" gives
             10 A (100 ampere-turns) to every slot: J_rel = Pcu_rel = A_ref / A. A variable gets an AT sweep too when its coil area moves
             by more than AT_AREA_SPAN = 3 % across its interval (the screen's area_span_rel): slot_half_deg, r_shoe and r_wall. slot_open
             and h_open do not touch the coil regions, so J and AT are the same numbers for them.
    noise    The reference moved by +-0.001, +-0.002 and +-0.003 (mm; degrees for slot_half_deg, where 0.001 deg is 1 micrometre at the
             slot wall) in each variable: 31 candidates that differ from the reference by micrometres, far too little to change the motor.
             What they still differ by is the mesh. The ladder repeats it with steps of 0.01, 0.03 and 0.1: how the noise of a finite
             difference, sigma / (sqrt(2) h) for a central difference with step h, gives way to curvature as h grows.
    fine     Re-scores on the fine mesh (its own Evaluator, evaluations_fine.csv): both ends of each sweep, the lowest-ripple point, the
             highest-torque point and the neighbours of the lowest-ripple point (plus, by --extra, two points on either side of a step
             that the coarse mesh showed): 23 candidates. The question is whether the order of the points, and so the choice of the best,
             survives a different mesh.
    map      11 x 11 values of two variables across both their feasible intervals, scoring J, one batch. The cells that design.screen
             refuses are written as rows too (with the reason), so the map shows where the room ends.

Run (from the Python to FEMM folder; every stage takes the solver slot once, and waits, printing a line every 30 s, if another run holds it)
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator screen      the intervals, no FEMM, 2 minutes: feasible_intervals.json, screen_scan.csv
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator j           5 sweeps under J (a variable's name after it runs only that one)
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator at          the sweeps under AT, for the variables whose coil area moves
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator noise       the micrometre perturbations (--steps 0.01,0.03,0.1 --label ladder: the ladder)
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator fine --neighbours --extra h_open=0.536664,0.602497
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator map slot_open slot_half_deg      (then lo_x hi_x lo_y hi_y to change the box)
    (uv run python -m optimization.sweeps_stator <stage> is the same.) A stage that is in its CSV already is not run again, so a stage that
    failed half way is repeated by running it again; delete its rows (not the file) to redo a finished one. Then
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures tables|shapes|noise|ladder|fine|map|derived|cost|figures
Writes out\\optimization\\sweeps_stator\\: evaluations.csv (J, noise, ladder and map rows), evaluations_AT.csv, evaluations_fine.csv, the two screen
files, logs\\ (what each stage printed) and the figures. Every row carries its stage label (the variable's name, "<name>_AT", "noise", "ladder", "fine", "map").

Numbers to expect (this PC, 2026-11-08, FEMM 5.0, coarse mesh unless said; FEMM gave the same numbers to the last digit in three sessions)
    intervals  slot_half_deg [1.5511, 3.7757] deg  low: open_vs_slot (the body must be 0.3 deg wider each side than the opening), high: tooth_min = 2.0 mm
               slot_open     [0.00001, 3.6281] mm  low: none (zero width, "invalid shapes"), high: open_vs_slot
               h_open        [0, 0.8010] mm        low: none (feasible even at 0, a knife edge), high: shoe_taper
               r_shoe        [46.3000, 53.6222] mm low: shoe_taper, high: the hidden layer split (S = B)
               r_wall        [48.6529, 58.3234] mm low: the hidden layer split, high: yoke = 3.0 mm
               each is one run, no gaps; coil layer area x0.558..x1.615 (slot_half_deg), x1 (slot_open, h_open), x0.404..x1.041 (r_shoe), x0.308..x1.275 (r_wall)
    reference  gamma 22.337 deg, T_mean 1.12173 N.m, ripple 19.619 % (0.2201 N.m), 8,892 nodes; fine mesh 22.457 deg, 1.12880 N.m, 18.448 %, 15,859 nodes
    size       345 candidates solved (7,935 FEMM solves: J 96, AT 58, noise 31, ladder 31, fine 23, map 106); 450 s of pool time (7.5 minutes;
               a sweep of 20 takes 21 to 26 s), 46 worker-minutes; the screen 132 s
    under J    slot_half_deg  T 0.577 -> 1.897 N.m, ripple 21.96 -> 12.10 %, both monotonic, best at the tooth limit (5.5 -> 16.2 A, Pcu_rel 0.55 -> 1.62)
               slot_open      T peaks at 1.185 (0.31 mm), 1.170 at 0.01, 1.014 at 3.63; ripple 8.23 % at 0.01, 10.7 % on a plateau at 0.2-0.9, 19.6 %
                              at 2.0, a hump of 24.4 % at 3.0, 22.5 % at 3.63: the best is the closed end
               h_open         flat: T 1.1216-1.1243, ripple 19.06-19.67 %, with steps of 0.3 point that the fine mesh does not repeat
               r_shoe         T 1.170 -> 0.412, ripple 18.41 -> 24.13 %: both best at the low end (shoe_taper)
               r_wall         T 0.311 -> 1.420, ripple 28.2 -> 18.97 % at 57.5 mm and then flat (19.14 at the yoke limit)
               T follows the coil area exactly (Spearman 1.00 for slot_half_deg, r_shoe, r_wall): J pays for slot area in current
    under AT   T stays between -1.6 % and +0.3 % of the reference (-2.8 % at r_wall's yoke limit) while Pcu_rel = J_rel runs 0.62 to 1.80
               (slot_half_deg), 0.96 to 2.48 (r_shoe), 0.78 to 3.25 (r_wall); the highest-AT torque of each is only +0.3 % and sits at
               1.7 to 2.8 times the reference current density: the trap design.py describes
    noise      micrometre steps (31 candidates): sigma 0.00056 N.m (0.05 %) and 0.0194 %-points pooled; 0.00001-0.0004 and 0.0012-0.0237 per
               variable once a straight line is taken out; step-scale roughness between neighbouring sweep points 0.00022 N.m and 0.059 points
    ladder     slot_open ripple slope at h = 0.001, 0.002, 0.003, 0.01, 0.03, 0.1: +6.1, -6.4, -0.4, +6.1, +7.6, +7.9 per mm (the sweep says +7.4)
    fine mesh  27 pairs: torque +0.47 % (sd 0.22), ripple -0.86 points (sd 0.59, from +0.93 to -1.86); the order of the points held in
               every sweep (Kendall tau 1.00) except r_wall's ripple (0.87: its two lowest points swap) and h_open (0.57: the coarse
               mesh does not see h_open, the fine mesh does)
    map        slot_open x slot_half_deg: 106 of 121 cells scored in 115 s, ripple 4.23-25.15 %, lowest at (0.01 mm, 3.7747 deg), 3 valleys above
               the noise, interaction share 3.3 % of the ripple's variance

Self-check (the answers follow).
  1. The J sweep of r_wall shows the torque rising all the way to the yoke limit. Why, and what has to be added to the black box's
     objectives so that an optimizer does not simply follow it?
  2. A finite difference of the ripple of slot_open with a step of 0.001 mm gives +6.1 points per mm, and with 0.002 mm gives -6.4. The
     sweep's slope is about +7.4. Which of the two is a gradient?
  3. Why are the sweeps' ends 0.001 inside the feasible interval and not on its end, and why is the interval found by bisection after a scan?
Answers.
  1. J gives every slot the reference current density, so a bigger slot gets more current, more ampere-turns, more torque (T follows the
     coil area, Spearman 1.00) and the same percentage ripple falls as the torque rises. The copper loss grows with the same slot (Pcu_rel
     = A / A_ref, 1.275 at the yoke limit) and nobody is charged for it. The copper loss has to be a third objective or a constraint
     (Pcu_rel <= 1, or the coil area no larger than the reference's); without it the optimum is the corner of the biggest slot.
  2. Neither. At steps of a micrometre the two copies of the design have different meshes, and the difference of two torques that differ
     by the mesh's noise (0.02 points of ripple) is divided by 2 x 0.001 mm: the noise of the slope is 0.0194 / (sqrt(2) x 0.001) = 14
     points per mm, twice the slope itself. A step of 0.03 mm or more (7.6, 7.9) agrees with the sweep.
  3. The sweep's candidates are cached and written by a key that rounds the values to 6 decimals, and a design that sits exactly on a limit
     must not be refused or accepted by that rounding. The scan only brackets an end to its own step (0.005 to 0.02); bisection finds it
     to 1e-5 in about ten screens, and the scan, not the bisection, is what would find a second island of feasible designs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

import femm_queue
from optimization import design

RUN = "sweeps_stator"
VARS = design.STATOR_VARS                      # ("slot_half_deg", "slot_open", "h_open", "r_shoe", "r_wall")
REF = design.REFERENCE
UNIT = {"slot_half_deg": "deg", "slot_open": "mm", "h_open": "mm", "r_shoe": "mm", "r_wall": "mm"}
CSV_J, CSV_AT, CSV_FINE = "evaluations.csv", "evaluations_AT.csv", "evaluations_fine.csv"

# The wide scans: (first, last, step) of each variable alone, the others at the reference. They run far past any end that could be wanted, so
# that the end found is the constraint's, not the scan's.
SCAN = {"slot_half_deg": (0.0, 5.5, 0.01), "slot_open": (0.0, 8.0, 0.01), "h_open": (0.0, 3.0, 0.005),
        "r_shoe": (46.0, 56.0, 0.01), "r_wall": (46.0, 64.0, 0.02)}
BISECT_TOL = 1e-5                              # an end is found to this much (mm or degrees)
MARGIN = 1e-3                                  # a sweep's end sits this far inside the constraint, so that rounding cannot flip it
SOFT_LOW = {"slot_open": 0.01, "h_open": 0.01}  # where a sweep starts when no constraint ends the low side: 10 micrometres, absurdly small on purpose
AT_AREA_SPAN = 0.03                            # "a few per cent": a variable whose coil area moves more than this across its interval also gets the AT sweep
UNIFORM, END_FRACTIONS = 13, (0.015, 0.035, 0.055)   # the sweep grid: 13 even steps, and 3 more points at these fractions of the interval from each end
MAP_N = 11                                     # the 2-D map is 11 x 11 = 121 candidates at most


def folder() -> Path:
    """out\\optimization\\sweeps_stator\\ (made if it is not there): every output of this script lives here."""
    return design.out_dir(RUN)


# ---------------------------------------------------------------------------------------------- step 1: the feasible intervals (no FEMM)
def layer_split(p) -> float:
    """r_mid - r_wall [mm]: positive means the slot cannot be cut into two coil layers of equal area inside its straight walls. motorgeo._slot
    puts the layers' boundary at r_mid (equal areas), and when r_mid is past r_wall the lower layer is drawn inside out. design.screen() sees
    that only as 'invalid shapes' with no number: it is a hidden constraint (G stays <= 0 all the way to it), so it is computed here."""
    a = math.radians(p.slot_half_deg)
    r_bottom = p.r_wall * math.tan(a)
    area = 0.5 * (p.r_wall ** 2 - p.r_shoe ** 2) * math.sin(2 * a) + r_bottom ** 2 / 2 * (math.pi + 2 * a + math.sin(2 * a))
    return math.sqrt(p.r_shoe ** 2 + area / math.sin(2 * a)) - p.r_wall


def probe(var: str, x: float) -> dict:
    """design.screen() of one value of one variable, never raising (a wide scan meets exceptions a narrow one cannot)."""
    try:
        return design.screen({var: float(x)})
    except Exception as error:
        return {"feasible": False, "reason": f"{type(error).__name__}: {error}", "G": dict.fromkeys(design.G_NAMES, 0.0), "metrics": {},
                "params": design.params({var: float(x)})}


def broken(s: dict) -> list[str]:
    """What stops a screened candidate: the constraints with G > 0, or, when the drawing failed first, the word 'drawing'."""
    names = [k for k, v in s["G"].items() if v > 0]
    if not names and s["reason"].startswith(("invalid shapes", "GeometryError")):
        names = ["drawing"]
    return names


def end_of_run(var: str, x_ok: float, x_bad: float) -> tuple[float, float]:
    """Bisect between a feasible and an infeasible value to BISECT_TOL: returns (last feasible, first infeasible)."""
    while abs(x_ok - x_bad) > BISECT_TOL:
        mid = (x_ok + x_bad) / 2
        if probe(var, mid)["feasible"]:
            x_ok = mid
        else:
            x_bad = mid
    return x_ok, x_bad


def scan_variable(var: str) -> tuple[dict, list[list]]:
    """The wide scan of one variable, then the bisection of both ends of the run that holds the reference value."""
    first, last, step = SCAN[var]
    xs = np.round(np.arange(first, last + step / 2, step), 6)
    shots = [probe(var, x) for x in xs]
    flags = [s["feasible"] for s in shots]
    runs, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        if not f and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(xs) - 1))
    ref = getattr(REF, var)
    main = next(r for r in runs if xs[r[0]] <= ref <= xs[r[1]])
    info = {"var": var, "unit": UNIT[var], "ref": ref, "scan": [first, last, step], "n_scanned": len(xs), "n_feasible": int(sum(flags)),
            "runs": [[float(xs[a]), float(xs[b])] for a, b in runs]}
    for side, i_ok, i_bad in (("lo", main[0], main[0] - 1), ("hi", main[1], main[1] + 1)):
        at_scan_limit = i_bad < 0 or i_bad >= len(xs)
        if at_scan_limit:
            x_end, x_bad, beyond = float(xs[i_ok]), None, None
        else:
            x_end, x_bad = end_of_run(var, float(xs[i_ok]), float(xs[i_bad]))
            beyond = probe(var, x_bad)
        at_end = probe(var, x_end)
        # "open": nothing in the constraint list ends the run, only the variable's own zero (a zero-width opening, a knife edge)
        is_open = at_scan_limit or (side == "lo" and x_end <= first + 1e-3)
        info[side] = {"x": x_end, "beyond_x": x_bad, "kind": "open" if is_open else "constraint",
                      "limit": [] if beyond is None else broken(beyond), "reason": "" if beyond is None else beyond["reason"],
                      "G_beyond": {} if beyond is None else {k: v for k, v in beyond["G"].items() if v > 0},
                      "G_at_end": at_end["G"], "layer_split_at_end": layer_split(at_end["params"]),
                      "tooth_min": at_end["metrics"].get("tooth_min"), "yoke": at_end["metrics"].get("yoke"),
                      "coil_layer_area_mm2": at_end["metrics"].get("coil_layer_area_mm2")}
    # coil area across the run: does the variable move the slot? (relative to the reference's)
    area = [shots[i]["metrics"]["coil_layer_area_mm2"] for i in range(main[0], main[1] + 1)]
    info["area_rel_min"], info["area_rel_max"] = min(area) / design.REFERENCE_COIL_AREA, max(area) / design.REFERENCE_COIL_AREA
    info["area_span_rel"] = info["area_rel_max"] - info["area_rel_min"]
    table = []
    for x, s in zip(xs, shots):
        m = s["metrics"]
        table.append([var, float(x), int(s["feasible"]), s["reason"][:90], *(round(s["G"][g], 5) for g in design.G_NAMES),
                      m.get("coil_layer_area_mm2", ""), m.get("tooth_min", ""), m.get("yoke", ""),
                      round(layer_split(s["params"]), 5) if s["G"]["shoe_taper"] <= 0 and s["G"]["open_vs_slot"] <= 0 else ""])
    return info, table


def stage_screen() -> None:
    """Step 1: scan every variable alone; write feasible_intervals.json and screen_scan.csv; print one block per variable."""
    started = time.perf_counter()
    out, rows = {}, []
    for var in VARS:
        t0 = time.perf_counter()
        info, table = scan_variable(var)
        out[var] = info
        rows += table
        lo, hi = info["lo"], info["hi"]
        print(f"{var} [{info['unit']}]: reference {info['ref']}; scanned {info['n_scanned']} values from {info['scan'][0]} to {info['scan'][1]} "
              f"in steps of {info['scan'][2]} ({time.perf_counter() - t0:.0f} s); feasible {info['n_feasible']}; runs {info['runs']}")
        for side, e in (("low ", lo), ("high", hi)):
            what = f"{e['kind']}, ended by {', '.join(e['limit']) or 'the scan limit'}"
            print(f"   {side} end {e['x']:.5f}  ({what}; just beyond: {e['reason'] or '-'}; G beyond {({k: round(v, 4) for k, v in e['G_beyond'].items()})})")
            print(f"            at the end: tooth_min {e['tooth_min']:.4f}, yoke {e['yoke']:.4f}, layer_split {e['layer_split_at_end']:+.5f}, "
                  f"coil layer {e['coil_layer_area_mm2']:.3f} mm2")
        print(f"   coil layer area across the interval: x{info['area_rel_min']:.3f} .. x{info['area_rel_max']:.3f} of the reference "
              f"(span {100 * info['area_span_rel']:.1f} %)")
    (folder() / "feasible_intervals.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    with open(folder() / "screen_scan.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["var", "x", "feasible", "reason", *(f"G_{g}" for g in design.G_NAMES), "coil_layer_area_mm2", "tooth_min", "yoke", "layer_split"])
        w.writerows(rows)
    print(f"{len(rows)} values screened in {time.perf_counter() - started:.0f} s; wrote feasible_intervals.json and screen_scan.csv")


def intervals() -> dict:
    """feasible_intervals.json, written by stage_screen."""
    return json.loads((folder() / "feasible_intervals.json").read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------------------------- the sweep grids
def sweep_ends(var: str, iv: dict) -> tuple[float, float]:
    """The first and last value of a sweep: the interval's ends pulled MARGIN inside their constraints; an end that no constraint ends
    (a zero-width opening, a knife-edge tip) starts at SOFT_LOW instead."""
    lo, hi = iv[var]["lo"], iv[var]["hi"]
    return (SOFT_LOW[var] if lo["kind"] == "open" else lo["x"] + MARGIN), (hi["x"] - MARGIN)


def grid_points(var: str, iv: dict) -> list[float]:
    """19 values across the sweep: 13 evenly spaced, and 3 more near each end (1.5, 3.5 and 5.5 % of the length from it), because an
    optimizer that is stopped by a constraint ends up at the end, and the last few per cent is where it would spend its time."""
    lo, hi = sweep_ends(var, iv)
    length = hi - lo
    xs = list(np.linspace(lo, hi, UNIFORM))
    for f in END_FRACTIONS:
        xs += [lo + f * length, hi - f * length]
    return sorted({round(float(x), 6) for x in xs})


# ---------------------------------------------------------------------------------------------- reading and showing rows
def stages_done(csv_name: str) -> set[str]:
    """The stage labels already in a CSV: a stage that is there is not run again (a re-run would double its rows)."""
    path = folder() / csv_name
    if not path.exists():
        return set()
    with open(path, newline="", encoding="utf-8") as f:
        return {r["stage"] for r in csv.DictReader(f)}


def read_rows(csv_name: str, stage: str | None = None) -> list[dict]:
    """The rows of a CSV as dicts of floats (strings stay strings), optionally only one stage's. Empty fields become nan."""
    with open(folder() / csv_name, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = []
    for r in rows:
        if stage is not None and r["stage"] != stage:
            continue
        d = {}
        for k, v in r.items():
            if k in ("run", "stage", "reason", "scoring", "overrides", "angle_T", "T"):
                d[k] = v
            elif k in ("feasible", "fine", "cached"):
                d[k] = v == "True"
            else:
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = float("nan")
        out.append(d)
    return out


def pcu_rel(row: dict) -> float:
    """The copper loss relative to the reference's, from an Evaluator row: (I / 10 A)^2 x (A_ref / A)."""
    return (row["amps"] / design.I_REF) ** 2 * design.REFERENCE_COIL_AREA / design.coil_area(row["metrics"])


def show(title: str, rows: list[dict], var: str | None, seconds: float) -> None:
    """One compact table per batch: the values the report quotes."""
    solved = [r for r in rows if r["feasible"] and not r.get("cached")]
    print(f"-- {title}: {len(rows)} candidates, {len(solved)} solved in {seconds:.0f} s")
    print(f"{'x':>9} {'amps':>7} {'T_mean':>8} {'rip %':>7} {'Pcu_rel':>7} {'gamma':>6} {'L_d mH':>7} {'L_q mH':>7} {'pf':>6}")
    for r in sorted(rows, key=lambda r: getattr(r["params"], var) if var else 0):
        x = getattr(r["params"], var) if var else float("nan")
        if r["feasible"]:
            print(f"{x:9.4f} {r['amps']:7.3f} {r['T_mean']:8.5f} {r['ripple_pct']:7.3f} {pcu_rel(r):7.3f} {r['gamma']:6.2f} "
                  f"{1e3 * r['L_d']:7.3f} {1e3 * r['L_q']:7.3f} {r['pf']:6.3f}")
        else:
            print(f"{x:9.4f}   not scored: {r['reason']}")


def pool():
    return femm_queue.femm_pool("optimization: stator sweeps", task_timeout_s=design.TASK_TIMEOUT_S)


# ---------------------------------------------------------------------------------------------- step 2: the one-at-a-time sweeps
def stage_sweeps(scoring: str, only: list[str] | None) -> None:
    """The one-at-a-time sweeps, one batch (one stage label) per variable. J: all five. AT: the variables whose coil area moves."""
    iv = intervals()
    if scoring == "J":
        variables, csv_name, tag = list(VARS), CSV_J, ""
    else:
        variables = [v for v in VARS if iv[v]["area_span_rel"] > AT_AREA_SPAN]
        csv_name, tag = CSV_AT, "_AT"
        print(f"AT sweeps for the variables whose coil area moves by more than {100 * AT_AREA_SPAN:.0f} %: {variables}")
    done = stages_done(csv_name)
    todo = [v for v in variables if v + tag not in done and (not only or v in only)]
    if not todo:
        print("nothing to run: every stage asked for is in the CSV already")
        return
    started = time.perf_counter()
    with pool() as p:
        ev = design.Evaluator(p, run=RUN, scoring=scoring, csv_name=csv_name)
        for v in todo:
            candidates = [{}] + [{v: x} for x in grid_points(v, iv)]
            t0 = time.perf_counter()
            rows = ev(candidates, stage=v + tag)
            show(f"{v}{tag} (scoring {scoring})", rows, v, time.perf_counter() - t0)
    print(f"{len(todo)} sweeps in {time.perf_counter() - started:.0f} s")


# ---------------------------------------------------------------------------------------------- step 3: noise, and coarse against fine
def stage_noise(steps: list[float], stage: str) -> None:
    """The reference with each variable moved by +-each step (micrometres: far too small to change the physics): what the scores
    still differ by is the remeshing noise of the coarse mesh. One batch, J scoring (the reference's J is its AT)."""
    if stage in stages_done(CSV_J):
        print(f"stage {stage} is in the CSV already")
        return
    candidates = [{}]
    for v in VARS:
        for h in steps:
            for sign in (+1, -1):
                candidates.append({v: round(getattr(REF, v) + sign * h, 6)})
    t0 = time.perf_counter()
    with pool() as p:
        ev = design.Evaluator(p, run=RUN, scoring="J", csv_name=CSV_J)
        rows = ev(candidates, stage=stage)
    for v in VARS:
        mine = [r for r in rows if r["overrides"] == {} or v in r["overrides"]]
        show(f"{stage} {v}", mine, v, 0.0)
    print(f"{len(candidates)} candidates in {time.perf_counter() - t0:.0f} s")


def fine_points(var: str, neighbours: bool) -> list[float]:
    """The sweep's points to re-score on the fine mesh: both ends, the lowest ripple, the highest torque, and (neighbours=True) the grid points
    on each side of the lowest ripple, to see whether the best point stays where it is."""
    rows = [r for r in read_rows(CSV_J, var) if r["feasible"]]
    rows = sorted(rows, key=lambda r: r[var])
    best_r = min(range(len(rows)), key=lambda i: rows[i]["ripple_pct"])
    best_t = max(range(len(rows)), key=lambda i: rows[i]["T_mean"])
    pick = {0, len(rows) - 1, best_r, best_t}
    if neighbours:
        pick |= {max(0, best_r - 1), min(len(rows) - 1, best_r + 1)}
    return [rows[i][var] for i in sorted(pick)]


def stage_fine(neighbours: bool, extra: dict[str, list[float]]) -> None:
    """Re-score points of every J sweep on the fine mesh (its own Evaluator and CSV), one batch. extra: more values to re-score,
    {variable: [value, ...]} (a step in a coarse sweep is worth a look on the fine mesh: is it the motor or the mesh?)."""
    if "fine" in stages_done(CSV_FINE):
        print("stage fine is in the CSV already")
        return
    candidates = [{}]
    for v in VARS:
        pts = sorted(set(fine_points(v, neighbours) + extra.get(v, [])))
        print(f"fine {v}: {pts}")
        candidates += [{v: x} for x in pts if abs(x - getattr(REF, v)) > 1e-9]
    t0 = time.perf_counter()
    with pool() as p:
        ev = design.Evaluator(p, run=RUN, scoring="J", fine=True, csv_name=CSV_FINE)
        rows = ev(candidates, stage="fine")
    for v in VARS:
        mine = [r for r in rows if r["overrides"] == {} or v in r["overrides"]]
        show(f"fine {v}", mine, v, 0.0)
    print(f"{len(candidates)} candidates in {time.perf_counter() - t0:.0f} s")


# ---------------------------------------------------------------------------------------------- step 4: the 2-D map
def stage_map(vx: str, vy: str, box: list[float] | None) -> None:
    """MAP_N x MAP_N values of two variables (J scoring), one batch. The cells the screen refuses are written as rows too (feasible False,
    with the reason), so the map can show where the room ends. box = lo_x hi_x lo_y hi_y, by default both variables' sweep ends."""
    if "map" in stages_done(CSV_J):
        print("stage map is in the CSV already")
        return
    iv = intervals()
    lo_x, hi_x = sweep_ends(vx, iv) if not box else box[:2]
    lo_y, hi_y = sweep_ends(vy, iv) if not box else box[2:]
    xs, ys = np.linspace(lo_x, hi_x, MAP_N), np.linspace(lo_y, hi_y, MAP_N)
    candidates = [{vx: round(float(x), 6), vy: round(float(y), 6)} for y in ys for x in xs]
    n_ok = sum(design.screen(c)["feasible"] for c in candidates)
    print(f"map {vx} [{lo_x:.4f}, {hi_x:.4f}] x {vy} [{lo_y:.4f}, {hi_y:.4f}]: {len(candidates)} cells, {n_ok} pass the screen")
    t0 = time.perf_counter()
    with pool() as p:
        ev = design.Evaluator(p, run=RUN, scoring="J", csv_name=CSV_J)
        rows = ev(candidates, stage="map")
    best = min((r for r in rows if r["feasible"]), key=lambda r: r["ripple_pct"])
    print(f"{n_ok} solved in {time.perf_counter() - t0:.0f} s; lowest ripple {best['ripple_pct']:.3f} % at {vx} {getattr(best['params'], vx):.4f}, "
          f"{vy} {getattr(best['params'], vy):.4f}")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="stator-tier sweeps of the geometry study")
    sub = ap.add_subparsers(dest="stage", required=True)
    sub.add_parser("screen")
    p = sub.add_parser("j")
    p.add_argument("only", nargs="*")
    p = sub.add_parser("at")
    p.add_argument("only", nargs="*")
    p = sub.add_parser("noise")
    p.add_argument("--steps", default="0.001,0.002,0.003")
    p.add_argument("--label", default="noise")
    p = sub.add_parser("fine")
    p.add_argument("--neighbours", action="store_true")
    p.add_argument("--extra", action="append", default=[], help="variable=value,value: more points for the fine mesh (repeatable)")
    p = sub.add_parser("map")
    p.add_argument("vx")
    p.add_argument("vy")
    p.add_argument("box", nargs="*", type=float)
    a = ap.parse_args(argv)
    if a.stage == "screen":
        stage_screen()
    elif a.stage == "j":
        stage_sweeps("J", a.only)
    elif a.stage == "at":
        stage_sweeps("AT", a.only)
    elif a.stage == "noise":
        stage_noise([float(s) for s in a.steps.split(",")], a.label)
    elif a.stage == "fine":
        stage_fine(a.neighbours, {e.split("=")[0]: [float(x) for x in e.split("=")[1].split(",")] for e in a.extra})
    elif a.stage == "map":
        stage_map(a.vx, a.vy, a.box or None)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
