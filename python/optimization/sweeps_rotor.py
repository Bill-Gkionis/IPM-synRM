"""optimization.sweeps_rotor: the rotor tier swept one variable at a time across everything the black box allows, so that the optimizer is chosen after looking at the landscape.

What it is, and why
    An optimizer is a bet on the shape of what it searches: a gradient method bets that the score is smooth, an evolutionary one that it has several valleys, a
    constrained one that the best design sits on a limit. This script looks first. It moves the rotor tier's variables one at a time through the black box
    (optimization.design; the stator, the winding and the magnets' size stay at the reference). B7 had shown that within +-1.5 % of the reference the mean torque
    follows the magnet flux and each dimension acts almost linearly; nobody had looked farther out, or at the ripple.
      mag_x       the magnet slot's centre along the d-axis [mm] (reference 39.4136): with mag_y, the radius of the V
      mag_y       the slot's distance from the d-axis [mm] (6.0): half the spacing of a pole's two magnets
      v_half_deg  the magnet's long side against the d-axis [deg] (60.0)
      bridge_mm   the mapped form of mag_x (design.decode): the outer bridge's thickness [mm] (reference 1.5585). The build guard (a bridge of at least 1.4 mm,
                  decision 17, a fixed design decision, never relaxed here) is then a plain lower bound, and "on the guard" means bridge_mm = 1.4.
    Six things come out of it: where each variable's feasible interval ends and which constraint ends it; whether the torque peaks inside the interval or
    keeps rising to the guard, and why; the shape of the torque and of the ripple along each variable; the noise floor of the scoring (is a gradient usable?);
    whether the coarse mesh orders designs as the fine one does; and, on one 2-D map, how the two variables that matter for the ripple interact.

The stages (one command each, below). Every FEMM stage opens its own femm_pool on the machine-wide solver slot: while another run holds the slot it prints
"waiting for the solver slot" and waits. Every stage is also safe to run twice: a candidate already in the CSV under the stage's label is not scored again.
    scan       step 1, Python only: design.screen() on a fine grid of each variable, the unbroken feasible runs, each end found by bisection and named by the
               constraint that is broken just beyond it (what_ends_it). Writes scan.json and scan_<name>.csv (every verdict, so anyone can look).
    sweeps     step 2: the reference (checked against the numbers the brief gives), then six sweeps of about 20 to 30 points: mag_x, mag_y and v_half_deg from the
               reference; mag_y and v_half_deg on the guard; and bridge_mm from the guard to the deep end of its interval (the radius sweep).
    noise      step 4a: the reference moved by hairs (+-0.001 to +-0.05 mm, +-0.001 to +-0.3 deg), one variable at a time. (sweeps --with-noise runs both in one session.)
    batch      the radius sweep again at four other settings of mag_y and v_half_deg (CROSS), the fine-mesh re-scores (run_fine: both ends, the best torque, the
               lowest ripple and the middle of every sweep) and the 2-D map (MAP), in one session: one wait for the slot, not three.
    valley     a line through the map's deepest valley (VALLEY), between its grid points. fine-map: the map's own front on the fine mesh.
    strip      the feasible room under the map's lowest row (STRIP): the centre rib's limit falls from 4.85 mm at 60 degrees to 3.48 at 30, so the map, which starts
               at mag_y 4.9, leaves a strip at the small angles; 18 points on it, at the rib's own limit and 0.4 and 0.8 mm above.
    analyse    steps 3 and 4 and the map's reading, from the CSVs alone (no FEMM): tables.txt and summary.json. sweeps_rotor_figures.py draws the pictures.

Method
    Feasible intervals. design.screen(candidate) says feasible or not and which G values are positive (G <= 0 is allowed). On a grid the unbroken runs of
    "feasible" are the intervals; each end is refined by 45 halvings between the last feasible grid point and the first infeasible one, and named by the G values
    that are positive 1e-4 beyond it. A scan limit is not a constraint: the report says so where an interval simply reaches the end of what was scanned.
    The guard trap. design.mag_x_for_bridge() bisects 60 times and returns the middle of the last interval, which lies on either side of the root by about 1e-15
    mm. At bridge_mm = 1.4 exactly, screen() then rejects about half the candidates (guard_bridge = +1.3e-15) although motorgeo.build accepts them: the
    scan "guard_mag_y_exact" (161 points) counts 74 of the 142 that would be feasible, "guard_v_half_deg_exact" 122 of 237. This script sits on the guard at
    1.4 + 1e-9 mm (GUARD_MM), where none is rejected, and the sweeps' ends sit END_MARGIN = 1e-5 inside a constraint's edge for the same reason.
    Sweep points. n points from lo to hi, x_k = lo + (hi - lo) [b (1 - cos(pi u_k))/2 + (1 - b) u_k], u_k = k/(n - 1), b = 0.7: a mix of a cosine spacing (dense at
    both ends, where the limits are) and an even one (no gap in the middle for a narrow valley to hide in). The reference's value is always one of the points.
    The torque split. The torque from the flux linkages is T = k (psi_d i_q - psi_q i_d); with psi_d = psi_pm + L_d i_d and psi_q = L_q i_q it is
    T = k psi_pm i_q + k (L_d - L_q) i_d i_q, the magnet's torque and the reluctance torque. k is read from the rows, not assumed: it is 6.0000 (3/2 x 4 pole
    pairs) to 6e-12 %. psi_pm is the no-load d-axis flux linkage, L_d and L_q are the apparent inductances of design.score().
    The leakage law. A steel strip beside a magnet that the magnet drives into saturation carries about B_sat x its thickness x the stack and no more, and that
    flux does not reach the air gap, so what is lost is linear in the strip's thickness: psi_pm = a - b t_bridge - c t_rib (leakage_fit; b for the two outer
    bridges of a V together, c for the centre rib). It is fitted on the sweeps from the reference (bridge up to 3.5 mm, rib up to 4.7 mm) and tested on every other row.
    The MTPA window. design.score() fits a parabola through the torque at 10, 25 and 40 degrees and keeps its top inside them: a design whose best angle is
    farther out is scored at the window's edge and its torque is a little low. 36 of the 361 scored rows are at an edge (the tables mark them *); every one has
    T <= 0.42 N.m, and the 309 rows at or above 0.5 N.m have angles from 14.9 to 36.3 deg. The ripple in per cent of such a torque says little.
    The noise floor. Scoring the same design twice gives the same numbers, but a design a hair away is meshed again, and the mesh changes the answer a little. Over
    the 7 points |h| <= 0.003 (the reference and its six perturbations) the true response is a straight line to far better than the noise, so what a fitted line
    misses is noise: sigma = sqrt(sum of squared misses / (n - 2)), pooled over the three variables. A central difference of step h has the error
    sigma / (sqrt(2) h) on a gradient b: it is usable (under 10 %) from h_10 = sigma / (0.1 sqrt(2) |b|) = 7.07 sigma / |b|.
    Shape. A turning point counts only if its prominence (scipy.signal.find_peaks: how far the curve falls from it, on both sides, before it reaches a point as good) is at least 3 sigma. The
    roughness is the rms miss of each interior point from the cubic through its two neighbours on each side, over what noise alone gives there: about 1 for a
    curve smooth to the noise at the spacing of its points, large where it has a valley or a kink at that scale.
    Coarse against fine. Of every two re-scored points whose coarse values differ by more than 3 sigma, how many does the fine mesh order the other way round
    (discordant_pairs)? Pairs closer than that are not an order the coarse mesh claims.
    The map. 8-neighbour local minima (local_minima_2d): a grid point lower than all its neighbours by 3 sigma. The map's front is design.nondominated() on
    (-T_mean, ripple_pct); its lowest ripple among the points with at least 80 % or 90 % of its best torque is the compromise usually asked for.

Run (from the Python to FEMM folder; about 10 minutes of one core for the scan, about 8 minutes of pool time for the rest, on a free slot)
    .venv/Scripts/python.exe -m optimization.sweeps_rotor scan
    .venv/Scripts/python.exe -m optimization.sweeps_rotor sweeps --with-noise
    .venv/Scripts/python.exe -m optimization.sweeps_rotor batch
    .venv/Scripts/python.exe -m optimization.sweeps_rotor valley
    .venv/Scripts/python.exe -m optimization.sweeps_rotor strip
    .venv/Scripts/python.exe -m optimization.sweeps_rotor fine-map
    .venv/Scripts/python.exe -m optimization.sweeps_rotor analyse          and then      -m optimization.sweeps_rotor_figures
Writes out\\optimization\\sweeps_rotor\\: evaluations.csv (every coarse row, the schema of design.FIELDS, with the stage label), evaluations_fine.csv, scan.json,
scan_<name>.csv, map.json, summary.json, tables.txt, the logs and the figures.

Numbers to expect (this PC, 2026-11-08, FEMM 5.0, coarse mesh unless said, scoring AT: 100 ampere-turns at each design's own MTPA angle)
    the reference  gamma 22.337 deg, T_mean 1.1217 N.m, ripple 0.22008 N.m (19.62 %); fine mesh 22.46 deg, 1.12880 N.m, 18.45 %: the brief's numbers
    size           367 coarse candidates (361 scored: 6 cells of the map are infeasible and cost no solve) and 30 fine ones; 6 pool sessions that solved something
    intervals      mag_x 24.6925 to 39.5760 (q_rib | guard_bridge; the reference is 0.162 below the top) and a spurious island 0 to 4.278 that screen() accepts
                   (no shaft check); mag_y 4.8511 to 6.6978 (centre_rib | guard_bridge); v_half_deg 56.6399 to 89 (guard_bridge | the end of the scan; nothing ends
                   it up to 135 degrees, the V turned inside out); bridge_mm 1.4 to 15.7401 (the guard | q_rib); on the guard mag_y 4.8511 to 11.4943 (centre_rib
                   | q_rib) and v_half_deg 0.5 to 135 with no end
    torque         rises to the guard along mag_x and bridge_mm (no turning point beyond 3 sigma; T 1.1217 at the reference, 1.1726 at bridge 1.4, 0.0016 at 15.7 mm);
                   rises as mag_y falls to the centre rib's limit (0.966 to 1.398 N.m from the guard end to the rib end; 1.535 on the guard); on the guard it rises with
                   v_half_deg to a plateau (1.173 at 60 deg, 1.354 at 86.7, 1.353 at 89: the 0.9 mN.m gain over the end is 1.6 sigma); from the reference, with mag_x
                   fixed, it falls with v_half_deg (1.129 at 56.64 to 0.743 at 89) because the bridge thickens from 1.4 to 3.3 mm
    the why        psi_pm = 33.71 - 5.152 t_bridge - 2.765 t_rib [mWb, mm]: R^2 0.9992 on the 50 fitting rows (rms miss 0.094 mWb), rms miss 0.63 mWb (3.6 %) on 168
                   others. From the reference's bridge 1.5585 to the guard, T +0.0522 N.m = magnet torque +0.0550 + reluctance torque -0.0028; psi_pm +0.925 mWb,
                   L_d -0.2 %, L_q +0.13 %, saliency 1.828 to 1.834, MTPA angle 22.32 to 21.78 deg. Slope at the guard -0.323 N.m per mm of bridge; at four other
                   settings -0.294 to -0.346: the answer does not depend on where mag_y and v_half_deg sit
    ripple         in per cent of the mean. Along the radius it only falls as the bridge thickens (20.97 % at the guard, 7.4 % at 3.5 mm); along mag_y and v_half_deg it
                   has valleys: on the guard mag_y has three interior valleys 1.8 to 1.9 mm apart (6.53, 8.42, 10.22) and v_half_deg two (36.1, 74.9)
    map            mag_y 4.9 to 7.5 against v_half_deg 30 to 89 on the guard: 111 of 117 points scored, T 0.361 to 1.630 N.m, ripple 6.01 to 25.03 %; the valleys run
                   as parallel diagonal bands (6 local minima); the front has 4 points, from T 1.630 at 19.55 % to T 1.435 at 6.01 % (54.6 deg, mag_y 4.9, on the
                   centre rib's limit); the line through it puts the bottom at 6.09 % at 54.0 deg (T 1.4248); on the fine mesh 5.56 % at T 1.4422. The strip under the
                   map (18 points): T 1.094 to 1.425, ripple 10.9 to 39.7 %, nothing that beats the front (the lowest ripple there is 10.86 % at T 1.371, 39 deg)
    noise          sigma 0.559 mN.m on the torque (0.05 %), 0.144 points on the ripple, 0.019 deg on the angle; the +-0.003 set spans 1.8 to 2.5 mN.m and 0.33 to 0.46
                   points; the exact reference sits 0.8 to 1.45 mN.m above the line through its six neighbours. The sweeps are 691 to 2612 sigma wide in torque and
                   79 to 167 sigma in ripple. Usable central-difference steps at the reference: torque 0.011 (mag_x), 0.016 (mag_y) mm, 1.1 deg (v_half_deg); ripple 0.19, 0.072
                   mm, 0.80 deg; on the guard 0.013 mm and 0.32 deg for the torque, 0.074 mm and 1.6 deg for the ripple
    fine mesh      21 designs of at least 0.5 N.m: T +0.60 % (+0.46 to +0.94), ripple -0.65 points (-1.27 to +0.34), MTPA angle +0.12 deg; discordant pairs: T 0 of 205,
                   ripple 1 of 197

Self-check (the answers follow).
  1. At bridge_mm = 1.4 exactly, half the designs of a scan on the guard are infeasible, yet motorgeo builds them. Why, and why does 1.4 + 1e-9 not have the problem?
  2. The ripple in per cent has its lowest value of the radius sweep, 1.46 %, at a bridge of 10 mm. Is that a finding about the motor?
  3. A central difference of the torque in mag_y with h = 0.003 mm: usable?
Answers.
  1. mag_x_for_bridge() returns the middle of its last bisection interval, so the mag_x it gives leaves a bridge that is 1e-15 mm over or under 1.4. Under, screen()
     computes guard_bridge = 1.4 - bridge = +1e-15 > 0 and calls the design infeasible; motorgeo.build compares the same sum another way round and accepts it.
     1e-9 mm is far above the bisection's 1e-15 and far below anything a drawing can show, so the sign of the margin no longer depends on the rounding.
  2. No. At 10 mm the torque is 0.0047 N.m (the magnets' flux is shunted away: psi_pm 0.066 mWb), so 1.46 % is a ripple of 69 micro-newton-metres over a mean of
     4.7 milli-newton-metres: a ratio of two numbers that are nearly nothing. An optimizer that minimises it will find such designs unless a torque floor stops it.
  3. No. The slope is -0.2115 N.m per mm, so the two points differ by 1.3 mN.m; the noise of their difference is sqrt(2) x 0.46 = 0.66 mN.m: a signal to noise of
     2 (the table's 2.1), and an error of about 50 % on the gradient. The step must be 0.0155 mm or more for 10 %.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

import femm_queue
from optimization import design

RUN = "sweeps_rotor"
REF = design.REFERENCE
GUARD_EXACT_MM = REF.bridge                 # 1.4 mm: decision 17, a manufacturing decision (never relaxed here)
GUARD_MM = GUARD_EXACT_MM + 1e-9            # "on the guard": see the note in the docstring (mag_x_for_bridge's rounding)
END_MARGIN = 1e-5                           # a sweep's end sits this far inside the edge, so that rounding cannot flip it


def folder() -> Path:
    """out\\optimization\\sweeps_rotor\\ (made if it is not there): every output of this script lives here."""
    return design.out_dir(RUN)


# ---------------------------------------------------------------------------------------------- step 1: the feasible intervals
# Each scan is one variable alone, the others at the reference (or, "guard_*", on the guard through bridge_mm). make(x) is the candidate.
SCANS = {
    "mag_x":            {"make": lambda x: {"mag_x": x}, "lo": 20.0, "hi": 40.0, "step": 0.01, "ref": REF.mag_x},
    "mag_x_low":        {"make": lambda x: {"mag_x": x}, "lo": 0.0, "hi": 20.0, "step": 0.05, "ref": None},
    "mag_y":            {"make": lambda x: {"mag_y": x}, "lo": 0.0, "hi": 16.0, "step": 0.01, "ref": REF.mag_y},
    "v_half_deg":       {"make": lambda x: {"v_half_deg": x}, "lo": 30.0, "hi": 89.0, "step": 0.05, "ref": REF.v_half_deg},
    "v_half_deg_over":  {"make": lambda x: {"v_half_deg": x}, "lo": 89.0, "hi": 135.0, "step": 0.1, "ref": None},
    "bridge":           {"make": lambda x: {"bridge_mm": x}, "lo": GUARD_MM, "hi": 20.0, "step": 0.01, "ref": design.bridge_of(REF)},
    "guard_mag_y":      {"make": lambda x: {"bridge_mm": GUARD_MM, "mag_y": x}, "lo": 0.0, "hi": 16.0, "step": 0.01, "ref": REF.mag_y},
    "guard_v_half_deg": {"make": lambda x: {"bridge_mm": GUARD_MM, "v_half_deg": x}, "lo": 30.0, "hi": 89.0, "step": 0.05, "ref": REF.v_half_deg},
    "guard_v_half_deg_low":  {"make": lambda x: {"bridge_mm": GUARD_MM, "v_half_deg": x}, "lo": 0.5, "hi": 30.0, "step": 0.1, "ref": None},
    "guard_v_half_deg_over": {"make": lambda x: {"bridge_mm": GUARD_MM, "v_half_deg": x}, "lo": 89.0, "hi": 135.0, "step": 0.1, "ref": None},
    # the same two on the guard with bridge_mm = 1.4 exactly, on a coarser grid: only to count what the rounding rejects
    "guard_mag_y_exact":      {"make": lambda x: {"bridge_mm": GUARD_EXACT_MM, "mag_y": x}, "lo": 4.0, "hi": 12.0, "step": 0.05, "ref": None, "count_only": True},
    "guard_v_half_deg_exact": {"make": lambda x: {"bridge_mm": GUARD_EXACT_MM, "v_half_deg": x}, "lo": 30.0, "hi": 89.0, "step": 0.25, "ref": None, "count_only": True},
}
SCAN_ORDER = list(SCANS)


def probe(candidate: dict) -> dict:
    """screen() of one candidate as a plain dict {"feasible", "reason", "G"}; it never raises. screen() catches only GeometryError, and a
    scan this wide can meet other exceptions (a square root of a negative number, say): a scan that stops at the first one finds nothing."""
    try:
        s = design.screen(design.decode(candidate))
    except Exception as error:
        return {"feasible": False, "reason": f"{type(error).__name__}: {error}", "G": {}}
    return {"feasible": bool(s["feasible"]), "reason": s["reason"], "G": s["G"]}


def positive(G: dict) -> list[str]:
    """The constraints that are broken (G > 0) in a G dict, in the order of design.G_NAMES."""
    return [g for g in design.G_NAMES if G.get(g, 0.0) > 0]


def grid(lo: float, hi: float, step: float) -> np.ndarray:
    """lo, lo + step, ... up to hi, computed as lo + k step (not by adding step again and again, which drifts)."""
    return lo + step * np.arange(int(round((hi - lo) / step)) + 1)


def contiguous_runs(flags) -> list[tuple[int, int]]:
    """[(first index, last index)] of every unbroken run of True in a list of booleans."""
    out, start = [], None
    for i, f in enumerate(list(flags) + [False]):
        if f and start is None:
            start = i
        if not f and start is not None:
            out.append((start, i - 1))
            start = None
    return out


def edge(make, x_in: float, x_out: float, rounds: int = 45) -> float:
    """The feasible end of an interval, found by bisection between a feasible x_in and an infeasible x_out (adjacent scan points). 45 halvings
    of a 0.01 step are far below a nanometre: the end is the constraint's own root, not a multiple of the scan step."""
    for _ in range(rounds):
        mid = (x_in + x_out) / 2
        if probe(make(mid))["feasible"]:
            x_in = mid
        else:
            x_out = mid
    return x_in


def what_ends_it(make, x_edge: float, direction: int) -> list[str]:
    """The constraint(s) broken just beyond an interval's end (1e-4 outward): their names, or the reason when none of the G values is positive."""
    r = probe(make(x_edge + direction * 1e-4))
    return positive(r["G"]) or [r["reason"] or "none"]


def scan_one(name: str, log) -> dict:
    """One scan: screen() at every grid point, the unbroken feasible runs, each end refined by bisection and labelled with the constraint that ends it.
    Returns the summary (what scan.json holds) and writes scan_<name>.csv (every point, so that anyone can look at the raw verdicts)."""
    spec = SCANS[name]
    make = spec["make"]
    xs = grid(spec["lo"], spec["hi"], spec["step"])
    started = time.perf_counter()
    rows = [probe(make(float(x))) for x in xs]
    flags = [r["feasible"] for r in rows]
    with open(folder() / f"scan_{name}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["x", "feasible", "reason", *(f"G_{g}" for g in design.G_NAMES)])
        for x, r in zip(xs, rows):
            w.writerow([f"{x:.6f}", r["feasible"], r["reason"], *(f"{r['G'].get(g, 0.0):.6g}" for g in design.G_NAMES)])
    intervals = []
    for i, j in ([] if spec.get("count_only") else contiguous_runs(flags)):
        lo_limit = hi_limit = "scan limit"
        lo, hi = float(xs[i]), float(xs[j])
        if i > 0:
            lo = edge(make, float(xs[i]), float(xs[i - 1]))
            lo_limit = ", ".join(what_ends_it(make, lo, -1))
        if j < len(xs) - 1:
            hi = edge(make, float(xs[j]), float(xs[j + 1]))
            hi_limit = ", ".join(what_ends_it(make, hi, +1))
        ref = spec["ref"]
        intervals.append({"lo": lo, "hi": hi, "lo_limit": lo_limit, "hi_limit": hi_limit,
                          "contains_reference": bool(ref is not None and lo - 1e-9 <= ref <= hi + 1e-9), "points": j - i + 1})
    # points that only the rounding of the guard rejects: infeasible, and the only broken constraint is guard_bridge by less than 1e-9
    flipped = sum(1 for r in rows if not r["feasible"] and positive(r["G"]) == ["guard_bridge"] and r["G"]["guard_bridge"] < 1e-9)
    out = {"scan": [float(xs[0]), float(xs[-1]), spec["step"]], "points": len(xs), "feasible_points": int(sum(flags)),
           "rejected_by_rounding": int(flipped), "intervals": intervals, "seconds": round(time.perf_counter() - started, 1)}
    log(f"scan {name}: {len(xs)} points {xs[0]:g}..{xs[-1]:g} step {spec['step']:g}, {sum(flags)} feasible, {flipped} rejected by guard rounding only, "
        f"{time.perf_counter() - started:.0f} s")
    for iv in intervals:
        log(f"    feasible {iv['lo']:.6f} .. {iv['hi']:.6f}  [{iv['lo_limit']} | {iv['hi_limit']}]{'  <- contains the reference' if iv['contains_reference'] else ''}")
    return out


def stage_scan(args) -> None:
    """Step 1: every scan in SCANS (or the ones named), Python only: no FEMM, no solver slot. Writes scan.json."""
    names = args.only or SCAN_ORDER
    path = folder() / "scan.json"
    summary = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    for name in names:
        summary[name] = scan_one(name, say)
        path.write_text(json.dumps(summary, indent=1), encoding="utf-8")


def say(text: str) -> None:
    """A printed line that reaches a log file at once (a background run is read while it works)."""
    print(text, flush=True)


# ---------------------------------------------------------------------------------------------- the FEMM stages: plumbing
def open_evaluator(pool, fine: bool = False):
    """An Evaluator that writes evaluations.csv (or evaluations_fine.csv). Every stage is a separate run of this script, so the row counter n
    would start again at 1 each time: it is set to the number of rows the file already holds, and goes on counting."""
    ev = design.Evaluator(pool, run=RUN, scoring="AT", fine=fine, csv_name="evaluations_fine.csv" if fine else "evaluations.csv")
    with open(ev.csv, newline="", encoding="utf-8") as f:
        ev.n = max(sum(1 for _ in f) - 1, 0)
    return ev


def held_in(csv_path: Path, stage: str) -> set[str]:
    """The overrides (as the JSON text the Evaluator writes) of the rows this CSV already holds under this stage label."""
    with open(csv_path, newline="", encoding="utf-8") as f:
        return {r["overrides"] for r in csv.DictReader(f) if r["stage"] == stage}


def run_stage(ev, stage: str, candidates: list[dict]) -> list[dict]:
    """Score the candidates under this stage label. Candidates the CSV already holds for the label are skipped, so a run that was stopped
    (or a stage run twice) goes on where it stopped and never writes a design twice under one label. Prints one line and the failures."""
    held = held_in(ev.csv, stage)
    todo = [c for c in candidates if json.dumps(design.decode(c)) not in held]
    if not todo:
        say(f"stage {stage}: all {len(candidates)} candidates are in {ev.csv.name} already")
        return []
    started = time.perf_counter()
    rows = ev(todo, stage=stage)
    solved = [r for r in rows if r["feasible"] and not r.get("cached")]
    failed = [r for r in rows if not r["feasible"]]
    line = f"stage {stage}: {len(todo)} candidates, {len(solved)} solved, {len(failed)} not scored, {time.perf_counter() - started:.0f} s"
    if solved:
        T, rip = [r["T_mean"] for r in solved], [r["ripple_pct"] for r in solved]
        line += f"; T_mean {min(T):.4f}..{max(T):.4f} N.m, ripple {min(rip):.2f}..{max(rip):.2f} %"
    say(line)
    for r in failed:
        say(f"    NOT SCORED {r['overrides']}: {r['reason']}")
    return rows


def pool_session(label: str):
    """The pool of this script: 8 workers on the solver slot. Waiting for the slot (another run holds it) is expected and prints a line every 30 s."""
    return femm_queue.femm_pool(f"optimization: {label}", task_timeout_s=design.TASK_TIMEOUT_S)


# ---------------------------------------------------------------------------------------------- step 2: the one-at-a-time sweeps
# What each sweep varies (var), what it holds apart from that (hold: entries of every candidate), the point that is always included (anchor:
# the reference's value), and how many points it gets: about 20, a few more where the interval is long. "guard_*" hold bridge_mm on the guard,
# so the mapping layer moves mag_x to keep the outer bridge at its minimum while the other variable moves. "bridge" is the radius sweep.
SWEEPS = {
    "mag_x":            {"var": "mag_x", "hold": {}, "n": 22, "anchor": REF.mag_x},
    "mag_y":            {"var": "mag_y", "hold": {}, "n": 20, "anchor": REF.mag_y},
    "v_half_deg":       {"var": "v_half_deg", "hold": {}, "n": 24, "anchor": REF.v_half_deg},
    "guard_mag_y":      {"var": "mag_y", "hold": {"bridge_mm": GUARD_MM}, "n": 20, "anchor": REF.mag_y},
    "guard_v_half_deg": {"var": "v_half_deg", "hold": {"bridge_mm": GUARD_MM}, "n": 24, "anchor": REF.v_half_deg},
    "bridge":           {"var": "bridge_mm", "hold": {}, "n": 15, "anchor": design.bridge_of(REF),
                         "extra": (1.45, 1.5, 1.6, 1.7, 1.8, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0)},   # more points where the question lives (the guard end) and where the torque falls away
}
SWEEP_ORDER = list(SWEEPS)
EXPECTED_COARSE = {"gamma": 22.34, "T_mean": 1.12173, "ripple_pp": 0.2201, "ripple_pct": 19.62}   # the reference, coarse mesh, as the brief gives it


def spaced(lo: float, hi: float, n: int, blend: float = 0.7) -> np.ndarray:
    """n points from lo to hi, denser toward both ends: a mix (blend) of a cosine spacing and an even one. The cosine alone would put its
    points 0.7 % of the length apart at the ends and 8 % apart in the middle (20 points); the mix keeps the ends dense (about 2 %) without leaving a
    gap in the middle that a narrow ripple valley could hide in."""
    u = np.linspace(0.0, 1.0, n)
    return lo + (hi - lo) * (blend * (1.0 - np.cos(np.pi * u)) / 2.0 + (1.0 - blend) * u)


def sweep_interval(scan: dict, name: str) -> tuple[float, float]:
    """The interval a sweep covers: the feasible interval of its scan that holds the reference, pulled END_MARGIN inside each end that a
    constraint makes (so the end point is feasible whatever the rounding does), but not at an end that is only the scan's own limit."""
    iv = next(i for i in scan[name]["intervals"] if i["contains_reference"])      # a scan and its sweep share a name
    lo = iv["lo"] if iv["lo_limit"] == "scan limit" else iv["lo"] + END_MARGIN
    hi = iv["hi"] if iv["hi_limit"] == "scan limit" else iv["hi"] - END_MARGIN
    return lo, hi


def sweep_points(name: str, lo: float, hi: float) -> list[float]:
    """The values of a sweep: spaced(lo, hi, n), the reference's value (an even point near it gives way), and the sweep's extra points."""
    spec = SWEEPS[name]
    pts = [float(x) for x in spaced(lo, hi, spec["n"])]
    gap = 0.3 * (hi - lo) / spec["n"]
    if lo < spec["anchor"] < hi:
        pts = [x for x in pts if x in (lo, hi) or abs(x - spec["anchor"]) > gap] + [float(spec["anchor"])]
    pts += [x for x in spec.get("extra", ()) if lo < x < hi]
    return sorted({round(x, 9) for x in pts})


def sweep_candidates(name: str, scan: dict) -> list[dict]:
    """The candidates of a sweep, in the order of the variable: {hold..., var: x}."""
    spec = SWEEPS[name]
    lo, hi = sweep_interval(scan, name)
    return [{**spec["hold"], spec["var"]: x} for x in sweep_points(name, lo, hi)]


def stage_sweeps(args) -> None:
    """Step 2: the reference first (a check of the whole chain against the numbers the brief gives), then every sweep of SWEEP_ORDER (or the
    ones named), all in one pool session. Needs scan.json."""
    scan = json.loads((folder() / "scan.json").read_text(encoding="utf-8"))
    names = args.only or SWEEP_ORDER
    with pool_session("rotor sweeps") as pool:
        ev = open_evaluator(pool)
        ref = run_stage(ev, "reference", [{}])
        if ref:
            r = ref[0]
            say("    the reference, coarse mesh: " + ", ".join(f"{k} {r[k]:.5g} (expected {v:g})" for k, v in EXPECTED_COARSE.items()))
        for name in names:
            run_stage(ev, name, sweep_candidates(name, scan))
        if args.with_noise:                       # one slot wait for both: the slot may be held by another sweep for minutes
            run_stage(ev, "noise", noise_candidates())


def parse_hold(items) -> dict:
    """['mag_y=5.2', 'bridge_mm=guard'] -> {'mag_y': 5.2, 'bridge_mm': GUARD_MM}: what a map given on the command line holds fixed ('guard' is bridge_mm on the guard)."""
    out = {}
    for item in items or []:
        k, v = item.split("=")
        out[k] = GUARD_MM if v == "guard" else float(v)
    return out


# ---------------------------------------------------------------------------------------------- step 4a: the noise floor
# Perturbations of the reference, in each variable's own unit (mm, mm, degrees). The three smallest are the brief's; the larger ones show where
# the signal (the real change of torque with the variable) climbs out of the noise, which is what a finite-difference step has to clear.
NOISE_H = {"mag_x": (0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05),
           "mag_y": (0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05),
           "v_half_deg": (0.001, 0.002, 0.003, 0.01, 0.03, 0.1, 0.3)}


def noise_candidates() -> list[dict]:
    """The reference and every perturbation of it, one variable at a time, up and down."""
    out = [{}]
    for var, steps in NOISE_H.items():
        base = getattr(REF, var)
        out += [{var: base + sign * h} for h in steps for sign in (-1, +1)]
    return out


def stage_noise(args) -> None:
    with pool_session("rotor sweeps: noise floor") as pool:
        run_stage(open_evaluator(pool), "noise", noise_candidates())


# ---------------------------------------------------------------------------------------------- step 4b: coarse against fine
def pick_for_fine(rows: list[dict]) -> list[dict]:
    """From one sweep's coarse rows (dicts of floats, sorted by x): the two ends, the best torque, the lowest ripple and the middle one.
    The ends are what the optimizer's bounds depend on, the two best points are where it will look, and the middle one says whether the
    mesh changes the shape between them. Returns the rows, each once."""
    order = sorted(rows, key=lambda r: r["x"])
    picks = [order[0], order[-1], max(order, key=lambda r: r["T_mean"]), min(order, key=lambda r: r["ripple_pct"]), order[len(order) // 2]]
    seen, out = set(), []
    for r in picks:
        if r["x"] not in seen:
            seen.add(r["x"])
            out.append(r)
    return out


def run_fine(pool) -> None:
    """Re-score points of every sweep on the fine mesh (its own Evaluator and CSV, evaluations_fine.csv): the same designs, found by their overrides in
    the coarse CSV, and the reference to check the fine mesh against the numbers the brief gives."""
    data = load_sweeps()
    ev = open_evaluator(pool, fine=True)
    run_stage(ev, "fine_reference", [{}])
    for name, rows in data.items():
        run_stage(ev, f"fine_{name}", [json.loads(r["overrides"]) for r in pick_for_fine(rows)])


def stage_fine(args) -> None:
    with pool_session("rotor sweeps: fine mesh") as pool:
        run_fine(pool)


# ---------------------------------------------------------------------------------------------- step 5: the 2-D map
# The pair is mag_y and v_half_deg, on the guard (bridge_mm held at its minimum, so mag_x follows). Why this pair, chosen from the one-at-a-time sweeps:
# the radius moves the ripple in one direction only (a thinner bridge, a larger ripple, a larger torque) and the torque keeps rising to the guard, so
# the optimizer will sit on the guard in the radius; mag_y and v_half_deg are the two variables whose ripple curves have several valleys, and on the
# guard mag_y has the widest ripple swing of all. The window is where the torque is worth having: mag_y from 4.9 (the centre rib's limit at 60 degrees)
# to 7.5 (the torque has fallen to 0.7 N.m there), v_half_deg over the whole scanned range, which holds both valleys the sweep showed (about 36 and 75).
# The corner of low mag_y and high v_half_deg is infeasible (the centre rib shrinks with both): screen() says so for nothing, and the map shows the edge.
MAP = {"x": "v_half_deg", "y": "mag_y", "xlo": 30.0, "xhi": 89.0, "nx": 13, "ylo": 4.9, "yhi": 7.5, "ny": 9, "hold": {"bridge_mm": GUARD_MM}}


def map_candidates(m: dict) -> tuple[dict, list[dict]]:
    """The map's spec (what map.json holds: the axes, the values and what is held) and its candidates, row by row of y."""
    xs = [float(v) for v in np.linspace(m["xlo"], m["xhi"], m["nx"])]
    ys = [float(v) for v in np.linspace(m["ylo"], m["yhi"], m["ny"])]
    spec = {"x": m["x"], "y": m["y"], "xs": xs, "ys": ys, "hold": m["hold"]}
    return spec, [{**m["hold"], m["x"]: x, m["y"]: y} for y in ys for x in xs]


def stage_map(args) -> None:
    """The 2-D map of MAP (or of the ranges given on the command line), under the stage label 'map'. map.json records which variables, which values and
    what was held, so that the analysis and the figure need nothing else."""
    m = dict(MAP)
    if args.x:
        m.update(x=args.x, y=args.y, xlo=args.xlo, xhi=args.xhi, ylo=args.ylo, yhi=args.yhi, nx=args.nx, ny=args.ny, hold=parse_hold(args.hold))
    spec, candidates = map_candidates(m)
    (folder() / "map.json").write_text(json.dumps(spec, indent=1), encoding="utf-8")
    say(f"map: {m['ny']} x {m['nx']} = {len(candidates)} points, {sum(design.screen(design.decode(c))['feasible'] for c in candidates)} of them feasible (the rest cost no solve)")
    with pool_session("rotor sweeps: 2-D map") as pool:
        run_stage(open_evaluator(pool), "map", candidates)


def stage_fine_map(args) -> None:
    """The map's own front (the grid points that no other grid point beats in both torque and ripple, by design.nondominated) on the fine mesh, under the
    label fine_map: does the best of the landscape, the deepest valley among them, survive a finer mesh? The ripple is where a mesh matters."""
    rows = [r for r in read_csv("evaluations.csv") if r["stage"] == "map" and r["feasible"]]
    keep = design.nondominated([design.objectives(r) for r in rows])
    front = [r for r, k in zip(rows, keep) if k]
    say(f"the map's front: {len(front)} of {len(rows)} grid points")
    with pool_session("rotor sweeps: the map's front on the fine mesh") as pool:
        run_stage(open_evaluator(pool, fine=True), "fine_map", [json.loads(r["overrides"]) for r in front])


# ---------------------------------------------------------------------------------------------- the map's deepest valley, between its grid points
# The map's lowest ripple sits on its lowest row (mag_y 4.9, the centre rib's limit) at v_half_deg 54.6, and the grid is 4.9 degrees apart there: the valley
# may be deeper between two columns. One line through it, 1.5 degrees apart, on the guard, up to 61.5 degrees: the centre rib reaches its 1.0 mm limit at
# about 61.75 degrees when mag_y is 4.9.
VALLEY = {"label": "valley", "var": "v_half_deg", "values": [48.0 + 1.5 * k for k in range(10)], "hold": {"bridge_mm": GUARD_MM, "mag_y": 4.9}}


def stage_valley(args) -> None:
    """The line through the map's deepest valley (VALLEY), under the stage label 'valley'."""
    with pool_session("rotor sweeps: the deepest valley") as pool:
        run_stage(open_evaluator(pool), VALLEY["label"], [{**VALLEY["hold"], VALLEY["var"]: v} for v in VALLEY["values"]])


# ---------------------------------------------------------------------------------------------- the strip the map left out
# The map starts at mag_y 4.9, the centre rib's limit at 60 degrees. The limit moves with the angle (screen() on the guard: 3.48 mm at 30 degrees, 4.28 at 45, 4.85
# at 60, 5.15 at 89), so below 60 degrees there is feasible room under the map's lowest row. Six angles, and at each the rib's own limit plus three offsets.
STRIP = {"label": "strip", "v_half_deg": [30.0, 33.0, 36.0, 39.0, 42.0, 45.0], "above_limit": [0.05, 0.4, 0.8]}


def min_mag_y(v_half_deg: float) -> float:
    """The smallest mag_y that is feasible on the guard at this V angle (the centre rib at its 1.0 mm limit): bisection on screen(), as the scan's edges are."""
    lo, hi = 1.0, 8.0                                         # infeasible, feasible
    for _ in range(40):
        mid = (lo + hi) / 2
        if probe({"bridge_mm": GUARD_MM, "mag_y": mid, "v_half_deg": v_half_deg})["feasible"]:
            hi = mid
        else:
            lo = mid
    return hi


def stage_strip(args) -> None:
    """The strip under the map's lowest row (STRIP), under the stage label 'strip'."""
    candidates = [{"bridge_mm": GUARD_MM, "v_half_deg": v, "mag_y": min_mag_y(v) + d} for v in STRIP["v_half_deg"] for d in STRIP["above_limit"]]
    with pool_session("rotor sweeps: the strip under the map") as pool:
        run_stage(open_evaluator(pool), STRIP["label"], candidates)


# ---------------------------------------------------------------------------------------------- the question: does the radius sweep depend on the other two?
# The radius sweep again, at other settings of mag_y and v_half_deg: a thin centre rib, a thick one, a narrow V, a flat V. The first bridge is the guard.
CROSS = {"bridge@mag_y=4.9,v_half_deg=60": {"mag_y": 4.9, "v_half_deg": 60.0}, "bridge@mag_y=7,v_half_deg=60": {"mag_y": 7.0, "v_half_deg": 60.0},
         "bridge@mag_y=6,v_half_deg=45": {"mag_y": 6.0, "v_half_deg": 45.0}, "bridge@mag_y=6,v_half_deg=85": {"mag_y": 6.0, "v_half_deg": 85.0}}
CROSS_BRIDGES = (GUARD_MM, 1.55, 1.7, 1.9, 2.2, 2.6, 3.0, 3.6, 4.5)


def stage_batch(args) -> None:
    """The radius sweeps at other settings (CROSS), the fine-mesh re-scores (step 4b) and the 2-D map (step 5), in ONE pool session: one wait for the
    solver slot, which another run may hold for minutes, instead of three. Needs the sweeps and scan.json."""
    spec, map_list = map_candidates(MAP)
    (folder() / "map.json").write_text(json.dumps(spec, indent=1), encoding="utf-8")
    with pool_session("rotor sweeps: radius at other settings, fine mesh, 2-D map") as pool:
        ev = open_evaluator(pool)
        for label, hold in CROSS.items():
            run_stage(ev, label, [{**hold, "bridge_mm": b} for b in CROSS_BRIDGES])
        say(f"map: {MAP['ny']} x {MAP['nx']} = {len(map_list)} points, {sum(design.screen(design.decode(c))['feasible'] for c in map_list)} of them feasible")
        run_stage(ev, "map", map_list)
        run_fine(pool)


# ---------------------------------------------------------------------------------------------- reading the CSVs back
NUMBERS = ("mag_x", "mag_y", "v_half_deg", "bridge_mm", "amps", "gamma", "T_mean", "ripple_pp", "ripple_pp_samples", "ripple_pct", "psi_pm",
           "psi_d", "psi_q", "L_d", "L_q", "saliency", "pf", "T_dq_mean", "nodes", "seconds")


def read_csv(name: str = "evaluations.csv") -> list[dict]:
    """Every row of a CSV of this run as a dict: the numbers as floats (nan where a cell is empty), feasible and cached as booleans."""
    path = folder() / name
    out = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            row = dict(r)
            for k in NUMBERS:
                row[k] = float(r[k]) if r.get(k, "") != "" else float("nan")
            row["feasible"] = r["feasible"] == "True"
            row["cached"] = r["cached"] == "True"
            row["centre_rib"] = design.LIMITS["centre_rib"] - float(r["G_centre_rib"])       # G is the limit less the value (<= 0 allowed)
            row["clipped"] = bool(row["gamma"] == row["gamma"] and (abs(row["gamma"] - design.GAMMAS[0]) < 0.005 or abs(row["gamma"] - design.GAMMAS[-1]) < 0.005))
            out.append(row)
    return out


def sweep_x(name: str, row: dict) -> float:
    """The x of a row in its sweep: the sweep's variable (for 'bridge' the bridge the candidate really has, which the CSV records)."""
    return row[SWEEPS[name]["var"]]


def load_sweeps(name_csv: str = "evaluations.csv", rows: list[dict] | None = None) -> dict[str, list[dict]]:
    """{sweep name: its scored rows sorted by x (each row has "x")}, from the CSV or from rows already read. A row the Evaluator took from its
    cache is kept (it is the same design, scored once already, and a point of this sweep too); a row that was not scored is left out."""
    rows = read_csv(name_csv) if rows is None else rows
    out = {}
    for name in SWEEP_ORDER:
        mine = [dict(r, x=sweep_x(name, r)) for r in rows if r["stage"] == name and r["feasible"]]
        out[name] = sorted(mine, key=lambda r: r["x"])
    return out


# ---------------------------------------------------------------------------------------------- steps 3 and 4: the analysis (no FEMM)
QUANTITIES = ("T_mean", "ripple_pct", "gamma", "psi_pm", "L_d", "L_q", "saliency", "pf")
SHOWN = {"psi_pm": 1e3, "L_d": 1e3, "L_q": 1e3}              # shown in mWb and mH, where the SI numbers are 0.0xx and 0.00x
TINY = 0.003                                                  # the brief's noise set: perturbations up to +-0.003 (mm, mm, degrees)
FD_ERROR = 0.10                                               # a finite-difference gradient is "usable" when its noise is under 10 % of it


def split_torque(rows: list[dict]) -> tuple[float, float]:
    """Add i_d, i_q, T_pm (the magnet's torque), T_rel (the rest) and share_pm to every row; return (k, the relative spread of k over the rows).
    The torque from the flux linkages is T = k (psi_d i_q - psi_q i_d), k = 3/2 pole pairs. Write psi_d = psi_pm + L_d i_d and psi_q = L_q i_q
    and it becomes T = k psi_pm i_q + k (L_d - L_q) i_d i_q: the magnet's torque and the reluctance torque. k is not assumed, it is read
    from the rows (if it were not one constant the split would mean nothing: the spread says so)."""
    ks = []
    for r in rows:
        r["i_d"], r["i_q"] = design.currents(r["amps"], r["gamma"])
        ks.append(r["T_dq_mean"] / (r["psi_d"] * r["i_q"] - r["psi_q"] * r["i_d"]))
    k = float(np.median(ks))
    for r in rows:
        r["T_pm"] = k * r["psi_pm"] * r["i_q"]
        r["T_rel"] = r["T_dq_mean"] - r["T_pm"]
        r["share_pm"] = r["T_pm"] / r["T_dq_mean"]
    return k, float((max(ks) - min(ks)) / k)


def cubic_roughness(x, y, noise: float) -> float:
    """How far the curve is from smooth at the spacing of its points, in units of the noise. Each interior point is predicted from the cubic through its
    two neighbours on each side (Lagrange weights w, so the points need not be evenly spaced); the miss of a point that lies on a smooth curve is
    tiny (the fourth derivative times the spacing to the fourth), and the miss of a point with noise on it is noise x sqrt(1 + sum w^2). The ratio of
    the rms miss to the rms of that expectation is about 1 for a curve smooth down to the noise, and large for a curve that has structure (a valley,
    a kink, a jump) at the scale of its points. nan if there is no noise floor to divide by."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    miss, expected = [], []
    for i in range(2, len(x) - 2):
        idx = [i - 2, i - 1, i + 1, i + 2]
        w = np.array([np.prod([x[i] - x[b] for b in idx if b != a]) / np.prod([x[a] - x[b] for b in idx if b != a]) for a in idx])
        miss.append(y[i] - float(w @ y[idx]))
        expected.append(noise * np.sqrt(1.0 + float((w ** 2).sum())))
    if not miss or not noise > 0:
        return float("nan")
    return float(np.sqrt(np.mean(np.square(miss))) / np.sqrt(np.mean(np.square(expected))))


def shape(x, y, kind: str, noise: float) -> dict:
    """What the curve y(x) (x increasing) looks like. kind is "max" for a torque (larger is better) or "min" for a ripple.
    noise is the scatter that scoring the same point again with a hair's difference gives (step 4a); a bump smaller than 3 noises is not a feature.
      best         the best point: its x, its value, and whether it is an end of the interval (a bound or a constraint) or inside
      gain_inside  how much better the best point is than the better of the two ends (0 when an end is the best)
      optima       interior local extrema of the kind sought whose prominence (how far the curve falls from it, on both sides, before it reaches a point as good) is at least
                   3 noises: the number of real valleys of a ripple, or real peaks of a torque; "other_extrema" are those of the opposite kind
      monotonic    no interior extremum of either kind at that prominence
      trend        the sign of last minus first
      roughness    cubic_roughness: about 1 when the curve is smooth to the noise at the spacing of the points, large when it has structure there"""
    from scipy.signal import find_peaks
    x, y = np.asarray(x, float), np.asarray(y, float)
    z = y if kind == "max" else -y                                  # larger z is better in both cases
    best = int(np.argmax(z))
    gain = float(z[best] - max(z[0], z[-1])) if 0 < best < len(z) - 1 else 0.0
    prom = max(3.0 * noise, 1e-12)
    good, _ = find_peaks(z, prominence=prom)
    bad, _ = find_peaks(-z, prominence=prom)
    return {"best_x": float(x[best]), "best_value": float(y[best]), "best_is_end": bool(best in (0, len(x) - 1)),
            "best_end": "lower end" if best == 0 else "upper end" if best == len(x) - 1 else "inside",
            "distance_from_nearest_end_pct": float(100.0 * min(x[best] - x[0], x[-1] - x[best]) / (x[-1] - x[0])),
            "gain_inside": gain, "optima": [[float(x[i]), float(y[i])] for i in good], "other_extrema": [[float(x[i]), float(y[i])] for i in bad],
            "monotonic": bool(len(good) == 0 and len(bad) == 0), "trend": "rising" if y[-1] > y[0] else "falling", "first": float(y[0]), "last": float(y[-1]),
            "span": float(y.max() - y.min()), "roughness": cubic_roughness(x, y, noise)}


def perturbed(row: dict) -> tuple[str, float] | None:
    """(variable, how far from the reference) of a row that differs from the reference in one rotor variable, or None for the reference itself."""
    off = [(v, row[v] - getattr(REF, v)) for v in NOISE_H if abs(row[v] - getattr(REF, v)) > 1e-9]
    return off[0] if len(off) == 1 else None


def noise_floor(rows: list[dict]) -> dict:
    """The scatter of the score under a hair's change of the geometry, and what it means for gradients.
    Per variable: the reference and its perturbations (stage "noise"). Over the brief's tiny set (|h| <= 0.003) the true response is a straight
    line to better than the noise (its curvature over 0.003 mm is nothing), so what a line fitted through those 7 points misses is noise:
    sigma = sqrt(sum of squared misses / (n - 2)); the same without the reference (6 points) says whether the exact reference is an outlier, and
    reference_minus_its_neighbours is how far it sits from the line through the six. The slope b comes from the fit through the whole ladder (|h| up to
    0.05 mm or 0.3 deg).
    A central difference of step h has the error sigma / (sqrt(2) h) on a gradient b, so it is usable (error under 10 %) from
    h_10 = sigma / (0.1 sqrt(2) |b|) = 7.07 sigma / |b|."""
    ref = next(r for r in rows if r["stage"] == "reference" and r["feasible"])
    out = {"reference": {q: ref[q] for q in QUANTITIES}, "variables": {}}
    pooled = {q: [] for q in QUANTITIES}
    for var in NOISE_H:
        mine = [(0.0, ref)] + [(h, r) for r in rows if r["stage"] == "noise" and r["feasible"] and not r["cached"]
                               for v_h in [perturbed(r)] if v_h and v_h[0] == var for h in [v_h[1]]]
        mine.sort(key=lambda t: t[0])
        h_all = np.array([h for h, _ in mine])
        tiny = np.abs(h_all) <= TINY + 1e-12
        entry = {"points": len(mine), "tiny_points": int(tiny.sum()), "quantities": {}}
        for q in QUANTITIES:
            y = np.array([r[q] for _, r in mine])
            a1, a0 = np.polyfit(h_all[tiny], y[tiny], 1)
            miss = y[tiny] - (a0 + a1 * h_all[tiny])
            sigma = float(np.sqrt(np.sum(miss ** 2) / (tiny.sum() - 2)))
            slope = float(np.polyfit(h_all, y, 1)[0])
            near = tiny & (h_all != 0)                                  # the same, without the reference: is the exact reference a mesh of its own?
            b1, b0 = np.polyfit(h_all[near], y[near], 1)
            sigma_near = float(np.sqrt(np.sum((y[near] - (b0 + b1 * h_all[near])) ** 2) / (near.sum() - 2)))
            ref_offset = float(y[h_all == 0][0] - b0)
            pooled[q].append(sigma ** 2)
            # signal to noise of a central difference at each step h: |y(+h) - y(-h)| over sqrt(2) sigma
            snr = {}
            for h in sorted({abs(t) for t in h_all if abs(t) > 0}):
                up = [r[q] for t, r in mine if abs(t - h) < 1e-12]
                dn = [r[q] for t, r in mine if abs(t + h) < 1e-12]
                if up and dn:
                    snr[f"{h:g}"] = float(abs(up[0] - dn[0]) / (np.sqrt(2.0) * sigma)) if sigma > 0 else float("inf")
            entry["quantities"][q] = {"sigma": sigma, "sigma_without_reference": sigma_near, "reference_minus_its_neighbours": ref_offset,
                                      "tiny_range": float(y[tiny].max() - y[tiny].min()), "slope": slope,
                                      "h_10": float(sigma / (FD_ERROR * np.sqrt(2.0) * abs(slope))) if slope else float("inf"), "snr_central_difference": snr}
        out["variables"][var] = entry
    out["pooled_sigma"] = {q: float(np.sqrt(np.mean(v))) for q, v in pooled.items()}
    return out


def discordant_pairs(a, b, resolved: float) -> tuple[int, int]:
    """(pairs the two scorings order differently, pairs compared): of every two points whose values in `a` differ by more than `resolved` (the noise the
    scoring cannot see through), how many does `b` put the other way round. 0 of m means the second mesh keeps every order the first one can resolve."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    pairs = [(i, j) for i in range(len(a)) for j in range(i + 1, len(a)) if abs(a[i] - a[j]) > resolved]
    return sum(1 for i, j in pairs if (a[i] - a[j]) * (b[i] - b[j]) < 0), len(pairs)


def coarse_vs_fine(coarse: list[dict], fine: list[dict], sig_T: float = 0.0, sig_R: float = 0.0) -> dict:
    """The fine-mesh re-scores held against the coarse scores of the same designs (found by their overrides). Per stage: the differences, the rank
    correlation of the coarse and the fine values over the re-scored points, whether the best points agree, and the discordant pairs among the pairs the
    coarse mesh resolves (values more than 3 noise floors apart: a difference smaller than that is not an order the coarse mesh can claim). "overall" pools
    every stage's points whose torque is at least 0.5 N.m (below that the ripple in per cent is a ratio of two small numbers)."""
    from scipy.stats import spearmanr
    index = {r["overrides"]: r for r in coarse if r["feasible"]}
    by_stage: dict[str, list] = {}
    for r in fine:
        if r["feasible"] and r["overrides"] in index:
            by_stage.setdefault(r["stage"], []).append({"coarse": index[r["overrides"]], "fine": r})
    out = {}
    for stage, pairs in by_stage.items():
        T_c, T_f = np.array([p["coarse"]["T_mean"] for p in pairs]), np.array([p["fine"]["T_mean"] for p in pairs])
        R_c, R_f = np.array([p["coarse"]["ripple_pct"] for p in pairs]), np.array([p["fine"]["ripple_pct"] for p in pairs])
        t_bad, t_all = discordant_pairs(T_c, T_f, 3 * sig_T)
        r_bad, r_all = discordant_pairs(R_c, R_f, 3 * sig_R)
        out[stage] = {
            "n": len(pairs), "T_diff_pct": [float(v) for v in 100 * (T_f / T_c - 1)], "ripple_diff_points": [float(v) for v in R_f - R_c],
            "gamma_diff_deg": [float(p["fine"]["gamma"] - p["coarse"]["gamma"]) for p in pairs],
            "T_rank_correlation": float(spearmanr(T_c, T_f)[0]) if len(pairs) > 2 else float("nan"),
            "ripple_rank_correlation": float(spearmanr(R_c, R_f)[0]) if len(pairs) > 2 else float("nan"),
            "T_discordant_pairs": [t_bad, t_all], "ripple_discordant_pairs": [r_bad, r_all],
            "best_T_agrees": bool(int(np.argmax(T_c)) == int(np.argmax(T_f))), "lowest_ripple_agrees": bool(int(np.argmin(R_c)) == int(np.argmin(R_f))),
            "rows": [{"overrides": p["fine"]["overrides"], "T_coarse": float(a), "T_fine": float(b), "ripple_coarse": float(c), "ripple_fine": float(d)}
                     for p, a, b, c, d in zip(pairs, T_c, T_f, R_c, R_f)]}
    pooled = {}
    for stage, pairs in by_stage.items():
        if stage == "fine_reference":
            continue
        for p in pairs:
            if p["coarse"]["T_mean"] >= 0.5:
                pooled[p["fine"]["overrides"]] = p
    if pooled:
        dT = np.array([100 * (p["fine"]["T_mean"] / p["coarse"]["T_mean"] - 1) for p in pooled.values()])
        dR = np.array([p["fine"]["ripple_pct"] - p["coarse"]["ripple_pct"] for p in pooled.values()])
        dg = np.array([p["fine"]["gamma"] - p["coarse"]["gamma"] for p in pooled.values()])
        t_bad, t_all = discordant_pairs([p["coarse"]["T_mean"] for p in pooled.values()], [p["fine"]["T_mean"] for p in pooled.values()], 3 * sig_T)
        r_bad, r_all = discordant_pairs([p["coarse"]["ripple_pct"] for p in pooled.values()], [p["fine"]["ripple_pct"] for p in pooled.values()], 3 * sig_R)
        out["overall"] = {"n": len(pooled), "T_diff_pct": {"mean": float(dT.mean()), "min": float(dT.min()), "max": float(dT.max())},
                          "ripple_diff_points": {"mean": float(dR.mean()), "min": float(dR.min()), "max": float(dR.max())},
                          "gamma_diff_deg": {"mean": float(dg.mean()), "min": float(dg.min()), "max": float(dg.max())},
                          "T_discordant_pairs": [t_bad, t_all], "ripple_discordant_pairs": [r_bad, r_all]}
    return out


def local_minima_2d(z: np.ndarray, prominence: float) -> list[tuple[int, int, float]]:
    """The grid points that are lower than all 8 neighbours (a grid edge has fewer neighbours) by at least `prominence` over the lowest of them:
    (row, column, value), the lowest first. A plateau of equal values would hide a minimum, which is why 'lower than ALL neighbours' is strict."""
    out = []
    ny, nx = z.shape
    for i in range(ny):
        for j in range(nx):
            if not np.isfinite(z[i, j]):
                continue
            around = [z[a, b] for a in range(max(0, i - 1), min(ny, i + 2)) for b in range(max(0, j - 1), min(nx, j + 2)) if (a, b) != (i, j)]
            around = [v for v in around if np.isfinite(v)]
            if around and z[i, j] < min(around) - prominence:
                out.append((i, j, float(z[i, j])))
    return sorted(out, key=lambda t: t[2])


def map_grids(rows: list[dict], spec: dict) -> dict:
    """The map's rows as 2-D arrays [y index, x index]: T_mean, ripple_pct, gamma, psi_pm, saliency (nan where a point was not scored)."""
    xs, ys = np.array(spec["xs"]), np.array(spec["ys"])
    grids = {q: np.full((len(ys), len(xs)), np.nan) for q in ("T_mean", "ripple_pct", "gamma", "psi_pm", "saliency", "bridge_mm", "mag_x")}
    for r in rows:
        if r["stage"] != "map" or not r["feasible"]:
            continue
        i, j = int(np.argmin(np.abs(ys - r[spec["y"]]))), int(np.argmin(np.abs(xs - r[spec["x"]])))
        for q in grids:
            grids[q][i, j] = r[q]
    return grids


def rep_points(rows: list[dict], anchor: float | None = None, n: int = 7) -> list[dict]:
    """About six to seven representative rows of a sweep for the report: both ends, the best torque, the lowest ripple, the point nearest the
    reference's value (the anchor), and the points nearest a third and two thirds of the way. Each once, in the order of x."""
    near = lambda x0: min(rows, key=lambda r: abs(r["x"] - x0))
    first, last = rows[0]["x"], rows[-1]["x"]
    picks = [rows[0], rows[-1], max(rows, key=lambda r: r["T_mean"]), min(rows, key=lambda r: r["ripple_pct"]),
             near(first + (last - first) / 3), near(first + 2 * (last - first) / 3)]
    if anchor is not None:
        picks.append(near(anchor))
    seen, out = set(), []
    for r in picks:
        if r["x"] not in seen:
            seen.add(r["x"])
            out.append(r)
    return sorted(out, key=lambda r: r["x"])[:n]


def row_line(r: dict, name: str) -> str:
    """One table line: x, T_mean, ripple, gamma, psi_pm [mWb], L_d, L_q [mH], saliency, pf, T_pm and T_rel [N.m] (the torque's split)."""
    return (f"{name:<16} {r['x']:>9.4f}  {r['T_mean']:>8.4f}  {r['ripple_pct']:>7.2f}  {r['gamma']:>6.2f}{'*' if r['clipped'] else ' '} {1e3 * r['psi_pm']:>8.3f}  "
            f"{1e3 * r['L_d']:>7.4f}  {1e3 * r['L_q']:>7.4f}  {r['saliency']:>6.3f}  {r['pf']:>6.3f}  {r['T_pm']:>7.4f}  {r['T_rel']:>7.4f}")


TABLE_HEAD = ("sweep            x            T_mean   ripple%  gamma   psi_pm     L_d      L_q     saliency  pf      T_pm    T_rel\n"
              "                 [mm|deg]     [N.m]    [%]      [deg]   [mWb]      [mH]     [mH]    L_q/L_d   -       [N.m]   [N.m]")


FIT_STAGES = ("reference", "mag_x", "mag_y", "v_half_deg")          # the sweeps from the reference: the fit's data; every other stage then tests it


FIT_STAGES = ("reference", "mag_x", "mag_y", "v_half_deg")          # the sweeps from the reference: the fit's data; every other stage then tests it


def leakage_fit(rows: list[dict], bridge_max: float = 3.5, rib_max: float = 4.7) -> dict:
    """The magnet flux linkage as the magnets' flux less what the saturated steel strips beside them shunt away: psi_pm = a - b (outer bridge) - c (centre rib).
    A strip of steel that the magnet drives into saturation carries a flux of about B_sat x its thickness x the stack whatever the magnet does, and that
    flux does not reach the air gap: so the loss is linear in the thickness, and b and c are the flux linkage lost per millimetre of each strip (b for the two
    outer bridges of a V together, c for the one centre rib). It is fitted by least squares on the sweeps from the reference (FIT_STAGES) with the bridge up to
    bridge_max and the rib up to rib_max mm (beyond them a strip no longer saturates fully and the line is not the law), and then tested on every other
    scored row inside the same limits: the guard sweeps, the other radius sweeps, the map. Returns the coefficients (mWb, mWb per mm), R^2 and the rms miss on
    the fit and on the test rows."""
    inside = [r for r in rows if r["feasible"] and r["bridge_mm"] <= bridge_max and r["centre_rib"] <= rib_max]
    use = [r for r in inside if r["stage"] in FIT_STAGES]
    test = [r for r in inside if r["stage"] not in FIT_STAGES]
    design_matrix = lambda rs: np.column_stack([np.ones(len(rs)), [r["bridge_mm"] for r in rs], [r["centre_rib"] for r in rs]])
    y = np.array([1e3 * r["psi_pm"] for r in use])
    coef, *_ = np.linalg.lstsq(design_matrix(use), y, rcond=None)
    miss = y - design_matrix(use) @ coef
    out = {"n": len(use), "bridge_max": bridge_max, "rib_max": rib_max, "a_mWb": float(coef[0]), "b_mWb_per_mm_outer_bridge": float(coef[1]),
           "c_mWb_per_mm_centre_rib": float(coef[2]), "r2": float(1 - (miss ** 2).sum() / ((y - y.mean()) ** 2).sum()), "rms_mWb": float(np.sqrt(np.mean(miss ** 2))),
           "largest_miss_mWb": float(np.abs(miss).max())}
    if test:
        y_t = np.array([1e3 * r["psi_pm"] for r in test])
        miss_t = y_t - design_matrix(test) @ coef
        out.update({"test_n": len(test), "test_rms_mWb": float(np.sqrt(np.mean(miss_t ** 2))), "test_largest_miss_mWb": float(np.abs(miss_t).max()),
                    "test_rms_over_mean_psi_pct": float(100 * np.sqrt(np.mean(miss_t ** 2)) / y_t.mean())})
    return out


def radius_split(rows: list[dict], at: list[float]) -> list[dict]:
    """The change of the torque from each bridge in `at` to the guard end of the sweep (its first row), split into the magnet's part and the rest, with the
    changes of everything the torque is made of: psi_pm, L_d, L_q, the saliency and the MTPA angle. Rows are the radius sweep sorted by bridge."""
    guard = rows[0]
    out = []
    for b in at:
        r = min(rows, key=lambda q: abs(q["x"] - b))
        out.append({"bridge_mm": r["x"], "dT": guard["T_mean"] - r["T_mean"], "dT_pm": guard["T_pm"] - r["T_pm"], "dT_rel": guard["T_rel"] - r["T_rel"],
                    "dpsi_pm_mWb": 1e3 * (guard["psi_pm"] - r["psi_pm"]), "dL_d_mH": 1e3 * (guard["L_d"] - r["L_d"]), "dL_q_mH": 1e3 * (guard["L_q"] - r["L_q"]),
                    "dsaliency": guard["saliency"] - r["saliency"], "dgamma_deg": guard["gamma"] - r["gamma"]})
    return out


def stage_analyse(args) -> None:
    """Steps 3 and 4 and the map's reading, from the CSVs alone (no FEMM): prints the tables and writes summary.json and tables.txt."""
    lines: list[str] = []

    def out(text: str = "") -> None:
        lines.append(text)
        say(text)

    scan = json.loads((folder() / "scan.json").read_text(encoding="utf-8"))
    coarse_all = read_csv("evaluations.csv")
    scored = [r for r in coarse_all if r["feasible"]]
    split_k, k_spread = split_torque(scored)
    summary: dict = {"k": split_k, "k_spread": k_spread, "scan": scan}
    fresh = [r for r in scored if not r["cached"]]
    out(f"coarse CSV: {len(coarse_all)} rows, {len(scored)} scored, {len(fresh)} solved afresh (the others came from the Evaluator's cache or were not feasible), "
        f"{sum(r['seconds'] for r in fresh):.0f} s of worker time in all")
    out(f"torque split: k = {split_k:.4f} (3/2 pole pairs = {1.5 * REF.pole_pairs:g}); spread of k over the {len(scored)} rows {100 * k_spread:.1e} %")
    # ---- step 1 again, for the record: the intervals and what ends them
    out("\nFEASIBLE INTERVALS (design.screen, Python only; others at the reference, 'guard_*' with bridge_mm on the guard)")
    for name in SCAN_ORDER:
        if name not in scan:
            continue
        sc = scan[name]
        head = f"  {name:<24} {sc['scan'][0]:g}..{sc['scan'][1]:g} step {sc['scan'][2]:g}: {sc['feasible_points']} of {sc['points']} points feasible"
        if sc["rejected_by_rounding"]:
            head += f"; {sc['rejected_by_rounding']} rejected by rounding at the guard"
        out(head)
        for iv in sc["intervals"]:
            out(f"      {iv['lo']:.5f} .. {iv['hi']:.5f}   ends: {iv['lo_limit']} | {iv['hi_limit']}{'   <- the reference is in it' if iv['contains_reference'] else ''}")
    # ---- the noise floor first: every threshold below uses it
    noise = noise_floor(coarse_all) if any(r["stage"] == "noise" for r in coarse_all) else None
    sig_T = noise["pooled_sigma"]["T_mean"] if noise else 0.0
    sig_R = noise["pooled_sigma"]["ripple_pct"] if noise else 0.0
    summary["noise"] = noise
    if noise:
        out("\nNOISE FLOOR (the reference and its perturbations; sigma from the +-0.003 set, 7 points a variable)")
        for q in ("T_mean", "ripple_pct", "gamma", "psi_pm", "L_d", "L_q"):
            each = ", ".join(f"{v}: {e['quantities'][q]['sigma'] * SHOWN.get(q, 1):.3g}" for v, e in noise["variables"].items())
            out(f"  {q:<11} pooled sigma {noise['pooled_sigma'][q] * SHOWN.get(q, 1):.6g}   ({each})")
        for var, e in noise["variables"].items():
            for q in ("T_mean", "ripple_pct"):
                d = e["quantities"][q]
                steps = ", ".join(f"{h}: {v:.1f}" for h, v in d["snr_central_difference"].items())
                unit, scale = ("mN.m", 1e3) if q == "T_mean" else ("points", 1.0)
                out(f"  {var:<10} {q:<10} the +-0.003 set spreads {d['tiny_range'] * scale:.3g} {unit} peak to peak (sigma {d['sigma'] * scale:.3g}, without the reference {d['sigma_without_reference'] * scale:.3g}; "
                    f"the reference sits {d['reference_minus_its_neighbours'] * scale:+.3g} {unit} from the line through its six neighbours); slope {d['slope']:+.4g} per unit; "
                    f"a central difference is usable (10 % error) from h = {d['h_10']:.3g}; its signal/noise at h = {steps}")
    # ---- the sweeps
    sweeps = load_sweeps(rows=coarse_all)
    summary["sweeps"] = {}
    out("\nSWEEPS (coarse mesh, 10 A, scoring AT; * = the MTPA angle is at the edge of its 10 to 40 degree window)")
    out(TABLE_HEAD)
    useful = 0.5 * next(r["T_mean"] for r in coarse_all if r["stage"] == "reference")     # ripple in per cent means little where the torque is nearly gone
    for name in SWEEP_ORDER:
        rows = sweeps[name]
        if len(rows) < 4:
            continue
        x = [r["x"] for r in rows]
        sT, sR = shape(x, [r["T_mean"] for r in rows], "max", sig_T), shape(x, [r["ripple_pct"] for r in rows], "min", sig_R)
        lo, hi = sweep_interval(scan, name)
        keep = [r for r in rows if r["T_mean"] >= useful]
        sRu = shape([r["x"] for r in keep], [r["ripple_pct"] for r in keep], "min", sig_R) if len(keep) > 4 else None
        summary["sweeps"][name] = {"var": SWEEPS[name]["var"], "n": len(rows), "interval_swept": [lo, hi], "T_mean": sT, "ripple_pct": sR, "ripple_pct_where_T_is_half_the_reference": sRu,
                                   "clipped_rows": sum(r["clipped"] for r in rows), "effect_over_noise_T": sT["span"] / sig_T if sig_T else None,
                                   "effect_over_noise_ripple": sR["span"] / sig_R if sig_R else None,
                                   "rows": [{k: r[k] for k in ("x", "mag_x", "mag_y", "v_half_deg", "bridge_mm", "centre_rib", "T_mean", "ripple_pct", "ripple_pp", "gamma", "clipped", "psi_pm", "psi_d", "psi_q",
                                                                "L_d", "L_q", "saliency", "pf", "T_pm", "T_rel", "share_pm", "i_d", "i_q")} for r in rows]}
        for r in rep_points(rows, SWEEPS[name]["anchor"]):
            out(row_line(r, name))
        out(f"   -> {name}: {len(rows)} points over {x[0]:.4f}..{x[-1]:.4f}, {sum(r['clipped'] for r in rows)} with the MTPA angle at its window's edge. "
            f"T_mean {sT['trend']}, {sT['span']:.4f} N.m wide ({sT['span'] / sig_T:.0f} noise floors): best {sT['best_value']:.4f} at {sT['best_x']:.4f} ({sT['best_end']}; gain over the better end "
            f"{sT['gain_inside']:.4f}); turning points beyond 3 noise: {len(sT['optima'])} peak, {len(sT['other_extrema'])} dip; roughness {sT['roughness']:.1f}")
        out(f"      ripple {sR['span']:.2f} points wide ({sR['span'] / sig_R:.0f} noise floors): lowest {sR['best_value']:.2f} % at {sR['best_x']:.4f} ({sR['best_end']}, {sR['distance_from_nearest_end_pct']:.1f} % "
            f"of the interval from the nearest end); valleys beyond 3 noise: {len(sR['optima'])} at {[round(v[0], 3) for v in sR['optima']]}; hills: {len(sR['other_extrema'])} at "
            f"{[round(v[0], 3) for v in sR['other_extrema']]}; roughness {sR['roughness']:.1f}")
        if sRu:
            out(f"      where T_mean is at least {useful:.3f} N.m: {len(keep)} points, ripple {min(r['ripple_pct'] for r in keep):.2f}..{max(r['ripple_pct'] for r in keep):.2f} %, "
                f"valleys {len(sRu['optima'])} at {[round(v[0], 3) for v in sRu['optima']]}, lowest at {sRu['best_x']:.4f} ({sRu['best_end']})")
    # ---- the gradients where the optimizer will work: on the guard, at the reference's value of the variable, from the two points on each side of it
    summary["guard_slopes"] = {}
    if noise:
        out("\nGRADIENTS ON THE GUARD (central slope over the neighbours of the reference's value; h_10 is the smallest central-difference step whose noise is under 10 % of it)")
        for name in ("guard_mag_y", "guard_v_half_deg"):
            rows = sweeps[name]
            i = min(range(1, len(rows) - 1), key=lambda k: abs(rows[k]["x"] - SWEEPS[name]["anchor"]))
            dx = rows[i + 1]["x"] - rows[i - 1]["x"]
            entry = {}
            for q, sig in (("T_mean", sig_T), ("ripple_pct", sig_R)):
                slope = (rows[i + 1][q] - rows[i - 1][q]) / dx
                entry[q] = {"slope": slope, "h_10": float(sig / (FD_ERROR * np.sqrt(2.0) * abs(slope)))}
            summary["guard_slopes"][name] = {"x": rows[i]["x"], "over": [rows[i - 1]["x"], rows[i + 1]["x"]], **entry}
            out(f"  {name:<17} at {rows[i]['x']:.3g} (over {rows[i - 1]['x']:.3f}..{rows[i + 1]['x']:.3f}): T_mean {entry['T_mean']['slope']:+.4f} per unit, h_10 {entry['T_mean']['h_10']:.3g}; "
                f"ripple {entry['ripple_pct']['slope']:+.3f} points per unit, h_10 {entry['ripple_pct']['h_10']:.3g}")
    # ---- the radius sweep: why the torque does what it does
    if sweeps["bridge"]:
        rows = sweeps["bridge"]
        fit = leakage_fit(coarse_all)
        summary["leakage_fit"] = fit
        out(f"\nTHE RADIUS SWEEP: psi_pm = {fit['a_mWb']:.2f} {fit['b_mWb_per_mm_outer_bridge']:+.3f} x (outer bridge) {fit['c_mWb_per_mm_centre_rib']:+.3f} x (centre rib), in mWb and mm, fitted on "
            f"{fit['n']} rows of the sweeps from the reference with the bridge up to {fit['bridge_max']:g} mm and the rib up to {fit['rib_max']:g} mm: R^2 {fit['r2']:.4f}, rms miss {fit['rms_mWb']:.3f} mWb, "
            f"largest miss {fit['largest_miss_mWb']:.2f} mWb")
        if "test_n" in fit:
            out(f"  tested on {fit['test_n']} other rows inside the same limits (guard sweeps, other radius sweeps, the map): rms miss {fit['test_rms_mWb']:.3f} mWb ({fit['test_rms_over_mean_psi_pct']:.1f} % of their mean psi_pm), "
                f"largest {fit['test_largest_miss_mWb']:.2f} mWb")
        split = radius_split(rows, [1.5585, 2.0, 2.5, 3.0, 3.5])
        summary["radius_split"] = split
        out("  from the bridge shown to the guard (1.4 mm), what the torque gains, and from where:")
        for d in split:
            out(f"    from {d['bridge_mm']:.3f} mm: T +{d['dT']:.4f} N.m = magnet {d['dT_pm']:+.4f} + reluctance {d['dT_rel']:+.4f}; psi_pm {d['dpsi_pm_mWb']:+.3f} mWb, "
                f"L_d {d['dL_d_mH']:+.4f} mH, L_q {d['dL_q_mH']:+.4f} mH, saliency {d['dsaliency']:+.4f}, MTPA angle {d['dgamma_deg']:+.2f} deg")
        near = rows[:3]
        slope = (near[2]["T_mean"] - near[0]["T_mean"]) / (near[2]["x"] - near[0]["x"])
        out(f"  the torque's slope at the guard: dT/d(bridge) = {slope:+.3f} N.m per mm over {near[0]['x']:.3f}..{near[2]['x']:.3f} mm; the noise floor of 3 sigma ({3 * sig_T * 1e3:.2f} mN.m) is "
            f"{3 * sig_T / abs(slope) * 1e3:.1f} micrometres of bridge")
        summary["slope_at_guard_N_m_per_mm"] = slope
    # ---- the radius sweeps at other values of the other two variables
    cross = sorted({r["stage"] for r in coarse_all if r["stage"].startswith("bridge@")})
    summary["cross"] = {}
    for stage in cross:
        rows = sorted([dict(r, x=r["bridge_mm"]) for r in coarse_all if r["stage"] == stage and r["feasible"]], key=lambda r: r["x"])
        if len(rows) < 3:
            continue
        sT = shape([r["x"] for r in rows], [r["T_mean"] for r in rows], "max", sig_T)
        sR = shape([r["x"] for r in rows], [r["ripple_pct"] for r in rows], "min", sig_R)
        fit_rows = [r for r in rows if r["T_mean"] > 0]
        slope = (rows[2]["T_mean"] - rows[0]["T_mean"]) / (rows[2]["x"] - rows[0]["x"])
        summary["cross"][stage] = {"n": len(rows), "T_mean": sT, "ripple_pct": sR, "slope_at_guard": slope,
                                   "rows": [{k: r[k] for k in ("x", "mag_x", "mag_y", "v_half_deg", "centre_rib", "T_mean", "ripple_pct", "gamma", "psi_pm", "L_d", "L_q", "saliency", "pf", "T_pm", "T_rel")} for r in fit_rows]}
        out(f"\n{stage}: {len(rows)} points, bridge {rows[0]['x']:.3f}..{rows[-1]['x']:.3f}: T_mean {sT['trend']}, best {sT['best_value']:.4f} at bridge {sT['best_x']:.3f} ({sT['best_end']}), "
            f"turning points {len(sT['optima'])} peak {len(sT['other_extrema'])} dip; slope at the guard {slope:+.3f} N.m per mm; ripple lowest {sR['best_value']:.2f} % at {sR['best_x']:.3f} ({sR['best_end']}), "
            f"valleys {len(sR['optima'])}")
        for r in rep_points(rows):
            out(row_line(r, stage[7:23]))
    # ---- the line through the deepest valley
    line = sorted([r for r in coarse_all if r["stage"] == VALLEY["label"] and r["feasible"]], key=lambda r: r["v_half_deg"])
    if line:
        low = min(line, key=lambda r: r["ripple_pct"])
        summary["valley"] = {"hold": VALLEY["hold"], "rows": [{k: r[k] for k in ("v_half_deg", "mag_y", "mag_x", "centre_rib", "T_mean", "ripple_pct", "gamma")} for r in line]}
        out(f"\nTHE DEEPEST VALLEY, a line through it (mag_y 4.9 on the guard, v_half_deg {line[0]['v_half_deg']:.1f} to {line[-1]['v_half_deg']:.1f} degrees, {len(line)} points scored): "
            f"lowest {low['ripple_pct']:.2f} % at {low['v_half_deg']:.1f} degrees, T {low['T_mean']:.4f} N.m, centre rib {low['centre_rib']:.3f} mm")
        out("  v_half_deg  rib[mm]  T_mean   ripple%   gamma   " + "   ".join(f"{r['v_half_deg']:.1f}: {r['ripple_pct']:.2f}" for r in line))
    # ---- the strip under the map
    strip = sorted([r for r in coarse_all if r["stage"] == STRIP["label"] and r["feasible"]], key=lambda r: (r["v_half_deg"], r["mag_y"]))
    if strip:
        summary["strip"] = [{k: r[k] for k in ("v_half_deg", "mag_y", "mag_x", "centre_rib", "T_mean", "ripple_pct", "gamma", "psi_pm")} for r in strip]
        out(f"\nTHE STRIP UNDER THE MAP (on the guard; mag_y from the centre rib's limit up; {len(strip)} points scored): T_mean {min(r['T_mean'] for r in strip):.4f}..{max(r['T_mean'] for r in strip):.4f}, "
            f"ripple {min(r['ripple_pct'] for r in strip):.2f}..{max(r['ripple_pct'] for r in strip):.2f} %")
        out("  v_half_deg  mag_y   rib[mm]  T_mean   ripple%  gamma")
        for r in strip:
            out(f"  {r['v_half_deg']:8.1f}  {r['mag_y']:6.3f}  {r['centre_rib']:6.3f}  {r['T_mean']:7.4f}  {r['ripple_pct']:7.2f}  {r['gamma']:6.2f}")
    # ---- coarse against fine
    fine_path = folder() / "evaluations_fine.csv"
    if fine_path.exists():
        fine_rows = read_csv("evaluations_fine.csv")
        split_torque([r for r in fine_rows if r["feasible"]])
        cf = coarse_vs_fine(coarse_all, fine_rows, sig_T, sig_R)
        summary["coarse_vs_fine"] = cf
        out("\nCOARSE AGAINST FINE (the same designs, in the order of the sweep's variable; the fine mesh finds its own MTPA angle)")
        for stage, c in cf.items():
            if stage == "overall":
                continue
            out(f"  {stage:<22} n {c['n']}: T diff % {', '.join(f'{v:+.2f}' for v in c['T_diff_pct'])}; ripple diff points {', '.join(f'{v:+.2f}' for v in c['ripple_diff_points'])}; "
                f"rank corr T {c['T_rank_correlation']:.2f}, ripple {c['ripple_rank_correlation']:.2f}; discordant pairs the coarse mesh resolves: T {c['T_discordant_pairs'][0]} of "
                f"{c['T_discordant_pairs'][1]}, ripple {c['ripple_discordant_pairs'][0]} of {c['ripple_discordant_pairs'][1]}; best T agrees {c['best_T_agrees']}, lowest ripple agrees {c['lowest_ripple_agrees']}")
        if "overall" in cf:
            o = cf["overall"]
            out(f"  overall, {o['n']} designs with at least 0.5 N.m: fine minus coarse: T {o['T_diff_pct']['mean']:+.2f} % (range {o['T_diff_pct']['min']:+.2f}..{o['T_diff_pct']['max']:+.2f}), "
                f"ripple {o['ripple_diff_points']['mean']:+.2f} points (range {o['ripple_diff_points']['min']:+.2f}..{o['ripple_diff_points']['max']:+.2f}), "
                f"MTPA angle {o['gamma_diff_deg']['mean']:+.2f} deg (range {o['gamma_diff_deg']['min']:+.2f}..{o['gamma_diff_deg']['max']:+.2f}); "
                f"discordant pairs: T {o['T_discordant_pairs'][0]} of {o['T_discordant_pairs'][1]}, ripple {o['ripple_discordant_pairs'][0]} of {o['ripple_discordant_pairs'][1]}")
        fine_ref = next((r for r in fine_rows if r["stage"] == "fine_reference" and r["feasible"]), None)
        if fine_ref:
            out(f"  the reference on the fine mesh: gamma {fine_ref['gamma']:.2f}, T_mean {fine_ref['T_mean']:.5f}, ripple {fine_ref['ripple_pct']:.2f} % (the brief: 22.46, 1.12880, 18.45 %)")
    # ---- the 2-D map
    map_path = folder() / "map.json"
    if map_path.exists() and any(r["stage"] == "map" for r in coarse_all):
        spec = json.loads(map_path.read_text(encoding="utf-8"))
        g = map_grids(coarse_all, spec)
        minima = local_minima_2d(g["ripple_pct"], 3.0 * sig_R)
        maxima = local_minima_2d(-g["T_mean"], 3.0 * sig_T)
        summary["map"] = {"spec": spec, "scored": int(np.isfinite(g["T_mean"]).sum()), "points": int(g["T_mean"].size),
                          "T_best": float(np.nanmax(g["T_mean"])), "ripple_best": float(np.nanmin(g["ripple_pct"])),
                          "ripple_minima": [{"x": spec["xs"][j], "y": spec["ys"][i], "ripple_pct": v, "T_mean": float(g["T_mean"][i, j])} for i, j, v in minima],
                          "T_maxima": [{"x": spec["xs"][j], "y": spec["ys"][i], "T_mean": -v, "ripple_pct": float(g["ripple_pct"][i, j])} for i, j, v in maxima],
                          "grids": {q: np.where(np.isfinite(v), v, None).tolist() for q, v in g.items()}}
        out(f"\n2-D MAP {spec['y']} (rows) against {spec['x']} (columns), held {spec['hold']}: {summary['map']['scored']} of {summary['map']['points']} points scored (the rest are infeasible); "
            f"T_mean {np.nanmin(g['T_mean']):.4f}..{np.nanmax(g['T_mean']):.4f}, ripple {np.nanmin(g['ripple_pct']):.2f}..{np.nanmax(g['ripple_pct']):.2f} %")
        out(f"  ripple valleys (a grid point lower than all its neighbours by 3 noise floors): {len(minima)}: " +
            "; ".join(f"({spec['x']} {spec['xs'][j]:.2f}, {spec['y']} {spec['ys'][i]:.3f}) {v:.2f} % at T {g['T_mean'][i, j]:.4f}" for i, j, v in minima))
        out(f"  torque peaks (a grid point higher than all its neighbours by 3 noise floors): {len(maxima)}: " +
            "; ".join(f"({spec['x']} {spec['xs'][j]:.2f}, {spec['y']} {spec['ys'][i]:.3f}) {-v:.4f} N.m at ripple {g['ripple_pct'][i, j]:.2f} %" for i, j, v in maxima))
        # the compromise that matters: the lowest ripple among the points that keep at least 90 % of the map's best torque, and the Pareto points
        T, R = g["T_mean"], g["ripple_pct"]
        ok = np.isfinite(T)
        pts = [(float(T[i, j]), float(R[i, j]), spec["xs"][j], spec["ys"][i]) for i in range(T.shape[0]) for j in range(T.shape[1]) if ok[i, j]]
        front = [p for p in pts if not any(q[0] >= p[0] and q[1] <= p[1] and (q[0] > p[0] or q[1] < p[1]) for q in pts)]
        summary["map"]["pareto"] = [{"T_mean": a, "ripple_pct": b, "x": c, "y": d} for a, b, c, d in sorted(front, reverse=True)]
        out(f"  the map's own front (no other grid point has more torque and less ripple): {len(front)} points: " +
            "; ".join(f"T {a:.3f} / {b:.2f} % at ({c:.1f}, {d:.2f})" for a, b, c, d in sorted(front, reverse=True)))
        for share in (0.9, 0.8):
            floor = share * float(np.nanmax(T))
            best = min((p for p in pts if p[0] >= floor), key=lambda p: p[1])
            out(f"  lowest ripple among the points with at least {int(share * 100)} % of the map's best torque ({floor:.3f} N.m): {best[1]:.2f} % at T {best[0]:.4f}, ({spec['x']} {best[2]:.2f}, {spec['y']} {best[3]:.3f})")
    (folder() / "summary.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf-8")
    (folder() / "tables.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------------------------- the command line
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m optimization.sweeps_rotor", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("scan", help="step 1: the feasible intervals (Python only)")
    p.add_argument("--only", nargs="*", help="scan names (default: all)")
    p = sub.add_parser("sweeps", help="step 2: the reference and the six one-at-a-time sweeps (FEMM)")
    p.add_argument("--only", nargs="*", choices=SWEEP_ORDER, help="sweep names (default: all)")
    p.add_argument("--with-noise", action="store_true", help="also run the noise-floor stage (step 4a) in the same pool session")
    sub.add_parser("analyse", help="steps 3 and 4: tables, shapes, noise floor, coarse against fine, the map (no FEMM)")
    sub.add_parser("noise", help="step 4a: tiny perturbations of the reference (FEMM)")
    sub.add_parser("fine", help="step 4b: points of every sweep on the fine mesh (FEMM)")
    sub.add_parser("strip", help="the feasible strip under the map's lowest row (FEMM)")
    sub.add_parser("valley", help="a line through the map's deepest valley (FEMM)")
    sub.add_parser("fine-map", help="the 2-D map's front on the fine mesh (FEMM)")
    sub.add_parser("batch", help="the radius sweeps at other settings, the fine-mesh re-scores and the 2-D map in one pool session (FEMM)")
    p = sub.add_parser("map", help="step 5: the 2-D map of MAP, or of the ranges given here (FEMM)")
    p.add_argument("--x")
    p.add_argument("--y")
    for k in ("xlo", "xhi", "ylo", "yhi"):
        p.add_argument(f"--{k}", type=float)
    p.add_argument("--nx", type=int, default=11)
    p.add_argument("--ny", type=int, default=11)
    p.add_argument("--hold", nargs="*", help="name=value pairs held fixed; bridge_mm=guard puts the outer bridge on the guard")
    args = parser.parse_args(argv)
    {"scan": stage_scan, "sweeps": stage_sweeps, "noise": stage_noise, "fine": stage_fine,
     "map": stage_map, "batch": stage_batch, "fine-map": stage_fine_map, "valley": stage_valley, "strip": stage_strip, "analyse": stage_analyse}[args.stage](args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
