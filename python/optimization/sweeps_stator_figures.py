"""optimization.sweeps_stator_figures: what the stator sweeps found, as numbers and as pictures. It reads the CSVs and starts no FEMM.

What it is, and why
    optimization/sweeps_stator.py runs the solver and writes rows; this script only reads them. It is kept apart so that a finding can be
    recomputed and a figure redrawn without waiting for the solver slot, and so that the way a landscape is described is written once and
    applied to every sweep alike (two sweeps described with two ways of counting valleys could not be compared).

The measures, each with the reason it is that one (sigma is a noise floor, SIGMA_K = 3: a difference below 3 sigma is not called real)
    two noise floors   "Micrometre": the noise stage's 31 candidates, the reference and +-0.001, 0.002, 0.003 on each variable. Their sample
                       standard deviation (n - 1) is what the mesh does to a design moved by micrometres, and so what a finite difference
                       of tiny step meets; the per-variable value is also given after a straight line (a real slope) is taken out.
                       "Step scale" (step_noise): the points of a sweep are tens to hundreds of micrometres apart and between them the mesh
                       changes more (the ladder shows it; h_open's ripple jumps by 0.3 points). For an interior point i with neighbours at
                       x[i-1] and x[i+1] let w = (x[i+1] - x[i]) / (x[i+1] - x[i-1]) and r[i] = y[i] - (w y[i-1] + (1 - w) y[i+1]), how far
                       the point lies from the straight line between its neighbours. Noise alone of size sigma gives Var(r) = sigma^2 (1 + w^2 +
                       (1 - w)^2), so sigma = 1.4826 x median|r - median r| / sqrt(1 + w^2 + (1 - w)^2) (the median absolute deviation, which
                       the mesh's rare jumps do not move). A bending curve adds to r, a flat one has only the mesh in it, so the figure is the
                       median over the three sweeps with the smallest swing. Every 'real or not' test uses this one: it is the larger.
    valleys            the prominence of a local minimum is the lower of its two shoulders above it, a shoulder being the highest value met
                       walking away until a lower value (or the end of the sweep) is reached; an end point has one shoulder, so a sweep that
                       falls to its end has a minimum there. Minima with prominence above 3 sigma count (the raw count is printed too).
                       On the 2-D map the same is found by a search that extends the path with the lowest highest-point-so-far (8 neighbours).
    kind of a sweep    monotonic if no interior extremum has prominence above 3 sigma; flat if the whole swing is below 3 sigma; and the
                       best point's place (an end, an end within 3 sigma, or interior). "Smooth" is held against a cubic fitted over the
                       sweep: its residual rms, as a share of the swing (curvature a cubic cannot follow, and the mesh) and in sigmas.
    finite difference  a central difference D(h) = (f(x + h) - f(x - h)) / (2 h) of two scores with independent noise sigma has the noise
                       sigma / (sqrt(2) h). Held to 10 % of the slope s it needs h = 10 sigma / (sqrt(2) |s|), printed as a step and as a share
                       of the interval; the ladder prints D(h) itself for h = 0.001 to 0.1.
    coarse vs fine     Kendall's tau of the order of a variable's re-scored points on the two meshes (+1 the same order), the offsets
                       (fine - coarse) and where the lowest ripple lies on each.
    the 2-D map        additive share = 1 - R^2 of the least-squares fit a(x) + b(y) over the cells that exist: 0 means the two variables act
                       separately, and the larger it is, the more they interact.
    derived            T / sqrt(Pcu_rel) under J, the first-order torque at the reference's copper loss (exact where torque goes with
                       the current, and torque here goes with the current and its square); the AT and J picks side by side; the room the
                       constraints leave (from screen_scan.csv).
    figures            one PNG per variable (torque, ripple and, where AT ran, the copper loss against the variable; J blue and AT orange,
                       the two categorical slots of the palette that were checked on this surface for both colour-vision floors; the infeasible
                       zones shaded and named; the reference drawn as a line; fine-mesh scores as open diamonds; a ring on the J sweep's best
                       point; a bar in the left margin the length of 2 sigma (step scale) in that panel) and one for the map (one blue
                       ramp, darker = larger, refused cells hatched). The tables printed here are the table view of every figure.

Run (from the Python to FEMM folder; no solver slot is needed)
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures tables     about six values of every sweep, J and AT side by side
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures shapes     the landscape of every sweep, and the variables ranked by the ripple's swing
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures noise      the noise floors, the signal against them, the gradient steps
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures ladder     the central differences at h = 0.001 .. 0.1
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures fine       coarse against fine
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures map        the 2-D map's findings
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures derived    T per root copper loss, J against AT picks, the room the limits leave
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures cost       candidates, solves and worker time
    .venv\\Scripts\\python.exe -m optimization.sweeps_stator_figures figures    the six PNGs, in out\\optimization\\sweeps_stator\\

Numbers to expect (the run of 2026-11-08 that sweeps_stator.py describes)
    noise floors   micrometre 0.00056 N.m and 0.0194 points pooled; step scale 0.00022 N.m and 0.0592 points
    ripple swing   slot_open 16.15 points (273 sigma), slot_half_deg 9.86 (166), r_wall 9.23 (156), r_shoe 5.71 (96), h_open 0.61 (10): the
                   map takes the first two (the two largest swings, and the two that the open_vs_slot constraint ties together)
    landscapes     under J the torque of slot_half_deg, r_shoe and r_wall is monotonic (a cubic leaves 0.1-1.0 % of the swing) and the ripple of
                   slot_half_deg and r_shoe too; the ripple of r_wall falls to 18.97 % at 57.5 mm and stays within 0.18 points to the yoke
                   limit (2 raw minima, none above 3 sigma); slot_open's torque has one interior peak (0.31 mm) and its ripple a hump at
                   3.0 mm between two end basins (3 raw minima, 2 above 3 sigma); h_open's torque and ripple are 12 and 10 sigma wide, a
                   cubic leaves 20 % of the swing, 5 raw ripple minima and 1 above 3 sigma (a one-point spike at 0.076 mm)
    gradients      ripple, step for a 10 % gradient (step noise): slot_open 0.056 mm (1.6 % of the interval), slot_half_deg 0.29 deg (13 %),
                   r_wall 0.86 mm (8.8 %), r_shoe 2.25 mm (31 %), h_open 0.24 mm (30 %); torque 0.0024 deg, 0.029 mm, 0.011 mm, 0.0145 mm
                   (0.1 to 0.8 %), and 0.87 mm (110 %) for h_open, whose torque has no slope above the noise
    fine mesh      the re-scored point with the lowest ripple is the same on both meshes in every sweep except h_open (fine 0.6025, coarse
                   0.7342) and r_wall (fine 57.79, coarse 57.52: one flat valley bottom); h_open on the fine mesh: 19.64 % at 0.01 mm,
                   18.45 % at 0.3, 18.2-18.3 % from 0.54 to 0.8 (the coarse mesh gives 19.6 % at both 0.01 and 0.3)
    equal loss     T / sqrt(Pcu_rel) is highest at slot_half_deg 3.6969 (1.4943; 1.4912 at the limit 3.7747), r_wall 57.984 (1.2709, within
                   0.5 % from 57.79 to 58.18), r_shoe 46.301 (1.1473), slot_open 0.3114 (1.1846); reference 1.1217
    room           tooth_min >= 3 mm up to slot_half_deg 3.160 (3.5 mm: 2.850); yoke >= 4 mm up to r_wall 57.360 (3.5 mm: 57.840); the opening may be
                   at most 2 r_si sin(slot_half_deg - 0.3 deg): 2.002 mm at slot_half_deg 1.5521, 3.628 mm at the reference's 2.57

Self-check (the answers follow).
  1. A sweep's ripple has 5 raw local minima and 1 above 3 sigma. How many valleys has the motor? What would change the answer?
  2. The map's interaction share is 3.3 %, yet its ripple has a second valley (8.8 % at slot_open 1.46 mm, on the widest slot). How can
     both be true?
  3. The noise bar in a figure is 2 x 0.059 points long for the ripple, not 2 x 0.019. Why?
Answers.
  1. The sweep can resolve one: the others are dips of less than 3 x 0.059 = 0.18 points, which is how far neighbouring points of a flat
     sweep differ through the mesh alone. A quieter mesh (the fine one, with its own noise measured) could confirm or drop them; counting
     them would be counting the mesh.
  2. The share is a share of variance over the whole map, which the two main effects fill (the ripple spans 4 to 25 %); a valley 1.1 points
     deep on one edge of the map is a small part of that variance. A low share says the two variables can be tuned one at a time to
     first order, not that the map is free of local structure.
  3. 0.019 is what the mesh does to a design moved by micrometres. A sweep moves by tens to hundreds of micrometres between points, and the
     mesh changes more over that distance (h_open's ripple steps by 0.3 points between neighbours): a bar of 0.019 would make the
     sweeps look cleaner than the black box can deliver.
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

from optimization import design
from optimization import sweeps_stator as ss

SIGMA_K = 3.0                                  # a difference counts as real when it is more than this many noise floors (sigma)
REF = design.REFERENCE


# ---------------------------------------------------------------------------------------------- loading
def column(rows: list[dict], name: str) -> np.ndarray:
    return np.array([r[name] for r in rows], float)


def load_sweep(var: str, scoring: str) -> list[dict]:
    """The scored rows of one sweep, sorted by the variable, one row per value (the reference appears once)."""
    csv_name, stage = (ss.CSV_J, var) if scoring == "J" else (ss.CSV_AT, var + "_AT")
    rows = [r for r in ss.read_rows(csv_name, stage) if r["feasible"]]
    seen, out = set(), []
    for r in sorted(rows, key=lambda r: r[var]):
        key = round(r[var], 6)
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


def load_noise(stage: str = "noise") -> tuple[dict, dict]:
    """(the reference row, {variable: its perturbed rows}) of a noise stage."""
    rows = ss.read_rows(ss.CSV_J, stage)
    ref = next(r for r in rows if json.loads(r["overrides"]) == {})
    clusters = {v: [] for v in ss.VARS}
    for r in rows:
        ov = json.loads(r["overrides"])
        if len(ov) == 1:
            clusters[next(iter(ov))].append(r)
    return ref, clusters


# ---------------------------------------------------------------------------------------------- shape of a sweep
def prominence(y: np.ndarray, i: int) -> float:
    """How deep the basin of the local minimum y[i] is: the lower of its two shoulders above it. A shoulder is the highest value met walking
    away from i until a value below y[i] (or the end of the sweep) is reached; an end point has one shoulder only (the other is infinite)."""
    n = len(y)
    left = right = math.inf
    if i > 0:
        top, j = y[i], i - 1
        while j >= 0 and y[j] >= y[i]:
            top, j = max(top, y[j]), j - 1
        left = top - y[i]
    if i < n - 1:
        top, j = y[i], i + 1
        while j < n and y[j] >= y[i]:
            top, j = max(top, y[j]), j + 1
        right = top - y[i]
    return min(left, right)


def local_minima(y: np.ndarray, tol: float = 0.0) -> list[int]:
    """Indices of the local minima of a sequence (the two end points included when they are lower than their one neighbour) whose
    prominence is above tol. tol = 0 is the raw count; tol = 3 sigma is the count a noise floor of sigma cannot make up."""
    n = len(y)
    idx = [i for i in range(n) if (i == 0 or y[i] < y[i - 1]) and (i == n - 1 or y[i] < y[i + 1])]
    return [i for i in idx if prominence(y, i) > tol]


def describe(x: np.ndarray, y: np.ndarray, sigma: float, best: str) -> dict:
    """The landscape of one sweep in numbers. best = "max" (torque) or "min" (ripple). All the 'real' tests use SIGMA_K x sigma."""
    tol = SIGMA_K * sigma
    n = len(y)
    sign = 1.0 if best == "min" else -1.0               # work with the quantity to MINIMISE
    z = sign * y
    i_best = int(np.argmin(z))
    frac = (x[i_best] - x[0]) / (x[-1] - x[0])
    valleys = [i for i in local_minima(z, tol) if 0 < i < n - 1]      # interior basins of z (the good side)
    ridges = [i for i in local_minima(-z, tol) if 0 < i < n - 1]      # interior humps of z (the bad side)
    where = "interior"
    if abs(z[0] - z[i_best]) <= tol and i_best != 0:
        where = "low end (within noise)"
    if abs(z[-1] - z[i_best]) <= tol and i_best != n - 1:
        where = "high end (within noise)" if where == "interior" else where + " and high end (within noise)"
    if i_best == 0:
        where = "low end"
    if i_best == n - 1:
        where = "high end"
    swing = float(y.max() - y.min())
    if swing <= tol:
        kind = "flat (swing within noise)"
    elif not valleys and not ridges:
        kind = "monotonic " + ("up" if y[-1] > y[0] else "down")
    else:
        good, bad = ("peak", "dip") if best == "max" else ("valley", "hump")
        kind = f"{len(valleys)} interior {good}(s) and {len(ridges)} {bad}(s) above the noise"
    poly = np.polynomial.Polynomial.fit(x, y, 3)        # a cubic: what "smooth" is held against
    res = y - poly(x)
    rms = float(np.sqrt(np.sum(res ** 2) / (n - 4)))
    return {"i_best": i_best, "x_best": float(x[i_best]), "y_best": float(y[i_best]), "frac": float(frac), "where": where, "kind": kind,
            "swing": swing, "n_min_raw": len(local_minima(y)), "n_min_sig": len(local_minima(y, tol)),
            "min_raw_at": [float(x[i]) for i in local_minima(y)], "min_sig_at": [float(x[i]) for i in local_minima(y, tol)],
            "valleys_at": [float(x[i]) for i in valleys], "ridges_at": [float(x[i]) for i in ridges],
            "cubic_rms": rms, "cubic_rms_over_range": rms / swing if swing else float("nan"),
            "cubic_rms_over_sigma": rms / sigma if sigma else float("nan")}


def slope_at(x: np.ndarray, y: np.ndarray, x0: float) -> float:
    """The central difference at x0 from the two sweep points on either side of it (x0 itself is one of the points when it is in the sweep)."""
    j = int(np.argmin(np.abs(x - x0)))
    lo, hi = max(0, j - 1), min(len(x) - 1, j + 1)
    return float((y[hi] - y[lo]) / (x[hi] - x[lo]))


# ---------------------------------------------------------------------------------------------- the noise floor
def noise_summary(stage: str = "noise") -> dict:
    """The remeshing noise from the reference and its tiny perturbations: sigma (n - 1) of torque and ripple, pooled over all clusters and
    per variable, and the same after a straight line is taken out of each cluster (a trend would not be noise)."""
    ref, clusters = load_noise(stage)
    out = {"per_var": {}}
    pooled = {"T_mean": [ref["T_mean"]], "ripple_pct": [ref["ripple_pct"]]}
    for v, rows in clusters.items():
        step = np.array([r[v] - getattr(REF, v) for r in rows] + [0.0])
        d = {}
        for m in ("T_mean", "ripple_pct"):
            y = np.array([r[m] for r in rows] + [ref[m]])
            pooled[m] += [r[m] for r in rows]
            line = np.polyfit(step, y, 1)
            d[m] = {"std": float(np.std(y, ddof=1)), "ptp": float(np.ptp(y)), "mean": float(np.mean(y)), "slope": float(line[0]),
                    "detrended_std": float(np.std(y - np.polyval(line, step), ddof=2))}
        out["per_var"][v] = {"n": len(step), "steps": sorted(step.tolist()), **d}
    out["pooled"] = {m: {"n": len(y), "std": float(np.std(y, ddof=1)), "mean": float(np.mean(y)), "ptp": float(np.ptp(y))} for m, y in pooled.items()}
    out["ref"] = {m: ref[m] for m in ("T_mean", "ripple_pct", "gamma")}
    return out


# ---------------------------------------------------------------------------------------------- the tables
def has_stage(csv_name: str, stage: str) -> bool:
    return stage in ss.stages_done(csv_name)


def rep_indices(x: np.ndarray, T: np.ndarray, rip: np.ndarray, x_ref: float) -> list[int]:
    """About six representative points of a sweep: both ends, the reference, the lowest ripple, the highest torque, and (until there are six)
    the grid points nearest a third, two thirds, a half, a sixth and five sixths of the way along."""
    pick = [0, len(x) - 1, int(np.argmin(np.abs(x - x_ref))), int(np.argmin(rip)), int(np.argmax(T))]
    for f in (1 / 3, 2 / 3, 1 / 2, 1 / 6, 5 / 6):
        if len(set(pick)) < 6:
            pick.append(int(np.argmin(np.abs(x - (x[0] + f * (x[-1] - x[0]))))))
    return sorted(set(pick))


def print_rows(rows: list[dict], var: str, idx: list[int]) -> None:
    print(f"{'value':>9} {'amps':>7} {'T_mean':>8} {'rip %':>7} {'Pcu_rel':>7} {'J_rel':>6} {'gamma':>6} {'L_d mH':>7} {'L_q mH':>7} {'pf':>6}")
    for i in idx:
        r = rows[i]
        print(f"{r[var]:9.4f} {r['amps']:7.3f} {r['T_mean']:8.5f} {r['ripple_pct']:7.3f} {r['Pcu_rel']:7.3f} {r['J_rel']:6.3f} {r['gamma']:6.2f} "
              f"{1e3 * r['L_d']:7.3f} {1e3 * r['L_q']:7.3f} {r['pf']:6.3f}")


def print_tables() -> None:
    """The report's compact tables: about six representative points of every sweep, and the AT scores at the same values."""
    for v in ss.VARS:
        J = load_sweep(v, "J")
        x, T, rip = column(J, v), column(J, "T_mean"), column(J, "ripple_pct")
        idx = rep_indices(x, T, rip, getattr(REF, v))
        print(f"== {v} [{ss.UNIT[v]}], {len(J)} points, scoring J (the reference's current density)")
        print_rows(J, v, idx)
        if has_stage(ss.CSV_AT, v + "_AT"):
            AT = {round(r[v], 6): r for r in load_sweep(v, "AT")}
            print(f"   the same values, scoring AT (100 ampere-turns whatever the slot)")
            print_rows([AT[round(x[i], 6)] for i in idx], v, list(range(len(idx))))


# ---------------------------------------------------------------------------------------------- the findings
def sigmas() -> tuple[float, float]:
    """(torque, ripple %) noise for every 'is this difference real?' test below: the step-scale roughness of step_noise(), not the micrometre
    noise of the noise stage. The points of a sweep lie tens to hundreds of micrometres apart, and between them the mesh changes more
    than in a step of a few micrometres (the ladder shows it), so the micrometre figure would call mesh jumps real."""
    return step_noise()


def print_shapes() -> dict:
    """The landscape of every sweep, for the torque and the ripple, J and AT: where the best is, monotonic or not, how many valleys,
    how smooth, and how big the swing is against the noise. Returns the numbers for the figures."""
    s_T, s_R = sigmas()
    micro = noise_summary()["pooled"] if has_stage(ss.CSV_J, "noise") else None
    print(f"noise floor used (step scale, from the flattest sweeps): torque {s_T:.5f} N.m, ripple {s_R:.4f} %-points; a difference is "
          f"'real' above {SIGMA_K:g} sigma" + (f" (the micrometre noise stage gave {micro['T_mean']['std']:.5f} N.m and {micro['ripple_pct']['std']:.4f} %-points)" if micro else ""))
    found = {}
    for scoring in ("J", "AT"):
        for v in ss.VARS:
            if scoring == "AT" and not has_stage(ss.CSV_AT, v + "_AT"):
                continue
            rows = load_sweep(v, scoring)
            x = column(rows, v)
            lo, hi = x[0], x[-1]
            for name, col, best, sig, unit in (("T_mean", "T_mean", "max", s_T, "N.m"), ("ripple", "ripple_pct", "min", s_R, "%")):
                d = describe(x, column(rows, col), sig, best)
                found[(v, scoring, name)] = d
                snr = d["swing"] / sig if sig else float("nan")
                print(f"{v:>13} {scoring:>2} {name:>7}: best ({best}) {d['y_best']:.5f} {unit} at {d['x_best']:.4f} (fraction {d['frac']:.2f} of "
                      f"{lo:.4f}..{hi:.4f}; {d['where']}); {d['kind']}; swing {d['swing']:.5f} = {snr:.0f} sigma; "
                      f"cubic residual {100 * d['cubic_rms_over_range']:.1f} % of the swing = {d['cubic_rms_over_sigma']:.1f} sigma")
                if name == "ripple":
                    print(f"{'':>26} local minima: {d['n_min_raw']} raw at {[round(a, 3) for a in d['min_raw_at']]}; "
                          f"{d['n_min_sig']} above {SIGMA_K:g} sigma at {[round(a, 3) for a in d['min_sig_at']]}; "
                          f"interior valleys {[round(a, 3) for a in d['valleys_at']]}, humps {[round(a, 3) for a in d['ridges_at']]}")
    print("ranking of the variables by the ripple's swing under J (the 2-D map's pair is taken from here):")
    ranked = sorted(ss.VARS, key=lambda v: -found[(v, "J", "ripple")]["swing"])
    for v in ranked:
        d, t = found[(v, "J", "ripple")], found[(v, "J", "T_mean")]
        print(f"   {v:>13}: ripple swing {d['swing']:.3f} %-points ({d['swing'] / s_R if s_R else float('nan'):.0f} sigma); torque swing {t['swing']:.4f} N.m")
    return found


# ---------------------------------------------------------------------------------------------- noise, gradients, steps
def print_noise() -> dict:
    """The noise stage and what it means for finite-difference gradients (the formulas are in the module docstring)."""
    n = noise_summary()
    s_T, s_R = n["pooled"]["T_mean"]["std"], n["pooled"]["ripple_pct"]["std"]
    print(f"reference: T_mean {n['ref']['T_mean']:.5f} N.m, ripple {n['ref']['ripple_pct']:.4f} %, gamma {n['ref']['gamma']:.3f} deg")
    print(f"pooled over the reference and its {n['pooled']['T_mean']['n'] - 1} perturbed copies: T_mean mean {n['pooled']['T_mean']['mean']:.5f} sigma "
          f"{s_T:.5f} N.m (peak-to-peak {n['pooled']['T_mean']['ptp']:.5f}); ripple mean {n['pooled']['ripple_pct']['mean']:.4f} sigma {s_R:.4f} "
          f"%-points (peak-to-peak {n['pooled']['ripple_pct']['ptp']:.4f})")
    ref, clusters = load_noise()
    for v in ss.VARS:
        d = n["per_var"][v]
        print(f"   {v:>13} steps {d['steps']}: T sigma {d['T_mean']['std']:.5f} (after a straight line {d['T_mean']['detrended_std']:.5f}), "
              f"ripple sigma {d['ripple_pct']['std']:.4f} (after a line {d['ripple_pct']['detrended_std']:.4f})")
        rows = sorted(clusters[v] + [ref], key=lambda r: r[v])
        print("                 " + "  ".join(f"{r[v] - getattr(REF, v):+.3f}:{r['T_mean']:.5f}/{r['ripple_pct']:.3f}" for r in rows))
    t_step, r_step = step_noise()
    print(f"step-scale roughness between neighbouring sweep points (the flattest sweeps): torque {t_step:.5f} N.m, ripple {r_step:.4f} %-points")
    print("signal against noise, and finite differences. The slope is the central difference at the reference between the two neighbouring sweep "
          "points; a central difference with step h has noise sigma / (sqrt(2) h), so a step of h = 10 sigma / (sqrt(2) |slope|) keeps it to 10 % of the slope:")
    out = {}
    for v in ss.VARS:
        rows = load_sweep(v, "J")
        x = column(rows, v)
        length = x[-1] - x[0]
        for name, col, sig_micro, sig_step in (("T_mean", "T_mean", s_T, t_step), ("ripple", "ripple_pct", s_R, r_step)):
            y = column(rows, col)
            s = slope_at(x, y, getattr(REF, v))
            h10 = {k: (10 * sg / (math.sqrt(2) * abs(s)) if s else float("inf")) for k, sg in (("micro", sig_micro), ("step", sig_step))}
            e1 = {k: sg / (math.sqrt(2) * 1e-3 * abs(s)) if s else float("inf") for k, sg in (("micro", sig_micro), ("step", sig_step))}
            out[(v, name)] = {"slope": s, "h10": h10, "noise_over_slope_at_1um": e1, "swing_over_sigma_step": float(np.ptp(y)) / sig_step}
            print(f"   {v:>13} {name:>7}: swing {np.ptp(y):.5f} = {np.ptp(y) / sig_step:5.0f} step-sigma; slope at the reference {s:+.5f} per {ss.UNIT[v]}; "
                  f"step for 10 % gradient noise: {h10['micro']:.4f} {ss.UNIT[v]} ({100 * h10['micro'] / length:.1f} % of the interval) with the micrometre noise, "
                  f"{h10['step']:.4f} ({100 * h10['step'] / length:.1f} %) with the step noise; noise / slope at a 0.001 step: {e1['micro']:.1f} and {e1['step']:.1f}")
    return out


def print_ladder() -> None:
    """Central differences of the reference at every step h that has both a plus and a minus copy (stages noise and ladder): the estimate
    should settle on the sweep's slope as h grows past the noise, until the curvature takes it away again."""
    stages = [s for s in ("noise", "ladder") if has_stage(ss.CSV_J, s)]
    rows = [r for s in stages for r in ss.read_rows(ss.CSV_J, s) if r["feasible"]]
    ref = next(r for r in rows if json.loads(r["overrides"]) == {})
    print("central difference (f(x0 + h) - f(x0 - h)) / (2 h) of the reference, per step h (torque in N.m per unit, ripple in %-points per unit):")
    for v in ss.VARS:
        by = {}
        for r in rows:
            ov = json.loads(r["overrides"])
            if len(ov) == 1 and v in ov:
                by[round(ov[v] - getattr(REF, v), 6)] = r
        line = []
        for h in sorted({abs(k) for k in by}):
            if h in by and -h in by:
                dT = (by[h]["T_mean"] - by[-h]["T_mean"]) / (2 * h)
                dR = (by[h]["ripple_pct"] - by[-h]["ripple_pct"]) / (2 * h)
                line.append(f"h {h:g}: T {dT:+.4f} rip {dR:+.3f}")
        print(f"   {v:>13} ({ss.UNIT[v]}): " + "; ".join(line))


# ---------------------------------------------------------------------------------------------- coarse against fine
def kendall(a: np.ndarray, b: np.ndarray) -> float:
    """Kendall's tau of two rankings: +1 the same order, -1 the opposite order (pairs of equal values count as 0)."""
    n, s = len(a), 0
    for i in range(n):
        for j in range(i + 1, n):
            s += np.sign(a[i] - a[j]) * np.sign(b[i] - b[j])
    return float(s / (n * (n - 1) / 2)) if n > 1 else float("nan")


def print_fine() -> dict:
    """Every fine-mesh row next to its coarse twin (same overrides): the offsets, and whether the order of the points is the same."""
    fine = [r for r in ss.read_rows(ss.CSV_FINE, "fine") if r["feasible"]]
    s_T, s_R = sigmas()
    out = {}
    all_dT, all_dR = [], []
    for v in ss.VARS:
        coarse = {round(r[v], 6): r for r in load_sweep(v, "J")}
        mine = sorted([r for r in fine if json.loads(r["overrides"]) == {} or v in json.loads(r["overrides"])], key=lambda r: r[v])
        print(f"== {v} [{ss.UNIT[v]}]: coarse against fine")
        print(f"{'value':>9} {'T coarse':>9} {'T fine':>9} {'dT %':>7} {'rip coarse':>10} {'rip fine':>9} {'d rip':>7} {'gamma c':>8} {'gamma f':>8} {'nodes':>6}")
        pairs = []
        for r in mine:
            c = coarse.get(round(r[v], 6))
            if c is None:
                continue
            pairs.append((c, r))
            print(f"{r[v]:9.4f} {c['T_mean']:9.5f} {r['T_mean']:9.5f} {100 * (r['T_mean'] / c['T_mean'] - 1):+7.2f} {c['ripple_pct']:10.3f} "
                  f"{r['ripple_pct']:9.3f} {r['ripple_pct'] - c['ripple_pct']:+7.3f} {c['gamma']:8.2f} {r['gamma']:8.2f} {r['nodes']:6.0f}")
        cT, fT = np.array([p[0]["T_mean"] for p in pairs]), np.array([p[1]["T_mean"] for p in pairs])
        cR, fR = np.array([p[0]["ripple_pct"] for p in pairs]), np.array([p[1]["ripple_pct"] for p in pairs])
        dT, dR = 100 * (fT / cT - 1), fR - cR
        all_dT += list(dT)
        all_dR += list(dR)
        best_c, best_f = int(np.argmin(cR)), int(np.argmin(fR))
        out[v] = {"tau_T": kendall(cT, fT), "tau_R": kendall(cR, fR), "dT_pct": dT.tolist(), "dR": dR.tolist(),
                  "best_ripple_coarse_x": float(pairs[best_c][0][v]), "best_ripple_fine_x": float(pairs[best_f][0][v])}
        print(f"   {len(pairs)} points; Kendall tau of the torque order {out[v]['tau_T']:+.2f}, of the ripple order {out[v]['tau_R']:+.2f}; "
              f"lowest ripple: coarse at {out[v]['best_ripple_coarse_x']:.4f}, fine at {out[v]['best_ripple_fine_x']:.4f}; "
              f"offset spread across the points: torque {np.ptp(dT):.2f} %, ripple {np.ptp(dR):.2f} %-points")
    print(f"all {len(all_dT)} pairs: torque offset mean {np.mean(all_dT):+.2f} %, sd {np.std(all_dT, ddof=1):.2f} %; "
          f"ripple offset mean {np.mean(all_dR):+.3f} %-points, sd {np.std(all_dR, ddof=1):.3f} (coarse sigma: torque {100 * s_T / 1.1217:.2f} %, ripple {s_R:.3f})")
    return out


# ---------------------------------------------------------------------------------------------- derived numbers for the report
def print_derived() -> None:
    """What J hides, what AT chooses, and the room the constraints leave: three small sets of numbers built from the sweeps and the scan."""
    print("torque per root copper loss, T / sqrt(Pcu_rel), under J (the first-order torque at the reference's copper loss; it is exact where T goes with the current):")
    for v in ss.VARS:
        rows = load_sweep(v, "J")
        x, T, P = column(rows, v), column(rows, "T_mean"), column(rows, "Pcu_rel")
        k = T / np.sqrt(P)
        i, r = int(np.argmax(k)), int(np.argmin(np.abs(x - getattr(REF, v))))
        print(f"   {v:>13}: highest {k[i]:.4f} at {x[i]:.4f}; reference {k[r]:.4f}; first point {k[0]:.4f}, last point {k[-1]:.4f}; "
              f"points within 0.5 % of the highest: {x[k >= 0.995 * k[i]].min():.4f}..{x[k >= 0.995 * k[i]].max():.4f}")
    print("what an optimizer would pick under AT (100 ampere-turns) against J, for the variables that move the slot:")
    for v in ss.VARS:
        if not has_stage(ss.CSV_AT, v + "_AT"):
            continue
        A, Jr = load_sweep(v, "AT"), load_sweep(v, "J")
        ref = A[int(np.argmin(np.abs(column(A, v) - getattr(REF, v))))]
        for rows, tag in ((A, "AT"), (Jr, "J")):
            T, R = column(rows, "T_mean"), column(rows, "ripple_pct")
            for what, k in (("highest torque", int(np.argmax(T))), ("lowest ripple", int(np.argmin(R)))):
                r = rows[k]
                print(f"   {v:>13} {tag:>2} {what:>14}: at {r[v]:.4f}: torque {r['T_mean']:.5f} ({100 * (r['T_mean'] / ref['T_mean'] - 1):+.2f} % on the reference), "
                      f"ripple {r['ripple_pct']:.3f} ({r['ripple_pct'] - ref['ripple_pct']:+.3f}), {r['amps']:.2f} A, J_rel {r['J_rel']:.3f}, Pcu_rel {r['Pcu_rel']:.3f}")
    print("room the constraints leave (from screen_scan.csv; the variable alone, the others at the reference):")
    import csv
    with open(ss.folder() / "screen_scan.csv", newline="", encoding="utf-8") as f:
        scan = [r for r in csv.DictReader(f) if r["feasible"] == "1"]
    for var, key, targets in (("slot_half_deg", "tooth_min", (4.0, 3.5, 3.0, 2.5, 2.0)), ("r_wall", "yoke", (5.0, 4.5, 4.0, 3.5, 3.0))):
        pts = sorted((float(r["x"]), float(r[key])) for r in scan if r["var"] == var)
        x, y = np.array([a for a, _ in pts]), np.array([b for _, b in pts])
        for t in targets:
            ok = x[y >= t]
            print(f"   {var} with {key} >= {t:g} mm: up to {ok.max():.3f} (reference {getattr(REF, var)}; reference {key} {y[int(np.argmin(np.abs(x - getattr(REF, var))))]:.3f} mm)")
    print("the opening's room against the slot (open_vs_slot): slot_open <= 2 r_si sin(slot_half_deg - 0.3 deg):")
    for sh in (1.5521, 1.8, 2.0, 2.57, 3.0, 3.5, 3.7747):
        print(f"   slot_half_deg {sh}: slot_open up to {2 * REF.r_si * math.sin(math.radians(sh - design.OPEN_MARGIN_DEG)):.3f} mm")


def print_cost() -> None:
    """How much FEMM the study used: candidates solved (rows that were not cached and passed the screen), solves, and the workers' own time."""
    total = 0
    for name in (ss.CSV_J, ss.CSV_AT, ss.CSV_FINE):
        rows = [r for r in ss.read_rows(name) if r["feasible"] and not r["cached"]]
        secs = float(np.nansum([r["seconds"] for r in rows]))
        by_stage = {}
        for r in rows:
            by_stage[r["stage"]] = by_stage.get(r["stage"], 0) + 1
        print(f"{name}: {len(rows)} candidates solved ({23 * len(rows)} FEMM solves), {secs / 60:.1f} worker-minutes; by stage {by_stage}")
        total += len(rows)
    print(f"total {total} candidates solved; rows written {sum(len(ss.read_rows(n)) for n in (ss.CSV_J, ss.CSV_AT, ss.CSV_FINE))} (including cached copies and the cells the screen refused)")


# ---------------------------------------------------------------------------------------------- the 2-D map
def load_map(stage: str = "map") -> dict:
    """The map stage as arrays: x and y values, torque, ripple, current, the screen's reason where a cell was refused. The first variable
    named in a row's overrides is x, the second y (stage_map writes them in that order)."""
    rows = ss.read_rows(ss.CSV_J, stage)
    vx, vy = list(json.loads(rows[0]["overrides"]))
    xs = sorted({round(r[vx], 6) for r in rows})
    ys = sorted({round(r[vy], 6) for r in rows})
    shape = (len(ys), len(xs))
    m = {"vx": vx, "vy": vy, "xs": np.array(xs), "ys": np.array(ys), "T": np.full(shape, np.nan), "R": np.full(shape, np.nan),
         "amps": np.full(shape, np.nan), "gamma": np.full(shape, np.nan), "ok": np.zeros(shape, bool), "reason": {}}
    for r in rows:
        i, j = xs.index(round(r[vx], 6)), ys.index(round(r[vy], 6))
        if r["feasible"]:
            m["ok"][j, i] = True
            m["T"][j, i], m["R"][j, i], m["amps"][j, i], m["gamma"][j, i] = r["T_mean"], r["ripple_pct"], r["amps"], r["gamma"]
        else:
            m["reason"][(j, i)] = r["reason"]
    return m


def grid_minima(Z: np.ndarray) -> list[tuple]:
    """The local minima of a 2-D array with holes (nan), against the 8 neighbours, each with its prominence: the lowest pass that has to be
    crossed to reach a cell with a lower value, less the cell's own value (inf for the lowest cell of all). The pass is found by a search that
    always extends the path with the lowest highest-point-so-far (a priority queue). Returns [((j, i), value, prominence)]."""
    import heapq
    ny, nx = Z.shape
    nb = [(dj, di) for dj in (-1, 0, 1) for di in (-1, 0, 1) if (dj, di) != (0, 0)]
    inside = lambda a, b: 0 <= a < ny and 0 <= b < nx and not np.isnan(Z[a, b])
    out = []
    for j in range(ny):
        for i in range(nx):
            if np.isnan(Z[j, i]) or any(inside(j + dj, i + di) and Z[j + dj, i + di] <= Z[j, i] for dj, di in nb):
                continue
            v0, best, heap, prom = Z[j, i], {(j, i): Z[j, i]}, [(Z[j, i], j, i)], math.inf
            while heap:
                key, a, b = heapq.heappop(heap)
                if best[(a, b)] < key:
                    continue
                if Z[a, b] < v0:
                    prom = key - v0
                    break
                for dj, di in nb:
                    c, d = a + dj, b + di
                    if inside(c, d) and max(key, Z[c, d]) < best.get((c, d), math.inf):
                        best[(c, d)] = max(key, Z[c, d])
                        heapq.heappush(heap, (best[(c, d)], c, d))
            out.append(((j, i), float(v0), float(prom)))
    return out


def additive_share(Z: np.ndarray) -> float:
    """The share of a map's variance that is NOT the sum of a function of x and a function of y: fit Z = a(x) + b(y) by least squares on the
    cells that exist and take 1 - R^2. 0 means the two variables act separately; the larger it is, the more they interact."""
    ny, nx = Z.shape
    cells = [(j, i) for j in range(ny) for i in range(nx) if not np.isnan(Z[j, i])]
    A = np.zeros((len(cells), nx + ny))
    for k, (j, i) in enumerate(cells):
        A[k, i] = A[k, nx + j] = 1.0
    y = np.array([Z[j, i] for j, i in cells])
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(np.sum((y - A @ coef) ** 2) / np.sum((y - y.mean()) ** 2))


def print_map() -> dict:
    """The map's findings: where the lowest ripple and the highest torque are, how many valleys, the valley floor row by row and column by
    column, how much of the variance is interaction, and the best cell once the opening is held to a buildable minimum."""
    m = load_map()
    s_T, s_R = sigmas()
    vx, vy, xs, ys = m["vx"], m["vy"], m["xs"], m["ys"]
    n_ok = int(m["ok"].sum())
    print(f"map {vx} (x, {len(xs)} values {xs[0]:.4f}..{xs[-1]:.4f}) by {vy} (y, {len(ys)} values {ys[0]:.4f}..{ys[-1]:.4f}): {n_ok} cells scored, "
          f"{m['ok'].size - n_ok} refused by the screen ({sorted(set(r.split(' +')[0] for r in m['reason'].values()))})")
    for name, Z in (("ripple %", m["R"]), ("T_mean N.m", m["T"])):
        j, i = np.unravel_index(np.nanargmin(Z), Z.shape)
        j2, i2 = np.unravel_index(np.nanargmax(Z), Z.shape)
        print(f"   {name}: lowest {Z[j, i]:.4f} at {vx} {xs[i]:.4f}, {vy} {ys[j]:.4f}; highest {Z[j2, i2]:.4f} at {vx} {xs[i2]:.4f}, {vy} {ys[j2]:.4f}; "
              f"range {np.nanmin(Z):.4f}..{np.nanmax(Z):.4f}")
    tol = SIGMA_K * step_noise()[1]
    mins = grid_minima(m["R"])
    print(f"   ripple valleys (local minima of the 8-neighbour grid): {len(mins)} raw; prominence above {SIGMA_K:g} x {step_noise()[1]:.3f} = {tol:.3f} %-points:")
    for (j, i), v, p in sorted(mins, key=lambda t: t[1]):
        print(f"      {vx} {xs[i]:.4f}, {vy} {ys[j]:.4f}: ripple {v:.3f} %, prominence {'inf (the lowest cell)' if math.isinf(p) else f'{p:.3f}'}"
              f"{'' if p > tol else '  (not above the noise)'}")
    print(f"   interaction share of the ripple's variance (what a(x) + b(y) cannot explain): {100 * additive_share(m['R']):.1f} %; "
          f"of the torque's: {100 * additive_share(m['T']):.1f} %")
    print(f"   ripple valley floor, for each {vy}: the {vx} with the lowest ripple")
    for j in range(len(ys) - 1, -1, -1):
        row = m["R"][j]
        i = int(np.nanargmin(row))
        i_hi = int(np.nanargmin(np.where(xs >= 1.4, row, np.nan))) if np.any((xs >= 1.4) & m["ok"][j]) else None
        print(f"      {vy} {ys[j]:.4f}: lowest {row[i]:.3f} % at {vx} {xs[i]:.4f}" +
              (f"; with {vx} >= 1.4 mm: {row[i_hi]:.3f} % at {xs[i_hi]:.4f}" if i_hi is not None else ""))
    print(f"   ripple valley floor, for each {vx}: the {vy} with the lowest ripple")
    for i in range(len(xs)):
        col = m["R"][:, i]
        j = int(np.nanargmin(col))
        print(f"      {vx} {xs[i]:.4f}: lowest {col[j]:.3f} % at {vy} {ys[j]:.4f} (torque {m['T'][j, i]:.4f} N.m, {m['amps'][j, i]:.2f} A)")
    mask = m["ok"] & (xs[None, :] >= 1.4)
    Zr = np.where(mask, m["R"], np.nan)
    j, i = np.unravel_index(np.nanargmin(Zr), Zr.shape)
    print(f"   best ripple with {vx} >= 1.4 mm: {Zr[j, i]:.3f} % at {vx} {xs[i]:.4f}, {vy} {ys[j]:.4f} (torque {m['T'][j, i]:.4f} N.m, {m['amps'][j, i]:.2f} A, "
          f"gamma {m['gamma'][j, i]:.1f})")
    return m


def step_noise() -> tuple[float, float]:
    """The roughness a sweep shows between neighbouring points (tens to hundreds of micrometres apart), as (torque, ripple): for each of the
    three sweeps with the smallest swing of the quantity, the robust spread of how far each point lies from the straight line between its
    two neighbours (scaled so that pure noise of sigma gives sigma); the median of the three. The flattest sweeps are used because a curve
    that bends adds to the spread, and a flat one has nothing but the mesh in it."""
    spreads = {"T_mean": [], "ripple_pct": []}
    for scoring in ("J", "AT"):
        for v in ss.VARS:
            if scoring == "AT" and not has_stage(ss.CSV_AT, v + "_AT"):
                continue
            rows = load_sweep(v, scoring)
            x = column(rows, v)
            for m, store in spreads.items():
                y = column(rows, m)
                r = []
                for i in range(1, len(x) - 1):
                    w = (x[i + 1] - x[i]) / (x[i + 1] - x[i - 1])
                    r.append((y[i] - (w * y[i - 1] + (1 - w) * y[i + 1])) / math.sqrt(1 + w * w + (1 - w) ** 2))
                r = np.array(r)
                store.append((float(np.ptp(y)), 1.4826 * float(np.median(np.abs(r - np.median(r))))))
    out = []
    for m in ("T_mean", "ripple_pct"):
        flattest = sorted(spreads[m])[:3]
        out.append(float(np.median([s for _, s in flattest])))
    return out[0], out[1]


# ---------------------------------------------------------------------------------------------- the figures
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
C_J, C_AT = "#2a78d6", "#eb6834"                    # categorical slots 1 and 2 of the palette (validated: light surface, CVD and normal-vision floors pass)
BLUES = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
LIMIT_TEXT = {"open_vs_slot": "opening reaches the\nslot wall (open_vs_slot)", "tooth_min": "tooth under 2 mm\n(tooth_min)",
              "yoke": "yoke under 3 mm\n(yoke)", "shoe_taper": "shoe too flat\n(shoe_taper)",
              "drawing": "the two coil layers cannot\nbe split (the drawing fails)", "wall_length": "walls too short\n(wall_length)"}
OPEN_TEXT = {"slot_open": "no limit:\nwidth 0", "h_open": "no limit:\nthickness 0"}
TITLE = {"slot_half_deg": "Half angle of the slot body", "slot_open": "Slot opening at the bore", "h_open": "Straight part of the opening (tooth tip height)",
         "r_shoe": "Radius where the shoe ends", "r_wall": "Radius where the slot walls end (slot depth)"}


def style() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 9, "axes.facecolor": SURFACE, "figure.facecolor": SURFACE,
                         "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.8,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID,
                         "grid.linewidth": 0.6, "grid.linestyle": "-", "savefig.facecolor": SURFACE})


def draw_variable(var: str, folder) -> str:
    """One figure per variable: torque, ripple and (where AT ran) the copper loss against the variable, J and AT, the feasible interval's ends
    shaded out and named, the reference marked, the fine-mesh re-scores as open diamonds, and the noise floor as a bar in the left margin."""
    style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    iv = ss.intervals()[var]
    J = load_sweep(var, "J")
    AT = load_sweep(var, "AT") if has_stage(ss.CSV_AT, var + "_AT") else None
    fine = [r for r in ss.read_rows(ss.CSV_FINE, "fine") if r["feasible"] and (json.loads(r["overrides"]) == {} or var in json.loads(r["overrides"]))] \
        if has_stage(ss.CSV_FINE, "fine") else []
    rows = 3 if AT else 2
    handles = [Line2D([0], [0], color=C_J, lw=1.0, marker="o", ms=4.0, mec=SURFACE, mew=1.0)]
    labels = ["J: reference current density held (I = 10 A x slot area / reference)"]
    if AT:
        handles.append(Line2D([0], [0], color=C_AT, lw=1.0, marker="o", ms=4.0, mec=SURFACE, mew=1.0))
        labels.append("AT: 100 ampere-turns held (10 A whatever the slot)")
    if fine:
        handles.append(Line2D([0], [0], ls="none", marker="D", ms=5.5, mfc=SURFACE, mec=C_J, mew=1.3))
        labels.append("J, fine mesh")
    handles.append(Line2D([0], [0], ls="none", mfc="none", mec=INK, marker="o", ms=8, mew=1.1))
    labels.append("highest torque and lowest ripple of the J sweep: where an optimizer would go")
    handles.append(Line2D([0], [0], ls="none", color=INK2, marker="|", ms=10, mew=1.3))
    labels.append("noise between neighbouring points (1 sigma, drawn to scale)")
    panel_h, gap, foot = 2.55, 0.16, 0.62                                           # inches: one panel, between panels, under the last one
    head = 0.80 + 0.205 * len(labels)                                                # inches above the first panel: title, subtitle, the legend's rows
    height = head + rows * panel_h + (rows - 1) * gap + foot
    fig, axes = plt.subplots(rows, 1, sharex=True, figsize=(7.6, height), dpi=150)
    fig.subplots_adjust(top=1 - head / height, bottom=foot / height, left=0.12, right=0.985, hspace=gap / panel_h)
    lo, hi = iv["lo"]["x"], iv["hi"]["x"]
    pad = 0.07 * (hi - lo)
    x_ref = getattr(REF, var)
    s_T, s_R = step_noise()
    panels = [("T_mean", "Mean torque [N.m]", s_T), ("ripple_pct", "Ripple [% of mean torque]", s_R)]
    if AT:
        panels.append(("Pcu_rel", "Copper loss, relative to the reference", 0.0))
    for ax, (col, label, sig) in zip(axes, panels):
        ax.axvspan(lo - pad, lo, color=GRID, alpha=0.55, lw=0)
        ax.axvspan(hi, hi + pad, color=GRID, alpha=0.55, lw=0)
        ax.axvline(x_ref, color=INK2, lw=0.8, alpha=0.55)
        for rowsx, color in ((AT, C_AT), (J, C_J)):                               # J on top: it is the study's scoring
            if rowsx:
                ax.plot(column(rowsx, var), column(rowsx, col), color=color, lw=1.0, marker="o", ms=4.0, mec=SURFACE, mew=1.0, zorder=3)
        if col != "Pcu_rel" and fine:
            ax.plot([r[var] for r in fine], [r[col] for r in fine], ls="none", marker="D", ms=5.5, mfc=SURFACE, mec=C_J, mew=1.3, zorder=4)
        if col == "Pcu_rel":
            ax.axhline(1.0, color=INK2, lw=0.8, alpha=0.55)
        if AT and col == "Pcu_rel":                                                 # two series: the legend, and a name at the end of each line where they part
            for rowsx, name in ((J, "J"), (AT, "AT")):
                ax.annotate(name, (column(rowsx, var)[-1], column(rowsx, col)[-1]), xytext=(-6, 6), textcoords="offset points", ha="right",
                            va="bottom", fontsize=8, color=INK2, zorder=6)
        if col in ("T_mean", "ripple_pct"):                                          # a ring on the extreme of the study's scoring: where an optimizer would go
            xj, yj = column(J, var), column(J, col)
            k = int(np.argmax(yj) if col == "T_mean" else np.argmin(yj))
            ax.plot([xj[k]], [yj[k]], ls="none", marker="o", ms=10, mfc="none", mec=INK, mew=1.1, zorder=6)
        ax.set_ylabel(label)
        ax.set_xlim(lo - pad, hi + pad)
        ax.tick_params(length=0)
        ax.grid(axis="x", visible=False)
    fig.canvas.draw()
    for ax, (col, label, sig) in zip(axes, panels):
        if sig:                                                                     # one noise bar per panel, in the empty zone left of the interval
            lo_y, hi_y = ax.get_ylim()
            ax.errorbar([lo - pad / 2], [lo_y + 0.16 * (hi_y - lo_y)], yerr=[sig], fmt="none", ecolor=INK2, elinewidth=1.2, capsize=2.5, capthick=1.2, zorder=5)
    for side, x_text, e in (("lo", lo - pad / 2, iv["lo"]), ("hi", hi + pad / 2, iv["hi"])):
        text = OPEN_TEXT[var] if e["kind"] == "open" else ", ".join(LIMIT_TEXT.get(k, k) for k in e["limit"])
        axes[0].text(x_text, 0.55, text, rotation=90, va="center", ha="center", fontsize=7.5, color=INK2, transform=axes[0].get_xaxis_transform())
    axes[0].text(x_ref, 1.02, "reference", fontsize=7.5, color=INK2, ha="center", va="bottom", transform=axes[0].get_xaxis_transform())
    axes[-1].set_xlabel(f"{var} [{ss.UNIT[var]}]")
    fig.text(0.012, 1 - 0.10 / height, f"{TITLE[var]}: torque and ripple", fontsize=12.5, fontweight="bold", color=INK, va="top")
    fig.text(0.012, 1 - 0.43 / height, "Band-free rotor, other variables at the reference; shaded = beyond the feasible interval", fontsize=8.5, color=INK2, va="top")
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.62 / height), frameon=False, fontsize=8, labelcolor=INK2, handlelength=2.2,
               labelspacing=0.35, borderaxespad=0.2)
    path = str(folder / f"stator_{var}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def draw_map(folder) -> str:
    """The 2-D map: ripple and torque over the pair, one blue ramp (darker = more), refused cells hatched, the reference ringed, the lowest ripple starred."""
    style()
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle
    m = load_map()
    vx, vy, xs, ys = m["vx"], m["vy"], m["xs"], m["ys"]
    cmap = LinearSegmentedColormap.from_list("blue", BLUES[1:])
    dx, dy = (xs[1] - xs[0]) / 2, (ys[1] - ys[0]) / 2
    xe, ye = np.append(xs - dx, xs[-1] + dx), np.append(ys - dy, ys[-1] + dy)
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.3), dpi=150, sharey=True)
    j_best, i_best = np.unravel_index(np.nanargmin(m["R"]), m["R"].shape)
    for ax, key, title in ((axes[0], "R", "Ripple [% of mean torque]"), (axes[1], "T", "Mean torque [N.m]")):
        Z = np.ma.masked_invalid(m[key])
        mesh = ax.pcolormesh(xe, ye, Z, cmap=cmap, shading="flat")
        for (j, i) in m["reason"]:
            ax.add_patch(Rectangle((xe[i], ye[j]), 2 * dx, 2 * dy, facecolor=GRID, edgecolor=AXIS, hatch="////", lw=0, zorder=2))
        ax.plot([getattr(REF, vx)], [getattr(REF, vy)], marker="o", ms=9, mfc="none", mec=INK, mew=1.6, zorder=6)
        ax.annotate("reference", (getattr(REF, vx), getattr(REF, vy)), xytext=(10, -14), textcoords="offset points", fontsize=8, color=INK, zorder=6)
        if key == "R":
            ax.plot([xs[i_best]], [ys[j_best]], marker="*", ms=13, mfc=INK, mec=SURFACE, mew=1.0, zorder=7)
            ax.annotate(f"lowest {m['R'][j_best, i_best]:.1f} %", (xs[i_best], ys[j_best]), xytext=(11, -4), textcoords="offset points", fontsize=8, color=INK, zorder=7)
            for (j, i), v, prom in grid_minima(m["R"]):                              # the other valleys that the noise cannot make up
                if prom > SIGMA_K * step_noise()[1] and (j, i) != (j_best, i_best):
                    ax.plot([xs[i]], [ys[j]], marker="D", ms=7, mfc="none", mec=INK, mew=1.3, zorder=7)
                    right = i > len(xs) - 3                                              # a label at the right edge would run under the colour bar
                    ax.annotate(f"valley {v:.1f} %", (xs[i], ys[j]), xytext=(6 if right else 0, -13), textcoords="offset points",
                                ha="right" if right else "center", fontsize=8, color=INK, zorder=7)
        else:
            j2, i2 = np.unravel_index(np.nanargmax(m["T"]), m["T"].shape)
            ax.plot([xs[i2]], [ys[j2]], marker="*", ms=13, mfc=SURFACE, mec=INK, mew=1.2, zorder=7)
            ax.annotate(f"highest {m['T'][j2, i2]:.2f}", (xs[i2], ys[j2]), xytext=(11, -12), textcoords="offset points", fontsize=8, color=SURFACE, zorder=7)
        ax.set_title(title, loc="left", fontsize=10.5, color=INK, fontweight="bold")
        ax.set_xlabel(f"{vx} [{ss.UNIT[vx]}]")
        ax.set_xlim(xe[0], xe[-1])
        ax.set_ylim(ye[0], ye[-1])
        ax.grid(False)
        ax.tick_params(length=0)
        cb = fig.colorbar(mesh, ax=ax, shrink=0.9, pad=0.02)
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=0, labelcolor=MUTED)
    axes[0].set_ylabel(f"{vy} [{ss.UNIT[vy]}]")
    fig.suptitle(f"{vx} against {vy}: ripple and torque, other variables at the reference", x=0.012, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    fig.text(0.012, 0.915, "J scoring, 11 x 11 values across both feasible intervals, each cell centred on its values; darker = larger; hatched = refused by "
             "the screen (open_vs_slot); diamond = a ripple valley the noise cannot make up", fontsize=8.5, color=INK2, va="top")
    fig.subplots_adjust(top=0.83, left=0.07, right=0.97, bottom=0.11, wspace=0.12)
    path = str(folder / f"stator_map_{vx}_{vy}.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def draw_all() -> None:
    folder = ss.folder()
    for v in ss.VARS:
        print("wrote", draw_variable(v, folder))
    if has_stage(ss.CSV_J, "map"):
        print("wrote", draw_map(folder))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "tables"
    {"tables": print_tables, "shapes": print_shapes, "noise": print_noise, "ladder": print_ladder, "fine": print_fine, "map": print_map,
     "derived": print_derived, "cost": print_cost, "figures": draw_all}[what]()
