"""A small but real NSGA-II run on the reference rotor: pymoo proposes, motorgeo draws, shapely inspects, FEMM measures.

Variables: magnet centre (mag_x, mag_y), half V angle, web. Objectives: highest mean torque at 10 A and the
candidate's OWN best current angle; lowest torque ripple. Constraints: the shapely rules (feasibility first).
Every evaluation becomes a CSV row. The best compromise (the "knee") then goes back to CAD as a DXF and a
SolidWorks equations file, and the DXF is read back and compared (regression test).

Run: uv run python optimize_demo.py   (from this folder; about 12 min with FEMM 4.2, 8 workers, coarse mesh; FEMM 5.0:
about 98 s; the numbers are in the README, Expected numbers)
This is a demonstration of the loop, not the real study: 8 candidates x 7 generations, coarse mesh.
"""
import csv
import math
import os
import time
import uuid
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import femm_session
import motorgeo as mg

OUT = os.path.abspath("out/opt")
NAMES = ["mag_y", "mag_x", "v_half_deg", "web"]
LOW, HIGH = [5.0, 38.9, 58.0, 1.0], [6.5, 39.6, 62.0, 1.8]
BASE = [6.0, 39.4136, 60.0, 1.5]                      # the reference design, seeded into the first generation
POP, GENS, WORKERS = 8, 7, 8                          # one candidate per worker; 8 workers keep the desktop responsive
I = 10.0                                              # A peak, x 10 turns = 100 A-turns
GAMMAS = (25.0, 45.0, 65.0)                           # 3 current angles, a parabola finds the best one
MEAN_POS = range(0, 60, 10)                           # 6 positions: their average cancels the ripple
RIPPLE_POS = range(0, 60, 4)                          # 15 positions over the 60 deg ripple period
LIMIT = 1.0                                           # mm, every rib, bridge and web


def evaluate(x):
    """One candidate -> one CSV row. Runs in a worker process with its own FEMM and its own files."""
    t0 = time.perf_counter()
    row = dict(x)
    p = mg.Params(**x)
    try:
        regions = mg.build(p)
    except mg.GeometryError as e:                     # can't even be drawn
        return {**row, "feasible": False, "reason": str(e), "G": [1.0, 1.0, 1.0, 1.0], "seconds": time.perf_counter() - t0}
    m, problems = mg.check(p, regions)
    G = [LIMIT - m["centre_rib"], LIMIT - m["q_rib"], LIMIT - m["outer_bridge"], LIMIT - m["web"]]   # <= 0 is allowed
    row.update({"centre_rib": m["centre_rib"], "q_rib": m["q_rib"], "outer_bridge": m["outer_bridge"], "web_mm": m["web"]})
    if problems:                                      # drawn, but too thin somewhere: no FEMM
        return {**row, "feasible": False, "reason": "; ".join(problems), "G": G, "seconds": time.perf_counter() - t0}

    import femm
    fem = os.path.join(OUT, "fem", uuid.uuid4().hex + ".fem")        # private file per candidate
    femm_session.open_femm()                          # this worker's own FEMM: own port, below-normal (desktop first)
    try:
        mg.to_femm(p, regions, fem)                   # coarse mesh: fine for ranking
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        cur = lambda g: (-I * math.sin(math.radians(g)), I * math.cos(math.radians(g)))
        T_g = [np.mean([model.solve(th, *cur(g))["T_age"] for th in MEAN_POS]) for g in GAMMAS]
        a2, a1, _ = np.polyfit(GAMMAS, T_g, 2)
        g_best = float(np.clip(-a1 / (2 * a2), 0, 90)) if a2 < 0 else GAMMAS[int(np.argmax(T_g))]
        T = [model.solve(th, *cur(g_best))["T_age"] for th in RIPPLE_POS]
    finally:
        femm_session.close_femm()
        for ext in (".fem", ".ans"):
            if os.path.exists(fem[:-4] + ext):
                os.remove(fem[:-4] + ext)
    return {**row, "feasible": True, "reason": "", "G": G, "gamma": g_best, "T_mean": float(np.mean(T)),
            "ripple": float(max(T) - min(T)), "seconds": time.perf_counter() - t0}


if __name__ == "__main__":
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.problem import Problem
    from pymoo.optimize import minimize

    os.makedirs(os.path.join(OUT, "fem"), exist_ok=True)
    log_path = os.path.join(OUT, "evaluations.csv")
    fields = ["gen", *NAMES, "feasible", "reason", "gamma", "T_mean", "ripple",
              "centre_rib", "q_rib", "outer_bridge", "web_mm", "seconds"]
    with open(log_path, "w", newline="") as f:
        csv.writer(f).writerow(fields)

    class Rotor(Problem):
        def __init__(self, pool):
            super().__init__(n_var=len(NAMES), n_obj=2, n_ieq_constr=4, xl=LOW, xu=HIGH)
            self.pool, self.gen = pool, 0

        def _evaluate(self, X, out, *args, **kwargs):
            rows = list(self.pool.map(evaluate, [dict(zip(NAMES, map(float, x))) for x in X]))
            with open(log_path, "a", newline="") as f:      # every evaluation, feasible or not, is a row
                w = csv.writer(f)
                for r in rows:
                    w.writerow([self.gen, *(r[n] for n in NAMES), r["feasible"], r["reason"], r.get("gamma"),
                                r.get("T_mean"), r.get("ripple"), r.get("centre_rib"), r.get("q_rib"),
                                r.get("outer_bridge"), r.get("web_mm"), round(r["seconds"], 1)])
            ok = [r["feasible"] for r in rows]
            print(f"generation {self.gen}: {sum(ok)}/{len(rows)} feasible, best torque "
                  f"{max([r['T_mean'] for r in rows if r['feasible']], default=float('nan')):.4f} N.m", flush=True)
            self.gen += 1
            # arrays, one row per candidate. A plain list would be read by pymoo as a list of COLUMNS and
            # silently scrambled (found 2026-10-28: the first demo run optimised mixed-up scores)
            out["F"] = np.array([[-r["T_mean"], r["ripple"]] if r["feasible"] else [0.0, 1.0] for r in rows])
            out["G"] = np.array([r["G"] for r in rows])

    rng = np.random.default_rng(7)
    X0 = np.vstack([BASE, rng.uniform(LOW, HIGH, size=(POP - 1, len(NAMES)))])     # seed with the reference design
    t = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        res = minimize(Rotor(pool), NSGA2(pop_size=POP, sampling=X0), ("n_gen", GENS), seed=7, verbose=False)
    print(f"{GENS} generations in {time.perf_counter() - t:.0f} s")

    # self-check: every Pareto design's scores must match its row in the CSV log
    import pandas as pd
    log = pd.read_csv(log_path)
    for x, f in zip(res.X, res.F):
        hit = log[np.isclose(log[NAMES].to_numpy(), x, atol=1e-9).all(axis=1)].iloc[0]
        assert np.isclose(-f[0], hit.T_mean) and np.isclose(f[1], hit.ripple), f"pymoo F {f} != log {hit.T_mean}, {hit.ripple}"
    print("self-check: pymoo's Pareto scores match the evaluation log")

    # the knee: closest to the ideal corner after scaling both objectives to 0..1
    F = res.F
    Fn = (F - F.min(axis=0)) / np.where(np.ptp(F, axis=0) > 0, np.ptp(F, axis=0), 1)
    knee = res.X[int(np.argmin(np.hypot(*Fn.T)))]
    kx = dict(zip(NAMES, map(float, knee)))
    print("Pareto set:")
    for x, f in sorted(zip(res.X, res.F), key=lambda z: z[1][0]):
        print("  " + ", ".join(f"{n} {v:.3f}" for n, v in zip(NAMES, x)) + f" -> T {-f[0]:.4f} N.m, ripple {f[1]:.4f} N.m")
    print("knee:", {k: round(v, 3) for k, v in kx.items()})

    # back to CAD: DXF + equations, then read the DXF back and compare (shapely)
    p = mg.Params(**kx)
    regions = mg.build(p)
    mg.to_dxf(p, regions, os.path.join(OUT, "knee_full.dxf"))
    mg.to_solidworks_equations(p, os.path.join(OUT, "knee_equations.txt"))
    rows = mg.compare_dxf(regions, os.path.join(OUT, "knee_full.dxf"))
    worst = max(rows, key=lambda r: abs(r[2] / r[1] - 1))
    print(f"knee DXF round trip: {len(rows)} regions, worst {worst[0]} {100 * (worst[2] / worst[1] - 1):+.3f} %")
    np.save(os.path.join(OUT, "pareto_X.npy"), res.X)
    np.save(os.path.join(OUT, "pareto_F.npy"), res.F)
