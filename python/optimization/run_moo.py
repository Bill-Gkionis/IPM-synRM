"""run_moo: a multi-objective optimizer run of the study, on design.py's black box, with every evaluation logged and the front read
from the whole log.

What it is
    pymoo proposes, design.Evaluator scores (on one femm_pool held for the whole run), and this file keeps the books:
      the problem    the variables and their bounds (a JSON file or one of the BOXES below), the two objectives of design.objectives()
                     (-mean torque, ripple in % of it), and the constraints: the screen's G values, one that says the row has no
                     scores (refused, or FEMM failed), and the torque floor (mean torque >= 0.8 N.m: the ripple in % of a torque
                     near zero is a ratio of two tiny numbers, and an optimizer minimising it would collect dead rotors; the rotor
                     sweeps found 1.46 % at 0.0047 N.m)
      the algorithm  NSGA-II (pymoo's default operators: simulated binary crossover and polynomial mutation, tournament selection),
                     seeded with designs given (the reference, a tier's knee), the rest of the first population from a Latin
                     hypercube over the box (so the start covers the box evenly instead of clumping)
      the books      every candidate is a row of out\\optimization\\<run>\\evaluations.csv as soon as its generation is back; after
                     each generation a line with the generation's best torque and least ripple, the front's size and its hypervolume
                     against a fixed reference point (so two runs can be compared on one scale)
      the answer     the non-dominated set of ALL feasible evaluations of the run (not only the last population: a good design an
                     optimizer stepped over is still in the log), its knee (design.knee), and both written to front.csv and knee.json

Why the front is read from the log and not from pymoo's result
    pymoo returns the last population's survivors. A design found in generation 3 and lost later through crowding is still a real,
    measured design. The log holds them all, so the front is computed from it (design.nondominated). The demo of 10-28 learned to
    check pymoo against the log the hard way (its first run optimised scrambled scores); here the log IS the result.

The hypervolume (how a run's progress is measured)
    The area of the objective plane that the front dominates, up to a reference point that every useful design beats: here
    (-0.8 N.m, 40 % ripple). Bigger is better. A front that moves outward, or fills a gap, grows it; noise does not.

Run (from the Python to FEMM folder)
    uv run python -m optimization.run_moo --run rotor_nsga2 --box rotor --scoring AT --pop 40 --gens 30 --seed 1
    options: --box NAME or --bounds FILE.json ({"name": [low, high], ...}); --seeds FILE.json (a list of candidate dicts);
             --fix JSON (overrides held for every candidate, e.g. a rotor tier's knee while the stator tier runs)
Writes out\\optimization\\<run>\\: evaluations.csv, progress.csv, front.csv, knee.json, front.png.

Self-check (the answers follow).
  1. A run's last population has 40 designs and its front 12. The log's front has 15. Where did the other 3 come from?
  2. Why does a candidate that has no scores (refused, or FEMM failed on it) get a constraint violation and not a bad objective value?
Answers.
  1. From earlier generations: designs that were non-dominated when measured but were later crowded out of the population (NSGA-II
     keeps the most spread-out members of a front when it must cut). They are real measurements, so the front keeps them.
  2. A bad objective value is a claim about the design ("little torque"), and the failure says nothing about the design's torque.
     A violated constraint says "not admissible" without inventing a number, and NSGA-II then ranks it below every feasible design.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

import femm_queue
from optimization import design

BOXES: dict[str, dict] = {}          # filled from out\optimization\boxes.json when the sweeps have set the bounds
HV_REF = (-0.8, 40.0)                # (-T [N.m], ripple [%]): every useful design beats it
TORQUE_FLOOR = 0.8                   # N.m: a constraint, not an objective. Without it the ripple in % rewards dead rotors (the rotor
                                     # sweeps: 1.46 % ripple at 0.0047 N.m, a ratio of two tiny numbers, not a smooth motor)


def hypervolume(F, ref=HV_REF) -> float:
    """The dominated area of a 2-objective front F (both minimised) up to ref: sort by the first objective, add rectangles."""
    F = np.asarray([f for f in F if f[0] < ref[0] and f[1] < ref[1]], float)
    if not len(F):
        return 0.0
    F = F[design.nondominated(F)]
    F = F[np.argsort(F[:, 0])]
    area, prev_y = 0.0, ref[1]
    for x, y in F:
        area += (ref[0] - x) * (prev_y - y)
        prev_y = y
    return float(area)


def latin_hypercube(n: int, lo, hi, rng) -> np.ndarray:
    """n points in the box [lo, hi]: each variable's range cut into n equal strata, one point in each, the strata paired at random."""
    d = len(lo)
    strata = np.array([rng.permutation(n) for _ in range(d)]).T
    return lo + (strata + rng.random((n, d))) / n * (hi - lo)


def load_box(args) -> dict:
    if args.bounds:
        return json.loads(Path(args.bounds).read_text(encoding="utf-8"))
    boxes = json.loads((design.out_dir("") / "boxes.json").read_text(encoding="utf-8")) if not BOXES else BOXES
    return boxes[args.box]


def main(argv=None) -> int:
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.core.problem import Problem
    from pymoo.optimize import minimize

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--box", default="")
    ap.add_argument("--bounds", default="")
    ap.add_argument("--scoring", default="P", choices=("AT", "J", "P"))
    ap.add_argument("--pop", type=int, default=40)
    ap.add_argument("--gens", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--seeds", default="", help="JSON file: a list of candidate dicts for the first population")
    ap.add_argument("--fix", default="{}", help="JSON: overrides held for every candidate")
    args = ap.parse_args(argv)

    box = load_box(args)
    names = list(box)
    xl, xu = np.array([box[n][0] for n in names], float), np.array([box[n][1] for n in names], float)
    fix = json.loads(args.fix)
    folder = design.out_dir(args.run)
    (folder / "setup.json").write_text(json.dumps({"args": vars(args), "box": box, "fix": fix, "hv_ref": HV_REF}, indent=2), encoding="utf-8")
    progress_path = folder / "progress.csv"
    with open(progress_path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow(["gen", "evaluated", "feasible", "best_T", "least_ripple_pct", "front_size", "hypervolume", "seconds"])

    seeds = json.loads(Path(args.seeds).read_text(encoding="utf-8")) if args.seeds else []
    rng = np.random.default_rng(args.seed)
    X0 = [[float(np.clip(s.get(n, (xl[i] + xu[i]) / 2), xl[i], xu[i])) for i, n in enumerate(names)] for s in seeds]
    lhs = latin_hypercube(args.pop - len(X0), xl, xu, rng) if args.pop > len(X0) else np.empty((0, len(names)))
    X0 = np.vstack([np.array(X0, float).reshape(-1, len(names)), lhs])
    log_F: list = []
    started = time.perf_counter()

    with femm_queue.femm_pool(f"optimization: {args.run}", task_timeout_s=design.TASK_TIMEOUT_S) as pool:
        ev = design.Evaluator(pool, run=args.run, scoring=args.scoring)

        class Study(Problem):
            def __init__(self):
                super().__init__(n_var=len(names), n_obj=2, n_ieq_constr=len(design.G_NAMES) + 2, xl=xl, xu=xu)
                self.gen = 0

            def _evaluate(self, X, out, *a, **k):
                rows = ev([{**fix, **dict(zip(names, map(float, x)))} for x in X], stage=args.run, gen=self.gen)
                F, G = [], []
                for r in rows:
                    ok = r["feasible"] and "T_mean" in r
                    F.append(list(design.objectives(r)) if ok else [0.0, 100.0])
                    G.append([r["G"][g] for g in design.G_NAMES] + [0.0 if ok else 1.0,          # any row without scores is inadmissible
                                                                      TORQUE_FLOOR - r["T_mean"] if ok else 0.0])
                    if ok and r["T_mean"] >= TORQUE_FLOOR:
                        log_F.append(design.objectives(r))
                out["F"], out["G"] = np.array(F, float), np.array(G, float)    # arrays: a list would be read as columns
                feas = [r for r in rows if r["feasible"] and "T_mean" in r and r["T_mean"] >= TORQUE_FLOOR]
                front = np.asarray(log_F)[design.nondominated(log_F)] if log_F else []
                line = [self.gen, len(rows), len(feas), max((r["T_mean"] for r in feas), default=float("nan")),
                        min((r["ripple_pct"] for r in feas), default=float("nan")), len(front), hypervolume(log_F),
                        round(time.perf_counter() - started, 1)]
                with open(progress_path, "a", newline="", encoding="utf-8") as f:
                    csv.writer(f).writerow(line)
                print(f"gen {self.gen:3d}: {len(feas):3d}/{len(rows)} feasible, best T {line[3]:.4f} N.m, least ripple {line[4]:.2f} %, "
                      f"front {line[5]}, HV {line[6]:.4f}, {line[7]:.0f} s", flush=True)
                self.gen += 1

        minimize(Study(), NSGA2(pop_size=args.pop, sampling=X0), ("n_gen", args.gens), seed=args.seed, verbose=False)

    write_front(folder, names)
    print(f"done in {time.perf_counter() - started:.0f} s; front and knee in {folder}")
    return 0


def write_front(folder: Path, names: list[str]) -> dict:
    """The non-dominated set of every feasible row of the run's log, and its knee: front.csv, knee.json, front.png."""
    import pandas as pd
    log = pd.read_csv(folder / "evaluations.csv")
    ok = log[(log.feasible == True) & log.T_mean.notna()].copy()      # noqa: E712 (pandas wants the comparison)
    ok = ok[ok.T_mean >= TORQUE_FLOOR]
    ok = ok.drop_duplicates(subset=["overrides"])
    F = np.c_[-ok.T_mean.to_numpy(float), ok.ripple_pct.to_numpy(float)]
    front = ok[design.nondominated(F)].sort_values("T_mean")
    front.to_csv(folder / "front.csv", index=False)
    Ff = np.c_[-front.T_mean.to_numpy(float), front.ripple_pct.to_numpy(float)]
    k = front.iloc[design.knee(Ff)]
    knee = {"overrides": json.loads(k.overrides), "T_mean": float(k.T_mean), "ripple_pct": float(k.ripple_pct),
            "ripple_pp": float(k.ripple_pp), "gamma": float(k.gamma), "pf": float(k.pf), "bridge_mm": float(k.bridge_mm),
            "front_size": int(len(front)), "feasible_evaluations": int(len(ok)), "hypervolume": hypervolume(F)}
    (folder / "knee.json").write_text(json.dumps(knee, indent=2), encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.scatter(ok.T_mean, ok.ripple_pct, s=8, color="#b8c2cc", label=f"every feasible design ({len(ok)})")
    ax.plot(front.T_mean, front.ripple_pct, "o-", ms=4, lw=1.2, color="#2f6fb5", label=f"the front ({len(front)})")
    ax.scatter([k.T_mean], [k.ripple_pct], s=110, facecolors="none", edgecolors="#c0392b", lw=1.6, label="the knee", zorder=4)
    ax.set_xlabel("mean torque at the scoring point [N·m]")
    ax.set_ylabel("ripple, peak-to-peak [% of the mean]")
    ax.set_title(f"{folder.name}: every design, the front and its knee", fontsize=10, loc="left")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(folder / "front.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return knee


if __name__ == "__main__":
    raise SystemExit(main())
