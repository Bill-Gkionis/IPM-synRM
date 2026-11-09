"""screening: Morris's method (SALib) on the study's real box, for the mean torque AND the ripple: which variables matter, and how
nonlinearly or jointly they act, before an optimizer spends its budget.

What it is, and how it differs from the bench's B7
    B7 proved the chain (Python -> SALib -> FEMM) on four rotor dimensions, in a box of +-1.5 % that the bridge guard allowed, on the
    mean torque at a fixed current angle. Two of its four dimensions are the magnet's size, which the study holds fixed (the magnets
    are bought). This run asks the study's question instead:
      the box      every variable of the joint problem, over the bounds the sweeps set (out\\optimization\\boxes.json), with the
                   rotor's depth given as bridge_mm (the mapping layer), so the guard cannot break a trajectory
      the outputs  the mean torque AND the ripple in % (design.objectives), each design at its own MTPA angle (design.score)
      the scoring  scoring="P": the stator's variables change the slot, so the current is set for the reference's copper loss
    and reads, for each output, mu* (the mean size of a variable's effect across its whole range), mu (its signed mean: + or - or
    both), sigma (the spread of the effect inside the box: curvature, or acting together with others) and sigma / mu*.

Morris's method in four lines (SALib 1.6.0's own definitions, read in its source for B7)
    1. Cut each variable's range into p = 4 levels. A trajectory starts at a random grid point and moves one variable at a time by a
       step of 2/3 of its range, every variable once, in a random order: k + 1 designs for k variables.
    2. Each move gives one elementary effect of the variable moved: EE = (y_after - y_before) / step (step in units of the range).
    3. r trajectories give r elementary effects per variable. mu* = mean |EE|, mu = mean EE, sigma = std(EE).
    4. A variable with small mu* can be held fixed. sigma near or above mu* says its effect depends on where the others are.
    Cost: r (k + 1) designs; with k = 8 and r = 12, 108 designs (about 2 minutes of FEMM on 8 workers).

A trajectory needs every one of its designs (an elementary effect is a difference), so every sample point is screened in Python
first, and the run stops before FEMM if any is refused: the box, not the sample, has to change.

Run (from the Python to FEMM folder)
    uv run python -m optimization.screening --box joint --r 12
Writes out\\optimization\\screening\\: evaluations.csv, morris.json (both analyses), morris_T.png and morris_ripple.png (B7's drawing).

Self-check (the answers follow).
  1. A variable has mu* 0.02 N.m and sigma 0.03 N.m for the torque. Can it be held fixed?
  2. Why is the screening run on the ripple as well, when the bench's B7 asked only the torque?
Answers.
  1. Not on that evidence alone. Its average effect is small, but sigma above mu* says the effect changes sign or size across the
     box, so somewhere it may matter a lot (with another variable at one end, say). Look at its elementary effects before fixing it.
  2. Because the ripple is the objective with the rugged landscape (§1.5 of the study note): a variable can do nothing to the mean
     torque and everything to the ripple (the slot opening, for one). Dropping it on the torque's evidence would lose the ripple's
     best lever.
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

import femm_queue
from optimization import design

LEVELS, SEED = 4, 1


def main(argv=None) -> int:
    from SALib.analyze import morris as morris_analyze
    from SALib.sample import morris as morris_sample

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--box", default="joint")
    ap.add_argument("--r", type=int, default=12)
    ap.add_argument("--scoring", default="P", choices=("AT", "J", "P"))
    ap.add_argument("--fix", default='{"h_open": 0.5}', help="overrides held for every design (the study's fixed choices)")
    ap.add_argument("--figures", action="store_true", help="only redraw the two figures from morris.json (no FEMM)")
    args = ap.parse_args(argv)
    if args.figures:
        saved = json.loads((design.out_dir("screening") / "morris.json").read_text(encoding="utf-8"))
        for key, label, unit in (("T_mean", "mean torque", "N.m"), ("ripple_pct", "ripple", "%")):
            title = (f"Morris screening (SALib) of the joint box, {len(saved['box'])} variables, r = {saved['r']}: {saved['designs']} designs, "
                     f"coarse mesh\noutput: {label} [{unit}] at the scoring point, each design at its own MTPA angle")
            draw(saved[key], title, unit, label, design.out_dir("screening") / f"morris_{'T' if key == 'T_mean' else 'ripple'}.png")
        return 0
    box = json.loads((design.out_dir("") / "boxes.json").read_text(encoding="utf-8"))[args.box]
    names = list(box)
    problem = {"num_vars": len(names), "names": names, "bounds": [box[n] for n in names]}
    X = morris_sample.sample(problem, args.r, num_levels=LEVELS, seed=SEED)
    fix = json.loads(args.fix)
    candidates = [{**fix, **dict(zip(names, map(float, x)))} for x in X]
    refused = [(c, s["reason"]) for c in candidates if not (s := design.screen(design.decode(c)))["feasible"]]
    print(f"Morris sample: r = {args.r}, k = {len(names)}: {len(X)} designs; refused in Python: {len(refused)}")
    if refused:
        for c, why in refused[:10]:
            print("  refused:", {k: round(v, 4) for k, v in c.items()}, "-", why)
        print("the box has to change: a trajectory needs every one of its designs")
        return 1
    with femm_queue.femm_pool("optimization: screening", task_timeout_s=design.TASK_TIMEOUT_S) as pool:
        rows = design.Evaluator(pool, run="screening", scoring=args.scoring)(candidates, stage=f"morris r={args.r}")
    failed = [r for r in rows if not (r["feasible"] and "T_mean" in r)]
    if failed:
        print(f"{len(failed)} designs have no scores (FEMM): {[r['reason'] for r in failed[:3]]}")
        return 1
    out = {"box": box, "r": args.r, "levels": LEVELS, "seed": SEED, "designs": len(X), "scoring": args.scoring}
    folder = design.out_dir("screening")
    for key, label, unit in (("T_mean", "mean torque", "N.m"), ("ripple_pct", "ripple", "%")):
        Y = np.array([r[key] for r in rows], float)
        res = morris_analyze.analyze(problem, X, Y, num_levels=LEVELS, seed=SEED, print_to_console=False)
        a = {"names": names, **{k: [float(v) for v in np.ma.filled(np.ma.asarray(res[k], dtype=float), np.nan)]
                                for k in ("mu", "mu_star", "sigma", "mu_star_conf")}}
        a["sigma_over_mu_star"] = [s / m if m else float("nan") for s, m in zip(a["sigma"], a["mu_star"])]
        a["range_of_Y"] = [float(Y.min()), float(Y.max())]
        out[key] = a
        order = np.argsort(-np.asarray(a["mu_star"]))
        print(f"\n{label} [{unit}], from {Y.min():.4f} to {Y.max():.4f} over the sample:")
        for i in order:
            print(f"  {names[i]:14s} mu* {a['mu_star'][i]:8.4f} (+-{a['mu_star_conf'][i]:.4f})  mu {a['mu'][i]:+8.4f}  "
                  f"sigma {a['sigma'][i]:8.4f}  sigma/mu* {a['sigma_over_mu_star'][i]:6.2f}")
        title = (f"Morris screening (SALib) of the joint box, {len(names)} variables, r = {args.r}: {len(X)} designs, coarse mesh\n"
                 f"output: {label} [{unit}] at the scoring point, each design at its own MTPA angle")
        draw(a, title, unit, label, folder / f"morris_{'T' if key == 'T_mean' else 'ripple'}.png")
    (folder / "morris.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {folder}")
    return 0


def draw(a: dict, title: str, unit: str, label: str, path) -> None:
    """B7's two Morris views (bench.b7_salib.draw_morris), with its words for the torque and N.m changed to this output's."""
    import matplotlib.pyplot as plt
    from bench.b7_salib import draw_morris
    fig = draw_morris({k: a[k] for k in ("names", "mu_star", "sigma", "mu_star_conf")}, title)
    for ax in fig.axes:
        for get, put in ((ax.get_xlabel, ax.set_xlabel), (ax.get_ylabel, ax.set_ylabel), (lambda: ax.get_title(loc="left"), lambda t: ax.set_title(t, fontsize=9.5, loc="left"))):   # B7 titles on the left
            text = get()
            if text:
                put(text.replace("N.m", unit).replace("the torque", f"the {label}").replace("of torque", f"of the {label}"))
    fig.savefig(path, dpi=130, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
