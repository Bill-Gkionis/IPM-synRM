"""finalize: the study's last step. The coarse fronts re-scored on the fine mesh, the knee read there, the independent run's knee held
against it, and the chosen design written out as a DXF of the whole motor.

What it does, in order
    1. the fronts      reads front.csv of each run named (the study's tiers, and the independent cross-check), and keeps their union's
                       non-dominated designs on the coarse mesh
    2. the re-score    up to --max of them on the fine mesh (all, if there are fewer; else half evenly spread along the front by
                       torque and half the nearest to the coarse knee, so the knee's neighbourhood is dense), plus the reference design
                       and each run's own knee, so the before and after, and the runs, are compared on the same mesh
    3. the knee        the fine front's knee (design.knee: both objectives scaled 0..1 over the fine front, nearest the ideal corner).
                       That is the chosen design.
    4. the drawing     export.write: the whole motor as a DXF, checked by reading it back, with a preview
Why re-score instead of optimizing again on the fine mesh: the sweeps measured the two meshes against each other (51 designs): the
fine mesh moves the torque by +0.5 % and the ripple by -0.7 to -0.9 points, and it ranks the designs the same way (Spearman 1.00). A
second search on the fine mesh would refine inside a ranking it does not change; re-scoring the front puts the final numbers on the
mesh the study trusts.

Run (from the Python to FEMM folder)
    uv run python -m optimization.finalize --runs joint_nsga2 rotor_nsga2 stator_grid --crosscheck crosscheck --max 24
Writes out\\optimization\\final\\: evaluations.csv (the fine rows), final.json, compare.csv, front_fine.png, and the DXF with its
preview (motor_optimized.dxf, motor_optimized.png).

Self-check (the answers follow).
  1. The coarse knee and the fine knee are different designs. Which one is the answer, and why can that happen when the ranking holds?
  2. Why is the independent run's knee re-scored on the fine mesh too, and not compared by its coarse numbers?
Answers.
  1. The fine knee. The ranking along the front holds, but the knee depends on the front's two ends (the scaling to 0..1), and the
     fine mesh lowers the ripple more on some designs than others (-0.86 +- 0.59 points in the stator sweeps): the ends move a
     little, and the nearest point to the corner can shift by a neighbour.
  2. Two designs compared on two meshes compare the meshes as much as the designs. On one mesh the difference is the design's.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import femm_queue
from optimization import design, export
from optimization.run_moo import TORQUE_FLOOR, hypervolume

REFERENCE = {}                                      # the design as drawn (h_open 0.3: the "before")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--crosscheck", default="crosscheck")
    ap.add_argument("--max", type=int, default=24)
    ap.add_argument("--name", default="motor_optimized")
    args = ap.parse_args(argv)
    root = design.out_dir("")
    fronts = []
    for run in [*args.runs, args.crosscheck]:
        path = root / run / "front.csv"
        if path.exists():
            f = pd.read_csv(path)
            f["source"] = run
            fronts.append(f)
        else:
            print(f"no front for {run} ({path})")
    union = pd.concat(fronts, ignore_index=True).drop_duplicates(subset=["overrides"])
    union = union[union.T_mean >= TORQUE_FLOOR]
    F = np.c_[-union.T_mean.to_numpy(float), union.ripple_pct.to_numpy(float)]
    front = union[design.nondominated(F)].sort_values("T_mean").reset_index(drop=True)
    Ff = np.c_[-front.T_mean.to_numpy(float), front.ripple_pct.to_numpy(float)]
    k = design.knee(Ff)
    if len(front) <= args.max:
        pick = list(range(len(front)))
    else:                                             # half spread along the front, half the nearest to the coarse knee
        span = np.where(np.ptp(Ff, axis=0) > 0, np.ptp(Ff, axis=0), 1.0)
        near = np.argsort(np.linalg.norm((Ff - Ff[k]) / span, axis=1))[: args.max // 2]
        spread = np.linspace(0, len(front) - 1, args.max - args.max // 2).round().astype(int)
        pick = sorted(set(near.tolist()) | set(spread.tolist()))
    chosen = front.iloc[pick].copy()
    knees = {}
    for run in [*args.runs, args.crosscheck]:
        kj = root / run / "knee.json"
        if kj.exists():
            knees[run] = json.loads(kj.read_text(encoding="utf-8"))
    candidates = [json.loads(o) for o in chosen.overrides] + [REFERENCE] + [kn["overrides"] for kn in knees.values()]
    print(f"union front: {len(front)} designs (coarse); re-scoring {len(chosen)} of them, the reference and {len(knees)} knees on the fine mesh")
    with femm_queue.femm_pool("optimization: finalize (fine)", task_timeout_s=design.TASK_TIMEOUT_S) as pool:
        rows = design.Evaluator(pool, run="final", scoring="P", fine=True)(candidates, stage="fine")
    folder = design.out_dir("final")
    fine = pd.read_csv(folder / "evaluations.csv").drop_duplicates(subset=["overrides"], keep="last")
    fine_ok = fine[(fine.feasible == True) & fine.T_mean.notna() & (fine.T_mean >= TORQUE_FLOOR)]     # noqa: E712
    on_front = fine_ok[fine_ok.overrides.isin([json.dumps(json.loads(o)) for o in chosen.overrides])]
    G = np.c_[-on_front.T_mean.to_numpy(float), on_front.ripple_pct.to_numpy(float)]
    fine_front = on_front[design.nondominated(G)].sort_values("T_mean").reset_index(drop=True)
    Gf = np.c_[-fine_front.T_mean.to_numpy(float), fine_front.ripple_pct.to_numpy(float)]
    best = fine_front.iloc[design.knee(Gf)]
    overrides = json.loads(best.overrides)
    p = design.params(overrides)
    lookup = {r: fine[fine.overrides == json.dumps(o)] for r, o in [("reference", REFERENCE)] + [(run, kn["overrides"]) for run, kn in knees.items()]}
    table = []
    for label, sub in [("chosen (fine knee)", fine_front.iloc[[design.knee(Gf)]])] + list(lookup.items()):
        if len(sub):
            r = sub.iloc[-1]
            table.append({"design": label, "T_mean": r.T_mean, "ripple_pct": r.ripple_pct, "ripple_pp": r.ripple_pp, "gamma": r.gamma,
                          "pf": r.pf, "psi_pm_mWb": 1e3 * r.psi_pm, "L_d_mH": 1e3 * r.L_d, "L_q_mH": 1e3 * r.L_q, "amps": r.amps,
                          "bridge_mm": r.bridge_mm, "overrides": r.overrides})
    pd.DataFrame(table).to_csv(folder / "compare.csv", index=False)
    final = {"overrides": overrides, "params": {v: getattr(p, v) for v in (*design.ALL_VARS, "h_open", "r_shoe")},
             "bridge_mm": design.bridge_of(p), "fine": {k: float(best[k]) for k in ("T_mean", "ripple_pct", "ripple_pp", "gamma", "pf", "amps")},
             "fine_front_size": int(len(fine_front)), "hv_fine_front": hypervolume(Gf), "union_front_coarse": int(len(front)),
             "sources": sorted(set(front.source)), "knees": knees, "compare": table}
    (folder / "final.json").write_text(json.dumps(final, indent=2, default=float), encoding="utf-8")
    for row in table:
        print(f"  {row['design']:24s} T {row['T_mean']:.4f} N.m  ripple {row['ripple_pct']:.2f} %  gamma {row['gamma']:.2f}  pf {row['pf']:.3f}  "
              f"psi_pm {row['psi_pm_mWb']:.2f} mWb  I {row['amps']:.2f} A")
    dxf = folder / f"{args.name}.dxf"
    export.write(overrides, dxf)
    check = export.compare(p, dxf)
    print("DXF read back:", "all checks pass" if check["ok"] else [c for c in check["checks"] if not c["ok"]])
    export.preview(dxf, dxf.with_suffix(".png"), f"the chosen design: {json.dumps({k: round(v, 3) for k, v in overrides.items()})}")
    _figure(fronts, fine_front, best, lookup, folder / "front_fine.png")
    print(f"wrote {folder}")
    return 0 if check["ok"] else 1


def _figure(fronts, fine_front, best, lookup, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colours = ["#2f6fb5", "#d08a2e", "#4c9f70", "#8e5ea2", "#7f8c8d"]
    fig, ax = plt.subplots(figsize=(8, 5.2))
    for f, c in zip(fronts, colours):
        ax.plot(f.T_mean, f.ripple_pct, "o-", ms=3, lw=0.9, color=c, alpha=0.75, label=f"{f.source.iloc[0]} (coarse)")
    ax.plot(fine_front.T_mean, fine_front.ripple_pct, "s-", ms=4, lw=1.4, color="#222222", label="the union front, fine mesh")
    ax.scatter([best.T_mean], [best.ripple_pct], s=150, facecolors="none", edgecolors="#c0392b", lw=1.8, zorder=5, label="the chosen design (fine knee)")
    for label, sub in lookup.items():
        if len(sub):
            r = sub.iloc[-1]
            ax.scatter([r.T_mean], [r.ripple_pct], marker="x" if label == "reference" else "D", s=50, zorder=5, label=f"{label} (fine)")
    ax.set_xlabel("mean torque at the scoring point, fixed copper loss [N·m]")
    ax.set_ylabel("ripple, peak-to-peak [% of the mean]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, frameon=False)
    ax.set_title("the fronts on the coarse mesh, the union re-scored on the fine mesh, and the chosen design", fontsize=10, loc="left")
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
