"""Figure for the optimizer demo: every rotor it tried, the Pareto front, the reference design, the knee, and the knee's shape.

Run: uv run python optimize_plot.py <output.png>   (from this folder; reads out/opt/evaluations.csv and pareto_X.npy)
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import motorgeo as mg

OUT = os.path.abspath("out/opt")
NAMES = ["mag_y", "mag_x", "v_half_deg", "web"]
INK, INK2, MUTED, GRID, AXIS, SURF, BLUE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#2a78d6"

log = pd.read_csv(os.path.join(OUT, "evaluations.csv"))
ok = log[log.feasible == True]                                   # noqa: E712 (pandas reads True/False)
X = np.load(os.path.join(OUT, "pareto_X.npy"))
par = pd.concat([ok[np.isclose(ok[NAMES].to_numpy(), x, atol=1e-9).all(axis=1)].head(1) for x in X])
F = np.column_stack([-par.T_mean, par.ripple])
Fn = (F - F.min(axis=0)) / np.where(np.ptp(F, axis=0) > 0, np.ptp(F, axis=0), 1)
knee = par.iloc[int(np.argmin(np.hypot(*Fn.T)))]
base = ok.iloc[0]                                                # generation 0, row 0 = the reference design (seeded)

fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.25, 1]}, facecolor=SURF)
ax = axs[0]
ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ax.spines.values():
    s.set_color(AXIS)
ax.tick_params(colors=MUTED)
ax.plot(ok.ripple, ok.T_mean, "o", ms=5, color=AXIS, label=f"every buildable rotor tried ({len(ok)})")
p = par.sort_values("ripple")
ax.plot(p.ripple, p.T_mean, "-o", ms=6, color=BLUE, lw=1.8, label=f"Pareto front ({len(par)})")
ax.plot(base.ripple, base.T_mean, "*", ms=15, color=INK, label="reference design (seeded into generation 0)")
ax.plot(knee.ripple, knee.T_mean, "o", ms=14, mfc="none", mec=BLUE, mew=2, label="knee: the design sent to CAD")
ax.set_xlabel("torque ripple, peak-to-peak [N·m]  (lower is better)", color=INK2)
ax.set_ylabel("mean torque at 10 A, own MTPA angle [N·m]  (higher is better)", color=INK2)
n_bad = int((log.feasible == False).sum())                     # noqa: E712
ax.set_title(f"{len(log)} rotors in {log.gen.max() + 1} generations ({n_bad} rejected before FEMM). "
             "Coarse mesh, demo only", loc="left", fontsize=10, color=INK)
ax.legend(frameon=False, fontsize=8, loc="lower right")

ax = axs[1]
ax.set_facecolor(SURF)
for row, style, lab in ((base, dict(color=INK, lw=1.0), "reference design"), (knee, dict(color=BLUE, lw=1.8), "knee design")):
    prm = mg.Params(**{n: float(row[n]) for n in NAMES})
    for r in mg.build(prm):
        if r.name in ("magnet_0_up", "magnet_0_dn", "barrier_0_up", "barrier_0_dn"):
            g = r.polygon()
            xy = np.array([mg.rot(q, -22.5) for q in g.exterior.coords])
            ax.plot(xy[:, 0], xy[:, 1], **style, label=lab if r.name == "magnet_0_up" else None)
t = np.radians(np.linspace(-25, 25, 200))
ax.plot(45 * np.cos(t), 45 * np.sin(t), color=AXIS, lw=0.8)
ax.axhline(0, color=MUTED, lw=0.7, ls="-.")
ax.set_xlim(32, 46); ax.set_ylim(-17, 17); ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2)
fmt = lambda row: ", ".join(f"{n} {row[n]:.2f}" for n in NAMES)
ax.set_title(f"One pole, d-axis horizontal\nreference: {fmt(base)}\nknee:      {fmt(knee)}", loc="left", fontsize=8.5,
             color=INK)
fig.tight_layout()
fig.savefig(sys.argv[1], dpi=150, facecolor=SURF)
print(f"wrote {sys.argv[1]} | base T {base.T_mean:.4f} ripple {base.ripple:.4f} | knee T {knee.T_mean:.4f} "
      f"ripple {knee.ripple:.4f} | {fmt(knee)}")
