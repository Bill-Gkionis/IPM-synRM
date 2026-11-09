"""The two figures of the project note MTPA, from run.py's results and report.py's summary.json.
  fig_angle.png        the best current angle against current: FEMM, and the formula with the 10 A constants and with
                       each current's own L_d and L_q
  fig_trajectory.png   the (i_d, i_q) plane: circles of equal current, lines of equal torque (FEMM), and the MTPA path,
                       where each torque line just touches a current circle
Colours: a data-visualization reference palette (slots 1 to 3, which its validator passes for all pairs) and its chart ink.

Run, from the code folder:  uv run python learn/mtpa/figures.py      (seconds; writes into out/learn/mtpa/)
"""
import json
import math
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / "dq_torque"))   # after this folder, so its own report.py wins
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import common as dq
import paths
import report

OUT = paths.path("out") / "learn" / "mtpa"
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
FEMM_C, F10_C, FOWN_C = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 10, "axes.edgecolor": AXIS,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURFACE,
                     "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.titlecolor": INK,
                     "axes.titlesize": 11, "legend.frameon": False, "legend.labelcolor": INK2})


def frame(ax):
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def note(ax, x, y, text, dx=6, dy=0, ha="left", va="center"):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, color=INK2, fontsize=9)


def fig_angle(S):
    rows = S["rows"]
    I = [float(k) for k in rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for key, col, name in (("FEMM angle [deg]", FEMM_C, "FEMM"),
                           ("formula, 10 A constants [deg]", F10_C, "formula, 10 A constants"),
                           ("formula, own constants [deg]", FOWN_C, "formula, each current's own L_d, L_q")):
        y = [rows[k][key] for k in rows]
        ax.plot(I, y, "-o", color=col, lw=1.8, ms=5.5, mec=SURFACE, mew=1.2, label=name)
        note(ax, I[-1], y[-1], name.split(",")[0] if key.startswith("FEMM") else
             ("10 A constants" if "10 A" in key else "own constants"))
    ax.axhline(45, color=AXIS, lw=0.8)
    note(ax, 1, 45, "45°: a pure synRM", dx=0, dy=5, va="bottom")
    frame(ax)
    ax.set_xlim(0, 49); ax.set_xticks(range(0, 41, 5)); ax.set_ylim(0, 50)
    ax.set_xlabel("current [A peak]")
    ax.set_ylabel("best current angle γ, from +q toward −d [°]")
    ax.set_title("Where MTPA puts the current, against current", loc="left")
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_angle.png", dpi=200)
    plt.close(fig)


def fig_trajectory(df, S):
    pts = df[(df.I > 0) & (df.gamma <= 60)].groupby(["I", "gamma"])[["id", "iq", "T_age"]].mean().reset_index()
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    th = np.radians(np.linspace(0, 90, 200))
    for I in sorted(pts.I.unique()):                                   # circles of equal current
        ax.plot(-I * np.sin(th), I * np.cos(th), color=GRID, lw=0.9, zorder=1)
        note(ax, 0, I, f"{I:g} A", dx=3, dy=0)
    xs = np.concatenate([pts.id.values, [0.0]]); ys = np.concatenate([pts.iq.values, [0.0]])
    zs = np.concatenate([pts.T_age.values, [0.0]])
    levels = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
    cs = ax.tricontour(xs, ys, zs, levels=[v for v in levels if v < zs.max()], colors=MUTED, linewidths=0.8)
    ax.clabel(cs, fmt=lambda v: f"{v:g} N·m", fontsize=8, colors=INK2)
    P = pd.DataFrame(S["points"])
    i_d = [-r.I * math.sin(math.radians(r.gamma)) for r in P.itertuples()]
    i_q = [r.I * math.cos(math.radians(r.gamma)) for r in P.itertuples()]
    ax.plot([0] + i_d, [0] + i_q, "-o", color=FEMM_C, lw=1.8, ms=5.5, mec=SURFACE, mew=1.2, zorder=4,
            label="MTPA, FEMM")
    Ig = np.linspace(0.5, 40, 80)
    g10 = [report.formula_angle(S["psi_pm"], S["dL10"], I) for I in Ig]
    ax.plot(-Ig * np.sin(np.radians(g10)), Ig * np.cos(np.radians(g10)), color=F10_C, lw=1.8, zorder=3,
            label="MTPA, formula with the 10 A constants")
    frame(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-30, 2); ax.set_ylim(0, 42)
    ax.set_xlabel("i_d [A]"); ax.set_ylabel("i_q [A]")
    ax.set_title("The MTPA path: where torque lines touch current circles", loc="left")
    ax.legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_trajectory.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    df = pd.read_csv(OUT / "solves.csv")
    with open(OUT / "summary.json", encoding="utf-8") as f:
        S = json.load(f)
    fig_angle(S)
    fig_trajectory(df, S)
    print(f"wrote fig_angle.png and fig_trajectory.png into {OUT}")
