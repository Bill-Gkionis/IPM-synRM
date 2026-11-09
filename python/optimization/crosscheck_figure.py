"""optimization.crosscheck_figure: the picture of the crosscheck run: every admissible design it scored, its front, its knee, and how the front grew.

What it shows
    Left, the plane of the two objectives (mean torque, ripple in per cent of the mean): every admissible design the run scored as a
    small grey dot (the context), the non-dominated front as a blue line with its designs on it, and the knee as one orange dot with its numbers.
    The corner of the axes is the reference point of the hypervolume (0.8 N.m, 40 %): the area the front dominates up to it is the number
    the run reports. Right, that hypervolume against the number of designs that reached FEMM, the run's own convergence record
    (progress.csv), in bands: the first block of space-filling designs, the rounds of the surrogate search, the probes (the surrogate's own optima and
    the ripple map), and the rounds that followed. The corner of the left panel's lower right can hold no design: one there would beat the front.
    The colours are the first two slots of the reference palette (blue, orange) on the light surface, checked with the palette validator
    (all pairs, light mode: every check passes); the grey is the neutral ink for context, never a series colour. Marks follow the mark
    specs: lines 2 px, the front's and knee's dots 8 and 12 px with a 2 px ring in the surface colour, the context dots small and
    translucent (there are a thousand of them), hairline solid grid, text in the text inks. front.csv is the table twin of the picture.

Run (from the Python to FEMM folder; reads evaluations.csv, front.csv, knee.json and progress.csv of out\\optimization\\crosscheck\\)
    .venv/Scripts/python.exe -m optimization.crosscheck_figure [--folder <another folder>]
Writes figure.png beside them.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from optimization import crosscheck as cc

DPI = 200
PX = 72.0 / DPI                                     # one pixel in points at this resolution
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"                 # categorical slots 1 and 2 (light mode)


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(1 * PX)
    ax.tick_params(colors=MUTED, labelcolor=INK2, length=3, width=1 * PX, labelsize=10)
    ax.grid(True, color=GRID, linewidth=1 * PX, linestyle="-")
    ax.set_axisbelow(True)


def draw(folder: Path) -> Path:
    adm = cc.read_scored(folder)
    front = pd.read_csv(folder / "front.csv").sort_values("T_mean")
    knee = json.loads((folder / "knee.json").read_text(encoding="utf-8"))
    prog = pd.read_csv(folder / "progress.csv")

    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=DPI, gridspec_kw={"width_ratios": [1.55, 1.0], "wspace": 0.22})
    fig.patch.set_facecolor(SURFACE)

    # ---- left: the objective plane
    style(ax)
    t, r = adm["T_mean"].to_numpy(), adm["ripple_pct"].to_numpy()
    shown = r <= cc.REF_POINT[1]
    ax.plot(t[shown], r[shown], linestyle="none", marker="o", markersize=5 * PX, markerfacecolor=MUTED, markeredgewidth=0,
            alpha=0.40, label=f"admissible designs ({len(adm)})", zorder=2)
    ax.plot(front["T_mean"], front["ripple_pct"], color=BLUE, linewidth=2 * PX, solid_capstyle="round", solid_joinstyle="round", zorder=3)
    ax.plot(front["T_mean"], front["ripple_pct"], linestyle="none", marker="o", markersize=8 * PX, markerfacecolor=BLUE,
            markeredgecolor=SURFACE, markeredgewidth=2 * PX, label=f"non-dominated front ({len(front)})", zorder=4)
    ax.plot([knee["T_mean"]], [knee["ripple_pct"]], linestyle="none", marker="o", markersize=12 * PX, markerfacecolor=ORANGE,
            markeredgecolor=SURFACE, markeredgewidth=2 * PX, label="knee", zorder=5)
    ax.set_xlim(cc.T_FLOOR - 0.03, float(t.max()) + 0.04 * (float(t.max()) - cc.T_FLOOR))
    ax.set_ylim(0, cc.REF_POINT[1])
    ax.set_xlabel("mean torque [N.m]", color=INK2, fontsize=11)
    ax.set_ylabel("torque ripple [% of the mean torque]", color=INK2, fontsize=11)
    # direct label, sparingly: the knee gets its name beside the dot, below and to the right where nothing can lie (a design there would dominate the
    # front); the ink text is not the series colour, the orange dot carries the identity. The numbers of the knee and of the two ends of the front
    # go in a small block in the empty band above the cloud (the legend takes its left corner), so no label can run off the axes or over the front.
    ax.annotate("knee", xy=(knee["T_mean"], knee["ripple_pct"]), xytext=(16, -16), textcoords="offset points", color=INK, fontsize=10, ha="left", va="top",
                arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 1 * PX, "shrinkA": 0, "shrinkB": 6})
    lo_r = front.loc[front["ripple_pct"].idxmin()]
    hi_t = front.loc[front["T_mean"].idxmax()]
    ax.text(0.985, 0.975,
            f"knee: {knee['T_mean']:.3f} N.m, {knee['ripple_pct']:.2f} %\n"
            f"(MTPA angle {knee['gamma_deg']:.1f} deg, power factor {knee['pf']:.3f})\n"
            f"lowest ripple: {lo_r['ripple_pct']:.2f} % at {lo_r['T_mean']:.3f} N.m\n"
            f"highest torque: {hi_t['T_mean']:.3f} N.m at {hi_t['ripple_pct']:.2f} %",
            transform=ax.transAxes, ha="right", va="top", color=INK2, fontsize=9.5, linespacing=1.5)
    leg = ax.legend(loc="upper left", frameon=False, fontsize=10, handletextpad=0.5, borderaxespad=0.2)
    for text in leg.get_texts():
        text.set_color(INK2)
    ax.set_title(f"Every admissible design, the front and its knee  (hypervolume {knee['hypervolume']:.2f} N.m·% to the corner)",
                 loc="left", color=INK, fontsize=11, pad=10)

    # ---- right: the front's growth
    style(bx)
    bx.plot(prog["reached"], prog["hv"], color=BLUE, linewidth=2 * PX, solid_capstyle="round", solid_joinstyle="round", zorder=3)
    bx.plot(prog["reached"].iloc[-1:], prog["hv"].iloc[-1:], linestyle="none", marker="o", markersize=8 * PX, markerfacecolor=BLUE,
            markeredgecolor=SURFACE, markeredgewidth=2 * PX, zorder=4)
    group = prog["stage"].map(lambda st: "design" if st == "design" else ("probes" if str(st).startswith("probe") else "rounds"))
    start = 0.0
    labels = {"design": "design", "rounds": "rounds", "probes": "probes"}
    for g, rows in prog.groupby((group != group.shift()).cumsum(), sort=False):
        name, end = group[rows.index[0]], float(rows["reached"].iloc[-1])
        bx.axvline(end, color=AXIS, linewidth=1 * PX, zorder=1)
        bx.text((start + end) / 2, 0.03, labels[name], color=MUTED, fontsize=9, va="bottom", ha="center", transform=bx.get_xaxis_transform())
        start = end
    bx.annotate(f"{prog['hv'].iloc[-1]:.2f}", xy=(prog["reached"].iloc[-1], prog["hv"].iloc[-1]), xytext=(0, -18), textcoords="offset points",
                color=INK, fontsize=10, ha="center")
    bx.set_xlabel("designs that reached FEMM", color=INK2, fontsize=11)
    bx.set_ylabel("hypervolume to (0.8 N.m, 40 %)  [N.m·%]", color=INK2, fontsize=11)
    bx.set_xlim(0, max(cc.BUDGET, float(prog["reached"].max())) * 1.02)
    bx.set_title("How the front grew", loc="left", color=INK, fontsize=11, pad=10)

    out = folder / "figure.png"
    fig.savefig(out, dpi=DPI, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m optimization.crosscheck_figure", description=__doc__.splitlines()[0])
    parser.add_argument("--folder", type=Path, default=None, help="the run's folder (default out\\optimization\\crosscheck)")
    args = parser.parse_args(argv)
    out = draw(args.folder or cc.out_folder())
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
