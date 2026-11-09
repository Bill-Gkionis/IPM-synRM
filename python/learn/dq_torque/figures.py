"""The three figures of the project note Torque in d-q - field, current and saliency.

It solves nothing: it reads run.py's results (solves.csv, probes.csv, fields/*.ans) and report.py's summary.json.
  fig_torque.png    torque against current angle at 10 A: magnets on and off; and the torque split into its three terms
  fig_bridges.png   |B| in one rotor pole, d axis up: magnets on with no current, on with 10 A on -d, off with 10 A on -d
  fig_saliency.png  the saliency L_q/L_d against current, magnets on and off
Colours: a data-visualization reference palette (categorical slots checked with its validator: all checks pass; the light
slots below 3:1 contrast are carried by direct labels and the note's tables) and its one-hue blue ramp for |B|.

Run, from the code folder:  uv run python learn/dq_torque/figures.py      (a few seconds)
Writes into out/learn/dq_torque/. The notes' copies are made by hand: no script writes
into the project notes.
"""
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

import common as c
import femfile
import motorgeo as mg
import report

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
ON, OFF = "#2a78d6", "#eb6834"                             # slots 1 and 2: magnets on, magnets off
MAG, QX, DX = "#1baf7a", "#eda100", "#e87ba4"              # slots 3 to 5: the three torque terms
BLUE_RAMP = LinearSegmentedColormap.from_list("blue", [SURFACE, "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                                                       "#256abf", "#184f95", "#104281", "#0d366b"])
plt.rcParams.update({"font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 10, "axes.edgecolor": AXIS,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURFACE,
                     "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.titlecolor": INK,
                     "axes.titlesize": 11, "legend.frameon": False, "legend.labelcolor": INK2})


def frame(ax):
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def label(ax, x, y, text, dx=4, dy=0, ha="left", va="center"):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, color=INK2, fontsize=9)


def fig_torque(df, S):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.3), sharey=True)
    for v, col, name in (("on", ON, "magnets on"), ("off", OFF, "magnets off")):
        gm = report.by_gamma(df, v)
        a.plot(gm.index, gm.T_age, color=col, lw=1.8, solid_capstyle="round", label=name)
    g, T = S["on"]["mtpa_gamma"], S["on"]["T_mtpa"]
    a.plot([g], [T], "o", ms=7, color=ON, mec=SURFACE, mew=1.5, zorder=5)
    label(a, g, T, f"MTPA: {T:.3f} N·m at {g:.1f}°", dx=8, dy=6, va="bottom")
    label(a, 20, 0.0, "magnets off: 0.0006 N·m at best", dx=0, dy=6, va="bottom")
    a.set_title("Torque against current angle, 10 A peak", loc="left")
    a.set_xlabel("current angle γ from +q toward −d [°]")
    a.set_ylabel("torque [N·m]  (FEMM depth 20 mm)")
    a.legend(loc="lower left")

    gm = report.by_gamma(df, "on")
    pm = S["on"]["psi_pm"]
    parts = pd.DataFrame([c.three_terms(report.PP, pm, r.psi_d, r.psi_q, r.id, r.iq) for r in gm.itertuples()],
                         index=gm.index)
    b.axhline(0, color=AXIS, lw=0.8)
    b.plot(gm.index, gm.T_age, color=ON, lw=1.8, label="total (FEMM)")
    b.plot(parts.index, parts["magnets x i_q"], color=MAG, lw=1.8, label="magnets × i_q")
    b.plot(parts.index, parts["q field x i_d"], color=QX, lw=1.8, label="iron's q field × i_d")
    b.plot(parts.index, parts["d field of i_d x i_q"], color=DX, lw=1.8, label="i_d's own d field × i_q")
    for key, gx, text in (("magnets x i_q", 20, "magnets × i_q"), ("q field x i_d", 33, "q field × i_d"),
                          ("d field of i_d x i_q", 45, "i_d's d field × i_q")):
        y = float(np.interp(gx, parts.index, parts[key]))
        label(b, gx, y, text, dx=0, dy=-7, ha="center", va="top")          # each just under its own line
    label(b, 8, float(np.interp(8, gm.index, gm.T_age)), "total", dx=0, dy=7, va="bottom")
    b.set_title("The same torque, as three field × current products", loc="left")
    b.set_xlabel("current angle γ from +q toward −d [°]")
    b.legend(loc="upper right", fontsize=9)
    for ax in (a, b):
        frame(ax)
        ax.set_xlim(0, 90)
        ax.set_xticks(range(0, 91, 15))
    a.set_ylim(-0.5, 1.3)
    fig.tight_layout()
    fig.savefig(c.OUT / "fig_torque.png", dpi=200)
    plt.close(fig)


def rotated(x, y, deg):
    t = math.radians(deg)
    return x * math.cos(t) - y * math.sin(t), x * math.sin(t) + y * math.cos(t)


def fig_bridges(S):
    p = c.params()
    regions = mg.build(p)
    turn = 90 - 90 / p.pole_pairs                     # rotate pole 0's d axis (22.5 deg as drawn) to point straight up
    win = (-19.5, 19.5, 27.5, 47.0)
    probes = pd.DataFrame(S["probes"])
    pts = c.probe_points(p, regions)
    panels = [("on", "no load", "Magnets on, no current"),
              ("on", "10 A on -d", "Magnets on, 10 A on −d"),
              ("off", "10 A on -d", "Magnets off, 10 A on −d")]
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.5))
    for ax, (v, state, title) in zip(axs, panels):
        f = femfile.read_ans(str(c.FIELDS / f"{v}_{state.replace(' ', '_').replace('-', 'minus')}.ans"))
        bx, by = femfile.element_b(f["xy"], f["A"], f["tri"])
        x, y = rotated(f["xy"][:, 0], f["xy"][:, 1], turn)
        cx, cy = x[f["tri"]].mean(axis=1), y[f["tri"]].mean(axis=1)
        keep = (cx > win[0] - 2) & (cx < win[1] + 2) & (cy > win[2] - 2) & (cy < win[3] + 2)
        tri = mtri.Triangulation(x, y, f["tri"], mask=~keep)
        pc = ax.tripcolor(tri, facecolors=np.hypot(bx, by), cmap=BLUE_RAMP, vmin=0, vmax=2.3, shading="flat",
                          rasterized=True)
        for r in regions:                              # the magnet slots: outlines only
            if r.material == "N45":
                ox, oy = r.polygon().exterior.xy
                ax.plot(*rotated(np.array(ox), np.array(oy), turn), color=INK2, lw=0.7)
        th = np.radians(np.linspace(0, 180, 400))
        for rad in (p.r_ro, p.r_si):
            ax.plot(rad * np.cos(th), rad * np.sin(th), color=MUTED, lw=0.6)
        for deg, name in ((90, "d"), (90 - 90 / p.pole_pairs, "q"), (90 + 90 / p.pole_pairs, "q")):
            t = math.radians(deg)                      # the axis line, labelled at its inner end, inside the panel
            ax.plot([27.0 * math.cos(t), 46.6 * math.cos(t)], [27.0 * math.sin(t), 46.6 * math.sin(t)], color=INK, lw=0.8)
            ax.annotate(name, (32.0 * math.cos(t), 32.0 * math.sin(t)), xytext=(5, 0), textcoords="offset points",
                        color=INK, fontsize=11, ha="left", va="center", fontweight="bold")
        mine = probes[(probes.variant == v) & (probes.state == state)].set_index("where").B_T
        for where, (px, py) in pts.items():
            rx, ry = rotated(px, py, turn)
            ax.plot([rx], [ry], "o", ms=5, color=INK, mec=SURFACE, mew=1.2, zorder=5)
            side = {"outer bridge, left": (-6, 0, "right"), "outer bridge, right": (6, 0, "left"),
                    "centre rib": (0, -9, "center")}[where]
            ax.annotate(f"{mine[where]:.2f} T", (rx, ry), xytext=side[:2], textcoords="offset points", ha=side[2],
                        va="center", color=INK, fontsize=9,
                        bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.85))
        ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3])
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(title, loc="left")
    cb = fig.colorbar(pc, ax=axs, shrink=0.65, pad=0.02, aspect=14)
    cb.set_label("|B| [T]", color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=MUTED)
    fig.savefig(c.OUT / "fig_bridges.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_saliency(S):
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    ax.axhline(1.0, color=AXIS, lw=0.8)
    label(ax, 0.5, 1.0, "1.0: no saliency (a round rotor)", dx=0, dy=-6, va="top")
    for v, col, name in (("on", ON, "magnets on"), ("off", OFF, "magnets off")):
        xs = [float(k) for k in S["saliency"][v]]
        ys = list(S["saliency"][v].values())
        ax.plot(xs, ys, "-o", color=col, lw=1.8, ms=6, mec=SURFACE, mew=1.2, label=name)
        label(ax, xs[-1], ys[-1], name, dx=8)
    frame(ax)
    ax.set_xlim(0, 52); ax.set_xticks([0, 5, 10, 20, 30, 40])
    ax.set_ylim(0.9, 2.0)
    ax.set_xlabel("current [A peak]  (10 turns per coil side)")
    ax.set_ylabel("saliency  L_q / L_d")
    ax.set_title("Saliency against current", loc="left")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(c.OUT / "fig_saliency.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    df = pd.read_csv(c.OUT / "solves.csv")
    with open(c.OUT / "summary.json", encoding="utf-8") as f:
        S = json.load(f)
    fig_torque(df, S)
    fig_bridges(S)
    fig_saliency(S)
    print(f"wrote fig_torque.png, fig_bridges.png and fig_saliency.png into {c.OUT}")
