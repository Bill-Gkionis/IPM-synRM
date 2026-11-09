"""optimization.thermal_figures: the four pictures of the thermal study (matplotlib PNGs drawn from what thermal.py computes; no solver).

What it is, and why
    A table of 24 sustainable-loss rows, eleven node temperatures in three coolings and two dozen heating curves is hard to read; these four pictures are what is
    looked at when choosing the cooling. thermal_study.main() builds the rows and calls make_all(rows); every picture has its numbers in a CSV beside it
    (sustainable.csv, time_constants.csv, overload.csv): the table view the data-visualization method asks for.

The pictures (each has a title, a note under it saying what is drawn, and a legend)
    fig_thermal_network.png         the network as a diagram: the eleven parts, every resistance in K/W for the 60 mm stack at the rated speed, where the
                                    losses go in, and the three outsides. Blue: stator iron; orange: winding; aqua: rotor; gray: housing and air.
    fig_thermal_node_temperatures.png  the temperature of every node at the rated point of each cooling (60 mm stack, class F, 3,000 rpm): the continuous
                                    current at which the first limit is reached. Dashed lines: the magnets' 80 C and class F's 155 C; the binding node is ringed.
    fig_thermal_sustainable_current.png  the continuous current density (left) and copper loss (right) against stack length for the three coolings. Solid:
                                    what the network allows (the smaller of the two limits). Dotted: what the winding alone would allow (class F). The gap
                                    between them is the cost of the magnets' 80 C.
    fig_thermal_transients.png      top: heating from cold at the continuous current, the winding's hot spot and the magnets (60 mm stack); bottom: the
                                    winding's hot spot after a step to 1.5x and 2x the continuous current from the hot state, and the magnets.

The look (the palette of optimization/winding_figures.py)
    Surface #fcfcfb, ink #0b0b0b and #52514e, hairline grid #e1e0d9, axes #c3c2b7; series blue #2a78d6, orange #eb6834, aqua #1baf7a; context gray #898781.
    The three coolings are blue (natural), orange (forced air), aqua (liquid) in every picture; the hot spot is orange and the magnets blue where the coolings
    are drawn apart. Aqua is 2.7:1 on the surface, so every aqua series also has a second channel (the marker shape) and is named in the legend and the CSV.

Run
    Not on its own: thermal_study.py calls make_all(rows) (uv run python -m optimization.thermal_study). Or: uv run python -m optimization.thermal_figures
Numbers to expect: four PNGs, 1430 px wide, drawn in a second or two.

Self-check (the answers follow).
  1. Why are the dotted and the solid lines of the third picture the same where the winding binds?
  2. Why is the heating curve of the magnets so much slower than the winding's in the liquid-cooled panel?
Answers.
  1. The solid line is the smaller of the winding's and the magnets' limits; where the winding is the smaller one, it is the dotted line.
  2. The winding's heat has a short way to the coolant (a thin liner, 36 teeth, the housing) and little mass; the rotor has a lot of mass (2 to 5 kg of steel) and
     can only lose heat through the air gap, the shaft and the end space, so it follows the stator with a lag of tens of minutes.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

import paths
from optimization import thermal
from optimization import winding_figures as wf
from optimization.thermal import CLASS_F_C, CLASS_H_C, MAGNET_LIMIT_C, Cooling, ThermalModel

SURFACE, INK, INK2, MUTED, GRID, AXIS = wf.SURFACE, wf.INK, wf.INK2, wf.MUTED, wf.GRID, wf.AXIS
BLUE, ORANGE, AQUA = wf.BLUE, wf.ORANGE, wf.AQUA
COOLINGS = (("natural", Cooling.natural(), BLUE, "o"), ("forced air", Cooling.forced_air(), ORANGE, "s"), ("liquid", Cooling.liquid(), AQUA, "^"))
OUT = lambda: paths.path("out") / "optimization" / "thermal"
STACK = 60.0


def _ax(ax, xlabel: str, ylabel: str) -> None:
    wf._style(ax, xlabel, ylabel)


# ---------------------------------------------------------------------------------------------------- 1. the network
LAYOUT = {"E": (3.0, 8.0), "A": (6.4, 6.0), "B": (10.2, 6.0), "W": (0.9, 6.0), "T": (0.9, 3.7), "Y": (3.6, 3.7), "H": (7.7, 3.7),
          "R": (3.4, 1.2), "M": (0.9, 1.2), "S": (6.6, 1.2), "G": (9.8, 1.2), "ambient": (12.6, 3.9), "coolant": (12.6, 2.4)}
FAMILY = {"W": ORANGE, "E": ORANGE, "T": BLUE, "Y": BLUE, "R": AQUA, "M": AQUA, "S": AQUA, "G": AQUA, "H": MUTED, "B": MUTED, "A": "#c3c2b7"}


def fig_network(path: Path) -> Path:
    m = ThermalModel(stack_mm=STACK, cooling=Cooling.natural())
    edges = m.edges(dt_outside=40.0)
    edges.append(("H", "coolant", dict(((a, b), r) for a, b, r in ThermalModel(stack_mm=STACK, cooling=Cooling.liquid()).edges(dt_outside=40.0))[("H", "coolant")]))
    ir = m.iron_placeholder()
    fig = wf._figure(14.3, 8.4, "The thermal network: eleven parts, the resistance of every path in K/W (60 mm stack, 3,000 rpm)",
                     "Each box is the mean temperature of a part; each line a thermal resistance in K/W (lower = easier for the heat). Orange: the winding (slot copper W, end turns E); "
                     "blue: stator iron (teeth T, yoke Y); aqua: rotor (steel R, magnets M, shaft S, bearings G); gray: housing H, end bells B, end-space air A. Losses go in at the nodes "
                     "named in the arrows (copper at W and E, iron at T, Y and R, magnet eddy loss at M, bearings at G, windage at T and R). The outside is natural convection here; the "
                     "other two coolings replace the last box.")
    ax = fig.add_axes([0.01, 0.12, 0.98, 0.82])
    ax.set_xlim(-0.8, 14.2)
    ax.set_ylim(0.0, 9.0)
    ax.axis("off")
    ax.set_facecolor(SURFACE)
    names = thermal.NODE_NAMES
    pos = dict(LAYOUT)
    for a, b, r in edges:
        if (a, b) == ("E", "B"):
            continue
        if b in ("ambient", "coolant"):
            continue
    drawn = []
    for a, b, r in edges:
        if b in ("ambient", "coolant"):
            xa, ya = pos[a]
            xb, yb = pos["ambient"] if b == "ambient" else pos["coolant"]
        else:
            (xa, ya), (xb, yb) = pos[a], pos[b]
        style = dict(color=INK2, linewidth=1.3, solid_capstyle="round", zorder=1)
        if b == "coolant":
            style = dict(color=MUTED, linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)
        if (a, b) == ("E", "B"):
            ax.plot([xa, 6.6, xb], [ya, 8.8, yb + 0.42], color=MUTED, linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)
            ax.text(6.6, 8.8, f"{r:.1f}", fontsize=8, color=INK2, ha="center", va="center", zorder=3, bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec=AXIS, lw=0.6))
            continue
        ax.plot([xa, xb], [ya, yb], **style)
        frac = {("H", "ambient"): 0.78, ("G", "B"): 0.3, ("A", "S"): 0.55, ("A", "R"): 0.3, ("B", "ambient"): 0.5, ("H", "coolant"): 0.78}.get((a, b), 0.5)
        xm, ym = xa + frac * (xb - xa), ya + frac * (yb - ya)
        label = f"{r:.2f}" if r < 10 else f"{r:.0f}"
        ax.text(xm, ym, label, fontsize=8.5, color=INK, ha="center", va="center", zorder=3, bbox=dict(boxstyle="round,pad=0.18", fc=SURFACE, ec=AXIS, lw=0.6))
        drawn.append((a, b))
    for key, (x, y) in pos.items():
        if key in ("ambient", "coolant"):
            continue
        ax.add_patch(FancyBboxPatch((x - 0.85, y - 0.42), 1.7, 0.84, boxstyle="round,pad=0.02,rounding_size=0.18", fc=FAMILY[key], ec=SURFACE, lw=2.0, zorder=2))
        ax.text(x, y + 0.08, key, ha="center", va="center", fontsize=13, fontweight="bold", color="#ffffff" if key not in ("A",) else INK, zorder=4)
        ax.text(x, y - 0.24, names[key], ha="center", va="center", fontsize=6.8, color="#ffffff" if key not in ("A",) else INK, zorder=4)
    for key, label in (("ambient", "ambient 40 C"), ("coolant", "coolant 40 C (jacket)")):
        x, y = pos[key]
        ax.plot([x, x], [y + 0.25, y - 0.25], color=INK, lw=2)
        for k in range(3):
            ax.plot([x - 0.35 + 0.1 * k, x + 0.35 - 0.1 * k], [y - 0.25 - 0.12 * k, y - 0.25 - 0.12 * k], color=INK, lw=2)
        ax.text(x, y + 0.5, label, ha="center", fontsize=8, color=INK2)
    # the three outsides
    rows = []
    for cname, c, col, mk in COOLINGS:
        mm = ThermalModel(stack_mm=STACK, cooling=c)
        e = {(a, b): r for a, b, r in mm.edges(dt_outside=40.0)}
        rows.append((cname, e.get(("H", "ambient")), e.get(("H", "coolant")), e[("B", "ambient")]))
    ax.text(12.6, 8.55, "the outside, by cooling [K/W]", ha="center", fontsize=8.5, color=INK, fontweight="bold")
    ax.text(12.6, 8.2, "H to air or coolant  |  B to air", ha="center", fontsize=7.5, color=INK2)
    for k, (cname, rh, rc, rb) in enumerate(rows):
        main = rh if rh is not None else rc
        ax.text(12.6, 7.8 - 0.42 * k, f"{cname}:  {main:.3g}  |  {rb:.3g}", ha="center", fontsize=8.2, color=INK2)
    # the sources
    cu = ir
    sources = [("W", "copper in the slots"), ("E", "copper in the end turns"), ("T", "iron + half the windage"), ("Y", "iron"), ("R", "rotor iron + half the windage"),
               ("M", "magnet eddy loss (input, 0)"), ("G", "bearings")]
    for key, text in sources:
        x, y = pos[key]
        dx, dy = (-0.6, 0.0)
        ax.annotate("", xy=(x - 0.82, y + 0.18), xytext=(x - 1.5, y + 0.18), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.6)) if key in ("W", "T", "M") else None
        if key in ("E",):
            ax.annotate("", xy=(x, y + 0.45), xytext=(x, y + 0.95), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.6))
            ax.text(x, y + 1.05, "P copper (ends)", ha="center", fontsize=7.5, color=ORANGE)
        if key in ("Y", "R", "G"):
            ax.annotate("", xy=(x, y - 0.45), xytext=(x, y - 0.95), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.6))
            ax.text(x, y - 1.05, {"Y": "P iron", "R": "P iron, windage", "G": "P bearings"}[key], ha="center", fontsize=7.5, color=ORANGE, va="top")
    ax.text(-0.75, 6.75, "P copper (slots)", fontsize=7.5, color=ORANGE)
    ax.text(-0.75, 4.45, "P iron, windage", fontsize=7.5, color=ORANGE)
    ax.text(-0.75, 1.9, "P magnets", fontsize=7.5, color=ORANGE)
    ax.text(0.9, 7.3, "winding: liner + slot", fontsize=7, color=INK2, ha="left")
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------------------------------- 2. node temperatures at the rated point
def fig_node_temperatures(path: Path) -> Path:
    order = [("hot spot (winding)", "hot"), ("W", "W"), ("E", "E"), ("T", "T"), ("Y", "Y"), ("H", "H"), ("B", "B"), ("A", "A"), ("R", "R"), ("M", "M"), ("S", "S"), ("G", "G")]
    fig = wf._figure(14.3, 5.6, "Temperature of every part at the rated point: the continuous current at which the first limit is reached (60 mm, class F, 3,000 rpm, 40 C ambient)",
                     "One dot per part (mean temperature of the node); the winding's hot spot is the higher of the slot's and the end turn's. Dashed vertical lines: the magnets' 80 C "
                     "(a hard limit) and class F's 155 C; dotted: the 40 C ambient or coolant. The ringed dot is the binding limit.")
    gs = fig.add_gridspec(1, 3, left=0.09, right=0.985, top=0.80, bottom=0.17, wspace=0.12)
    for k, (cname, c, col, mk) in enumerate(COOLINGS):
        ax = fig.add_subplot(gs[0, k])
        m = ThermalModel(stack_mm=STACK, cooling=c)
        s = m.sustainable()
        T = dict(s["T"])
        vals = [s["hot_spot"] if key == "hot" else T[key] for _, key in order]
        ys = np.arange(len(order))[::-1]
        _ax(ax, "temperature [C]", "")
        ax.set_facecolor(SURFACE)
        ax.set_yticks(ys)
        ax.set_yticklabels([lab if key == "hot" else f"{key}  {thermal.NODE_NAMES[key]}" for lab, key in order], fontsize=8, color=INK2)
        if k > 0:
            ax.set_yticklabels([])
        ax.axvline(MAGNET_LIMIT_C, color=INK2, lw=1.2, linestyle=(0, (4, 3)))
        ax.axvline(CLASS_F_C, color=INK2, lw=1.2, linestyle=(0, (4, 3)))
        ax.axvline(c.ambient_c, color=MUTED, lw=1.0, linestyle=(0, (1, 2)))
        ax.hlines(ys, c.ambient_c, vals, color=GRID, lw=1.2, zorder=1)
        ax.scatter(vals, ys, s=64, color=col, edgecolor=SURFACE, linewidth=1.6, marker=mk, zorder=3)
        bind_key = "hot" if s["binding"] == "winding" else "M"
        yb = ys[[key for _, key in order].index(bind_key)]
        ax.scatter([vals[[key for _, key in order].index(bind_key)]], [yb], s=230, facecolor="none", edgecolor=INK, linewidth=1.6, zorder=4)
        ax.set_xlim(30, 165)
        ax.text(0.0, 1.10, f"{cname}", transform=ax.transAxes, fontsize=11, fontweight="bold", color=INK)
        ax.text(0.0, 1.03, f"J = {s['J_rms']:.2f} A/mm2, copper {s['P_cu']:.0f} W, binding: {s['binding']}", transform=ax.transAxes, fontsize=8.5, color=INK2)
        ax.text(MAGNET_LIMIT_C + 1.5, ys[-1] - 0.75, "80 C", fontsize=7.5, color=INK2)
        ax.text(CLASS_F_C - 1.5, ys[-1] - 0.75, "155 C", fontsize=7.5, color=INK2, ha="right")
        ax.set_ylim(-1.2, len(order) - 0.3)
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------------------------------- 3. sustainable current against stack length
def fig_sustainable(path: Path, rows: list[dict]) -> Path:
    fig = wf._figure(14.3, 5.4, "Sustainable current density and copper loss against stack length (class F, 3,000 rpm, copper fill 0.45)",
                     "Solid: what the network allows (the smaller of the winding's and the magnets' limits). Dotted: what the winding alone would allow (class F, 155 C). Where the two "
                     "differ the magnets bind (N45, 80 C). Iron loss: the provisional placeholder at 200 Hz; magnet eddy loss: 0. Class H gives the same solid line in every case shown.")
    gs = fig.add_gridspec(1, 2, left=0.065, right=0.985, top=0.84, bottom=0.17, wspace=0.17)
    axj, axp = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    _ax(axj, "stack length [mm] (rows of magnets x 20 mm)", "continuous current density [A/mm2 rms]")
    _ax(axp, "stack length [mm]", "copper loss at the limit [W]")
    for cname, c, col, mk in COOLINGS:
        rr = sorted([r for r in rows if r["cooling"] == cname and r["class"] == "F"], key=lambda r: r["stack_mm"])
        x = [r["stack_mm"] for r in rr]
        axj.plot(x, [r["J_rms_A_mm2"] for r in rr], color=col, lw=2.2, marker=mk, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.5, label=f"{cname}: allowed")
        axj.plot(x, [r["J_winding_limit"] for r in rr], color=col, lw=1.6, linestyle=(0, (1.5, 2.5)), marker=mk, markersize=7, markerfacecolor=SURFACE, markeredgecolor=col,
                 label=f"{cname}: winding alone")
        axp.plot(x, [r["P_cu_W"] for r in rr], color=col, lw=2.2, marker=mk, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.5)
        axj.text(x[-1] + 1.5, rr[-1]["J_rms_A_mm2"], cname, color=INK2, fontsize=9, va="center")
        axp.text(x[-1] + 1.5, rr[-1]["P_cu_W"], cname, color=INK2, fontsize=9, va="center")
    for ax in (axj, axp):
        ax.set_xticks([20, 40, 60, 80])
        ax.set_xlim(15, 96)
        ax.set_ylim(bottom=0)
    leg = axp.legend(handles=axj.get_legend_handles_labels()[0], labels=axj.get_legend_handles_labels()[1], loc="upper left", frameon=False, fontsize=8.2, ncol=1)
    for t in leg.get_texts():
        t.set_color(INK2)
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------------------------------- 4. transients
def fig_transients(path: Path, rows: list[dict]) -> Path:
    fig = wf._figure(14.3, 8.6, "Heating curves at the continuous current, and overloads from the hot state (60 mm, class F, 3,000 rpm)",
                     "Top: from cold at the continuous current of each cooling, the winding's hot spot (orange) and the magnets (blue) against time; dashed: the magnets' 80 C and class F's 155 C. "
                     "Bottom: from the hot state (the steady state at the continuous current) the hot spot after a step to 1.5x (orange, dashed) and 2x (orange) the current, and the magnets "
                     "at 2x (blue, dots); dashed lines: class F 155 C, class H 180 C, magnets 80 C. A curve that stays under a line does not reach that limit in the time shown.")
    gs = fig.add_gridspec(2, 3, left=0.065, right=0.985, top=0.90, bottom=0.115, wspace=0.18, hspace=0.34)
    xmax_top = {"natural": 240.0, "forced air": 120.0, "liquid": 80.0}
    xmax_bot = {"natural": 120.0, "forced air": 22.0, "liquid": 1.5}
    for k, (cname, c, col, mk) in enumerate(COOLINGS):
        m = ThermalModel(stack_mm=STACK, cooling=c)
        s = m.sustainable()
        j = s["J_rms"]
        ax = fig.add_subplot(gs[0, k])
        _ax(ax, "time [min]", "temperature [C]" if k == 0 else "")
        tr = m.transient(xmax_top[cname] * 60.0, j, start="cold", dt_max_s=20.0)
        ax.plot(tr.t / 60, tr.hot_spot, color=ORANGE, lw=2.2, label="hot spot")
        ax.plot(tr.t / 60, tr.magnet, color=BLUE, lw=2.2, label="magnets")
        ax.axhline(MAGNET_LIMIT_C, color=INK2, lw=1.1, linestyle=(0, (4, 3)))
        ax.axhline(CLASS_F_C, color=INK2, lw=1.1, linestyle=(0, (4, 3)))
        ax.text(xmax_top[cname] * 0.99, MAGNET_LIMIT_C + 2, "80 C", ha="right", fontsize=7.5, color=INK2)
        ax.text(xmax_top[cname] * 0.99, CLASS_F_C + 2, "155 C", ha="right", fontsize=7.5, color=INK2)
        ax.set_xlim(0, xmax_top[cname])
        ax.set_ylim(35, 170)
        ax.text(0.02, 1.02, f"{cname}: J = {j:.2f} A/mm2 from cold", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK)
        if k == 0:
            leg = ax.legend(loc="center right", frameon=False, fontsize=8.5)
            for t in leg.get_texts():
                t.set_color(INK2)
        ax = fig.add_subplot(gs[1, k])
        _ax(ax, "time [min]", "temperature [C]" if k == 0 else "")
        for factor, ls in ((1.5, (0, (5, 2))), (2.0, "-")):
            o = m.overload(factor, j, start="hot")
            tr = o["transient"]
            ax.plot(tr.t / 60, tr.hot_spot, color=ORANGE, lw=2.2, linestyle=ls, label=f"hot spot, {factor:g}x ({factor * j:.1f} A/mm2)")
            if factor == 2.0:
                ax.plot(tr.t / 60, tr.magnet, color=BLUE, lw=2.0, marker="o", markevery=max(1, len(tr.t) // 12), markersize=5, markeredgecolor=SURFACE, label="magnets, 2x")
        ax.axhline(MAGNET_LIMIT_C, color=INK2, lw=1.1, linestyle=(0, (4, 3)))
        ax.axhline(CLASS_F_C, color=INK2, lw=1.1, linestyle=(0, (4, 3)))
        ax.axhline(CLASS_H_C, color=INK2, lw=1.1, linestyle=(0, (4, 3)))
        for y, t in ((MAGNET_LIMIT_C, "80 C"), (CLASS_F_C, "F 155 C"), (CLASS_H_C, "H 180 C")):
            ax.text(xmax_bot[cname] * 0.01, y + 2, t, ha="left", fontsize=7.5, color=INK2)
        ax.set_xlim(0, xmax_bot[cname])
        ax.set_ylim(60, 200)
        ax.text(0.02, 1.02, f"{cname}: overload from hot", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK)
        leg = ax.legend(loc="lower right", frameon=False, fontsize=8)
        for t in leg.get_texts():
            t.set_color(INK2)
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return path


def make_all(rows: list[dict]) -> list[Path]:
    out = OUT()
    out.mkdir(parents=True, exist_ok=True)
    return [fig_network(out / "fig_thermal_network.png"), fig_node_temperatures(out / "fig_thermal_node_temperatures.png"),
            fig_sustainable(out / "fig_thermal_sustainable_current.png", rows), fig_transients(out / "fig_thermal_transients.png", rows)]


if __name__ == "__main__":
    from optimization import thermal_study
    print(*make_all(thermal_study.sustainable_table()), sep="\n")
