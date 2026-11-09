"""optimization.winding_figures: the four pictures of the winding study (matplotlib PNGs drawn from what winding.py computed; no FEMM).

What it is, and why
    A table of 136 front designs and 20 optimizer runs is hard to read; these four pictures are what the study looks at. winding.main() builds the context (the model, the enumeration,
    the grader, the runs, the example winding) and calls make_all(ctx); every picture has its numbers in a CSV beside it (front_exact.csv, nsga_runs.csv, nsga_history.csv): the 'table
    view' the data-visualization method asks for.

The pictures (1430 px wide; each has its title, a one-paragraph note under it that says what is drawn, and a legend)
    fig_exact_front.png   The exact front on two projections that share the loss axis, peak current against loss (left) and corner speed against loss (right), over the other allowed
                          designs in gray (those under 68 W; the rest is off the plot). Both x axes are logarithmic: current and corner speed both go as a / N, so the 120 distinct points
                          are spread by N / a instead of crowding at one end. The example winding is the diamond; the 3000 rpm reference is a hairline. Two panels, not one with two y scales:
                          a second axis would invent a relation between two quantities that only share a parameter.
    fig_nsga_fronts.png   NSGA-II (left) and NSGA-III (right) over the exact front, the same two projections, for the run of median hypervolume ratio of each algorithm (its final
                          population's non-dominated designs). A blue point with a coloured centre was found exactly, a blue point alone was missed, a coloured point alone is a design
                          the exact front beats (a near miss: the same current, a few per cent more loss).
    fig_hv_history.png    Left, the hypervolume ratio of the population against the generation, every seed as a thin line and the mean of each algorithm thick; right, the share of the
                          exact front's points present, in the population (solid) and in the archive (dashed): the grades of the table, generation by generation.
    fig_loss_split.png    The copper loss against the copper fill for three strand sizes (0.2, 0.5, 1.0 mm), and the AC loss alone, from the model's closed forms: the three strands lie on
                          one another in the left panel (DC loss, which only the fill sets), and in the right panel the AC loss is 0.002 ... 0.05 W against 57 W.

The look (the data-visualization method, light surface)
    Palette: surface #fcfcfb, ink #0b0b0b and #52514e, hairline solid grid #e1e0d9, axes #c3c2b7; the exact front is categorical slot 1 (blue #2a78d6), NSGA-II slot 2 (orange #eb6834), NSGA-III
    slot 3 (aqua #1baf7a); everything that is not the subject is gray. The method's validator (validate_palette.py --pairs all, because these are scatter forms) passes the first three slots:
    worst colour-blind Delta E 9.2 (target 8), worst normal-vision 24.0 (floor 15), and warns that aqua is 2.74:1 on the surface, so aqua always has a second channel (triangles) and the
    series are named in titles, legends and the CSV tables. Marks: markers 8 px with a 2 px ring in the surface colour, lines 2 px (seeds 1 px), no dashed grid, text in the ink colours
    and never in a series colour, a legend whenever there are two series, direct labels only where they help (the example, the 3000 rpm line, the end of each curve).

Run
    Not on its own: winding.py --label <label> calls make_all(ctx) and writes the PNGs next to its other outputs, out\\optimization\\winding\\<label>\\ (--no-figures skips them).

Numbers to expect (the reference run, 2026-11-08): four PNGs of 1430 px width: fig_exact_front.png 1430 x 650, fig_nsga_fronts.png 1430 x 1014, fig_hv_history.png 1430 x 611,
fig_loss_split.png 1430 x 611, each 90 ... 180 kB; drawing the four takes about a second.

Self-check (the answers follow).
  1. Why are the current and the speed axes logarithmic?
  2. In fig_nsga_fronts, what does a blue point with an orange centre mean, and what does an orange point with no blue one under it?
  3. Why does fig_loss_split have two panels, and not the AC loss on a second y axis?
Answers.
  1. Both go as a / N, and the front has 120 values of it from 0.25 to 60 (a factor of 240): on a linear axis most of the points would sit within a few per cent of one end.
  2. Blue under orange: the optimizer found that exact point. Blue alone: it missed it. Orange alone: the optimizer found a design the exact front beats, usually the same current as a
     blue point with more loss.
  3. The DC loss is 57 W and the AC loss 0.05 W: on one axis the second would be a line on the floor, and on two axes (one scaled to fill the plot) the reader would be told that the two
     compete, when one is a thousandth of the other. Two panels, each with its own honest scale, say it as it is.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# ---- the look: the data-visualization method's reference palette on its light surface. Categorical slots 1 - 3 (blue, orange, aqua) validate for scatter forms
# ---- (all-pairs: worst CVD Delta E 9.2, worst normal-vision Delta E 24.0); aqua is 2.7:1 on the surface, so it always has a second channel (marker shape)
# ---- and the data are also in CSV tables (the relief the palette asks for).
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"        # exact front, NSGA-II, NSGA-III (slots 1, 2, 3)
CONTEXT = "#898781"                                          # the de-emphasis gray: everything that is not the subject
PT = 0.75                                                    # points per CSS pixel: marks are specified in pixels (marker >= 8 px, line 2 px, ring 2 px)
FONT = ["Segoe UI", "DejaVu Sans"]
ALG_COLOUR = {"NSGA-II": ORANGE, "NSGA-III": AQUA}
ALG_MARKER = {"NSGA-II": "o", "NSGA-III": "^"}
LOSS_TOP = 68.0                                              # the loss axis of the front's pictures ends here [W]; designs above are off the plot


def _style(ax, xlabel: str, ylabel: str, logx: bool = False) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3, color=AXIS)
    ax.grid(True, color=GRID, linewidth=1 * PT, linestyle="-")
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    if logx:
        ax.set_xscale("log")


def _figure(width: float, height: float, title: str, note: str):
    plt.rcParams["font.family"] = FONT
    fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=11, color=INK, fontweight="bold")
    fig.text(0.012, 0.012, note, ha="left", va="bottom", fontsize=7.6, color=INK2, wrap=True)
    return fig


def _legend(ax, handles, loc="upper left") -> None:
    leg = ax.legend(handles=handles, loc=loc, frameon=False, fontsize=8.5, handletextpad=0.5, borderaxespad=0.2)
    for text in leg.get_texts():
        text.set_color(INK2)


def _dot(color, marker="o", size=6.0, label=None, edge=SURFACE) -> Line2D:
    """A legend handle for a scatter mark: its colour with the surface ring."""
    return Line2D([], [], marker=marker, linestyle="", markersize=size, markerfacecolor=color, markeredgecolor=edge, markeredgewidth=2 * PT, label=label)


def _line(color, label) -> Line2D:
    return Line2D([], [], color=color, linewidth=2 * PT, label=label)


def _save(fig, path: Path) -> Path:
    path = Path(path)
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------------------------- 1. the exact front
def fig_exact_front(ctx: dict, path: Path) -> Path:
    """The exact front on two projections that share the loss axis: (peak current, loss) and (corner speed, loss), over the other allowed designs in gray."""
    model, en, ex = ctx["model"], ctx["en"], ctx["example"]
    F = en["e"]["F"]
    allowed = np.flatnonzero(en["feasible"])
    fr = en["front"]
    fig = _figure(11.0, 5.0, "The exact front of the winding study: loss against current, and against corner speed",
                  f"Every combination of N, wire, n_s and a was evaluated ({en['total']:,}; {en['n_feasible']:,} allowed). The front is the {len(fr)} non-dominated designs "
                  f"({ctx['distinct']} distinct points), one for each N / a. Gray: the other allowed designs under {LOSS_TOP:g} W. Geometry '{ctx['label']}', {model.v_dc:g} V bus, "
                  f"copper at {model.t_cu:g} deg C. Table view: front_exact.csv.")
    j, row = ex["index"], ex["row"]
    for k, (col, xlabel, sign) in enumerate(((1, "peak phase current [A]   (the inverter's size)", 1.0), (2, "corner speed at the rated torque [rpm]   (48 V bus)", -1.0))):
        ax = fig.add_subplot(1, 2, k + 1)
        _style(ax, xlabel, "copper loss at 3000 rpm, DC + AC [W]" if k == 0 else "", logx=True)
        ax.set_position([0.07 + 0.495 * k, 0.19, 0.415, 0.67])
        x_all, y_all = sign * F[allowed, col], F[allowed, 0]
        sel = (y_all <= LOSS_TOP) & (x_all > 0)
        ax.scatter(x_all[sel], y_all[sel], s=(5 * PT) ** 2, color=CONTEXT, alpha=0.30, linewidths=0, zorder=1, rasterized=True)
        ax.scatter(sign * F[fr, col], F[fr, 0], s=(8 * PT) ** 2, color=BLUE, edgecolor=SURFACE, linewidth=2 * PT, zorder=3)
        ax.scatter([sign * F[j, col]], [F[j, 0]], s=(10 * PT) ** 2, marker="D", color=ORANGE, edgecolor=SURFACE, linewidth=2 * PT, zorder=4)
        ax.set_ylim(55.8, LOSS_TOP)
        if k == 1:
            ax.axvline(model.rpm_ref, color=AXIS, linewidth=1 * PT, zorder=2)
            ax.text(model.rpm_ref * 1.07, 67.2, f"{model.rpm_ref:g} rpm reference", color=INK2, fontsize=8, va="center")
        ax.annotate(f"example: N = {row['N']}, a = {row['a']}, {row['ns']} x {row['d_mm']:g} mm", xy=(sign * F[j, col], F[j, 0]), xytext=(0.5, 0.62), textcoords="axes fraction",
                    fontsize=8.2, color=INK2, ha="center", bbox={"boxstyle": "round,pad=0.3", "fc": SURFACE, "ec": "none", "alpha": 0.94},
                    arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 1 * PT})
        if k == 0:
            _legend(ax, [_dot(BLUE, label="exact front"), _dot(ORANGE, "D", 7, "example winding"), _dot(CONTEXT, size=4, label="other allowed designs", edge="none")])
    return _save(fig, path)


# ---------------------------------------------------------------------------------------------- 2. NSGA-II and NSGA-III over the exact front
def _median_run(runs: list[dict], name: str) -> dict:
    mine = sorted((r for r in runs if r["algorithm"] == name), key=lambda r: r["records"][-1]["pop"]["hv_ratio"])
    return mine[len(mine) // 2]


def fig_nsga_fronts(ctx: dict, path: Path) -> Path:
    """NSGA-II (left) and NSGA-III (right) over the exact front, on the two projections: the run of median hypervolume ratio of each algorithm."""
    grader, runs = ctx["grader"], ctx["runs"]
    F_ex = grader.front
    n_seeds = len({r["seed"] for r in runs})
    last = runs[0]["records"][-1]
    fig = _figure(11.0, 7.8, "NSGA-II and NSGA-III over the exact front",
                  f"Blue: the exact front's {len(F_ex)} distinct points. Coloured: the non-dominated designs of the final population of the run of median hypervolume ratio "
                  f"(of {n_seeds} seeds; population {ctx['pop']}, {last['n_eval']:,} evaluations). A blue point with a coloured centre was found exactly, a blue point alone was missed, "
                  f"a coloured point alone is a design the exact front beats. Table view: nsga_runs.csv, nsga_history.csv.")
    for c, name in enumerate(("NSGA-II", "NSGA-III")):
        run = _median_run(runs, name)
        rec = run["records"][-1]
        keys = {grader.key(f) for f in run["pop_F"]}
        n_hit = sum(grader.key(f) in keys for f in F_ex)
        for r_, (col, xlabel, sign) in enumerate(((1, "peak phase current [A]", 1.0), (2, "corner speed at the rated torque [rpm]", -1.0))):
            ax = fig.add_subplot(2, 2, 2 * r_ + c + 1)
            _style(ax, xlabel, "copper loss at 3000 rpm [W]" if c == 0 else "", logx=True)
            ax.scatter(sign * F_ex[:, col], F_ex[:, 0], s=(8 * PT) ** 2, color=BLUE, edgecolor=SURFACE, linewidth=2 * PT, zorder=2)
            pf = run["pop_F"]
            ax.scatter(sign * pf[:, col], pf[:, 0], s=(5 * PT) ** 2, marker=ALG_MARKER[name], color=ALG_COLOUR[name], edgecolor=SURFACE, linewidth=1 * PT, zorder=3)
            ax.set_ylim(55.8, LOSS_TOP)
            if r_ == 0:
                ax.set_title(f"{name}, seed {run['seed']}: HV ratio {rec['pop']['hv_ratio']:.3f}, {n_hit} of {len(F_ex)} exact points", loc="left", fontsize=9.5, color=INK,
                             fontweight="bold")
                _legend(ax, [_dot(BLUE, label="exact front"), _dot(ALG_COLOUR[name], ALG_MARKER[name], 5.5, f"{name} final population")], loc="upper right")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.12, wspace=0.12, hspace=0.30)
    return _save(fig, path)


# ---------------------------------------------------------------------------------------------- 3. how the grades grow with the generations
def fig_history(ctx: dict, path: Path) -> Path:
    """Hypervolume ratio (left) and share of the exact points found (right) against the generation, every seed and the mean of each algorithm."""
    runs, grader = ctx["runs"], ctx["grader"]
    n_seeds = len({r["seed"] for r in runs})
    gens = np.array([x["gen"] for x in runs[0]["records"]])
    fig = _figure(11.0, 4.7, "How the optimizers approach the exact front, generation by generation",
                  f"Thin lines: the {n_seeds} seeds of each algorithm; thick: their mean. Left: hypervolume of the population's non-dominated designs over the exact front's (1 = the "
                  f"front). Right: share of the exact front's {len(grader.keys)} distinct points present in the population (solid) and in everything the run has evaluated, its archive "
                  f"(dashed). A population of {ctx['pop']} cannot hold more than {100 * min(1.0, ctx['pop'] / len(grader.keys)):.0f} % of them. Table view: nsga_history.csv.")
    for k, (metric, label, ylim) in enumerate((("hv_ratio", "hypervolume ratio", (0.0, 1.02)), ("hit_share", "share of the exact front found", (0.0, 0.75)))):
        ax = fig.add_subplot(1, 2, k + 1)
        _style(ax, "generation   (100 evaluations each)", label)
        ax.set_position([0.065 + 0.49 * k, 0.19, 0.42, 0.66])
        for name in ("NSGA-II", "NSGA-III"):
            mine = [r for r in runs if r["algorithm"] == name]
            pop = np.array([[x["pop"][metric] for x in r["records"]] for r in mine])
            for y in pop:
                ax.plot(gens, y, color=ALG_COLOUR[name], linewidth=1 * PT, alpha=0.40, zorder=2)
            ax.plot(gens, pop.mean(axis=0), color=ALG_COLOUR[name], linewidth=2 * PT, zorder=3)
            if metric == "hit_share":
                arch = np.array([[x["archive"][metric] for x in r["records"]] for r in mine])
                ax.plot(gens, arch.mean(axis=0), color=ALG_COLOUR[name], linewidth=2 * PT, linestyle="--", zorder=3)
            ax.text(gens[-1] + 1.5, pop.mean(axis=0)[-1], name, color=INK2, fontsize=8.5, va="center")
        ax.axvline(50, color=AXIS, linewidth=1 * PT, zorder=1)
        ax.set_xlim(1, gens[-1] + 18)
        ax.set_ylim(*ylim)
        if k == 0:
            _legend(ax, [_line(ORANGE, "NSGA-II"), _line(AQUA, "NSGA-III")], loc="lower right")
    return _save(fig, path)


# ---------------------------------------------------------------------------------------------- 4. where the loss goes: fill, not strand size
def fig_loss_split(ctx: dict, path: Path) -> Path:
    """The loss against the fill for three wire sizes, and the AC part alone on its own scale (a dual axis would invent a relation: two panels instead)."""
    model = ctx["model"]
    k = np.linspace(0.03, 0.45, 200)
    wires = [(0.2e-3, BLUE), (0.5e-3, ORANGE), (1.0e-3, AQUA)]
    fig = _figure(11.0, 4.7, "Where the copper loss goes: the fill sets it, the strand size hardly matters",
                  f"At a given fill k the DC loss is 1.5 rho l_mt c AT^2 / (k A_layer), the same for every N, a and wire. The AC (proximity) loss is "
                  f"k A_layer n_sides L sigma d^2 omega^2 G / 16 with G = {model.G:.3e} T^2/rad^2 from FEMM, at {model.rpm_ref:g} rpm. Closed forms of the model, not extra solves.")
    ax1, ax2 = fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2)
    ax1.set_position([0.065, 0.19, 0.42, 0.66])
    ax2.set_position([0.565, 0.19, 0.40, 0.66])
    _style(ax1, "copper fill of a coil layer k", "copper loss at 3000 rpm, DC + AC [W]")
    _style(ax2, "copper fill of a coil layer k", "AC (proximity) loss alone [W]")
    p_dc = model.dc_loss_at_fill(k)
    for d, col in wires:
        p_ac = model.n_sides * model.stack * k * model.a_layer * model.sigma * d ** 2 * model.w_ref ** 2 * model.G / 16.0
        ax1.plot(k, p_dc + p_ac, color=col, linewidth=2 * PT, zorder=3)
        ax2.plot(k, p_ac, color=col, linewidth=2 * PT, zorder=3)
        ax2.text(0.455, p_ac[-1], f"{d * 1e3:g} mm", color=INK2, fontsize=8.5, va="center")
    ax1.set_ylim(0, 450)
    ax1.axvline(0.45, color=AXIS, linewidth=1 * PT, zorder=1)
    ax1.text(0.44, 300, "fill limit 0.45", color=INK2, fontsize=8.3, rotation=90, va="center", ha="right")
    ax1.text(0.20, p_dc[np.argmin(abs(k - 0.20))] * 1.07, "the three wire sizes lie on one another", color=INK2, fontsize=8.3)
    ax2.set_xlim(0.03, 0.52)
    ax2.set_ylim(0, None)
    _legend(ax2, [_line(c, f"strand {d * 1e3:g} mm") for d, c in wires], loc="upper left")
    return _save(fig, path)


def make_all(ctx: dict) -> list[Path]:
    """Draw every figure into ctx['folder'] and return the paths."""
    folder = Path(ctx["folder"])
    ctx = {**ctx, "distinct": len(ctx["grader"].front)}
    return [fig_exact_front(ctx, folder / "fig_exact_front.png"), fig_nsga_fronts(ctx, folder / "fig_nsga_fronts.png"), fig_history(ctx, folder / "fig_hv_history.png"),
            fig_loss_split(ctx, folder / "fig_loss_split.png")]
