"""optimization.sweeps_rotor_figures: the pictures of the rotor tier's sweeps, drawn from the CSVs of optimization.sweeps_rotor and nothing else.

What it draws, and why
    The tables say what happened at a handful of points; a picture says what the curve does between them, which is what decides the optimizer. One figure for
    each sweep group, and four more for the questions the sweeps were run to answer:
      sweep_mag_x, sweep_mag_y, sweep_v_half_deg, sweep_guard_mag_y, sweep_guard_v_half_deg, sweep_bridge
                          mean torque above, ripple below, against the variable. The feasible interval's ends are the edge of a light wash that carries the name of
                          the constraint that ends it (or a line, where the end is only the end of the scan); a ring marks the reference's value, a black ring the
                          best point; the fine-mesh re-scores are orange diamonds over the blue coarse line.
      physics_bridge      the radius sweep, nine small panels: the torque, the magnet's torque k psi_pm i_q, the reluctance torque, psi_pm, L_d, L_q, the saliency, the
                          MTPA angle and the power factor, against the outer bridge: what the torque is made of as the V moves toward the guard.
      radius_other_holds  the radius sweep again at four other settings of mag_y and v_half_deg, each over the reference setting in grey: does the shape depend on them?
      noise_floor         each variable's perturbations of the reference: the change of torque and of ripple against the size of the change, with a band of 3 sigma.
      map                 the 2-D map: torque (blue) and ripple (orange) side by side, the ripple's local minima ringed, the infeasible corner named.

How it follows the data-viz rules (the palette is the reference one, light surface)
    One axis to a panel: torque and ripple are two panels on one x, never two scales on one panel, because the way two scales are laid against each other is arbitrary
    and invents a relation. Colour follows the thing, not the rank: blue is the coarse mesh in every figure, orange the fine mesh, grey the context. The first three
    slots of the palette validate in every pair (worst pair 9.2 against the 8 target, normal vision 24.0 against 15); only the first two are used. The map has two
    sequential scales on one figure, so the second takes the next slot's hue (orange, a ramp derived at the blue ramp's lightness steps) and each panel has its own scale
    bar. Lines are 2 px, markers 8 px with a 2 px ring in the surface colour, the grid a solid hairline. A marker is drawn only where it is not within 2.2 % of the axis of
    the last one (spread): where a sweep is dense, the rings would eat the line and make it look dashed. Labels name the extremes and the reference, never every point; text
    wears the ink colours, not the data colour. The radius figures put the guard on the right: the V moves outward to the right, as in the question asked.
    Not done: a table view and hover (these are PNGs; sweeps_rotor's tables.txt and summary.json carry every number).

Run (after optimization.sweeps_rotor analyse; no FEMM)
    .venv/Scripts/python.exe -m optimization.sweeps_rotor_figures
Writes ten PNGs of 85 to 160 kB (130 dpi) to out\\optimization\\sweeps_rotor\\ and prints their paths.

Numbers to expect (this PC, 2026-11-08): the six sweep figures draw 22, 20, 24, 21, 25 and 30 points; the map is 13 columns by 9 rows with 6 blank cells; the ripple panels
leave out the fine-mesh points of torque under 0.05 N.m; the noise figure's band is plus and minus 3 x 0.559 mN.m on the torque and 3 x 0.144 points on the ripple.

Self-check (the answers follow).
  1. Why are the torque and the ripple on two stacked panels and not on one panel with a scale on each side?
  2. The ripple panel of sweep_bridge does not show the fine-mesh point at a bridge of 15.74 mm. Why not?
  3. The map's torque is blue and its ripple orange. Why not two blue ramps?
Answers.
  1. With two scales the reader cannot tell what is a coincidence of where the scales were set: moving one axis moves the curve across the other. On stacked panels each
     measure has its one honest axis and the shared x still lines a torque up with the ripple at the same value of the variable.
  2. Its torque is 0.0017 N.m, and the ripple in per cent of that is a ratio of nearly nothing: the fine mesh says 18 % where the coarse one says 3 %, and neither is a fact about the
     motor. Points under 0.05 N.m are left out of the ripple panel (the torque panel still shows them).
  3. Two scales on one figure in the same hue would read as one scale. The second sequential ramp takes another hue, so a dark blue cell and a dark orange cell can never be
     mistaken for the same quantity.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from optimization import sweeps_rotor as sr

# ---------------------------------------------------------------------------------------------- the look (data-visualization palette, light surface)
SURFACE, INK, INK2, MUTED, GRID, AXIS, WASH = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#eeede8"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"      # categorical slots 1 to 3: the only three that validate in every pair of a small-multiples chart
GREY = "#b9b8b0"                                          # the context series of an emphasis chart
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]      # sequential, steps 100 to 700 of the palette's blue
ORANGE_RAMP = ["#f8d8cc", "#edb29d", "#e08a6b", "#cf6238", "#ab4821", "#853412", "#5e2208"]    # the second sequential ramp: slot 2's hue at the blue ramp's lightness steps
LINE, MARKER, RING = 1.5, 6.0, 1.5                         # points: 2 px line, 8 px marker, 2 px surface ring (96 px to the inch)

FAMILY = "Segoe UI" if any(f.name == "Segoe UI" for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS, "axes.linewidth": 0.75,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.75, "grid.linestyle": "-", "axes.axisbelow": True,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.5, "legend.frameon": False,
    "font.family": FAMILY, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.major.width": 0.75, "ytick.major.width": 0.75, "figure.dpi": 100,
})
BLUE_MAP = LinearSegmentedColormap.from_list("blue_seq", BLUE_RAMP)
ORANGE_MAP = LinearSegmentedColormap.from_list("orange_seq", ORANGE_RAMP)

LABELS = {"mag_x": "mag_x [mm]  (slot centre along the d-axis)", "mag_y": "mag_y [mm]  (slot centre from the d-axis)",
          "v_half_deg": "v_half_deg [deg]  (magnet against the d-axis)", "bridge_mm": "outer bridge [mm]  (the radius of the V)"}
TITLES = {"mag_x": "mag_x alone: the V moves out along the d-axis, the others at the reference",
          "mag_y": "mag_y alone: the two magnets of a V move apart, the others at the reference",
          "v_half_deg": "v_half_deg alone: the V opens, the others at the reference",
          "guard_mag_y": "mag_y on the guard: the outer bridge held at 1.4 mm, mag_x follows",
          "guard_v_half_deg": "v_half_deg on the guard: the outer bridge held at 1.4 mm, mag_x follows",
          "bridge": "the radius sweep: the outer bridge from the guard inward, mag_y and v_half_deg at the reference"}


def out_path(name: str) -> Path:
    return sr.folder() / name


def finish(fig, name: str) -> Path:
    path = out_path(name)
    fig.savefig(path, dpi=130, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    return path


def style_axes(ax) -> None:
    ax.tick_params(length=3, width=0.75)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.75)


def spread(x, gap_fraction: float = 0.022) -> np.ndarray:
    """A mask over the points that are at least gap_fraction of the axis span from the last marker drawn: where the sweep is dense (the guard end) the
    rings of neighbouring markers would eat the line and make it look dashed, so there the line alone carries the points."""
    x = np.asarray(x, float)
    keep, last = np.zeros(len(x), bool), None
    for i in np.argsort(x):
        if last is None or abs(x[i] - last) >= gap_fraction * (x.max() - x.min()):
            keep[i], last = True, x[i]
    return keep


def ring_markers(ax, x, y, color, size=MARKER, **kw):
    """Dots with the 2 px ring in the surface colour, so that they stay legible where they cross a line, on the points that are not crowded."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = spread(x)
    return ax.plot(x[m], y[m], "o", color=color, markersize=size, markeredgecolor=SURFACE, markeredgewidth=RING, zorder=3, **kw)


LIMIT_NOTES = {"guard_bridge": "guard_bridge: outer bridge under 1.4 mm", "centre_rib": "centre_rib: rib between the magnets under 1.0 mm",
               "q_rib": "q_rib: steel between poles under 1.0 mm"}


def headline(fig, title: str, sub: str) -> None:
    fig.text(0.005, 0.995, title, ha="left", va="top", fontsize=11, color=INK, fontweight="bold")
    fig.text(0.005, 0.955, sub, ha="left", va="top", fontsize=8.5, color=INK2)


# ---------------------------------------------------------------------------------------------- one sweep: torque and ripple against the variable
INVERTED = ("bridge",)                 # the radius sweeps are drawn with the guard on the right: the V moves outward to the right, as in the question asked


def sweep_figure(name: str, rows: list[dict], fine: list[dict], iv: dict, anchor: float | None, noise: dict | None) -> Path:
    """Two stacked panels on one x axis (never two scales on one panel): mean torque above, ripple below. The feasible interval's ends are the edge
    of a light wash with the constraint's name (or a line, where the end is only the scan's own limit); a ring at the reference's value; the best points
    labelled and nothing else; the fine-mesh re-scores over the coarse line, with a legend because there are then two series."""
    var = sr.SWEEPS[name]["var"]
    x = np.array([r["x"] for r in rows])
    fig, (ax_t, ax_r) = plt.subplots(2, 1, sharex=True, figsize=(9.2, 6.6), gridspec_kw={"hspace": 0.1, "top": 0.86, "bottom": 0.1, "left": 0.09, "right": 0.985})
    pad = 0.06 * (x[-1] - x[0])
    xlo, xhi = x[0] - pad, x[-1] + pad
    ends = ((iv["lo"], iv["lo_limit"], xlo), (iv["hi"], iv["hi_limit"], xhi))
    for ax in (ax_t, ax_r):
        style_axes(ax)
        ax.set_xlim(*((xhi, xlo) if name in INVERTED else (xlo, xhi)))
        for edge, limit, far in ends:
            if limit == "scan limit":
                ax.axvline(edge, color=AXIS, lw=0.75, zorder=1)
            else:
                ax.axvspan(*sorted((edge, far)), color=WASH, lw=0, zorder=0)
    clipped = np.array([bool(r["clipped"]) for r in rows])
    for ax, key, ylabel, kind in ((ax_t, "T_mean", "mean torque  [N.m]", "max"), (ax_r, "ripple_pct", "torque ripple  [% of the mean]", "min")):
        y = np.array([r[key] for r in rows])
        ax.plot(x, y, "-", color=BLUE, lw=LINE, zorder=2, label="coarse mesh (every point)")
        ring_markers(ax, x[~clipped], y[~clipped], BLUE)
        if clipped.any():
            m = spread(x[clipped])
            ax.plot(x[clipped][m], y[clipped][m], "o", mfc=SURFACE, mec=BLUE, mew=1.2, ms=MARKER, zorder=3, label="MTPA angle at the edge of its 10 to 40 deg window")
        shown = [r for r in fine if key == "T_mean" or r["T_mean"] >= 0.05]       # the ripple in per cent of a torque near zero says nothing
        if shown:
            ax.plot([r["x"] for r in shown], [r[key] for r in shown], "D", color=ORANGE, ms=MARKER - 0.5, mec=SURFACE, mew=RING, zorder=4, label="fine mesh (re-scored)")
        best = int(np.argmax(y) if kind == "max" else np.argmin(y))
        ax.plot(x[best], y[best], "o", mfc="none", mec=INK, mew=1.1, ms=MARKER + 5, zorder=5)
        right = (x[best] > x[0] + 0.5 * (x[-1] - x[0])) != (name in INVERTED)         # which side of the panel the point is on, as drawn
        text = f"highest {y[best]:.4f} N.m at {x[best]:.4g}" if kind == "max" else f"lowest {y[best]:.2f} % at {x[best]:.4g}"
        ax.annotate(text, (x[best], y[best]), xytext=(-14 if right else 14, 12), textcoords="offset points", ha="right" if right else "left", va="bottom",
                    fontsize=8.5, color=INK2)
        if anchor is not None and x[0] <= anchor <= x[-1]:
            i = int(np.argmin(np.abs(x - anchor)))
            ax.plot(x[i], y[i], "o", mfc="none", mec=INK2, mew=0.9, ms=MARKER + 2.5, zorder=5)
            if ax is ax_t:
                side_right = (x[i] > x[0] + 0.5 * (x[-1] - x[0])) != (name in INVERTED)
                ax.annotate("reference", (x[i], y[i]), xytext=(-10 if side_right else 10, -16), textcoords="offset points", ha="right" if side_right else "left",
                            fontsize=8, color=INK2)
        ax.set_ylabel(ylabel)
        pad_y = 0.08 * (y.max() - y.min() or 1.0)
        ax.set_ylim(y.min() - pad_y, y.max() + 2.2 * pad_y)
    ax_r.set_xlabel(LABELS[var] + ("  -  the V moves outward to the right" if name in INVERTED else ""))
    for edge, limit, far in ends:                                       # the name of each limit stands in the margin beside it, clear of the data
        if limit == "scan limit":
            label = "the guard: bridge_mm 1.4" if name == "bridge" else "end of the scanned range"
        else:
            label = LIMIT_NOTES.get(limit, limit)
        ax_t.text((edge + far) / 2, 0.97, label, transform=ax_t.get_xaxis_transform(), rotation=90, ha="center", va="top", fontsize=7.5, color=MUTED)
    handles, labels = ax_t.get_legend_handles_labels()
    if len(labels) > 1:
        ax_t.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=len(labels), handletextpad=0.4, columnspacing=1.6)
    sub = f"FEMM 5.0, coarse mesh, 10 A at each point's own MTPA angle, scoring AT, {len(rows)} points"
    if noise:
        sub += f"; noise floor {1e3 * noise['pooled_sigma']['T_mean']:.2g} mN.m on the torque, {noise['pooled_sigma']['ripple_pct']:.2g} points on the ripple"
    headline(fig, TITLES[name], sub)
    return finish(fig, f"sweep_{name}.png")


# ---------------------------------------------------------------------------------------------- the physics of the radius sweep
PHYSICS_PANELS = (("T_mean", "mean torque  [N.m]", 1.0), ("T_pm", "magnet torque  k psi_pm i_q  [N.m]", 1.0), ("T_rel", "reluctance torque  [N.m]", 1.0),
                  ("psi_pm", "magnet flux linkage psi_pm  [mWb]", 1e3), ("L_d", "L_d  [mH]", 1e3), ("L_q", "L_q  [mH]", 1e3),
                  ("saliency", "saliency  L_q / L_d", 1.0), ("gamma", "MTPA angle  [deg]", 1.0), ("pf", "power factor", 1.0))


def physics_figure(rows: list[dict], xmax: float = 8.0) -> Path:
    """The torque and what it is made of, against the outer bridge: nine small panels, one series each, so each carries its own scale and no
    panel needs a legend (its title names it)."""
    rows = [r for r in rows if r["x"] <= xmax]
    x = np.array([r["x"] for r in rows])
    fig, axes = plt.subplots(3, 3, sharex=True, figsize=(10.4, 8.0), gridspec_kw={"hspace": 0.16, "wspace": 0.32, "top": 0.9, "bottom": 0.07, "left": 0.07, "right": 0.985})
    for ax, (key, label, scale) in zip(axes.flat, PHYSICS_PANELS):
        style_axes(ax)
        y = np.array([r[key] for r in rows]) * scale
        ax.plot(x, y, "-", color=BLUE, lw=LINE, zorder=2)
        ring_markers(ax, x, y, BLUE, size=MARKER - 1)
        ax.set_title(label, fontsize=9, loc="left", color=INK2)
        ax.set_xlim(x.max() + 0.15, x.min() - 0.15)                     # the guard on the right: the V moves outward to the right
    for ax in axes[-1]:
        ax.set_xlabel("outer bridge [mm]  (guard 1.4)")
    headline(fig, "What the torque is made of, as the V moves outward toward the guard",
             f"the radius sweep (mag_y and v_half_deg at the reference), outer bridge 1.4 to {xmax:g} mm, coarse mesh, 10 A at each point's own MTPA angle; the guard is at the right")
    return finish(fig, "physics_bridge.png")


# ---------------------------------------------------------------------------------------------- the radius sweeps at other values of the other two variables
def radius_other_holds(base: list[dict], cross: dict, xmax: float = 5.0) -> Path | None:
    """Does the shape of the radius sweep depend on where mag_y and v_half_deg sit? One column for each other setting: its own curve in blue, the
    reference setting's in grey for comparison (emphasis: one series is the point, the rest is context). Torque above, ripple below, the guard at the right."""
    if not cross:
        return None
    names = list(cross)
    fig, axes = plt.subplots(2, len(names), sharex=True, figsize=(2.9 * len(names) + 1.2, 6.2), squeeze=False,
                             gridspec_kw={"hspace": 0.1, "wspace": 0.3, "top": 0.84, "bottom": 0.1, "left": 0.08, "right": 0.985})
    b = [r for r in base if r["x"] <= xmax]
    for col, name in enumerate(names):
        rows = [r for r in cross[name]["rows"] if r["x"] <= xmax]
        for row, (key, ylabel) in enumerate((("T_mean", "mean torque  [N.m]"), ("ripple_pct", "torque ripple  [% of the mean]"))):
            ax = axes[row][col]
            style_axes(ax)
            ax.plot([r["x"] for r in b], [r[key] for r in b], "-", color=GREY, lw=LINE, zorder=1)
            ax.plot([r["x"] for r in rows], [r[key] for r in rows], "-", color=BLUE, lw=LINE, zorder=2)
            ring_markers(ax, [r["x"] for r in rows], [r[key] for r in rows], BLUE, size=MARKER - 1)
            ax.set_xlim(xmax + 0.1, 1.3)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == 0:
                ax.set_title(name.replace("bridge@", "").replace(",", ",  ").replace("=", " "), fontsize=9, loc="left", color=INK2)
            if row == 1:
                ax.set_xlabel("outer bridge [mm]")
    from matplotlib.lines import Line2D
    fig.legend([Line2D([0], [0], color=BLUE, lw=LINE), Line2D([0], [0], color=GREY, lw=LINE)], ["the setting named above each column", "the reference setting (mag_y 6, v_half_deg 60)"],
               loc="upper right", bbox_to_anchor=(0.985, 0.935), ncol=2, handlelength=1.6)
    headline(fig, "The radius sweep at other settings of mag_y and v_half_deg",
             "the guard (1.4 mm) is at the right; the V moves outward to the right; coarse mesh, 10 A at each point's own MTPA angle, scoring AT")
    return finish(fig, "radius_other_holds.png")


# ---------------------------------------------------------------------------------------------- the noise floor
def noise_figure(rows: list[dict], noise: dict) -> Path:
    """What a hair's change of the geometry does to the score: each variable's perturbations of the reference, the change of torque (above) and of ripple
    (below) from the reference, on a symmetric-log axis (the ladder runs from 0.001 to 0.05 mm), with a band of 3 sigma around zero."""
    ref = next(r for r in rows if r["stage"] == "reference" and r["feasible"])
    fig, axes = plt.subplots(2, 3, figsize=(10.4, 6.2), gridspec_kw={"hspace": 0.2, "wspace": 0.32, "top": 0.86, "bottom": 0.1, "left": 0.09, "right": 0.985})
    units = {"mag_x": "mm", "mag_y": "mm", "v_half_deg": "deg"}
    for col, var in enumerate(sr.NOISE_H):
        pts = sorted([(h, r) for r in rows if r["stage"] == "noise" and r["feasible"] and not r["cached"] for v, h in [sr.perturbed(r) or (None, None)] if v == var], key=lambda t: t[0])
        for row, (key, ylabel, scale) in enumerate((("T_mean", "change of mean torque  [mN.m]", 1e3), ("ripple_pct", "change of ripple  [points]", 1.0))):
            ax = axes[row][col]
            style_axes(ax)
            sigma = noise["variables"][var]["quantities"][key]["sigma"] * scale
            ax.axhspan(-3 * sigma, 3 * sigma, color=WASH, lw=0, zorder=0)
            ax.axhline(0, color=AXIS, lw=0.75, zorder=1)
            h = np.array([p[0] for p in pts])
            y = (np.array([p[1][key] for p in pts]) - ref[key]) * scale
            ax.set_xscale("symlog", linthresh=sr.NOISE_H[var][0])
            ax.plot(h, y, "o", color=BLUE, ms=MARKER - 0.5, mec=SURFACE, mew=RING, zorder=3)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == 1:
                ax.set_xlabel(f"{var} minus its reference value [{units[var]}]")
            if row == 0:
                ax.set_title(var, fontsize=9, loc="left", color=INK2)
    headline(fig, "The noise floor: the score under a hair's change of the geometry",
             "coarse mesh; each point is the reference with one variable moved; the band is plus and minus 3 sigma, sigma from the points within 0.003")
    return finish(fig, "noise_floor.png")


# ---------------------------------------------------------------------------------------------- the 2-D map
def map_figure(spec: dict, grids: dict, minima: list[dict], anchor: tuple[float, float] | None) -> Path:
    """The two outputs of the 2-D map side by side, each its own one-hue ramp (torque blue, ripple orange) with its own scale: the ripple's local
    minima ringed, the highest torque marked."""
    xs, ys = np.array(spec["xs"]), np.array(spec["ys"])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0), gridspec_kw={"wspace": 0.28, "top": 0.86, "bottom": 0.12, "left": 0.07, "right": 0.985})
    dx, dy = (xs[1] - xs[0]) / 2, (ys[1] - ys[0]) / 2
    ex, ey = np.append(xs - dx, xs[-1] + dx), np.append(ys - dy, ys[-1] + dy)
    for ax, key, cmap, label in ((axes[0], "T_mean", BLUE_MAP, "mean torque  [N.m]"), (axes[1], "ripple_pct", ORANGE_MAP, "torque ripple  [% of the mean]")):
        z = np.array(grids[key], float)
        mesh = ax.pcolormesh(ex, ey, np.ma.masked_invalid(z), cmap=cmap, shading="flat", edgecolors=SURFACE, linewidth=1.0)
        ax.set_xlim(ex[0], ex[-1])
        ax.set_ylim(ey[0], ey[-1])
        ax.grid(False)
        ax.set_xlabel(LABELS.get(spec["x"], spec["x"]))
        ax.set_ylabel(LABELS.get(spec["y"], spec["y"]))
        cb = fig.colorbar(mesh, ax=ax, fraction=0.05, pad=0.03)
        cb.set_label(label, color=INK2)
        cb.outline.set_edgecolor(AXIS)
        cb.ax.tick_params(colors=MUTED, labelsize=8.5)
        for s in ax.spines.values():
            s.set_edgecolor(AXIS)
        if key == "T_mean":
            i, j = np.unravel_index(np.nanargmax(z), z.shape)
            ax.plot(xs[j], ys[i], "o", mfc="none", mec=INK, mew=1.4, ms=11)
            ax.annotate(f"highest {z[i, j]:.3f} N.m", (xs[j], ys[i]), xytext=(-12, 12) if xs[j] > xs.mean() else (12, 12), textcoords="offset points", fontsize=8.5, color=INK,
                        va="bottom", ha="right" if xs[j] > xs.mean() else "left", bbox={"boxstyle": "round,pad=0.2", "fc": SURFACE, "ec": "none", "alpha": 0.85})
        else:
            for m in minima:
                ax.plot(m["x"], m["y"], "o", mfc="none", mec=INK, mew=1.4, ms=11)
            if minima:
                low = min(minima, key=lambda m: m["ripple_pct"])
                ax.annotate(f"lowest {low['ripple_pct']:.2f} %", (low["x"], low["y"]), xytext=(-12, 12) if low["x"] > xs.mean() else (12, 12), textcoords="offset points",
                            fontsize=8.5, color=INK, va="bottom", ha="right" if low["x"] > xs.mean() else "left",
                            bbox={"boxstyle": "round,pad=0.2", "fc": SURFACE, "ec": "none", "alpha": 0.85})
        blank = np.argwhere(~np.isfinite(z))
        if len(blank):                                                    # the cells that are not scored: infeasible, and named
            ax.text(xs[blank[:, 1]].mean(), ys[blank[:, 0]].mean(), "infeasible:\ncentre_rib\nunder 1.0 mm", ha="center", va="center", fontsize=8, color=MUTED)
        if anchor:
            ax.plot(*anchor, "x", color=INK, ms=7, mew=1.4)
    hold = ", ".join(f"{k} {v:g}" for k, v in spec["hold"].items()) or "nothing"
    headline(fig, f"The 2-D map: {spec['y']} against {spec['x']}",
             f"{len(xs)} x {len(ys)} points, held: {hold}; rings: local minima of the ripple (3 noise floors below every neighbour) and the highest torque; coarse mesh, 10 A, scoring AT")
    return finish(fig, "map.png")


# ---------------------------------------------------------------------------------------------- all of them
def main() -> int:
    coarse = sr.read_csv("evaluations.csv")
    sr.split_torque([r for r in coarse if r["feasible"]])
    sweeps = sr.load_sweeps(rows=coarse)
    scan = json.loads((sr.folder() / "scan.json").read_text(encoding="utf-8"))
    noise = sr.noise_floor(coarse) if any(r["stage"] == "noise" for r in coarse) else None
    fine_path = sr.folder() / "evaluations_fine.csv"
    fine_all = sr.read_csv("evaluations_fine.csv") if fine_path.exists() else []
    made = []
    for name in sr.SWEEP_ORDER:
        rows = sweeps[name]
        if len(rows) < 4:
            continue
        iv = next(i for i in scan[name]["intervals"] if i["contains_reference"])
        fine = [dict(r, x=sr.sweep_x(name, r)) for r in fine_all if r["stage"] == f"fine_{name}" and r["feasible"]]
        made.append(sweep_figure(name, rows, fine, iv, sr.SWEEPS[name]["anchor"], noise))
    if len(sweeps["bridge"]) > 3:
        made.append(physics_figure(sweeps["bridge"]))
    cross = {}
    for stage in sorted({r["stage"] for r in coarse if r["stage"].startswith("bridge@")}):
        rows = sorted([dict(r, x=r["bridge_mm"]) for r in coarse if r["stage"] == stage and r["feasible"]], key=lambda r: r["x"])
        if len(rows) > 3:
            cross[stage] = {"rows": rows}
    figure = radius_other_holds(sweeps["bridge"], cross)
    if figure:
        made.append(figure)
    if noise:
        made.append(noise_figure(coarse, noise))
    map_path = sr.folder() / "map.json"
    if map_path.exists() and any(r["stage"] == "map" for r in coarse):
        spec = json.loads(map_path.read_text(encoding="utf-8"))
        grids = sr.map_grids(coarse, spec)
        sig_R = noise["pooled_sigma"]["ripple_pct"] if noise else 0.0
        minima = [{"x": spec["xs"][j], "y": spec["ys"][i], "ripple_pct": v} for i, j, v in sr.local_minima_2d(grids["ripple_pct"], 3.0 * sig_R)]
        made.append(map_figure(spec, {k: np.where(np.isfinite(v), v, np.nan) for k, v in grids.items()}, minima, None))
    for path in made:
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
