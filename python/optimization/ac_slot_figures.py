"""optimization.ac_slot_figures: the pictures of the slot studies (matplotlib PNGs drawn from what optimization/ac_slot_study.py computed; no FEMM).

What it is, and why
    The slot studies give hundreds of numbers; these figures are what they say. Every figure has its title, a note under it that says what is drawn, units on its axes, and its numbers
    in a table beside it (sweep.csv, circ.json, ...): the 'table view'.

The pictures (written by ac_slot_study.py into out/optimization/copper_ac/slot/)
    fig_slot_rac.png       R_ac / R_dc of the slot against frequency, for the two windings, layers in phase and 60 degrees apart (left: the whole range, log-log in the extra part
                           R_ac / R_dc - 1, with the machine's frequencies marked; right: the ratio itself on a linear axis up to the PWM frequency).
    fig_slot_collapse.png  (R_ac / R_dc - 1) / (d k / delta^2)^2 against d / delta for the two windings: if the extra loss followed one law of the strand and the fill the points would lie on one
                           horizontal line; where they leave it, something else (the steel, the strands' self-shielding) acts.
    fig_slot_formula.png   the winding study's formula against FEMM, strand by strand: the ratio of FEMM's excess loss to the formula against frequency (total, and the strands next to the opening)
                           and, at 20 kHz, every strand of one layout as a point (formula against FEMM).
    fig_slot_maps.png      the current density inside the strands (top: a zoom to the opening) and the loss of every strand over its DC loss (bottom: the whole slot) at three cases.
    fig_slot_circ.png      the loss added by connecting the strands of a turn in parallel, over the DC loss, against frequency, for four placements of the strands in the turns.
    The look follows optimization/winding_figures.py: surface #fcfcfb, ink #0b0b0b and #52514e, grid #e1e0d9, axes #c3c2b7; series blue #2a78d6, orange #eb6834, aqua #1baf7a; context gray #898781.

Run
    Not on its own: ac_slot_study.py calls these after its stages. Each takes a second or two (the maps a few seconds: 80 thousand triangles).

Self-check (the answers follow).
  1. Why is the extra loss R_ac / R_dc - 1 drawn on a logarithmic axis, and not the ratio?
  2. In the collapse figure the points of the two windings sit on one line at low d / delta. What does that say about the slot, and what would a spread say?
Answers.
  1. At the fundamental the extra loss is 1e-4 of the DC loss and at the switching frequency three times it: seven orders of magnitude across the plot; on a linear axis the first would be a
     line on the floor. The ratio itself (right panel) is the number the loss model uses, the extra loss on a log scale is how it grows (as f^2 while the strands are thin).
  2. That one geometry constant K [m^2] describes the slot's field at the strands: the extra loss of any winding is K (d k / delta^2)^2. A spread at low d / delta would say that the field depends
     on how the copper is arranged, not just on how much there is: it would not collapse.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.patches import Circle

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
PT = 0.75
FONT = ["Segoe UI", "DejaVu Sans"]
CMAP = LinearSegmentedColormap.from_list("surface_blue_ink", [SURFACE, "#9ec4ef", BLUE, "#0b2a52"])
WCOLOUR = {"W1": BLUE, "W2": ORANGE}


def _style(ax, xlabel: str, ylabel: str) -> None:
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


def _figure(width: float, height: float, title: str, note: str):
    plt.rcParams["font.family"] = FONT
    fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=11, color=INK, fontweight="bold")
    fig.text(0.012, 0.012, note, ha="left", va="bottom", fontsize=7.6, color=INK2, wrap=True)
    return fig


def _legend(ax, loc="best", **kw) -> None:
    leg = ax.legend(loc=loc, frameon=False, fontsize=8.2, handletextpad=0.5, borderaxespad=0.2, **kw)
    for text in leg.get_texts():
        text.set_color(INK2)


def _save(fig, path) -> None:
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def _colorbar(fig, ax, mappable, label: str | None = None, ticks=None, gap: float = 0.007, width: float = 0.008):
    """A colour bar as tall as the axes' DRAWN box (an equal-aspect axes is smaller than the cell it was given), just to its right. `ticks` (optional) are written as plain numbers."""
    ax.apply_aspect()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + gap, pos.y0, width, pos.height])
    cb = fig.colorbar(mappable, cax=cax)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    cb.outline.set_edgecolor(AXIS)
    if ticks is not None:
        cb.set_ticks(list(ticks))
        cb.ax.set_yticklabels([f"{t:g}" for t in ticks])
        cb.ax.minorticks_off()
    if label:
        cb.set_label(label, fontsize=8.5, color=INK2)
    return cb


def _plain_log_ticks(vmin: float, vmax: float) -> list[float]:
    """1, 2, 5, 10, 20, 50 ... inside [vmin, vmax] (for a logarithmic colour bar)."""
    return [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000) if vmin <= t <= vmax]


KEY = (("200 Hz\n3,000 rpm", 200.0), ("1.09 kHz\n16,330 rpm", 1090.0), ("10 kHz", 1e4), ("20 kHz", 2e4))


def _mean_by_f(rows, winding, case, key):
    fs = sorted({r["f_Hz"] for r in rows if r["winding"] == winding and r["case"] == case})
    mean, lo, hi = [], [], []
    for f in fs:
        vals = [r[key] for r in rows if r["winding"] == winding and r["case"] == case and r["f_Hz"] == f]
        mean.append(np.mean(vals))
        lo.append(np.min(vals))
        hi.append(np.max(vals))
    return np.array(fs), np.array(mean), np.array(lo), np.array(hi)


def fig_rac(an: dict, path: Path) -> None:
    rows = an["rows"]
    fig = _figure(11.0, 5.2, "R_ac / R_dc of one slot of the optimized stator, strand by strand in FEMM",
                  "Both windings carry the rated 118.6 ampere-turns a coil side (18.25 A peak a phase) with every strand sharing the current equally (series circuits); copper at 100 degC; "
                  "'in phase' = both layers of the slot carry one phase (6 of the 9 slots), '60 deg' = two phases (3 of 9). Band: the 3 random layouts. FEMM 5.0 time-harmonic, 416 strands "
                  "(W1) or 104 (W2) in a 3-slot-pitch cell. Table view: sweep.csv.")
    ax = fig.add_subplot(1, 2, 1)
    _style(ax, "frequency [Hz]", "R_ac / R_dc - 1   (the extra loss over the DC loss)")
    ax.set_position([0.07, 0.22, 0.42, 0.66])
    for wn, lab in (("W1", "W1: 26 turns x 8 x 0.315 mm"), ("W2", "W2: 13 turns x 4 x 0.630 mm")):
        fs, mean, lo, hi = _mean_by_f(rows, wn, "same", "R")
        ax.fill_between(fs, lo - 1, hi - 1, color=WCOLOUR[wn], alpha=0.18, linewidth=0)
        ax.plot(fs, mean - 1, color=WCOLOUR[wn], lw=2 * PT, marker="o", ms=4.5, mec=SURFACE, mew=1.2, label=f"{lab}, layers in phase")
        fd, md, _, _ = _mean_by_f(rows, wn, "diff", "R")
        ax.plot(fd, md - 1, color=WCOLOUR[wn], lw=1.6 * PT, ls="--", label=f"{lab}, layers 60 deg apart")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(3e-8, 1e3)
    for text, f in KEY:
        ax.axvline(f, color=AXIS, lw=1 * PT, zorder=0)
        ax.text(f * 1.06, 3.5e-8 * 1.6, text, fontsize=7.2, color=INK2, va="bottom")
    _legend(ax, "upper left")
    ax = fig.add_subplot(1, 2, 2)
    _style(ax, "frequency [Hz]", "R_ac / R_dc")
    ax.set_position([0.58, 0.22, 0.40, 0.66])
    for wn in ("W1", "W2"):
        fs, mean, lo, hi = _mean_by_f(rows, wn, "same", "R")
        sel = fs <= 6.5e4
        ax.fill_between(fs[sel], lo[sel], hi[sel], color=WCOLOUR[wn], alpha=0.18, linewidth=0)
        ax.plot(fs[sel], mean[sel], color=WCOLOUR[wn], lw=2 * PT, marker="o", ms=4.5, mec=SURFACE, mew=1.2, label=f"{wn}, layers in phase")
        fd, md, _, _ = _mean_by_f(rows, wn, "diff", "R")
        ax.plot(fd[fd <= 6.5e4], md[fd <= 6.5e4], color=WCOLOUR[wn], lw=1.6 * PT, ls="--", label=f"{wn}, layers 60 deg apart")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e3, 6.5e4)
    ax.set_ylim(0.97, 200)
    ax.set_yticks([1, 2, 5, 10, 20, 50, 100])
    ax.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"])
    for text, f in KEY[2:]:
        ax.axvline(f, color=AXIS, lw=1 * PT, zorder=0)
        ax.text(f * 1.04, 1.02, text, fontsize=7.2, color=INK2, va="bottom")
    _legend(ax, "upper left")
    _save(fig, path)


def fig_collapse(an: dict, path: Path) -> None:
    rows = an["rows"]
    fig = _figure(9.0, 5.0, "One law for the slot's extra loss? (R_ac / R_dc - 1) / (d k / delta^2)^2 against d / delta",
                  "K = (R_ac / R_dc - 1) / (d k / delta^2)^2 [mm^2] with the strand diameter d, the copper fill k (0.448 for both) and the skin depth delta. Where K is flat the extra loss is the slot's field "
                  "acting on strands thin against the skin depth; it falls where the strands shield themselves and the steel's eddy currents act. Layers in phase, mean of 3 layouts.")
    ax = fig.add_subplot(1, 1, 1)
    _style(ax, "d / delta  (strand diameter over skin depth at that frequency)", "K [mm^2]")
    ax.set_position([0.09, 0.22, 0.88, 0.66])
    for wn, lab in (("W1", "W1: 0.315 mm strands"), ("W2", "W2: 0.630 mm strands")):
        for case, ls, extra in (("same", "-", "layers in phase"), ("diff", "--", "layers 60 deg apart")):
            sel = [r for r in rows if r["winding"] == wn and r["case"] == case]
            xs = sorted({r["d_over_delta"] for r in sel})
            ys = [1e6 * np.mean([r["K"] for r in sel if r["d_over_delta"] == x]) for x in xs]
            m = np.array(xs) >= 0.012
            ax.plot(np.array(xs)[m], np.array(ys)[m], ls=ls, color=WCOLOUR[wn], lw=2 * PT, marker="o", ms=4.5, mec=SURFACE, mew=1.2, label=f"{lab}, {extra}")
    ax.set_xscale("log")
    _legend(ax, "lower left")
    _save(fig, path)


def fig_formula(an: dict, res: dict, path: Path) -> None:
    rows = an["rows"]
    tables = an["tables"]
    fig = _figure(11.0, 5.4, "The winding study's formula against FEMM, strand by strand",
                  "Left: FEMM's loss above DC over the formula pi sigma w^2 d^4 L |B|^2 / 128 with B the field of the SAME strands in a magnetostatic solve (W1, W2) and, dashed, with the lone wire's exact "
                  "proximity and skin factors. Right: every strand of one W1 layout at 20 kHz (d / delta = 0.59), FEMM against the formula, coloured by depth below the coil's top edge. Layers in phase.")
    ax = fig.add_subplot(1, 2, 1)
    _style(ax, "frequency [Hz]", "FEMM excess loss / formula")
    ax.set_position([0.07, 0.24, 0.42, 0.63])
    for wn, lab in (("W1", "W1 (0.315 mm)"), ("W2", "W2 (0.630 mm)")):
        sel = [r for r in rows if r["winding"] == wn and r["case"] == "same"]
        fs = sorted({r["f_Hz"] for r in sel})
        mean = lambda key: np.array([np.mean([r[key] for r in sel if r["f_Hz"] == f]) for f in fs])
        near = np.array([np.mean([r["groups"][0]["ratio"] for r in sel if r["f_Hz"] == f]) for f in fs])
        fa = np.array(fs)
        use = fa >= 500
        ax.plot(fa[use], mean("ratio_formula")[use], color=WCOLOUR[wn], lw=2 * PT, marker="o", ms=4.5, mec=SURFACE, mew=1.2, label=f"{lab}: all strands")
        ax.plot(fa[use], near[use], color=WCOLOUR[wn], lw=1.4 * PT, marker="^", ms=4.5, mec=SURFACE, mew=1.0, label=f"{lab}: the strands within 1 mm of the opening side")
        ax.plot(fa[use], mean("ratio_exact")[use], color=WCOLOUR[wn], lw=1.4 * PT, ls="--", label=f"{lab}: with the lone wire's factors")
    ax.axhline(1.0, color=AXIS, lw=1.2 * PT)
    ax.set_xscale("log")
    _legend(ax, "lower left")
    ax = fig.add_subplot(1, 2, 2)
    _style(ax, "formula [mW a strand]", "FEMM excess loss [mW a strand]")
    ax.set_position([0.58, 0.24, 0.32, 0.63])
    tab = tables[("W1", "same", 1, 20000.0)]
    exc = (tab["p_femm"] - tab["p_dc"]) * 1e3
    form = tab["p_formula"] * 1e3
    top = float(tab["xy"][:, 1].min())
    depth = tab["xy"][:, 1] - top
    sc = ax.scatter(form, exc, c=depth, cmap=CMAP, s=(5 * PT) ** 2 * 2, edgecolors="none", norm=Normalize(0, float(depth.max())))
    lim = [min(form.min(), exc.min()) * 0.8, max(form.max(), exc.max()) * 1.2]
    ax.plot(lim, lim, color=AXIS, lw=1.2 * PT)
    ax.set_xscale("log")
    ax.set_yscale("log")
    _colorbar(fig, ax, sc, label="depth below the coil's top edge [mm]", gap=0.012)
    _save(fig, path)


def _strand_polys(field: dict, ac: dict) -> tuple[np.ndarray, np.ndarray]:
    """Triangle vertices (m, 3, 2) and |J| over the DC density of the strand's own current, per triangle (the mean of the three nodal values)."""
    strand = field["strand"]
    j_dc = np.abs(ac["current_A"][strand]) / ac["area_m2"][strand]
    rel = np.abs(field["J_nodes"]).mean(axis=1) / j_dc
    return field["tri_xy"], rel


def fig_maps(res: dict, path: Path) -> None:
    from optimization import ac_slot_study as study
    cases = list(study.MAP_CASES)
    fig = _figure(13.0, 8.8, "Current density inside the strands, and the loss of each strand, at the rated current",
                  "Top row: |J| over the DC current density of the same strand, in a window around the slot opening (the coil's top edge is at the window's lower third), per triangle of the FEMM mesh. "
                  "Bottom row: every strand of the slot, coloured by its loss over its own DC loss. 118.6 ampere-turns a coil side, layers in phase, copper at 100 degC. The three cases: W1 at "
                  "20 kHz (d / delta = 0.59), W2 at 20 kHz (1.18) and W2 at 100 kHz (2.6).")
    for k, (wn, f) in enumerate(cases):
        ac = res[f"map_{wn}_{f:g}"]
        d_mm = study.WINDINGS[wn]["d"]
        verts, rel = _strand_polys(ac["field"], ac)
        x0, y0 = 0.0, float(ac["xy"][:, 1].min())
        ax = fig.add_axes([0.055 + 0.318 * k, 0.57, 0.24, 0.32])
        ax.set_anchor("W")
        ax.set_facecolor(SURFACE)
        sel = (np.abs(verts[:, :, 0]).max(axis=1) < 1.5) & (verts[:, :, 1].max(axis=1) < y0 + 2.6)
        pc = PolyCollection(verts[sel], array=rel[sel], cmap=CMAP, norm=Normalize(float(np.percentile(rel, 1)), float(np.percentile(rel, 99.5))), edgecolors="none")
        ax.add_collection(pc)
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(y0 - 0.15, y0 + 2.6)
        ax.set_aspect("equal")
        ax.set_title(f"{wn} ({d_mm} mm), {f / 1e3:g} kHz: |J| / J_dc", fontsize=9.5, color=INK, loc="left", fontweight="bold")
        ax.set_xlabel("x across the slot [mm]", fontsize=9, color=INK2)
        ax.set_ylabel("y, radial [mm]", fontsize=9, color=INK2)
        ax.tick_params(colors=INK2, labelsize=8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(AXIS)
        _colorbar(fig, ax, pc)
        # the second row: the slot with every strand coloured by its loss ratio
        ax2 = fig.add_axes([0.055 + 0.318 * k, 0.085, 0.24, 0.40])
        ax2.set_anchor("W")
        ax2.set_facecolor(SURFACE)
        p_dc = study.ac_slot.dc_loss_W(ac, study.layers_for(study.WINDINGS[wn], "same"))
        ratio = ac["loss_W"] / p_dc
        patches = PolyCollection([], array=[])
        norm = LogNorm(vmin=max(1.0, float(ratio.min())), vmax=float(ratio.max()))
        circles = [Circle((x, y), d_mm / 2) for x, y in ac["xy"]]
        from matplotlib.collections import PatchCollection
        pcs = PatchCollection(circles, cmap=CMAP, norm=norm, edgecolors="none")
        pcs.set_array(ratio)
        ax2.add_collection(pcs)
        ax2.set_xlim(-3.4, 3.4)
        ax2.set_ylim(float(ac["xy"][:, 1].min()) - 0.4, float(ac["xy"][:, 1].max()) + 0.4)
        ax2.set_aspect("equal")
        ax2.set_title(f"{wn}, {f / 1e3:g} kHz: strand loss / its DC loss", fontsize=9.5, color=INK, loc="left", fontweight="bold")
        ax2.set_xlabel("x across the slot [mm]", fontsize=9, color=INK2)
        ax2.set_ylabel("y, radial [mm]", fontsize=9, color=INK2)
        ax2.tick_params(colors=INK2, labelsize=8)
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax2.spines[side].set_color(AXIS)
        _colorbar(fig, ax2, pcs, ticks=_plain_log_ticks(norm.vmin, norm.vmax))
    _save(fig, path)


def fig_circ(rows: list[dict], path: Path) -> None:
    fig = _figure(11.0, 5.2, "What connecting the strands of a turn in parallel costs: R_ac / R_dc of the slot with the strands in parallel circuits",
                  "R_ac / R_dc - 1 of the slot (the extra loss over the DC loss) when the strands of a turn are connected in parallel (one circuit a turn, the turn's current shared as the strands' positions in "
                  "the field make it), against the same strands with equal sharing (dashed gray): the difference is the circulating-current loss. Both windings at the rated current, layers in phase, copper at "
                  "100 degC. Placement of the strands in the turns: line = every turn spans the slot's depth (the worst), random, cluster = compact blobs, row = side by side at one depth. "
                  "Table view: circ.json.")
    cols = {"line": ORANGE, "random": BLUE, "cluster": AQUA, "row": INK2}
    marks = {"line": "o", "random": "s", "cluster": "^", "row": "D"}
    for k, wn in enumerate(("W1", "W2")):
        ax = fig.add_subplot(1, 2, k + 1)
        _style(ax, "frequency [Hz]", "R_ac / R_dc - 1" if k == 0 else "")
        ax.set_position([0.07 + 0.50 * k, 0.26, 0.41, 0.62])
        for arr in ("line", "random", "cluster", "row"):
            sel = sorted([r for r in rows if r["winding"] == wn and r["arrangement"] == arr], key=lambda r: r["f_Hz"])
            ax.plot([r["f_Hz"] for r in sel], [r["R_parallel"] - 1.0 for r in sel], color=cols[arr], lw=2 * PT, marker=marks[arr], ms=5, mec=SURFACE, mew=1.2, label=arr)
        sel = sorted([r for r in rows if r["winding"] == wn and r["arrangement"] == "line"], key=lambda r: r["f_Hz"])
        ax.plot([r["f_Hz"] for r in sel], [r["R_series"] - 1.0 for r in sel], color=MUTED, lw=1.5 * PT, ls="--", label="equal sharing (no circulating currents)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title({"W1": "W1: 26 turns x 8 strands of 0.315 mm", "W2": "W2: 13 turns x 4 strands of 0.630 mm"}[wn], fontsize=9.5, color=INK, loc="left", fontweight="bold")
        _legend(ax, "lower right")
    _save(fig, path)
