"""optimization.maps_figures: the figures of phase 1 (matplotlib PNGs, no FEMM): the flux maps, the torque map with the MTPA path, cross-saturation, the inductances, the iron loss.

What it is, and why
    maps_report.py builds flux_map.npz and iron_map.npz and says whether they can be trusted; this file draws what they say. Each picture has its numbers in a CSV beside it
    (fig_<name>.csv: the "table view" that a picture alone cannot give). Written to be read on their own: every panel has a title that says what it shows, and a note under the
    figure says what is drawn and what to look at.

The pictures (written to out\\optimization\\maps\\chosen\\figures\\; 130 dpi)
    fig_flux_maps.png        psi_d and psi_q over the current plane (the motoring half, i_q >= 0; the table holds both halves), with the MTPA path and the circles of 1, 2 and 3 times
                             the rated current. psi_d falls through zero as the d current grows: the magnets' flux is cancelled; psi_q grows with i_q and flattens: saturation.
    fig_torque_map.png       the torque with the MTPA path of the maps and the path of the constant-parameter formula (the rated point's constants): the paths part above the rated current.
    fig_cross_saturation.png psi_d against i_d for several i_q (the fan is the cross-saturation: the q current lowers the d flux), and the cross terms L_dq and L_qd against i_q
                             (two derivatives of two different maps that must be equal: Maxwell's reciprocity, drawn as a line and as marks that fall on it).
    fig_inductances.png      the incremental L_dd and L_qq over the plane, and the apparent and incremental inductances along the MTPA path.
    fig_iron_loss.png        the iron loss of teeth, yoke and rotor against the current along the MTPA path at three speeds (3,000 rpm, 8,000 rpm and the top speed 16,330 rpm),
                             for the CAD's 0.3 mm lamination, with the datasheet's 0.356 mm total as the thin dashed line.
    fig_iron_fit.png         the Bertotti fit against the AK Steel M-19 datasheet (loss against frequency at four flux densities) and its error, with the cross-validation folds.
    fig_positions.png        how many rotor positions the mean needs (the error of a K-position mean against the 60-position one) and the ripple at 15 and 30 positions.
    fig_mirror.png           how far the mirror symmetry i_q -> -i_q holds across the plane (the reason the negative half was solved, not mirrored).

The look (the data-visualization method, the palette of optimization/winding_figures.py)
    surface #fcfcfb, ink #0b0b0b and #52514e, hairline grid #e1e0d9, axes #c3c2b7; series blue #2a78d6, orange #eb6834, aqua #1baf7a (aqua always with a second channel: its own marker
    shape or line style, because it is 2.7:1 on the surface), context gray #898781. Maps use one sequential blue ramp or a blue-orange diverging ramp with a neutral midpoint at zero;
    lines are 2 px, markers at least 8 px with a 2 px surface ring; text is never in a series colour; one axis per panel (no dual axes: quantities of different units get their own panel,
    or are made dimensionless against their declared limit).

Run (from the Python to FEMM folder; after maps_report's flux, cross and loss sections; seconds)
    uv run python -m optimization.maps_figures            all the figures
    uv run python -m optimization.maps_figures flux_maps iron_loss      some of them

Self-check (the answers follow).
  1. In fig_cross_saturation the lines of psi_d against i_d fan out as i_q grows. What would a motor without cross-saturation show?
  2. Why is the iron loss drawn in three panels, one per speed, and not on one set of axes?
Answers.
  1. The same curve for every i_q: psi_d would depend on i_d only, and the fan would collapse into one line.
  2. The loss grows with the speed to a power between 1 and 2, so the 16,330 rpm loss is about 20 times the 3,000 rpm one; on one axis the slow-speed curves would be flat lines
     on the floor. Three panels, each with its own scale and its speed in its title, tell the truth about each.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

from optimization import iron_loss, maps, maps_run

OUT = maps_run.OUT
FIG = OUT / "figures"
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
CONTEXT = "#898781"
BLUE_RAMP = ["#9ec5f2", "#5d9be6", "#2a78d6", "#14478a"]            # light to dark: an ordered series (i_q = 0, 1, 2, 3 x rated)
PT = 0.75
FONT = ["Segoe UI", "DejaVu Sans"]
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f1f6fd", "#b9d3f2", "#6aa3e8", "#2a78d6", "#154c92"])
DIV = LinearSegmentedColormap.from_list("div_orange_blue", ["#c4491a", "#eb6834", "#f6cdb8", "#ecebe6", "#bcd6f3", "#2a78d6", "#154c92"])
RPM_LIST = (3000.0, 8000.0, 16330.0)
PP = 4
TOP = 0.88                                                           # the panels' top edge (figure fraction)


def _style(ax, xlabel: str, ylabel: str, grid: bool = True) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3, color=AXIS)
    if grid:
        ax.grid(True, color=GRID, linewidth=1 * PT, linestyle="-")
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)


def _figure(width: float, height: float, title: str, note: str):
    """A figure with its title (top left) and its note (bottom left, wrapped); fig._bottom is the bottom edge the panels must keep so that the note and the axis labels fit."""
    plt.rcParams["font.family"] = FONT
    fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=11, color=INK, fontweight="bold")
    fig.text(0.012, 0.012, note, ha="left", va="bottom", fontsize=7.6, color=INK2, wrap=True)
    n_lines = math.ceil(len(note) / (width * 15.0))                                   # about 15 characters of 7.6 pt text per inch of width
    fig._bottom = min(0.40, (0.012 * height + n_lines * 0.135 + 0.62) / height)       # the note's lines (0.135 in each) and 0.62 in for tick labels and the axis label
    return fig


def _panel(fig, k: int, n: int, left: float = 0.055, right: float = 0.985, gap: float = 0.07, cbar: float = 0.0) -> list[float]:
    """The position of panel k of n in a row, above the figure's reserved bottom margin; cbar reserves room at each panel's right for a colour bar (figure fraction)."""
    w = (right - left - gap * (n - 1)) / n
    return [left + k * (w + gap), fig._bottom, w - cbar, TOP - fig._bottom]


def _cbar(fig, ax, mappable, label: str | None = None, width: float = 0.009, pad: float = 0.008, fmt: str | None = None):
    ax.apply_aspect()                                              # an equal-aspect axes is smaller than the box it was given: put the bar beside the axes as drawn
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + pad, pos.y0, width, pos.height])
    cb = fig.colorbar(mappable, cax=cax)
    if fmt:
        cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter(fmt))
    cb.ax.tick_params(labelsize=7.5, colors=INK2)
    cb.outline.set_edgecolor(AXIS)
    if label:
        cb.set_label(label, color=INK2, fontsize=8.5)
    return cb


def _save(fig, name: str, rows: list[dict] | None = None) -> Path:
    FIG.mkdir(parents=True, exist_ok=True)
    path = FIG / f"fig_{name}.png"
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    if rows:
        with open(FIG / f"fig_{name}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), restval="")       # rows of different panels may have different columns
            w.writeheader()
            w.writerows(rows)
    return path


def _dot(color, marker="o", size=6.0, label=None) -> Line2D:
    return Line2D([], [], marker=marker, linestyle="", markersize=size, markerfacecolor=color, markeredgecolor=SURFACE, markeredgewidth=2 * PT, label=label)


def _line(color, label, style="-", marker=None) -> Line2D:
    return Line2D([], [], color=color, linewidth=2 * PT, linestyle=style, marker=marker, markersize=6, markeredgecolor=SURFACE, label=label)


def _legend(ax, handles, loc="best", ncol=1) -> None:
    leg = ax.legend(handles=handles, loc=loc, frameon=False, fontsize=8.3, handletextpad=0.5, borderaxespad=0.2, ncol=ncol)
    for t in leg.get_texts():
        t.set_color(INK2)


def _circles(ax, I_R: float, mults=(1, 2, 3), label=True) -> None:
    g = np.radians(np.linspace(0, 90, 120))
    for m in mults:
        ax.plot(-m * I_R * np.sin(g), m * I_R * np.cos(g), color=CONTEXT, linewidth=1 * PT, linestyle=(0, (4, 3)), zorder=3)
        if label:
            ax.text(-m * I_R * math.sin(math.radians(88)), m * I_R * math.cos(math.radians(88)) + 0.7, f"{m}x", color=INK2, fontsize=7.5, ha="right", va="bottom")


def _plane(fm: maps.FluxMaps, nx: int = 161, ny: int = 121):
    ids = np.linspace(fm.range_d[0], fm.range_d[1], nx)
    iqs = np.linspace(0.0, fm.range_q, ny)
    A, B = np.meshgrid(ids, iqs, indexing="ij")
    return ids, iqs, A, B


# ============================================================================================== the pictures
def fig_flux_maps(fm: maps.FluxMaps) -> Path:
    I_R = fm.rated
    ids, iqs, A, B = _plane(fm)
    psi_d, psi_q = fm.psi_dq(A, B)
    path = fm.mtpa(np.linspace(0.02, 3.0, 60) * I_R)
    zero = ids[np.argmin(np.abs(psi_d[:, 0]))]
    fig = _figure(13.4, 5.6, "The flux linkages of the optimized motor over the current plane (motoring half, i_q >= 0)",
                  f"FEMM 5.0, fine mesh, 2D 20 mm, the 10-turn placeholder winding (scale to another winding with s = turns / (10 paths): psi x s, currents / s). Currents are peak, in the "
                  f"amplitude-invariant d-q frame; d is the magnet axis. Each point is the mean over 30 rotor positions of a ripple period. Dashed arcs: constant current, 1, 2 and 3 times the "
                  f"rated {I_R:.2f} A; orange: the MTPA path of the maps. Left, psi_d crosses zero near i_d = {zero:.0f} A on the d axis (the current cancels the magnets' "
                  f"{1e3 * fm.psi_pm0:.1f} mWb); right, psi_q rises with i_q and flattens as the steel saturates. Table view: fig_flux_maps.csv.")
    for k, (Z, title, cmap) in enumerate(((1e3 * psi_d, "psi_d [mWb]", DIV), (1e3 * psi_q, "psi_q [mWb]", SEQ))):
        ax = fig.add_subplot(1, 2, k + 1)
        _style(ax, "i_d [A]", "i_q [A]", grid=False)
        ax.set_position(_panel(fig, k, 2, left=0.055, right=0.975, gap=0.08, cbar=0.03))
        levels = np.arange(np.floor(Z.min() / 5) * 5, np.ceil(Z.max() / 5) * 5 + 0.1, 5.0)          # 5 mWb steps; 0 is a level of psi_d
        norm = TwoSlopeNorm(vcenter=0.0, vmin=levels[0], vmax=levels[-1]) if k == 0 else None       # the diverging ramp's neutral midpoint is zero flux
        cs = ax.contourf(A, B, Z, levels=levels, cmap=cmap, norm=norm)
        lines = ax.contour(A, B, Z, levels=levels[::2], colors=INK2, linewidths=0.5)
        ax.clabel(lines, fmt="%.0f", fontsize=7, inline=True)
        ax.plot(path[1], path[2], color=ORANGE, linewidth=2.0, zorder=4)
        _circles(ax, I_R)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(fm.range_d[0], fm.range_d[1])
        ax.set_ylim(0, fm.range_q)
        ax.set_title(title, loc="left", fontsize=9.5, color=INK, pad=4)
        _cbar(fig, ax, cs, "mWb")
    rows = [{"i_d_A": float(a), "i_q_A": float(b), "psi_d_mWb": float(1e3 * pd_), "psi_q_mWb": float(1e3 * pq_)}
            for a, b, pd_, pq_ in zip(A[::10, ::10].ravel(), B[::10, ::10].ravel(), psi_d[::10, ::10].ravel(), psi_q[::10, ::10].ravel())]
    return _save(fig, "flux_maps", rows)


def fig_torque_map(fm: maps.FluxMaps) -> Path:
    I_R = fm.rated
    ids, iqs, A, B = _plane(fm)
    T = fm.torque(A, B)
    cur = np.linspace(0.02, 3.0, 60) * I_R
    g_dq, d_dq, q_dq, T_dq = fm.mtpa(cur)
    # the constant-parameter formula with the constants of the rated point
    g1, d1, q1, T1 = fm.mtpa(I_R)
    L1 = fm.inductances(d1, q1)
    psi_pm, Ld, Lq = fm.psi_pm0, L1["L_d_app"], L1["L_q_app"]
    dL = Lq - Ld
    s = (-psi_pm + np.sqrt(psi_pm ** 2 + 8 * dL ** 2 * cur ** 2)) / (4 * dL * cur)
    g_f = np.degrees(np.arcsin(np.minimum(s, 1.0)))
    d_f, q_f = -cur * np.sin(np.radians(g_f)), cur * np.cos(np.radians(g_f))
    fig = _figure(8.6, 7.4, "Torque of the optimized motor over the current plane, and where the best angle goes (MTPA)",
                  "Filled contours: the mean torque of FEMM's air-gap integral (N.m). Orange: the MTPA path of the flux maps (the angle of the largest torque at each current). Aqua with triangles: the "
                  f"constant-parameter formula's path with the constants of the rated point (psi_PM {1e3 * psi_pm:.2f} mWb, L_d {1e3 * Ld:.3f} mH, L_q {1e3 * Lq:.3f} mH): it follows the maps "
                  "up to about the rated current and then bends away from them. Dashed arcs: 1, 2 and 3 times the rated current. Table view: fig_torque_map.csv.")
    ax = fig.add_subplot(1, 1, 1)
    _style(ax, "i_d [A]", "i_q [A]", grid=False)
    pos = _panel(fig, 0, 1, left=0.085, right=0.90)
    ax.set_position(pos)
    cs = ax.contourf(A, B, T, levels=np.arange(np.floor(T.min() * 2) / 2, np.ceil(T.max()) + 0.01, 0.5), cmap=SEQ)
    lines = ax.contour(A, B, T, levels=np.arange(0, np.ceil(T.max()) + 0.01, 1.0), colors=INK2, linewidths=0.5)
    ax.clabel(lines, fmt="%.0f", fontsize=7, inline=True)
    _circles(ax, I_R)
    ax.plot(d_dq, q_dq, color=ORANGE, linewidth=2.0, zorder=5)
    ax.plot(d_f, q_f, color=AQUA, linewidth=2.0, zorder=5, linestyle="-")
    ax.plot(d_f[::6], q_f[::6], linestyle="", marker="^", markersize=7, markerfacecolor=AQUA, markeredgecolor=SURFACE, markeredgewidth=2 * PT, zorder=6)
    for m in (1, 2, 3):
        gm, dm, qm, Tm = fm.mtpa(m * I_R)
        ax.plot([dm], [qm], marker="o", markersize=8, markerfacecolor=ORANGE, markeredgecolor=SURFACE, markeredgewidth=2 * PT, zorder=7)
        ax.annotate(f"{m}x: {gm:.1f} deg, {Tm:.2f} N.m", (dm, qm), xytext=(-9, 7), textcoords="offset points", fontsize=7.8, color=INK, ha="right")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(fm.range_d[0], fm.range_d[1])
    ax.set_ylim(0, fm.range_q)
    _legend(ax, [_line(ORANGE, "MTPA path of the maps", marker="o"), _line(AQUA, "constant-parameter formula (rated constants)", marker="^")], loc="lower right")
    _cbar(fig, ax, cs, "N.m")
    rows = [{"I_A": float(c), "gamma_maps_deg": float(g), "T_maps_Nm": float(t), "gamma_formula_deg": float(gf)} for c, g, t, gf in zip(cur, g_dq, T_dq, g_f)]
    return _save(fig, "torque_map", rows)


def fig_cross_saturation(fm: maps.FluxMaps) -> Path:
    I_R = fm.rated
    fig = _figure(12.6, 5.6, "Cross-saturation: the q current lowers the d flux, and the two cross terms are one number",
                  f"Left: psi_d against i_d at four values of i_q (FEMM maps, bicubic spline). With no cross-saturation the four lines would be one. The dashed gray line is the no-load flux "
                  f"{1e3 * fm.psi_pm0:.1f} mWb; at i_d = 0 the q current alone takes it down to {1e3 * float(fm.psi_dq(0.0, 3 * I_R)[0]):.1f} mWb at 3 times the rated current. Right: the cross terms L_dq = "
                  "d psi_d / d i_q (lines) and L_qd = d psi_q / d i_d (marks) against i_q at three d currents: two derivatives of two different maps that Maxwell's reciprocity says "
                  "are equal; the marks fall on the lines. Table view: fig_cross_saturation.csv.")
    ax = fig.add_subplot(1, 2, 1)
    _style(ax, "i_d [A]", "psi_d [mWb]")
    ax.set_position(_panel(fig, 0, 2, left=0.065, right=0.975, gap=0.09))
    x = np.linspace(fm.range_d[0], fm.range_d[1], 200)
    rows = []
    for m, col in zip((0, 1, 2, 3), BLUE_RAMP):
        y = 1e3 * fm.psi_dq(x, np.full_like(x, m * I_R))[0]
        ax.plot(x, y, color=col, linewidth=2 * PT)
        ax.annotate(f"i_q = {m} x rated" if m else "i_q = 0", (x[-1], y[-1]), xytext=(4, 2 if m != 3 else -7), textcoords="offset points", fontsize=8, color=INK2, va="center")
        rows += [{"panel": "psi_d", "i_q_x_rated": m, "i_d_A": float(a), "psi_d_mWb": float(b)} for a, b in zip(x[::10], y[::10])]
    ax.axhline(1e3 * fm.psi_pm0, color=CONTEXT, linewidth=1 * PT, linestyle=(0, (4, 3)))
    ax.axvline(0, color=AXIS, linewidth=1 * PT)
    ax.set_xlim(fm.range_d[0], fm.range_d[1] + 12)
    ax = fig.add_subplot(1, 2, 2)
    _style(ax, "i_q [A]", "cross inductance [mH]")
    ax.set_position(_panel(fig, 1, 2, left=0.065, right=0.975, gap=0.09))
    yq = np.linspace(0, fm.range_q, 120)
    for m, col, mk in zip((0, -1, -2), (BLUE, ORANGE, AQUA), ("o", "s", "^")):
        i_d = m * I_R
        L = fm.inductances(np.full_like(yq, i_d), yq)
        ax.plot(yq, 1e3 * L["L_dq"], color=col, linewidth=2 * PT)
        sel = slice(5, None, 12)
        ax.plot(yq[sel], 1e3 * L["L_qd"][sel], linestyle="", marker=mk, markersize=7, markerfacecolor=SURFACE, markeredgecolor=col, markeredgewidth=1.6)
        rows += [{"panel": "cross", "i_d_x_rated": m, "i_q_A": float(a), "L_dq_mH": float(1e3 * b), "L_qd_mH": float(1e3 * c)} for a, b, c in zip(yq[::10], L["L_dq"][::10], L["L_qd"][::10])]
    ax.axhline(0, color=AXIS, linewidth=1 * PT)
    _legend(ax, [_line(BLUE, "i_d = 0"), _line(ORANGE, "i_d = -1 x rated"), _line(AQUA, "i_d = -2 x rated"),
                 Line2D([], [], marker="o", linestyle="", markerfacecolor=SURFACE, markeredgecolor=INK2, markersize=7, label="marks: L_qd (colour of its line)")], loc="upper right")
    return _save(fig, "cross_saturation", rows)


def fig_inductances(fm: maps.FluxMaps) -> Path:
    I_R = fm.rated
    ids, iqs, A, B = _plane(fm, 121, 91)
    L = fm.inductances(A, B)
    cur = np.linspace(0.1, 3.0, 60) * I_R
    g, d, q, T = fm.mtpa(cur)
    Lp = fm.inductances(d, q)
    k1 = int(np.argmin(np.abs(cur - I_R)))
    fig = _figure(13.4, 5.6, "The inductances of the optimized motor: where they fall, and the apparent against the incremental",
                  "Left and middle: the incremental inductances L_dd = d psi_d / d i_d and L_qq = d psi_q / d i_q over the current plane (mH), with the MTPA path in orange. Right: along the MTPA path, "
                  "the apparent inductances L_d,app = (psi_d - psi_PM0) / i_d and L_q,app = psi_q / i_q (solid) and the incremental ones (dashed): the incremental inductance, the one in "
                  "v = R i + L di/dt, falls faster than the apparent one because it is the slope of a curve that is bending over. d is the magnet axis: L_d < L_q (saliency "
                  f"{float(Lp['L_q_app'][k1] / Lp['L_d_app'][k1]):.2f} at the rated current). Table view: fig_inductances.csv.")
    xs = (0.050, 0.325, 0.640)
    for k, (Z, title) in enumerate(((1e3 * L["L_dd"], "L_dd, incremental [mH]"), (1e3 * L["L_qq"], "L_qq, incremental [mH]"))):
        ax = fig.add_subplot(1, 3, k + 1)
        _style(ax, "i_d [A]", "i_q [A]", grid=False)
        ax.set_position([xs[k], fig._bottom, 0.19, TOP - fig._bottom])
        lv = np.linspace(np.floor(Z.min() * 10) / 10, np.ceil(Z.max() * 10) / 10, 13)
        cs = ax.contourf(A, B, Z, levels=lv, cmap=SEQ)
        lines = ax.contour(A, B, Z, levels=lv[::2], colors=INK2, linewidths=0.5)
        ax.clabel(lines, fmt="%.1f", fontsize=7, inline=True)
        ax.plot(d, q, color=ORANGE, linewidth=2.0, zorder=4)
        _circles(ax, I_R, label=False)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(fm.range_d[0], fm.range_d[1])
        ax.set_ylim(0, fm.range_q)
        ax.set_title(title, loc="left", fontsize=9.5, color=INK, pad=4)
        _cbar(fig, ax, cs, fmt="%.1f")
    ax = fig.add_subplot(1, 3, 3)
    _style(ax, "current magnitude [x rated]", "inductance [mH]")
    ax.set_position([xs[2], fig._bottom, 0.985 - xs[2], TOP - fig._bottom])
    x = cur / I_R
    series = (("L_d_app", BLUE, "-", "o"), ("L_dd", BLUE, (0, (5, 3)), None), ("L_q_app", ORANGE, "-", "s"), ("L_qq", ORANGE, (0, (5, 3)), None))
    for key, col, style, mk in series:
        ax.plot(x, 1e3 * Lp[key], color=col, linewidth=2 * PT, linestyle=style)
        if mk:
            ax.plot(x[::6], 1e3 * Lp[key][::6], linestyle="", marker=mk, markersize=7, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=2 * PT)
    rows = [{"I_x_rated": float(a), "gamma_deg": float(b), **{k: float(1e3 * Lp[k][i]) for k in ("L_d_app", "L_q_app", "L_dd", "L_qq", "L_dq", "L_qd")}}
            for i, (a, b) in enumerate(zip(x, g)) if i % 3 == 0]
    ax.set_xlim(0, 3.05)
    ax.set_ylim(0, None)
    _legend(ax, [_line(BLUE, "L_d apparent", marker="o"), _line(BLUE, "L_d incremental", style=(0, (5, 3))), _line(ORANGE, "L_q apparent", marker="s"),
                 _line(ORANGE, "L_q incremental", style=(0, (5, 3)))], loc="lower left")
    return _save(fig, "inductances", rows)


def fig_iron_loss(fm: maps.FluxMaps, imap: iron_loss.IronLossMap) -> Path:
    I_R = fm.rated
    cur = np.linspace(0.0, 3.0, 61) * I_R
    g, d, q, T = fm.mtpa(np.maximum(cur, 1e-6))
    with open(iron_loss.FIT_JSON, encoding="utf-8") as f:
        rms = json.load(f)["fit"]["rms"]
    fig = _figure(13.4, 5.6, "Iron loss of the optimized motor along the MTPA path, at three speeds",
                  "Whole machine, 20 mm stack, 10-turn placeholder winding (the iron loss depends on the ampere-turns, so another winding only rescales the current axis: s times the current). Solid, 2 px: "
                  "the CAD's 0.3 mm lamination (the headline: the classical term scaled by (0.3 / 0.356)^2 = 0.710, hysteresis and excess left as fitted); thin dashed gray: the total at the datasheet's "
                  f"0.356 mm. Bertotti model fitted to AK Steel M-19 29 gauge (rms error {100 * rms:.1f} % on the 112 datasheet points), flux density "
                  "waveforms of every steel triangle from FEMM over one electrical period. Each panel has its own scale: the loss grows with the speed as f^1 ... f^2. Table view: fig_iron_loss.csv.")
    rows = []
    for k, rpm in enumerate(RPM_LIST):
        ax = fig.add_subplot(1, 3, k + 1)
        _style(ax, "current magnitude [x rated]", "iron loss [W]")
        ax.set_position(_panel(fig, k, 3, left=0.055, right=0.985, gap=0.055))
        r = imap(d, q, rpm, lamination_mm=0.3)
        r356 = imap(d, q, rpm, lamination_mm=0.356)
        x = cur / I_R
        for name, col, mk in (("teeth", BLUE, "o"), ("yoke", ORANGE, "s"), ("rotor", AQUA, "^")):
            ax.plot(x, r[name], color=col, linewidth=2 * PT)
            ax.plot(x[::10], r[name][::10], linestyle="", marker=mk, markersize=7, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=2 * PT)
        ax.plot(x, r["total"], color=INK, linewidth=2 * PT)
        ax.plot(x, r356["total"], color=CONTEXT, linewidth=1 * PT, linestyle=(0, (4, 3)))
        ax.set_xlim(0, 3.0)
        ax.set_ylim(0, None)
        ax.set_title(f"{rpm:,.0f} rpm  (f = {PP * rpm / 60:.0f} Hz)", loc="left", fontsize=9.5, color=INK, pad=4)
        if k == 0:
            _legend(ax, [_line(BLUE, "teeth", marker="o"), _line(ORANGE, "yoke", marker="s"), _line(AQUA, "rotor", marker="^"), _line(INK, "total, 0.3 mm"),
                         _line(CONTEXT, "total, 0.356 mm", style=(0, (4, 3)))], loc="upper left")
        rows += [{"rpm": rpm, "I_x_rated": float(a), "teeth_W": float(b), "yoke_W": float(c), "rotor_W": float(e), "total_0p3mm_W": float(f_), "total_0p356mm_W": float(h)}
                 for a, b, c, e, f_, h in zip(x[::5], r["teeth"][::5], r["yoke"][::5], r["rotor"][::5], r["total"][::5], r356["total"][::5])]
    return _save(fig, "iron_loss", rows)


def fig_iron_fit() -> Path:
    with open(iron_loss.FIT_JSON, encoding="utf-8") as fh:
        js = json.load(fh)
    fit = js["fit"]
    f, B, P = iron_loss.datasheet_points()
    fig = _figure(13.4, 5.8, "The iron-loss model against the M-19 datasheet: the fit and its cross-validation",
                  f"AK Steel Di-Max M-19, 29 gauge (0.356 mm), typical core loss in W/kg (EMERF Lamination Steels, MIT 6.685 excerpt, p. 5). Left: loss against frequency at four flux densities, "
                  f"dots the datasheet, lines the model k_h f B^a + k_e f^2 B^2 + k_x f^1.5 B^1.5 with k_h = {fit['kh']:.4f}, a = {fit['alpha']:.3f}, k_x = {fit['kx']:.3e} fitted and "
                  f"k_e = {fit['ke']:.3e} computed from the conductivity and the thickness (no fit). Middle and right: the relative error of the model at the 112 points (gray dots: the points it was fitted on) "
                  "and the cross-validation: refit with a whole flux density (middle, blue) or a whole frequency (right, orange) left out; the mark is the rms error at the left-out points, the bar runs to the "
                  f"largest; a circle is an interior fold (judged), a cross an edge fold (the model extrapolates: shown, not judged). Declared limits for the interior folds: rms <= {100 * iron_loss.TOL_RMS:.0f} %, "
                  f"largest <= {100 * iron_loss.TOL_MAX:.0f} % (dashed gray). Table view: fig_iron_fit.csv.")
    ax = fig.add_subplot(1, 3, 1)
    _style(ax, "frequency [Hz]", "loss [W/kg]")
    ax.set_position(_panel(fig, 0, 3, left=0.055, right=0.985, gap=0.06))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ff = np.geomspace(40, 2400, 80)
    for b, col, mk in ((0.4, BLUE_RAMP[0], "o"), (1.0, BLUE_RAMP[1], "s"), (1.5, BLUE_RAMP[2], "^"), (1.7, BLUE_RAMP[3], "D")):
        sel = np.isclose(B, b)
        ax.plot(ff, iron_loss.predict(fit, ff, b), color=col, linewidth=2 * PT)
        ax.plot(f[sel], P[sel], linestyle="", marker=mk, markersize=6, markerfacecolor=INK, markeredgecolor=SURFACE, markeredgewidth=1 * PT)
        ax.annotate(f"{b} T", (f[sel].max(), P[sel][np.argmax(f[sel])]), xytext=(6, 0), textcoords="offset points", fontsize=8, color=INK2, va="center")
    ax.set_xlim(40, 3200)
    err = iron_loss.predict(fit, f, B) / P - 1
    for k, (kind, xdata, col, xlabel, logx) in enumerate((("B", B, BLUE, "flux density B [T]", False), ("f", f, ORANGE, "frequency [Hz]", True))):
        ax = fig.add_subplot(1, 3, k + 2)
        _style(ax, xlabel, "relative error of the model [%]")
        ax.set_position(_panel(fig, k + 1, 3, left=0.055, right=0.985, gap=0.06))
        if logx:
            ax.set_xscale("log")
        ax.plot(xdata, 100 * err, linestyle="", marker="o", markersize=4, markerfacecolor=CONTEXT, markeredgecolor=SURFACE, markeredgewidth=0.5, alpha=0.7)
        for row in js["cross_validation"]["folds"]:
            if row["fold"] != kind:
                continue
            x = row["left_out"]
            ax.plot([x, x], [-100 * row["max"], 100 * row["max"]], color=col, linewidth=1.5 * PT, alpha=0.8)
            ax.plot([x], [100 * row["rms"]], marker="o" if row["interior"] else "X", markersize=7, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=2 * PT, linestyle="")
        ax.axhline(0, color=AXIS, linewidth=1 * PT)
        for lim in (iron_loss.TOL_RMS, -iron_loss.TOL_RMS, iron_loss.TOL_MAX, -iron_loss.TOL_MAX):
            ax.axhline(100 * lim, color=CONTEXT, linewidth=1 * PT, linestyle=(0, (4, 3)))
        ax.set_ylim(-26, 26)
    rows = [{"f_Hz": float(a), "B_T": float(b), "datasheet_W_per_kg": float(c), "model_W_per_kg": float(c * (1 + e)), "relative_error": float(e)} for a, b, c, e in zip(f, B, P, err)]
    return _save(fig, "iron_fit", rows)


def fig_positions() -> Path:
    with open(OUT / "report_flux.json", encoding="utf-8") as fh:
        rep = json.load(fh)
    from optimization import maps_report
    table = rep["positions"]["table"]
    Ks = sorted({r["K"] for r in table})
    I_R = maps_run.Context().I_R
    fig = _figure(13.4, 5.8, "How many rotor positions does the mean need? (60 positions, 1 degree apart, as the reference)",
                  "At eight points from no load to the corner of the grid (i_d, i_q) = (-3, 3) x rated, the mean of K equally spaced positions per 60 degrees against the 60-position mean. Left: the "
                  "worst error of psi (relative to |psi|), against its declared limit of 1e-3. Middle: the same for the torque, against its declared limit (2e-3 of T + 2 mN.m, drawn as 1 for every "
                  "point: the error is plotted as a fraction of the limit that belongs to its own point). Right: the Fourier-rebuilt peak-to-peak ripple at K = 15 (the project's standard) and K = 30 "
                  "(the maps') against the 60-position value at each point: K = 15 reads the rated point's ripple 2.4 % low, because the waveform has harmonics beyond the 7th of the ripple period. "
                  "Six positions give the means to 2e-4 of |psi|; the maps use 30. Table view: fig_positions.csv.")
    ax = fig.add_subplot(1, 3, 1)
    _style(ax, "positions per 60 degrees, K", "worst psi error / its limit")
    ax.set_position(_panel(fig, 0, 3, left=0.055, right=0.985, gap=0.065))
    ax.set_yscale("log")
    ax.set_xscale("log")
    e_psi = [max(r["psi_err"] for r in table if r["K"] == K) / maps_report.K_TOL_PSI for K in Ks]
    ax.plot(Ks, e_psi, color=BLUE, linewidth=2 * PT, marker="o", markersize=7, markerfacecolor=BLUE, markeredgecolor=SURFACE, markeredgewidth=2 * PT)
    ax.axhline(1.0, color=CONTEXT, linewidth=1 * PT, linestyle=(0, (4, 3)))
    ax.set_xticks(Ks)
    ax.set_xticklabels([str(k) for k in Ks])
    ax.annotate("declared limit", (Ks[-1], 1.0), xytext=(-4, 4), textcoords="offset points", fontsize=8, color=INK2, ha="right")
    ax = fig.add_subplot(1, 3, 2)
    _style(ax, "positions per 60 degrees, K", "worst T error / its limit")
    ax.set_position(_panel(fig, 1, 3, left=0.055, right=0.985, gap=0.065))
    ax.set_yscale("log")
    ax.set_xscale("log")
    e_T = [max(r["T_err_Nm"] / (maps_report.K_TOL_T * abs(r["T"]) + maps_report.K_TOL_T_ABS) for r in table if r["K"] == K) for K in Ks]
    ax.plot(Ks, e_T, color=ORANGE, linewidth=2 * PT, marker="s", markersize=7, markerfacecolor=ORANGE, markeredgecolor=SURFACE, markeredgewidth=2 * PT)
    ax.axhline(1.0, color=CONTEXT, linewidth=1 * PT, linestyle=(0, (4, 3)))
    ax.set_xticks(Ks)
    ax.set_xticklabels([str(k) for k in Ks])
    ax.annotate("declared limit", (Ks[-1], 1.0), xytext=(-4, 4), textcoords="offset points", fontsize=8, color=INK2, ha="right")
    ax = fig.add_subplot(1, 3, 3)
    _style(ax, "point (i_d, i_q) in rated currents", "torque ripple peak-to-peak [mN.m]")
    ax.set_position(_panel(fig, 2, 3, left=0.055, right=0.985, gap=0.065))
    pts = sorted({r["point"] for r in table})
    xs = np.arange(len(pts))

    def get(K, key):
        return [1e3 * next(r[key] for r in table if r["K"] == K and r["point"] == p) for p in pts]

    ax.plot(xs, get(15, "ripple60"), linewidth=0, marker="D", markersize=9, markerfacecolor=SURFACE, markeredgecolor=INK, markeredgewidth=1.6)
    ax.plot(xs - 0.14, get(15, "ripple_K"), linestyle="", marker="o", markersize=6, markerfacecolor=BLUE, markeredgecolor=SURFACE, markeredgewidth=2 * PT)
    ax.plot(xs + 0.14, get(30, "ripple_K"), linestyle="", marker="s", markersize=6, markerfacecolor=ORANGE, markeredgecolor=SURFACE, markeredgewidth=2 * PT)
    labels = []
    for p in pts:
        r = next(r for r in table if r["point"] == p and r["K"] == 15)
        labels.append(f"({r['i_d'] / I_R:.1f}, {r['i_q'] / I_R:.1f})")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=7.5, rotation=30)
    _legend(ax, [Line2D([], [], marker="D", linestyle="", markerfacecolor=SURFACE, markeredgecolor=INK, markersize=8, label="60 positions (reference)"),
                 _dot(BLUE, "o", label="15 positions (the project's standard)"), _dot(ORANGE, "s", label="30 positions (the maps')")], loc="upper left")
    rows = [{"point": r["point"], "i_d_A": r["i_d"], "i_q_A": r["i_q"], "K": r["K"], "psi_err": r["psi_err"], "T_err_Nm": r["T_err_Nm"], "ripple60_Nm": r.get("ripple60", ""),
             "ripple_K_Nm": r.get("ripple_K", "")} for r in table]
    return _save(fig, "positions", rows)


def fig_mirror() -> Path:
    from optimization import maps_report
    built = np.load(OUT / "flux_map.npz")
    ids, iqs = built["id_A"], built["iq_A"]
    I_R = float(built["rated_current_A"])
    jp = [j for j, v in enumerate(iqs) if v > 1e-9]
    jm = [int(np.argmin(np.abs(iqs + iqs[j]))) for j in jp]
    mag = (np.hypot(ids[:, None], iqs[jp][None, :]) / I_R).ravel()
    psi = np.hypot(built["psi_d"][:, jp], built["psi_q"][:, jp])
    e_q = (100 * np.abs(built["psi_q"][:, jm] + built["psi_q"][:, jp]) / psi).ravel()
    e_d = (100 * np.abs(built["psi_d"][:, jm] - built["psi_d"][:, jp]) / psi).ravel()
    T_sum = (built["T_age"][:, jm] + built["T_age"][:, jp]).ravel()
    T_rel = (100 * np.abs(T_sum) / np.maximum(np.abs(built["T_age"][:, jp]).ravel(), 0.05))
    med = lambda y, lo, hi: float(np.median(y[(mag > lo) & (mag <= hi)]))
    fig = _figure(13.4, 5.4, "How far the mirror symmetry i_q -> -i_q holds: the reason the negative half of the maps was solved, not mirrored",
                  "Each dot is one pair of grid points (i_d, +i_q) and (i_d, -i_q), both solved by FEMM (30 and 15 positions). If the rotor and winding were mirror-symmetric about the d axis, "
                  "psi_d(-i_q) = psi_d(i_q), psi_q(-i_q) = -psi_q(i_q) and T(-i_q) = -T(i_q) would hold at every point. The medians are "
                  f"{med(e_q, 0, 1):.2f} % (psi_q) and {med(T_rel, 0, 1):.2f} % (T) up to the rated current, {med(e_q, 2, 3):.2f} % and {med(T_rel, 2, 3):.2f} % at 2 ... 3 times it; single points reach 2 % in psi_q and "
                  "3 % in T (the 36-slot, 8-pole winding is not mirror-symmetric about the d axis, and saturation makes it show). Dashed gray: the limits declared for the check (0.5 % of |psi|, 1 % of T). "
                  "Table view: fig_mirror.csv.")
    for k, (y, ylab, lim) in enumerate(((e_q, "psi_q asymmetry [% of |psi|]", maps_report.TOL_MIRROR_PSI * 100), (e_d, "psi_d asymmetry [% of |psi|]", maps_report.TOL_MIRROR_PSI * 100),
                                          (T_rel, "torque asymmetry [% of T]", maps_report.TOL_MIRROR_T_REL * 100))):
        ax = fig.add_subplot(1, 3, k + 1)
        _style(ax, "current magnitude [x rated]", ylab)
        ax.set_position(_panel(fig, k, 3, left=0.055, right=0.985, gap=0.065))
        ax.plot(mag, y, linestyle="", marker="o", markersize=5, markerfacecolor=BLUE, markeredgecolor=SURFACE, markeredgewidth=0.6, alpha=0.8)
        ax.axhline(lim, color=CONTEXT, linewidth=1 * PT, linestyle=(0, (4, 3)))
        ax.set_xlim(0, 4.3)
        ax.set_ylim(0, None)
    rows = [{"current_x_rated": float(a), "psi_q_asym_pct": float(b), "psi_d_asym_pct": float(c), "T_sum_Nm": float(d), "T_asym_pct": float(e)}
            for a, b, c, d, e in zip(mag[::3], e_q[::3], e_d[::3], T_sum[::3], T_rel[::3])]
    return _save(fig, "mirror", rows)


FIGURES = ("flux_maps", "torque_map", "cross_saturation", "inductances", "iron_loss", "iron_fit", "positions", "mirror")


def main(argv: list[str] | None = None) -> int:
    which = (sys.argv[1:] if argv is None else argv) or list(FIGURES)
    fm = maps.FluxMaps.load()
    imap = iron_loss.IronLossMap() if iron_loss.IRON_MAP.exists() else None
    made = []
    for name in which:
        if name == "flux_maps":
            made.append(fig_flux_maps(fm))
        elif name == "torque_map":
            made.append(fig_torque_map(fm))
        elif name == "cross_saturation":
            made.append(fig_cross_saturation(fm))
        elif name == "inductances":
            made.append(fig_inductances(fm))
        elif name == "iron_loss":
            made.append(fig_iron_loss(fm, imap))
        elif name == "iron_fit":
            made.append(fig_iron_fit())
        elif name == "positions":
            made.append(fig_positions())
        elif name == "mirror":
            made.append(fig_mirror())
        else:
            raise SystemExit(f"unknown figure {name!r}: {', '.join(FIGURES)}")
    for p in made:
        print(f"wrote {p} ({p.stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
