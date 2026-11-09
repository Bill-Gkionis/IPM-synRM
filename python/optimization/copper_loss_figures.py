"""optimization.copper_loss_figures: the pictures of the copper-loss answer, drawn from the fast model (optimization/copper_loss.py). No FEMM: a few seconds.

What it is, and why
    optimization/copper_loss_study.py prints the tables that answer 'does the AC loss of the copper matter in this motor, where, and for which winding'. These are the same numbers as
    pictures, for the report and for looking at. Every figure has its title (the case), a note under it saying what is drawn and with which assumptions, units on its axes, and the
    numbers of its table in the study's printout. The look follows optimization/ac_slot_figures.py: surface #fcfcfb, ink #0b0b0b and #52514e, grid #e1e0d9, axes #c3c2b7; series blue
    #2a78d6, orange #eb6834, aqua #1baf7a; context gray #898781.

The pictures (written into out/optimization/copper_ac/final/)
    fig_final_speed.png  the parts of the AC loss over the winding's DC loss, 1,000 to 16,330 rpm, for W1 (26 turns x 8 strands of 0.315 mm, 4 paths) and W2 (13 x 4 x 0.630 mm, 2 paths): the
                         machine field's proximity loss, the skin effect, the PWM ripple, and what the circulating currents add for strands in parallel placed at random or side by side. Rated
                         ampere-turns at the MTPA angle at every speed (the inverter's voltage limit is marked: beyond it the real machine runs in field weakening).
    fig_final_map.png    the proximity + skin loss over the DC loss for the strand diameter and the speed at a constant copper fill, at the rated and at twice the rated ampere-turns.
    fig_final_pwm.png    the PWM ripple current and its loss for the six windings of the winding study's inverter table, at 10 and 20 kHz on 48 V and 20 kHz on 60 V.
    fig_final_circ.png   the circulating-current loss of strands in parallel against the speed for the four placements of the strands in the turns, for W1 and W2.

Run
    uv run python -m optimization.copper_loss_figures        (about 2 s; also called by optimization.copper_loss_study)

Self-check (the answers follow).
  1. The loss in the map rises towards the top right (thick strands, high speed) as d^2 times the speed squared. Where would the rise stop, and why?
  2. The ripple loss of a winding in the PWM figure is highest for the windings with the lowest inductance. Why does a higher bus voltage not help?
Answers.
  1. When d / delta is no longer small the strand shields itself (ac_theory.proximity_ratio falls below 1: 0.93 at d / delta 1.8, 0.78 at 2.5) and the iron screens the field of the high harmonics; the
     model carries the first (exact lone-wire factors) and not the second (see copper_loss.py: p_prox is an upper bound).
  2. The ripple is the volt-second error of the inverter, which is proportional to the bus voltage: V_dc / (L f_sw). A higher bus makes the ripple larger, not smaller (compare the 48 V and 60 V bars);
     only a larger inductance (more series turns) or a higher switching frequency lowers it.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: F401  (selects the backend before the figure helpers are used)
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter

import paths
from optimization import ac_slot_figures as fg
from optimization import copper_loss as cl

OUT = paths.path("out") / "optimization" / "copper_ac" / "final"
AT_RATED = 118.6
STACK_MM = 20.0
FILL_MAP = 0.448
SPEED_MAX = 16330.0
W1 = {"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367}
W2 = {"N": 13, "a": 2, "n_s": 4, "d_mm": 0.630, "d_over_mm": 0.7042}
INVERTER_TABLE = [                                      # (current limit [A], N, a, n_s, d_mm) from the winding study's section 9.5
    (10, 20, 2, 13, 0.280), (15, 27, 4, 19, 0.200), (20, 20, 4, 13, 0.280), (30, 14, 4, 29, 0.224), (50, 8, 4, 26, 0.315), (100, 4, 4, 13, 0.630)]
MACHINE = paths.path("out") / "optimization" / "copper_ac" / "machine"


def mtpa_angle(at: float) -> float:
    """The MTPA current angle [deg] at `at` ampere-turns a coil side: interpolated over the FEMM points of optimization/ac_machine.py (12.33, 23.02, 37.28 and 43.37 deg at 0.5, 1, 2 and 3 x the
    rated 118.6 ampere-turns)."""
    import json
    pts = [q for q in json.loads((MACHINE / "points.json").read_text(encoding="utf-8")) if q["group"] == "mtpa"]
    return float(np.interp(at, [q["at"] for q in pts], [q["gamma"] for q in pts]))


def d_over_for(d_mm: float) -> float:
    """The overall (enamelled) diameter [mm] of a strand of bare diameter d_mm: interpolated over the winding study's wire catalogue points."""
    return float(np.interp(d_mm, [0.2, 0.28, 0.315, 0.63], [0.239, 0.329, 0.367, 0.7042]))


def i_for(w: dict, at: float = AT_RATED) -> float:
    """The peak phase current [A] that gives `at` ampere-turns a coil side to the winding: N i / a = at."""
    return at * w["a"] / w["N"]


def run(w: dict, speeds, gamma: float, at: float = AT_RATED, f_sw: float = 20000.0, v_dc: float = 48.0, arrangement: str = "transposed") -> dict:
    """The model on an array of speeds [rpm] for the winding w at `at` ampere-turns and the current angle gamma: the arrays of copper_loss.losses_arrays."""
    return cl.losses_arrays(w["N"], w["a"], w["n_s"], w["d_mm"], w["d_over_mm"], i_for(w, at), gamma, np.asarray(speeds, float), f_sw, v_dc, 100.0, STACK_MM, arrangement=arrangement)


def corner_rpm(w: dict, gamma: float, at: float = AT_RATED, v_dc: float = 48.0) -> float | None:
    """The lowest speed [rpm] (steps of 20) at which the winding's voltage need at the MTPA point exceeds what the bus can make (modulation index 1), or None below 16,330 rpm."""
    sp = np.arange(1000.0, SPEED_MAX + 1.0, 20.0)
    over = run(w, sp, gamma, at, 20000.0, v_dc)["overmodulated"]
    hit = np.flatnonzero(over)
    return float(sp[hit[0]]) if len(hit) else None


def fig_speed(path: Path) -> dict:
    gamma = mtpa_angle(AT_RATED)
    speeds = np.linspace(1000.0, SPEED_MAX, 80)
    fig = fg._figure(12.0, 5.8, "AC loss of the copper over the winding's DC loss, 1,000 to 16,330 rpm: W1 with 0.315 mm strands and W2 with 0.63 mm strands",
                     f"118.6 ampere-turns a coil side at the MTPA angle ({gamma:.1f} deg) at every speed, copper at 100 degC, 20 mm stack, 48 V bus, 20 kHz. Proximity: the machine's own field (rotor and "
                     "armature harmonics, FEMM magnetostatic solves of the sector; an upper bound, the iron's screening of the high harmonics is not in it). Skin: the fundamental current in its own strand. PWM: the "
                     "ripple current's rms (space-vector PWM) times the slot's AC resistance at the carrier harmonics. Circulating: strands of a turn in PARALLEL circuits, placed at random or side by side "
                     "(a transposed winding has none). The dotted vertical line marks where the winding's voltage need exceeds the 48 V bus at this current and angle (field weakening beyond). "
                     "Source: optimization/copper_loss_figures.py (copper_loss.py).")
    out = {}
    for k, (title, w) in enumerate((("W1: 26 turns x 8 strands of 0.315 mm, 4 paths", W1), ("W2: 13 turns x 4 strands of 0.630 mm, 2 paths", W2))):
        ax = fig.add_axes([0.065 + 0.49 * k, 0.24, 0.42, 0.64])
        fg._style(ax, "speed [rpm]", "loss over the DC loss of the winding [%]" if k == 0 else "")
        base = run(w, speeds, gamma)
        dc = base["p_dc_slot"] + base["p_dc_end"]
        series = [("proximity, the machine's field", 100 * base["p_prox"] / dc, fg.BLUE, "-"),
                  ("skin effect of the fundamental", 100 * base["p_skin"] / dc, fg.MUTED, ":"),
                  ("PWM ripple, 20 kHz", 100 * base["p_pwm"] / dc, fg.AQUA, "-")]
        for arr, col, ls in (("random", fg.ORANGE, "-"), ("row", fg.ORANGE, "--")):
            r = run(w, speeds, gamma, arrangement=arr)
            series.append((f"circulating, strands {'at random' if arr == 'random' else 'in a row'} (parallel)", 100 * r["p_circ"] / dc, col, ls))
        for label, y, col, ls in series:
            ax.semilogy(speeds, np.maximum(y, 1e-12), color=col, lw=2 * fg.PT, ls=ls, label=label)
        corner = corner_rpm(w, gamma)
        if corner is not None:
            ax.axvline(corner, color=fg.INK2, lw=1.2 * fg.PT, ls=":")
            ax.text(corner + 150, 1.5e-7 * 4, "voltage limit", rotation=90, fontsize=7.8, color=fg.INK2, va="bottom")
        ax.set_ylim(1e-7, 100.0)
        ax.set_xlim(speeds[0], speeds[-1])
        ax.set_title(title, fontsize=9.5, color=fg.INK, loc="left", fontweight="bold")
        if k == 0:
            fg._legend(ax, "lower right")
        out[title] = {"speeds": speeds, "series": {s[0]: s[1] for s in series}, "corner_rpm": corner}
    fg._save(fig, path)
    return out


def fig_map(path: Path) -> dict:
    ds = np.geomspace(0.2, 1.4, 56)
    speeds = np.linspace(1000.0, SPEED_MAX, 56)
    fig = fg._figure(12.0, 5.4, "Proximity + skin loss over the DC loss for the strand diameter and the speed, at a constant copper fill of 0.448",
                     "13 turns a coil side and 2 paths, the number of strands in hand following the diameter (continuous) so that the copper area is the same, 20 mm stack, copper at 100 degC, the MTPA angle at each "
                     "current. Left: the rated 118.6 ampere-turns a coil side; right: twice that. Contour lines at 0.1, 1 and 5 per cent. The machine's field is that of the FEMM points (an upper bound at the high "
                     "harmonics); no PWM, no circulating currents (a transposed winding). Source: optimization/copper_loss_figures.py (copper_loss.py).")
    out = {"d_mm": ds, "rpm": speeds}
    n_s = FILL_MAP * 36.218 / (13 * math.pi * ds ** 2 / 4.0)
    norm = LogNorm(vmin=1e-3, vmax=30.0)
    for k, fac in enumerate((1.0, 2.0)):
        at = fac * AT_RATED
        r = cl.losses_arrays(13, 2, n_s[:, None], ds[:, None], ds[:, None] * 1.12, at * 2 / 13, mtpa_angle(at), speeds[None, :], 0.0, 48.0, 100.0, STACK_MM)
        share = 100.0 * (r["p_prox"] + r["p_skin"]) / (r["p_dc_slot"] + r["p_dc_end"])
        out[f"share_{fac:g}x"] = share
        ax = fig.add_axes([0.065 + 0.49 * k, 0.24, 0.38, 0.64])
        fg._style(ax, "speed [rpm]", "strand diameter [mm]" if k == 0 else "")
        ax.grid(False)
        mesh = ax.pcolormesh(speeds, ds, np.maximum(share, 1e-3), norm=norm, cmap=fg.CMAP, shading="auto")
        cs = ax.contour(speeds, ds, share, levels=[0.1, 1.0, 5.0], colors=[fg.INK], linewidths=0.9)
        ax.clabel(cs, fmt=lambda v: f"{v:g} %", fontsize=7.8)
        ax.set_yscale("log")
        ax.set_yticks([0.2, 0.28, 0.4, 0.63, 1.0, 1.4])
        ax.set_yticklabels(["0.2", "0.28", "0.4", "0.63", "1.0", "1.4"])
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_title(f"{fac:g} x rated ({at:.0f} ampere-turns a coil side, MTPA angle {mtpa_angle(at):.1f} deg)", fontsize=9.5, color=fg.INK, loc="left", fontweight="bold")
        cb = fig.colorbar(mesh, ax=ax, fraction=0.05, pad=0.02)
        cb.ax.tick_params(labelsize=7.8, colors=fg.INK2)
        cb.set_label("[% of the DC loss]", color=fg.INK2, fontsize=8.5)
    fg._save(fig, path)
    return out


def fig_pwm(path: Path) -> list[dict]:
    gamma = mtpa_angle(AT_RATED)
    cases = (("10 kHz, 48 V", 10000.0, 48.0, fg.BLUE), ("20 kHz, 48 V", 20000.0, 48.0, fg.ORANGE), ("20 kHz, 60 V", 20000.0, 60.0, fg.AQUA))
    fig = fg._figure(12.0, 5.4, "PWM ripple of the six windings of the inverter table: the ripple current and what it costs, at 3,000 rpm and the rated ampere-turns",
                     "Each winding is the winding study's choice for its inverter limit (turns a coil side, parallel paths, strands in hand x strand diameter), at 118.6 ampere-turns a coil side, MTPA angle, "
                     "3,000 rpm, copper at 100 degC, 20 mm stack, strands transposed. The ripple's rms is that of space-vector PWM with the d-axis inductance of the 2-D model scaled with the series turns squared "
                     "(an upper bound); its loss is that rms squared times the slot's AC resistance at the carrier's harmonics, over the winding's DC loss. Source: optimization/copper_loss_figures.py.")
    rows = []
    labels = []
    for lim, n, a, n_s, d in INVERTER_TABLE:
        w = {"N": n, "a": a, "n_s": n_s, "d_mm": d, "d_over_mm": d_over_for(d)}
        labels.append(f"{lim} A\n{n} t, {a} p\n{n_s} x {d:g}")
        cells = {}
        for name, f_sw, v, _ in cases:
            r = run(w, 3000.0, gamma, AT_RATED, f_sw, v)
            cells[name] = (float(r["ripple_rms_a"]), 100.0 * float(r["p_pwm"]) / float(r["p_dc_slot"] + r["p_dc_end"]), float(r["l_phase_h"]))
        rows.append({"limit": lim, "winding": w, "i_peak": i_for(w), "cells": cells})
    x = np.arange(len(rows))
    for j, (ylabel, idx, ylog) in enumerate((("ripple current, rms [A]", 0, False), ("ripple loss over the winding's DC loss [%]", 1, True))):
        ax = fig.add_axes([0.065 + 0.49 * j, 0.27, 0.42, 0.61])
        fg._style(ax, "inverter current limit; turns a coil side, paths; strands x diameter [mm]", ylabel)
        for i, (name, _, _, col) in enumerate(cases):
            ax.bar(x + (i - 1) * 0.26, [r["cells"][name][idx] for r in rows], width=0.24, color=col, label=name)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7.8)
        if ylog:
            ax.set_yscale("log")
        if j == 0:
            fg._legend(ax, "upper left")
    fg._save(fig, path)
    return rows


def fig_circ(path: Path) -> dict:
    gamma = mtpa_angle(AT_RATED)
    speeds = np.linspace(1000.0, SPEED_MAX, 80)
    fig = fg._figure(12.0, 5.6, "What strands in parallel cost: the circulating-current loss over the DC loss for four placements of the strands of a turn",
                     f"118.6 ampere-turns a coil side at the MTPA angle ({gamma:.1f} deg), copper at 100 degC, 20 mm stack. The strands of one turn are connected in parallel; line = every turn spans the depth of the "
                     "slot (the worst), random, cluster = compact blobs, row = side by side at one depth. From the FEMM slot model with its circulating-current law (optimization/ac_slot_study.py, circ stage), "
                     "the slot's own field and the machine's A spread; the strands' ends are assumed to join at the terminals, where the whole turn is in parallel. A transposed winding has none. "
                     "Source: optimization/copper_loss_figures.py (copper_loss.py).")
    cols = {"line": fg.ORANGE, "random": fg.BLUE, "cluster": fg.AQUA, "row": fg.INK2}
    out = {}
    for k, (title, w) in enumerate((("W1: 26 turns x 8 strands of 0.315 mm, 4 paths", W1), ("W2: 13 turns x 4 strands of 0.630 mm, 2 paths", W2))):
        ax = fig.add_axes([0.065 + 0.49 * k, 0.24, 0.42, 0.64])
        fg._style(ax, "speed [rpm]", "circulating-current loss over the DC loss [%]" if k == 0 else "")
        for arr in ("line", "random", "cluster", "row"):
            r = run(w, speeds, gamma, arrangement=arr)
            y = 100 * r["p_circ"] / (r["p_dc_slot"] + r["p_dc_end"])
            ax.semilogy(speeds, np.maximum(y, 1e-12), color=cols[arr], lw=2 * fg.PT, label=arr)
            out[(title, arr)] = (speeds, y)
        ax.set_ylim(1e-4, 100.0)
        ax.set_xlim(speeds[0], speeds[-1])
        ax.set_title(title, fontsize=9.5, color=fg.INK, loc="left", fontweight="bold")
        if k == 0:
            fg._legend(ax, "upper left")
    fg._save(fig, path)
    return out


def main() -> int:
    import time
    OUT.mkdir(parents=True, exist_ok=True)
    cl.tables()
    t0 = time.perf_counter()
    for name, fn in (("fig_final_speed.png", fig_speed), ("fig_final_map.png", fig_map), ("fig_final_pwm.png", fig_pwm), ("fig_final_circ.png", fig_circ)):
        fn(OUT / name)
        print(f"written {OUT / name}")
    print(f"done in {time.perf_counter() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
