"""learn/proximity/02_two_wires.py: two neighbouring wires: a current in one makes a field, the field changing in time drives eddy currents in the other, and the eddy currents pull the
current of both toward one side.

The idea (the explanation for the proximity-effect note, from the field of one wire)
    1. A wire carrying a current I makes, at a distance s, the field B = mu0 I / (2 pi s) (Ampere's law: the circle of radius s encloses I; B points round the wire, across a neighbour that
       lies beside it).
    2. With I alternating at frequency f, that field alternates and the neighbour sits in a changing field. A wire in a field B(t) across it, carrying no net current of its own, loses
       power (the proximity effect, 01_one_wire.py): at low frequency pi sigma w^2 B^2 d^4 / 128 per metre, with B the PEAK field at the neighbour. So the neighbour of a live wire heats up
       although it carries nothing. Case (a) of this script checks that chain: the field of wire 1 at wire 2 against FEMM's loss in wire 2.
    3. If the neighbour also carries current, the two fields add. For currents in the SAME direction the fields cancel between the wires and add outside, and the induced electric field pushes each
       wire's current toward the side FACING AWAY from the other; for OPPOSITE currents (a go-and-return pair) the fields add between the wires and the current of each is pulled toward the side facing
       the other. Either way the current density is no longer even and the loss rises, the more the closer the wires are: R_ac / R_dc against the spacing (the right panel).
    4. To first order the rise is the sum of the two passive-wire losses of case (a): each wire of the pair sits in the field I / (2 pi s) of the other, whichever way the other's current flows, and the
       eddy-current pattern it answers with (odd across the wire) does not mix, in the loss, with the even pattern of its own current. The right panel checks it: the dotted line.

What the script does
    A 0.63 mm copper wire pair (100 degC, d / delta = 2 at 57.9 kHz: thick enough for the crowding to be seen) in FEMM's time-harmonic solver through the project's solver slot.
      (a) wire 1 carries 1 A (peak), wire 2 a series circuit of 0 A: the eddy-current loss in wire 2 against the formula with B = mu0 I / (2 pi s), for centre distances s = 1.1 d ... 8 d;
      (b) both carry 1 A in the same direction, and (c) in opposite directions: |J| over the DC density for s = 1.3 d (the pictures) and R_ac / R_dc of the pair for s = 1.1 d ... 8 d.
    Writes out/learn/proximity/two_wires.png and prints the table.

Run (from the Python to FEMM folder)
    uv run python learn/proximity/02_two_wires.py            (about 4 s with the solver slot free, 24 solves)

Numbers to expect (FEMM 5.0.20261101, mesh d / 14 in the wires, 2 degree arcs, the box 60 mm in radius; the run prints them)
    (a) FEMM's loss in the passive wire over the formula, s / d = 1.1, 1.2, 1.5, 2, 3, 5, 8:  1.0018, 0.9911, 0.9785, 0.9661, 0.9569, 0.9732, 0.9911 (33.8 uW at 1.1 d, 0.63 uW at 8 d, for 20 mm of wire).
        So the formula is right to a few per cent and the loss falls as 1 / s^2 (the slope between 2 d and 8 d is 1.98). Do not read more into the third digit: a finer mesh (d / 28, 1 degree arcs)
        and a box four times larger (the numbers 'mesh_over_d', 'arc_deg' and 'radius_mm' in pair_task) move single points by up to 2 per cent (0.95 to 0.99 for s = 1.5 d ... 5 d over four such
        variants, tried once, no trend), and the passive wire's loss is a small number made from a field that the mesh has to carry across the gap.
    (b), (c) R_ac / R_dc of the pair, s / d = 1.1 ... 8:  same direction 1.0673 ... 1.0213, opposite directions 1.0668 ... 1.0212, against the lone wire's 1.0205 (the table has all seven spacings).
        The first-order prediction (lone wire + the passive wire's loss of (a) over a wire's DC loss) is 1.0669 ... 1.0214: within 0.09 per cent of DC of both at every spacing, and the mean of the two to 1.5e-4.

Self-check (the answers follow).
  1. In case (a) the loss in the passive wire falls with the spacing as 1 / s^2. Why that exponent?
  2. The same-direction and the opposite-direction pair differ at s = 2 d by 1.03463 against 1.03316. Which loses more, and why is the difference so small?
Answers.
  1. The field at the neighbour is mu0 I / (2 pi s), which falls as 1 / s, and the proximity loss goes as the square of the field: 1 / s^2.
  2. The same-direction pair, by 0.15 per cent of its DC loss (about a tenth of the extra loss over the lone wire). To first order the two are the same: each wire sits in the field I / (2 pi s) of its
     neighbour whichever way the neighbour's current flows, and the eddy currents that field drives (odd across the wire) do not mix with the wire's own current distribution (even), so the pair's extra
     loss is twice the passive wire's loss of (a), the same in (b) and (c). What is left is second order (how each wire's eddy currents disturb the other's current); this script does not separate it.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import prox_common as pc
from optimization import ac_theory

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import NullFormatter

D_MM = 0.63
D_OVER_DELTA = 2.0
SPACINGS = (1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0)         # centre distance over d
PICTURE_S = 1.3


def pair_task(tag: str, s_over_d: float, i1: complex, i2: complex) -> dict:
    s = s_over_d * D_MM
    return {"tag": tag, "f": pc.frequency_for(D_OVER_DELTA, D_MM), "radius_mm": 60.0, "B": None, "mesh_over_d": 1 / 14, "arc_deg": 2.0, "fem_dir": str(pc.OUT / "fem"),
            "wires": [{"x": -s / 2, "y": 0.0, "d": D_MM, "circuit": "w1"}, {"x": s / 2, "y": 0.0, "d": D_MM, "circuit": "w2"}],
            "circuits": {"w1": (i1, "series"), "w2": (i2, "series")}}


def main() -> int:
    pc.OUT.mkdir(parents=True, exist_ok=True)
    tasks = []
    for s in SPACINGS:
        tasks += [pair_task(f"two_a_{s:g}", s, 1.0, 0.0), pair_task(f"two_b_{s:g}", s, 1.0, 1.0), pair_task(f"two_c_{s:g}", s, 1.0, -1.0)]
    tasks += [pair_task(f"two_a_{PICTURE_S:g}", PICTURE_S, 1.0, 0.0), pair_task(f"two_b_{PICTURE_S:g}", PICTURE_S, 1.0, 1.0), pair_task(f"two_c_{PICTURE_S:g}", PICTURE_S, 1.0, -1.0)]
    results = pc.solve(tasks, "learn/proximity: two wires", workers=4)
    by = {t["tag"]: r for t, r in zip(tasks, results)}
    f = pc.frequency_for(D_OVER_DELTA, D_MM)
    rho, d, L = pc.RHO, D_MM * 1e-3, 20e-3
    mu0 = ac_theory.MU0
    print(f"{D_MM} mm wires, d / delta = {D_OVER_DELTA}, f = {f:.0f} Hz, copper at 100 degC; currents 1 A peak")
    print("(a) the passive wire: FEMM's loss against the formula with B = mu0 I / (2 pi s) at its centre")
    print("    s / d    FEMM loss [uW]   formula [uW]   FEMM / formula")
    rows_a = []
    for s in SPACINGS:
        r = by[f"two_a_{s:g}"]
        femm_p = r["loss_W"][1]
        b = mu0 * 1.0 / (2 * math.pi * s * D_MM * 1e-3)
        form = float(ac_theory.proximity_exact(d, f, b, rho)) * L
        rows_a.append((s, femm_p, form))
        print(f"  {s:6.1f}   {1e6 * femm_p:14.5f}   {1e6 * form:12.5f}   {femm_p / form:10.4f}")
    lone = float(ac_theory.skin_ratio(d, f, rho))
    print("(b), (c): R_ac / R_dc of the pair (both wires); the first-order prediction is the lone wire's R_ac / R_dc plus the passive wire's loss of (a) over a wire's DC loss")
    print("    s / d    same direction   opposite directions   first-order prediction   (lone wire alone: %.5f)" % lone)
    rows_bc = []
    for s in SPACINGS:
        vals = []
        for case in ("b", "c"):
            r = by[f"two_{case}_{s:g}"]
            r_dc = rho * L / np.array(r["area_m2"])
            vals.append(sum(r["loss_W"]) / sum(0.5 * 1.0 ** 2 * r_dc))
        ra = by[f"two_a_{s:g}"]
        pred = lone + ra["loss_W"][1] / (0.5 * 1.0 ** 2 * rho * L / ra["area_m2"][1])        # each of the two wires adds the passive wire's loss; both over the pair's DC loss
        rows_bc.append((s, *vals, pred))
        print(f"  {s:6.1f}   {vals[0]:12.5f}   {vals[1]:18.5f}   {pred:18.5f}")
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    cmap = LinearSegmentedColormap.from_list("sb", [pc.SURFACE, "#9ec4ef", pc.BLUE, "#0b2a52"])
    fig = plt.figure(figsize=(12.0, 8.4), facecolor=pc.SURFACE)
    fig.text(0.012, 0.985, "Two neighbouring wires: a current in one, eddy currents in the other, and the crowding when both carry current", ha="left", va="top", fontsize=11,
             color=pc.INK, fontweight="bold")
    fig.text(0.012, 0.012, f"Two {D_MM} mm copper wires at 100 degC, {f / 1e3:.1f} kHz (d / delta = {D_OVER_DELTA}), centres {PICTURE_S} d apart in the pictures; FEMM 5.0 time-harmonic, 1 A peak a wire. "
             "Pictures: |J| over the DC density of 1 A in one wire. Source: learn/proximity/02_two_wires.py.", ha="left", va="bottom", fontsize=7.6, color=pc.INK2, wrap=True)
    titles = {"a": "(a) wire 1 carries 1 A, wire 2 nothing", "b": "(b) both carry 1 A, same direction", "c": "(c) 1 A and -1 A (go and return)"}
    for k, case in enumerate(("a", "b", "c")):
        r = by[f"two_{case}_{PICTURE_S:g}"]
        j_dc = 1.0 / np.array(r["area_m2"])[0]
        rel = np.abs(r["J"]).mean(axis=1) / j_dc
        ax = fig.add_subplot(2, 3, k + 1)
        coll = pc.plot_wires(ax, r, rel, cmap, Normalize(0, float(rel.max())), lim=None)
        s = PICTURE_S * D_MM
        ax.set_xlim(-s / 2 - D_MM * 0.62, s / 2 + D_MM * 0.62)
        ax.set_ylim(-D_MM * 0.62, D_MM * 0.62)
        ax.set_title(titles[case], fontsize=9.5, color=pc.INK, loc="left", fontweight="bold")
        ax.set_xlabel("x [mm]", fontsize=9, color=pc.INK2)
        cb = fig.colorbar(coll, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=7.5, colors=pc.INK2)
    ax = fig.add_subplot(2, 2, 3)
    pc.style_axes(ax, "centre distance over wire diameter  s / d", "loss in the passive wire [uW]  (case a)")
    ax.loglog([r[0] for r in rows_a], [1e6 * r[1] for r in rows_a], "o-", color=pc.BLUE, lw=2.0, ms=5, mec=pc.SURFACE, label="FEMM")
    ax.loglog([r[0] for r in rows_a], [1e6 * r[2] for r in rows_a], "--", color=pc.ORANGE, lw=2.0, label="formula with B = mu0 I / (2 pi s)")
    leg = ax.legend(frameon=False, fontsize=8.2)
    for t in leg.get_texts():
        t.set_color(pc.INK2)
    ax = fig.add_subplot(2, 2, 4)
    pc.style_axes(ax, "centre distance over wire diameter  s / d", "R_ac / R_dc of the pair")
    ax.semilogx([r[0] for r in rows_bc], [r[1] for r in rows_bc], "o-", color=pc.BLUE, lw=2.0, ms=5, mec=pc.SURFACE, label="same direction (b)")
    ax.semilogx([r[0] for r in rows_bc], [r[2] for r in rows_bc], "s-", color=pc.ORANGE, lw=2.0, ms=5, mec=pc.SURFACE, label="opposite directions (c)")
    ax.semilogx([r[0] for r in rows_bc], [r[3] for r in rows_bc], ":", color=pc.INK2, lw=2.0, label="first order: lone wire + the passive wire's loss of (a)")
    ax.axhline(lone, color=pc.GRAY, lw=1.2, ls="--", label="a lone wire")
    leg = ax.legend(frameon=False, fontsize=8.2)
    for t in leg.get_texts():
        t.set_color(pc.INK2)
    for ax in fig.axes[-2:]:
        ax.set_xticks([1.1, 2, 3, 5, 8])
        ax.set_xticklabels(["1.1", "2", "3", "5", "8"])
        ax.xaxis.set_minor_formatter(NullFormatter())
    fig.subplots_adjust(left=0.07, right=0.955, top=0.93, bottom=0.10, wspace=0.28, hspace=0.34)
    path = pc.OUT / "two_wires.png"
    fig.savefig(path, dpi=130, facecolor=pc.SURFACE)
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
