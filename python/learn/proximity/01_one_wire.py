"""learn/proximity/01_one_wire.py: the current density inside ONE round wire as the frequency rises: the skin effect (its own current) and the proximity effect (a field from outside).

The idea, from the smallest piece (the explanation for the proximity-effect note)
    1. A current in a wire makes a magnetic field around it and inside it. At DC the current spreads evenly over the wire: J = I / S everywhere.
    2. With the current alternating, the field inside the wire alternates too, and a changing field drives an electric field (Faraday: curl E = -dB/dt). Inside the wire that electric
       field pushes current. In the wire's own field the push is against the current at the centre and with it near the surface: the current is pushed OUT to the surface. That is the SKIN
       EFFECT, and its scale is the skin depth delta = sqrt(rho / (pi f mu0)): 5.4 mm at 200 Hz, 0.54 mm at 20 kHz in copper at 100 degC.
    3. A wire with NO net current in a field made by something else (a neighbour, the slot's walls) feels the same Faraday push: the field B (along y, say) changes in time, and an electric
       field E = w B x appears along the wire, proportional to the distance x from the line through the wire's axis that runs along the field. The wire answers with an eddy current: upward on
       one side of that line and downward on the other, so the net current is zero but the loss is not. That is the PROXIMITY EFFECT. At low frequency the current density is
       J = j w sigma B x (it grows linearly across the wire: the wire is 'rolled' by the field) and the loss per metre is the integral of sigma |E|^2 / 2 over the section, pi sigma w^2 B^2 d^4 / 128;
       at high frequency the eddy currents' own field shields the inside, and the current sits in a layer of thickness delta on the two sides facing across the field.
    4. Both are measured by one number, d / delta, the wire's diameter over the skin depth: nothing happens at d / delta = 0.1, the effects are felt at 1 and strong at 3 and above.

What the script does
    For a 0.63 mm copper wire at 100 degC it solves four frequencies (d / delta = 0.3, 1, 2, 4) twice in FEMM's time-harmonic solver (the project's solver slot, about 4 s in all):
      the skin effect:      the wire alone in a series circuit of 1 A (peak);
      the proximity effect: the wire in a series circuit of 0 A (no net current) in a uniform field of 1 T (peak) along y, brought in as A = -B x on the outer boundary (100 radii out).
    and draws, for each, the current density in the wire's cross-section (|J| over the DC density for the skin effect; |J| over its low-frequency scale w sigma B a for the proximity effect:
    1 at the edge at low frequency), and below them the profile across the wire: FEMM's triangles against the exact Bessel-function solution (optimization/ac_theory.py).
    It prints, for each frequency, R_ac / R_dc and the proximity loss over the low-frequency formula pi sigma w^2 B^2 d^4 / 128 from FEMM and from the exact solution.

Run (from the Python to FEMM folder)
    uv run python learn/proximity/01_one_wire.py            (about 5 s of FEMM, 8 solves)
Writes out/learn/proximity/one_wire.png and prints the table.

Numbers to expect (FEMM 5.0.20261101, mesh d / 14 in the wire and 2 degree arcs; the run prints them):
    d / delta      0.3         1          2          4
    R_ac / R_dc    1.000011    1.001295   1.020414   1.264971     exact 1.000011, 1.001301, 1.020492, 1.264643 (FEMM within 4e-4)
    proximity / low-frequency formula
                   0.99995     0.99320    0.90167    0.37945      exact 0.99994, 0.99289, 0.89764, 0.36931 (FEMM within 4e-4 up to d / delta = 1, 0.4 % above at 2, 2.7 % above at 4)
    The proximity error at d / delta = 4 is the mesh: the skin layer is d / 4 thick and the elements d / 14. optimization/ac_calibration.py (convergence stage) takes it apart; the demo keeps the
    cheap mesh because it is for pictures.

Self-check (the answers follow).
  1. At d / delta = 4 the picture shows the current density at the surface 1.5 times the DC value and 0.86 times at the centre, and the resistance has risen by only a quarter. How do the two fit?
  2. In the proximity picture the current density is zero along the vertical line through the wire's centre. Why?
Answers.
  1. The resistance is the area average of |J|^2 over the square of the DC density, for the same total current. Most of the section carries |J| near J_dc; the 1.5 J_dc is reached only in the thin
     ring at the surface (the profile at the bottom right of the figure), and averaging |J|^2 over the section gives 1.265. The picture shows J; the resistance averages J^2 over the area.
  2. The electric field E = w B x is proportional to the distance x from that line: by symmetry the push is up on one side of it and down on the other, so on the line itself there is no push
     (and no current), and the two sides' currents close through each other: no net current.
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

D_MM = 0.63
RATIOS = (0.3, 1.0, 2.0, 4.0)
B_PEAK = 1.0                      # T: the problem is linear, so the value is only a scale


def make_tasks() -> list[dict]:
    tasks = []
    for r in RATIOS:
        f = pc.frequency_for(r, D_MM)
        wire = {"x": 0.0, "y": 0.0, "d": D_MM, "circuit": "w"}
        tasks.append({"tag": f"one_skin_{r:g}", "f": f, "wires": [wire], "circuits": {"w": (1.0, "series")}, "B": None, "radius_mm": 20.0 * D_MM, "mesh_over_d": 1 / 14, "arc_deg": 2.0,
                      "fem_dir": str(pc.OUT / "fem")})
        tasks.append({"tag": f"one_prox_{r:g}", "f": f, "wires": [wire], "circuits": {"w": (0.0, "series")}, "B": B_PEAK, "radius_mm": 50.0 * D_MM, "mesh_over_d": 1 / 14, "arc_deg": 2.0,
                      "fem_dir": str(pc.OUT / "fem")})
    return tasks


def main() -> int:
    pc.OUT.mkdir(parents=True, exist_ok=True)
    tasks = make_tasks()
    results = pc.solve(tasks, "learn/proximity: one wire", workers=4)
    by = {t["tag"]: r for t, r in zip(tasks, results)}
    rho, d, L = pc.RHO, D_MM * 1e-3, 20e-3
    print("  d/delta   f [Hz]    R_ac/R_dc: FEMM   exact      proximity loss / low-frequency formula: FEMM   exact")
    for r in RATIOS:
        f = pc.frequency_for(r, D_MM)
        s, p = by[f"one_skin_{r:g}"], by[f"one_prox_{r:g}"]
        r_dc = rho * L / s["area_m2"][0]
        femm_skin = s["loss_W"][0] / (0.5 * r_dc)
        femm_prox = p["loss_W"][0] / (float(ac_theory.proximity_lowfreq(d, f, B_PEAK, rho)) * L)
        print(f"  {r:6.1f}  {f:9.0f}        {femm_skin:9.6f}      {float(ac_theory.skin_ratio(d, f, rho)):9.6f}               {femm_prox:9.6f}                     {float(ac_theory.proximity_ratio(d, f, rho)):9.6f}")
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    cmap = LinearSegmentedColormap.from_list("sb", [pc.SURFACE, "#9ec4ef", pc.BLUE, "#0b2a52"])
    fig = plt.figure(figsize=(12.0, 9.0), facecolor=pc.SURFACE)
    fig.text(0.012, 0.985, "One round wire as the frequency rises: where the current goes", ha="left", va="top", fontsize=11, color=pc.INK, fontweight="bold")
    fig.text(0.012, 0.012, f"A {D_MM} mm copper wire at 100 degC, FEMM 5.0 time-harmonic. Columns: d / delta = {', '.join(f'{r:g}' for r in RATIOS)}. Top: the wire alone with its own current (skin effect), |J| over the DC density. "
             "Middle: no net current in a uniform field (proximity effect), |J| over w sigma B a (1 at the edge at low frequency). Bottom: along the wire's x axis, FEMM's triangles (dots) and the exact solution (lines). "
             "Source: learn/proximity/01_one_wire.py.", ha="left", va="bottom", fontsize=7.6, color=pc.INK2, wrap=True)
    a = D_MM / 2.0                    # the wire's radius in mm (the pictures' axes are in mm)
    a_m = a * 1e-3                    # and in metres (the current densities are in A/m^2)
    for k, r in enumerate(RATIOS):
        f = pc.frequency_for(r, D_MM)
        s, p = by[f"one_skin_{r:g}"], by[f"one_prox_{r:g}"]
        j_dc = abs(s["current_A"][0]) / s["area_m2"][0]
        rel_s = np.abs(s["J"]).mean(axis=1) / j_dc
        scale_p = 2.0 * math.pi * f / rho * B_PEAK * a_m               # w sigma B a [A/m^2]: the low-frequency current density at the wire's edge
        rel_p = np.abs(p["J"]).mean(axis=1) / scale_p
        ax = fig.add_subplot(3, 4, k + 1)
        coll = pc.plot_wires(ax, s, rel_s, cmap, Normalize(0, max(1.0, float(rel_s.max()))), lim=a * 1.05)
        ax.set_title(f"d / delta = {r:g}: skin, |J| / J_dc  (max {rel_s.max():.2f})", fontsize=9, color=pc.INK, loc="left", fontweight="bold")
        cb = fig.colorbar(coll, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=7.5, colors=pc.INK2)
        ax = fig.add_subplot(3, 4, 4 + k + 1)
        coll = pc.plot_wires(ax, p, rel_p, cmap, Normalize(0, max(1.0, float(rel_p.max()))), lim=a * 1.05)
        ax.set_title(f"proximity, |J| / (w sigma B a)  (max {rel_p.max():.2f})", fontsize=9, color=pc.INK, loc="left", fontweight="bold")
        cb = fig.colorbar(coll, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=7.5, colors=pc.INK2)
        ax = fig.add_subplot(3, 4, 8 + k + 1)
        pc.style_axes(ax, "x / a  (across the wire, along the field's normal)", "|J| relative" if k == 0 else "")
        xs = np.linspace(0, 1, 200)
        ax.plot(xs, np.abs(ac_theory.skin_profile(xs, d, f, rho)), color=pc.BLUE, lw=2.0, label="skin, exact")
        cen = s["tri"].mean(axis=1)
        sel = np.abs(cen[:, 1]) < a * 0.08
        ax.plot(np.abs(cen[sel, 0]) / a, np.abs(s["J"]).mean(axis=1)[sel] / j_dc, ".", color=pc.BLUE, ms=3, alpha=0.5)
        jp = np.abs(ac_theory.proximity_current_density(xs, np.zeros_like(xs), d, f, B_PEAK, rho)) / scale_p
        ax.plot(xs, jp, color=pc.ORANGE, lw=2.0, label="proximity, exact")
        cenp = p["tri"].mean(axis=1)
        selp = np.abs(cenp[:, 1]) < a * 0.08
        ax.plot(np.abs(cenp[selp, 0]) / a, np.abs(p["J"]).mean(axis=1)[selp] / scale_p, ".", color=pc.ORANGE, ms=3, alpha=0.5)
        if k == 0:
            leg = ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(0.0, 0.92))
            for t in leg.get_texts():
                t.set_color(pc.INK2)
    fig.subplots_adjust(left=0.05, right=0.96, top=0.93, bottom=0.10, wspace=0.28, hspace=0.32)
    path = pc.OUT / "one_wire.png"
    fig.savefig(path, dpi=130, facecolor=pc.SURFACE)
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
