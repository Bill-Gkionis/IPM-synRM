"""learn/proximity/04_rac_vs_frequency.py: R_ac / R_dc of a slot against frequency, and the formula that predicts it from ONE static solve.

The idea (the explanation for the proximity-effect note, from the formula to the winding)
    1. The loss of a strand of diameter d in a transverse field of peak B is pi sigma w^2 B^2 d^4 L / 128 (01_one_wire.py), while its DC loss at the current i is i^2 rho L / (2 S) with S = pi d^2 / 4.
    2. In a slot, B at a strand is proportional to the current: B = mu0 beta N n_s i, with N n_s i the current of the whole layer and beta [1/m] a number that depends only on WHERE the strand is in the slot
       (a strand-free solve of the slot's field gives it). Put B in the ratio of the two losses and the strand's d and the winding's turns and strands collapse into the copper fill k = N n_s S / A_layer and
       the skin depth delta (delta^2 = 2 / (w mu0 sigma)):
           (R_ac - R_dc) / R_dc  =  K (d k / delta^2)^2,       K = A_layer^2 <beta^2> / 4  a constant of the slot [m^2]
       so ONE magnetostatic solve of the slot gives K, and then the AC resistance of any winding at any frequency follows. (The derivation is in the docstring of optimization/ac_slot_study.py.)
    3. The law holds while the strands are thin against the skin depth (the lone wire's own shielding, optimization/ac_theory.proximity_ratio, takes over above d / delta of 1) and while the steel does
       not matter (its laminations' eddy currents reduce the field at high frequency).

What the script does
    For the winding study's example winding (26 turns x 8 strands of 0.315 mm, 4 paths, in the optimized stator's slot at the rated 118.6 ampere-turns a coil side; 416 strands) it solves
      the slot in FEMM's time-harmonic solver at twelve frequencies from 100 Hz to 200 kHz, twice: with the laminated stator steel (0.356 mm sheets) and with lossless steel (no eddy currents in it);
      the same slot once in the magnetostatic solver, from which K is made (R_ac / R_dc - 1 at 200 Hz over (d k / delta^2)^2: the law's constant, and the field of every strand);
    and draws R_ac / R_dc against frequency: FEMM's two curves, the law's curve (K from the static solve, times the lone wire's exact proximity factor, plus its skin effect), and the lone strand's skin
    effect alone. It prints the table. The project's solver slot is used (about 15 s of FEMM).

Run (from the Python to FEMM folder)
    uv run python learn/proximity/04_rac_vs_frequency.py
Writes out/learn/proximity/rac_vs_frequency.png and prints the table.

Numbers to expect (FEMM 5.0.20261101, 60,523 nodes, the random layout of seed 1; the run prints them): K = 15.889 mm^2 (fill 0.4476). With the steel's eddy currents switched off, the law's extra loss
(R_ac / R_dc - 1) is within 0.2 % of FEMM's up to 40 kHz (d / delta 0.83), 0.4 % at 60 kHz, 1.3 % at 100 kHz and 4.9 % at 200 kHz (d / delta 1.86; the law is above FEMM there). With the laminated
steel (0.356 mm sheets) FEMM's extra loss is 4.1 % below the law's at 5 kHz, 11 % at 20 kHz and 39 % at 200 kHz: the steel screens the field. R_ac / R_dc, laminated steel: 1.0004 at 200 Hz, 1.011 at
1,090 Hz, 1.230 at 5 kHz, 1.888 at 10 kHz, 4.397 at 20 kHz, 28.1 at 60 kHz. The lone strand's skin effect is 1.0002 at 20 kHz and 1.015 at 200 kHz. About 18 s with the solver slot free (25 solves).

Self-check (the answers follow).
  1. The law says (R_ac - R_dc) / R_dc goes as d^2 at a given fill. A designer asks: thinner strands, lower loss? What is the answer, and what is its cost?
  2. Why does the curve of the lone strand lie so far below the slot's?
Answers.
  1. Yes: the extra loss falls as the square of the strand's diameter, so 0.2 mm strands lose 2.5 times less extra than 0.315 mm ones, 6 times less than 0.5 mm, at the same copper fill. The cost is the fill: a
     thinner wire has relatively more enamel, so fewer strands fit in the slot (the winding study's wire catalogue and its 0.45 fill limit), and a winding of many strands is hard to wind.
  2. The lone strand has only its own field to fight (d / delta is 0.59 at 20 kHz: its skin effect is 1.0002). In the slot, each strand sits in the field of the whole slot's current, 400 strands' worth:
     62 mT near the opening against 0.72 mT at a strand's own surface (mu0 i / (2 pi a) at 0.570 A), about 90 times as much, and the loss goes as the square of the field.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))          # the code folder, wherever this is run from
import numpy as np

import femm_queue
import paths
import prox_common as pc
from optimization import ac_slot, ac_slot_study as study, ac_theory, design

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

N, N_S, D_MM, D_OVER_MM = 26, 8, 0.315, 0.367
AT_RATED = 118.6
FREQS = (100.0, 200.0, 500.0, 1090.0, 2000.0, 5000.0, 10000.0, 20000.0, 40000.0, 60000.0, 100000.0, 200000.0)
LOSSLESS = {"steel": {"mu_r": 4416.0, "sigma_MS": 0.0, "d_lam_mm": 0.0, "lam_fill": 1.0}}
A_LAYER_M2 = 36.218e-6


def main() -> int:
    pc.OUT.mkdir(parents=True, exist_ok=True)
    ov = json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"]
    i_s = AT_RATED / (N * N_S)
    fem = str(pc.OUT / "fem")
    tasks = [ac_slot.task_ms("rac_ms", fem, ov, n_turns=N, n_s=N_S, d_mm=D_MM, d_over_mm=D_OVER_MM, seed=1)]
    for tag, opt in (("lam", {}), ("lossless", LOSSLESS)):
        for f in FREQS:
            tasks.append(ac_slot.task_ac(f"rac_{tag}_{f:g}", fem, ov, n_turns=N, n_s=N_S, d_mm=D_MM, d_over_mm=D_OVER_MM, f=f, i_layers=(i_s, i_s), seed=1, options=opt))
    with femm_queue.femm_pool("learn/proximity: R_ac/R_dc against frequency", workers=8, task_timeout_s=300) as pool:
        res = pool.map(ac_slot.run_task, tasks)
    bad = [r for r in res if r.get("status") != "ok"]
    if bad:
        raise RuntimeError(str(bad[0].get("error")))
    by = {t["tag"]: r for t, r in zip(tasks, res)}
    ms = by["rac_ms"]
    rho = pc.RHO
    fill = N * N_S * math.pi * (D_MM * 1e-3) ** 2 / 4.0 / A_LAYER_M2
    # K from the static solve: the formula's total at one frequency over the DC loss, over (d k / delta^2)^2 (the d^4 f^2 scaling makes it frequency independent)
    p = design.params(ov)
    tab = study.strand_table(by["rac_lam_200"], ms, None, (i_s, i_s), 200.0, D_MM)
    delta200 = float(ac_theory.skin_depth(200.0, rho))
    K = float(tab["p_formula"].sum() / tab["p_dc"].sum()) / ((D_MM * 1e-3 * fill / delta200 ** 2) ** 2)
    print(f"the slot constant K from one static solve: {1e6 * K:.3f} mm^2  (fill {fill:.4f}, {len(ms['xy'])} strands, {ms['nodes']} nodes)")
    print("   f [Hz]   d/delta   FEMM laminated   FEMM lossless steel   the law (K, exact lone-wire factors)   lone strand's skin effect")
    rows = []
    for f in FREQS:
        delta = float(ac_theory.skin_depth(f, rho))
        d = D_MM * 1e-3
        law = 1.0 + (float(ac_theory.skin_ratio(d, f, rho)) - 1.0) + K * (d * fill / delta ** 2) ** 2 * float(ac_theory.proximity_ratio(d, f, rho))
        lam = by[f"rac_lam_{f:g}"]["loss_W"].sum() / ac_slot.dc_loss_W(by[f"rac_lam_{f:g}"], (i_s, i_s)).sum()
        ll = by[f"rac_lossless_{f:g}"]["loss_W"].sum() / ac_slot.dc_loss_W(by[f"rac_lossless_{f:g}"], (i_s, i_s)).sum()
        rows.append((f, d / delta, lam, ll, law, float(ac_theory.skin_ratio(d, f, rho))))
        print(f"  {f:8g}   {d / delta:6.3f}   {lam:13.5f}   {ll:17.5f}   {law:30.5f}   {rows[-1][5]:.6f}")
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    fig = plt.figure(figsize=(9.5, 5.6), facecolor=pc.SURFACE)
    fig.text(0.012, 0.985, "R_ac / R_dc of one slot against frequency: FEMM, and the law made from ONE static solve", ha="left", va="top", fontsize=11, color=pc.INK, fontweight="bold")
    fig.text(0.012, 0.012, "26 turns x 8 strands of 0.315 mm in the optimized stator's slot at 118.6 ampere-turns a coil side (416 strands), copper at 100 degC. The law: R_ac / R_dc = 1 + (lone strand's skin effect - 1) + "
             f"K (d k / delta^2)^2 x the lone wire's exact proximity factor, K = {1e6 * K:.2f} mm^2 from the magnetostatic solve. Source: learn/proximity/04_rac_vs_frequency.py.",
             ha="left", va="bottom", fontsize=7.6, color=pc.INK2, wrap=True)
    ax = fig.add_axes([0.08, 0.17, 0.89, 0.72])
    pc.style_axes(ax, "frequency [Hz]", "R_ac / R_dc  (the extra loss over the DC loss, plus 1)")
    f_ = np.array([r[0] for r in rows])
    ax.loglog(f_, [r[2] - 1 + 1e-9 for r in rows], "o", color=pc.ORANGE, ms=6, mec=pc.SURFACE, mew=1.2, label="FEMM, laminated steel (0.356 mm)")
    ax.loglog(f_, [r[3] - 1 + 1e-9 for r in rows], "s", color=pc.BLUE, ms=6, mec=pc.SURFACE, mew=1.2, label="FEMM, steel without eddy currents")
    ax.loglog(f_, [r[4] - 1 + 1e-9 for r in rows], "-", color=pc.INK2, lw=1.6, label="the law with K from the static solve")
    ax.loglog(f_, [r[5] - 1 + 1e-9 for r in rows], "--", color=pc.GRAY, lw=1.6, label="a lone strand's skin effect")
    ax.set_ylabel("R_ac / R_dc - 1", fontsize=9, color=pc.INK2)
    leg = ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    for t in leg.get_texts():
        t.set_color(pc.INK2)
    path = pc.OUT / "rac_vs_frequency.png"
    fig.savefig(path, dpi=130, facecolor=pc.SURFACE)
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
