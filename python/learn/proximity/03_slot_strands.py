"""learn/proximity/03_slot_strands.py: the strands of a real slot: those next to the opening lose most, the deep ones almost nothing, and the fundamental frequency hardly matters.

The idea (the explanation for the proximity-effect note, from the wire to the slot)
    A slot of the stator holds hundreds of strands. The current of all of them makes the slot's leakage field: it points across the slot (from one tooth wall to the other), and its strength at a given
    depth is Ampere's law applied to all the current BELOW that depth: B(y) = mu0 I(below y) / (slot width). So it is zero at the bottom of the slot, grows with the height above it, and is largest
    next to the opening, where it is also joined by the field of the rotor entering through the opening. Each strand is a wire in that field: its proximity loss is pi sigma w^2 d^4 |B|^2 / 128 per metre
    (01_one_wire.py), the strongest where the field is strongest. And it goes as the square of the frequency: at 200 Hz nothing, at the inverter's switching frequency something.

What the script does
    The winding study's example (optimization section 9.5): 26 turns a coil side of 8 strands of 0.315 mm each, 4 parallel paths, in the optimized stator's slot at the rated 118.6 ampere-turns a coil
    side (0.570 A peak in each of the 416 strands), strands packed at random in the two coil layers. FEMM 5.0 solves the slot (with the neighbouring teeth, the yoke, the opening and the air gap, 60 thousand
    nodes) at 200 Hz (the fundamental at 3,000 rpm), 5 kHz and 20 kHz (the switching frequency) and once in the static solver (the field of the slot at each strand), in the project's solver slot,
    about 7 s in all. It draws the slot with every strand coloured by its loss over
    its own DC loss, the current density inside the strands at the top and in the depth of the slot (the picture of crowding), and the loss of every strand against its depth. It prints the slot's R_ac / R_dc
    and where the extra loss is made. Everything about the model is in optimization/ac_slot.py; the studies of it (the frequency sweep, the formula strand by strand, circulating currents) are in
    optimization/ac_slot_study.py.

Run (from the Python to FEMM folder)
    uv run python learn/proximity/03_slot_strands.py
Writes out/learn/proximity/slot_strands.png and prints the table.

Numbers to expect (FEMM 5.0.20261101, 60,523 nodes, the random layout of seed 1; the run prints them): R_ac / R_dc 1.00038 at 200 Hz, 1.2298 at 5 kHz and 4.3967 at 20 kHz (the same-phase layers of
6 of the 9 slots of the sector; the three random layouts of optimization/ac_slot_study.py give 4.17 to 4.40 at 20 kHz, mean 4.27, and seed 1 is the highest); 63 % of the extra loss is made in the first 3 mm below the top of the coil (32 % in the first millimetre) and 9 % in layer 1, at every frequency; the slot's own field
at the strands, peak: 62 mT in the first millimetre, 47 mT at 1 to 3 mm, 33 mT at 3 to 6.5 mm and 12 mT on average in layer 1 (the largest 108 mT, the smallest 0.29 mT). About 7 s in all with the
solver slot free (four solves: three frequencies and the static field).

Self-check (the answers follow).
  1. At 200 Hz the slot's R_ac / R_dc is 1.00038. Is the proximity effect then negligible in this motor?
  2. Why does the picture of the current density at the top of the slot show the crowding inside each strand, but the strands deep in the slot show none?
Answers.
  1. At the machine's fundamental, in the slot's OWN field, yes. The proximity loss of the field from the rotor, with its harmonics (the slot-passing field is at 9 f_e) and the circulating currents between strands
     of one turn are separate matters, treated in ac_machine.py and the circulating-current stage; and the PWM ripple at 10 or 20 kHz meets an R_ac / R_dc of 1.9 and 4.4 (04_rac_vs_frequency.py).
  2. The crowding follows the field across the strand, and the field at the strand follows the slot's depth: 62 mT in the first millimetre below the top of the coil, 12 mT on average in layer 1, 0.3 mT
     for the deepest strand; the loss goes as the square of the field.
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
from optimization import ac_slot, ac_theory, design

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.patches import Circle

N, A_PATHS, N_S, D_MM, D_OVER_MM = 26, 4, 8, 0.315, 0.367
AT_RATED = 118.6
FREQS = (200.0, 5000.0, 20000.0)


def main() -> int:
    pc.OUT.mkdir(parents=True, exist_ok=True)
    ov = json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"]
    p = design.params(ov)
    i_s = AT_RATED / (N * N_S)                                         # A peak a strand
    tasks = [ac_slot.task_ac(f"slot_{f:g}", str(pc.OUT / "fem"), ov, n_turns=N, n_s=N_S, d_mm=D_MM, d_over_mm=D_OVER_MM, f=f, i_layers=(i_s, i_s), seed=1, want_field=True) for f in FREQS]
    tasks.append(ac_slot.task_ms("slot_static", str(pc.OUT / "fem"), ov, n_turns=N, n_s=N_S, d_mm=D_MM, d_over_mm=D_OVER_MM, seed=1))        # the same strands, magnetostatic: B at each
    with femm_queue.femm_pool("learn/proximity: slot strands", workers=4, task_timeout_s=300) as pool:
        res = pool.map(ac_slot.run_task, tasks)
    bad = [r for r in res if r.get("status") != "ok"]
    if bad:
        raise RuntimeError(str(bad[0].get("error")))
    ms = res.pop()                                                        # the static solve (one solve with layer 0 at 1 A a strand, one with layer 1)
    assert np.allclose(ms["xy"], res[0]["xy"]), "the static and the AC runs must have the same strand layout"
    b_peak = i_s * np.hypot(np.sum(ms["Bx_mean"], axis=0), np.sum(ms["By_mean"], axis=0))       # T, peak, at the centre of each strand, both layers at i_s in phase
    top = ac_slot.layer_polygons(p)[0].bounds[1]
    print(f"26 turns x 8 strands x 0.315 mm, 4 paths, {AT_RATED} ampere-turns a coil side, {i_s:.3f} A peak a strand, {res[0]['nodes']} nodes, copper at 100 degC")
    print("   f [Hz]   d/delta   R_ac/R_dc    share of the extra loss by depth below the top of the coil: 0-1 mm   1-3 mm   3-6.5 mm   layer 1")
    for f, r in zip(FREQS, res):
        p_dc = ac_slot.dc_loss_W(r, (i_s, i_s))
        extra = r["loss_W"] - p_dc
        depth = r["xy"][:, 1] - top
        shares = [float(extra[(depth >= lo) & (depth < hi)].sum() / extra.sum()) for lo, hi in ((0, 1), (1, 3), (3, 6.55), (6.55, 99))]
        print(f"  {f:8g}   {D_MM * 1e-3 / float(ac_theory.skin_depth(f, pc.RHO)):6.3f}   {r['loss_W'].sum() / p_dc.sum():9.5f}                                            "
              + "   ".join(f"{100 * s:6.1f} %" for s in shares))
    depth0 = res[0]["xy"][:, 1] - top
    bands = ((0, 1, "0-1 mm"), (1, 3, "1-3 mm"), (3, 6.55, "3-6.5 mm"), (6.55, 99, "layer 1"))
    print("the slot's own field at the strands (peak, at the strand centres, both layers in phase), mean |B| by depth below the top of the coil:")
    print("   " + "   ".join(f"{name}: {1e3 * b_peak[(depth0 >= lo) & (depth0 < hi)].mean():5.1f} mT" for lo, hi, name in bands))
    print(f"   largest |B| of any strand {1e3 * b_peak.max():.1f} mT, smallest {1e3 * b_peak.min():.2f} mT")
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    cmap = LinearSegmentedColormap.from_list("sb", [pc.SURFACE, "#9ec4ef", pc.BLUE, "#0b2a52"])
    fig = plt.figure(figsize=(12.0, 8.6), facecolor=pc.SURFACE)
    fig.text(0.012, 0.985, "The strands of a slot at 20 kHz: the loss is made next to the opening", ha="left", va="top", fontsize=11, color=pc.INK, fontweight="bold")
    fig.text(0.012, 0.012, "26 turns x 8 strands of 0.315 mm in the optimized stator's slot, 118.6 ampere-turns a coil side, copper at 100 degC, FEMM 5.0 time-harmonic. Left: every strand coloured by its loss over its "
             "own DC loss at 20 kHz. Middle: |J| over the DC density inside the strands, a window at the top of the coil and one 4 mm deeper. Right: the loss of each strand over its DC loss against its depth "
             "below the top of the coil, at three frequencies. Source: learn/proximity/03_slot_strands.py.", ha="left", va="bottom", fontsize=7.6, color=pc.INK2, wrap=True)
    r = res[2]
    p_dc = ac_slot.dc_loss_W(r, (i_s, i_s))
    ratio = r["loss_W"] / p_dc
    ax = fig.add_subplot(1, 3, 1)
    pcs = PatchCollection([Circle((x, y), D_MM / 2) for x, y in r["xy"]], cmap=cmap, norm=LogNorm(vmin=max(1.0, ratio.min()), vmax=ratio.max()), edgecolors="none")
    pcs.set_array(ratio)
    ax.add_collection(pcs)
    ax.set_xlim(-3.4, 3.4)
    ax.set_ylim(r["xy"][:, 1].min() - 0.3, r["xy"][:, 1].max() + 0.3)
    ax.set_aspect("equal")
    ax.set_title("loss of each strand / its DC loss", fontsize=9.5, color=pc.INK, loc="left", fontweight="bold")
    ax.set_xlabel("x across the slot [mm]", fontsize=9, color=pc.INK2)
    ax.set_ylabel("y, radial [mm]  (the opening is at the bottom of the picture)", fontsize=9, color=pc.INK2)
    ax.tick_params(colors=pc.INK2, labelsize=8)
    cb = fig.colorbar(pcs, ax=ax, fraction=0.046, pad=0.02)
    cb.set_ticks([2, 5, 10, 20, 30])
    cb.ax.set_yticklabels(["2", "5", "10", "20", "30"])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=7.5, colors=pc.INK2)
    field = r["field"]
    j_dc = np.abs(r["current_A"][field["strand"]]) / r["area_m2"][field["strand"]]
    rel = np.abs(field["J_nodes"]).mean(axis=1) / j_dc
    verts = field["tri_xy"]
    for k, (y0, label) in enumerate(((top, "a window at the top of the coil"), (top + 4.0, "a window 4 mm deeper"))):
        ax2 = fig.add_axes([0.38 + 0.0, 0.54 - 0.45 * k, 0.26, 0.40])
        sel = (np.abs(verts[:, :, 0]).max(axis=1) < 1.3) & (verts[:, :, 1].max(axis=1) < y0 + 2.4) & (verts[:, :, 1].min(axis=1) > y0 - 0.1)
        pcl = PolyCollection(verts[sel], array=rel[sel], cmap=cmap, norm=Normalize(float(np.percentile(rel, 1)), float(np.percentile(rel, 99.5))), edgecolors="none")
        ax2.add_collection(pcl)
        ax2.set_xlim(-1.3, 1.3)
        ax2.set_ylim(y0 - 0.1, y0 + 2.4)
        ax2.set_aspect("equal")
        ax2.set_title(f"|J| / J_dc, {label}", fontsize=9, color=pc.INK, loc="left", fontweight="bold")
        ax2.tick_params(colors=pc.INK2, labelsize=7.5)
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax2.spines[side].set_color(pc.AXIS)
        cb2 = fig.colorbar(pcl, ax=ax2, fraction=0.046, pad=0.02)
        cb2.ax.tick_params(labelsize=7.5, colors=pc.INK2)
    ax3 = fig.add_axes([0.71, 0.12, 0.27, 0.78])
    pc.style_axes(ax3, "depth below the top of the coil [mm]", "strand loss / its DC loss - 1")
    for f, rr, col in zip(FREQS, res, (pc.GRAY, pc.AQUA, pc.ORANGE)):
        pdc = ac_slot.dc_loss_W(rr, (i_s, i_s))
        ax3.semilogy(rr["xy"][:, 1] - top, np.maximum(rr["loss_W"] / pdc - 1.0, 1e-7), ".", ms=4, color=col, alpha=0.8, label=f"{f / 1e3:g} kHz")
    leg = ax3.legend(frameon=False, fontsize=8.5, markerscale=2.5)
    for t in leg.get_texts():
        t.set_color(pc.INK2)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.93, bottom=0.12, wspace=0.30)
    path = pc.OUT / "slot_strands.png"
    fig.savefig(path, dpi=130, facecolor=pc.SURFACE)
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
