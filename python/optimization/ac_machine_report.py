"""optimization.ac_machine_report: the analysis of the machine's field at the copper (no FEMM): the check against the first study, the table of the operating points, the fits the fast
model uses, the figure.

What it is, and why
    optimization/ac_machine.py solves the sector at the operating points and keeps the raw fields. This module turns them into the numbers of the report, and into the tables that
    optimization/copper_loss.py packages: for every operating point the harmonic spectrum G_k of the mean square dB/d(theta_e) over the coil area (and per coil layer and per depth in the slot),
    the spectrum H_k of the vector potential's spatial variance inside a layer (the circulating currents' source), the flux linkages and the torque. Then:
    1. THE CHECK against the first study. Its femm.json (out/optimization/winding/optimized/) holds G for the same geometry, the same fine mesh and the same 48 positions at 100 ampere-turns and
       its own MTPA angle: the 'first study' point of ac_machine must give the same G to the last digits (FEMM 5.0 is deterministic here), or something in this pipeline is wrong.
    2. THE POINTS: ampere-turns, angle, i_d, i_q, torque, flux linkages, G and its harmonic split, per layer.
    3. THE FITS. G_k, H_k, psi_d, psi_q and the torque are fitted by a quadratic polynomial in (i_d, i_q) (6 coefficients each: 1, i_d, i_q, i_d^2, i_d i_q, i_q^2) by least squares over the
       operating points. A leave-one-out error says how well a point is predicted by the fit of the others: the fit is an interpolation over the region of the points, not an extrapolation.
    4. THE PICTURE: where in the slot the field changes fastest (G by depth) at no load, at the rated MTPA point, at 3 x rated and in field weakening, and the harmonic spectrum.

Run
    Not on its own: ac_machine.py --stage analyse calls report(points). Writes out/optimization/copper_ac/machine/summary.json, tables.json (what copper_loss.py packages) and
    fig_machine_field.png, and prints the tables.

Self-check (the answers follow).
  1. The 'first study' point must reproduce its G to the last digits. Why is that a check and not a coincidence of two codes that share their input?
  2. Why fit G_k, a quadratic quantity, with a quadratic polynomial in the currents?
Answers.
  1. The two share the geometry, the mesh and the currents, but not the code that reads the field: ac_machine.field_task is its own function, the Fourier analysis here is its own, and the
     comparison is made on G, a number made from every triangle and every position. A mistake in the position set, the triangle selection, the area weighting or the harmonic weights would show.
  2. If the field in the copper were linear in the currents, B_k = B_k0 + i_d B_kd + i_q B_kq at every point and harmonic, then |B_k|^2 summed over the area is exactly a quadratic form in
     (1, i_d, i_q). The iron is not linear (saturation), so the fit is checked point by point; where the field is mostly the rotor's or mostly the armature's, the quadratic form is nearly exact.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import paths
from optimization import ac_machine as am

FIRST_STUDY = paths.path("out") / "optimization" / "winding" / "optimized" / "femm.json"


def _first_study() -> dict | None:
    if not FIRST_STUDY.exists():
        return None
    return json.loads(FIRST_STUDY.read_text(encoding="utf-8"))


def report(points: list[dict], log=print) -> dict:
    stats = []
    for k in range(len(points)):
        raw = am.load_point(k)
        st = am.harmonic_stats(raw)
        st["name"], st["amps"], st["gamma"], st["i_d"], st["i_q"], st["at"] = points[k]["name"], points[k]["amps"], points[k]["gamma"], points[k]["i_d"], points[k]["i_q"], points[k]["at"]
        st["nodes"], st["coil_triangles"] = raw["nodes"], raw["coil_triangles"]
        stats.append(st)
    out: dict = {"points": points, "stats": stats}
    # ---- the check against the first study
    first = _first_study()
    idx = next(i for i, q in enumerate(points) if q["group"] == "reference")
    if first is not None:
        g_ref = first["statistics"]["48"]["G"]
        gk_ref = np.array(first["statistics"]["48"]["G_k"])
        g_new = stats[idx]["G"]
        log(f"\nTHE CHECK against the first study (winding/optimized/femm.json): G at 100 ampere-turns and {points[idx]['gamma']:.4f} deg")
        log(f"  G = {g_new:.10e} T^2/rad^2 here, {g_ref:.10e} there: relative difference {g_new / g_ref - 1:+.2e}")
        gk_new = np.array(stats[idx]["G_k"])
        log(f"  the harmonics: largest relative difference {np.max(np.abs(gk_new[1:] - gk_ref[1:]) / np.maximum(gk_ref[1:], 1e-9 * g_ref)):.2e} over k = 1 ... {len(gk_ref) - 1}")
        psi_d_ref = float(np.mean(first["psi_d_A"]))
        log(f"  mean psi_d of the 24 positions of the first study's set A: {psi_d_ref * 1e3:.6f} mWb there; here (the same 24 positions): "
            f"{np.mean(am.load_point(idx)['psi_d'][::2]) * 1e3:.6f} mWb")
        out["check_first_study"] = {"G_here": g_new, "G_there": g_ref, "relative": g_new / g_ref - 1}
    # ---- the points
    log("\nTHE OPERATING POINTS (currents at the 10-turn model; torque and flux linkages of the whole machine; G over the coil area)")
    log("  name                          A pk   AT    gamma    i_d     i_q   T_mean  ripple%   psi_d    psi_q    G [1e-4]   top layer  bottom   k=1 %  5 %  7 %  9 %  rest %")
    for st in stats:
        gk = np.array(st["G_k"])
        share = lambda k: 100 * gk[k] / st["G"] if st["G"] > 0 else 0.0
        rest = 100 - sum(share(k) for k in (1, 5, 7, 9)) if st["G"] > 0 else 0.0
        top, bot = float(np.sum(st["G_k_layer"][0])), float(np.sum(st["G_k_layer"][1]))
        ripple = f"{100 * st['T_ripple_pp'] / st['T_mean']:7.2f}" if abs(st["T_mean"]) > 1e-3 else "      -"           # a ratio to a torque of zero (no load) means nothing
        log(f"  {st['name']:28s} {st['amps']:6.2f} {st['at']:6.1f} {st['gamma']:6.2f} {st['i_d']:7.3f} {st['i_q']:7.3f} {st['T_mean']:7.4f} {ripple} "
            f"{st['psi_d'] * 1e3:8.3f} {st['psi_q'] * 1e3:8.3f}  {1e4 * st['G']:8.3f}  {1e4 * top:8.3f} {1e4 * bot:8.3f}   {share(1):5.1f} {share(5):4.1f} {share(7):4.1f} {share(9):4.1f} {rest:5.1f}")
    # ---- the fits
    gk = np.array([st["G_k"] for st in stats])                 # (points, 25)
    gk_l0 = np.array([st["G_k_layer"][0] for st in stats])
    gk_l1 = np.array([st["G_k_layer"][1] for st in stats])
    hk = np.array([st["H_k"] for st in stats])
    psi = np.array([[st["psi_d"], st["psi_q"], st["T_mean"]] for st in stats])
    fit_g = am.fit_split(points, gk)
    fit_h = am.fit_split(points, hk)
    fit_psi = am.fit_quadratic(points, psi)
    names = [st["name"] for st in stats]
    log("\nTHE FITS  V_g = V_g(no load) + AT^2 g_tot(AT, gamma) a_g(AT, gamma) per harmonic group, for G (the field) and H (the vector potential's spread); quadratic in (i_d, i_q) for psi and T")
    log(f"  groups of harmonics k: {[g if len(g) == 1 else f'{g[0]}-{g[-1]}' for g in am.GROUPS]}; representative k of each group: {np.round(fit_g['k_eff'], 2).tolist()}")
    log(f"  total G: rms residual {100 * fit_g['rms_total_rel']:.2f} % of the value; leave-one-out (a point predicted by the fit of the others):")
    for j, i in enumerate(fit_g["loaded"]):
        log(f"     leaving out {names[i]:28s}: predicted {1e4 * (stats[i]['G'] + fit_g['loo_total'][j]):8.3f} (true {1e4 * stats[i]['G']:8.3f}) e-4, error {100 * fit_g['loo_total_rel'][j]:+5.1f} %; "
            f"the harmonic groups' errors in per cent of the total: " + " ".join(f"{100 * x:+5.1f}" for x in fit_g["loo_group_rel"][j]))
    log(f"  H (A's spread): total rms residual {100 * fit_h['rms_total_rel']:.2f} %, leave-one-out error up to {100 * np.abs(fit_h['loo_total_rel']).max():.1f} %")
    log(f"  flux linkages and torque, quadratic fit: largest leave-one-out error psi_d {1e3 * np.abs(fit_psi['loo'][:, 0]).max():.2f} mWb (of {1e3 * np.ptp(psi[:, 0]):.1f}), "
        f"psi_q {1e3 * np.abs(fit_psi['loo'][:, 1]).max():.2f} mWb (of {1e3 * np.ptp(psi[:, 1]):.1f}), T {np.abs(fit_psi['loo'][:, 2]).max():.2f} N.m (of {np.ptp(psi[:, 2]):.2f})")
    ser = lambda d: {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in d.items()}
    out["fits"] = {"G": ser(fit_g), "H": ser(fit_h), "psi_T": {"coef": fit_psi["coef"].tolist(), "rms": fit_psi["rms"].tolist(), "loo": fit_psi["loo"].tolist()}}
    out["ranges"] = {"i_d": [min(q["i_d"] for q in points), max(q["i_d"] for q in points)], "i_q": [min(q["i_q"] for q in points), max(q["i_q"] for q in points)],
                     "amps_max": max(q["amps"] for q in points), "gamma_max": max(q["gamma"] for q in points), "gamma_min_loaded": min(q["gamma"] for q in points if q["at"] > 0),
                     "at_max": max(q["at"] for q in points)}
    depth = np.array(stats[0]["depth_mm"])
    out["depth_mm"] = depth.tolist()
    (am.OUT / "summary.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    pack = lambda f: {"v0": f["v0"].tolist(), "g_coef": f["g_coef"].tolist(), "a_coef": f["a_coef"].tolist(), "k_eff": f["k_eff"].tolist()}
    tables = {"i_ref_turns": am.I_REF_TURNS, "at_rated": am.AT_RATED, "groups": [g if len(g) == 1 else [g[0], g[-1]] for g in am.GROUPS], "ranges": out["ranges"],
              "G": pack(fit_g), "H": pack(fit_h), "psi_T_coef": fit_psi["coef"].tolist(),
              "points": [{"name": q["name"], "i_d": q["i_d"], "i_q": q["i_q"], "at": q["at"], "gamma": q["gamma"]} for q in points],
              "G_k_points": gk.tolist(), "H_k_points": hk.tolist(), "psi_T_points": psi.tolist()}
    (am.OUT / "tables.json").write_text(json.dumps(tables, indent=1), encoding="utf-8")
    figure(stats, am.OUT / "fig_machine_field.png")
    return out


def figure(stats: list[dict], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    surface, ink, ink2, grid, axis = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.0), facecolor=surface)
    pick = {"no load": (gray, "no load (the rotor's field alone)"), "MTPA 1 x rated": (blue, "rated, MTPA (118.6 ampere-turns)"), "MTPA 3 x rated": (orange, "3 x rated, MTPA"),
            "FW 1 x rated, 70 deg": (aqua, "rated current, field weakening 70 deg")}
    ax = axes[0]
    for st in stats:
        if st["name"] in pick:
            col, lab = pick[st["name"]]
            ax.plot(np.array(st["depth_mm"]) - 45.8, np.array(st["G_depth"]) * 1e4, color=col, lw=2.0, marker="o", ms=3.5, mec=surface, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("depth below the bore [mm]  (the opening is at 0, the coil starts at 0.9 mm, the layers meet at 7.5 mm)", fontsize=8.5, color=ink2)
    ax.set_ylabel("G = mean square of dB/d(theta_e) [1e-4 T^2/rad^2]", fontsize=8.5, color=ink2)
    ax.set_title("Where in the slot the field changes fastest", fontsize=10, color=ink, loc="left", fontweight="bold")
    ax = axes[1]
    for st in stats:
        if st["name"] in pick:
            col, lab = pick[st["name"]]
            k = np.arange(1, len(st["G_k"]))
            ax.plot(k, np.array(st["G_k"][1:]) * 1e4, color=col, lw=1.8, marker="o", ms=3.5, mec=surface, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("harmonic k of the electrical angle  (k = 9: the slot-passing harmonic, 36 slots / 4 pole pairs)", fontsize=8.5, color=ink2)
    ax.set_ylabel("G_k [1e-4 T^2/rad^2]", fontsize=8.5, color=ink2)
    ax.set_ylim(1e-8, 3e3)                                         # headroom for the legend above the data
    ax.set_title("The harmonics the strands see (G_k weighs k^2)", fontsize=10, color=ink, loc="left", fontweight="bold")
    for ax in axes:
        ax.set_facecolor(surface)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(axis)
        ax.tick_params(colors=ink2, labelsize=8.5)
        ax.grid(True, color=grid, lw=0.75)
        ax.set_axisbelow(True)
        leg = ax.legend(frameon=False, fontsize=8.2, loc="best" if ax is axes[0] else "upper right")
        for t in leg.get_texts():
            t.set_color(ink2)
    fig.text(0.012, 0.01, "FEMM 5.0 magnetostatic, fine mesh, the sector's 6,176 coil triangles over 48 rotor positions of one electrical period; area-weighted means; G is the mean square of dB/d(theta_e), "
             "G_k its share at harmonic k of the electrical angle (the sum over k is G). Source: optimization/ac_machine.py.", fontsize=7.5, color=ink2, wrap=True)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, dpi=130, facecolor=surface)
    plt.close(fig)
