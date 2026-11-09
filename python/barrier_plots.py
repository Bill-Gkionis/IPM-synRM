"""Figures and numbers for the air-band study, from barrier_study.py's CSVs and .ans files (no FEMM needed).

Run: uv run python barrier_plots.py <output folder>   (from this folder; writes PNGs + summary.csv there)
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.tri import Triangulation

import femfile
import motorgeo as mg

SRC = os.path.abspath("out/barrier")
DST = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "out/barrier/figs")
PREFIX = "Air band study - "
VARIANTS = ["with band", "no band"]
COL = {"with band": "#2a78d6", "no band": "#eb6834"}          # categorical slots 1 and 2
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUES = LinearSegmentedColormap.from_list("blues", [SURF, "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
                                                    "#184f95", "#0d366b"])
P = mg.Params()
PP = P.pole_pairs
RPM = 1000.0
W_E = PP * RPM * 2 * math.pi / 60                               # electrical rad/s at 1000 rpm

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlesize": 10, "axes.titlecolor": INK,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "legend.frameon": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF, "lines.linewidth": 1.8})


def save(fig, name):
    path = os.path.join(DST, PREFIX + name + ".png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def harmonics(y):
    """Amplitudes of harmonics 0..N/2 of one period sampled evenly (last point not repeated)."""
    c = np.fft.rfft(y) / len(y)
    c[1:] *= 2
    return np.abs(c), c


# ----------------------------------------------------------------------------- the numbers
def load():
    """Read the study's CSVs and boil them down to one dict of numbers (and curves) per variant."""
    ph1 = pd.read_csv(os.path.join(SRC, "phase1.csv"))
    ph2 = pd.read_csv(os.path.join(SRC, "phase2.csv"))
    S = {v: summarise(ph1, ph2, v) for v in VARIANTS}
    summary = pd.DataFrame({v: {k: x for k, x in S[v].items() if not k.startswith("_")} for v in VARIANTS})
    summary["change_pct"] = 100 * (summary["no band"] / summary["with band"] - 1)
    summary.to_csv(os.path.join(SRC, "summary.csv"))
    pd.set_option("display.width", 160)
    print(summary.round(4).to_string())
    return ph1, ph2, S


def summarise(ph1, ph2, v):
    s = {}
    nl = ph1[(ph1.study == "noload") & (ph1.variant == v)].sort_values("theta_e")
    amp, _ = harmonics(nl.psi_a.values)
    s["psi_pm_mWb"] = amp[1] * 1e3                               # fundamental of phase-A flux linkage
    s["psi_d0_mWb"] = nl.psi_d.mean() * 1e3
    X = np.fft.rfft(nl.psi_a.values)
    e = np.fft.irfft(1j * np.arange(len(X)) * X, len(nl)) * W_E  # e = d(psi)/dt = w_e d(psi)/d(theta_e)
    s["_e"], s["_theta"] = e, nl.theta_e.values
    e_amp, _ = harmonics(e)
    s["_e_amp"] = e_amp
    s["emf_peak_V_1000rpm"] = e_amp[1]
    s["emf_thd_pct"] = 100 * math.sqrt((e_amp[2:] ** 2).sum()) / e_amp[1]
    s["ke_Vs_per_rad"] = e_amp[1] / (RPM * 2 * math.pi / 60)
    cg = ph1[(ph1.study == "cogging") & (ph1.variant == v)].sort_values("theta_e")
    s["cogging_pp_mNm"] = (cg.T_age.max() - cg.T_age.min()) * 1e3
    g = ph1[(ph1.study == "gamma") & (ph1.variant == v)]
    gm = g.groupby(["I", "gamma"])[["T_age", "T_dq", "psi_d", "psi_q", "id", "iq"]].mean().reset_index()
    s["_gamma"] = gm
    for I in sorted(gm.I.unique()):
        row = gm[gm.I == I]
        k = int(np.argmax(row.T_age.values)); k = min(max(k, 1), len(row) - 2)
        a2, a1, a0 = np.polyfit(row.gamma.values[k - 1:k + 2], row.T_age.values[k - 1:k + 2], 2)
        g_opt = -a1 / (2 * a2)
        s[f"mtpa_gamma_{I:g}A"] = g_opt
        s[f"T_mtpa_{I:g}A"] = a0 + a1 * g_opt + a2 * g_opt ** 2
        s[f"T_q_{I:g}A"] = float(row[row.gamma == 0].T_age.iloc[0])
        d = row[np.isclose(row.gamma, 90)].iloc[0]
        q = row[np.isclose(row.gamma, 0)].iloc[0]
        s[f"Ld_{I:g}A_mH"] = (d.psi_d - s["psi_d0_mWb"] / 1e3) / d.id * 1e3
        s[f"Lq_{I:g}A_mH"] = q.psi_q / q.iq * 1e3
        s[f"saliency_{I:g}A"] = s[f"Lq_{I:g}A_mH"] / s[f"Ld_{I:g}A_mH"]
    rp = ph2[(ph2.study == "ripple") & (ph2.variant == v)].sort_values("theta_e")
    T = rp.T_age.values[:-1]                                     # 0..59: one 60-deg period, 60 is the repeat
    s["_ripple"] = rp
    s["T_mean_mtpa_10A"] = T.mean()
    s["ripple_pp_Nm"] = T.max() - T.min()
    s["ripple_pct"] = 100 * (T.max() - T.min()) / T.mean()
    r_amp, _ = harmonics(T)
    s["_r_amp"] = r_amp
    gap = pd.read_csv(os.path.join(SRC, f"gap_noload_{v.replace(' ', '_')}.csv"))
    b_amp, _ = harmonics(gap.br.values)
    s["_gap"], s["_b_amp"] = gap, b_amp
    s["B1_gap_T"] = b_amp[1]
    s["B_gap_thd_pct"] = 100 * math.sqrt((b_amp[2:60] ** 2).sum()) / b_amp[1]
    s["T_pm_10A_Nm"] = 1.5 * PP * s["psi_pm_mWb"] / 1e3 * 10 * math.cos(math.radians(s["mtpa_gamma_10A"]))
    return s


# ----------------------------------------------------------------------------- field pictures from the .ans
def field_mesh(path, rotor_deg):
    """Triangulation of the solution, rotor turned to its real angle, tiled over 270 deg (periodic), pole 0 centred."""
    sol = femfile.read_ans(path)
    xy, A, tri = sol["xy"].copy(), sol["A"], sol["tri"]
    r_split = P.r_ro + P.airgap / 2                              # middle of the unmeshed AGE band
    rotor = np.zeros(len(xy), bool)
    rotor[np.unique(tri[np.hypot(*xy[tri].mean(axis=1).T) < r_split])] = True
    rot = lambda pts, deg: pts @ np.array([[math.cos(math.radians(deg)), math.sin(math.radians(deg))],
                                           [-math.sin(math.radians(deg)), math.cos(math.radians(deg))]])
    xy[rotor] = rot(xy[rotor], rotor_deg)
    centre = 22.5 + rotor_deg                                    # pole 0 (north) d-axis, stator frame
    tiles = [rot(xy, 90 * k) for k in (-1, 0, 1)]                # periodic: the field repeats every 90 deg
    allxy = rot(np.vstack(tiles), 90 - centre)                   # turn the picture: pole 0 points up
    alltri = np.vstack([tri + k * len(xy) for k in range(3)])
    allA = np.concatenate([A] * 3)
    bx, by = femfile.element_b(xy, A, tri)
    allB = np.concatenate([np.hypot(bx, by)] * 3)
    ang = np.degrees(np.arctan2(*allxy[alltri].mean(axis=1).T[::-1])) - 90
    return Triangulation(allxy[:, 0], allxy[:, 1], alltri, mask=np.abs(ang) > 45), allA, allB, centre


def outlines(ax, v, rotor_deg, centre, zoom=False):
    """Region outlines; a steel-filled band is drawn dotted, to show where the band would be."""
    regions = mg.build(mg.Params(barrier="air"))
    for e, r in mg.unique_edges(regions):
        if r.name.startswith("gap_") or r.material == "Coil" and zoom:
            continue
        pts = np.array(e.points() + [e.b])
        if r.side == "rotor":
            pts = np.array([mg.rot(q, rotor_deg) for q in pts])
        style = (dict(color=INK2, lw=0.9, ls=(0, (1, 1.5))) if v == "no band" and r.name.startswith("barrier") else
                 dict(color=INK, lw=1.3, ls="-") if r.name.startswith("barrier") else
                 dict(color="#d1495b", lw=1.2, ls="-") if r.material == "N45" else
                 dict(color=INK2, lw=0.5, ls="-"))
        for k in (-1, 0, 1):
            q = np.array([mg.rot(xy, 90 * k + 90 - centre) for xy in pts])
            ang = np.degrees(np.arctan2(q[:, 1], q[:, 0])) - 90
            if np.all(np.abs(ang) > 46):
                continue
            ax.plot(q[:, 0], q[:, 1], **style)


def field_panel(ax, path, rotor_deg, v, levels, norm, title, zoom=None):
    t, A, B, centre = field_mesh(path, rotor_deg)
    pc = ax.tripcolor(t, facecolors=B, cmap=BLUES, norm=norm, shading="flat", rasterized=True)
    ax.tricontour(t, A, levels=levels, colors="#3a3a38", linewidths=0.45 if zoom is None else 0.55,
                  negative_linestyles="solid", alpha=0.8)
    outlines(ax, v, rotor_deg, centre, zoom=zoom is not None)
    ax.set_aspect("equal"); ax.grid(False)
    ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_visible(False)
    if zoom:
        ax.set_xlim(*zoom[:2]); ax.set_ylim(*zoom[2:])
    else:
        ax.set_xlim(-46, 46); ax.set_ylim(20, 65)
    ax.set_title(title, loc="left", fontsize=9.5)
    return pc


if __name__ == "__main__":
    os.makedirs(DST, exist_ok=True)
    ph1, ph2, S = load()
    snap = pd.concat([ph1, ph2]).dropna(subset=["ans"]).set_index("ans").rotor_deg      # saved .ans -> rotor angle
    fpath = lambda name: os.path.join(SRC, "fields", name + ".ans")

    # 1. field maps: no load and MTPA, both variants, same colour scale and same flux-line spacing
    norm = Normalize(0, 2.2)
    lv = np.linspace(-0.0072, 0.0072, 49)                        # equal steps of A: each gap = same flux
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 7.2), layout="constrained")
    for j, v in enumerate(VARIANTS):
        tag = v.replace(" ", "_")
        pc = field_panel(axs[0, j], fpath(f"noload_{tag}"), snap[f"noload_{tag}"], v, lv, norm,
                         f"{v}: no current (magnets only)")
        field_panel(axs[1, j], fpath(f"mtpa_{tag}"), snap[f"mtpa_{tag}"], v, lv, norm,
                    f"{v}: 10 A at its best current angle ({S[v]['mtpa_gamma_10A']:.0f}°)")
    cb = fig.colorbar(pc, ax=axs, shrink=0.6, pad=0.02)
    cb.set_label("|B| [T]", color=INK2); cb.outline.set_visible(False)
    fig.suptitle("Flux lines (grey) and flux density (blue) on a 90° window around one north pole. Every gap between two "
                 "neighbouring lines\ncarries the same flux, so more lines = more flux. Red: magnets. Black: the air band. "
                 "Dotted: where the band was.", x=0.01, ha="left", fontsize=10.5, color=INK)
    save(fig, "field maps")

    # 2. close-up of one magnet at no load: where its flux goes
    lv2 = np.linspace(-0.0072, 0.0072, 97)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.4), layout="constrained")
    for ax, v in zip(axs, VARIANTS):
        tag = v.replace(" ", "_")
        pc = field_panel(ax, fpath(f"noload_{tag}"), snap[f"noload_{tag}"], v, lv2, norm,
                         f"{v}: ψ(magnet → winding) = {S[v]['psi_pm_mWb']:.1f} mWb", zoom=(-16, 3, 30, 47))
    cb = fig.colorbar(pc, ax=axs, shrink=0.75, pad=0.02)
    cb.set_label("|B| [T]", color=INK2); cb.outline.set_visible(False)
    fig.suptitle("Close-up, no current: the left magnet of the pole. Lines that close around the magnet without "
                 "reaching the stator are leakage.", x=0.01, ha="left", fontsize=10.5, color=INK)
    save(fig, "magnet close-up")

    # 3. no-load flux linkage, back-EMF and its harmonics
    fig, axs = plt.subplots(3, 1, figsize=(8, 9.5), gridspec_kw={"height_ratios": [1, 1, 0.9]})
    for v in VARIANTS:
        nl = ph1[(ph1.study == "noload") & (ph1.variant == v)].sort_values("theta_e")
        axs[0].plot(nl.theta_e, nl.psi_a * 1e3, color=COL[v], label=v)
        axs[1].plot(S[v]["_theta"], S[v]["_e"], color=COL[v], label=v)
    for ax, yl, loc in zip(axs[:2], ("phase-A flux linkage [mWb]", f"phase-A back-EMF at {RPM:.0f} rpm [V]"),
                           ("lower left", "upper left")):
        ax.set_ylabel(yl); ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 60)); ax.legend(loc=loc)
    axs[1].set_xlabel("rotor position [electrical degrees]")
    axs[0].set_title("No current: the magnets alone, as the rotor turns", loc="left")
    orders = np.arange(1, 20, 2)
    w = 0.38
    for k, v in enumerate(VARIANTS):
        axs[2].bar(orders + (k - 0.5) * w, S[v]["_e_amp"][orders], width=w - 0.04, color=COL[v], label=v)
    axs[2].set_xticks(orders); axs[2].set_xlabel("harmonic (multiples of the electrical frequency)")
    axs[2].set_ylabel("back-EMF amplitude [V]"); axs[2].legend(loc="upper right")
    axs[2].set_title(f"Back-EMF harmonics. THD: with band {S['with band']['emf_thd_pct']:.1f} %, "
                     f"no band {S['no band']['emf_thd_pct']:.1f} %", loc="left")
    fig.tight_layout()
    save(fig, "no-load flux linkage and back-EMF")

    # 4. air-gap flux density at no load + spectrum
    fig, axs = plt.subplots(2, 1, figsize=(8, 6.6), gridspec_kw={"height_ratios": [1.2, 1]})
    for v in VARIANTS:
        gp = S[v]["_gap"]
        axs[0].plot(gp.deg, gp.br, color=COL[v], lw=1.2, label=v)
    axs[0].axhline(0, color=AXIS, lw=0.8)
    axs[0].set_xlim(0, 90); axs[0].set_xlabel("position around the gap, stator frame [mechanical degrees]")
    axs[0].set_ylabel("radial B in the gap [T]"); axs[0].legend(loc="lower left")
    axs[0].set_title("No current: radial flux density in the middle of the air gap (one N and one S pole)", loc="left")
    orders = np.arange(1, 40, 2)
    for k, v in enumerate(VARIANTS):
        axs[1].bar(orders + (k - 0.5) * 0.8, S[v]["_b_amp"][orders], width=0.72, color=COL[v], label=v)
    axs[1].set_xticks(orders[::2]); axs[1].set_xlabel("harmonic (per electrical cycle)")
    axs[1].set_ylabel("amplitude [T]"); axs[1].legend(loc="upper right")
    axs[1].set_title(f"Fundamental: with band {S['with band']['B1_gap_T']:.3f} T, no band {S['no band']['B1_gap_T']:.3f} T. "
                     "The 17th/19th are the slot harmonics", loc="left")
    fig.tight_layout()
    save(fig, "air-gap flux density")

    # 5. cogging torque
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for v in VARIANTS:
        cg = ph1[(ph1.study == "cogging") & (ph1.variant == v)].sort_values("theta_e")
        ax.plot(cg.theta_e, cg.T_age * 1e3, "-o", ms=3.5, color=COL[v],
                label=f"{v}: {S[v]['cogging_pp_mNm']:.1f} mN·m peak-to-peak")
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.set_xlim(0, 20); ax.set_xlabel("rotor position [electrical degrees] (one cogging period = 20°)")
    ax.set_ylabel("cogging torque [mN·m]"); ax.legend(loc="upper right")
    ax.set_title("No current: cogging torque (the magnets pulling on the stator teeth)", loc="left")
    fig.tight_layout()
    save(fig, "cogging")

    # 6. torque vs current angle, one panel per current level
    levels_I = sorted(S["with band"]["_gamma"].I.unique())
    fig, axs = plt.subplots(1, len(levels_I), figsize=(14, 3.9), sharex=True)
    for ax, I in zip(axs, levels_I):
        for v in VARIANTS:
            gm = S[v]["_gamma"]; row = gm[gm.I == I]
            ax.plot(row.gamma, row.T_age, color=COL[v], label=v)
            ax.plot(S[v][f"mtpa_gamma_{I:g}A"], S[v][f"T_mtpa_{I:g}A"], "o", ms=6, color=COL[v], mec=SURF, mew=1.5)
        ax.set_title(f"{I:g} A peak = {I * P.turns:g} A-turns", loc="left")
        ax.set_xlabel("current angle γ [°] (0 = all q)")
        ax.set_xticks(range(0, 91, 30)); ax.set_ylim(bottom=0)
    axs[0].set_ylabel("mean torque [N·m]"); axs[0].legend(loc="lower center")
    fig.suptitle("Torque vs current angle (mean over 6 rotor positions). Dots: the best angle (MTPA) for each design",
                 x=0.01, ha="left", fontsize=10.5, color=INK)
    fig.tight_layout()
    save(fig, "torque vs current angle")

    # 7. where the torque comes from at 10 A
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
    for ax, v in zip(axs, VARIANTS):
        gm = S[v]["_gamma"]; row = gm[gm.I == 10].sort_values("gamma")
        t_pm = 1.5 * PP * S[v]["psi_pm_mWb"] / 1e3 * row.iq
        ax.plot(row.gamma, row.T_age, color=COL[v], label="total (FEMM)")
        ax.plot(row.gamma, t_pm, color=COL[v], ls="--", lw=1.4, label="magnet: 1.5·p·ψ_PM·i_q")
        ax.plot(row.gamma, row.T_age - t_pm, color=COL[v], ls=":", lw=1.6, label="the rest: reluctance")
        ax.set_title(f"{v}, 10 A", loc="left"); ax.set_xlabel("current angle γ [°]"); ax.set_xticks(range(0, 91, 15))
        ax.axhline(0, color=AXIS, lw=0.8); ax.legend(loc="upper right")
    axs[0].set_ylabel("torque [N·m]")
    fig.suptitle("Splitting the torque: magnet part (from the no-load flux linkage) and the reluctance remainder",
                 x=0.01, ha="left", fontsize=10.5, color=INK)
    fig.tight_layout()
    save(fig, "torque split")

    # 8. torque ripple at MTPA
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.8), gridspec_kw={"width_ratios": [1.6, 1]})
    for v in VARIANTS:
        rp = S[v]["_ripple"]
        axs[0].plot(rp.theta_e, rp.T_age, color=COL[v],
                    label=f"{v}: mean {S[v]['T_mean_mtpa_10A']:.3f} N·m, ripple {S[v]['ripple_pct']:.1f} %")
    axs[0].set_xlim(0, 60); axs[0].set_xlabel("rotor position [electrical degrees] (one full ripple period = 60°)")
    axs[0].set_ylabel("torque [N·m]"); axs[0].set_ylim(0.6, 1.27); axs[0].legend(loc="lower left")
    axs[0].set_title("10 A at each design's MTPA angle", loc="left")
    orders = np.array([6, 12, 18, 24, 36])
    for k, v in enumerate(VARIANTS):
        axs[1].bar(np.arange(len(orders)) + (k - 0.5) * 0.38, S[v]["_r_amp"][orders // 6], width=0.34,
                   color=COL[v], label=v)
    axs[1].set_xticks(range(len(orders)), [f"{o}th" for o in orders])
    axs[1].set_ylabel("amplitude [N·m]"); axs[1].legend(loc="upper right")
    axs[1].set_title("ripple harmonics (per electrical cycle)", loc="left")
    fig.tight_layout()
    save(fig, "torque ripple")

    # 9. inductances and saliency vs current
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for v in VARIANTS:
        Ld = [S[v][f"Ld_{I:g}A_mH"] for I in levels_I]
        Lq = [S[v][f"Lq_{I:g}A_mH"] for I in levels_I]
        axs[0].plot(levels_I, Lq, "-o", ms=4, color=COL[v], label=f"{v}: L_q")
        axs[0].plot(levels_I, Ld, "--s", ms=4, color=COL[v], label=f"{v}: L_d")
        axs[1].plot(levels_I, [q / d for q, d in zip(Lq, Ld)], "-o", ms=4, color=COL[v], label=v)
    axs[0].set_ylabel("apparent inductance [mH]"); axs[0].set_ylim(bottom=0); axs[0].legend(loc="lower left", ncol=2)
    axs[1].set_ylabel("saliency L_q / L_d"); axs[1].legend(loc="lower left"); axs[1].set_ylim(bottom=1)
    for ax in axs:
        ax.set_xlabel("current [A peak] (x 10 turns = A-turns)"); ax.set_xticks(levels_I)
    axs[0].set_title("Inductances fall as the steel saturates", loc="left")
    axs[1].set_title("Saliency: what makes reluctance torque", loc="left")
    fig.tight_layout()
    save(fig, "inductances")

    # 10. the 90 deg slice against the hand-drawn full 360 deg model (validate_full.py)
    vpath = os.path.abspath("out/validate/validate.csv")
    if os.path.exists(vpath):
        val = pd.read_csv(vpath)
        style = {"motorgeo 90": dict(color=COL["with band"], lw=1.8), "hand-drawn 360": dict(color="#1baf7a", lw=1.8)}
        name = {"motorgeo 90": "motorgeo, 90° slice (fine mesh, 27.9k nodes)", "hand-drawn 360": "the hand-drawn .FEM, full 360°"}
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))
        for model, g in val.groupby("model"):
            for I, gi in g.groupby("I"):
                gi = gi.sort_values("theta_e")
                axs[0].plot(gi.theta_e, gi.psi_a * 1e3, ls="-" if I == 0 else "--", marker="o", ms=3.5,
                            label=f"{name[model]}, {'no current' if I == 0 else f'{I:g} A'}", **style[model])
                if I > 0:
                    axs[1].plot(gi.theta_e, gi["T"], marker="o", ms=3.5, label=name[model], **style[model])
        axs[0].set_ylabel("phase-A flux linkage [mWb]"); axs[0].legend(loc="upper center", fontsize=7.5)
        axs[1].set_ylabel("torque [N·m]"); axs[1].legend(loc="upper left", fontsize=7.5)
        axs[1].axhline(0, color=AXIS, lw=0.8)
        for ax in axs:
            ax.set_xlabel("rotor north-pole d-axis minus phase-A axis [electrical degrees]")
            ax.set_xlim(0, 345); ax.set_xticks(range(0, 361, 60))
        axs[0].set_title("Flux linkage: stator field fixed (i_a = I, i_b = i_c = −I/2), rotor turning", loc="left")
        axs[1].set_title("Torque at 10 A, same test", loc="left")
        fig.tight_layout()
        save(fig, "slice vs full model")
