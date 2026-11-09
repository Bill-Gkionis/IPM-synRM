"""Power factor and field harmonics of the reference rotor, from FEMM output already on disk (no new solves).

Sources: out/barrier (air-band study: phase1/phase2 CSVs, air-gap B samples), out/barrier/reference_open.csv,
out/flux_explorer/states.csv. Prints the numbers for the project notes and writes its figures.
Run: uv run python pf_harmonics.py <figure folder>   (from this folder)
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import motorgeo as mg

SRC = os.path.abspath("out/barrier")
DST = sys.argv[1] if len(sys.argv) > 1 else os.path.abspath("out/pf_harmonics")
PREFIX = "Power factor and harmonics - "
P = mg.Params()
PP = P.pole_pairs
W_E = PP * 1000 * 2 * math.pi / 60          # electrical rad/s at 1000 rpm
COL = {"with band": "#2a78d6", "no band": "#eb6834"}
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlesize": 10, "axes.titlecolor": INK, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "lines.linewidth": 1.8})


def save(fig, name):
    os.makedirs(DST, exist_ok=True)
    fig.savefig(os.path.join(DST, PREFIX + name + ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def spectrum(y):
    """Amplitude of each harmonic of one evenly sampled period (index = harmonic order)."""
    return np.abs(np.fft.rfft(y)) / len(y) * np.r_[1, np.full(len(y) // 2, 2.0)][: len(y) // 2 + 1]


def thd(a, keep=None):
    orders = np.arange(2, len(a)) if keep is None else np.array([n for n in range(2, len(a)) if keep(n)])
    return 100 * math.sqrt((a[orders] ** 2).sum()) / a[1]


def pf(psi_d, psi_q, i_d, i_q):
    """Internal power factor: angle between voltage (= w x psi, turned 90 deg) and current; R and iron loss ignored."""
    return (psi_d * i_q - psi_q * i_d) / (np.hypot(psi_d, psi_q) * np.hypot(i_d, i_q))


# ================================================================ 1. power factor from psi_d, psi_q
ph1 = pd.read_csv(os.path.join(SRC, "phase1.csv"))
ph2 = pd.read_csv(os.path.join(SRC, "phase2.csv"))
g = ph1[ph1.study == "gamma"].groupby(["variant", "I", "gamma"])[["psi_d", "psi_q", "id", "iq", "T_age"]].mean().reset_index()
g["pf"] = pf(g.psi_d, g.psi_q, g.id, g.iq)
g["T_check"] = 1.5 * PP * g.pf * np.hypot(g.psi_d, g.psi_q) * np.hypot(g.id, g.iq)   # T = 1.5 p |psi||i| PF
print("check T = 1.5 p |psi| |i| PF against FEMM's torque: max diff",
      f"{(g.T_check - g.T_age).abs().max():.4f} N.m (should be ~0: same identity as the dq torque)")

ref = pd.read_csv(os.path.join(SRC, "reference_open.csv"))
r = ref[ref.study == "gamma"].groupby("gamma")[["psi_d", "psi_q", "id", "iq", "T_age"]].mean().reset_index()
r["pf"] = pf(r.psi_d, r.psi_q, r.id, r.iq)

rows = []
for v in ("with band", "no band"):
    for I in (5.0, 10.0, 15.0, 20.0):
        s = g[(g.variant == v) & (g.I == I)].sort_values("gamma")
        k = int(np.argmax(s.T_age.values)); k = min(max(k, 1), len(s) - 2)
        a2, a1, a0 = np.polyfit(s.gamma.values[k - 1:k + 2], s.T_age.values[k - 1:k + 2], 2)
        gm = -a1 / (2 * a2)
        at = lambda col: float(np.interp(gm, s.gamma, s[col]))
        rows.append({"design": v, "I": I, "gamma_mtpa": gm, "T": a0 + a1 * gm + a2 * gm ** 2,
                     "psi_d_mWb": 1e3 * at("psi_d"), "psi_q_mWb": 1e3 * at("psi_q"),
                     "pf_mtpa": pf(at("psi_d"), at("psi_q"), -I * math.sin(math.radians(gm)), I * math.cos(math.radians(gm))),
                     "pf_q_axis": float(s[s.gamma == 0].pf.iloc[0])})
k = int(np.argmax(r.T_age.values)); a2, a1, a0 = np.polyfit(r.gamma.values[k - 1:k + 2], r.T_age.values[k - 1:k + 2], 2)
gm = -a1 / (2 * a2)
rows.append({"design": "magnet inside band", "I": 10.0, "gamma_mtpa": gm, "T": a0 + a1 * gm + a2 * gm ** 2,
             "psi_d_mWb": 1e3 * float(np.interp(gm, r.gamma, r.psi_d)), "psi_q_mWb": 1e3 * float(np.interp(gm, r.gamma, r.psi_q)),
             "pf_mtpa": float(np.interp(gm, r.gamma, r.pf)), "pf_q_axis": float(r[r.gamma == 0].pf.iloc[0])})
mt = pd.DataFrame(rows)
print("\nPower factor at MTPA (and with all-q current):")
print(mt.round(3).to_string(index=False))

# ---- time-domain cross-check at MTPA 10 A: rebuild a full period of psi_a from 0..60 deg of three phases
print("\nTime-domain cross-check at MTPA, 10 A, 1000 rpm (voltage = d psi/dt, no R):")
xcheck = {}
for v in ("with band", "no band"):
    rp = ph2[(ph2.study == "ripple") & (ph2.variant == v)].sort_values("theta_e").iloc[:60]   # 0..59 deg
    pa, pb, pc = rp.psi_a.values, rp.psi_b.values, rp.psi_c.values
    psi = np.concatenate([pa, -pb, pc, -pa, pb, -pc])          # psi_a(th+60)=-psi_b(th), (th+120)=psi_c(th), ...
    th = np.radians(np.arange(360))
    X = np.fft.rfft(psi)
    e = np.fft.irfft(1j * np.arange(len(X)) * X, 360) * W_E   # e = d psi / dt
    i_d, i_q = rp.id.iloc[0], rp.iq.iloc[0]
    i = i_d * np.cos(th) - i_q * np.sin(th)
    P3 = 3 * np.mean(e * i)
    S3 = 3 * np.sqrt(np.mean(e ** 2)) * np.sqrt(np.mean(i ** 2))
    E1 = np.fft.rfft(e)[1]; I1 = np.fft.rfft(i)[1]
    pf_disp = math.cos(np.angle(E1) - np.angle(I1))
    xcheck[v] = (P3 / S3, pf_disp, e, i)
    print(f"  {v}: true PF (P/S, harmonics included) {P3 / S3:.4f} | fundamental-only PF {pf_disp:.4f} | "
          f"power {P3:.2f} W = torque x speed {np.mean(rp.T_age.values) * 1000 * 2 * math.pi / 60:.2f} W")

# ---- every rotor of the flux explorer: PF at its own MTPA, 10 A
fx = pd.read_csv(os.path.join("out", "flux_explorer", "states.csv"))
from flux_explorer_sweep import BASE, SLIDERS, geometries
gx = {gid: x for gid, x in geometries()}
pos = fx.groupby(["gid", "band", "I", "gamma"])[["T_age", "psi_d", "psi_q"]].mean().reset_index()
expl = []
for (gid, band), s in pos[pos.I == 10].groupby(["gid", "band"]):
    s = s.sort_values("gamma")                                  # 15-deg samples; the peak is refined by a parabola
    k = int(np.argmax(s.T_age.values)); k = min(max(k, 1), len(s) - 2)
    a2, a1, a0 = np.polyfit(s.gamma.values[k - 1:k + 2], s.T_age.values[k - 1:k + 2], 2)
    gm = float(np.clip(-a1 / (2 * a2), 0, 90))
    pd_, pq_ = np.interp(gm, s.gamma, s.psi_d), np.interp(gm, s.gamma, s.psi_q)
    key = next((k2 for k2 in SLIDERS if gx[gid][k2] != BASE[k2]), "zero change (reference design)")
    expl.append({"gid": gid, "band": band, "changed": key, "value": gx[gid].get(key, ""), "gamma": gm,
                 "T": a0 + a1 * gm + a2 * gm ** 2, "pf": pf(pd_, pq_, -10 * math.sin(math.radians(gm)), 10 * math.cos(math.radians(gm)))})
ex = pd.DataFrame(expl)
print("\nPF at 10 A MTPA for every rotor of the flux explorer (coarse mesh):")
print(ex.sort_values(["band", "changed", "value"]).round(3).to_string(index=False))

# ================================================================ 2. harmonics
# ---- winding factors and the stator MMF, straight from motorgeo's winding table
table = mg.winding(P)
cond = [(mg.slot_angle(P, k), q, s) for k, lay in table.items() for q, s in lay]
kw = {}
for n in range(1, 26):
    a = [s * np.exp(1j * n * PP * math.radians(ang)) for ang, q, s in cond if q == "A"]
    kw[n] = abs(sum(a)) / len(a)
th = np.radians(np.arange(0, 360, 0.05))
F = np.zeros_like(th)
for ang, q, s in cond:                                          # i_a = 1, i_b = i_c = -1/2: phase A's axis
    F -= s * {"A": 1.0, "B": -0.5, "C": -0.5}[q] * (th > math.radians(ang))
F -= F.mean()
mmf = np.abs(np.fft.rfft(F)) / len(F) * 2
mmf_e = {m / PP: mmf[m] / mmf[PP] for m in range(1, 120) if mmf[m] / mmf[PP] > 0.005}
print("\nwinding factors kw_n:", ", ".join(f"{n}: {kw[n]:.3f}" for n in range(1, 26) if n % 2))
print("stator MMF harmonics (electrical order: % of fundamental):",
      ", ".join(f"{o:g}: {100 * a:.1f}" for o, a in sorted(mmf_e.items()) if o <= 30))

# ---- air-gap B: no load, MTPA, all-q (all orders, even and odd)
gap = {}
for v in ("with_band", "no_band"):
    for case in ("noload", "mtpa", "qaxis"):
        d = pd.read_csv(os.path.join(SRC, f"gap_{case}_{v}.csv"))
        gap[(v, case)] = (d.deg.values, d.br.values, spectrum(d.br.values))
for v in ("with_band", "no_band"):
    for case in ("noload", "mtpa"):
        a = gap[(v, case)][2]
        top = sorted(range(2, 60), key=lambda n: -a[n])[:7]
        print(f"gap B {v:9s} {case:6s}: B1 {a[1]:.3f} T, THD {thd(a[:60]):.1f} %, even-order share "
              f"{100 * math.sqrt((a[2:60:2] ** 2).sum()) / a[1]:.1f} %, largest: " + ", ".join(f"{n}:{a[n]:.3f}" for n in sorted(top)))

# ---- back-EMF per phase vs line to line
emf = {}
for v in ("with band", "no band"):
    nl = ph1[(ph1.study == "noload") & (ph1.variant == v)].sort_values("theta_e")
    def e_of(psi):
        X = np.fft.rfft(psi)
        return np.fft.irfft(1j * np.arange(len(X)) * X, len(psi)) * W_E
    ea, eb = e_of(nl.psi_a.values), e_of(nl.psi_b.values)
    sa, sl = spectrum(ea), spectrum(ea - eb)
    emf[v] = (nl.theta_e.values, ea, ea - eb, sa, sl)
    print(f"EMF {v:9s}: phase peak {sa[1]:.2f} V THD {thd(sa):.1f} % | line-to-line peak {sl[1]:.2f} V (x{sl[1] / sa[1]:.3f}) "
          f"THD {thd(sl):.1f} % | 3rd: phase {sa[3]:.3f} V, line {sl[3]:.4f} V | 17th/19th phase {sa[17]:.3f}/{sa[19]:.3f} V")

# ================================================================ 3. figures
# PF vs current angle
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
for ax, I in zip(axs, (10.0, 20.0)):
    for v in ("with band", "no band"):
        s = g[(g.variant == v) & (g.I == I)].sort_values("gamma")
        ax.plot(s.gamma, s.pf, color=COL[v], label=v)
        m = mt[(mt.design == v) & (mt.I == I)].iloc[0]
        ax.plot(m.gamma_mtpa, m.pf_mtpa, "o", ms=6, color=COL[v], mec=SURF, mew=1.5)
        ax.annotate(f"MTPA {m.pf_mtpa:.2f}", (m.gamma_mtpa, m.pf_mtpa), textcoords="offset points", xytext=(6, -12),
                    fontsize=8, color=INK2)
    if I == 10:
        ax.plot(r.gamma, r.pf, color=MUTED, ls="--", lw=1.4, label="magnet inside band (reference)")
    ax.set_title(f"{I:g} A peak ({I * P.turns:g} A-turns)", loc="left")
    ax.set_xlabel("current angle γ [°] (0 = all q, 90 = all −d)"); ax.set_xticks(range(0, 91, 15)); ax.set_ylim(0, 1)
axs[0].set_ylabel("internal power factor cos φ"); axs[0].legend(loc="lower center")
fig.suptitle("Power factor from FEMM flux linkages (fine mesh, mean of 6 rotor positions). No resistance, no end windings",
             x=0.01, ha="left", fontsize=10.5, color=INK)
fig.tight_layout(); save(fig, "PF vs current angle")

# dq phasor diagrams at MTPA, 10 A
fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
for ax, v in zip(axs, ("with band", "no band")):
    m = mt[(mt.design == v) & (mt.I == 10.0)].iloc[0]
    gm = math.radians(m.gamma_mtpa); i_d, i_q = -10 * math.sin(gm), 10 * math.cos(gm)
    psd, psq = m.psi_d_mWb, m.psi_q_mWb
    k = 10 / 25                                                  # draw mWb on the same scale as amps: 25 mWb = 10 units
    vec = {"current i": (i_d, i_q, INK), "flux linkage ψ": (psd * k, psq * k, COL[v]), "voltage v = jωψ": (-psq * k, psd * k, "#e34948")}
    for name, (x, y, c) in vec.items():
        ax.annotate("", xy=(x, y), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=c, lw=2))
        ax.text(x * 1.06, y * 1.06 + 0.3, name, color=c, fontsize=9, ha="center")
    ax.annotate("", xy=(psd * k, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="-", color=COL[v], lw=1, ls=":"))
    ax.axhline(0, color=AXIS, lw=0.8); ax.axvline(0, color=AXIS, lw=0.8)
    ax.text(10.6, 0.25, "d (magnet)", color=MUTED, fontsize=8, ha="right"); ax.text(0.25, 10.6, "q", color=MUTED, fontsize=8)
    ax.set_xlim(-10.8, 10.8); ax.set_ylim(-1.5, 11.2); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    phi = math.degrees(math.acos(m.pf_mtpa))
    ax.set_title(f"{v}, MTPA 10 A: φ = {phi:.0f}°, cos φ = {m.pf_mtpa:.2f}\n"
                 f"ψd {psd:.1f}, ψq {psq:.1f} mWb · id {i_d:.1f}, iq {i_q:.1f} A", loc="left", fontsize=9.5)
fig.suptitle("The dq picture: voltage leads flux linkage by 90°; power factor is the cosine of the angle between voltage "
             "and current", x=0.01, ha="left", fontsize=10.5, color=INK)
fig.tight_layout(); save(fig, "dq phasors")

# air-gap B spectra, no load vs load (with band)
fig, axs = plt.subplots(2, 1, figsize=(9, 6.4), gridspec_kw={"height_ratios": [1.1, 1]})
for case, c, lab in (("noload", "#2a78d6", "no current"), ("mtpa", "#e34948", "10 A at MTPA")):
    deg, br, a = gap[("with_band", case)]
    axs[0].plot(deg, br, color=c, lw=1.2, label=lab)
    orders = np.arange(1, 41)
    axs[1].bar(orders + (-0.2 if case == "noload" else 0.2), a[orders], width=0.38, color=c, label=lab)
axs[0].axhline(0, color=AXIS, lw=0.8); axs[0].set_xlim(0, 90); axs[0].legend(loc="lower left")
axs[0].set_xlabel("position around the gap [mechanical °] (one N and one S pole)"); axs[0].set_ylabel("radial B [T]")
axs[0].set_title("Air-gap radial flux density, the reference design (with band)", loc="left")
axs[1].set_xticks(range(1, 41, 2)); axs[1].set_xlabel("harmonic order (per electrical cycle)")
axs[1].set_ylabel("amplitude [T]"); axs[1].legend(loc="upper right")
axs[1].set_title("Its harmonics: all orders, even ones included", loc="left")
fig.tight_layout(); save(fig, "air-gap B spectrum")

# winding factors + MMF
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
n_odd = [n for n in range(1, 26) if n % 2]
axs[0].bar(n_odd, [kw[n] for n in n_odd], width=1.2, color="#2a78d6")
for n in (1, 17, 19):
    axs[0].annotate(f"{kw[n]:.3f}", (n, kw[n]), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8, color=INK2)
axs[0].set_xticks(n_odd); axs[0].set_ylim(0, 1.05); axs[0].set_xlabel("harmonic order n")
axs[0].set_ylabel("winding factor k_w,n"); axs[0].set_title("What the winding lets through (phase A)", loc="left")
ords = [o for o in sorted(mmf_e) if o <= 30]
axs[1].bar(ords, [100 * mmf_e[o] for o in ords], width=0.6, color="#eb6834")
axs[1].set_xlabel("harmonic order (electrical)"); axs[1].set_ylabel("% of fundamental")
axs[1].set_title("Stator MMF with i_a = I, i_b = i_c = −I/2 (current sheet steps)", loc="left")
fig.tight_layout(); save(fig, "winding factors and MMF")

# EMF phase vs line
fig, axs = plt.subplots(1, 2, figsize=(11, 3.7), gridspec_kw={"width_ratios": [1.3, 1]})
thv, ea, eab, sa, sl = emf["with band"]
axs[0].plot(thv, ea, color="#2a78d6", label=f"phase A (THD {thd(sa):.1f} %)")
axs[0].plot(thv, eab, color="#eb6834", label=f"line A−B (THD {thd(sl):.1f} %)")
axs[0].set_xlim(0, 360); axs[0].set_xticks(range(0, 361, 60)); axs[0].set_xlabel("rotor position [electrical °]")
axs[0].set_ylabel("back-EMF at 1000 rpm [V]"); axs[0].legend(loc="upper left"); axs[0].axhline(0, color=AXIS, lw=0.8)
axs[0].set_title("Reference design, no current", loc="left")
orders = np.arange(1, 22, 2)
axs[1].bar(orders - 0.22, sa[orders] / sa[1] * 100, width=0.42, color="#2a78d6", label="phase")
axs[1].bar(orders + 0.22, sl[orders] / sl[1] * 100, width=0.42, color="#eb6834", label="line to line")
axs[1].set_xticks(orders); axs[1].set_ylim(0, 15); axs[1].set_ylabel("% of fundamental")
axs[1].set_xlabel("harmonic order"); axs[1].legend(loc="upper right")
axs[1].set_title("The 3rd, 9th, 15th… vanish between lines", loc="left")
fig.tight_layout(); save(fig, "back-EMF phase vs line")
print("\nfigures in", DST)
