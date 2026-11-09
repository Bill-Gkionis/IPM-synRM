"""magnetic.py: every number and figure of the note "01 Magnetic".

The note climbs one chain, and this script computes each step of it:

  A. a current makes a field             -> a wire, a loop, an electron's orbit, its spin
  B. many aligned moments                -> iron's magnetisation, then silicon steel's
  C. what FEMM knows about M-19          -> its B-H curve as J = B - mu0 H and mu_r(B),
                                            checked against FEMM's other M-19 and the maker's curve
  D. hysteresis and its loss             -> the coursework's Steinmetz exercise refitted,
                                            next to the steel maker's loss figures
  E. temperature and the Curie point     -> Weiss's molecular field: m = tanh(m / t)
  F. the magnet                          -> N45's straight line, (BH)max, the knee when hot,
                                            and where our motor's magnets work

Run from the Python to FEMM folder, after fields.py and current.py (this script reads their JSON):
    $env:MOTOR_FIGURES_DIR = "out/materials/figures"      # keeps the figures out of the docs figures folder
    uv run python materials/magnetic.py
Outputs: out/materials/magnetic.json, magnetic_m19_curve.csv, and five figures (also copied into the folder
MOTOR_FIGURES_DIR names, and into the repository's docs figures folder
only if it says `docs`: see common.save()).

FEMM's library: section C reads FEMM 5.0's library (common.MATLIB) for the ASM "M-19" curve. FEMM 5.0 has no
"M-19 Steel" any more, so that block is FEMM 4.2's, which motorgeo.M19_STEEL_FEMM42 keeps (common.m19_steel_block()).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import common as cm

N = {}                                   # every number the note quotes, saved as JSON at the end
plt = cm.style()

# ============================================================================ A. a current makes a field
# The smallest piece: a long straight wire. Ampere's law around a circle of radius r gives
# B * 2 pi r = mu0 * I, so B = mu0 I / (2 pi r).
I, r = 10.0, 0.01                         # our FEMM test current (10 A) at 1 cm from the wire
B_wire = cm.MU0 * I / (2 * math.pi * r)
B_earth = 50e-6                           # Earth's field is roughly 25-65 uT; 50 uT is the usual round number
N["A1_wire"] = {"I_A": I, "r_m": r, "B_T": B_wire, "times_earth": B_wire / B_earth}

# Bend the wire into a loop of radius R: at its centre B = mu0 I / (2R); seen from far away it acts
# like a small magnet of moment m = I * A (current times area, in A m^2).
R = 0.01
N["A2_loop"] = {"R_m": R, "B_centre_T": cm.MU0 * I / (2 * R), "m_Am2": I * math.pi * R ** 2}

# An electron going round a nucleus is such a loop. In Bohr's hydrogen atom the orbit has radius a0
# and speed v = alpha * c. Its current is the charge passing per second: I = e v / (2 pi a0).
v = cm.ALPHA * cm.C
I_orbit = cm.E * v / (2 * math.pi * cm.A0)
m_orbit = I_orbit * math.pi * cm.A0 ** 2               # = (e / 2 me) * L with L = me v a0 = hbar
L_orbit = cm.ME * v * cm.A0
N["A3_orbit"] = {"v_m_per_s": v, "I_A": I_orbit, "m_Am2": m_orbit, "L_over_hbar": L_orbit / cm.HBAR,
                 "mu_B": cm.MU_B, "m_over_muB": m_orbit / cm.MU_B}

# ============================================================================ B. many aligned moments
# Iron: bcc, 2 atoms per cube of side a. Each atom carries about 2.22 Bohr magnetons in the metal
# (not 4: the 3d electrons are shared in a band). Line them all up: M = n * m.
a_fe = 0.28665e-9                          # m, bcc iron at room temperature        [CHECK report A, A12]
n_fe = 2 / a_fe ** 3                       # atoms per m^3
m_fe = 2.22 * cm.MU_B                      # A m^2 per atom at 0 K                    [CHECK report A, A12]
Ms0 = n_fe * m_fe
Ms_rt = 1.71e6                             # A/m at room temperature (Kittel, Table 11.? )  [CHECK report A, A12]
N["B1_iron"] = {"a_m": a_fe, "n_per_m3": n_fe, "moment_muB": 2.22, "Ms_0K_A_per_m": Ms0,
                "mu0Ms_0K_T": cm.MU0 * Ms0, "Ms_RT_A_per_m": Ms_rt, "mu0Ms_RT_T": cm.MU0 * Ms_rt}

# Silicon steel: Si atoms take iron's place in the lattice and carry no moment, so they dilute it.
# M-19 has about 3 % silicon by weight (plus a little aluminium).          [CHECK report A, A13/A17]
wt_si = 3.2
x_si = (wt_si / 28.086) / (wt_si / 28.086 + (100 - wt_si) / 55.845)       # atomic fraction
Js_est = cm.MU0 * Ms_rt * (1 - x_si)
N["B2_silicon_steel"] = {"wt_pct_Si": wt_si, "atomic_fraction_Si": x_si, "Js_estimate_T": Js_est}

# ============================================================================ C. what FEMM knows about M-19
blocks = cm.read_fem_blocks()
m19 = blocks["M-19 Steel"]
B, H = m19["bh"][:, 0], m19["bh"][:, 1]
J = B - cm.MU0 * H                           # polarisation: what the material adds on top of mu0 H
with np.errstate(divide="ignore", invalid="ignore"):
    mur = np.where(H > 0, B / (cm.MU0 * H), np.nan)
k_max = np.nanargmax(mur)

def h_at(b, BB=B, HH=H):
    """H [A/m] at a given B on a curve, interpolated linearly in log(H) (the curve's natural scale)."""
    ok = HH > 0
    return float(np.exp(np.interp(b, BB[ok], np.log(HH[ok]))))

points = [0.5, 1.0, 1.3, 1.5, 1.7, 1.8, 2.0]
N["C1_femm_M19"] = {
    "points": int(len(B)), "B_max_T": float(B[-1]), "H_max_A_per_m": float(H[-1]),
    "J_at_top_T": float(J[-1]), "sigma_MS_per_m": m19["Sigma"], "d_lam_mm": m19["d_lam"],
    "lam_fill": m19["LamFill"], "mu_x_linear": m19["Mu_x"],
    "mur_max": float(mur[k_max]), "mur_max_at_B_T": float(B[k_max]),
    "H_at_B": {f"{b:.1f}T": h_at(b) for b in points},
    "mur_at_B": {f"{b:.1f}T": b / (cm.MU0 * h_at(b)) for b in points},
}
np.savetxt(cm.OUT / "magnetic_m19_curve.csv", np.column_stack([B, H, J, mur]), delimiter=",",
           header="B_T,H_A_per_m,J_T,mu_r", comments="", fmt="%.6g")

# FEMM's library holds a second, independent M-19 curve, digitised from the ASM Metals Handbook.
# (FEMM 4.2's library also held "M-19 Steel", our curve. FEMM 5.0's does not, so that block comes from motorgeo.)
lib = cm.read_matlib()
asm = next(b for b in lib if b["BlockName"] == "M-19" and "Metals Handbook" in b["folder"])
lib_m19 = cm.m19_steel_block()
assert np.allclose(lib_m19["bh"], m19["bh"]), "our .FEM's M-19 is not FEMM 4.2's library M-19 Steel"
Ba, Ha = asm["bh"][:, 0], asm["bh"][:, 1]
N["C2_femm_two_curves"] = {
    "ours_is_library_folder": lib_m19["folder"], "asm_folder": asm["folder"], "asm_sigma_MS_per_m": asm["Sigma"],
    "H_ratio_asm_over_ours": {f"{b:.1f}T": h_at(b, Ba, Ha) / h_at(b) for b in [0.5, 1.0, 1.3, 1.5, 1.7]},
}

# The steel maker's own DC curve, when the cross-check has transcribed it (research/m19_dc_maker.csv:
# columns B_T,H_A_per_m). Until then the figure shows FEMM's two curves only.
maker_csv = Path(__file__).parent / "research" / "m19_dc_maker.csv"
maker = np.loadtxt(maker_csv, delimiter=",", skiprows=1) if maker_csv.exists() else None
if maker is not None:
    N["C3_maker_vs_femm"] = {f"{b:.2f}T": {"H_maker": float(h), "H_femm": h_at(b), "ratio_femm_over_maker": h_at(b) / h}
                             for b, h in maker}

# --- figure: the B-H curve
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.semilogx(H[1:], B[1:], color=cm.CAT[0], label="FEMM 'M-19 Steel' (our model)")
ax.semilogx(Ha[1:], Ba[1:], color=cm.CAT[1], ls="--", lw=1.5, label="FEMM 'M-19', ASM Metals Handbook")
if maker is not None:
    ax.semilogx(maker[:, 1], maker[:, 0], "o", color=cm.CAT[2], ms=6, label="steel maker's DC curve")
ax.semilogx(H[1:], J[1:], color=cm.MUTED, lw=1.2, label="J = B - mu0 H (our model)")
ax.axhline(2.0, color=cm.AXIS, lw=0.8)
cm.label_end(ax, 12, 2.0, "J levels off near 2.0 T: saturation", dy=6)
ax.set_xlim(10, 3e5); ax.set_ylim(0, 2.45)
ax.set_xlabel("H, field strength [A/m] (log scale)"); ax.set_ylabel("B, flux density [T]")
ax.set_title("M-19: our FEMM curve is the maker's; FEMM's other one is not")
ax.legend(loc="lower right")
cm.save(fig, "materials-m19-bh.png"); plt.close(fig)

# --- figure: relative permeability against B
fig, ax = plt.subplots(figsize=(7.2, 4.0))
ok = np.isfinite(mur)
ax.semilogy(B[ok], mur[ok], color=cm.CAT[0])
ax.plot(B[k_max], mur[k_max], "o", color=cm.CAT[0])
cm.label_end(ax, B[k_max], mur[k_max], f"peak {mur[k_max]:,.0f} at {B[k_max]:.2f} T", dy=10)
for b in (1.3, 1.8, 2.0):
    mu_b = b / (cm.MU0 * h_at(b))
    ax.plot(b, mu_b, "o", color=cm.INK2, ms=4)
    cm.label_end(ax, b, mu_b, f"{mu_b:,.0f} at {b:.1f} T", dx=-6, dy=-10, ha="right")
ax.set_xlim(0, 2.3); ax.set_ylim(20, 2e4)
ax.set_xlabel("B [T]"); ax.set_ylabel("mu_r = B / (mu0 H) (log scale)")
ax.set_title("How much M-19 multiplies the field, as it fills up")
cm.save(fig, "materials-m19-permeability.png"); plt.close(fig)

# ============================================================================ D. hysteresis and its loss
# The coursework's Exercise 3 (lecture 10, P80-P84): a soft-iron core magnetised at 50 Hz.
# Steinmetz: energy per cycle E_h = k_h * Bmax^n  ->  ln E_h = n ln Bmax + ln k_h, a straight line.
ex3 = np.array([[40, 0.06, 3], [60, 0.17, 16], [80, 0.36, 60], [160, 0.94, 270], [800, 1.51, 580]])
n_fit, lnk_fit = np.polyfit(np.log(ex3[:, 1]), np.log(ex3[:, 2]), 1)
k_fit = math.exp(lnk_fit)
P_ex3 = k_fit * 1.1 ** n_fit * 60
N["D1_steinmetz_ex3"] = {"n_least_squares": n_fit, "kh_least_squares": k_fit, "P_1p1T_60Hz_W_per_m3": P_ex3,
                         "course_n": 1.67, "course_kh": 317.35, "course_P_W_per_m3": 22326.2}

# The steel maker's figures, per cycle and per cubic metre, for comparison (Cleveland-Cliffs DI-MAX,
# June 2023: table 5 p.7 and table 8 p.8). 1 lb = 0.45359237 kg; density 7650 kg/m^3.
LB, RHO = 0.45359237, 7650.0
m15_h_1T = 0.0068 / LB * RHO                       # M-15 hysteresis loss at 1.0 T [J/m^3 per cycle]
m19_total_1p5_60 = {g: wlb / LB for g, wlb in {"29ga_0.356mm": 1.55, "26ga_0.470mm": 1.65, "24ga_0.635mm": 2.00}.items()}
N["D2_maker_losses"] = {
    "M15_hysteresis_1T_J_per_m3_cycle": m15_h_1T,
    "ex3_soft_iron_hysteresis_1T_J_per_m3_cycle": k_fit,
    "M19_max_total_W_per_kg_1p5T_60Hz": m19_total_1p5_60,
    "M19_29ga_total_J_per_m3_cycle_1p5T": m19_total_1p5_60["29ga_0.356mm"] / 60 * RHO,
    "ex3_soft_iron_hysteresis_1p5T_J_per_m3_cycle": k_fit * 1.5 ** n_fit,
}

# M-19's own hysteresis per cycle, separated from the maker's typical loss by current.py (run it first).
import json
sep = json.loads((cm.OUT / "current.json").read_text())["F2_loss_separation_29ga_typical"]
m19_h = {1.0: sep["1.0T"]["hysteresis_J_per_m3_cycle"], 1.5: sep["1.5T"]["hysteresis_J_per_m3_cycle"]}
N["D2_maker_losses"]["M19_29ga_hysteresis_J_per_m3_cycle (current.py, F2)"] = {f"{b}T": h for b, h in m19_h.items()}
N["D2_maker_losses"]["M19_over_ex3_soft_iron"] = {f"{b}T": h / (k_fit * b ** n_fit) for b, h in m19_h.items()}

fig, ax = plt.subplots(figsize=(6.6, 4.2))
bb = np.linspace(0.05, 1.6, 100)
ax.loglog(ex3[:, 1], ex3[:, 2], "o", color=cm.CAT[0], ms=6, label="coursework Ex.3 data (soft iron)")
ax.loglog(bb, k_fit * bb ** n_fit, color=cm.CAT[0], lw=1.4, label=f"Steinmetz fit: {k_fit:.0f} * B^{n_fit:.2f}")
ax.loglog(list(m19_h), list(m19_h.values()), "D", color=cm.CAT[1], ms=6,
          label="M-19 (29 ga), separated from the maker's typical loss")
ax.set_xlabel("Bmax [T] (log)"); ax.set_ylabel("energy lost per cycle [J/m^3] (log)")
ax.set_title("Hysteresis: the loop's area, cycle after cycle")
ax.legend(loc="upper left")
cm.save(fig, "materials-steinmetz.png"); plt.close(fig)

# ============================================================================ E. temperature: the molecular field
# A moment in a field B has two states (spin up or down) whose populations follow Boltzmann's factor;
# the average is m = tanh(mu B / kT). Weiss: the field that aligns each moment comes mostly from its
# neighbours and is proportional to M itself, so m = tanh(m * Tc / T). Solve it by iteration.
def weiss(t, classical=False):
    """Spontaneous magnetisation m = M/Ms at reduced temperature t = T/Tc (spin 1/2, or classical)."""
    if t >= 1:
        return 0.0
    m = 1.0
    for _ in range(2000):
        x = m / t if not classical else 3 * m / t
        m_new = math.tanh(x) if not classical else (1 / math.tanh(x) - 1 / x)
        if abs(m_new - m) < 1e-12:
            break
        m = m_new
    return m

def coeff(Tc_K, T1=20.0, T2=80.0, classical=False):
    """Average % change of m per kelvin between T1 and T2 (deg C), as datasheets quote it."""
    m1, m2 = weiss((T1 + 273.15) / Tc_K, classical), weiss((T2 + 273.15) / Tc_K, classical)
    return 100 * (m2 - m1) / m1 / (T2 - T1)

Tc_fe, Tc_ndfeb = 770 + 273.15, 310 + 273.15          # iron (coursework L10 P44); N45 (Arnold)
N["E1_weiss"] = {
    "iron_coeff_pct_per_K_spin_half": coeff(Tc_fe), "iron_coeff_pct_per_K_classical": coeff(Tc_fe, classical=True),
    "N45_coeff_pct_per_K_spin_half": coeff(Tc_ndfeb), "N45_coeff_pct_per_K_classical": coeff(Tc_ndfeb, classical=True),
    "N45_datasheet_pct_per_K": -0.12,
    "t_motor_iron": [(20 + 273.15) / Tc_fe, (120 + 273.15) / Tc_fe],
    "t_motor_N45": [(20 + 273.15) / Tc_ndfeb, (120 + 273.15) / Tc_ndfeb],
    "m_N45_20C_120C": [weiss((20 + 273.15) / Tc_ndfeb), weiss((120 + 273.15) / Tc_ndfeb)],
    "m_iron_20C": weiss((20 + 273.15) / Tc_fe),
}

# The coursework's Exercise 4 (lecture 10, P85-P89): nickel above its Curie point follows
# chi = C / (T - Tc). Rearranged, T = C * (1/chi) + Tc: a straight line in 1/chi.
T4 = np.array([800, 900, 1000, 1100, 1200.0]); chi4 = np.array([3.30, 2.10, 1.55, 1.20, 1.00]) * 1e-5
C4, Tc4 = np.polyfit(1 / chi4, T4, 1)
N["E2_curie_weiss_ex4"] = {"C": C4, "Tc_K": Tc4, "Tc_C": Tc4 - 273.15, "course_Tc_K": 628.03, "course_C": 5.7e-3}

fig, ax = plt.subplots(figsize=(6.8, 4.2))
tt = np.linspace(0.005, 1.2, 240)
ax.plot(tt, [weiss(t) for t in tt], color=cm.CAT[0], label="Weiss, spin 1/2: m = tanh(m/t)")
ax.plot(tt, [weiss(t, True) for t in tt], color=cm.CAT[1], ls="--", lw=1.5, label="Weiss, classical moments")
for name, Tc in (("M-19 (iron)", Tc_fe), ("N45", Tc_ndfeb)):
    t1, t2 = (20 + 273.15) / Tc, (120 + 273.15) / Tc
    ax.axvspan(t1, t2, color=cm.AXIS, alpha=0.35, lw=0)
    ax.text((t1 + t2) / 2, 0.08, f"{name}\nat 20-120 °C", ha="center", fontsize=8.5, color=cm.INK2)
ax.set_xlim(0, 1.2); ax.set_ylim(0, 1.05)
ax.set_xlabel("t = T / Tc (temperature over the Curie point)"); ax.set_ylabel("m = M / Ms")
ax.set_title("Heat against alignment: the motor's two magnetic materials")
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.18))
cm.save(fig, "materials-weiss-curie.png"); plt.close(fig)

# ============================================================================ F. the magnet: N45
n45 = blocks["N45"]
mur_rec, Hc_femm = n45["Mu_x"], n45["H_c"]
Br_femm = cm.MU0 * mur_rec * Hc_femm                 # FEMM's straight line crosses H = 0 here
# Arnold Magnetic Technologies, N45 datasheet Rev. 210607
Br_nom, Br_min, Br_max = 1.35, 1.32, 1.38            # T
HcJ_min, HcB_nom = 955e3, 955e3                      # A/m
BHmax_nom, BHmax_min, BHmax_max = 350e3, 334e3, 366e3
a_Br, a_HcJ = -0.12e-2, -0.62e-2                     # per K, measured 20-80 deg C
BHmax_line = Br_nom ** 2 / (4 * cm.MU0 * mur_rec)    # the best B*H on a straight recoil line is at B = Br/2
BHmax_femm = Br_femm ** 2 / (4 * cm.MU0 * mur_rec)   # the same for FEMM's line (Br 1.359 T): like for like

def Br_at(T):  return Br_nom * (1 + a_Br * (T - 20))
def HcJ_at(T): return HcJ_min * (1 + a_HcJ * (T - 20))
def knee_B(T):  return Br_at(T) - cm.MU0 * mur_rec * HcJ_at(T)    # B where the intrinsic curve drops

temps = [20, 60, 80, 100, 120]
N["F1_N45_femm_vs_datasheet"] = {
    "mur_rec_femm": mur_rec, "Hc_femm_A_per_m": Hc_femm, "Br_femm_T": Br_femm,
    "Br_datasheet_T": [Br_min, Br_nom, Br_max], "sigma_femm_MS_per_m": n45["Sigma"],
    "resistivity_datasheet_uOhm_cm_parallel": 150, "sigma_datasheet_MS_per_m": 1 / 150e-8 / 1e6,
    "BHmax_straight_line_kJ_per_m3": BHmax_line / 1e3, "BHmax_datasheet_kJ_per_m3": [BHmax_min / 1e3, BHmax_nom / 1e3, BHmax_max / 1e3],
}
N["F2_N45_hot"] = {f"{T}C": {"Br_T": Br_at(T), "HcJ_kA_per_m": HcJ_at(T) / 1e3, "knee_B_T": knee_B(T)} for T in temps}
# 100 and 120 deg C lie past Arnold's measured 20-80 deg C: there the knee is a straight-line extrapolation, and
# the makers disagree on how fast HcJ falls (report A, A21): Eclipse -0.60, Arnold -0.62, DuraMag -0.69 %/K.
N["F2_knee_spread_extrapolated"] = {
    f"{T}C": {name: Br_at(T) - cm.MU0 * mur_rec * HcJ_min * (1 + a * (T - 20))
              for name, a in {"Eclipse -0.60": -0.60e-2, "Arnold -0.62": -0.62e-2, "DuraMag -0.69": -0.69e-2}.items()}
    for T in (100, 120)}

# Where our magnets work. fields.py measured B along each magnet's own direction in the solved sectors of
# the air-band study. With no current, the median sits on the load line of the circuit itself: on the
# recoil line that fixes H, and B / (mu0 |H|) is the permeance coefficient Pc, the load line's slope.
# With 10 A at the best angle, the coils' field moves the point a little; those are shown as markers.
import json
fields = json.loads((cm.OUT / "fields.json").read_text())      # run fields.py first
def operating(Bm):
    Hm = (Bm - Br_femm) / (cm.MU0 * mur_rec)
    return Hm, Bm / (cm.MU0 * abs(Hm))

# The textbook circuit: magnet thickness lm in series with the gap g, flux concentrated from the two
# magnets of a pole (2 x 9.08 mm wide) into one pole arc at the gap.
lm, g = 1.57e-3, 0.8e-3
pole_arc = 2 * math.pi * 45.4e-3 / 8
Pc_simple = (lm / g) * (pole_arc / (2 * 9.08e-3))
op = {}
for name, case in (("band", "with_band"), ("filled", "no_band")):
    Bm = fields[f"noload_{case}"]["B_parallel_median_T"]
    B10 = fields[f"mtpa_{case}"]["B_parallel_median_T"]
    Hm, Pc = operating(Bm)
    hot = {}
    for T in temps:           # the same load line, hot: the line meets the hotter recoil line lower down
        Hh = -Br_at(T) / (cm.MU0 * (Pc + mur_rec))
        hot[f"{T}C"] = {"B_T": -cm.MU0 * Pc * Hh, "H_kA_per_m": Hh / 1e3,
                        "margin_to_knee_kA_per_m": (HcJ_at(T) + Hh) / 1e3}
    op[name] = {"B_T": Bm, "H_kA_per_m": Hm / 1e3, "Pc": Pc, "BH_kJ_per_m3": Bm * abs(Hm) / 1e3,
                "B_over_Br": Bm / Br_femm, "BH_over_BHmax": Bm * abs(Hm) / BHmax_femm, "hot": hot,
                "at_10A_B_T": B10, "at_10A_H_kA_per_m": operating(B10)[0] / 1e3}
N["F3_operating_points"] = {"Pc_simple_circuit": Pc_simple, **op}

# FEMM fills each magnet's slot (1.77 x 9.28 mm) with N45; the real magnet is 1.57 x 9.08 mm.
Pc_slot = (1.77e-3 / g) * (pole_arc / (2 * 9.28e-3))
Pc_real = op["band"]["Pc"] * Pc_simple / Pc_slot                 # FEMM's no-load Pc, scaled to the real magnet
N["F3_operating_points"]["Pc_simple_slot_filling_magnet"] = Pc_slot
N["F3_operating_points"]["real_magnet_noload_band"] = {"Pc": Pc_real, "B_T": Br_femm * Pc_real / (Pc_real + mur_rec)}

# The weakest triangle (10 A, best angle, cold field) against the knee. As the magnet heats, its B falls with
# Br (same circuit), while the knee climbs: they meet at the crossing temperature.
weak = fields["mtpa_with_band"]["weakest_element"]["B_parallel_T"]
cdf = fields["mtpa_with_band"]["area_pct_below_B"]
tb, pb = np.array([float(k) for k in cdf]), np.array(list(cdf.values()))
lo, hi = 20.0, 200.0
for _ in range(60):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if weak * Br_at(mid) / Br_at(20) > knee_B(mid) else (lo, mid)
N["F4_weakest_spot"] = {
    "B_min_cold_10A_T": weak, "element": fields["mtpa_with_band"]["weakest_element"],
    "crossing_C_if_B_falls_with_Br": lo,
    "area_pct_below_knee_cold_field": {f"{T}C": float(np.interp(knee_B(T), tb, pb)) for T in temps},
    "area_pct_below_120C_knee_by_maker": {k: float(np.interp(b, tb, pb))
                                          for k, b in N["F2_knee_spread_extrapolated"]["120C"].items()},
}

# --- figure: the second quadrant, hot and cold
fig, ax = plt.subplots(figsize=(7.4, 4.8))
ramp = [cm.BLUE_RAMP[k] for k in (250, 350, 450, 550, 650)]
Hgrid = np.linspace(-1300e3, 0, 600)
for T, col in zip(temps, ramp):
    br, hcj = Br_at(T), HcJ_at(T)
    Jcurve = np.where(Hgrid > -hcj, br + cm.MU0 * (mur_rec - 1) * Hgrid, np.nan)   # intrinsic: flat, then drops
    Bline = br + cm.MU0 * mur_rec * Hgrid
    Bcurve = np.where(Hgrid > -hcj, Bline, np.nan)
    ax.plot(Hgrid / 1e3, Bcurve, color=col, lw=1.8)
    ax.plot([-hcj / 1e3, -hcj / 1e3], [knee_B(T), 0], color=col, lw=1.0, ls=":")
    ax.plot(-hcj / 1e3, knee_B(T), "o", color=col, ms=5)
    cm.label_end(ax, -hcj / 1e3, knee_B(T), f"{T} °C knee", dx=-6, dy=8 if T == 20 else 0, ha="right")
# FEMM's N45 is the 20 deg C line carried on past the knee, all the way to B = 0 at its Hc
Hf = np.linspace(-Hc_femm, -HcJ_at(20), 20)
ax.plot(Hf / 1e3, Br_femm + cm.MU0 * mur_rec * Hf, color=cm.MUTED, lw=1.2, ls="--")
ax.text(-1285, 0.025, "grey dashes: FEMM's N45 runs on to B = 0", fontsize=8.5, color=cm.INK2)
Hl = np.linspace(-500e3, 0, 50)
for name, colour in (("band", cm.CAT[1]), ("filled", cm.CAT[2])):
    Pc = op[name]["Pc"]
    label = (f"our magnets (air band): {op[name]['B_T']:.2f} T, Pc {Pc:.1f}" if name == "band" else
             f"band filled with steel: {op[name]['B_T']:.2f} T, Pc {Pc:.1f}")
    ax.plot(Hl / 1e3, -cm.MU0 * Pc * Hl, color=colour, lw=1.2, ls="--")
    ax.plot([op[name]["H_kA_per_m"]], [op[name]["B_T"]], color=colour, lw=1.2, ls="--", marker="D", ms=6,
            label=label)
    ax.plot([op[name]["at_10A_H_kA_per_m"]], [op[name]["at_10A_B_T"]], "D", ms=6, mfc=cm.SURFACE, mec=colour, mew=1.4)
ax.plot([], [], "D", ms=6, mfc=cm.SURFACE, mec=cm.INK2, mew=1.4, label="the same magnets with 10 A at the best angle")
ax.legend(loc="upper left", title="load lines, no current, cold (FEMM medians)", title_fontsize=9, alignment="left")
ax.text(-1280, 1.08, "recoil lines at 20, 60, 80, 100 and 120 °C,\nlight to dark; each ends at its knee",
        fontsize=8.5, color=cm.INK2, va="top")
ax.set_xlim(-1300, 0); ax.set_ylim(0, 1.5)
ax.set_xlabel("H inside the magnet [kA/m]"); ax.set_ylabel("B [T]")
ax.set_title("N45, second quadrant: the knee climbs as the magnet heats")
cm.save(fig, "materials-n45-demag.png"); plt.close(fig)

# ============================================================================ write and print
cm.write_numbers("magnetic", N)
c = N["C1_femm_M19"]
print(f"wire 10 A at 1 cm: {B_wire * 1e3:.3f} mT ({B_wire / B_earth:.1f}x Earth)")
print(f"Bohr orbit: I = {I_orbit * 1e3:.3f} mA, m = {m_orbit:.4e} A m^2 = {m_orbit / cm.MU_B:.4f} mu_B, L = {L_orbit / cm.HBAR:.4f} hbar")
print(f"iron: n = {n_fe:.4e}/m^3, Ms(0K) = {Ms0:.4e} A/m, mu0 Ms = {cm.MU0 * Ms0:.3f} T (0 K), {cm.MU0 * Ms_rt:.3f} T (RT)")
print(f"Fe-{wt_si}%Si: x_Si = {x_si:.4f}, Js estimate = {Js_est:.3f} T;  FEMM M-19 J at {B[-1]} T = {J[-1]:.3f} T")
print(f"M-19 mu_r max = {c['mur_max']:.0f} at {c['mur_max_at_B_T']:.2f} T;  mu_r at points: "
      + ", ".join(f"{k} {v:.0f}" for k, v in c["mur_at_B"].items()))
print("FEMM ASM M-19 / ours, H ratio: " + ", ".join(f"{k} {v:.2f}" for k, v in N["C2_femm_two_curves"]["H_ratio_asm_over_ours"].items()))
print(f"Steinmetz Ex3 (least squares): n = {n_fit:.3f}, kh = {k_fit:.1f}, P(1.1 T, 60 Hz) = {P_ex3:.0f} W/m^3")
print("M-19 hysteresis per cycle: " + ", ".join(f"{b} T {h:.0f} J/m^3 ({h / (k_fit * b ** n_fit):.2f} of Ex.3's soft iron)"
                                                 for b, h in m19_h.items()))
print(f"M-15 hysteresis at 1 T = {m15_h_1T:.1f} J/m^3/cycle; M-19 29ga total limit at 1.5 T/60 Hz = "
      f"{N['D2_maker_losses']['M19_29ga_total_J_per_m3_cycle_1p5T']:.1f} J/m^3/cycle")
e1 = N["E1_weiss"]
print(f"Weiss 20->80 C: iron {e1['iron_coeff_pct_per_K_spin_half']:.4f} %/K (spin 1/2); "
      f"N45 {e1['N45_coeff_pct_per_K_spin_half']:.3f} (spin 1/2) / {e1['N45_coeff_pct_per_K_classical']:.3f} (classical) %/K vs -0.12")
print(f"Curie-Weiss Ex4 refit: Tc = {Tc4:.1f} K, C = {C4:.3e}")
f1 = N["F1_N45_femm_vs_datasheet"]
print(f"N45: FEMM Br = {Br_femm:.4f} T; straight-line (BH)max = {BHmax_line / 1e3:.0f} kJ/m^3 (datasheet 334-366)")
for T in temps:
    h = N["F2_N45_hot"][f"{T}C"]
    print(f"  {T:>3} C: Br {h['Br_T']:.3f} T, HcJ {h['HcJ_kA_per_m']:.0f} kA/m, knee at B = {h['knee_B_T']:.3f} T")
print(f"Pc simple circuit = {Pc_simple:.2f}")
for name in ("band", "filled"):
    o = op[name]
    print(f"  {name}: B {o['B_T']:.2f} T, H {o['H_kA_per_m']:.0f} kA/m, Pc {o['Pc']:.2f}, BH/BHmax {o['BH_over_BHmax']:.2f}; "
          + "margin to knee: " + ", ".join(f"{k} {v['margin_to_knee_kA_per_m']:.0f}" for k, v in o["hot"].items()))
for T, v in N["F2_knee_spread_extrapolated"].items():
    print(f"  knee at {T} (extrapolated): " + ", ".join(f"{k} {b:.3f} T" for k, b in v.items()))
w = N["E1_weiss"]
print(f"Weiss m: N45 {w['m_N45_20C_120C'][0]:.3f} at 20 C, {w['m_N45_20C_120C'][1]:.3f} at 120 C; iron {w['m_iron_20C']:.4f} at 20 C")
o = N["F3_operating_points"]
print(f"Pc by hand: real magnet {o['Pc_simple_circuit']:.2f}, slot-filling {o['Pc_simple_slot_filling_magnet']:.2f}; real magnet at no load: "
      f"Pc {o['real_magnet_noload_band']['Pc']:.2f}, B {o['real_magnet_noload_band']['B_T']:.3f} T; BH/BHmax (band) {o['band']['BH_over_BHmax']:.3f}")
f4 = N["F4_weakest_spot"]
print(f"weakest element {f4['B_min_cold_10A_T']:.3f} T (area {f4['element']['area_mm2']:.5f} mm^2 at {f4['element']['centroid_mm']}); "
      f"meets the knee at {f4['crossing_C_if_B_falls_with_Br']:.1f} C; area below the knee, cold field: "
      + ", ".join(f"{k} {v:.3f} %" for k, v in f4["area_pct_below_knee_cold_field"].items()))
print("area below the 120 C knee, by maker: " + ", ".join(f"{k} {v:.3f} %" for k, v in N["F4_weakest_spot"]["area_pct_below_120C_knee_by_maker"].items()))
